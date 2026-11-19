#!/usr/bin/env python3
"""
joint_genotype.py — all-sites joint genotyping over the callable CDS targets.

One chromosome per invocation (SLURM array); ``--merge`` concatenates.

The one thing this script exists to get right
---------------------------------------------
**Invariant sites are emitted, not just variants.** ``bcftools call -v`` gives a
variants-only VCF, in which a monomorphic site and an uncallable site are
indistinguishable — both are simply absent. piN/piS is a ratio of per-site
rates, so its denominator is the number of sites actually *observed*, and a
variants-only VCF silently replaces that denominator with a guess. The failure
mode is not noise: regions with poorer callability get smaller denominators,
which inflates their apparent diversity. Eliminated regions have poorer
callability. That is the artifact this whole pipeline is built to avoid, and it
walks straight back in through a missing ``-v`` flag.

The output is therefore large. Restricting to callable CDS targets (about
14 Mb, not 276 Mb) is what keeps it tractable.

Depth assumptions
-----------------
Input BAMs are expected to be depth-matched (``depth_match.py``), so that
retained and eliminated targets carry equal expected depth in every worm.
``--min-mq`` must match the value used there, or the depth the caller sees is
not the depth that was matched.

``-B`` (no BAQ) is deliberate: BAQ down-weights bases near indels, and its
effect scales with local divergence, which differs between region classes.
Indel realignment quality is handled by the caller's own model instead.

Example
-------
    python joint_genotype.py --bams bam_matched/*.matched.bam \\
        --genome $DB_DIR/ascaris/v50/AG.v50.fa \\
        --targets-bed targets/callable_cds_targets.bed \\
        --chrom chr01 --outdir vcf --threads 8

    python joint_genotype.py --merge --outdir vcf --threads 8
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path


def require(tool: str) -> str:
    p = shutil.which(tool)
    if p is None:
        logging.error("%s not found on PATH", tool)
        sys.exit(2)
    return p


def sh(cmd: str, desc: str) -> None:
    logging.info("%s", desc)
    logging.debug("CMD: %s", cmd)
    proc = subprocess.run(["bash", "-o", "pipefail", "-c", cmd],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        logging.error("%s FAILED (rc=%d)\n%s", desc, proc.returncode,
                      proc.stderr.strip()[-3000:])
        sys.exit(proc.returncode)
    if proc.stderr.strip():
        for line in proc.stderr.strip().splitlines()[-10:]:
            logging.debug("  | %s", line)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="All-sites joint genotyping over callable CDS targets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--merge", action="store_true",
                   help="Concatenate per-chromosome BCFs and index; ignores "
                        "--chrom and the calling options")

    p.add_argument("--bams", nargs="+", type=Path,
                   help="Depth-matched BAMs (one per worm)")
    p.add_argument("--genome", type=Path, help="Germline reference FASTA")
    p.add_argument("--targets-bed", type=Path,
                   help="callable_cds_targets.bed")
    p.add_argument("--chrom", help="Chromosome to call (array mode)")

    p.add_argument("--min-mq", type=int, default=30,
                   help="MAPQ floor. MUST match depth_match.py --min-mq.")
    p.add_argument("--min-bq", type=int, default=20, help="Base-quality floor")
    p.add_argument("--max-depth", type=int, default=200,
                   help="bcftools mpileup -d, per file per position. A cap, "
                        "not a target: it stops collapsed repeats from "
                        "dominating runtime.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def do_merge(args, bcftools: str) -> int:
    outdir = args.outdir.expanduser().resolve()
    raw = outdir / "raw"
    parts = sorted(raw.glob("*.bcf"),
                   key=lambda p: (len(p.stem), p.stem))  # chr01..chr19 natural
    if not parts:
        logging.error("No per-chromosome BCFs in %s", raw)
        return 1
    logging.info("Concatenating %d part(s): %s", len(parts),
                 ", ".join(p.stem for p in parts))
    merged = outdir / "allsites.bcf"
    lst = outdir / "parts.txt"
    lst.write_text("\n".join(str(p) for p in parts) + "\n")
    sh(f"{bcftools} concat -f {lst} -O b -o {merged} --threads {args.threads}",
       "bcftools concat")
    sh(f"{bcftools} index -f --threads {args.threads} {merged}", "index")

    # Site counts are the sanity check worth having in the log: an all-sites
    # file should be overwhelmingly invariant.
    sh(f"{bcftools} index -n {merged} > {outdir}/n_sites.txt", "count sites")
    n = (outdir / "n_sites.txt").read_text().strip()
    sh(f"{bcftools} view -H -v snps,indels {merged} | wc -l > {outdir}/n_variant.txt",
       "count variant sites")
    v = (outdir / "n_variant.txt").read_text().strip()
    logging.info("allsites.bcf: %s total sites, %s variant (%.2f%%)",
                 n, v, 100 * int(v) / max(int(n), 1))
    if int(v) == int(n):
        logging.error("Every site is variant — the caller was run in "
                      "variants-only mode. piN/piS needs the invariant sites; "
                      "re-run without -v.")
        return 1
    logging.info("Merged: %s", merged)
    return 0


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    bcftools = require("bcftools")

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if args.merge:
        return do_merge(args, bcftools)

    missing = [n for n in ("bams", "genome", "targets_bed", "chrom")
               if getattr(args, n) is None]
    if missing:
        logging.error("Calling mode needs: %s", ", ".join("--" + m.replace("_", "-")
                                                          for m in missing))
        return 2

    raw = outdir / "raw"
    raw.mkdir(exist_ok=True)
    out_bcf = raw / f"{args.chrom}.bcf"
    if out_bcf.exists():
        logging.info("%s exists; nothing to do", out_bcf)
        return 0

    bams = [b.expanduser().resolve() for b in args.bams]
    for b in bams:
        if not b.exists():
            logging.error("Missing BAM: %s", b)
            return 2
        if not (Path(str(b) + ".bai").exists() or b.with_suffix(".bai").exists()):
            logging.error("Missing index for %s", b)
            return 2
    logging.info("%d sample(s), chromosome %s", len(bams), args.chrom)

    bamlist = raw / f"bams_{args.chrom}.txt"
    bamlist.write_text("\n".join(str(b) for b in bams) + "\n")

    # Per-chromosome slice of the target BED. bcftools -R on a whole-genome BED
    # would make every array task scan every interval.
    targets = args.targets_bed.expanduser()
    chrom_bed = raw / f"targets_{args.chrom}.bed"
    n = 0
    with targets.open() as fh, chrom_bed.open("w") as out:
        for line in fh:
            if line.split("\t", 1)[0] == args.chrom:
                out.write(line)
                n += 1
    if n == 0:
        logging.warning("No target intervals on %s — writing nothing", args.chrom)
        return 0
    logging.info("%d target interval(s) on %s", n, args.chrom)

    # -a AD,DP: needed downstream for per-sample depth filtering and allele
    #   balance. FORMAT/DP is what a later setGT step keys on.
    # -B: BAQ off, see module docstring.
    # call -m without -v: ALL sites, variant and invariant. This is the point.
    cmd = (
        f"{bcftools} mpileup -f {args.genome.expanduser()} -b {bamlist} "
        f"-R {chrom_bed} -q {args.min_mq} -Q {args.min_bq} -d {args.max_depth} "
        # Indels are NOT skipped (-I omitted deliberately): the segregating-LoF
        # analysis needs coding indels. They are split off downstream — piN/piS
        # uses SNPs plus invariant sites only, since an indel has no meaningful
        # per-site denominator.
        f"-a FORMAT/AD,FORMAT/DP -B "
        f"--threads {max(1, args.threads // 2)} -O u "
        f"| {bcftools} call -m -f GQ -O b "
        f"--threads {max(1, args.threads // 2)} -o {out_bcf}"
    )
    sh(cmd, f"mpileup + call ({args.chrom}, all sites)")
    sh(f"{bcftools} index -f --threads {args.threads} {out_bcf}", "index")

    sh(f"{bcftools} index -n {out_bcf}", "count")
    logging.info("Done: %s", out_bcf)
    return 0


if __name__ == "__main__":
    sys.exit(main())
