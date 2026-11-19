#!/usr/bin/env python3
"""
mappability_mask.py — build a low-mappability exclusion mask with genmap, and
report how much of each region class it removes.

Why
---
Eliminated sequence is repeat-rich. Short reads that cannot be placed uniquely
there produce spurious coverage and, worse, spurious heterozygotes: paralogous
copies collapse onto one locus and their fixed differences are called as
segregating sites. That inflates diversity in exactly the gene class the paper
is about, in exactly the direction that would manufacture the result.

The audit already showed that eliminated coverage in these worms is real rather
than mismapped (somatic controls read ~0.01-0.05 against a spermatid control of
1.12). This mask is the second line of defence, applied per site rather than
per region, before any diversity statistic is computed.

Method
------
genmap computes, for every position, the number of times its k-mer occurs in
the genome allowing E mismatches. Mappability is 1/frequency, so 1.0 means the
k-mer is unique. Positions below ``--min-mappability`` are masked.

k should match the read length actually used (150 for these libraries). E=2
allows for real divergence between *A. lumbricoides* reads and the *A. suum*
reference, which matters here: with E=0 the mask would be blind to paralogs
that are distinguishable only under an unrealistically strict match.

The masked fraction PER REGION CLASS is the number to read. If the mask removes
70% of eliminated sequence and 5% of retained, every downstream comparison is
made on very unequal footing and the per-callable-site normalisation has to be
airtight.

Example
-------
    python mappability_mask.py --genome $DB_DIR/ascaris/v50/AG.v50.fa \\
        --eliminated-bed $DB_DIR/ascaris/v50/AG.v50.eliminated_strict.bed \\
        --kmer 150 --errors 2 --outdir $POPGEN_DIR/mask --threads 16
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path


def require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        logging.error("%s not found on PATH. Install with:", tool)
        logging.error("  conda install -n popgen "
                      "-c conda-forge -c bioconda -y genmap")
        sys.exit(2)
    return path


def run(cmd: list[str], desc: str) -> None:
    logging.info("%s", desc)
    logging.debug("CMD: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logging.error("%s FAILED (rc=%d)", desc, proc.returncode)
        logging.error("stderr tail:\n%s", proc.stderr.strip()[-4000:])
        sys.exit(proc.returncode)


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Read a BED into merged {chrom: [(start, end)]}."""
    raw: dict[str, list[tuple[int, int]]] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    out: dict[str, list[tuple[int, int]]] = {}
    for c, ivs in raw.items():
        ivs.sort()
        merged: list[tuple[int, int]] = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        out[c] = merged
    return out


def total_bp(ivs: dict[str, list[tuple[int, int]]]) -> int:
    return sum(e - s for v in ivs.values() for s, e in v)


def intersect_bp(a: dict[str, list[tuple[int, int]]],
                 b: dict[str, list[tuple[int, int]]]) -> int:
    """Total bp of intersection between two merged interval sets."""
    total = 0
    for chrom, av in a.items():
        bv = b.get(chrom)
        if not bv:
            continue
        i = j = 0
        while i < len(av) and j < len(bv):
            lo = max(av[i][0], bv[j][0])
            hi = min(av[i][1], bv[j][1])
            if hi > lo:
                total += hi - lo
            if av[i][1] < bv[j][1]:
                i += 1
            else:
                j += 1
    return total


def bedgraph_to_mask(bedgraph: Path, dest: Path, threshold: float) -> int:
    """
    Write intervals with mappability < ``threshold`` to ``dest``, merging
    adjacent runs. Returns the number of masked bp.

    Streamed line by line: the bedgraph for a 276 Mb genome is large and there
    is no reason to hold it in memory.
    """
    masked = 0
    n_iv = 0
    with bedgraph.open() as fh, dest.open("w") as out:
        cur_chrom, cur_start, cur_end = None, None, None
        for line in fh:
            if not line.strip():
                continue
            chrom, start, end, value = line.split("\t")[:4]
            if float(value) >= threshold:
                continue
            start, end = int(start), int(end)
            if cur_chrom == chrom and start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                if cur_chrom is not None:
                    out.write(f"{cur_chrom}\t{cur_start}\t{cur_end}\n")
                    masked += cur_end - cur_start
                    n_iv += 1
                cur_chrom, cur_start, cur_end = chrom, start, end
        if cur_chrom is not None:
            out.write(f"{cur_chrom}\t{cur_start}\t{cur_end}\n")
            masked += cur_end - cur_start
            n_iv += 1
    logging.info("Mask: %d intervals, %.2f Mb", n_iv, masked / 1e6)
    return masked


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build a genmap low-mappability mask and report its effect "
                    "per region class.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--genome", required=True, type=Path,
                   help="Germline reference FASTA (AG.v50.fa)")
    p.add_argument("--eliminated-bed", required=True, type=Path,
                   help="Eliminated regions, germline coordinates")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--kmer", type=int, default=150,
                   help="k-mer length; match the read length (150 here)")
    p.add_argument("--errors", type=int, default=2,
                   help="Mismatches allowed. 2 accounts for real A. lumbricoides "
                        "vs A. suum divergence; 0 would under-mask paralogs.")
    p.add_argument("--min-mappability", type=float, default=1.0,
                   help="Positions below this are masked (1.0 = keep only "
                        "uniquely mappable positions)")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--index-dir", type=Path,
                   help="genmap index directory (default <outdir>/genmap_index). "
                        "Reused if it already exists — indexing is the slow step.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    genmap = require("genmap")
    genome = args.genome.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    index_dir = (args.index_dir or (outdir / "genmap_index")).expanduser().resolve()

    # genmap refuses to write into an existing index directory, so only build
    # when it is absent. Indexing 276 Mb takes ~20 min and ~10 GB.
    if index_dir.exists() and any(index_dir.iterdir()):
        logging.info("Reusing genmap index: %s", index_dir)
    else:
        if index_dir.exists():
            index_dir.rmdir()
        run([genmap, "index", "-F", str(genome), "-I", str(index_dir)],
            f"genmap index ({genome.name})")

    prefix = outdir / f"mappability_k{args.kmer}_e{args.errors}"
    bedgraph = Path(f"{prefix}.bedgraph")
    if bedgraph.exists():
        logging.info("Reusing existing bedgraph: %s", bedgraph)
    else:
        run([genmap, "map", "-K", str(args.kmer), "-E", str(args.errors),
             "-I", str(index_dir), "-O", str(prefix),
             "-T", str(args.threads), "-bg"],
            f"genmap map (K={args.kmer}, E={args.errors})")

    mask_bed = outdir / f"lowmap_k{args.kmer}_e{args.errors}.bed"
    bedgraph_to_mask(bedgraph, mask_bed, args.min_mappability)

    # --- effect per region class ------------------------------------------
    elim = read_bed(args.eliminated_bed.expanduser())
    mask = read_bed(mask_bed)

    elim_bp = total_bp(elim)
    mask_bp = total_bp(mask)
    elim_masked = intersect_bp(elim, mask)

    fai = Path(str(genome) + ".fai")
    if not fai.exists():
        logging.warning("No .fai beside the genome; cannot compute the retained "
                        "denominator. Run `samtools faidx %s`.", genome)
        return 0
    genome_bp = sum(int(l.split("\t")[1]) for l in fai.read_text().splitlines() if l.strip())
    retained_bp = genome_bp - elim_bp
    retained_masked = mask_bp - elim_masked

    logging.info("=" * 70)
    logging.info("Mask effect by region class")
    logging.info("  retained   : %8.2f Mb, %6.2f Mb masked (%.1f%%)",
                 retained_bp / 1e6, retained_masked / 1e6,
                 100 * retained_masked / retained_bp)
    logging.info("  eliminated : %8.2f Mb, %6.2f Mb masked (%.1f%%)",
                 elim_bp / 1e6, elim_masked / 1e6,
                 100 * elim_masked / elim_bp)
    ratio = ((elim_masked / elim_bp) / (retained_masked / retained_bp)
             if retained_masked else float("nan"))
    logging.info("  eliminated is %.2fx more masked than retained", ratio)
    logging.info("")
    logging.info("Read this before trusting any downstream diversity estimate:")
    logging.info("  * a large asymmetry does NOT invalidate the comparison, but")
    logging.info("    every statistic must then be per CALLABLE site, never per bp;")
    logging.info("  * if eliminated survivorship is very low, check how many")
    logging.info("    eliminated GENES retain enough callable CDS to be tested.")
    logging.info("Mask written: %s", mask_bed)
    logging.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
