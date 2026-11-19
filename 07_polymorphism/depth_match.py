#!/usr/bin/env python3
"""
depth_match.py — equalise sequencing depth across region classes and across
worms before genotyping. This is the step that removes the analysis's main
source of bias, and it has to happen before a single genotype is called.

The problem
-----------
Eliminated regions are covered at 64-112% of retained depth across the ten
selected worms (median 78%, on the callable CDS targets), because only the
gonad fraction of each prep carries them (see the germline fraction audit).
Unequal depth is not merely unequal power: heterozygote calling is
depth-dependent, so a site at 10x has a higher chance of being
miscalled homozygous than the same site at 25x. Left alone, that depresses
apparent diversity in eliminated regions specifically — which is the exact
comparison the paper turns on, biased in the exact direction that would look
like a result.

Worse, the deficit differs per worm (germline fraction 0.34 to 0.72), so the
bias is not even a constant offset that might cancel in a ratio.

The fix
-------
Downsample reads so that **every worm has the same expected depth in both
retained and eliminated targets**. Retained reads are thinned harder than
eliminated ones within a worm, and worms with more gonad are thinned harder
than worms with less. After this, genotype-calling power is uniform by
construction and no downstream filter has to model the asymmetry.

This is the same reasoning as downsampling before ICE rather than trusting CPM
to fix depth differences between Hi-C libraries: normalisation applied after a
depth-dependent transformation cannot undo it.

Cost: throwing away reads. The target depth is set by the worst-covered class
of the worst-covered worm, so choose ``--target-depth`` explicitly once you
have seen the audit of the full-depth BAMs, rather than letting one bad library
drag the whole cohort down. Worms that cannot reach the target are reported and
skipped rather than silently left un-thinned.

Changes 2026-10-12
------------------
* The achieved depth is now always measured afresh on the matched BAM. Before,
  mean_depth() reused any mosdepth summary in ``tmp/`` with the right name, so a
  summary left by an earlier build survived a rebuild at another target. That
  is how SRR31692077's table came to say 8.00x and 8.02x for a BAM that
  mosdepth measures at 13.00x in both classes.
* A cached summary of an INPUT BAM is reused only if it is newer than both the
  BAM and the BED, and the cache name carries the MAPQ floor.
* The shared target BEDs are rewritten only when their content changes, and
  then atomically. Every array task writes them into the same --outdir, and an
  in-place rewrite could truncate the file under a sibling task reading it.

Example
-------
    python depth_match.py --bam bam_full/SRR31692074.bam \\
        --targets-bed targets/callable_cds_targets.bed \\
        --eliminated-bed $DB_DIR/ascaris/v50/AG.v50.eliminated_strict.bed \\
        --target-depth 13 --outdir bam_matched --threads 8
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path


def require(tool: str) -> str:
    p = shutil.which(tool)
    if p is None:
        logging.error("%s not found on PATH (activate the popgen env and set "
                      "PATH to its bin)", tool)
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


# ---------------------------------------------------------------------------
# Interval helpers (kept local so this script has no project imports)
# ---------------------------------------------------------------------------

def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    raw: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    out = {}
    for c, ivs in raw.items():
        ivs.sort()
        merged = [ivs[0]]
        for s, e in ivs[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        out[c] = merged
    return out


def overlap_bp(iv: tuple[int, int], ivs: list[tuple[int, int]]) -> int:
    s, e = iv
    total = 0
    for bs, be in ivs:
        if be <= s:
            continue
        if bs >= e:
            break
        total += min(e, be) - max(s, bs)
    return total


def split_targets(targets: dict[str, list[tuple[int, int]]],
                  elim: dict[str, list[tuple[int, int]]],
                  outdir: Path, min_frac: float
                  ) -> tuple[Path, Path, int, int]:
    """
    Split the target intervals into retained-only and eliminated-only BEDs.

    Intervals straddling a boundary are DROPPED, not assigned. A read pair
    overlapping both files would otherwise be extracted twice and thinned at
    two different rates, producing a depth spike exactly at the CBR — the one
    place in the genome where the paper makes a quantitative claim.
    """
    ret_bed = outdir / "targets_retained.bed"
    eli_bed = outdir / "targets_eliminated.bed"
    ret_lines, eli_lines = [], []
    n_drop = 0
    for chrom in sorted(targets):
        e_ivs = elim.get(chrom, [])
        for s, e in targets[chrom]:
            frac = overlap_bp((s, e), e_ivs) / (e - s)
            if frac >= min_frac:
                eli_lines.append(f"{chrom}\t{s}\t{e}\n")
            elif frac <= 1 - min_frac:
                ret_lines.append(f"{chrom}\t{s}\t{e}\n")
            else:
                n_drop += 1
    for path, lines in ((ret_bed, ret_lines), (eli_bed, eli_lines)):
        write_if_changed(path, "".join(lines))
    logging.info("Targets split: %d retained, %d eliminated, %d boundary "
                 "intervals dropped", len(ret_lines), len(eli_lines), n_drop)
    return ret_bed, eli_bed, len(ret_lines), len(eli_lines)


def write_if_changed(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` only if it differs, via a temporary file and
    an atomic rename.

    All array tasks share one --outdir, so they all write these BEDs. Rewriting
    in place truncates the file under any sibling that is reading it; a rename
    swaps it in whole, and leaving an unchanged file alone also keeps its mtime,
    which mean_depth() uses to decide whether a cached summary is still valid.
    """
    if path.exists() and path.read_text() == text:
        return
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# Depth
# ---------------------------------------------------------------------------

def mean_depth(bam: Path, bed: Path, tmp: Path, min_mq: int,
               threads: int, reuse: bool = True) -> float:
    """
    Mean depth over the intervals in ``bed``, from mosdepth's summary.

    The 'total' row of mosdepth's region summary is the coverage-weighted mean
    across all listed regions, which is what the downsampling fraction needs —
    a median over intervals would be pulled around by short ones.

    A cached summary is reused only when ``reuse`` is set and the summary is
    newer than both ``bam`` and ``bed``; the cache name carries the MAPQ floor.
    Before 2026-10-12 any summary with the right name was reused, so the
    achieved depth of a rebuilt matched BAM could be read from the previous
    build (SRR31692077: 8.00x / 8.02x reported, 13.00x / 13.00x on disk).
    """
    prefix = tmp / f"md_{bam.stem}_{bed.stem}_q{min_mq}"
    summary = Path(f"{prefix}.mosdepth.summary.txt")
    valid = (reuse and summary.exists()
             and summary.stat().st_mtime > max(bam.stat().st_mtime,
                                               bed.stat().st_mtime))
    if valid:
        logging.info("Reusing %s (newer than %s and %s)", summary.name,
                     bam.name, bed.name)
    else:
        summary.unlink(missing_ok=True)
        sh(f"mosdepth -n -x -t {threads} -Q {min_mq} --by {bed} {prefix} {bam}",
           f"mosdepth {bam.name} over {bed.name}")
    for line in summary.read_text().splitlines():
        f = line.split("\t")
        if f[0] == "total_region":
            return float(f[3])
    logging.error("No total_region row in %s", summary)
    sys.exit(1)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Downsample a BAM to a common depth in retained and "
                    "eliminated targets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--bam", required=True, type=Path,
                   help="Full-depth, coordinate-sorted, indexed BAM")
    p.add_argument("--targets-bed", required=True, type=Path,
                   help="callable_cds_targets.bed from make_callable_targets.py")
    p.add_argument("--eliminated-bed", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--target-depth", type=float, required=True,
                   help="Depth to match BOTH classes to, in every worm. Set it "
                        "from the audit of the full-depth BAMs; do not guess.")
    p.add_argument("--min-mq", type=int, default=30,
                   help="MAPQ floor, used both for measuring depth here and "
                        "for calling later — they must agree or the achieved "
                        "depth will not be the depth the caller sees")
    p.add_argument("--min-frac", type=float, default=0.9,
                   help="Overlap needed to call a target interval eliminated "
                        "or retained; straddling intervals are dropped")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--seed", type=int, default=42,
                   help="Subsampling seed. Fixed so the cohort is reproducible; "
                        "vary it only for a deliberate sensitivity run.")
    p.add_argument("--tolerance", type=float, default=0.90,
                   help="Skip the BAM if a class sits below this multiple of "
                        "--target-depth (it cannot be matched by thinning)")
    p.add_argument("--force", action="store_true")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    samtools = require("samtools")
    require("mosdepth")

    bam = args.bam.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    tmp = outdir / "tmp"
    tmp.mkdir(exist_ok=True)

    stem = bam.stem
    out_bam = outdir / f"{stem}.matched.bam"
    if out_bam.exists() and not args.force:
        logging.info("%s exists; nothing to do (--force to redo)", out_bam)
        return 0

    ret_bed, eli_bed, n_ret, n_eli = split_targets(
        read_bed(args.targets_bed.expanduser()),
        read_bed(args.eliminated_bed.expanduser()),
        outdir, args.min_frac)
    if n_eli == 0:
        logging.error("No eliminated target intervals — check the BEDs share a "
                      "coordinate system (both must be germline v50)")
        return 2

    d_ret = mean_depth(bam, ret_bed, tmp, args.min_mq, args.threads)
    d_eli = mean_depth(bam, eli_bed, tmp, args.min_mq, args.threads)
    logging.info("%s: retained %.2fx, eliminated %.2fx (ratio %.3f)",
                 stem, d_ret, d_eli, d_eli / d_ret if d_ret else float("nan"))

    floor = args.tolerance * args.target_depth
    if d_eli < floor or d_ret < floor:
        logging.warning("=" * 70)
        logging.warning("%s cannot reach target depth %.1fx (retained %.2f, "
                        "eliminated %.2f). SKIPPING — thinning cannot add "
                        "reads, and leaving it un-thinned would reintroduce "
                        "exactly the asymmetry this step removes.",
                        stem, args.target_depth, d_ret, d_eli)
        logging.warning("Either lower --target-depth for the whole cohort or "
                        "drop this worm. Do not mix.")
        logging.warning("=" * 70)
        (outdir / f"{stem}.SKIPPED").write_text(
            f"retained={d_ret:.3f}\neliminated={d_eli:.3f}\n"
            f"target={args.target_depth}\n")
        return 3

    f_ret = min(1.0, args.target_depth / d_ret)
    f_eli = min(1.0, args.target_depth / d_eli)
    logging.info("Subsampling fractions: retained %.4f, eliminated %.4f",
                 f_ret, f_eli)

    part_r = tmp / f"{stem}.ret.bam"
    part_e = tmp / f"{stem}.eli.bam"
    # -M with -L restricts to reads whose alignment actually overlaps a target
    # interval, rather than the looser default; keeps the two parts disjoint.
    sh(f"{samtools} view -@ {args.threads} -b -M -L {ret_bed} -q {args.min_mq} "
       f"--subsample {f_ret:.6f} --subsample-seed {args.seed} "
       f"-o {part_r} {bam}", "extract + thin retained targets")
    sh(f"{samtools} view -@ {args.threads} -b -M -L {eli_bed} -q {args.min_mq} "
       f"--subsample {f_eli:.6f} --subsample-seed {args.seed} "
       f"-o {part_e} {bam}", "extract + thin eliminated targets")

    sh(f"{samtools} merge -@ {args.threads} -f -o - {part_r} {part_e} "
       f"| {samtools} sort -@ {args.threads} -m 1G -T {tmp}/{stem}.srt "
       f"-o {out_bam}", "merge + sort")
    sh(f"{samtools} index -@ {args.threads} {out_bam}", "index")

    # Verify rather than assume: the achieved depths are what the caller sees.
    # Always measured afresh (reuse=False): this BAM was just written, so any
    # summary of the same name describes an earlier build.
    a_ret = mean_depth(out_bam, ret_bed, tmp, args.min_mq, args.threads,
                       reuse=False)
    a_eli = mean_depth(out_bam, eli_bed, tmp, args.min_mq, args.threads,
                       reuse=False)
    logging.info("Achieved: retained %.2fx, eliminated %.2fx (ratio %.3f)",
                 a_ret, a_eli, a_eli / a_ret if a_ret else float("nan"))
    if abs(a_ret - a_eli) > 0.15 * args.target_depth:
        logging.warning("Achieved depths differ by more than 15%% of target — "
                        "check for boundary intervals or duplicate reads")

    (outdir / f"{stem}.depth_match.tsv").write_text(
        "sample\tdepth_retained_in\tdepth_eliminated_in\tfrac_retained\t"
        "frac_eliminated\tdepth_retained_out\tdepth_eliminated_out\ttarget\n"
        f"{stem}\t{d_ret:.4f}\t{d_eli:.4f}\t{f_ret:.6f}\t{f_eli:.6f}\t"
        f"{a_ret:.4f}\t{a_eli:.4f}\t{args.target_depth}\n")

    for p in (part_r, part_e):
        p.unlink(missing_ok=True)
    logging.info("Done: %s", out_bam)
    return 0


if __name__ == "__main__":
    sys.exit(main())
