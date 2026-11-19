#!/usr/bin/env python3
"""
germline_fraction_audit.py — the go/no-go test for a polymorphism analysis of
eliminated genes in whole-worm Ascaris WGS.

The question
------------
Eliminated sequence is present only in germline nuclei. Landeryou et al. (2025,
PRJNA1195031) extracted DNA from *whole* adult worms, so each library is a
mixture of somatic cells (eliminated regions absent) and gonad (eliminated
regions present at 2 copies per nucleus). The usable depth over eliminated
regions is therefore some unknown fraction of the retained-region depth.

If that fraction is high (say >0.3), eliminated genes are genotypable and the
whole pi_N/pi_S and MK-test programme is on. If it is ~0.02, it is off, and it
is far better to learn that from three pilot BAMs than from a completed
pipeline whose eliminated-gene diversity estimates are made of noise.

What it measures, per sample
----------------------------
  * median depth in retained windows (autosomes only)
  * median depth in eliminated windows
  * germline_fraction = eliminated / retained median depth
  * callability: fraction of eliminated windows at depth >= 5, >= 10, >= 20
  * X:A depth ratio, which separates females (~1.0) from males (~0.5)

Sex matters for interpretation: a gravid female's uterus is packed with
pre-elimination embryos, so females should show a *higher* germline fraction
than males, and their eliminated-region reads carry paternal alleles from one
or more sires. Do not pool the two sexes before checking this panel.

Depth is measured with mosdepth if available (fast, MAPQ-filtered) and falls
back to ``samtools bedcov`` otherwise.

Example
-------
    python germline_fraction_audit.py \
        --bams bam/*.bam \
        --eliminated-bed $DB_DIR/ascaris/v50/AG.v50.eliminated_strict.bed \
        --fai $DB_DIR/ascaris/v50/AG.v50.fa.fai \
        --window 10000 --min-mapq 20 --outdir audit_20260818
"""

from __future__ import annotations

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- figure style (paper convention) ---------------------------------------
# retained = blue/left, eliminated = red/right, everywhere.
try:                                    # lib/paper_colors.py if importable
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
    from paper_colors import RETAINED, ELIMINATED       # type: ignore
except Exception:                       # standalone fallback
    RETAINED, ELIMINATED = "#2C6FBB", "#C0392B"

plt.rcParams.update({
    "font.family": "Arial",
    "svg.fonttype": "none",             # keep text editable in Illustrator
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 15,
    "figure.dpi": 150,
})


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------

def read_fai(path: Path) -> dict[str, int]:
    """Return {sequence_name: length} from a samtools .fai index."""
    lengths: dict[str, int] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.split("\t")
            lengths[f[0]] = int(f[1])
    logging.info("Reference: %d sequences, %.1f Mb",
                 len(lengths), sum(lengths.values()) / 1e6)
    return lengths


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Read a BED into {chrom: [(start, end), ...]}, merged and sorted."""
    raw: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))

    merged: dict[str, list[tuple[int, int]]] = {}
    total = 0
    for chrom, ivs in raw.items():
        ivs.sort()
        out: list[tuple[int, int]] = []
        for s, e in ivs:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
        total += sum(e - s for s, e in out)
    logging.info("%s: %d intervals on %d sequences, %.2f Mb",
                 path.name, sum(len(v) for v in merged.values()),
                 len(merged), total / 1e6)
    return merged


def make_windows(lengths: dict[str, int], window: int) -> pd.DataFrame:
    """Tile every sequence into fixed windows; drop the ragged terminal window."""
    rows = []
    for chrom, ln in lengths.items():
        for start in range(0, ln - window + 1, window):
            rows.append((chrom, start, start + window))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    logging.info("Windows: %d of %d bp", len(df), window)
    return df


def overlap_fraction(df: pd.DataFrame,
                     intervals: dict[str, list[tuple[int, int]]]) -> np.ndarray:
    """
    Fraction of each window covered by ``intervals``.

    Windows are already sorted within a chromosome, so a single merged sweep
    per chromosome is O(n + m) rather than the O(n*m) a naive per-window search
    would cost on 30,000 windows x 100+ tracts.
    """
    frac = np.zeros(len(df), dtype=float)
    for chrom, sub in df.groupby("chrom", sort=False):
        ivs = intervals.get(chrom)
        if not ivs:
            continue
        idx = sub.index.to_numpy()
        starts = sub["start"].to_numpy()
        ends = sub["end"].to_numpy()
        j = 0
        for k in range(len(starts)):
            ws, we = starts[k], ends[k]
            while j < len(ivs) and ivs[j][1] <= ws:
                j += 1
            covered = 0
            jj = j
            while jj < len(ivs) and ivs[jj][0] < we:
                covered += min(we, ivs[jj][1]) - max(ws, ivs[jj][0])
                jj += 1
            frac[idx[k]] = covered / (we - ws)
    return frac


# ---------------------------------------------------------------------------
# Depth
# ---------------------------------------------------------------------------

def depth_mosdepth(bam: Path, windows_bed: Path, outdir: Path,
                   min_mapq: int, threads: int) -> pd.DataFrame:
    """Per-window mean depth via mosdepth --by."""
    prefix = outdir / f"mosdepth_{bam.stem}"
    regions = Path(f"{prefix}.regions.bed.gz")
    if not regions.exists():
        cmd = ["mosdepth", "-n", "-x", "-t", str(threads),
               "-Q", str(min_mapq), "--by", str(windows_bed),
               str(prefix), str(bam)]
        logging.info("mosdepth %s", bam.name)
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    df = pd.read_csv(regions, sep="\t", header=None,
                     names=["chrom", "start", "end", "depth"])
    return df


def depth_bedcov(bam: Path, windows_bed: Path, outdir: Path,
                 min_mapq: int) -> pd.DataFrame:
    """
    Per-window mean depth via ``samtools bedcov`` (fallback when mosdepth is
    absent). bedcov reports the *sum* of per-base depth, so divide by width.
    """
    out = outdir / f"bedcov_{bam.stem}.tsv"
    if not out.exists():
        logging.info("samtools bedcov %s (slower than mosdepth)", bam.name)
        with out.open("w") as fh:
            subprocess.run(["samtools", "bedcov", "-Q", str(min_mapq),
                            str(windows_bed), str(bam)],
                           check=True, stdout=fh)
    df = pd.read_csv(out, sep="\t", header=None,
                     names=["chrom", "start", "end", "cov_sum"])
    df["depth"] = df["cov_sum"] / (df["end"] - df["start"])
    return df[["chrom", "start", "end", "depth"]]


# ---------------------------------------------------------------------------
# Per-sample statistics
# ---------------------------------------------------------------------------

def sample_stats(df: pd.DataFrame, x_prefix: str,
                 hi_pct: float) -> dict:
    """
    Summarise one sample's window depths.

    ``hi_pct`` trims the top percentile of retained-window depth before taking
    medians. Collapsed repeats pile reads into a handful of windows and are
    much denser in eliminated sequence; without the trim the germline fraction
    is biased upward by exactly the regions least suitable for genotyping.
    """
    auto = df[~df["chrom"].str.startswith(x_prefix)]
    ret = auto[auto["cls"] == "retained"]["depth"]
    elim = auto[auto["cls"] == "eliminated"]["depth"]

    if ret.empty:
        raise ValueError("no retained autosomal windows survived filtering")

    hi = np.percentile(ret, hi_pct)
    ret_t = ret[ret <= hi]
    elim_t = elim[elim <= hi] if not elim.empty else elim

    d_ret = float(np.median(ret_t))
    d_elim = float(np.median(elim_t)) if not elim_t.empty else np.nan

    x_ret = df[df["chrom"].str.startswith(x_prefix) & (df["cls"] == "retained")]["depth"]
    x_ret = x_ret[x_ret <= hi]

    return {
        "n_retained_windows": int(ret_t.size),
        "n_eliminated_windows": int(elim_t.size),
        "depth_retained_median": d_ret,
        "depth_eliminated_median": d_elim,
        "germline_fraction": d_elim / d_ret if d_ret > 0 else np.nan,
        "elim_frac_ge5": float((elim >= 5).mean()) if elim.size else np.nan,
        "elim_frac_ge10": float((elim >= 10).mean()) if elim.size else np.nan,
        "elim_frac_ge20": float((elim >= 20).mean()) if elim.size else np.nan,
        "x_to_autosome": float(np.median(x_ret) / d_ret) if x_ret.size and d_ret > 0 else np.nan,
        "depth_p99_retained": float(hi),
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def save(fig, outdir: Path, stem: str) -> None:
    """Export PNG/SVG/PDF. SVG is the one Illustrator should open — Arial is
    not installed on Isaac compute nodes, so PNG/PDF render in a fallback."""
    for ext in ("png", "svg", "pdf"):
        fig.savefig(outdir / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    logging.info("Figure: %s.{png,svg,pdf}", stem)


def plot_depth_distributions(per_window: pd.DataFrame, outdir: Path,
                             x_prefix: str) -> None:
    """Pooled autosomal window-depth distributions, retained left / eliminated right."""
    auto = per_window[~per_window["chrom"].str.startswith(x_prefix)]
    ret = auto[auto["cls"] == "retained"]["depth"].to_numpy()
    elim = auto[auto["cls"] == "eliminated"]["depth"].to_numpy()

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    parts = ax.violinplot([ret, elim], positions=[0, 1], widths=0.75,
                          showmedians=True, showextrema=False)
    for body, colour in zip(parts["bodies"], (RETAINED, ELIMINATED)):
        body.set_facecolor(colour)
        body.set_alpha(0.75)
        body.set_edgecolor("black")
    parts["cmedians"].set_color("black")

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Retained\n(n = {ret.size:,})",
                        f"Eliminated\n(n = {elim.size:,})"])
    ax.set_ylabel("Window depth (x)")
    ax.set_ylim(bottom=0)
    save(fig, outdir, "window_depth_by_region")


def plot_germline_fraction(stats: pd.DataFrame, outdir: Path) -> None:
    """Per-sample germline fraction, ordered, with the usability threshold drawn."""
    s = stats.sort_values("germline_fraction")
    fig, ax = plt.subplots(figsize=(max(6.0, 0.28 * len(s)), 5.0))
    ax.bar(range(len(s)), s["germline_fraction"], color=ELIMINATED,
           edgecolor="black", linewidth=0.4)
    ax.axhline(0.30, ls="--", lw=1.5, color="grey")
    ax.text(0.01, 0.31, "0.30", transform=ax.get_yaxis_transform(),
            va="bottom", fontsize=13, color="grey")
    ax.set_xticks(range(len(s)))
    ax.set_xticklabels(s["sample"], rotation=90, fontsize=9)
    ax.set_ylabel("Eliminated / retained depth")
    ax.set_xlabel(f"Sample (n = {len(s)})")
    save(fig, outdir, "germline_fraction_by_sample")


def plot_sex(stats: pd.DataFrame, outdir: Path) -> None:
    """X:autosome depth against germline fraction — sex, and its consequence."""
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(stats["x_to_autosome"], stats["germline_fraction"],
               s=60, color=ELIMINATED, edgecolor="black", linewidth=0.5,
               rasterized=True)
    ax.axvline(0.75, ls="--", lw=1.5, color="grey")
    ax.set_xlabel("X : autosome depth")
    ax.set_ylabel("Eliminated / retained depth")
    ax.text(0.98, 0.02, "males <- | -> females", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=13, color="grey")
    save(fig, outdir, "sex_vs_germline_fraction")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Audit germline (eliminated-region) coverage in whole-worm WGS.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bams", nargs="+", required=True, type=Path,
                   help="Sorted, indexed BAMs mapped to the germline assembly")
    p.add_argument("--eliminated-bed", required=True, type=Path,
                   help="Eliminated regions in germline coordinates "
                        "(AG.v50.eliminated_strict.bed)")
    p.add_argument("--fai", required=True, type=Path,
                   help="samtools .fai for the germline FASTA")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--window", type=int, default=10000,
                   help="Window size in bp")
    p.add_argument("--min-frac", type=float, default=0.9,
                   help="A window is called eliminated at >= this overlap and "
                        "retained at <= 1 - this; the rest are boundary windows "
                        "and are dropped")
    p.add_argument("--min-mapq", type=int, default=20,
                   help="Minimum MAPQ. Eliminated sequence is repeat-rich, so "
                        "counting MAPQ-0 reads would manufacture coverage that "
                        "cannot be genotyped.")
    p.add_argument("--mask-bed", type=Path,
                   help="Optional low-mappability / blacklist BED to exclude")
    p.add_argument("--x-prefix", default="chrX",
                   help="Prefix identifying X sequences (chrX1..chrX5 in v50)")
    p.add_argument("--hi-pct", type=float, default=99.0,
                   help="Trim retained-window depth above this percentile")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    lengths = read_fai(args.fai.expanduser())
    elim = read_bed(args.eliminated_bed.expanduser())

    windows = make_windows(lengths, args.window)
    frac = overlap_fraction(windows, elim)
    windows["elim_frac"] = frac
    windows["cls"] = np.where(frac >= args.min_frac, "eliminated",
                     np.where(frac <= 1 - args.min_frac, "retained", "boundary"))

    if args.mask_bed:
        mask = read_bed(args.mask_bed.expanduser())
        mfrac = overlap_fraction(windows, mask)
        n_before = len(windows)
        windows = windows[mfrac < 0.5].reset_index(drop=True)
        logging.info("Mask removed %d of %d windows", n_before - len(windows), n_before)

    counts = windows["cls"].value_counts()
    logging.info("Window classes: %s", counts.to_dict())
    if counts.get("eliminated", 0) < 50:
        logging.warning("Only %d eliminated windows at %d bp — consider a "
                        "smaller --window", counts.get("eliminated", 0), args.window)

    keep = windows[windows["cls"] != "boundary"].reset_index(drop=True)
    windows_bed = outdir / f"windows_{args.window}.bed"
    keep[["chrom", "start", "end"]].to_csv(windows_bed, sep="\t",
                                           header=False, index=False)

    use_mosdepth = shutil.which("mosdepth") is not None
    if not use_mosdepth:
        logging.warning("mosdepth not found; falling back to samtools bedcov")

    rows, pooled = [], []
    for bam in args.bams:
        bam = bam.expanduser().resolve()
        if not bam.exists():
            logging.warning("Missing BAM, skipped: %s", bam)
            continue
        depth = (depth_mosdepth(bam, windows_bed, outdir, args.min_mapq, args.threads)
                 if use_mosdepth
                 else depth_bedcov(bam, windows_bed, outdir, args.min_mapq))

        merged = keep.merge(depth, on=["chrom", "start", "end"], how="inner")
        if len(merged) != len(keep):
            logging.warning("%s: %d of %d windows matched the depth table",
                            bam.name, len(merged), len(keep))

        # Drop windows that are zero in this sample AND sit in an assembly gap:
        # a genuine zero in eliminated sequence is a result, not an artefact,
        # so only all-N windows are removed, and that is handled by requiring
        # the retained median to be positive rather than by dropping here.
        stat = sample_stats(merged, args.x_prefix, args.hi_pct)
        stat["sample"] = bam.stem
        rows.append(stat)
        merged["sample"] = bam.stem
        pooled.append(merged)

        logging.info("%-28s ret=%.1fx elim=%.1fx  germline_fraction=%.3f  X:A=%.2f",
                     bam.stem, stat["depth_retained_median"],
                     stat["depth_eliminated_median"], stat["germline_fraction"],
                     stat["x_to_autosome"])

    if not rows:
        logging.error("No BAMs processed")
        return 1

    stats = pd.DataFrame(rows)[[
        "sample", "depth_retained_median", "depth_eliminated_median",
        "germline_fraction", "elim_frac_ge5", "elim_frac_ge10", "elim_frac_ge20",
        "x_to_autosome", "n_retained_windows", "n_eliminated_windows",
        "depth_p99_retained",
    ]]
    stats.to_csv(outdir / "germline_fraction_summary.tsv", sep="\t", index=False)

    per_window = pd.concat(pooled, ignore_index=True)
    per_window.to_csv(outdir / "window_depth.tsv.gz", sep="\t", index=False)

    plot_depth_distributions(per_window, outdir, args.x_prefix)
    plot_germline_fraction(stats, outdir)
    if stats["x_to_autosome"].notna().any():
        plot_sex(stats, outdir)

    # --- verdict -----------------------------------------------------------
    # Report how many samples PASS, not the cohort median. Gonad carryover is
    # all-or-nothing per worm, so in a cohort where most preps are clean the
    # median is ~0 while a usable minority sits well above threshold. A median
    # here would say "abandon" when the correct action is "select".
    gf = stats["germline_fraction"]
    n_good = int((gf >= 0.30).sum())
    n_marginal = int(((gf >= 0.10) & (gf < 0.30)).sum())
    logging.info("=" * 70)
    logging.info("Germline fraction across %d sample(s): median %.3f, max %.3f",
                 len(stats), float(gf.median()), float(gf.max()))
    logging.info("  >= 0.30 (usable)   : %d", n_good)
    logging.info("  0.10-0.30 (marginal): %d", n_marginal)
    logging.info("  < 0.10 (somatic)   : %d", int((gf < 0.10).sum()))
    if n_good >= 6:
        logging.info("VERDICT: %d sample(s) carry enough germline DNA. Proceed, "
                     "restricted to those.", n_good)
    elif n_good + n_marginal >= 6:
        logging.info("VERDICT: marginal. Only %d clearly usable; including the "
                     "%d marginal samples is possible but needs a per-class "
                     "callable mask and depth matching.", n_good, n_marginal)
    else:
        logging.info("VERDICT: too few samples with germline DNA (%d usable). "
                     "Do not build the polymorphism analysis on this dataset — "
                     "ask Jianbin for gonad-derived WGS instead.", n_good)
    logging.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
