#!/usr/bin/env python3
"""
make_callable_targets.py — build the callable CDS target set, and report how
many genes of each region class survive it.

This runs BEFORE genotyping, for two reasons.

**It bounds the compute.** An all-sites VCF over 276 Mb x 10 samples is large
and slow. Every statistic we want (piN/piS, MK, LoF load, 4-fold neutral
control) lives inside coding sequence, so genotyping is restricted to CDS.

**It answers the question that can invalidate the analysis.** Eliminated
sequence is repeat-rich and will be masked harder than retained sequence. If
masking removes most eliminated CDS, the contrast is not between two gene
classes but between a well-observed class and a badly-observed one. That has
to be quantified here, in gene counts, not discovered later as a strange
result. The per-class survivorship table at the end is the output that matters.

What it does
------------
  1. Read CDS features from the miniprot GFF, grouped by mRNA (the same key,
     ``gene1``, that ``kaks_annotated.tsv`` uses, so the tables join).
  2. Assign each gene to retained / eliminated by the region BED.
  3. Restrict to autosomes — X has a different Ne and this cohort is mixed-sex.
  4. Subtract the low-mappability mask.
  5. Drop genes with too little callable CDS left to estimate anything.

Note on gene models: 18,588 miniprot models merge into ~15,652 loci, so model
counts over-count genes by roughly 16%. Keys are kept at model level anyway
because that is what the Ka/Ks table is keyed on; the redundancy affects counts
quoted in text, not the join.

Example
-------
    python make_callable_targets.py \\
        --gff ../As_miniprot.gff \\
        --eliminated-bed $DB_DIR/ascaris/v50/AG.v50.eliminated_strict.bed \\
        --fai $DB_DIR/ascaris/v50/AG.v50.fa.fai \\
        --mask-bed $POPGEN_DIR/mask/lowmap_k150_e2.bed \\
        --outdir $POPGEN_DIR/targets
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
    from paper_colors import RETAINED, ELIMINATED       # type: ignore
except Exception:
    RETAINED, ELIMINATED = "#2C6FBB", "#C0392B"

plt.rcParams.update({
    "font.family": "Arial", "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelsize": 20, "xtick.labelsize": 16, "ytick.labelsize": 16,
    "legend.fontsize": 15, "figure.dpi": 150,
})

Interval = tuple[int, int]


# ---------------------------------------------------------------------------
# Interval helpers
# ---------------------------------------------------------------------------

def merge(ivs: list[Interval]) -> list[Interval]:
    """Merge overlapping/adjacent intervals. Input need not be sorted."""
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [ivs[0]]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def subtract(a: list[Interval], b: list[Interval]) -> list[Interval]:
    """
    a minus b, both merged and sorted. Linear sweep, not O(n*m): the mask can
    hold millions of intervals and a per-gene rescan of it would dominate.
    """
    if not b:
        return list(a)
    out: list[Interval] = []
    j = 0
    for s, e in a:
        cur = s
        while j > 0 and b[j - 1][1] > cur:
            j -= 1
        k = j
        while k < len(b) and b[k][0] < e:
            bs, be = b[k]
            if be <= cur:
                k += 1
                continue
            if bs > cur:
                out.append((cur, min(bs, e)))
            cur = max(cur, be)
            if cur >= e:
                break
            k += 1
        if cur < e:
            out.append((cur, e))
    return out


def total(ivs: list[Interval]) -> int:
    return sum(e - s for s, e in ivs)


def contains(ivs: list[Interval], pos: int) -> bool:
    """Membership of a point in a merged interval list (binary search)."""
    lo, hi = 0, len(ivs) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid]
        if pos < s:
            hi = mid - 1
        elif pos >= e:
            lo = mid + 1
        else:
            return True
    return False


def read_bed(path: Path) -> dict[str, list[Interval]]:
    raw: dict[str, list[Interval]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return {c: merge(v) for c, v in raw.items()}


# ---------------------------------------------------------------------------
# GFF
# ---------------------------------------------------------------------------

def read_cds(gff: Path) -> dict[str, dict]:
    """
    Collect CDS intervals per mRNA from a miniprot GFF.

    Returns {mrna_id: {chrom, strand, cds: [(start, end), ...]}} in 0-based
    half-open coordinates. GFF is 1-based inclusive, hence the start-1.
    """
    attr_id = re.compile(r"(?:^|;)ID=([^;]+)")
    attr_parent = re.compile(r"(?:^|;)Parent=([^;]+)")

    genes: dict[str, dict] = {}
    n_cds = 0
    opener = gzip.open if str(gff).endswith(".gz") else open
    with opener(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                m = attr_id.search(f[8])
                if m:
                    genes.setdefault(m.group(1), {"chrom": f[0], "strand": f[6],
                                                  "cds": []})
            elif f[2] == "CDS":
                m = attr_parent.search(f[8])
                if not m:
                    continue
                gid = m.group(1)
                g = genes.setdefault(gid, {"chrom": f[0], "strand": f[6],
                                           "cds": []})
                g["cds"].append((int(f[3]) - 1, int(f[4])))
                n_cds += 1

    genes = {k: v for k, v in genes.items() if v["cds"]}
    for g in genes.values():
        g["cds"] = merge(g["cds"])       # collapse isoform/frame duplication
    logging.info("GFF: %d mRNA(s) with CDS, %d CDS features", len(genes), n_cds)
    return genes


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build callable CDS targets and report per-class survivorship.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--gff", required=True, type=Path,
                   help="miniprot GFF in germline (v50) coordinates")
    p.add_argument("--eliminated-bed", required=True, type=Path)
    p.add_argument("--fai", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--mask-bed", type=Path,
                   help="Low-mappability BED from mappability_mask.py. Omit to "
                        "measure survivorship of the unmasked CDS set.")
    p.add_argument("--min-callable-cds", type=int, default=300,
                   help="Drop genes with less callable CDS than this (bp). "
                        "300 bp = 100 codons; below that a per-gene rate is "
                        "noise, though genes are still pooled by class.")
    p.add_argument("--region-rule", default="midpoint",
                   choices=["midpoint", "majority"],
                   help="How a gene is assigned to a region. 'midpoint' matches "
                        "annotate_regions.py in the Ka/Ks pipeline; use it "
                        "unless you have a reason not to.")
    p.add_argument("--autosome-regex", default=r"^chr\d+$",
                   help="Sequences to keep. Default excludes chrX*, chrM, chrUn*.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    keep_seq = re.compile(args.autosome_regex)
    lengths = {l.split("\t")[0]: int(l.split("\t")[1])
               for l in Path(args.fai).expanduser().read_text().splitlines()
               if l.strip()}
    autosomes = sorted(c for c in lengths if keep_seq.match(c))
    logging.info("Autosomes kept (%s): %d sequences, %.1f Mb",
                 args.autosome_regex, len(autosomes),
                 sum(lengths[c] for c in autosomes) / 1e6)

    elim = read_bed(args.eliminated_bed.expanduser())
    mask = read_bed(args.mask_bed.expanduser()) if args.mask_bed else {}
    if args.mask_bed:
        logging.info("Mask: %.2f Mb over %d sequences",
                     sum(total(v) for v in mask.values()) / 1e6, len(mask))
    else:
        logging.warning("No --mask-bed given; reporting UNMASKED survivorship")

    genes = read_cds(args.gff.expanduser())

    rows = []
    target_ivs: dict[str, list[Interval]] = {}
    for gid, g in genes.items():
        chrom = g["chrom"]
        if not keep_seq.match(chrom):
            continue
        cds = g["cds"]
        cds_bp = total(cds)
        if cds_bp == 0:
            continue

        e_ivs = elim.get(chrom, [])
        if args.region_rule == "midpoint":
            mid = (cds[0][0] + cds[-1][1]) // 2
            region = "eliminated" if contains(e_ivs, mid) else "retained"
        else:
            overlap = cds_bp - total(subtract(cds, e_ivs))
            region = "eliminated" if overlap * 2 > cds_bp else "retained"

        callable_ivs = subtract(cds, mask.get(chrom, []))
        callable_bp = total(callable_ivs)

        rows.append({
            "gene": gid, "chrom": chrom, "start": cds[0][0], "end": cds[-1][1],
            "region": region, "cds_bp": cds_bp, "callable_bp": callable_bp,
            "callable_frac": callable_bp / cds_bp,
            "pass": callable_bp >= args.min_callable_cds,
        })
        if callable_bp >= args.min_callable_cds:
            target_ivs.setdefault(chrom, []).extend(callable_ivs)

    df = pd.DataFrame(rows)
    if df.empty:
        logging.error("No genes on the kept sequences — check --autosome-regex "
                      "against the GFF's sequence names")
        return 1

    df.to_csv(outdir / "gene_callable.tsv", sep="\t", index=False)

    # BED of everything genotyping should look at
    bed = outdir / "callable_cds_targets.bed"
    with bed.open("w") as fh:
        for chrom in sorted(target_ivs):
            for s, e in merge(target_ivs[chrom]):
                fh.write(f"{chrom}\t{s}\t{e}\n")

    # --- survivorship ------------------------------------------------------
    logging.info("=" * 70)
    logging.info("Callable CDS survivorship by region class (autosomes only)")
    logging.info("%-12s %7s %7s %10s %10s %8s", "class", "genes", "pass",
                 "CDS Mb", "callable", "median f")
    summary = []
    for region in ("retained", "eliminated"):
        sub = df[df["region"] == region]
        if sub.empty:
            continue
        logging.info("%-12s %7d %7d %10.2f %10.2f %8.3f",
                     region, len(sub), int(sub["pass"].sum()),
                     sub["cds_bp"].sum() / 1e6, sub["callable_bp"].sum() / 1e6,
                     sub["callable_frac"].median())
        summary.append({
            "region": region, "genes": len(sub), "genes_pass": int(sub["pass"].sum()),
            "cds_bp": int(sub["cds_bp"].sum()),
            "callable_bp": int(sub["callable_bp"].sum()),
            "median_callable_frac": float(sub["callable_frac"].median()),
            "pass_rate": float(sub["pass"].mean()),
        })
    pd.DataFrame(summary).to_csv(outdir / "callable_summary.tsv",
                                 sep="\t", index=False)

    ret = df[df["region"] == "retained"]
    eli = df[df["region"] == "eliminated"]
    if not eli.empty and not ret.empty:
        logging.info("")
        logging.info("Eliminated genes retained at >= %d bp callable CDS: "
                     "%d of %d (%.1f%%)", args.min_callable_cds,
                     int(eli["pass"].sum()), len(eli), 100 * eli["pass"].mean())
        logging.info("Retained   genes retained at >= %d bp callable CDS: "
                     "%d of %d (%.1f%%)", args.min_callable_cds,
                     int(ret["pass"].sum()), len(ret), 100 * ret["pass"].mean())
        loss_ratio = ((1 - eli["callable_frac"].median()) /
                      max(1 - ret["callable_frac"].median(), 1e-9))
        logging.info("Eliminated CDS loses %.1fx as much sequence to masking "
                     "as retained", loss_ratio)
        if eli["pass"].mean() < 0.5:
            logging.warning("Fewer than half of eliminated genes survive. The "
                            "contrast is now between a well-observed and a "
                            "poorly-observed class — report the surviving gene "
                            "counts prominently and check whether survivors are "
                            "a biased subset (e.g. by GC or expression).")
    logging.info("Targets: %s", bed)
    logging.info("=" * 70)

    # --- figure ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    parts = ax.violinplot([ret["callable_frac"], eli["callable_frac"]],
                          positions=[0, 1], widths=0.75,
                          showmedians=True, showextrema=False)
    for body, colour in zip(parts["bodies"], (RETAINED, ELIMINATED)):
        body.set_facecolor(colour)
        body.set_alpha(0.75)
        body.set_edgecolor("black")
    parts["cmedians"].set_color("black")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Retained\n(n = {len(ret):,})",
                        f"Eliminated\n(n = {len(eli):,})"])
    ax.set_ylabel("Callable fraction of CDS")
    ax.set_ylim(0, 1.02)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(outdir / f"callable_fraction_by_region.{ext}",
                    bbox_inches="tight")
    plt.close(fig)
    logging.info("Figure: callable_fraction_by_region.{png,svg,pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
