#!/usr/bin/env python3
"""
het_by_region.py — per-sample heterozygosity split by region class.

Why this exists
---------------
The duplicate check reported per-sample heterozygosity spanning 0.215 to 0.534
across ten libraries, and the pairwise table showed something sharper: IBS0 is
depressed against EVERY partner for two libraries rather than for particular
pairs. Kinship is pairwise, so no pedigree explains a worm that looks close to
all nine others. A sample that rarely emits a confident homozygous-alternate
call does explain it, and that is what excess apparent heterozygosity looks
like from the outside.

Genome-wide excess heterozygosity is a nuisance. Excess heterozygosity
CONCENTRATED IN ELIMINATED REGIONS is a threat to the paper's result, because
collapsed paralogy is already known to be ~12x worse there, and paralogous
copies differ at fixed amino-acid positions — so the spurious variants they
create are enriched for nonsynonymous changes and inflate piN specifically in
the eliminated class.

This script separates the two. It asks, for each worm, whether its
heterozygosity excess is a property of the library or a property of the region
class, which is the difference between a cohort footnote and a confounded
headline.

What is measured, per sample and region class
---------------------------------------------
* het, hom-ref and hom-alt rates over called genotypes
* the eliminated:retained het ratio, and each sample's deviation from the
  cohort median of that ratio -- the actual quantity of interest, since it is
  internally controlled: a library that is globally het-inflated has a high het
  rate in both classes and a NORMAL ratio
* mean allele balance at heterozygous calls, if the callset carries AD. Real
  heterozygotes centre on 0.5; collapsed paralogs skew low because the
  "alternate allele" is a minority second copy

Region assignment reuses the CDS intervals from the same GFF and
gene_callable.tsv that pi_by_region.py uses, so a site is retained or
eliminated here exactly when it is there. It builds only a flat interval map,
not the full codon index, so it is far lighter than pi_by_region.py.

Example
-------
    python het_by_region.py \\
        --vcf vcf/filtered_miss4/snps_invariant.bcf \\
        --gff As_miniprot.gff \\
        --gene-callable targets/gene_callable.tsv \\
        --outdir results/het_by_region
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pi_by_region import read_cds

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
    from paper_colors import RETAINED, ELIMINATED       # type: ignore
except Exception:
    RETAINED, ELIMINATED = "#2C6FBB", "#C0392B"

plt.rcParams.update({
    "font.family": "Arial", "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelsize": 17, "xtick.labelsize": 11, "ytick.labelsize": 13,
    "legend.fontsize": 13, "figure.dpi": 150,
})

REGIONS = ("retained", "eliminated")
RCODE = {"retained": 0, "eliminated": 1}


class RegionMap:
    """
    Position -> region class, as sorted per-chromosome interval arrays.

    Only CDS intervals of analysable genes are loaded, matching the sites
    pi_by_region.py actually uses. Overlapping intervals from two genes of
    different class would be ambiguous, so they are dropped and counted rather
    than silently assigned to whichever sorted first.
    """

    def __init__(self):
        self.starts: dict[str, np.ndarray] = {}
        self.ends: dict[str, np.ndarray] = {}
        self.code: dict[str, np.ndarray] = {}

    def build(self, genes: dict, region_of: dict) -> int:
        per_chrom: dict[str, list] = {}
        for gene, g in genes.items():
            r = region_of.get(gene)
            if r not in RCODE:
                continue
            for s, e in g["cds"]:
                per_chrom.setdefault(g["chrom"], []).append((s, e, RCODE[r]))
        n_conflict = 0
        for chrom, ivs in per_chrom.items():
            ivs.sort()
            merged: list[list] = []
            for s, e, c in ivs:
                if merged and s < merged[-1][1]:
                    if merged[-1][2] != c:
                        # two classes claim the same bases: unusable
                        n_conflict += 1
                        merged[-1][1] = min(merged[-1][1], s)
                        continue
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e, c])
            merged = [m for m in merged if m[1] > m[0]]
            self.starts[chrom] = np.array([m[0] for m in merged], np.int64)
            self.ends[chrom] = np.array([m[1] for m in merged], np.int64)
            self.code[chrom] = np.array([m[2] for m in merged], np.int8)
        total = sum(len(v) for v in self.starts.values())
        logging.info("Region map: %d interval(s) on %d contig(s)%s", total,
                     len(self.starts),
                     f"; {n_conflict} class-conflicting overlap(s) trimmed"
                     if n_conflict else "")
        return total

    def lookup(self, chrom: str, pos0: int) -> int:
        st = self.starts.get(chrom)
        if st is None or not st.size:
            return -1
        i = int(np.searchsorted(st, pos0, "right")) - 1
        if i < 0:
            return -1
        return int(self.code[chrom][i]) if pos0 < self.ends[chrom][i] else -1


def vcf_samples(vcf: Path) -> list[str]:
    r = subprocess.run(["bash", "-c", f"bcftools query -l {vcf}"],
                       capture_output=True, text=True, check=True)
    return [s for s in r.stdout.split() if s]


def has_ad(vcf: Path) -> bool:
    r = subprocess.run(["bash", "-c", f"bcftools view -h {vcf} | grep -c "
                                      f"'^##FORMAT=<ID=AD,' || true"],
                       capture_output=True, text=True)
    return r.stdout.strip() not in ("", "0")


def stream(vcf: Path, rmap: RegionMap, n: int, want_ad: bool,
           variants_only: bool) -> dict:
    """
    One pass, tallying genotype classes per sample per region.

    ``variants_only`` restricts to segregating sites. That is the right
    denominator for a het RATE that is comparable to the duplicate check, and
    it also keeps the stream ~100x shorter than walking every invariant site.
    """
    fields = "%CHROM\\t%POS[\\t%GT" + ("\\t%AD" if want_ad else "") + "]\\n"
    filt = " -i 'N_ALT=1 && TYPE=\"snp\"'" if variants_only else ""
    cmd = f"bcftools query -f '{fields}'{filt} {vcf}"
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 22)

    per = 2 if want_ad else 1
    # keyed by the integer region code that RegionMap.lookup returns, not by
    # region name -- these two must not drift apart
    counts = {c: {k: np.zeros(n, np.int64)
                  for k in ("called", "het", "hom_ref", "hom_alt", "missing")}
              for c in RCODE.values()}
    ab_sum = {c: np.zeros(n) for c in RCODE.values()}
    ab_n = {c: np.zeros(n, np.int64) for c in RCODE.values()}

    n_lines = n_mapped = 0
    for line in proc.stdout:
        n_lines += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 + per * n:
            continue
        code = rmap.lookup(f[0], int(f[1]) - 1)
        if code < 0:
            continue
        n_mapped += 1
        c = counts[code]
        for k in range(n):
            gt = f[2 + per * k].replace("|", "/")
            if gt in ("./.", ".", ""):
                c["missing"][k] += 1
                continue
            a = gt.split("/")
            if len(a) != 2 or "." in a:
                c["missing"][k] += 1
                continue
            c["called"][k] += 1
            if a[0] != a[1]:
                c["het"][k] += 1
                if want_ad:
                    ad = f[3 + per * k].split(",")
                    if len(ad) >= 2:
                        try:
                            ref_d, alt_d = int(ad[0]), int(ad[1])
                        except ValueError:
                            continue
                        if ref_d + alt_d > 0:
                            ab_sum[code][k] += alt_d / (ref_d + alt_d)
                            ab_n[code][k] += 1
            elif a[0] == "0":
                c["hom_ref"][k] += 1
            else:
                c["hom_alt"][k] += 1

    proc.wait()
    if proc.returncode != 0:
        logging.error("bcftools query failed:\n%s", proc.stderr.read()[-2000:])
        sys.exit(proc.returncode)
    logging.info("Streamed %d record(s); %d fell in classified CDS", n_lines,
                 n_mapped)
    if n_mapped == 0:
        logging.error("No sites mapped into the region map. Check that the "
                      "GFF contig names match the VCF.")
        sys.exit(1)
    return {"counts": counts, "ab_sum": ab_sum, "ab_n": ab_n}


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Per-sample heterozygosity by region class.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf", required=True, type=Path)
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--gene-callable", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--all-sites", action="store_true",
                   help="Include invariant sites in the denominator. Slower "
                        "and gives a per-base rather than per-SNP het rate; "
                        "the default (variants only) is what the duplicate "
                        "check reported, so leave it off to compare with that.")
    p.add_argument("--no-allele-balance", action="store_true",
                   help="Skip AD parsing even if the callset carries it")
    p.add_argument("--ratio-tolerance", type=float, default=0.25,
                   help="Flag a sample whose eliminated:retained het ratio "
                        "differs from the cohort median by more than this "
                        "fraction. The RATIO is the quantity that matters: a "
                        "globally het-inflated library has a normal ratio and "
                        "cannot bias the class contrast, while a library with "
                        "a high ratio is inflating one class only.")
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

    samples = vcf_samples(args.vcf.expanduser())
    n = len(samples)
    logging.info("Samples: %d", n)

    gc = pd.read_csv(args.gene_callable.expanduser(), sep="\t")
    gc = gc[gc["pass"].astype(bool) & gc.region.isin(REGIONS)]
    region_of = dict(zip(gc.gene, gc.region))
    genes = read_cds(args.gff.expanduser(), set(gc.gene))
    rmap = RegionMap()
    if rmap.build(genes, region_of) == 0:
        logging.error("Empty region map")
        return 1

    want_ad = (not args.no_allele_balance) and has_ad(args.vcf.expanduser())
    logging.info("Allele balance: %s", "on (FORMAT/AD present)" if want_ad
                 else "off")

    res = stream(args.vcf.expanduser(), rmap, n, want_ad, not args.all_sites)
    counts, ab_sum, ab_n = res["counts"], res["ab_sum"], res["ab_n"]

    rows = []
    for k, s in enumerate(samples):
        row = {"sample": s}
        for name, code in RCODE.items():
            c = counts[code]
            called = c["called"][k]
            row[f"called_{name}"] = int(called)
            row[f"missing_{name}"] = int(c["missing"][k])
            row[f"het_{name}"] = c["het"][k] / called if called else np.nan
            row[f"homalt_{name}"] = (c["hom_alt"][k] / called if called
                                     else np.nan)
            row[f"ab_{name}"] = (ab_sum[code][k] / ab_n[code][k]
                                 if ab_n[code][k] else np.nan)
        row["het_ratio_eli_ret"] = (row["het_eliminated"] / row["het_retained"]
                                    if row["het_retained"] else np.nan)
        rows.append(row)
    df = pd.DataFrame(rows)

    med = float(df.het_ratio_eli_ret.median())
    df["ratio_vs_cohort"] = df.het_ratio_eli_ret / med
    df = df.sort_values("het_ratio_eli_ret", ascending=False)
    df.to_csv(outdir / "het_by_region.tsv", sep="\t", index=False)

    logging.info("=" * 82)
    logging.info("%-14s %9s %9s %9s %9s %8s", "sample", "het_ret", "het_eli",
                 "eli/ret", "vs cohort", "AB_eli")
    for _, r in df.iterrows():
        logging.info("%-14s %9.4f %9.4f %9.3f %9.2f %8s", r["sample"],
                     r.het_retained, r.het_eliminated, r.het_ratio_eli_ret,
                     r.ratio_vs_cohort,
                     f"{r.ab_eliminated:.3f}" if r.ab_eliminated ==
                     r.ab_eliminated else "n/a")

    logging.info("")
    logging.info("Cohort median eliminated:retained het ratio = %.3f", med)
    logging.info("Genome-wide het rate spread: %.4f - %.4f (%.1fx)",
                 df.het_retained.min(), df.het_retained.max(),
                 df.het_retained.max() / df.het_retained.min()
                 if df.het_retained.min() else np.nan)

    off = df[(df.ratio_vs_cohort - 1).abs() > args.ratio_tolerance]
    logging.info("")
    if off.empty:
        logging.info("VERDICT: every worm's eliminated:retained het ratio is "
                     "within %.0f%% of the cohort median. Whatever drives the "
                     "%.1fx spread in absolute heterozygosity acts on both "
                     "region classes alike, so it cannot generate the class "
                     "contrast — it inflates numerator and denominator "
                     "together. Report the spread as a library-quality "
                     "observation, not a confound.",
                     100 * args.ratio_tolerance,
                     df.het_retained.max() / df.het_retained.min()
                     if df.het_retained.min() else float("nan"))
    else:
        logging.warning("VERDICT: %d sample(s) carry a class-SPECIFIC het "
                        "excess or deficit:", len(off))
        for _, r in off.iterrows():
            logging.warning("  %-14s eli/ret = %.3f (%.2fx the cohort median "
                            "%.3f)", r["sample"], r.het_ratio_eli_ret,
                            r.ratio_vs_cohort, med)
        logging.warning("This is the confound that matters. Re-run "
                        "jackknife_individuals.py with --drop-set naming these "
                        "libraries and confirm the piN/piS fold survives "
                        "without them before quoting it.")

    if want_ad:
        logging.info("")
        logging.info("Allele balance at heterozygous calls (0.5 = clean "
                     "diploid het; low = the alternate allele is a minority "
                     "read stack, the collapsed-paralog signature):")
        logging.info("  retained   %.4f - %.4f", df.ab_retained.min(),
                     df.ab_retained.max())
        logging.info("  eliminated %.4f - %.4f", df.ab_eliminated.min(),
                     df.ab_eliminated.max())
        gap = (df.ab_retained - df.ab_eliminated).median()
        logging.info("  median within-sample retained - eliminated = %+.4f",
                     gap)
        if gap > 0.02:
            logging.info("  Eliminated hets sit further from 0.5 in the same "
                         "libraries, which is collapsed paralogy rather than a "
                         "per-library problem — already handled by the "
                         "allele-balance screen in filter_genotypes.py.")

    # --- figure -------------------------------------------------------------
    d = df.sort_values("het_retained")
    x = np.arange(len(d))
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.8))
    w = 0.38
    axes[0].bar(x - w / 2, d.het_retained, w, color=RETAINED,
                edgecolor="black", linewidth=0.5, label="Retained")
    axes[0].bar(x + w / 2, d.het_eliminated, w, color=ELIMINATED,
                edgecolor="black", linewidth=0.5, label="Eliminated")
    axes[0].set_ylabel("Heterozygous fraction")
    axes[0].legend(frameon=False)

    axes[1].axhline(med, color="grey", ls="--", lw=1.5,
                    label=f"Cohort median ({med:.2f})")
    axes[1].axhspan(med * (1 - args.ratio_tolerance),
                    med * (1 + args.ratio_tolerance), color="grey", alpha=0.15)
    axes[1].scatter(x, d.het_ratio_eli_ret, s=90, color=ELIMINATED,
                    edgecolor="black", linewidth=0.6, zorder=3)
    axes[1].set_ylabel("Eliminated : retained het")
    axes[1].legend(frameon=False)

    for a in axes:
        a.set_xticks(x)
        a.set_xticklabels(d["sample"], rotation=90)
        a.set_ylim(bottom=0)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(outdir / f"S_het_by_region.{ext}", bbox_inches="tight")
    plt.close(fig)
    logging.info("")
    logging.info("Figure: S_het_by_region.{png,svg,pdf}")
    logging.info("Export from the SVG — Arial is absent on Isaac nodes.")
    logging.info("=" * 82)
    return 0


if __name__ == "__main__":
    sys.exit(main())
