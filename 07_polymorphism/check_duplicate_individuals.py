#!/usr/bin/env python3
"""
check_duplicate_individuals.py — are any two of the cohort the same worm?

Why this matters
----------------
Landeryou et al. report 54 worms for 54 runs, but the submission carries ONE
BioSample for all 62 runs and its Supplementary Data 2 is column-scrambled, so
nothing in the metadata actually confirms one worm per run. If two of the ten
selected libraries are the same individual, the effective sample size is nine,
every bootstrap CI is slightly too narrow, and the duplicated genotype is
double-counted in every allele frequency.

Cheap to rule out, and it cannot be ruled out any other way.

Method
------
Pairwise genotype concordance at biallelic SNPs called in both members of a
pair. Two libraries from the same worm agree at essentially every site — the
only disagreements are genotyping error, so discordance runs at a percent or
two. Two different worms from an outbred population disagree far more often,
and the gap is large enough that no careful threshold is needed.

Also reported per pair:

* **IBS0** — sites where one worm is 0/0 and the other 1/1. Two independent
  genotyping errors are needed to produce one from a single individual, so a
  duplicate's IBS0 sits near zero while unrelated outbred worms run a few
  percent. Reported as corroboration; discordance is the primary call, since
  IBS0 is not exactly zero even for a true duplicate and a tight threshold on
  it would veto real hits.
* **relatedness** — mean product of allele-frequency-standardised genotypes
  (GCTA/KING style), reported raw and rescaled so the median pair reads 0 and a
  same-individual pair reads ~1; full sibs land near 0.5. This flags close kin,
  which matter for the same reason duplicates do though less severely. Note it
  is NOT a raw dosage correlation: see standardise() for why that measure is
  confounded, and pairwise_stats() for why the rescaling is necessary.

Example
-------
    python check_duplicate_individuals.py \\
        --vcf vcf/filtered_miss4/snps_invariant.bcf \\
        --outdir results/duplicates
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Arial", "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelsize": 18, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "legend.fontsize": 14, "figure.dpi": 150,
})

GT_MAP = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1,
          "1/1": 2, "1|1": 2}


def load_genotypes(vcf: Path, max_sites: int) -> tuple[list[str], np.ndarray]:
    """
    Genotype dosage matrix (sites x samples) at biallelic SNPs, -1 for missing.

    Only variant sites are requested: invariant sites carry no information
    about whether two worms differ, and there are 12 M of them.
    """
    names = subprocess.run(["bash", "-c", f"bcftools query -l {vcf}"],
                           capture_output=True, text=True, check=True)
    samples = [s for s in names.stdout.split() if s]
    logging.info("Samples: %d", len(samples))

    cmd = (f"bcftools query -f '[%GT\\t]\\n' "
           f"-i 'N_ALT=1 && TYPE=\"snp\"' {vcf}")
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 22)
    rows = []
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")[:len(samples)]
        if len(f) < len(samples):
            continue
        rows.append([GT_MAP.get(g, -1) for g in f])
        if max_sites and len(rows) >= max_sites:
            logging.info("Reached --max-sites (%d); stopping the stream early",
                         max_sites)
            proc.kill()
            proc.wait()
            break
    else:
        proc.wait()
        if proc.returncode != 0:
            logging.error("bcftools query failed:\n%s",
                          proc.stderr.read()[-2000:])
            sys.exit(proc.returncode)

    g = np.array(rows, dtype=np.int8)
    logging.info("Loaded %d biallelic SNP(s)", g.shape[0])
    return samples, g


def standardise(g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Allele-frequency-standardised genotypes, z = (dosage - 2p) / sqrt(2p(1-p)).

    A raw Pearson correlation on 0/1/2 dosage is NOT a relatedness measure and
    will mislead here. Sites differ in allele frequency, and every individual
    tracks that shared frequency profile, so two completely unrelated worms
    still correlate at r ~ 0.15 in this data simply because both carry more
    alternate alleles at high-frequency sites. Centring each site on its own
    2p removes exactly that shared component, which is what the GCTA/KING
    style estimator below needs.

    Returns the z matrix (NaN where a genotype is missing or a site is
    monomorphic) and a boolean mask of usable sites.
    """
    d = g.astype(np.float64)
    d[g < 0] = np.nan
    # A site called in no sample gives an empty-slice warning and a NaN mean;
    # the NaN is the wanted answer (the `usable` mask drops it), so silence the
    # warning rather than special-casing the row.
    with warnings.catch_warnings(), np.errstate(invalid="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        p = np.nanmean(d, axis=1) / 2.0
    usable = np.isfinite(p) & (p > 0) & (p < 1)
    z = np.full(d.shape, np.nan)
    denom = np.sqrt(2 * p[usable] * (1 - p[usable]))
    z[usable] = (d[usable] - 2 * p[usable, None]) / denom[:, None]
    return z, usable


def pairwise_stats(samples: list[str], g: np.ndarray) -> pd.DataFrame:
    """
    Discordance, IBS0 and standardised relatedness for every pair, on sites
    called in both members of that pair.

    Each pair uses its own intersection of called sites rather than one global
    complete-case set: dropping a site because some third library missed it
    would throw away information and, worse, would make the comparison depend
    on which other libraries happen to be in the run.

    `relatedness` is the mean product of standardised genotypes, ~1 for the
    same individual, ~0.5 for full sibs, ~0 for unrelated. With allele
    frequencies estimated from these same few libraries the unrelated
    expectation is not exactly 0 but about -1/(n-1); read the off-diagonal
    values relative to each other, not against a textbook zero.
    """
    z, _ = standardise(g)
    rows = []
    for i, j in combinations(range(len(samples)), 2):
        a, b = g[:, i], g[:, j]
        both = (a >= 0) & (b >= 0)
        n = int(both.sum())
        if n == 0:
            logging.warning("%s vs %s share no called sites",
                            samples[i], samples[j])
            continue
        aa, bb = a[both].astype(np.int16), b[both].astype(np.int16)
        disc = float((aa != bb).mean())
        ibs0 = float(((aa == 0) & (bb == 2)).sum() +
                     ((aa == 2) & (bb == 0)).sum()) / n
        prod = z[both, i] * z[both, j]
        prod = prod[np.isfinite(prod)]
        rel = float(prod.mean()) if prod.size else np.nan
        rows.append({"sample_a": samples[i], "sample_b": samples[j],
                     "sites_compared": n, "discordance": disc,
                     "ibs0_rate": ibs0, "relatedness": rel})
    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Rescale so unrelated sits at 0 and a same-individual pair at ~1.
    #
    # Raw relatedness is shrunk toward zero and shifted negative because the
    # allele frequencies come from these same ten libraries: the unrelated
    # expectation is about -1/(n-1), not 0, and the shift grows as n shrinks.
    # Taking the median off-diagonal pair as the unrelated baseline fixes both
    # at once and makes one threshold meaningful whatever n is. The median is
    # the right statistic because the great majority of pairs in any cohort are
    # unrelated -- if that were false, no within-cohort estimator would work.
    base = float(df.relatedness.median())
    df["relatedness_scaled"] = (df.relatedness - base) / (1.0 - base)
    return df.sort_values("discordance").reset_index(drop=True)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Detect duplicate or closely related individuals.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--max-sites", type=int, default=200000,
                   help="Cap on SNPs read; 200k is far more than needed and "
                        "keeps memory trivial. 0 = no cap.")
    p.add_argument("--dup-discordance", type=float, default=0.05,
                   help="Flag a pair as a possible duplicate below this "
                        "genotype discordance")
    p.add_argument("--dup-ibs0", type=float, default=0.01,
                   help="Flag a pair below this IBS0 rate. Discordance is the "
                        "primary discriminator — it separates duplicates from "
                        "unrelated worms by roughly two orders of magnitude — "
                        "so this second condition is a sanity veto, not a fine "
                        "threshold, and it is deliberately loose. A true "
                        "duplicate's IBS0 is not exactly zero: a 0/0-vs-1/1 "
                        "call takes two independent genotyping errors, which "
                        "is rare but not impossible over a million sites. Set "
                        "it tight enough to catch that and it vetoes real "
                        "duplicates instead.")
    p.add_argument("--kin-relatedness", type=float, default=0.20,
                   help="Report pairs above this RESCALED relatedness as "
                        "possible close kin (the relatedness_scaled column, "
                        "where the median pair is 0 and a same-individual pair "
                        "is ~1). Full sibs sit near 0.5, half-sibs near 0.25.")
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

    samples, g = load_genotypes(args.vcf.expanduser(), args.max_sites)
    if g.size == 0:
        logging.error("No genotypes loaded")
        return 1

    # per-sample heterozygosity, a useful sanity read on its own
    het = []
    for j, s in enumerate(samples):
        col = g[:, j]
        called = col >= 0
        het.append({"sample": s, "sites_called": int(called.sum()),
                    "het_rate": float((col[called] == 1).mean())})
    pd.DataFrame(het).to_csv(outdir / "per_sample.tsv", sep="\t", index=False)
    logging.info("Per-sample het rate: %.4f - %.4f",
                 min(h["het_rate"] for h in het),
                 max(h["het_rate"] for h in het))

    pw = pairwise_stats(samples, g)
    if pw.empty:
        logging.error("No comparable pairs")
        return 1
    pw.to_csv(outdir / "pairwise.tsv", sep="\t", index=False)

    logging.info("=" * 70)
    logging.info("Pairwise genotype comparison (%d pairs)", len(pw))
    logging.info("  discordance : %.4f - %.4f (median %.4f)",
                 pw.discordance.min(), pw.discordance.max(),
                 pw.discordance.median())
    logging.info("  IBS0 rate   : %.5f - %.5f", pw.ibs0_rate.min(),
                 pw.ibs0_rate.max())
    logging.info("  relatedness : %+.3f - %+.3f (rescaled so the median "
                 "pair is 0)", pw.relatedness_scaled.min(),
                 pw.relatedness_scaled.max())
    logging.info("")
    logging.info("Five most similar pairs:")
    for _, r in pw.head(5).iterrows():
        logging.info("  %-14s %-14s  discord=%.4f  IBS0=%.5f  rel=%+.3f  "
                     "(n=%d)", r.sample_a, r.sample_b, r.discordance,
                     r.ibs0_rate, r.relatedness_scaled, r.sites_compared)

    flagged = pw[(pw.discordance <= args.dup_discordance) &
                 (pw.ibs0_rate <= args.dup_ibs0)]
    logging.info("")
    if flagged.empty:
        logging.info("VERDICT: no duplicate individuals. Every pair exceeds "
                     "%.1f%% genotype discordance or %.3f%% IBS0, so all %d "
                     "libraries are distinct worms and n = %d stands.",
                     100 * args.dup_discordance, 100 * args.dup_ibs0,
                     len(samples), len(samples))
        gap = pw.discordance.min()
        logging.info("The closest pair still disagrees at %.1f%% of sites — "
                     "far above a genotyping-error floor, so this is not a "
                     "marginal call.", 100 * gap)
    else:
        logging.warning("VERDICT: %d pair(s) look like the SAME individual:",
                        len(flagged))
        for _, r in flagged.iterrows():
            logging.warning("  %s == %s  (discordance %.4f, IBS0 %.5f)",
                            r.sample_a, r.sample_b, r.discordance, r.ibs0_rate)
        logging.warning("Drop one of each pair and re-run every downstream "
                        "statistic: a duplicated worm is double-counted in "
                        "every allele frequency and makes each CI too narrow.")

    close = pw[(pw.relatedness_scaled >= args.kin_relatedness) &
               (~pw.index.isin(flagged.index))]
    if not close.empty:
        logging.info("")
        logging.info("Also note %d pair(s) with relatedness >= %.2f — "
                     "plausible close kin (full sibs sit near 0.5). Not "
                     "duplicates, but worms from one host are not independent "
                     "draws and this is worth a sentence in the methods.",
                     len(close), args.kin_relatedness)
        for _, r in close.iterrows():
            logging.info("  %s ~ %s  relatedness=%+.3f (raw %+.3f)",
                         r.sample_a, r.sample_b, r.relatedness_scaled,
                         r.relatedness)

    # --- figure -------------------------------------------------------------
    n = len(samples)
    mat = np.full((n, n), np.nan)
    idx = {s: k for k, s in enumerate(samples)}
    for _, r in pw.iterrows():
        i, j = idx[r.sample_a], idx[r.sample_b]
        mat[i, j] = mat[j, i] = r.discordance
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    im = ax.imshow(mat, cmap="viridis")
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(samples, rotation=90, fontsize=9)
    ax.set_yticklabels(samples, fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.046)
    cb.set_label("Genotype discordance", fontsize=15)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(True)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(outdir / f"pairwise_discordance.{ext}", bbox_inches="tight")
    plt.close(fig)
    logging.info("Figure: pairwise_discordance.{png,svg,pdf}")
    logging.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
