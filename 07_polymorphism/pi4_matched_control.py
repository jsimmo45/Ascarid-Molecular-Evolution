#!/usr/bin/env python3
"""
pi4_matched_control.py — is elevated pi0/pi4 in eliminated genes explained by
their lower neutral diversity (i.e. by reduced local Ne)?

The question
------------
Eliminated genes show pi0/pi4 ~2x higher than retained genes, AND pi4 ~3x lower.
Since dS is unchanged between Ascaris and Parascaris, the mutation rate is the
same in both classes, so pi = 4*Ne*mu with 3x lower pi4 points at ~3x lower Ne.
Reduced Ne makes selection less efficient (Ne*s smaller) and elevates pi0/pi4
on its own, with no change in how constrained the genes are.

The earlier neighbourhood control (retained genes near a CBR: 41% lower piS,
piN/piS unchanged) argues against this, but 41% is not 300% — it does not span
the range in question. This script does.

The trap, and why the obvious implementation is wrong
-----------------------------------------------------
The tempting approach is to rank individual retained genes by their own pi4 and
compare pi0/pi4 in the low-pi4 ones. **That manufactures the result.** A gene
has ~165 four-fold sites here, so at pi4 ~0.0065 it carries roughly ONE
four-fold SNP. Ranking genes on such an estimate selects, at the bottom, genes
whose pi4 is low *by sampling accident* — and pi4 is the DENOMINATOR of
pi0/pi4, so those genes get an inflated ratio automatically. Regression to the
mean would produce exactly the pattern being tested for.

This script therefore aggregates genes into **positionally contiguous windows**
with a floor on four-fold-site content, so each window's pi4 is estimated from
thousands of sites rather than one SNP, and windows are defined by genomic
position — never by the statistic being tested.

What it produces
----------------
1. Windows of both classes on shared axes: pi0/pi4 versus pi4. If the two
   classes form ONE relationship, reduced Ne explains the difference. If
   eliminated windows sit systematically ABOVE retained windows at matched pi4,
   it does not.
2. A direct matched comparison over whatever pi4 range the two classes share.
3. An OLS fit of pi0/pi4 on log(pi4) across retained windows, with the
   eliminated mean overlaid — so the comparison is stated even if the ranges
   barely overlap.

Example
-------
    python pi4_matched_control.py \\
        --per-gene results/pi_clustfilt/per_gene.tsv \\
        --gene-callable targets/gene_callable.tsv \\
        --outdir results/pi4_matched
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

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


def make_windows(df: pd.DataFrame, min_s4: float) -> pd.DataFrame:
    """
    Group positionally contiguous genes into windows, each holding at least
    ``min_s4`` four-fold sites.

    Windows are defined by POSITION, never by pi4 — see the module docstring.
    A trailing window below the floor is merged into the previous one rather
    than reported with an unreliable estimate.
    """
    rows = []
    for chrom, sub in df.sort_values(["chrom", "start"]).groupby("chrom",
                                                                sort=True):
        acc = []
        s4 = 0.0
        for _, g in sub.iterrows():
            acc.append(g)
            s4 += g.S4_sites
            if s4 >= min_s4:
                rows.append(_summarise(acc, chrom))
                acc, s4 = [], 0.0
        if acc:
            if rows and rows[-1]["chrom"] == chrom:
                merged = rows.pop()
                rows.append(_merge(merged, _summarise(acc, chrom)))
            elif s4 >= min_s4 / 2:
                rows.append(_summarise(acc, chrom))
    return pd.DataFrame(rows)


def _summarise(acc, chrom) -> dict:
    d = pd.DataFrame(acc)
    return {"chrom": chrom, "start": int(d.start.min()), "end": int(d.end.max()),
            "genes": len(d), "S4_sites": d.S4_sites.sum(),
            "N0_sites": d.N0_sites.sum(),
            "pi_4_sum": d.pi_4.sum(), "pi_0_sum": d.pi_0.sum()}


def _merge(a: dict, b: dict) -> dict:
    return {"chrom": a["chrom"], "start": min(a["start"], b["start"]),
            "end": max(a["end"], b["end"]), "genes": a["genes"] + b["genes"],
            "S4_sites": a["S4_sites"] + b["S4_sites"],
            "N0_sites": a["N0_sites"] + b["N0_sites"],
            "pi_4_sum": a["pi_4_sum"] + b["pi_4_sum"],
            "pi_0_sum": a["pi_0_sum"] + b["pi_0_sum"]}


def finish(w: pd.DataFrame) -> pd.DataFrame:
    w = w[(w.S4_sites > 0) & (w.N0_sites > 0)].copy()
    w["pi4"] = w.pi_4_sum / w.S4_sites
    w["pi0"] = w.pi_0_sum / w.N0_sites
    w = w[w.pi4 > 0].copy()
    w["pi0_pi4"] = w.pi0 / w.pi4
    return w


def pooled_ratio(w: pd.DataFrame) -> float:
    S4, N0 = w.S4_sites.sum(), w.N0_sites.sum()
    if S4 <= 0 or N0 <= 0:
        return np.nan
    pi4 = w.pi_4_sum.sum() / S4
    return (w.pi_0_sum.sum() / N0) / pi4 if pi4 > 0 else np.nan


def boot_dist(w: pd.DataFrame, fn, n_boot: int, seed: int = 0) -> np.ndarray:
    """Bootstrap replicates over windows, returned rather than summarised.

    boot() below reduces these to percentiles; the trend test in TEST 2 needs
    the standard error, which cannot be recovered from a percentile interval
    without assuming normality.
    """
    rng = np.random.default_rng(seed)
    i = np.arange(len(w))
    v = [fn(w.iloc[rng.choice(i, len(i), replace=True)]) for _ in range(n_boot)]
    return np.array([x for x in v if x == x])


def boot(w: pd.DataFrame, fn, n_boot: int, seed: int = 0):
    rng = np.random.default_rng(seed)
    i = np.arange(len(w))
    v = [fn(w.iloc[rng.choice(i, len(i), replace=True)]) for _ in range(n_boot)]
    v = np.array([x for x in v if x == x])
    return (np.nan, np.nan) if v.size == 0 else tuple(np.percentile(v, [2.5, 97.5]))


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Test whether elevated pi0/pi4 in eliminated genes is "
                    "explained by their lower pi4 (reduced Ne).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--per-gene", required=True, type=Path,
                   help="per_gene.tsv from pi_by_region.py (needs pi_4, "
                        "S4_sites, pi_0, N0_sites)")
    p.add_argument("--gene-callable", required=True, type=Path,
                   help="gene_callable.tsv — supplies chrom/start/end")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-s4", type=float, default=10000,
                   help="Minimum four-fold sites per window. Higher = less "
                        "noisy pi4, fewer windows. Eliminated genes hold only "
                        "~70k four-fold sites in total, so much above 10000 "
                        "leaves too few eliminated windows to compare.")
    p.add_argument("--n-boot", type=int, default=2000)
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

    pg = pd.read_csv(args.per_gene.expanduser(), sep="\t")
    need = {"pi_4", "S4_sites", "pi_0", "N0_sites", "region", "gene"}
    missing = need - set(pg.columns)
    if missing:
        logging.error("per_gene.tsv lacks %s — re-run pi_by_region.py with the "
                      "version that reports pi0/pi4", ", ".join(sorted(missing)))
        return 2
    gc = pd.read_csv(args.gene_callable.expanduser(),
                     sep="\t")[["gene", "chrom", "start", "end"]]
    df = pg.merge(gc, on="gene", how="inner")
    df = df[df.region.isin(["retained", "eliminated"])]
    logging.info("Genes: %s", df.region.value_counts().to_dict())

    wins = {}
    for region in ("retained", "eliminated"):
        w = finish(make_windows(df[df.region == region], args.min_s4))
        w["region"] = region
        wins[region] = w
        logging.info("%-11s %d window(s), median %d genes and %.0f four-fold "
                     "sites each; pi4 range %.5f-%.5f",
                     region, len(w), int(w.genes.median()),
                     w.S4_sites.median(), w.pi4.min(), w.pi4.max())
    ret, eli = wins["retained"], wins["eliminated"]
    allw = pd.concat([ret, eli], ignore_index=True)
    allw.to_csv(outdir / "windows.tsv", sep="\t", index=False)

    if len(eli) < 3 or len(ret) < 10:
        logging.error("Too few windows to compare; lower --min-s4")
        return 1

    logging.info("=" * 70)

    # --- 1. do the classes form one relationship? --------------------------
    logging.info("TEST 1 — matched-pi4 comparison over the shared range")
    lo, hi = eli.pi4.min(), eli.pi4.max()
    ret_m = ret[(ret.pi4 >= lo) & (ret.pi4 <= hi)]
    logging.info("  eliminated pi4 range: %.5f - %.5f", lo, hi)
    logging.info("  retained windows in that range: %d of %d",
                 len(ret_m), len(ret))

    rows = []
    if len(ret_m) >= 5:
        r_val = pooled_ratio(ret_m)
        r_ci = boot(ret_m, pooled_ratio, args.n_boot)
        e_val = pooled_ratio(eli)
        e_ci = boot(eli, pooled_ratio, args.n_boot)
        logging.info("  retained   (pi4-matched) pi0/pi4 = %.4f [%.4f, %.4f] "
                     "  mean pi4 %.5f", r_val, r_ci[0], r_ci[1],
                     ret_m.pi_4_sum.sum() / ret_m.S4_sites.sum())
        logging.info("  eliminated               pi0/pi4 = %.4f [%.4f, %.4f] "
                     "  mean pi4 %.5f", e_val, e_ci[0], e_ci[1],
                     eli.pi_4_sum.sum() / eli.S4_sites.sum())
        logging.info("  matched fold = %.3f", e_val / r_val)
        rows += [{"test": "matched", "region": "retained", "windows": len(ret_m),
                  "pi0_pi4": r_val, "lo": r_ci[0], "hi": r_ci[1]},
                 {"test": "matched", "region": "eliminated", "windows": len(eli),
                  "pi0_pi4": e_val, "lo": e_ci[0], "hi": e_ci[1]}]
        if r_ci[1] < e_ci[0]:
            logging.info("  READ: at MATCHED neutral diversity the eliminated "
                         "class is still higher, intervals disjoint. Reduced "
                         "Ne does not explain the contrast.")
        elif e_ci[1] < r_ci[0]:
            logging.info("  READ: eliminated is LOWER at matched pi4 — "
                         "unexpected; inspect before interpreting.")
        else:
            logging.info("  READ: intervals overlap at matched pi4. Reduced Ne "
                         "is a sufficient explanation over this range; report "
                         "it as such rather than claiming otherwise.")
    else:
        logging.warning("  Only %d retained window(s) overlap the eliminated "
                        "pi4 range — the direct match is underpowered; rely on "
                        "TEST 2.", len(ret_m))

    # --- 2. the relationship across retained windows ------------------------
    logging.info("")
    logging.info("TEST 2 — pi0/pi4 versus pi4 across retained windows")
    x = np.log(ret.pi4.to_numpy())
    y = ret.pi0_pi4.to_numpy()
    rho, p = spearmanr(ret.pi4, ret.pi0_pi4)
    slope, intercept = np.polyfit(x, y, 1)
    logging.info("  Spearman(pi4, pi0/pi4) = %+.3f, p = %.3g (n = %d windows)",
                 rho, p, len(ret))
    logging.info("  OLS pi0/pi4 = %.4f %+.4f * ln(pi4)", intercept, slope)

    # A typical retained window, for the decomposition below: the trend
    # evaluated at the retained class's own pooled pi4, not the mean of the
    # per-window ratios, so it is on the same footing as the two other terms.
    r_pi4 = ret.pi_4_sum.sum() / ret.S4_sites.sum()
    r_typ = intercept + slope * np.log(r_pi4)

    e_pi4 = eli.pi_4_sum.sum() / eli.S4_sites.sum()
    e_ratio = pooled_ratio(eli)
    pred = intercept + slope * np.log(e_pi4)
    resid = y - (intercept + slope * x)
    sd = resid.std(ddof=2)

    # How far the eliminated class sits above the trend, in units of the
    # scatter of INDIVIDUAL retained windows. Descriptive only -- it says the
    # eliminated class is more extreme than most single retained windows, which
    # is a statement about the spread, not about whether the class mean differs
    # from the trend. Reporting it is useful; testing with it is not, and an
    # earlier version of this script did exactly that and drew the opposite
    # conclusion from TEST 1 on the same data.
    z_scatter = (e_ratio - pred) / sd if sd > 0 else np.nan

    # The test. e_ratio is a POOLED estimate over all eliminated windows, so it
    # must be compared against its own standard error, not against the spread of
    # single windows. Both terms of the difference carry uncertainty:
    #   - the eliminated pooled ratio, from a bootstrap over its windows;
    #   - the fitted prediction at that pi4, s * sqrt(1/n + (x0 - xbar)^2 / Sxx).
    e_boot = boot_dist(eli, pooled_ratio, args.n_boot)
    se_eli = float(np.std(e_boot, ddof=1)) if e_boot.size > 1 else np.nan
    x0 = float(np.log(e_pi4))
    sxx = float(((x - x.mean()) ** 2).sum())
    se_pred = (sd * np.sqrt(1.0 / len(x) + (x0 - x.mean()) ** 2 / sxx)
               if sxx > 0 else np.nan)
    se_diff = float(np.sqrt(se_eli ** 2 + se_pred ** 2))
    z = (e_ratio - pred) / se_diff if se_diff > 0 else np.nan
    p_trend = float(2 * norm.sf(abs(z))) if z == z else np.nan

    logging.info("  At the eliminated mean pi4 (%.5f) the retained "
                 "relationship predicts pi0/pi4 = %.4f (SE %.4f)",
                 e_pi4, pred, se_pred)
    logging.info("  Eliminated observed = %.4f (SE %.4f)", e_ratio, se_eli)
    logging.info("  Difference = %+.4f, %.2f SE above the trend, P = %.3g",
                 e_ratio - pred, z, p_trend)
    logging.info("  (For scale: that is %.2f SD of the scatter of individual "
                 "retained windows -- descriptive, not the test.)", z_scatter)
    rows.append({"test": "trend", "region": "eliminated", "windows": len(eli),
                 "pi0_pi4": e_ratio, "predicted_from_retained": pred,
                 "se_eliminated": se_eli, "se_prediction": se_pred,
                 "z_vs_retained_trend": z, "p_vs_retained_trend": p_trend,
                 "z_vs_window_scatter": z_scatter,
                 "spearman_rho": rho, "spearman_p": p})

    if rho < 0 and p < 0.05:
        logging.info("  Note: pi0/pi4 DOES decline with pi4 among retained "
                     "windows, i.e. lower-Ne retained regions are genuinely "
                     "less efficiently selected. The question is whether that "
                     "trend is enough to account for the eliminated class.")
    if z == z and z > 0 and p_trend < 0.05:
        share = ((pred - r_typ) / (e_ratio - r_typ)
                 if r_typ == r_typ and e_ratio != r_typ else np.nan)
        logging.info("  READ: eliminated sits significantly above the retained "
                     "trend extrapolated to its own pi4 (P = %.3g). Reduced Ne "
                     "does not account for it.", p_trend)
        if share == share:
            logging.info("  Decomposition: of the total elevation above a "
                         "typical retained window, the trend explains %.0f%% "
                         "and %.0f%% remains.", 100 * share, 100 * (1 - share))
    elif z == z and p_trend >= 0.05:
        logging.info("  READ: eliminated is not distinguishable from the "
                     "retained trend at its own pi4 (P = %.3g). Over this "
                     "range reduced Ne is a sufficient explanation and the "
                     "constraint claim must be softened accordingly.", p_trend)
    elif z == z and z < 0:
        logging.info("  READ: eliminated sits BELOW the retained trend — "
                     "unexpected; inspect before interpreting.")

    pd.DataFrame(rows).to_csv(outdir / "pi4_matched_summary.tsv", sep="\t",
                              index=False)

    # --- figure -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    ax.scatter(ret.pi4, ret.pi0_pi4, s=26, color=RETAINED, alpha=0.55,
               edgecolor="none", rasterized=True,
               label=f"Retained windows (n = {len(ret)})")
    ax.scatter(eli.pi4, eli.pi0_pi4, s=70, color=ELIMINATED,
               edgecolor="black", linewidth=0.5,
               label=f"Eliminated windows (n = {len(eli)})")
    xs = np.linspace(ret.pi4.min(), ret.pi4.max(), 100)
    ax.plot(xs, intercept + slope * np.log(xs), color=RETAINED, linewidth=2,
            linestyle="--", label="Retained trend")
    ax.axvspan(lo, hi, color=ELIMINATED, alpha=0.08, zorder=0)
    ax.set_xlabel(r"$\pi_4$  (neutral diversity)")
    ax.set_ylabel(r"$\pi_0/\pi_4$")
    ax.set_xscale("log")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc="upper right")
    for ext in ("png", "svg", "pdf"):
        fig.savefig(outdir / f"pi0pi4_vs_pi4.{ext}", bbox_inches="tight")
    plt.close(fig)
    logging.info("Figure: pi0pi4_vs_pi4.{png,svg,pdf}")
    logging.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
