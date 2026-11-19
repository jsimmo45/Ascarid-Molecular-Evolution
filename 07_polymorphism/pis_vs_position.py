#!/usr/bin/env python3
"""
pis_vs_position.py — is the depressed synonymous diversity in eliminated genes a
property of the sequence, or of the neighbourhood it sits in?

The problem this addresses
--------------------------
piN/piS is ~1.9x higher in eliminated genes, replicating the dN/dS result. But
piS is also ~3x LOWER there. That opens a competing explanation for the whole
finding: if eliminated regions have reduced local effective population size —
plausible if they are repeat-rich, recombination-poor, and under strong linked
selection — then selection is less efficient there because Ne*s is smaller, with
no change in how constrained those genes actually are. Reduced Ne predicts
elevated piN/piS AND reduced piS, which is exactly the observed pattern.

Both stories are "reduced efficacy of selection", but they are different claims:
one is about the genes, the other about the address.

The decisive test
-----------------
**Retained genes flanking eliminated tracts share the neighbourhood but not the
class.** If linked selection in a low-recombination environment is depressing
piS, those flanking retained genes should be depressed too, grading with
distance from the boundary. If piS depression is intrinsic to eliminated
sequence, they should look like ordinary retained genes right up to the CBR,
and piS should step at the boundary exactly as dN/dS does.

That is the same step-versus-gradient logic as Figure 4A/4D in the Ka/Ks paper,
applied to a different quantity, and it discriminates the two hypotheses
cleanly.

Supporting tests
----------------
* piS versus signed distance to the nearest CBR, binned, both classes.
* Within eliminated genes, piS versus tract size — under linked selection,
  bigger tracts have more interior and should be more depressed.
* Within each class, piS interior versus CBR-proximal.

Note on what piS means here: it is diversity per synonymous site pooled over
genes in a bin, sum(pi_syn)/sum(S_sites), not a mean of per-gene ratios. With
20 haplotypes most genes carry very few synonymous SNPs, so per-gene piS is
mostly zeros and a mean of it is meaningless.

Intervals
---------
Every interval is a 95% percentile interval from a gene bootstrap: genes are
resampled with replacement --n-boot times and the pooled statistic recomputed.
Results quotes two contrasts, and both resample each group independently:

* near-minus-far piN/piS: retained genes within --proximal-kb of a CBR minus
  retained genes beyond it. "P(boot > 0)" is the fraction of resamples above 0.
* eliminated / near-retained piN/piS: eliminated genes against the CBR-proximal
  retained genes, i.e. matched on neighbourhood. "P(boot > 1)" is the fraction
  above 1.

Each group draws from its own random stream, spawned from --seed, so a group's
resamples do not depend on which other groups were run before it.

Outputs (in --outdir)
---------------------
    per_gene_position.tsv     genes analysed, with signed CBR distance and
                              tract size
    pis_position_summary.tsv  one row per group or test; 'proximity' rows carry
                              piS and piN/piS with intervals, 'contrast' rows
                              the two contrasts above
    pis_distance_profile.tsv  binned piS profile (TEST 2)
    pis_vs_cbr_distance.{png,svg,pdf}

Example
-------
    python pis_vs_position.py \\
        --per-gene results/pi_filtered_miss4/per_gene.tsv \\
        --kaks ../kaks_regions/kaks_annotated.tsv \\
        --gene-callable targets/gene_callable.tsv \\
        --eliminated-bed $DB_DIR/ascaris/v50/AG.v50.eliminated_strict.bed \\
        --outdir results/pis_position
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
    from paper_colors import RETAINED, ELIMINATED       # type: ignore
except Exception:
    RETAINED, ELIMINATED = "#2C6FBB", "#C0392B"

plt.rcParams.update({
    # Arial where installed, DejaVu Sans (bundled with matplotlib) otherwise --
    # the same list paper_colors uses. A bare "Arial" made matplotlib log a
    # findfont warning for every text element on nodes without Arial: ~600
    # lines per run, burying the results.
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelsize": 20, "xtick.labelsize": 16, "ytick.labelsize": 16,
    "legend.fontsize": 15, "figure.dpi": 150,
})


def pooled_pis(df: pd.DataFrame) -> float:
    S = df.S_sites.sum()
    return df.pi_syn.sum() / S if S > 0 else np.nan


def pooled_pin_pis(df: pd.DataFrame) -> float:
    S, N = df.S_sites.sum(), df.N_sites.sum()
    if S <= 0 or N <= 0:
        return np.nan
    piS = df.pi_syn.sum() / S
    return (df.pi_non.sum() / N) / piS if piS > 0 else np.nan


def boot_ci(df: pd.DataFrame, fn, n_boot: int, seed: int = 0):
    rng = np.random.default_rng(seed)
    idx = np.arange(len(df))
    vals = []
    for _ in range(n_boot):
        v = fn(df.iloc[rng.choice(idx, len(idx), replace=True)])
        if v == v:
            vals.append(v)
    if not vals:
        return (np.nan, np.nan)
    return tuple(np.percentile(vals, [2.5, 97.5]))


def boot_pin_pis(df: pd.DataFrame, n_boot: int, seed) -> np.ndarray:
    """Gene-bootstrap replicates of pooled piN/piS for one group.

    Same statistic as pooled_pin_pis(), computed on arrays rather than on
    DataFrame slices, which matters for the 6,649-gene group at 2,000
    resamples. Sums skip missing values, as pandas does there. A resample with no sites or no synonymous
    diversity has no defined ratio and stays NaN; the caller reports how many
    replicates were defined. ``seed`` is anything np.random.default_rng takes.
    """
    pi_n = df.pi_non.to_numpy(dtype=float)
    pi_s = df.pi_syn.to_numpy(dtype=float)
    s_sites = df.S_sites.to_numpy(dtype=float)
    n_sites = df.N_sites.to_numpy(dtype=float)
    n = len(df)
    reps = np.full(n_boot, np.nan)
    if n == 0:
        return reps
    rng = np.random.default_rng(seed)
    for b in range(n_boot):
        i = rng.integers(0, n, n)
        S, N = np.nansum(s_sites[i]), np.nansum(n_sites[i])
        if S <= 0 or N <= 0:
            continue
        piS = np.nansum(pi_s[i]) / S
        if piS > 0:
            reps[b] = (np.nansum(pi_n[i]) / N) / piS
    return reps


def pct_ci(reps: np.ndarray) -> tuple[float, float]:
    """95% percentile interval over the defined replicates."""
    ok = reps[np.isfinite(reps)]
    if ok.size == 0:
        return (np.nan, np.nan)
    lo, hi = np.percentile(ok, [2.5, 97.5])
    return (float(lo), float(hi))


def summarise_contrast(estimate: float, reps: np.ndarray, null: float) -> dict:
    """Point estimate, 95% interval, fraction of replicates above ``null``,
    and the number of defined replicates."""
    ok = reps[np.isfinite(reps)]
    lo, hi = pct_ci(reps)
    return {"estimate": estimate, "ci_lo": lo, "ci_hi": hi, "null": null,
            "frac_above_null": float(np.mean(ok > null)) if ok.size else np.nan,
            "n_boot_valid": int(ok.size)}


def read_bed_tracts(path: Path) -> dict[str, list[tuple[int, int]]]:
    raw: dict[str, list[tuple[int, int]]] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
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


def tract_size(tracts: dict, chrom: str, pos: int) -> float:
    for s, e in tracts.get(chrom, []):
        if s <= pos < e:
            return e - s
    return np.nan


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Test whether depressed piS is positional (linked selection) "
                    "or intrinsic to eliminated sequence.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--per-gene", required=True, type=Path,
                   help="per_gene.tsv from pi_by_region.py")
    p.add_argument("--kaks", required=True, type=Path,
                   help="kaks_annotated.tsv — supplies min_cbr_dist1")
    p.add_argument("--gene-callable", required=True, type=Path,
                   help="gene_callable.tsv — supplies coordinates")
    p.add_argument("--eliminated-bed", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--proximal-kb", type=float, default=150,
                   help="Distance defining CBR-proximal, in kb. 150 matches the "
                        "window used for the dN/dS step in Figure 4A.")
    p.add_argument("--bins-kb", default="0,50,100,150,250,400,700,1000",
                   help="Comma-separated bin edges in kb for the distance profile")
    p.add_argument("--min-s-sites", type=float, default=2000,
                   help="Minimum pooled synonymous sites for a bin to be drawn")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0,
                   help="Base seed for every bootstrap. 0 reproduces the piS "
                        "intervals of earlier runs exactly.")
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

    df = pd.read_csv(args.per_gene.expanduser(), sep="\t")
    gc = pd.read_csv(args.gene_callable.expanduser(),
                     sep="\t")[["gene", "chrom", "start", "end"]]
    k = pd.read_csv(args.kaks.expanduser(), sep="\t")[
        ["gene1", "min_cbr_dist1"]].rename(columns={"gene1": "gene"})

    df = df.merge(gc, on="gene", how="left").merge(k, on="gene", how="left")
    df = df[df.region.isin(["retained", "eliminated"])].copy()
    df = df[df.min_cbr_dist1.notna()].copy()
    logging.info("Genes with a CBR distance: %d (%s)", len(df),
                 df.region.value_counts().to_dict())

    # Signed distance: eliminated positive/right, matching the paper's Figure 4
    # convention after the 2026-08-17 sign flip.
    df["signed_kb"] = np.where(df.region == "eliminated",
                               df.min_cbr_dist1 / 1e3,
                               -df.min_cbr_dist1 / 1e3)

    tracts = read_bed_tracts(args.eliminated_bed.expanduser())
    mid = ((df.start + df.end) // 2).astype(int)
    df["tract_kb"] = [tract_size(tracts, c, m) / 1e3
                      for c, m in zip(df.chrom, mid)]

    df.to_csv(outdir / "per_gene_position.tsv", sep="\t", index=False)

    prox = args.proximal_kb
    results = []

    # --- THE decisive test -------------------------------------------------
    logging.info("=" * 70)
    logging.info("TEST 1 — do retained genes NEAR a CBR look depressed?")
    logging.info("  (they share the neighbourhood but not the class)")
    ret = df[df.region == "retained"]
    eli = df[df.region == "eliminated"]
    r_prox = ret[ret.min_cbr_dist1 <= prox * 1e3]
    r_far = ret[ret.min_cbr_dist1 > prox * 1e3]

    # One independent stream per group for the piN/piS replicates. The piS
    # intervals keep boot_ci() and its single seed, so they match earlier runs.
    group_seed = dict(zip(("near", "far", "eli"),
                          np.random.SeedSequence(args.seed).spawn(3)))
    ratio_reps: dict[str, np.ndarray] = {}
    for key, label, sub in (("near", "retained, <= %.0f kb" % prox, r_prox),
                            ("far", "retained, > %.0f kb" % prox, r_far),
                            ("eli", "eliminated (all)", eli)):
        if sub.empty:
            continue
        v = pooled_pis(sub)
        lo, hi = boot_ci(sub, pooled_pis, args.n_boot, seed=args.seed)
        rr = pooled_pin_pis(sub)
        ratio_reps[key] = boot_pin_pis(sub, args.n_boot, group_seed[key])
        rr_lo, rr_hi = pct_ci(ratio_reps[key])
        logging.info("  %-24s n=%5d  piS=%.5f [%.5f, %.5f]  "
                     "piN/piS=%.4f [%.4f, %.4f]",
                     label, len(sub), v, lo, hi, rr, rr_lo, rr_hi)
        results.append({"test": "proximity", "group": label, "genes": len(sub),
                        "piS": v, "piS_lo": lo, "piS_hi": hi, "piN_piS": rr,
                        "piN_piS_lo": rr_lo, "piN_piS_hi": rr_hi,
                        "S_sites": sub.S_sites.sum()})

    if not r_prox.empty and not r_far.empty:
        ratio = pooled_pis(r_prox) / pooled_pis(r_far)
        logging.info("")
        logging.info("  retained proximal / retained far piS = %.3f", ratio)
        logging.info("  eliminated / retained far piS        = %.3f",
                     pooled_pis(eli) / pooled_pis(r_far))
        logging.info("")
        if ratio > 0.85:
            logging.info("  READ: retained genes beside a CBR have essentially "
                         "normal piS. The depression is INTRINSIC to eliminated "
                         "sequence, not a property of the neighbourhood — piS "
                         "steps at the boundary, as dN/dS does. Linked selection "
                         "in a shared low-recombination environment is not a "
                         "sufficient explanation.")
        elif ratio < 0.7:
            logging.info("  READ: retained genes beside a CBR are ALSO "
                         "depressed. The neighbourhood matters, so reduced local "
                         "Ne is a live explanation for the elevated piN/piS and "
                         "must be addressed directly — the elevated ratio may "
                         "reflect where these genes sit rather than what they do.")
        else:
            logging.info("  READ: intermediate. Report the number and treat both "
                         "explanations as partially in play.")

    # --- the two contrasts quoted in Results -------------------------------
    # Groups are resampled independently, so each replicate pairs one
    # resample of one group with one of the other.
    contrasts = []
    if {"near", "far", "eli"} <= set(ratio_reps):
        b_near, b_far, b_eli = (ratio_reps["near"], ratio_reps["far"],
                                ratio_reps["eli"])
        est_near, est_far = pooled_pin_pis(r_prox), pooled_pin_pis(r_far)
        est_eli = pooled_pin_pis(eli)

        diff = summarise_contrast(est_near - est_far, b_near - b_far, 0.0)
        ok = np.isfinite(b_eli) & np.isfinite(b_near) & (b_near > 0)
        fold_reps = np.full(args.n_boot, np.nan)
        fold_reps[ok] = b_eli[ok] / b_near[ok]
        fold = summarise_contrast(
            est_eli / est_near if est_near > 0 else np.nan, fold_reps, 1.0)

        logging.info("")
        logging.info("  Contrasts: gene bootstrap, %d resamples per group, "
                     "groups resampled independently, seed %d", args.n_boot,
                     args.seed)
        logging.info("  piN/piS near - far retained   %+.4f [%+.4f, %+.4f]  "
                     "P(boot > 0) = %.2f  (%d/%d defined)",
                     diff["estimate"], diff["ci_lo"], diff["ci_hi"],
                     diff["frac_above_null"], diff["n_boot_valid"], args.n_boot)
        logging.info("  piN/piS eliminated / near     %.3f [%.3f, %.3f]  "
                     "P(boot > 1) = %.3f  (%d/%d defined)",
                     fold["estimate"], fold["ci_lo"], fold["ci_hi"],
                     fold["frac_above_null"], fold["n_boot_valid"], args.n_boot)
        spans0 = diff["ci_lo"] <= 0 <= diff["ci_hi"]
        logging.info("  READ: the near-minus-far interval %s 0, so piN/piS "
                     "beside a CBR is %s that of interior retained genes.",
                     "includes" if spans0 else "EXCLUDES",
                     "indistinguishable from" if spans0 else "different from")
        contrasts = [
            {"test": "contrast",
             "group": f"retained <= {prox:.0f} kb minus retained > {prox:.0f} kb",
             "quantity": "piN/piS difference", **diff,
             "genes_a": len(r_prox), "genes_b": len(r_far)},
            {"test": "contrast",
             "group": f"eliminated / retained <= {prox:.0f} kb",
             "quantity": "piN/piS fold", **fold,
             "genes_a": len(eli), "genes_b": len(r_prox)},
        ]
    else:
        logging.warning("  Contrasts skipped: a group is empty (%s present)",
                        ", ".join(sorted(ratio_reps)) or "none")

    # --- distance profile ---------------------------------------------------
    logging.info("")
    logging.info("TEST 2 — piS versus distance to the nearest CBR")
    edges = [float(x) for x in args.bins_kb.split(",")]
    rows = []
    for region, sub_all in (("retained", ret), ("eliminated", eli)):
        for lo_kb, hi_kb in zip(edges[:-1], edges[1:]):
            sub = sub_all[(sub_all.min_cbr_dist1 / 1e3 >= lo_kb) &
                          (sub_all.min_cbr_dist1 / 1e3 < hi_kb)]
            if sub.empty or sub.S_sites.sum() < args.min_s_sites:
                continue
            v = pooled_pis(sub)
            lo, hi = boot_ci(sub, pooled_pis, args.n_boot, seed=args.seed)
            centre = (lo_kb + hi_kb) / 2
            rows.append({"region": region, "bin_lo_kb": lo_kb,
                         "bin_hi_kb": hi_kb,
                         "centre_kb": centre if region == "eliminated" else -centre,
                         "genes": len(sub), "S_sites": sub.S_sites.sum(),
                         "piS": v, "piS_lo": lo, "piS_hi": hi,
                         "piN_piS": pooled_pin_pis(sub)})
            logging.info("  %-11s %5.0f-%5.0f kb  n=%5d  piS=%.5f [%.5f, %.5f]",
                         region, lo_kb, hi_kb, len(sub), v, lo, hi)
    profile = pd.DataFrame(rows)
    profile.to_csv(outdir / "pis_distance_profile.tsv", sep="\t", index=False)

    # Correlation within each class: a gradient supports linked selection,
    # flatness supports an intrinsic property.
    logging.info("")
    for region, sub in (("retained", ret), ("eliminated", eli)):
        s = sub[sub.S_sites > 0]
        if len(s) < 20:
            continue
        rho, p = spearmanr(s.min_cbr_dist1, s.pi_syn / s.S_sites)
        logging.info("  Spearman(piS, distance) within %-11s rho=%+.3f p=%.3g "
                     "(n=%d)", region, rho, p, len(s))
        results.append({"test": "spearman_distance", "group": region,
                        "genes": len(s), "rho": rho, "p": p})

    # --- tract size ---------------------------------------------------------
    logging.info("")
    logging.info("TEST 3 — within eliminated genes, piS versus tract size")
    e = eli[eli.tract_kb.notna() & (eli.S_sites > 0)]
    if len(e) >= 20:
        rho, p = spearmanr(e.tract_kb, e.pi_syn / e.S_sites)
        logging.info("  Spearman(piS, tract size) rho=%+.3f p=%.3g (n=%d)",
                     rho, p, len(e))
        logging.info("  Negative rho = bigger tracts more depressed = the "
                     "linked-selection prediction.")
        results.append({"test": "spearman_tract", "group": "eliminated",
                        "genes": len(e), "rho": rho, "p": p})
        med = e.tract_kb.median()
        for label, sub in ((f"tract <= {med:.0f} kb", e[e.tract_kb <= med]),
                           (f"tract >  {med:.0f} kb", e[e.tract_kb > med])):
            v = pooled_pis(sub)
            lo, hi = boot_ci(sub, pooled_pis, args.n_boot, seed=args.seed)
            logging.info("  %-18s n=%4d  piS=%.5f [%.5f, %.5f]", label,
                         len(sub), v, lo, hi)
            results.append({"test": "tract_split", "group": label,
                            "genes": len(sub), "piS": v, "piS_lo": lo,
                            "piS_hi": hi})
    else:
        logging.warning("  Too few eliminated genes with a tract assignment")

    # Contrast rows last, so the columns of the earlier rows keep their order.
    # Count columns are nullable integers: the contrast rows leave 'genes'
    # empty, which would otherwise turn every count into a float (233.0).
    summary = pd.DataFrame(results + contrasts)
    for col in ("genes", "n_boot_valid", "genes_a", "genes_b"):
        if col in summary:
            summary[col] = summary[col].astype("Int64")
    summary.to_csv(outdir / "pis_position_summary.tsv", sep="\t", index=False)

    # --- figure -------------------------------------------------------------
    if not profile.empty:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        for region, colour in (("retained", RETAINED),
                               ("eliminated", ELIMINATED)):
            sub = profile[profile.region == region].sort_values("centre_kb")
            if sub.empty:
                continue
            ax.errorbar(sub.centre_kb, sub.piS,
                        yerr=[sub.piS - sub.piS_lo, sub.piS_hi - sub.piS],
                        fmt="o-", color=colour, ecolor=colour, capsize=4,
                        linewidth=2, markersize=8,
                        label=f"{region.capitalize()} (n = {sub.genes.sum():,})")
        ax.axvline(0, ls="--", lw=1.5, color="grey")
        ax.set_xlabel("Distance from nearest CBR (kb)")
        ax.set_ylabel(r"$\pi_S$")
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False)
        for ext in ("png", "svg", "pdf"):
            fig.savefig(outdir / f"pis_vs_cbr_distance.{ext}",
                        bbox_inches="tight")
        plt.close(fig)
        logging.info("Figure: pis_vs_cbr_distance.{png,svg,pdf}")

    logging.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
