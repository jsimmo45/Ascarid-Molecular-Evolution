#!/usr/bin/env python3
"""
sweep_min_s4.py — how small can the pi4-matched windows get before the match
stops meaning anything?

The problem
-----------
pi4_matched_control.py at --min-s4 5000 leaves only **8 retained windows** in
the eliminated pi4 range. That is the weakest number in the figure and the
first thing a referee will poke at. Lowering --min-s4 makes more, smaller
windows, which spreads retained pi4 wider and puts more retained windows into
the eliminated range — so the window count can be improved simply by asking
for it.

Which is exactly why a naive sweep gives the wrong answer.

Smaller windows estimate pi4 from fewer sites, and pi4 is the DENOMINATOR of
pi0/pi4. A window whose pi4 is low by sampling accident gets an inflated
pi0/pi4 automatically, and it is also, for the same accident, more likely to be
selected into the "matched" low-pi4 range. That is the regression-to-the-mean
trap pi4_matched_control.py's whole window design exists to avoid, reappearing
one level up: not ranking GENES on their own pi4, but selecting WINDOWS on
theirs. Push --min-s4 low enough and the matched comparison becomes a
measurement of sampling noise that happens to point the right way.

So "more windows" is not the objective. The objective is the smallest --min-s4
at which the matched comparison is still estimating what it claims to.

The guard: split-half validation
--------------------------------
Within each window, the constituent genes are split into two halves (balanced
on four-fold-site content). pi4 is estimated from half A and pi0/pi4 from
half B. The two estimates share no SNPs, so:

* **Selection on pi4_A is independent of the outcome pi0/pi4_B.** A window that
  looks low-pi4 by accident in half A carries no matching accident in half B,
  so regression to the mean cannot transfer between them.
* The naive statistic uses the same SNPs in both, and therefore carries an
  automatic negative correlation between pi4 and pi0/pi4 from sampling noise
  alone, whether or not any real relationship exists.

Comparing the two is the diagnostic. Where windows are large enough, the naive
and split-half estimates agree. Where they are too small, the naive one stays
confidently negative while the split-half one collapses toward zero — the gap
between them IS the artifact, measured rather than argued about.

Both the TEST 2 trend and the matched fold get a split-half twin here, because
both are exposed to the same trap.

A second, separate artifact the sweep also exposes
--------------------------------------------------
The matched range is [min, max] of eliminated pi4. Smaller windows make
eliminated pi4 noisier, which drags that minimum down and the maximum up, so
the retained comparator quietly absorbs windows sitting at higher pi4 — and
therefore lower pi0/pi4 — which inflates the fold. **The split-half twin does
not protect against this**, since it corrects selection on the denominator, not
the width of the range that selection runs over, and both estimators read the
same range.

So the sweep reports fold drift separately. On a synthetic set built with a
known 1.75x class offset on top of a known Ne gradient, the matched fold read
1.90 at --min-s4 1000 and 1.70 at 15000: shrinking windows inflates it. If the
real data drifts the same way, the large-window value is the conservative one
and the extra windows bought by shrinking come with a biased fold attached.

Outputs
-------
sweep_min_s4.tsv, one row per --min-s4 tried, and a two-panel figure
(supplementary fig. S15; file stem S18_min_s4_sweep): matched fold against
window size, and the naive-versus-split-half divergence that says where to
stop. provenance.json records the run and what the figure drew.

Figure fixes, 2026-10-12 (figure only; the sweep is unchanged)
--------------------------------------------------------------
* Axis titles name the quantity, in plain text: "Minimum four-fold sites per
  window" (was the flag, "--min-s4"), "Eliminated / retained π₀/π₄",
  "Spearman ρ(π₄, π₀/π₄)" (was mathtext).
* The matched fold is drawn in ink (was region red; a fold is a metric).
* The split-half twin is drawn at every setting that has one, not only where
  the naive fold exists: the 10,000-site twin the text quotes was missing.
* Two count rows: naive matched windows (ink) and split-half mean (grey).
* Both panels share an x-range and ticks taken from the settings drawn; a
  setting with nothing to draw is logged, not given an empty tick.
* --replot-from redraws the figure from an earlier sweep_min_s4.tsv in
  seconds, without inputs and without rerunning the sweep.
* The log ends with legend notes for every mark.

Split counts, 2026-10-13 (two new columns; every existing value unchanged)
--------------------------------------------------------------------------
The split-half fold averages only the splits whose half A selects at least 5
retained windows, and the split-half rho only those with at least 10 usable
windows, so the number of splits behind each mean varies with window size. The
table now records them, splithalf_fold_splits and splithalf_rho_splits, and the
legend notes give the fold's count per setting. A table from an earlier run
lacks both columns; --replot-from then says so instead of guessing.

Example
-------
    python sweep_min_s4.py \\
        --per-gene results/pi_clustfilt/per_gene.tsv \\
        --gene-callable targets/gene_callable.tsv \\
        --values 1500,2000,3000,5000,7500,10000,15000 \\
        --outdir results/pi4_sweep

    # figure only, from an existing run:
    python sweep_min_s4.py --outdir results \\
        --replot-from results/pi4_sweep_20260822_pi_filtered_miss4 \\
        --tag pi_filtered_miss4_replot
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from pi4_matched_control import boot, finish, make_windows, pooled_ratio

# style_boot locates lib/paper_colors.py and warns if it cannot, rather than
# silently substituting a different palette.
from style_boot import paper_colors

if paper_colors is None:                # pragma: no cover
    raise SystemExit(
        "sweep_min_s4.py builds a journal supplement and needs the full house "
        "style, not just the colours. Set PAPER_COLORS_DIR to the directory "
        "holding paper_colors.py.")

apply_style = paper_colors.apply_style
audit = paper_colors.audit
figure = paper_colors.figure
headroom = paper_colors.headroom
make_outdir = paper_colors.make_outdir
panel_letter = paper_colors.panel_letter
pc_save = paper_colors.save
FONT = paper_colors.FONT
LINEWIDTH = paper_colors.LINEWIDTH
INK = paper_colors.NEUTRAL["ink"]
INK_LIGHT = paper_colors.NEUTRAL["ink_light"]
RULE = paper_colors.NEUTRAL["rule"]

#: Supplementary fig. S15 (the window sweep; file stem S18_min_s4_sweep kept
#: from an earlier numbering): two panels side by side at double-column width.
FIG_W = "double"
FIG_H_MM = 68.0

REGIONS = ("retained", "eliminated")


def _compact(v: float) -> str:
    """1500 -> '1.5k', 10000 -> '10k'. Shorter than 10^3 notation and exact."""
    if v >= 1000:
        t = v / 1000.0
        return f"{t:.0f}k" if abs(t - round(t)) < 1e-9 else f"{t:g}k"
    return f"{v:g}"


def assign_windows(df: pd.DataFrame, min_s4: float) -> pd.DataFrame:
    """
    Same contiguous grouping as pi4_matched_control.make_windows, but returning
    a per-GENE window label instead of an aggregate.

    Gene identity has to survive for the split-half check to be possible at
    all. The grouping rule is replicated exactly — including the trailing-window
    merge and the drop below min_s4/2 — and `--self-test` asserts that
    aggregating by these labels reproduces make_windows() column for column, so
    the two cannot quietly diverge.

    Genes in a dropped trailing window get window id -1.
    """
    out = []
    for chrom, sub in df.sort_values(["chrom", "start"]).groupby("chrom",
                                                                sort=True):
        sub = sub.copy()
        ids = np.full(len(sub), -1, dtype=np.int64)
        wid = 0
        acc: list[int] = []
        s4 = 0.0
        last_closed = -1
        for i, (_, g) in enumerate(zip(range(len(sub)), sub.itertuples())):
            acc.append(i)
            s4 += getattr(g, "S4_sites")
            if s4 >= min_s4:
                for j in acc:
                    ids[j] = wid
                last_closed = wid
                wid += 1
                acc, s4 = [], 0.0
        if acc:
            if last_closed >= 0:
                # trailing genes merge into the previous window on this contig
                for j in acc:
                    ids[j] = last_closed
            elif s4 >= min_s4 / 2:
                for j in acc:
                    ids[j] = wid
                wid += 1
            # else: left at -1, dropped, matching make_windows
        # zero-padded so a lexicographic sort of the label is also positional
        # order; unpadded, "chr01:10" sorts before "chr01:2"
        sub["window"] = [f"{chrom}:{k:06d}" if k >= 0 else "" for k in ids]
        out.append(sub)
    res = pd.concat(out, ignore_index=True)
    return res[res.window != ""].copy()


def agg(d: pd.DataFrame) -> pd.DataFrame:
    """Window-level pi4 and pi0/pi4 from a gene table carrying `window`."""
    w = d.groupby("window", sort=True).agg(
        genes=("gene", "size"), S4_sites=("S4_sites", "sum"),
        N0_sites=("N0_sites", "sum"), pi_4_sum=("pi_4", "sum"),
        pi_0_sum=("pi_0", "sum")).reset_index()
    w = w[(w.S4_sites > 0) & (w.N0_sites > 0)].copy()
    w["pi4"] = w.pi_4_sum / w.S4_sites
    w["pi0"] = w.pi_0_sum / w.N0_sites
    w = w[w.pi4 > 0].copy()
    w["pi0_pi4"] = w.pi0 / w.pi4
    return w


def split_halves(d: pd.DataFrame, rng) -> pd.DataFrame:
    """
    Label each gene A or B within its window, balancing four-fold content.

    Genes are shuffled and then dealt to whichever half currently holds fewer
    four-fold sites. Random alone would leave the halves unequal in a window of
    a dozen genes, and an unequal split makes pi4_A and pi0/pi4_B differ in
    precision for a reason that has nothing to do with the biology.
    """
    lab = np.empty(len(d), dtype=object)
    pos = {g: i for i, g in enumerate(d.index)}
    for _, sub in d.groupby("window", sort=False):
        idx = np.array(sub.index.to_numpy(), copy=True)
        rng.shuffle(idx)
        tot = {"A": 0.0, "B": 0.0}
        for gi in idx:
            side = "A" if tot["A"] <= tot["B"] else "B"
            tot[side] += float(d.at[gi, "S4_sites"])
            lab[pos[gi]] = side
    d = d.copy()
    d["half"] = lab
    return d


def halves_windows(d: pd.DataFrame, rng):
    """
    (pi4 from half A, pi0/pi4 from half B) per window, on windows where both
    halves carry usable sites.
    """
    d = split_halves(d, rng)
    a = agg(d[d.half == "A"])[["window", "pi4", "S4_sites"]].rename(
        columns={"pi4": "pi4_A", "S4_sites": "S4_A"})
    b = agg(d[d.half == "B"])[["window", "pi0_pi4", "pi_4_sum", "pi_0_sum",
                               "S4_sites", "N0_sites"]].rename(
        columns={"pi0_pi4": "pi0_pi4_B"})
    return a.merge(b, on="window", how="inner")


def evaluate(df: pd.DataFrame, min_s4: float, n_boot: int, n_splits: int,
             seed: int) -> dict:
    """One --min-s4 setting: naive matched fold, and its split-half twin."""
    genes = {r: assign_windows(df[df.region == r], min_s4) for r in REGIONS}
    wins = {r: agg(genes[r]) for r in REGIONS}
    ret, eli = wins["retained"], wins["eliminated"]
    row = {"min_s4": min_s4, "n_ret_windows": len(ret),
           "n_eli_windows": len(eli),
           "median_s4_per_ret_window": float(ret.S4_sites.median())
           if len(ret) else np.nan,
           "median_genes_per_ret_window": float(ret.genes.median())
           if len(ret) else np.nan}
    if len(ret) < 10 or len(eli) < 3:
        logging.warning("min_s4=%g: too few windows (%d retained, %d "
                        "eliminated) — skipped", min_s4, len(ret), len(eli))
        return row

    # --- naive matched comparison, exactly as pi4_matched_control does it ---
    lo, hi = eli.pi4.min(), eli.pi4.max()
    ret_m = ret[(ret.pi4 >= lo) & (ret.pi4 <= hi)]
    row["n_ret_matched"] = len(ret_m)
    row["eli_pi4_lo"], row["eli_pi4_hi"] = lo, hi
    if len(ret_m) >= 5:
        r_val, e_val = pooled_ratio(ret_m), pooled_ratio(eli)
        r_ci, e_ci = boot(ret_m, pooled_ratio, n_boot), boot(eli, pooled_ratio,
                                                             n_boot)
        row.update({"ret_matched_pi0pi4": r_val, "ret_lo": r_ci[0],
                    "ret_hi": r_ci[1], "eli_pi0pi4": e_val, "eli_lo": e_ci[0],
                    "eli_hi": e_ci[1], "matched_fold": e_val / r_val,
                    "fold_lo": e_ci[0] / r_ci[1], "fold_hi": e_ci[1] / r_ci[0],
                    "disjoint": bool(r_ci[1] < e_ci[0])})
    row["naive_spearman"] = float(spearmanr(ret.pi4, ret.pi0_pi4)[0])

    # --- split-half twins ---------------------------------------------------
    # Averaged over several random splits: one split is itself a coin toss, and
    # a single unlucky one would be read as evidence about window size.
    rho_sh, fold_sh, nmatch_sh = [], [], []
    for s in range(n_splits):
        rng = np.random.default_rng(seed + s)
        h_ret = halves_windows(genes["retained"], rng)
        h_eli = halves_windows(genes["eliminated"], rng)
        if len(h_ret) < 10 or len(h_eli) < 3:
            continue
        ok = h_ret[(h_ret.pi4_A > 0) & h_ret.pi0_pi4_B.notna()]
        if len(ok) >= 10:
            rho_sh.append(float(spearmanr(ok.pi4_A, ok.pi0_pi4_B)[0]))
        # matched fold with SELECTION on half A and MEASUREMENT on half B
        lo_a, hi_a = h_eli.pi4_A.min(), h_eli.pi4_A.max()
        sel = h_ret[(h_ret.pi4_A >= lo_a) & (h_ret.pi4_A <= hi_a)]
        if len(sel) >= 5:
            r_b = pooled_ratio(sel)
            e_b = pooled_ratio(h_eli)
            if r_b == r_b and e_b == e_b and r_b > 0:
                fold_sh.append(e_b / r_b)
                nmatch_sh.append(len(sel))
    if rho_sh:
        row["splithalf_spearman"] = float(np.mean(rho_sh))
        row["splithalf_spearman_sd"] = float(np.std(rho_sh))
    if fold_sh:
        row["splithalf_fold"] = float(np.mean(fold_sh))
        row["splithalf_fold_sd"] = float(np.std(fold_sh))
        row["splithalf_n_matched"] = float(np.mean(nmatch_sh))
    if "naive_spearman" in row and "splithalf_spearman" in row:
        row["spearman_gap"] = row["naive_spearman"] - row["splithalf_spearman"]
    # How many of the n_splits splits each split-half mean rests on. Added
    # last, so every earlier column keeps its position in the table.
    row["splithalf_fold_splits"] = len(fold_sh)
    row["splithalf_rho_splits"] = len(rho_sh)
    return row


def self_test(df: pd.DataFrame, min_s4: float) -> None:
    """assign_windows must reproduce make_windows exactly."""
    for region in REGIONS:
        sub = df[df.region == region]
        ref = finish(make_windows(sub, min_s4)).sort_values(
            ["chrom", "start"]).reset_index(drop=True)
        mine = agg(assign_windows(sub, min_s4)).sort_values(
            "window").reset_index(drop=True)
        assert len(ref) == len(mine), (
            f"{region}: {len(ref)} windows from make_windows, {len(mine)} from "
            f"assign_windows at min_s4={min_s4}")
        for col in ("genes", "S4_sites", "N0_sites", "pi_4_sum", "pi_0_sum"):
            a, b = ref[col].to_numpy(float), mine[col].to_numpy(float)
            assert np.allclose(a, b), f"{region}: {col} differs"
    logging.info("Self-test passed: window assignment reproduces "
                 "pi4_matched_control.make_windows at min_s4=%g", min_s4)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sweep --min-s4 for the pi4-matched control, with a "
                    "split-half guard against regression to the mean "
                    "(supplementary fig. S15).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--per-gene", type=Path,
                   help="per_gene.tsv (required unless --replot-from)")
    p.add_argument("--gene-callable", type=Path,
                   help="gene_callable.tsv (required unless --replot-from)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--replot-from", type=Path, default=None,
                   help="Redraw the figure from an earlier run's "
                        "sweep_min_s4.tsv (its folder, or the file itself) "
                        "without rerunning the sweep. The table is copied into "
                        "the new output folder, so the folder stands alone.")
    p.add_argument("--values", default="1500,2000,3000,5000,7500,10000,15000",
                   help="Comma-separated --min-s4 values to try")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--n-splits", type=int, default=25,
                   help="Random split-half replicates averaged per setting")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--tag", default=None,
                   help="Suffix for the dated output directory")
    p.add_argument("--gap-tolerance", type=float, default=0.15,
                   help="Largest acceptable gap between the naive and "
                        "split-half Spearman. Above this the naive statistic "
                        "is being carried by shared sampling noise and the "
                        "window size is too small to trust.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = p.parse_args(argv)
    if args.replot_from is None and (args.per_gene is None
                                     or args.gene_callable is None):
        p.error("--per-gene and --gene-callable are required unless "
                "--replot-from is given")
    return args


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

#: Axis titles. Plain text with Unicode subscripts, as dN/dS: no mathtext.
X_LABEL = "Minimum four-fold sites per window"
A_LABEL = "Eliminated / retained π₀/π₄"
B_LABEL = "Spearman ρ(π₄, π₀/π₄)"

#: Fraction of panel A's height kept free under the lowest drawn value (or
#: the reference line at 1) for the two rows of window counts.
COUNT_BAND = 0.24
#: Heights of the two count rows, in axes fraction: naive above split-half,
#: the order of the legend.
COUNT_ROWS = {"naive": 0.12, "splithalf": 0.035}


def _col(sw: pd.DataFrame, name: str) -> pd.Series:
    """A column of the sweep table, or all-NaN if that run never wrote it."""
    return sw[name] if name in sw else pd.Series(np.nan, index=sw.index)


def text_problems(fig, tol_px: float = 1.0) -> list[str]:
    """Drawn text overlapping other drawn text or running past the canvas
    edge. audit() checks neither, and both save without complaint."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    texts = list(fig.texts)
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks():
                loc = tick.get_loc()
                if loc is not None and lo - 1e-9 <= loc <= hi + 1e-9:
                    texts += [t for t in (tick.label1, tick.label2)
                              if t.get_visible()]
            texts.append(axis.label)
        texts += list(ax.texts)
        leg = ax.get_legend()
        if leg is not None:
            texts += list(leg.get_texts())
    texts = [t for t in texts if t.get_visible() and (t.get_text() or "").strip()]
    boxes = [t.get_window_extent(rend) for t in texts]
    out = []
    for t, b in zip(texts, boxes):
        if b.x0 < -tol_px or b.y0 < -tol_px or b.x1 > W + tol_px \
                or b.y1 > H + tol_px:
            out.append(f"text {t.get_text()[:24]!r} runs past the canvas edge")
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = boxes[i], boxes[j]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) > tol_px
                    and min(a.y1, b.y1) - max(a.y0, b.y0) > tol_px):
                out.append(f"text {texts[i].get_text()[:16]!r} overlaps "
                           f"{texts[j].get_text()[:16]!r}")
    return out


def plot_sweep(sw: pd.DataFrame, outdir: Path, stem: str):
    """
    Supplementary fig. S15 (stem kept from an earlier numbering).

    A  matched fold of pi0/pi4 (eliminated / retained) against window size,
       with its bootstrap interval, and the split-half twin at EVERY setting
       that has one -- the twin needs fewer matched windows than the naive
       comparison, so it can exist where the naive fold does not (it is not
       restricted to the naive fold's rows, which hid the 10,000-site point
       the text quotes). Two rows of counts along the bottom: retained windows
       inside the eliminated pi4 range, naive (ink) over split-half mean
       (grey). Both series in ink tones: a fold change is a metric, not a
       region.
    B  naive and split-half Spearman rho between pi4 and pi0/pi4.

    Both panels share an x-range and ticks, taken from the settings drawn in
    either panel; a setting with nothing to draw (too few windows) is listed
    in the log rather than given an empty tick. Returns (files, problems,
    facts for the log and provenance).
    """
    fig, axes = figure(FIG_W, height_mm=FIG_H_MM, ncols=2)
    facts: dict = {}
    mf = sw[_col(sw, "matched_fold").notna()]
    sh = sw[_col(sw, "splithalf_fold").notna()]
    ax = axes[0]
    a_handles = []                       # legend order: matched, then twin
    if len(mf) or len(sh):
        ax.axhline(1.0, color=RULE, ls=":", lw=LINEWIDTH["rule"])
        if len(mf):
            a_handles.append(ax.errorbar(mf.min_s4, mf.matched_fold,
                        yerr=[np.clip(mf.matched_fold - mf.fold_lo, 0, None),
                              np.clip(mf.fold_hi - mf.matched_fold, 0, None)],
                        fmt="o-", color=INK, ms=3.0, lw=LINEWIDTH["data"],
                        capsize=2.0, elinewidth=LINEWIDTH["thin"],
                        label="Matched fold"))
        if len(sh):
            a_handles += ax.plot(sh.min_s4, sh.splithalf_fold, "s--",
                                 color=INK_LIGHT, ms=3.0, lw=LINEWIDTH["data"],
                                 label="Split-half twin")
        # Limits from every value drawn (see the 2026-08 notes: a point outside
        # the limits reports a window extent past the panel, which audit()
        # then flags against the legend), with 1 always in view, and a band
        # kept free at the bottom for the two count rows.
        drawn = [float(v) for frame, cols in
                 ((mf, ("fold_lo", "fold_hi", "matched_fold")),
                  (sh, ("splithalf_fold",)))
                 for c in cols if c in frame for v in frame[c] if v == v]
        lo_v, hi_v = min(1.0, min(drawn)), max(drawn)
        top = hi_v + 0.06 * (hi_v - lo_v)
        bottom = (lo_v - COUNT_BAND * top) / (1.0 - COUNT_BAND)
        ax.set_ylim(bottom, top)
        for frame, key, row, colour in (
                (mf, "n_ret_matched", "naive", INK),
                (sh, "splithalf_n_matched", "splithalf", INK_LIGHT)):
            for _, r in frame.iterrows():
                v = r.get(key, np.nan)
                if v == v:
                    # (data x, axes-fraction y): a fixed row under the data.
                    ax.annotate(f"{int(round(float(v))):,}",
                                (r.min_s4, COUNT_ROWS[row]),
                                xycoords=ax.get_xaxis_transform(), ha="center",
                                va="bottom", fontsize=FONT["annotation"],
                                color=colour)
        ax.set_xscale("log")
        ax.set_xlabel(X_LABEL)
        ax.set_ylabel(A_LABEL)
        facts["A_matched_fold_at"] = [float(v) for v in mf.min_s4]
        facts["A_splithalf_fold_at"] = [float(v) for v in sh.min_s4]
        facts["A_splithalf_only_at"] = sorted(set(facts["A_splithalf_fold_at"])
                                              - set(facts["A_matched_fold_at"]))
        facts["A_ylim"] = [round(bottom, 4), round(top, 4)]

    e = sw[_col(sw, "naive_spearman").notna()]
    if len(e):
        ax = axes[1]
        ax.axhline(0, color=RULE, ls=":", lw=LINEWIDTH["rule"])
        ax.plot(e.min_s4, e.naive_spearman, "o-", color=INK, ms=3.0,
                lw=LINEWIDTH["data"], label="Naive (shared SNPs)")
        if "splithalf_spearman" in e:
            ax.errorbar(e.min_s4, e.splithalf_spearman,
                        yerr=e.get("splithalf_spearman_sd"),
                        fmt="s--", color=INK_LIGHT, ms=3.0,
                        lw=LINEWIDTH["data"], capsize=2.0,
                        elinewidth=LINEWIDTH["thin"],
                        label="Split-half (independent)")
        ax.set_xscale("log")
        ax.set_xlabel(X_LABEL)
        ax.set_ylabel(B_LABEL)
        rho_drawn = [v for c in ("naive_spearman", "splithalf_spearman")
                     if c in e for v in e[c] if v == v]
        if "splithalf_spearman_sd" in e:
            m_, s_ = e.splithalf_spearman.fillna(0), \
                e.splithalf_spearman_sd.fillna(0)
            rho_drawn += list(m_ + s_) + list(m_ - s_)
        if rho_drawn:
            ax.set_ylim(min(rho_drawn) - 0.05, max(0.0, max(rho_drawn)) + 0.05)
        headroom(ax, 0.10)

    # One x-range and one set of ticks for both panels, from the settings
    # drawn in either, so a vertical line through both meets the same setting.
    drawn_x = sorted({float(v) for v in
                      list(facts.get("A_matched_fold_at", []))
                      + list(facts.get("A_splithalf_fold_at", []))
                      + list(e.min_s4)})
    facts["not_drawn"] = sorted(float(v) for v in sw.min_s4.dropna()
                                if float(v) not in set(drawn_x))
    if drawn_x:
        for a in axes:
            a.set_xlim(drawn_x[0] / 1.25, drawn_x[-1] * 1.25)
            # Ticks at the settings swept, not at powers of ten: the default
            # log locator labels minor decades too, and at 88 mm "3 x 10^3"
            # and "4 x 10^3" overprint each other.
            a.set_xticks(drawn_x)
            a.set_xticklabels([_compact(v) for v in drawn_x])
            a.xaxis.set_minor_formatter(NullFormatter())
            a.tick_params(axis="x", which="minor", length=0)

    # Legends ABOVE each panel: an overlap is impossible by construction.
    for a, handles in zip(axes, (a_handles, None)):
        if a.get_legend_handles_labels()[0]:
            extra = {"handles": handles} if handles else {}
            a.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2,
                     frameon=False, fontsize=FONT["legend"], handlelength=1.4,
                     handletextpad=0.5, columnspacing=1.2, borderaxespad=0.0,
                     **extra)
    for a, letter in zip(axes, "AB"):
        panel_letter(a, letter, dy=1.20)
    problems = list(audit(fig)) + text_problems(fig)
    if problems:
        logging.warning("S15 breaks %d rule(s): %s", len(problems),
                        "; ".join(problems))
    written = pc_save(fig, str(outdir), stem)
    plt.close(fig)
    return written, problems, facts


def cell(r, key: str, spec: str = ".3f") -> str:
    """Format one result field, or '-' when that setting produced none."""
    v = r.get(key, np.nan)
    if v is None or (isinstance(v, float) and v != v):
        return "-"
    try:
        return format(int(v) if spec == "d" else float(v), spec)
    except (TypeError, ValueError):
        return "-"


def split_count_notes(sw: pd.DataFrame, drawn_at, n_splits) -> list[str]:
    """
    Legend-note lines giving, per setting whose split-half fold is drawn, how
    many splits its mean rests on, and the same for rho where any setting used
    fewer than all of them. A table written before 2026-10-13 has no count
    columns, and the lines say so rather than implying every split was used.
    """
    of = f" of {int(n_splits)}" if n_splits else ""
    if "splithalf_fold_splits" not in sw:
        return ["Splits behind each split-half mean: not recorded in this "
                "table (written before 2026-10-13); rerun the sweep to record "
                "them."]
    drawn = {float(v) for v in drawn_at}
    lines = []
    sub = sw[sw.min_s4.astype(float).isin(drawn)]
    if len(sub):
        lines.append("Splits behind each split-half fold%s: %s." % (
            of, "; ".join(f"{_compact(float(r.min_s4))} "
                          f"{cell(r, 'splithalf_fold_splits', 'd')}"
                          for _, r in sub.iterrows())))
    if "splithalf_rho_splits" in sw and "splithalf_spearman" in sw \
            and n_splits and sw.splithalf_spearman.notna().any():
        shown = sw[sw.splithalf_spearman.notna()]
        rho = shown[shown.splithalf_rho_splits < int(n_splits)]
        if len(rho):
            lines.append("Panel B, split-half rho from fewer than all %d "
                         "splits (at least 10 usable windows needed): %s." % (
                             int(n_splits),
                             "; ".join(f"{_compact(float(r.min_s4))} "
                                       f"{cell(r, 'splithalf_rho_splits', 'd')}"
                                       for _, r in rho.iterrows())))
        else:
            lines.append("Panel B, split-half rho: every setting drawn used "
                         "all %d splits." % int(n_splits))
    return lines


def report(sw: pd.DataFrame, gap_tolerance: float) -> None:
    """The results table, the verdict and the fold drift, to the log."""
    logging.info("=" * 88)
    logging.info("%9s %7s %7s %8s %9s %9s %8s %8s", "min_s4", "ret_w",
                 "eli_w", "matched", "fold", "SH fold", "rho", "SH rho")
    for _, r in sw.iterrows():
        logging.info("%9g %7s %7s %8s %9s %9s %8s %8s", r.min_s4,
                     cell(r, "n_ret_windows", "d"),
                     cell(r, "n_eli_windows", "d"),
                     cell(r, "n_ret_matched", "d"),
                     cell(r, "matched_fold"), cell(r, "splithalf_fold"),
                     cell(r, "naive_spearman", "+.3f"),
                     cell(r, "splithalf_spearman", "+.3f"))
    if _col(sw, "matched_fold").isna().all():
        logging.warning("")
        logging.warning("No setting produced >=5 retained windows inside the "
                        "eliminated pi4 range, so there is no matched "
                        "comparison to report at any window size. That is a "
                        "statement about the data, not a failure: the two "
                        "classes barely share a pi4 range. Fall back on the "
                        "fitted-trend comparison in pi4_matched_control.py "
                        "(TEST 2) and say so explicitly.")

    ok = sw[sw.spearman_gap.abs() <= gap_tolerance] \
        if "spearman_gap" in sw else sw.iloc[0:0]
    logging.info("")
    if ok.empty:
        logging.warning("VERDICT: at EVERY window size tried, the naive and "
                        "split-half trends disagree by more than %.2f. The "
                        "naive statistic is being carried by shared sampling "
                        "noise across this whole range — keep the largest "
                        "window size and do not lower it to buy window count.",
                        gap_tolerance)
    else:
        best = ok.loc[ok.min_s4.idxmin()]
        cur = sw[sw.min_s4 == 5000]
        logging.info("VERDICT: the naive and split-half trends agree (gap "
                     "<= %.2f) down to min_s4 = %g, giving %s matched "
                     "retained window(s) against %s at the current setting "
                     "of 5000.", gap_tolerance, best.min_s4,
                     cell(best, "n_ret_matched", "d"),
                     cell(cur.iloc[0], "n_ret_matched", "d") if len(cur)
                     else "?")
        if best.get("matched_fold", np.nan) == best.get("matched_fold",
                                                        np.nan):
            disjoint = str(best.get("disjoint")).strip().lower() == "true"
            logging.info("  At that setting the matched fold is %s "
                         "[%s, %s]%s; the split-half twin, which selects "
                         "windows on data it does not then measure, gives %s. "
                         "Agreement between those two is what licenses the "
                         "smaller window — not the window count itself.",
                         cell(best, "matched_fold"), cell(best, "fold_lo"),
                         cell(best, "fold_hi"),
                         " (CIs disjoint)" if disjoint else "",
                         cell(best, "splithalf_fold"))
    # --- the second artifact: range widening --------------------------------
    # The matched range is [min, max] of eliminated pi4. Smaller windows make
    # eliminated pi4 noisier, which pushes that min DOWN and that max UP, so
    # the retained comparator absorbs higher-pi4 windows whose pi0/pi4 is
    # lower -- inflating the fold. The split-half twin does NOT protect against
    # this: it corrects selection on the denominator, not the width of the
    # window the selection runs over, and both estimators read the same range.
    # Validated on a synthetic set with a known 1.75x class offset, where the
    # fold ran 1.90 at min_s4 1000 and 1.70 at 15000.
    d2 = sw[_col(sw, "matched_fold").notna()]
    if len(d2) >= 3:
        drift = float(d2.matched_fold.iloc[0] - d2.matched_fold.iloc[-1])
        logging.info("")
        logging.info("FOLD DRIFT across the sweep: %.3f at min_s4=%g down to "
                     "%.3f at min_s4=%g (change %+.3f).",
                     d2.matched_fold.iloc[0], d2.min_s4.iloc[0],
                     d2.matched_fold.iloc[-1], d2.min_s4.iloc[-1], -drift)
        if drift > 0.05:
            logging.warning("The fold RISES as windows shrink. That is the "
                            "range-widening artifact, not extra signal — a "
                            "noisier eliminated pi4 stretches the matched "
                            "range and lets in retained windows that sit at "
                            "higher pi4 and therefore lower pi0/pi4. Quote a "
                            "LARGE-window value: it is the conservative end, "
                            "and any gain in window count bought by shrinking "
                            "comes with an inflated fold attached.")
        elif drift < -0.05:
            logging.info("The fold FALLS as windows shrink, i.e. the current "
                         "large-window setting is if anything generous. Report "
                         "the smaller-window value as the conservative one.")
        else:
            logging.info("The fold is stable across window size (|change| <= "
                         "0.05), so the estimate does not depend on this "
                         "choice and the window count can be taken at the "
                         "setting that passes the split-half check.")
    sh = sw[_col(sw, "splithalf_fold").notna()]
    if len(sh):
        logging.info("SPLIT-HALF FOLD over every setting that has one "
                     "(%s): mean %.3f, SD %.3f (n-1) / %.3f (n), range "
                     "%.3f-%.3f.",
                     ", ".join(_compact(v) for v in sh.min_s4),
                     sh.splithalf_fold.mean(), sh.splithalf_fold.std(ddof=1)
                     if len(sh) > 1 else 0.0, sh.splithalf_fold.std(ddof=0),
                     sh.splithalf_fold.min(), sh.splithalf_fold.max())

    logging.info("")
    logging.info("HOW TO READ THIS: do not simply take the setting with the "
                 "most matched windows. Smaller windows ALWAYS produce more "
                 "overlap, because a noisier pi4 spreads retained windows into "
                 "the eliminated range — the extra windows are the artifact, "
                 "not extra evidence. The split-half columns are the check: "
                 "where they track the naive ones, the windows are big enough; "
                 "where they collapse toward zero while the naive stays "
                 "strong, they are not.")


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    for noisy in ("fontTools", "matplotlib.font_manager", "matplotlib"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    apply_style()
    outdir = Path(make_outdir("pi4_sweep",
                              base=str(args.outdir.expanduser().resolve()),
                              tag=args.tag))
    logging.info("Output: %s", outdir)

    prov: dict = {"argv": sys.argv, "gap_tolerance": args.gap_tolerance}
    if args.replot_from is not None:
        # --- figure only, from an earlier run's table ----------------------
        src = args.replot_from.expanduser()
        if src.is_dir():
            src = src / "sweep_min_s4.tsv"
        if not src.is_file():
            logging.error("--replot-from: no sweep_min_s4.tsv at %s", src)
            return 2
        sw = pd.read_csv(src, sep="\t")
        if "min_s4" not in sw:
            logging.error("%s has no min_s4 column; is it a sweep table?", src)
            return 2
        logging.info("Replot from %s (%d settings); the sweep is NOT rerun",
                     src, len(sw))
        if src.resolve() != (outdir / "sweep_min_s4.tsv").resolve():
            shutil.copyfile(src, outdir / "sweep_min_s4.tsv")
        old = src.parent / "provenance.json"
        try:
            old_prov = json.loads(old.read_text()) if old.is_file() else {}
        except (OSError, ValueError):
            old_prov = {}
        prov.update({"replot_from": str(src.resolve()),
                     "replot_from_sha256": hashlib.sha256(
                         src.read_bytes()).hexdigest(),
                     "per_gene": old_prov.get("per_gene"),
                     "values": [float(v) for v in sw.min_s4],
                     "n_boot": old_prov.get("n_boot"),
                     "n_splits": old_prov.get("n_splits")})
    else:
        pg = pd.read_csv(args.per_gene.expanduser(), sep="\t")
        gc = pd.read_csv(args.gene_callable.expanduser(),
                         sep="\t")[["gene", "chrom", "start", "end"]]
        df = pg.merge(gc, on="gene", how="inner")
        df = df[df.region.isin(REGIONS)]
        logging.info("Genes: %s", df.region.value_counts().to_dict())

        values = [float(v) for v in str(args.values).split(",") if v.strip()]
        self_test(df, values[len(values) // 2])

        rows = []
        for v in values:
            logging.info("min_s4 = %g", v)
            rows.append(evaluate(df, v, args.n_boot, args.n_splits, args.seed))
        sw = pd.DataFrame(rows)
        sw.to_csv(outdir / "sweep_min_s4.tsv", sep="\t", index=False)
        prov.update({"per_gene": str(args.per_gene), "values": values,
                     "n_boot": args.n_boot, "n_splits": args.n_splits})
    prov["rows"] = sw.to_dict("records")
    # Written now and again after the figure, so a failed plot still leaves
    # the run documented.
    (outdir / "provenance.json").write_text(json.dumps(prov, indent=2,
                                                       default=str))

    report(sw, args.gap_tolerance)

    # --- figure -------------------------------------------------------------
    written, problems, facts = plot_sweep(sw, outdir, "S18_min_s4_sweep")
    logging.info("Figure: %s", ", ".join(Path(w).name for w in written))
    logging.info("Export from the SVG — Arial is absent on Isaac nodes.")
    nb, ns = prov.get("n_boot"), prov.get("n_splits")
    logging.info("")
    logging.info("LEGEND NOTES (supplementary fig. S15)")
    logging.info("  A: circles = matched fold of pi0/pi4, eliminated / "
                 "retained windows inside the eliminated pi4 range, with "
                 "%s-replicate bootstrap 95%% intervals; squares = split-half "
                 "twin (pi4 from half A selects, pi0/pi4 from half B measures), "
                 "mean over the splits (of %s random splits) whose half A "
                 "selects at least 5 retained windows; dotted line = no "
                 "difference.", nb if nb else "?", ns if ns else "?")
    for line in split_count_notes(sw, facts.get("A_splithalf_fold_at", []),
                                  ns):
        logging.info("     %s", line)
    logging.info("     Numbers along the bottom = retained windows inside the "
                 "eliminated pi4 range: naive (ink, upper row) and split-half "
                 "mean (grey, lower row).")
    if facts.get("A_splithalf_only_at"):
        logging.info("     Split-half twin only (the naive comparison had fewer "
                     "than 5 matched windows) at: %s.",
                     ", ".join(_compact(v) for v in facts["A_splithalf_only_at"]))
    logging.info("  B: Spearman rho between pi4 and pi0/pi4 across retained "
                 "windows; naive = both from the same SNPs; split-half = mean "
                 "+/- SD over the splits; dotted line = 0.")
    if facts.get("not_drawn"):
        logging.info("  Settings in the table with nothing to draw (too few "
                     "windows): %s.",
                     ", ".join(_compact(v) for v in facts["not_drawn"]))
    prov.update({"figure": [Path(w).name for w in written],
                 "figure_facts": facts, "audit": problems,
                 "split_count_notes": split_count_notes(
                     sw, facts.get("A_splithalf_fold_at", []), ns)})
    (outdir / "provenance.json").write_text(json.dumps(prov, indent=2,
                                                       default=str))
    logging.info("=" * 88)
    return 0


if __name__ == "__main__":
    sys.exit(main())
