#!/usr/bin/env python3
"""
figure_polymorphism.py — assemble the polymorphism figure and its supplements.

Main figure, four panels, reading left to right as an argument:

  A  Every estimator agrees      — forest of eliminated/retained fold across
                                   dN/dS, piN/piS, pi0/pi4, and the pi4-matched
                                   contrast, all on one axis with CIs.
  B  It survives the Ne control  — pi0/pi4 for every gene in each class (the
                                   left pair) and for the pi4-matched windows
                                   (the right pair, hatched). All four bars are
                                   the same estimator; see panel_b().
  C  The mechanism               — pi0/pi4 against pi4 across windows, the
                                   retained trend drawn, eliminated overlaid.
  D  It is not adaptation        — DoS by class, MAF-filtered.

Panel A is a forest rather than grouped bars on purpose. dN/dS is quoted as a
median of per-gene ratios (the Ka/Ks paper's convention) while piN/piS is a
pooled ratio of sums; the two are different estimators on different scales and
putting their raw values on a shared axis would invite a false comparison. The
FOLD between classes is the quantity that is comparable, so that is what the
panel shows.

Paper conventions enforced: retained blue and on the left, eliminated red and
on the right; Arial; no top or right spines; n embedded in tick labels;
PNG/SVG/PDF with svg.fonttype='none' so Illustrator keeps the text live.
Export from the SVG — Arial is not installed on Isaac compute nodes.


CHANGES, 2026-10-14 (figure only; every number is computed as before)
--------------------------------------------------------------------
  * Panel B drew the pi4-matched retained windows (0.340) beside the
    eliminated GENES (0.520), so it read as a 1.53-fold matched contrast; the
    matched comparison in panel A and in Results is against the 9 eliminated
    WINDOWS (0.489; 1.44-fold). The eliminated windows now have a bar: the
    outer pair is every gene (panel A's 1.82-fold "no NG weighting" row), the
    inner pair the matched windows (1.44-fold), hatched.
  * The hatch is drawn as line segments clipped to the bar (hatch_bar), not
    with matplotlib's hatch=, which writes an SVG <pattern>. Illustrator
    dropped that pattern fill: the matched bar of the 2026-10-08 build reached
    the assembled figure as an empty outline.
  * No mathtext. pi_N, pi_0, pi_4 and the 10^k ticks of panel C were mathtext,
    which the SVG carries in DejaVu Sans with 4.9 pt subscripts and exponents
    (under the 6 pt floor); Arial has no subscript digits to replace them with.
    Labels are now plain text in the label font, as dN/dS is: "πN/πS",
    "π0/π4", "π4 (neutral diversity)". Panel C's log axis reads 0.001, 0.01.
  * "win." is "windows"; panel D's label reads "DoS (MAF ≥ 0.1)" (the MK row
    is still selected by its table key, "MAF>=0.1").
  * --inputs-from RUN_DIR rebuilds from the inputs a previous run recorded in
    its provenance.json, and checks panel A against the values recorded there
    (the bootstrap is seeded, so the same inputs must reproduce them exactly).
  * provenance.json records the panel B bars.


CHANGES, 2026-10-29 (figure only; every number is computed as before)
--------------------------------------------------------------------
  * Panel B bar order, after Jianbin Wang's review: retained genes, eliminated
    genes, then the hatched retained and eliminated pi4-matched windows, with
    a gap of PAIR_GAP bar pitches between the two pairs. The nested order
    (genes outside, windows inside) made the reader work out which bars were
    compared; now each comparison is a pair with retained left and eliminated
    right, as in every other figure. Bar labels name the class and n
    ("Retained" over "9,587 genes"), and each pair is named once beneath it
    ("All genes", "π4-matched windows", minor tick labels), which replaces the
    old "Retained, π4-matched (13 windows)" label. Values and hatching
    unchanged.


Example
-------
    # rebuild the 2026-10-08 figure from the inputs it recorded
    python figure_polymorphism.py \\
        --inputs-from results/figures/figure_polymorphism_20261008_pi_filtered_miss4 \\
        --outdir results/figures --tag pi_filtered_miss4

    # or name every input
    python figure_polymorphism.py \\
        --per-gene results/pi_filtered_miss4/per_gene.tsv \\
        --pi-summary results/pi_filtered_miss4/pi_summary.tsv \\
        --mk-summary results/pi_filtered_miss4/mk_summary.tsv \\
        --matched-summary results/pi4_matched_5k_miss4/pi4_matched_summary.tsv \\
        --windows results/pi4_matched_5k_miss4/windows.tsv \\
        --outdir results/figures --tag pi_filtered_miss4
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import transforms as mtransforms
from matplotlib.collections import LineCollection
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FuncFormatter, NullFormatter

# style_boot locates lib/paper_colors.py and warns if it cannot, rather than
# silently substituting a different palette.
from style_boot import RETAINED, ELIMINATED, paper_colors

if paper_colors is None:                # pragma: no cover
    raise SystemExit(
        "figure_polymorphism.py builds a journal figure and needs the full "
        "house style, not just the colours. Set PAPER_COLORS_DIR to the "
        "directory holding paper_colors.py.")

apply_style = paper_colors.apply_style
audit = paper_colors.audit
figure = paper_colors.figure
headroom = paper_colors.headroom
make_outdir = paper_colors.make_outdir
panel_letter = paper_colors.panel_letter
pc_save = paper_colors.save
FONT = paper_colors.FONT
LINEWIDTH = paper_colors.LINEWIDTH
NEUTRAL = paper_colors.NEUTRAL["ink_light"]
INK = paper_colors.NEUTRAL["ink"]
RULE = paper_colors.NEUTRAL["rule"]

#: Window bars in panel B are a sub-category of their region, so they keep the
#: region hue and are distinguished by texture (paper_colors: hue = region,
#: texture = subclass). The texture is drawn by hatch_bar(), not hatch=.
HATCH_SPACING_MM = 1.2

#: Panel B draws two comparisons side by side (every gene, then the
#: pi4-matched windows); the second pair starts this many bar pitches further
#: right, so the space between the pairs groups them on the page.
PAIR_GAP = 0.5

#: Panel B names each pair once, centred under it, below the bar labels.
GROUP_LABEL = {"genes": "All genes", "windows": "π4-matched windows"}

#: Figure 6 (was 7) geometry. Two rows of two panels at double-column width;
#: well inside the 200 mm practical cap that leaves room for the caption.
FIG_W = "double"
FIG_H_MM = 128.0

#: Windows panel C could not draw (pi4 = 0). Recorded rather than discarded:
#: an omitted point is indistinguishable from missing data on the page.
PANEL_C_DROPPED = 0

#: Inputs a run records and --inputs-from can reuse. --local-divergence is
#: left out on purpose: it builds a supplement that is no longer used
#: (supplementary fig. S16 comes from reference_bias_check.py), so a rebuild
#: should not fail on it. Give it explicitly to build that supplement.
REUSABLE_INPUTS = ("per_gene", "pi_summary", "mk_summary", "matched_summary",
                   "windows", "kaks")
REQUIRED_INPUTS = ("per_gene", "pi_summary", "mk_summary")
DEFAULT_N_BOOT = 2000
DEFAULT_MK_SET = "MAF>=0.1"


def save(fig, outdir, stem: str) -> None:
    """Save through paper_colors so timestamps are suppressed and runs diff."""
    written = pc_save(fig, str(outdir), stem)
    plt.close(fig)
    logging.info("  %s", ", ".join(Path(w).name for w in written))


def boot_fold(sub_r, sub_e, fn, n_boot, seed=0):
    """Bootstrap the eliminated/retained fold, resampling genes within class."""
    rng = np.random.default_rng(seed)
    ir, ie = np.arange(len(sub_r)), np.arange(len(sub_e))
    out = []
    for _ in range(n_boot):
        a = fn(sub_e.iloc[rng.choice(ie, len(ie), replace=True)])
        b = fn(sub_r.iloc[rng.choice(ir, len(ir), replace=True)])
        if a == a and b == b and b != 0:
            out.append(a / b)
    return (np.nan, np.nan) if not out else tuple(np.percentile(out, [2.5, 97.5]))


def pooled_pin_pis(d):
    S, N = d.S_sites.sum(), d.N_sites.sum()
    if S <= 0 or N <= 0:
        return np.nan
    piS = d.pi_syn.sum() / S
    return (d.pi_non.sum() / N) / piS if piS > 0 else np.nan


def pooled_pi0_pi4(d):
    S4, N0 = d.S4_sites.sum(), d.N0_sites.sum()
    if S4 <= 0 or N0 <= 0:
        return np.nan
    pi4 = d.pi_4.sum() / S4
    return (d.pi_0.sum() / N0) / pi4 if pi4 > 0 else np.nan


def decimal_tick(v, _pos=None):
    """Label a decade of a log axis as a plain number: 0.001, 0.01, 0.1, 1.

    The default log formatter writes 10^k in mathtext (DejaVu Sans, exponent
    at 70% of the tick size). Anything between decades is left unlabelled.
    """
    if not np.isfinite(v) or v <= 0:
        return ""
    k = np.log10(v)
    if abs(k - round(k)) > 1e-9:
        return ""
    k = int(round(k))
    return f"{10.0 ** k:.{max(0, -k)}f}"


def hatch_bar(ax, bar, spacing_mm=HATCH_SPACING_MM, color="black", lw=None):
    """Draw 45-degree hatching inside one bar as line segments clipped to it.

    Call on a laid-out figure (hatch_window_bars() draws first). The segments
    are laid out in inches from the bar's lower-left corner, so they run at 45
    degrees on the page whatever the axes' aspect; the corner is resolved at
    draw time, so a later layout pass carries them with the bar. A margin
    around the bar covers any small change in its size. Plain paths with a
    clip path survive the SVG import into Illustrator, where a matplotlib
    hatch (an SVG <pattern>) did not.
    """
    lw = LINEWIDTH["hairline"] if lw is None else lw
    fig = ax.figure
    (x0, y0), (x1, y1) = ax.transData.transform(
        [(bar.get_x(), bar.get_y()),
         (bar.get_x() + bar.get_width(), bar.get_y() + bar.get_height())])
    w_in = abs(x1 - x0) / fig.dpi
    h_in = abs(y1 - y0) / fig.dpi
    pad = 0.1                                   # inches beyond every edge
    step = spacing_mm / 25.4
    length = h_in + 2 * pad
    offsets = np.arange(-length - pad, w_in + pad + step, step)
    segments = [[(o, -pad), (o + length, -pad + length)] for o in offsets]
    anchor = mtransforms.ScaledTranslation(bar.get_x(), bar.get_y(),
                                           ax.transData)
    lines = LineCollection(segments, colors=color, linewidths=lw,
                           transform=fig.dpi_scale_trans + anchor,
                           zorder=bar.get_zorder() + 0.5)
    lines.set_capstyle("butt")
    ax.add_collection(lines, autolim=False)
    # Clip AFTER add_collection(), and to the bar's path rather than the bar:
    # set_clip_path(<Rectangle>) stores a clip box, not a clip path, and
    # add_collection() then replaces it with the axes box, which leaves the
    # hatch running across the whole panel.
    lines.set_clip_path(bar.get_path(), bar.get_transform())
    return lines


def hatch_window_bars(fig) -> int:
    """Hatch every bar panel_b() marked as a window bar; returns the count.

    Run once the figure is complete, so the hatch is laid out on the solved
    geometry. The hatch sits inside the axes and does not move the layout.
    """
    fig.canvas.draw()
    n = 0
    for ax in fig.axes:
        for patch in list(ax.patches):
            if getattr(patch, "_window_bar", False):
                hatch_bar(ax, patch)
                n += 1
    return n


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

def gene_dnds(pg, kaks):
    """
    Per-gene dN/dS restricted to the genes in this polymorphism analysis.

    pi_by_region.py already merges Ka and Ks into per_gene.tsv when it is given
    --kaks, so that is the first choice: it guarantees the dN/dS genes are
    exactly the piN/piS genes. The standalone --kaks table is only a fallback
    for a per_gene.tsv built without it.
    """
    if {"Ka", "Ks"}.issubset(pg.columns):
        d = pg[pg.Ka.notna() & pg.Ks.notna() & (pg.Ks > 0)].copy()
        if len(d):
            d["dnds"] = d.Ka / d.Ks
            return d
    if kaks is None:
        return None
    k = kaks.rename(columns={"gene1": "gene"})
    if "Ka_Ks" in k.columns:
        k = k[["gene", "Ka_Ks"]].rename(columns={"Ka_Ks": "dnds"})
    elif {"Ka", "Ks"}.issubset(k.columns):
        k = k[["gene", "Ka", "Ks"]]
        k = k[k.Ks > 0]
        k["dnds"] = k.Ka / k.Ks
        k = k[["gene", "dnds"]]
    else:
        logging.warning("--kaks table has neither Ka_Ks nor Ka+Ks; "
                        "skipping the dN/dS row of panel A")
        return None
    d = pg.merge(k, on="gene", how="inner")
    d = d[d.dnds.notna()]
    return d if len(d) else None


def panel_a(ax, pg, kaks, matched, n_boot, rows=None):
    """
    Forest of eliminated/retained fold across estimators.

    `rows` lets a caller pass back the list returned by an earlier call so the
    individual-panel export does not repeat several thousand bootstraps.
    """
    if rows is None:
        ret = pg[pg.region == "retained"]
        eli = pg[pg.region == "eliminated"]
        rows = []

        # dN/dS on exactly the genes in this analysis, median per gene (the
        # Ka/Ks paper's convention), bootstrapped over genes.
        k = gene_dnds(pg, kaks)
        if k is not None:
            kr, ke = k[k.region == "retained"], k[k.region == "eliminated"]
            if len(kr) and len(ke):
                fold = ke.dnds.median() / kr.dnds.median()
                lo, hi = boot_fold(kr, ke, lambda d: d.dnds.median(), n_boot)
                rows.append(("dN/dS\n(divergence)", fold, lo, hi, len(ke),
                             INK))

        fold = pooled_pin_pis(eli) / pooled_pin_pis(ret)
        lo, hi = boot_fold(ret, eli, pooled_pin_pis, n_boot)
        rows.append(("πN/πS", fold, lo, hi, len(eli), INK))

        if {"S4_sites", "N0_sites"}.issubset(pg.columns):
            fold = pooled_pi0_pi4(eli) / pooled_pi0_pi4(ret)
            lo, hi = boot_fold(ret, eli, pooled_pi0_pi4, n_boot)
            rows.append(("π0/π4\n(no NG weighting)", fold, lo, hi,
                         len(eli), INK))

        # pi4-matched contrast, from pi4_matched_control.py
        if matched is not None:
            m = matched[matched.test == "matched"]
            if len(m) == 2:
                r = m[m.region == "retained"].iloc[0]
                e = m[m.region == "eliminated"].iloc[0]
                # Fold CI approximated from the two independent window
                # bootstraps; deliberately conservative (it is wider than a
                # joint resample would give).
                rows.append(("π0/π4\n(π4-matched)",
                             e.pi0_pi4 / r.pi0_pi4,
                             e["lo"] / r["hi"], e["hi"] / r["lo"],
                             int(e.windows), INK))

    y = np.arange(len(rows))[::-1]
    for (label, f, lo, hi, n, c), yy in zip(rows, y):
        if lo == lo and hi == hi:
            ax.plot([lo, hi], [yy, yy], color=c,
                    linewidth=LINEWIDTH["data"], solid_capstyle="butt")
        ax.plot([f], [yy], "o", color=c, markersize=4.5,
                markeredgecolor="black", markeredgewidth=LINEWIDTH["hairline"],
                zorder=3)
    ax.axvline(1.0, ls="--", lw=LINEWIDTH["rule"], color=RULE)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel("Eliminated / retained fold")
    finite = [v for r in rows for v in (r[1], r[2], r[3]) if v == v]
    if finite:
        span = max(finite) - min(min(finite), 1.0)
        ax.set_xlim(min(min(finite), 1.0) - 0.08 * span,
                    max(finite) + 0.08 * span)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    return rows


def panel_b(ax, pi_summary, matched):
    """
    The Ne control, every bar on ONE estimator (pi0/pi4). Two contrasts, side
    by side:

      left pair   every gene in each class (pi_summary.tsv) -- the 1.82-fold
                  "no NG weighting" row of panel A
      right pair  the pi4-matched comparison of pi4_matched_control.py: the
                  retained windows whose pi4 falls in the range of the
                  eliminated windows, against those eliminated windows -- the
                  1.44-fold row of panel A. Hatched: windows, not genes.

    Within each pair, retained on the left and eliminated on the right; the
    pairs are PAIR_GAP bar pitches apart (Jianbin Wang's order, 2026-10-29;
    until then the pairs were nested, genes outside and windows inside).
    Until 2026-10-14 the eliminated-windows bar was missing, so the matched
    retained windows stood beside the eliminated GENES and the panel read as a
    1.53-fold contrast.

    Everything here is pi0/pi4, not piN/piS. The matched contrast is only
    available as pi0/pi4 (pi4_matched_control.py matches windows on pi4), and
    mixing a pi0/pi4 bar into a piN/piS axis would be the same false comparison
    panel A is built to avoid. pi0/pi4 is also the stronger estimator for this
    panel: it carries no Nei-Gojobori weighting.

    Window bars are marked (`_window_bar`); hatch_window_bars() hatches them
    once the figure is laid out. Returns the bars drawn, as dicts.
    """
    r = pi_summary[pi_summary.region == "retained"].iloc[0]
    e = pi_summary[pi_summary.region == "eliminated"].iloc[0]
    # Left pair: every gene in each class.
    bars = [dict(label=f"Retained\n{int(r.genes):,} genes", region="retained",
                 unit="genes", n=int(r.genes), value=float(r.pi0_pi4),
                 lo=float(r.pi0_pi4_lo), hi=float(r.pi0_pi4_hi)),
            dict(label=f"Eliminated\n{int(e.genes):,} genes",
                 region="eliminated", unit="genes", n=int(e.genes),
                 value=float(e.pi0_pi4), lo=float(e.pi0_pi4_lo),
                 hi=float(e.pi0_pi4_hi))]

    # Right pair: the pi4-matched windows, retained then eliminated.
    if matched is not None:
        m = matched[matched.test == "matched"]
        rm, em = m[m.region == "retained"], m[m.region == "eliminated"]
        if len(rm) == 1 and len(em) == 1:
            rm, em = rm.iloc[0], em.iloc[0]
            bars.append(dict(
                label=f"Retained\n{int(rm.windows)} windows",
                region="retained", unit="windows", n=int(rm.windows),
                value=float(rm.pi0_pi4), lo=float(rm["lo"]),
                hi=float(rm["hi"])))
            bars.append(dict(
                label=f"Eliminated\n{int(em.windows)} windows",
                region="eliminated", unit="windows", n=int(em.windows),
                value=float(em.pi0_pi4), lo=float(em["lo"]),
                hi=float(em["hi"])))
        else:
            logging.warning("Panel B: the matched summary has %d retained and "
                            "%d eliminated 'matched' rows (need one each); "
                            "drawing the gene bars only", len(rm), len(em))

    vals = np.array([b["value"] for b in bars], float)
    los = np.array([b["lo"] for b in bars], float)
    his = np.array([b["hi"] for b in bars], float)
    cols = [RETAINED if b["region"] == "retained" else ELIMINATED for b in bars]
    # One bar pitch apart within a pair; the window pair is shifted right by
    # PAIR_GAP so the two comparisons read as separate pairs.
    x = np.arange(len(bars), dtype=float)
    x[np.array([b["unit"] == "windows" for b in bars], dtype=bool)] += PAIR_GAP
    drawn = ax.bar(x, vals, color=cols, edgecolor="black",
                   linewidth=LINEWIDTH["thin"], width=0.62)
    for patch, b in zip(drawn, bars):
        patch._window_bar = b["unit"] == "windows"
    ax.errorbar(x, vals,
                yerr=[np.clip(vals - los, 0, None), np.clip(his - vals, 0, None)],
                fmt="none", ecolor="black", capsize=2.0,
                linewidth=LINEWIDTH["thin"])
    ax.set_xticks(x)
    ax.set_xticklabels([b["label"] for b in bars])
    # Each pair is named once, centred under it ("All genes", "pi4-matched
    # windows"). Minor tick labels, so the layout reserves their room like any
    # tick label; they sit below the two-line bar labels: major tick length +
    # major pad + the bar labels' height + a 3 pt gap, enough that the pair
    # name does not read as a third line of the bar label above it.
    units = [u for u in ("genes", "windows") if any(b["unit"] == u for b in bars)]
    centres = [float(np.mean(x[np.array([b["unit"] == u for b in bars])]))
               for u in units]
    size = FontProperties(size=plt.rcParams["xtick.labelsize"]).get_size_in_points()
    n_lines = max(b["label"].count("\n") + 1 for b in bars)
    pad = (plt.rcParams["xtick.major.size"] + plt.rcParams["xtick.major.pad"]
           + n_lines * size * 1.2 + 3.0)
    ax.set_xticks(centres, minor=True)
    ax.set_xticklabels([GROUP_LABEL[u] for u in units], minor=True)
    ax.tick_params(axis="x", which="minor", length=0, pad=pad, labelsize=size)
    ax.set_ylabel("π0/π4")
    ax.set_ylim(bottom=0)
    return bars


def panel_c(ax, windows):
    """pi0/pi4 versus pi4 across windows, with the retained trend."""
    # A log x-axis and a log-scale fit both require pi4 > 0. Windows with
    # pi4 == 0 carry no information about the trend anyway, so drop them here
    # rather than letting log() produce -inf and take polyfit down with it.
    ok = (windows.pi4 > 0) & windows.pi0_pi4.notna()
    global PANEL_C_DROPPED
    dropped = int((~ok).sum())
    PANEL_C_DROPPED = dropped
    if dropped:
        logging.warning("Panel C: dropping %d window(s) with pi4 = 0 or "
                        "undefined pi0/pi4", dropped)
    windows = windows[ok]
    ret = windows[windows.region == "retained"]
    eli = windows[windows.region == "eliminated"]
    ax.scatter(ret.pi4, ret.pi0_pi4, s=3.0, color=RETAINED, alpha=0.5,
               edgecolor="none", rasterized=True,
               label=f"Retained ({len(ret)} windows)")
    ax.scatter(eli.pi4, eli.pi0_pi4, s=12.0, color=ELIMINATED,
               edgecolor="black", linewidth=LINEWIDTH["hairline"],
               label=f"Eliminated ({len(eli)} windows)")
    if len(ret) >= 3:
        x = np.log(ret.pi4.to_numpy())
        slope, intercept = np.polyfit(x, ret.pi0_pi4.to_numpy(), 1)
        xs = np.geomspace(ret.pi4.min(), ret.pi4.max(), 100)
        ax.plot(xs, intercept + slope * np.log(xs), color=RETAINED,
                lw=LINEWIDTH["data"], ls="--")
    ax.set_xscale("log")
    # Plain-number decades (0.001, 0.01), not mathtext 10^k.
    ax.xaxis.set_major_formatter(FuncFormatter(decimal_tick))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("π4 (neutral diversity)")
    ax.set_ylabel("π0/π4")
    # The trend descends, so data occupy upper-left (low pi4, high ratio) and
    # lower-right. headroom() reserves a blank strip along the top FIRST, and
    # the legend goes into it. Never loc="best": it re-picks at draw time and
    # moves between runs. A legend OUTSIDE the axes is the house default, but
    # constrained_layout resolves an over-wide outside legend by shrinking the
    # panel, and at 88 mm this one would.
    top = float(np.nanmax(windows.pi0_pi4)) if len(windows) else 1.0
    ax.set_ylim(0, top)
    headroom(ax, 0.30)
    ax.legend(loc="upper right", frameon=False, fontsize=FONT["legend"],
              handlelength=1.2, handletextpad=0.4, borderaxespad=0.3,
              scatterpoints=1)


def panel_d(ax, mk, variant_set):
    """DoS by class for one variant set."""
    m = mk[mk.variant_set == variant_set]
    if m.empty:
        ax.text(0.5, 0.5, f"no MK rows for {variant_set}", ha="center",
                transform=ax.transAxes)
        return
    vals, los, his, labels, cols = [], [], [], [], []
    for region, c in (("retained", RETAINED), ("eliminated", ELIMINATED)):
        row = m[m.region == region]
        if row.empty:
            continue
        row = row.iloc[0]
        vals.append(row.DoS); los.append(row.DoS_lo); his.append(row.DoS_hi)
        labels.append(f"{region.capitalize()}\n(n = {int(row.genes):,})")
        cols.append(c)
    if not vals:
        ax.text(0.5, 0.5, "no MK rows", ha="center", transform=ax.transAxes)
        return
    vals = np.asarray(vals, float)
    los = np.asarray(los, float)
    his = np.asarray(his, float)
    x = np.arange(len(vals))
    ax.bar(x, vals, color=cols, edgecolor="black",
           linewidth=LINEWIDTH["thin"], width=0.55)
    ax.errorbar(x, vals,
                yerr=[np.clip(vals - los, 0, None), np.clip(his - vals, 0, None)],
                fmt="none", ecolor="black", capsize=2.0,
                linewidth=LINEWIDTH["thin"])
    ax.axhline(0, color=RULE, lw=LINEWIDTH["rule"])
    ax.set_xticks(x); ax.set_xticklabels(labels)
    # The table key is "MAF>=0.1"; the label reads "MAF ≥ 0.1".
    ax.set_ylabel(f"DoS ({variant_set.replace('>=', ' ≥ ')})")


def supp_local_divergence(ld, outdir):
    """Supplementary: piN/piS by local-divergence stratum."""
    fig, axes = figure("double", height_mm=62.0, ncols=2)
    order = list(dict.fromkeys(ld.stratum))
    for region, colour in (("retained", RETAINED), ("eliminated", ELIMINATED)):
        s = ld[ld.region == region].set_index("stratum").reindex(order).dropna(
            subset=["piS"])
        xs = [order.index(i) for i in s.index]
        axes[0].plot(xs, s.piS, "o-", color=colour,
                     lw=LINEWIDTH["data"], ms=3.0, label=region.capitalize())
        axes[1].plot(xs, s.piN_piS, "o-", color=colour,
                     lw=LINEWIDTH["data"], ms=3.0, label=region.capitalize())
    for a in axes:
        a.set_xticks(range(len(order)))
        a.set_xticklabels(order, rotation=45, ha="right")
        a.set_xlabel("Variants within ±150 bp")
        a.set_ylim(bottom=0)
    axes[0].set_ylabel("πS")
    axes[1].set_ylabel("πN/πS")
    # constrained_layout is on via apply_style(); tight_layout() fights it and
    # matplotlib warns that the layout silently changed engine.
    headroom(axes[1], 0.20)
    axes[1].legend(loc="upper right", frameon=False, fontsize=FONT["legend"])
    for a, letter in zip(axes, "AB"):
        panel_letter(a, letter)
    audit(fig)
    save(fig, outdir, "S_pi_by_local_divergence")


# ---------------------------------------------------------------------------
# Rebuilding from a previous run
# ---------------------------------------------------------------------------

def load_previous(path: Path):
    """(path of provenance.json, its contents) for --inputs-from."""
    path = path.expanduser()
    prov = path / "provenance.json" if path.is_dir() else path
    if not prov.is_file():
        raise SystemExit(f"--inputs-from: no provenance.json at {prov}")
    try:
        return prov.resolve(), json.loads(prov.read_text())
    except json.JSONDecodeError as exc:
        raise SystemExit(f"--inputs-from: {prov} is not valid JSON ({exc})")


def resolve_recorded(recorded: str, prov: Path) -> Path:
    """
    An input path as a previous run recorded it, made absolute.

    A relative path was relative to the directory that run started in, which
    provenance.json does not record. Try the current directory, then the run
    directory and each of its ancestors (a run started at the top of the
    project wrote results/... paths). Exits if none exists.
    """
    p = Path(recorded).expanduser()
    if p.is_absolute():
        if not p.exists():
            raise SystemExit(f"--inputs-from: recorded input not found: {p}")
        return p
    candidates = [Path.cwd() / p, prov.parent / p] + [a / p for a in
                                                       prov.parent.parents]
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise SystemExit(
        f"--inputs-from: cannot find the recorded input {recorded!r}; looked "
        f"in the current directory and every ancestor of {prov.parent}. "
        "Give it explicitly.")


def fill_from_previous(args, prov: Path, previous: dict) -> None:
    """Fill every input not given on the command line from the previous run."""
    recorded = previous.get("inputs") or {}
    for key in REUSABLE_INPUTS:
        flag = "--" + key.replace("_", "-")
        if getattr(args, key) is not None:
            logging.info("  %-18s %s (given)", flag, getattr(args, key))
        elif key in recorded:
            setattr(args, key, resolve_recorded(recorded[key], prov))
            logging.info("  %-18s %s (from the previous run)", flag,
                         getattr(args, key))
    if args.n_boot is None and previous.get("n_boot") is not None:
        args.n_boot = int(previous["n_boot"])
        logging.info("  %-18s %d (from the previous run)", "--n-boot",
                     args.n_boot)
    if args.mk_variant_set is None and previous.get("mk_variant_set"):
        args.mk_variant_set = str(previous["mk_variant_set"])
        logging.info("  %-18s %s (from the previous run)", "--mk-variant-set",
                     args.mk_variant_set)


def compare_panel_a(rows, previous_rows, rtol=1e-9):
    """
    Check panel A against the values a previous run recorded, row by row.

    Returns (same, largest relative difference, notes). Rows are matched by
    position, because the row labels changed on 2026-10-14 (mathtext to
    plain text) while the estimators and their order did not.
    """
    notes = []
    if len(rows) != len(previous_rows):
        return False, np.nan, [f"{len(rows)} rows now, {len(previous_rows)} "
                               "in the previous run"]
    worst, same = 0.0, True
    for new, old in zip(rows, previous_rows):
        name = new[0].replace("\n", " ")
        for i, key in ((1, "fold"), (2, "lo"), (3, "hi")):
            a, b = float(new[i]), float(old.get(key, np.nan))
            if np.isnan(a) and np.isnan(b):
                continue
            if np.isnan(a) or np.isnan(b):
                same = False
                notes.append(f"{name}: {key} {a} now, {b} before")
                continue
            rel = abs(a - b) / max(abs(b), 1e-300)
            worst = max(worst, rel)
            if rel > rtol:
                same = False
                notes.append(f"{name}: {key} {a:.6f} now, {b:.6f} before")
        if int(new[4]) != int(old.get("n", -1)):
            same = False
            notes.append(f"{name}: n {int(new[4])} now, {old.get('n')} before")
    return same, worst, notes


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assemble the polymorphism figure and supplements.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--inputs-from", type=Path,
                   help="A previous figure_polymorphism run directory (or its "
                        "provenance.json). Its recorded inputs, --n-boot and "
                        "--mk-variant-set fill any not given here, and panel "
                        "A is checked against the values it recorded.")
    p.add_argument("--per-gene", type=Path,
                   help="pi_by_region per_gene.tsv (required unless "
                        "--inputs-from supplies it)")
    p.add_argument("--pi-summary", type=Path,
                   help="pi_by_region pi_summary.tsv (required unless "
                        "--inputs-from supplies it)")
    p.add_argument("--mk-summary", type=Path,
                   help="pi_by_region mk_summary.tsv (required unless "
                        "--inputs-from supplies it)")
    p.add_argument("--matched-summary", type=Path,
                   help="pi4_matched_control pi4_matched_summary.tsv")
    p.add_argument("--windows", type=Path,
                   help="pi4_matched_control windows.tsv")
    p.add_argument("--kaks", type=Path)
    p.add_argument("--local-divergence", type=Path,
                   help="Builds the S_pi_by_local_divergence supplement; "
                        "never taken from --inputs-from")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--mk-variant-set", default=None,
                   help="Which MK row to draw in panel D (default "
                        f"{DEFAULT_MK_SET!r}). Use the MAF-filtered one: the "
                        "unfiltered alpha comparison is partly circular (see "
                        "the MK notes).")
    p.add_argument("--n-boot", type=int, default=None,
                   help=f"bootstrap replicates (default {DEFAULT_N_BOOT})")
    p.add_argument("--tag", default=None,
                   help="Suffix for the dated output directory, so a rerun "
                        "cannot overwrite the run a figure came from")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    # matplotlib subsets fonts through fontTools on every PDF save, which emits
    # a few hundred INFO lines and buries the panel-A values.
    for noisy in ("fontTools", "matplotlib.font_manager", "matplotlib"):
        logging.getLogger(noisy).setLevel(logging.WARNING)

    prov_prev, previous = None, None
    if args.inputs_from is not None:
        prov_prev, previous = load_previous(args.inputs_from)
        logging.info("Inputs from the previous run: %s", prov_prev)
        if previous.get("argv"):
            logging.info("  its command: %s", " ".join(map(str,
                                                           previous["argv"])))
        fill_from_previous(args, prov_prev, previous)
    missing = [k for k in REQUIRED_INPUTS if getattr(args, k) is None]
    if missing:
        raise SystemExit("missing " + ", ".join("--" + k.replace("_", "-")
                                                for k in missing)
                         + " (give them, or --inputs-from a previous run)")
    if args.n_boot is None:
        args.n_boot = DEFAULT_N_BOOT
    if args.mk_variant_set is None:
        args.mk_variant_set = DEFAULT_MK_SET

    apply_style()

    outdir = Path(make_outdir("figure_polymorphism",
                              base=str(args.outdir.expanduser().resolve()),
                              tag=args.tag))
    logging.info("Output: %s", outdir)

    pg = pd.read_csv(args.per_gene.expanduser(), sep="\t")
    pi_summary = pd.read_csv(args.pi_summary.expanduser(), sep="\t")
    mk = pd.read_csv(args.mk_summary.expanduser(), sep="\t")
    matched = (pd.read_csv(args.matched_summary.expanduser(), sep="\t")
               if args.matched_summary else None)
    windows = (pd.read_csv(args.windows.expanduser(), sep="\t")
               if args.windows else None)
    # Read the whole Ka/Ks table: gene_dnds() picks whichever of Ka_Ks or
    # Ka+Ks the file actually carries, so a column-name change upstream is a
    # warning here rather than a KeyError.
    kaks = pd.read_csv(args.kaks.expanduser(), sep="\t") if args.kaks else None

    logging.info("Building main figure at %s width, %.0f mm high",
                 FIG_W, FIG_H_MM)
    fig, axes = figure(FIG_W, height_mm=FIG_H_MM, nrows=2, ncols=2)
    (axA, axB), (axC, axD) = axes

    rows = panel_a(axA, pg, kaks, matched, args.n_boot)
    bars = panel_b(axB, pi_summary, matched)
    if windows is not None:
        panel_c(axC, windows)
    else:
        axC.set_visible(False)
    panel_d(axD, mk, args.mk_variant_set)

    for ax, letter in ((axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")):
        if ax.get_visible():
            panel_letter(ax, letter)
    logging.info("Panel B: %d window bar(s) hatched", hatch_window_bars(fig))

    problems = audit(fig)
    if problems:
        logging.warning("Figure 6 breaks %d journal rule(s) — listed above. "
                        "Fix before assembling.", len(problems))
    save(fig, outdir, "Figure_polymorphism")

    # individual panels, for rebuilding one at a time in Illustrator
    for name, fn in (("panelA_fold_forest",
                      lambda a: panel_a(a, pg, kaks, matched, args.n_boot,
                                        rows=rows)),
                     ("panelB_pi0pi4_matched",
                      lambda a: panel_b(a, pi_summary, matched)),
                     ("panelD_dos", lambda a: panel_d(a, mk, args.mk_variant_set))):
        f, a = figure("single", height_mm=62.0)
        fn(a)
        hatch_window_bars(f)
        audit(f)
        save(f, outdir, name)
    if windows is not None:
        f, a = figure("single", height_mm=62.0)
        panel_c(a, windows)
        audit(f)
        save(f, outdir, "panelC_pi0pi4_vs_pi4")

    if args.local_divergence:
        logging.info("Building supplements")
        supp_local_divergence(
            pd.read_csv(args.local_divergence.expanduser(), sep="\t"), outdir)

    logging.info("")
    logging.info("Panel A values (eliminated / retained fold):")
    for label, f, lo, hi, n, _ in rows:
        logging.info("  %-28s %.3f [%.3f, %.3f]",
                     label.replace("\n", " "), f, lo, hi)
    logging.info("Panel B values (pi0/pi4):")
    for b in bars:
        logging.info("  %-34s %.4f [%.4f, %.4f]",
                     b["label"].replace("\n", " "), b["value"], b["lo"],
                     b["hi"])
    by = {(b["region"], b["unit"]): b["value"] for b in bars}
    for unit in ("genes", "windows"):
        if ("retained", unit) in by and ("eliminated", unit) in by:
            logging.info("  eliminated / retained, %-7s %.3f", unit,
                         by[("eliminated", unit)] / by[("retained", unit)])

    reproduces = None
    if previous is not None:
        if previous.get("rows"):
            reproduces, worst, notes = compare_panel_a(rows, previous["rows"])
            if reproduces:
                logging.info("Panel A reproduces the previous run (%d rows; "
                             "largest relative difference %.1e).", len(rows),
                             worst)
            else:
                logging.warning("Panel A DIFFERS from the previous run:")
                for n_ in notes:
                    logging.warning("  %s", n_)
        else:
            logging.warning("The previous run recorded no panel A rows; "
                            "nothing to compare.")

    (outdir / "provenance.json").write_text(json.dumps({
        "argv": sys.argv,
        "inputs_from": str(prov_prev) if prov_prev else None,
        "inputs": {k: str(getattr(args, k)) for k in
                   ("per_gene", "pi_summary", "mk_summary", "matched_summary",
                    "windows", "kaks", "local_divergence")
                   if getattr(args, k, None)},
        "rows": [{"estimator": r[0].replace("\n", " "), "fold": r[1],
                  "lo": r[2], "hi": r[3], "n": r[4]} for r in rows],
        "panel_a_reproduces_previous": reproduces,
        "panel_b": [{"bar": b["label"].replace("\n", " "),
                     "region": b["region"], "unit": b["unit"], "n": b["n"],
                     "pi0_pi4": b["value"], "lo": b["lo"], "hi": b["hi"]}
                    for b in bars],
        "panel_c_windows_dropped": PANEL_C_DROPPED,
        "mk_variant_set": args.mk_variant_set,
        "n_boot": args.n_boot,
    }, indent=2, default=str))

    logging.info("")
    logging.info("Panel C dropped %d window(s) with pi4 = 0 — quote this in "
                 "the caption; a reader cannot see an omitted point.",
                 PANEL_C_DROPPED)
    logging.info("Export from the SVG — Arial is absent on Isaac compute "
                 "nodes, so PNG/PDF render in a fallback font.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
