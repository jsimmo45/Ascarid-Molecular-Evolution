#!/usr/bin/env python3
"""
figure_cohort_integrity.py — supplementary fig. S17: the result depends on no
single worm and no library batch. (S19 before the 2026-10-12 renumbering; the
output stem keeps S19_cohort_integrity so the assembly manifest stays valid.)

Changes 2026-10-31 (S17 at 185 mm, like every other figure)
------------------------------------------------------------
Every main and supplementary figure is 185 mm wide except S17, which was laid
out by hand at 136 mm because the 185 mm build wasted space. That build used
two equal columns and started D's three-entry legend at D's left spine. The
legend is wider than D, so constrained_layout reserved room for it down the
whole right side: B and D came out 48 mm wide with a 31 mm blank strip to
their right, and the legend still ran past the edge (audit() does not look
for that, and reported the figure clean).

* D's legend ends at D's right edge and runs left over the column of
  held-out names, as in the hand layout of 2026-10-30.
* The columns are sized so A's square matrix fills the left column (A over
  C): fit_left_column() measures, after a draw, the height row 1 gives its
  axes and the width the left column gives its axes, and sets the width ratio
  that makes them equal (two passes on the real data). B and D take the rest.
* Panel letters follow the hand layout: one height per row, A and C at the
  figure's left edge, B and D at the left of their own labels
  (place_letters()). They are placed after the layout is solved and frozen.
* outside_canvas() checks that no label, legend or letter runs past the edge,
  and the log and provenance report what it finds.
* Statistics and the three input checks are unchanged.

What this assembles, and why it is one figure
---------------------------------------------
Three scripts each answer part of one question — is the cohort what it claims to
be, and does the headline survive its composition?

    A  check_duplicate_individuals.py  are any two libraries the same worm?
    B  het_by_region.py                does heterozygosity vary between worms?
    C  het_by_region.py                does it vary BETWEEN REGION CLASSES?
    D  jackknife_individuals.py        does the fold survive dropping worms?

Separately each is a loose end. Together they are an argument with a shape: the
cohort is ten distinct worms (A); they differ in library quality (B); that
difference is genome-wide rather than class-specific, so it cannot generate a
contrast between classes (C); and dropping any worm, or any batch, leaves the
fold where it was (D). C is the hinge — B on its own looks alarming and means
nothing.

This is an ASSEMBLER: it reads the three summary tables and draws, rather than
recomputing anything. The source scripts keep their own working figures for
analysis; this builds the journal panel. Nothing here can disagree with them,
because nothing here recalculates.

Example
-------
    python figure_cohort_integrity.py \\
        --pairwise  results/duplicates_filtered_miss4/pairwise.tsv \\
        --het       results/het_by_region_filtered_miss4/het_by_region.tsv \\
        --jackknife results/jackknife_filtered_miss4/jackknife_summary.tsv \\
        --outdir    results/figures
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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerTuple

from style_boot import RETAINED, ELIMINATED, paper_colors

if paper_colors is None:                # pragma: no cover
    raise SystemExit(
        "figure_cohort_integrity.py builds a journal supplement and needs the "
        "full house style. Set PAPER_COLORS_DIR to the directory holding "
        "paper_colors.py.")

apply_style = paper_colors.apply_style
audit = paper_colors.audit
figure = paper_colors.figure
headroom = paper_colors.headroom
make_outdir = paper_colors.make_outdir
pc_save = paper_colors.save
FONT = paper_colors.FONT
LINEWIDTH = paper_colors.LINEWIDTH
MM_PER_IN = paper_colors.MM_PER_IN
INK = paper_colors.NEUTRAL["ink"]
RULE = paper_colors.NEUTRAL["rule"]

#: One legend placement for every panel: above the axes, left-aligned. Shared
#: so panel C, which builds its own legend for a composite key, cannot drift
#: from the others. D alone is right-aligned (see main()), because its legend
#: is wider than its axes.
LEGEND_ABOVE = dict(loc="lower left", bbox_to_anchor=(0.0, 1.02),
                    frameon=False, fontsize=FONT["legend"], handlelength=1.4,
                    handletextpad=0.5, columnspacing=1.0, borderaxespad=0.0)

#: 185 mm, the width of every other main and supplementary figure.
FIG_W = "double"
FIG_H_MM = 132.0

#: Starting share of the axes width for the left column (A over C). Only the
#: first guess: fit_left_column() solves the share at which A's square matrix
#: fills its column, and B and D take the remainder.
LEFT_SHARE_START = 0.38
#: Stop adjusting once the left column and row 1 agree to within this.
FIT_TOL_MM = 0.2

#: Discordance is a MAGNITUDE, not an identity, so it takes the house
#: sequential ramp rather than a categorical map. Using viridis or similar would
#: put blue and red into a panel where those hues are reserved for region.
#: pair_ramp() is the house sequential ramp. getattr rather than a direct
#: attribute so an older paper_colors without it degrades to the same six
#: values rather than failing to import the whole figure.
_RAMP = getattr(paper_colors, "pair_ramp", None)
_RAMP_HEX = _RAMP(6) if _RAMP else ["#CFC7E2", "#B3A6D4", "#9784C4",
                                    "#7B62B0", "#5F4396", "#432873"]
if _RAMP is None:
    logging.warning("paper_colors has no pair_ramp(); using the literal house "
                    "sequential values for the panel A colour scale")
DISCORD_CMAP = LinearSegmentedColormap.from_list(
    "pair_ramp", list(reversed(_RAMP_HEX)))


def short(name: str, n: int = 4) -> str:
    """Final n characters of a run accession.

    Ten full SRR accessions on both axes of an 88 mm matrix leaves no room for
    the matrix. Within this cohort the last four digits are unique; the caption
    states the convention.
    """
    return str(name)[-n:]


def panel_a(ax, pw: pd.DataFrame, dup_threshold: float, names) -> dict:
    """Pairwise genotype discordance, as a matrix."""
    idx = {s: i for i, s in enumerate(names)}
    m = np.full((len(names), len(names)), np.nan)
    for r in pw.itertuples():
        i, j = idx[r.sample_a], idx[r.sample_b]
        m[i, j] = m[j, i] = r.discordance
    im = ax.imshow(m, cmap=DISCORD_CMAP, aspect="equal")
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels([short(s) for s in names], rotation=90)
    ax.set_yticklabels([short(s) for s in names])
    ax.tick_params(length=1.5, pad=1.5)
    # The diagonal is undefined, not zero. Left blank rather than filled with 0,
    # which would read as a perfect match and is exactly the wrong impression in
    # a panel about duplicate detection.
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(True)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("Genotype discordance", fontsize=FONT["axis_label"])
    cb.ax.tick_params(labelsize=FONT["tick"], length=1.5)
    cb.outline.set_linewidth(LINEWIDTH["spine"])
    lo = float(pw.discordance.min())
    return {"n_samples": len(names), "min_discordance": lo,
            "flagged": int((pw.discordance <= dup_threshold).sum())}


def panel_bc(axB, axC, het: pd.DataFrame, tol: float, names) -> dict:
    """Heterozygosity per worm by class, and the class ratio."""
    # Accession order, matching A and D, NOT sorted by heterozygosity. Sorting
    # B and C by their own values makes each panel prettier and makes the figure
    # unreadable: a reader tracking one worm across four panels has to find it
    # again in every one. Accession order also puts the sequencing batches
    # adjacent, which is what panel D's grouped drops are about.
    d = het.set_index("sample").reindex(names).reset_index()
    x = np.arange(len(d), dtype=float)
    w = 0.38
    axB.bar(x - w / 2, d.het_retained, w, color=RETAINED, edgecolor="black",
            linewidth=LINEWIDTH["thin"], label="Retained")
    axB.bar(x + w / 2, d.het_eliminated, w, color=ELIMINATED,
            edgecolor="black", linewidth=LINEWIDTH["thin"], label="Eliminated")
    axB.set_ylabel("Heterozygous fraction")
    axB.set_ylim(0, float(max(d.het_retained.max(), d.het_eliminated.max())))
    headroom(axB, 0.10)

    med = float(d.het_ratio_eli_ret.median())
    band = axC.axhspan(med * (1 - tol), med * (1 + tol), color=RULE,
                       alpha=0.12, lw=0, zorder=1)
    line = axC.axhline(med, color=RULE, ls="--", lw=LINEWIDTH["rule"],
                       zorder=2)
    # The ratio belongs to a worm, not to a region, so it is drawn in ink as in
    # panel D. Region red here would read as "eliminated values".
    axC.scatter(x, d.het_ratio_eli_ret, s=12.0, color=INK,
                edgecolor="black", linewidth=LINEWIDTH["hairline"], zorder=3)
    axC.set_ylabel("Eliminated / retained\nheterozygous fraction")
    axC.set_ylim(0, float(d.het_ratio_eli_ret.max()))
    headroom(axC, 0.10)
    # The band is drawn, so it needs a key. One composite entry -- the dashed
    # line over a grey swatch, as plotted -- rather than a second entry: a
    # two-entry row is wider than the panel, and constrained_layout answers
    # that by shrinking all four panels (about 47 -> 37 mm in testing).
    axC.legend([(band, line)],
               [f"Cohort median ({med:.2f}) ± {100 * tol:.0f}%"],
               handler_map={tuple: HandlerTuple(ndivide=1)}, ncol=1,
               **LEGEND_ABOVE)

    for a in (axB, axC):
        a.set_xticks(x)
        a.set_xticklabels([short(s) for s in d["sample"]], rotation=90)
        a.set_xlim(-0.7, len(d) - 0.3)
        a.tick_params(axis="x", length=1.5, pad=1.5)
    off = d[(d.het_ratio_eli_ret / med - 1).abs() > tol]
    return {"cohort_median_ratio": med,
            "het_spread_fold": float(d.het_retained.max() / d.het_retained.min()),
            "outside_band": [str(s) for s in off["sample"]]}


def panel_d(ax, jk: pd.DataFrame, names) -> dict:
    """Leave-one-out and grouped-drop folds against the full cohort."""
    full = jk[jk.dropped == "(none)"]
    if full.empty:
        ax.text(0.5, 0.5, "no full-cohort row", ha="center",
                transform=ax.transAxes)
        return {}
    full_fold = float(full.fold_piN_piS.iloc[0])
    d = jk[jk.dropped != "(none)"].copy()
    # Grouped drops are a different question from leave-one-out — they are the
    # only thing here that can see a SHARED artifact — so they are separated
    # rather than sorted in among the singles.
    d["grouped"] = d.dropped.str.contains(r"\+", regex=True)
    rank = {s: i for i, s in enumerate(names)}
    singles = (d[~d.grouped].assign(_o=lambda x: x.dropped.map(rank))
               .sort_values("_o").drop(columns="_o"))
    d = pd.concat([singles, d[d.grouped]], ignore_index=True)
    y = np.arange(len(d))[::-1]

    ax.axvline(full_fold, color=INK, ls="--", lw=LINEWIDTH["rule"],
               label=f"Full cohort ({full_fold:.2f})")
    # Single vs grouped is not a REGION distinction, so it must not borrow the
    # region hues — the same rule that keeps panel A off viridis. Both are drawn
    # in ink and separated by marker shape, which also survives greyscale.
    for grouped, marker, lab in ((False, "o", "One worm held out"),
                                 (True, "s", "Group held out")):
        sub = d[d.grouped == grouped]
        if sub.empty:
            continue
        ys = y[d.grouped.to_numpy() == grouped]
        ax.scatter(sub.fold_piN_piS, ys, s=12.0, marker=marker, color=INK,
                   edgecolor="black", linewidth=LINEWIDTH["hairline"],
                   zorder=3, label=lab)
    ax.set_yticks(y)
    ax.set_yticklabels([("+".join(short(p) for p in s.split("+"))
                         if "+" in s else short(s)) for s in d.dropped])
    # Plain text, as dN/dS: no mathtext italics (2026-10-09 rule).
    ax.set_xlabel("Eliminated / retained \u03c0N/\u03c0S")
    ax.set_ylabel("Worm or group held out")
    ax.tick_params(axis="y", length=1.5, pad=1.5)
    # Scale to the SPREAD of the folds, not to 1. This panel answers "does the
    # result move when a worm is dropped", and the answer lives in a range about
    # 0.05 wide; forcing 1.0 into view compresses every point into the last
    # tenth of the panel and shows nothing. The reference line is drawn only
    # when it is close enough to be informative rather than to be a margin.
    lo = float(min(d.fold_piN_piS.min(), full_fold))
    hi = float(max(d.fold_piN_piS.max(), full_fold))
    pad = max(0.30 * (hi - lo), 0.02)
    if lo - pad <= 1.0:
        ax.axvline(1.0, color=RULE, ls=":", lw=LINEWIDTH["rule"])
        lo = min(lo, 1.0)
    else:
        logging.info("D  x-axis is zoomed to the jackknife spread; the "
                     "no-difference line at 1.0 is off-scale (every value is "
                     "%.2f or above). Say so in the caption.", lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(-0.8, len(d) - 0.2)
    return {"full_fold": full_fold,
            "loo_min": float(d[~d.grouped].fold_piN_piS.min()),
            "loo_max": float(d[~d.grouped].fold_piN_piS.max()),
            "grouped": {r.dropped: float(r.fold_piN_piS)
                        for r in d[d.grouped].itertuples()}}


def fit_left_column(fig, axA, axB, axC, axD, tol_mm: float = FIT_TOL_MM,
                    max_pass: int = 6) -> dict:
    """Size the left column so panel A's square matrix fills it.

    A uses aspect="equal", so constrained_layout fits the largest square that
    its cell allows and centres it, leaving the difference empty: beside the
    matrix when the column is wider than the row is tall, above and below it
    when the column is narrower. The square fills its cell exactly when the
    left column's axes width equals row 1's axes height. Both are read off
    panels that fill their cells: C (left column, auto aspect) gives the
    column width, and B (row 1, auto aspect) gives the row height.

    The axes widths of the two columns follow the gridspec width ratios, so
    the share that makes them equal is row height / total axes width. The
    margins move a little when the share changes (A's colorbar is sized from
    its parent), hence the loop; it settles in one or two passes.

    Returns the share, the size of A's matrix, the passes taken and the
    remaining mismatch, for the log and the provenance file.
    """
    gs = axA.get_subplotspec().get_gridspec()
    ratios = gs.get_width_ratios() or [1.0, 1.0]
    share = ratios[0] / sum(ratios)
    w_mm, h_mm = (float(v) * MM_PER_IN for v in fig.get_size_inches())
    gap = float("nan")
    for n_pass in range(1, max_pass + 1):
        fig.canvas.draw()
        pos_b, pos_c, pos_d = (ax.get_position() for ax in (axB, axC, axD))
        row_h = pos_b.height * h_mm
        col_w = pos_c.width * w_mm
        axes_w = (pos_c.width + pos_d.width) * w_mm
        gap = col_w - row_h            # > 0: empty band beside A; < 0: A short
        logging.debug("fit pass %d: share %.4f, column %.2f mm, row %.2f mm",
                      n_pass, share, col_w, row_h)
        if abs(gap) <= tol_mm:
            break
        share = min(max(row_h / axes_w, 0.20), 0.60)
        gs.set_width_ratios([share, 1.0 - share])
    else:
        logging.warning("Left column did not settle in %d passes (A and its "
                        "column still differ by %.1f mm); the figure is "
                        "drawn, but check panel A.", max_pass, gap)
    pos_a = axA.get_position()
    return {"left_share": round(share, 4), "passes": n_pass,
            "a_matrix_mm": round(min(pos_a.width * w_mm, pos_a.height * h_mm), 1),
            "residual_mm": round(gap, 2)}


def place_letters(fig, panels: dict) -> None:
    """Bold panel letters: one height per row, at the left of each panel.

    Placed in figure coordinates after the layout is frozen, so they sit where
    the panels actually are rather than at a fraction of an axes box, which put
    B a different distance above its axes than A. Each row's letters share one
    top, the top of the tallest thing in the row (B's legend in row 1; the C
    and D legends in row 2). A and C share the figure's left edge, as in the
    hand layout. B and D each sit at the left of their own labels: aligning
    them would carry B's letter out into the gap beside A's colorbar, because
    D's held-out group names are much wider than B's tick labels.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    to_fig = fig.transFigure.inverted()
    box = {k: ax.get_tightbbox(rend) for k, ax in panels.items()}
    x = {k: box[k].x0 for k in panels}
    x["A"] = x["C"] = min(box["A"].x0, box["C"].x0)
    for row in ("AB", "CD"):
        top = max(box[k].y1 for k in row)
        for k in row:
            fx, fy = to_fig.transform((x[k], top))
            fig.text(fx, fy, k, fontsize=FONT["panel_letter"],
                     fontweight="bold", ha="left", va="top")


def outside_canvas(fig, tol_mm: float = 0.2) -> list[str]:
    """Labels, tick labels, legends and letters that run past the canvas edge.

    audit() checks type size, width and legend-on-data overlap but not this,
    and the equal-column build put D's legend off the right edge while audit()
    reported the figure clean. Axis extents come from Axis.get_tightbbox, which
    counts only the tick labels actually drawn.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    px_per_mm = fig.dpi / MM_PER_IN
    w_px, h_px = fig.bbox.width, fig.bbox.height
    items = []
    for ax in fig.axes:
        name = ax.get_label() or "axes"
        for axis in (ax.xaxis, ax.yaxis):
            bb = axis.get_tightbbox(rend)
            if bb is not None and bb.width > 0:
                items.append((f"{name} {axis.axis_name}-axis labels", bb))
        leg = ax.get_legend()
        if leg is not None and leg.get_visible():
            items.append((f"{name} legend", leg.get_window_extent(rend)))
        for t in ax.texts:
            if t.get_visible() and t.get_text().strip():
                items.append((f"{name} text '{t.get_text()[:20]}'",
                              t.get_window_extent(rend)))
    for t in fig.texts:
        if t.get_visible() and t.get_text().strip():
            items.append((f"letter {t.get_text()}", t.get_window_extent(rend)))
    found = []
    for what, bb in items:
        over = max(-bb.x0, -bb.y0, bb.x1 - w_px, bb.y1 - h_px) / px_per_mm
        if over > tol_mm:
            found.append(f"{what} run {over:.1f} mm past the edge of the canvas")
    return found


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Supplementary fig. S17 (stem S19) — cohort integrity.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--pairwise", required=True, type=Path,
                   help="pairwise.tsv from check_duplicate_individuals.py")
    p.add_argument("--het", required=True, type=Path,
                   help="het_by_region.tsv from het_by_region.py")
    p.add_argument("--jackknife", required=True, type=Path,
                   help="jackknife_summary.tsv from jackknife_individuals.py")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--dup-discordance", type=float, default=0.05,
                   help="Match check_duplicate_individuals.py, for the count "
                        "of flagged pairs reported in the log")
    p.add_argument("--ratio-tolerance", type=float, default=0.25,
                   help="Match het_by_region.py; sets the shaded band in C")
    p.add_argument("--tag", default=None)
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    for noisy in ("fontTools", "matplotlib.font_manager", "matplotlib"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    apply_style()

    outdir = Path(make_outdir("s19_cohort_integrity",
                              base=str(args.outdir.expanduser().resolve()),
                              tag=args.tag))
    logging.info("Output: %s", outdir)

    pw = pd.read_csv(args.pairwise.expanduser(), sep="\t")
    het = pd.read_csv(args.het.expanduser(), sep="\t")
    jk = pd.read_csv(args.jackknife.expanduser(), sep="\t")

    # The three tables must describe the same cohort. They come from separate
    # runs and separate output directories, so mixing a stale one in is the
    # obvious failure and it would be invisible in the drawn figure.
    s_pw = set(pw.sample_a) | set(pw.sample_b)
    s_het = set(het["sample"])
    s_jk = {s for s in jk.dropped if s not in ("(none)",) and "+" not in s}
    if not (s_pw == s_het == s_jk):
        logging.error("The three inputs describe different cohorts:\n"
                      "  pairwise  %d: %s\n  het       %d: %s\n"
                      "  jackknife %d: %s",
                      len(s_pw), sorted(s_pw), len(s_het), sorted(s_het),
                      len(s_jk), sorted(s_jk))
        logging.error("Re-run the missing one against the same callset rather "
                      "than assembling a figure from mismatched runs.")
        return 1
    logging.info("Cohort agrees across all three inputs: %d worms", len(s_pw))

    fig, axes = figure(FIG_W, height_mm=FIG_H_MM, nrows=2, ncols=2,
                       gridspec_kw={"width_ratios": [LEFT_SHARE_START,
                                                     1.0 - LEFT_SHARE_START]})
    (axA, axB), (axC, axD) = axes
    for ax, name in ((axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")):
        ax.set_label(name)              # names the panel in outside_canvas()
    order = sorted(s_pw)
    a_info = panel_a(axA, pw, args.dup_discordance, order)
    bc_info = panel_bc(axB, axC, het, args.ratio_tolerance, order)
    d_info = panel_d(axD, jk, order)

    # Every legend goes ABOVE its axes. Placing them in a corner means guessing
    # which corner the data leaves free, and that guess depends on the data:
    # ordering B and C by accession moved the tall bars to the left, and zooming
    # D's x-axis filled the panel. Both guesses were wrong on the real data and
    # audit() caught both. Above the axes, an overlap cannot happen.
    # Panel C builds its own legend in panel_bc, because its key is a
    # composite swatch; it uses the same LEGEND_ABOVE settings.
    # D's legend ends at D's right edge and runs left over the column of
    # held-out names, as in the hand layout of 2026-10-30. Started at D's left
    # spine it is wider than D, and constrained_layout answers that by taking
    # a right margin from the whole column, which shrinks B and D and leaves
    # the left column too narrow for A's matrix (it stopped A at 47 of its
    # 51.5 mm in testing).
    if axB.get_legend_handles_labels()[0]:
        axB.legend(ncol=2, **LEGEND_ABOVE)
    if axD.get_legend_handles_labels()[0]:
        axD.legend(ncol=3, **dict(LEGEND_ABOVE, loc="lower right",
                                  bbox_to_anchor=(1.0, 1.02)))

    # Size the left column to A's matrix, then freeze the solved layout so the
    # letters (placed in figure coordinates) cannot move the panels.
    layout = fit_left_column(fig, axA, axB, axC, axD)
    fig.set_layout_engine("none")
    logging.info("Layout: 185 mm; A's matrix %.1f mm square, left column "
                 "%.1f%% of the axes width (%d pass(es), residual %.2f mm)",
                 layout["a_matrix_mm"], 100 * layout["left_share"],
                 layout["passes"], layout["residual_mm"])
    place_letters(fig, {"A": axA, "B": axB, "C": axC, "D": axD})

    problems = audit(fig)
    overflow = outside_canvas(fig)
    for item in overflow:
        logging.warning("Off the canvas: %s", item)
    problems += overflow
    if problems:
        logging.warning("S17 breaks %d journal rule(s) — listed above.",
                        len(problems))
    written = pc_save(fig, str(outdir), "S19_cohort_integrity")
    plt.close(fig)

    (outdir / "provenance.json").write_text(json.dumps({
        "argv": sys.argv,
        "inputs": {"pairwise": str(args.pairwise), "het": str(args.het),
                   "jackknife": str(args.jackknife)},
        "panel_a": a_info, "panel_bc": bc_info, "panel_d": d_info,
        "layout": layout,
        "audit_problems": problems,
    }, indent=2, default=str))

    logging.info("=" * 74)
    logging.info("A  closest pair disagrees at %.1f%% of genotypes; %d pair(s) "
                 "below the duplicate threshold",
                 100 * a_info["min_discordance"], a_info["flagged"])
    logging.info("B  heterozygosity spans %.1fx across the cohort",
                 bc_info["het_spread_fold"])
    logging.info("C  cohort median eliminated:retained ratio %.3f; %d worm(s) "
                 "outside the +/-%.0f%% band%s",
                 bc_info["cohort_median_ratio"], len(bc_info["outside_band"]),
                 100 * args.ratio_tolerance,
                 (": " + ", ".join(bc_info["outside_band"]))
                 if bc_info["outside_band"] else "")
    if d_info:
        logging.info("D  leave-one-out %.3f-%.3f against %.3f for the full "
                     "cohort", d_info["loo_min"], d_info["loo_max"],
                     d_info["full_fold"])
        for k, v in d_info.get("grouped", {}).items():
            logging.info("   grouped drop %-28s %.3f", k, v)
    logging.info("")
    logging.info("Panel C is the one that carries B. A 2.5x spread in absolute "
                 "heterozygosity looks alarming on its own; what decides "
                 "whether it can bias a between-class contrast is whether the "
                 "RATIO moves, and that is panel C.")
    logging.info("Figure: %s", ", ".join(Path(w).name for w in written))
    logging.info("Export from the SVG — Arial is absent on Isaac nodes.")
    logging.info("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
