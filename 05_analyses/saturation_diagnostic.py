#!/usr/bin/env python3
"""
saturation_diagnostic.py - is the Toxocara signal an artefact of dS saturation?

Supplementary figure S11 (S12 before the 2026-10-12 renumbering; the output
stem saturation_diagnostic is unchanged). The manuscript uses *Toxocara* as the
outgroup and
reports results from comparisons whose median dS is 1.39-1.43, against a rough
reliability ceiling of dS ~ 1.5 for pairwise yn00 estimates. That invites one
obvious objection, and this figure answers it before it is raised.

Three panels:

    A   cumulative dS per comparison, with the 1.0 and 1.5 marks drawn, so the
        fraction of pairs past each threshold can be read off directly
    B   median dN/dS against dS across all pairs pooled, which is where a
        saturation artefact shows itself: as dS estimates inflate near
        saturation, omega is pushed down, so a declining trend at high dS is the
        signature rather than a biological finding. Above the axes: the
        within-comparison Spearman rho of each group of pairs (the pooled rho
        goes to the summary, since it cannot be read on its own; see below)
    C   the retained-versus-eliminated dN/dS fold change in *Toxocara*
        recomputed at successively stricter dS cut-offs, with bootstrap
        intervals, the comparison named above the axes (--annotated-label) and
        the eliminated n over each interval. A flat line means the contrast
        does not depend on the saturated tail

**Why the diagnostic is needed at all.** Absolute dN/dS is not comparable
between comparisons of different divergence depth: the *Toxocara* pairs sit near
0.10 and the ascarid-ascarid pairs near 0.18-0.21, which read side by side as
*Toxocara* genes being the more constrained. Panel B is where that is tested,
and it has to be read carefully.

**The pooled correlation in panel B cannot be interpreted on its own.** The six
comparisons cluster in dS -- the ascarid pairs near 0.2-0.37, the Toxocara pairs
near 1.4, with little in between -- so a pooled negative rho may be nothing more
than two clusters at different heights, which is comparison identity rather than
any effect of dS. On synthetic data built with omega deliberately independent of
dS *inside* every comparison, the pooled rho still came out at -0.45. The
per-comparison correlations, written to `omega_vs_ds_within.tsv` and printed in
the summary, are what decide it: if omega also falls with dS within the Toxocara
comparisons, which span the saturation threshold internally, the effect is real;
if those are flat, the pooled trend is between-comparison structure and says
nothing about saturation.

**What a null result looks like.** Panel C is the one that carries the argument.
If the fold change is stable as the cut-off tightens, the eliminated-gene result
does not rest on saturated pairs. Watch the intervals as well as the point
estimates: tightening the cut-off discards genes, so a flat line with intervals
that widen to uselessness is not evidence of stability. The n at each cut-off is
printed on the panel for that reason.

Colours in A and B are the same sequential ramp used for the six comparisons in
Figure 1C, assigned in ascending median dS, so a comparison keeps its colour
across both figures.

Figure fixes, 2026-10-12 (figure only; every table is unchanged): full pair
names with × in the legend, in two rows under the x labels (ascarid pairs,
then outgroup pairs) with the blank band above it gone; panel B's statistics
above the axes with true minus signs, so nothing crosses the dS = 1.5 rule;
panel C's markers in ink (a fold change is a metric, not a region), "n = 333"
labels and the comparison named; panel letters by default (--no-panel-letters);
the panels laid out from their measured text. The summary ends with legend
notes for what each mark means.

Usage
-----
    python saturation_diagnostic.py \\
      --as-pu   kaks_results/kaks_filtered.tsv \\
      --as-bsc  kaks_As_Bsc_v2/kaks_filtered.tsv \\
      --pu-bsc  kaks_Pu_Bsc_v2/kaks_filtered.tsv \\
      --tox-as  kaks_Tox_As/kaks_filtered.tsv \\
      --tox-pu  kaks_Tox_Pu/kaks_filtered.tsv \\
      --tox-bsc kaks_Tox_Bsc_v2/kaks_filtered.tsv \\
      --annotated kaks_Tox_As_annotated/kaks_annotated.tsv \\
      --region-col tox_region

Panels A and B need only the six pairwise tables. Panel C needs --annotated and
is skipped without it.

Requires numpy, pandas, scipy and matplotlib, plus paper_colors.py and
rate_tables.py, which are located automatically in lib/.
"""

import argparse
import hashlib
import json
import logging
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation
from scipy import stats

_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, 'lib'), os.path.join(os.path.dirname(_HERE), 'lib')):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, _cand)

from paper_colors import (COL, FONT, LINEWIDTH, NEUTRAL, apply_style, audit,
                          mm, pair_ramp, save)
from rate_tables import read_rate_table

apply_style()

EXTS = ('.svg', '.pdf', '.png')

#: Thresholds reported in the summary. 1.5 is the operative one: pairwise dS
#: beyond it is generally treated as unreliable, and the Toxocara comparisons
#: sit just under it.
THRESHOLDS = (1.0, 1.5, 2.0)


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return 'n/a'
    if p <= 0:
        # Underflowed to exactly zero. log10(0) is -inf and int() of that
        # raises, so a large-n test with a strong effect would crash the run
        # at the formatting step. Report the floor honestly instead.
        return '< 10^-300'
    if p >= 1e-3:
        return f'{p:.3f}'
    exp = int(np.floor(np.log10(p)))
    return f'{p / 10 ** exp:.2f} x 10^{exp}'


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def exceedance(ds, thresholds=THRESHOLDS):
    """Fraction of pairs above each dS threshold."""
    ds = np.asarray(ds, dtype=float)
    ds = ds[np.isfinite(ds)]
    return {t: float((ds > t).mean()) for t in thresholds}


def omega_vs_ds(df_all, n_bins=18, min_per_bin=50):
    """Median dN/dS in dS bins, pooled over every comparison.

    Pooled rather than per-comparison on purpose: the question is whether omega
    depends on dS at all, and each comparison alone spans too narrow a dS range
    to answer it. Bins are quantiles of the pooled dS so every bin carries
    comparable weight.
    """
    d = df_all.dropna(subset=['dS', 'dN_dS'])
    d = d[(d['dS'] > 0) & np.isfinite(d['dN_dS'])]
    if len(d) < n_bins * min_per_bin:
        n_bins = max(4, len(d) // max(min_per_bin, 1))
    edges = np.unique(np.quantile(d['dS'], np.linspace(0, 1, n_bins + 1)))
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = d[(d['dS'] >= lo) & (d['dS'] < hi)]
        if len(sel) < min_per_bin:
            continue
        rows.append({'ds_mid': float(np.median(sel['dS'])),
                     'omega_med': float(np.median(sel['dN_dS'])),
                     'omega_q1': float(sel['dN_dS'].quantile(0.25)),
                     'omega_q3': float(sel['dN_dS'].quantile(0.75)),
                     'n': int(len(sel))})
    out = pd.DataFrame(rows)
    # Spearman on the per-gene values, not on the binned medians: binning first
    # would inflate the correlation by discarding within-bin scatter.
    rho, p = stats.spearmanr(d['dS'], d['dN_dS'])
    return out, float(rho), float(p)


def omega_ds_within(per_pair):
    """Spearman of dN/dS on dS *inside* each comparison.

    The pooled correlation is not sufficient on its own. Comparisons cluster in
    dS -- the ascarid pairs near 0.2-0.37, the Toxocara pairs near 1.4, with
    little between -- so a pooled negative rho could be nothing more than two
    clusters at different heights, which is comparison identity rather than any
    effect of dS. If omega also falls with dS *within* a comparison, where every
    pair shares one divergence and one alignment pipeline, that reading is
    closed off. The Toxocara comparisons are the informative ones: they span the
    saturation threshold internally.
    """
    rows = []
    for name, df in per_pair.items():
        d = df.dropna(subset=['dS', 'dN_dS'])
        d = d[(d['dS'] > 0) & np.isfinite(d['dN_dS'])]
        if len(d) < 100:
            continue
        rho, p = stats.spearmanr(d['dS'], d['dN_dS'])
        rows.append({'comparison': name, 'n': int(len(d)),
                     'median_dS': float(d['dS'].median()),
                     'rho': float(rho), 'p': float(p)})
    return pd.DataFrame(rows)


def fold_by_cutoff(df, region_col, cutoffs, rng, n_boot=2000):
    """Retained-vs-eliminated dN/dS fold change at successive dS cut-offs.

    The interval is a percentile bootstrap over genes, resampled within each
    region so the group sizes are preserved.
    """
    rows, skipped = [], []
    d0 = df.dropna(subset=['dS', 'dN_dS', region_col])
    for cut in cutoffs:
        d = d0 if cut is None else d0[d0['dS'] <= cut]
        ret = d.loc[d[region_col] == 'retained', 'dN_dS'].to_numpy()
        eli = d.loc[d[region_col] == 'eliminated', 'dN_dS'].to_numpy()
        if len(ret) < 20 or len(eli) < 20:
            logging.warning('dS cutoff %s leaves %d retained / %d eliminated -- '
                            'skipped', cut, len(ret), len(eli))
            skipped.append((cut, len(ret), len(eli)))
            continue
        fold = float(np.median(eli) / np.median(ret))
        boot = np.empty(n_boot)
        for b in range(n_boot):
            r = ret[rng.integers(0, len(ret), len(ret))]
            e = eli[rng.integers(0, len(eli), len(eli))]
            boot[b] = np.median(e) / np.median(r)
        _, p = stats.mannwhitneyu(ret, eli, alternative='two-sided')
        rows.append({'cutoff': cut, 'fold': fold,
                     'lo': float(np.percentile(boot, 2.5)),
                     'hi': float(np.percentile(boot, 97.5)),
                     'n_ret': int(len(ret)), 'n_elim': int(len(eli)),
                     'med_ret': float(np.median(ret)),
                     'med_elim': float(np.median(eli)), 'p': float(p)})
    return pd.DataFrame(rows), skipped


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

#: Display names of the six comparisons (the internal keys, 'As x Pu' etc.,
#: stay in the tables). Species-only text joined by ×, so house_text()
#: italicises the names.
DISPLAY = {'As x Pu': 'Ascaris × P. univalens',
           'As x Bsc': 'Ascaris × B. schroederi',
           'Pu x Bsc': 'P. univalens × B. schroederi',
           'Tox x As': 'T. canis × Ascaris',
           'Tox x Pu': 'T. canis × P. univalens',
           'Tox x Bsc': 'T. canis × B. schroederi'}

MINUS = '−'


def signed(v, nd=2):
    """'+0.06' / '−0.24', with a true minus."""
    return f'{v:+.{nd}f}'.replace('-', MINUS)


def rho_range(within, outgroup_prefix='Tox'):
    """Within-comparison Spearman rho, as a range per group of pairs:
    {'ascarid': (lo, hi, n), 'outgroup': (lo, hi, n)}; a group with no
    comparison (fewer than 100 usable pairs each) is left out."""
    out = {}
    if within is None or not len(within):
        return out
    is_out = within['comparison'].str.startswith(outgroup_prefix)
    for name, sel in (('ascarid', within[~is_out]), ('outgroup', within[is_out])):
        if len(sel):
            out[name] = (float(sel['rho'].min()), float(sel['rho'].max()),
                         int(len(sel)))
    return out


def range_text(lo, hi):
    return signed(lo) if abs(hi - lo) < 0.005 else f'{signed(lo)} to {signed(hi)}'


def drawn_texts(fig):
    """Every text object actually drawn (tick labels outside the view limits
    are returned by findobj but never drawn, so they are excluded)."""
    out = list(fig.texts)
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            tol = 1e-9 * max(abs(hi - lo), 1e-300)
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                loc = tick.get_loc()
                if loc is None or not (lo - tol <= loc <= hi + tol):
                    continue
                out += [lab for lab in (tick.label1, tick.label2)
                        if lab.get_visible()]
            out.append(axis.label)
            out.append(axis.get_offset_text())
        out += list(ax.texts)
    for leg in fig.legends:
        out += list(leg.get_texts())
    return [t for t in out if t.get_visible() and (t.get_text() or '').strip()]


def text_problems(fig, tol_px=1.0):
    """Drawn text running past the canvas edge, or overlapping other drawn
    text. Neither is checked by audit(), and both save without complaint."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    out = []
    texts = drawn_texts(fig)
    boxes = [t.get_window_extent(rend) for t in texts]
    for t, bb in zip(texts, boxes):
        if (bb.x0 < -tol_px or bb.y0 < -tol_px or bb.x1 > W + tol_px
                or bb.y1 > H + tol_px):
            out.append(f'text {t.get_text()[:28]!r} runs past the canvas edge')
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = boxes[i], boxes[j]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) > tol_px and
                    min(a.y1, b.y1) - max(a.y0, b.y0) > tol_px):
                out.append(f'text {texts[i].get_text()[:20]!r} overlaps '
                           f'{texts[j].get_text()[:20]!r}')
    return out


def place_row(fig, axes, shares, gaps_mm, pad_mm=0.8, iterations=4,
              min_top_px=0.0, bottom_extra_px=0.0):
    """Lay out panels side by side from their measured text, so the row fills
    the canvas whatever font the node renders with. At least `min_top_px` is
    kept above the axes (the panel-letter row) and `bottom_extra_px` below
    everything the axes draw (the legend)."""
    W, H = fig.bbox.width, fig.bbox.height
    pad = mm(pad_mm) * fig.dpi
    gaps = [mm(g) * fig.dpi for g in gaps_mm]
    shares = np.asarray(shares, float) / np.sum(shares)
    for _ in range(iterations):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        tight = [ax.get_tightbbox(r) for ax in axes]
        boxes = [ax.get_window_extent(r) for ax in axes]
        lefts = [max(0.0, b.x0 - t.x0) for t, b in zip(tight, boxes)]
        rights = [max(0.0, t.x1 - b.x1) for t, b in zip(tight, boxes)]
        top = max(max(0.0, t.y1 - b.y1) for t, b in zip(tight, boxes))
        bottom = max(max(0.0, b.y0 - t.y0) for t, b in zip(tight, boxes))
        top, bottom = max(top, min_top_px), bottom + bottom_extra_px
        avail = W - 2 * pad - sum(lefts) - sum(rights) - sum(gaps)
        widths = avail * shares
        x = pad
        y0, hh = pad + bottom, H - 2 * pad - bottom - top
        for i, ax in enumerate(axes):
            x += lefts[i]
            ax.set_position([x / W, y0 / H, widths[i] / W, hh / H])
            x += widths[i] + rights[i] + (gaps[i] if i < len(gaps) else 0.0)


def above_axes(ax, s, x=0.0, ha='left', dy_pt=3.0, **kw):
    """Text just above the axes, at axes fraction `x`."""
    tr = ax.transAxes + ScaledTranslation(0, dy_pt / 72,
                                          ax.figure.dpi_scale_trans)
    opts = dict(ha=ha, va='bottom', fontsize=FONT['annotation'],
                color=NEUTRAL['ink'], clip_on=False, linespacing=1.25)
    opts.update(kw)
    return ax.text(x, 1.0, s, transform=tr, **opts)


def keep_inside_x(ax, texts, gap_pt=1.5, iterations=3):
    """Widen the x-limits symmetrically until every text in `texts` sits
    inside the axes; the n labels over the first and last intervals would
    otherwise cross the spines. Returns how far the limits were widened, in
    data units per side."""
    fig = ax.figure
    gap = gap_pt / 72 * fig.dpi
    total = 0.0
    for _ in range(iterations):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        bb = ax.get_window_extent(r)
        ext = [t.get_window_extent(r) for t in texts]
        over = max([bb.x0 + gap - e.x0 for e in ext]
                   + [e.x1 - (bb.x1 - gap) for e in ext] + [0.0])
        if over <= 0.5:
            break
        lo, hi = ax.get_xlim()
        d = over * (hi - lo) / bb.width
        ax.set_xlim(lo - d, hi + d)
        total += d
    return total


def n_labels_fit(ax, texts, gap_pt=1.5):
    """True if no two of the n labels touch (they sit at one height)."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes = sorted((t.get_window_extent(r) for t in texts), key=lambda b: b.x0)
    gap = gap_pt / 72 * fig.dpi
    return all(b.x0 - a.x1 >= gap for a, b in zip(boxes, boxes[1:]))


def plot_saturation(per_pair, order, ramp, binned, within, folds, outdir,
                    stem, width='double', height_mm=62.0, letters=True,
                    ds_view=3.0, annotated_label=None):
    """Supplementary fig. S11 (S12 before renumbering).

    A cumulative dS per comparison; B median dN/dS per pooled dS bin (line),
    its interquartile range (band) and each comparison's median (points),
    with the within-comparison Spearman rho of each group of pairs above the
    axes (the pooled rho cannot be read on its own -- see the module notes --
    and goes to the summary); C the eliminated/retained dN/dS fold change by
    dS cut-off for the --annotated comparison, named above the axes, with the
    eliminated n over each interval.

    The panels are laid out from their measured text (place_row), the legend
    of the six comparisons sits under the x labels in two rows (ascarid pairs,
    then outgroup pairs), and the panel letters sit level with each panel's
    left edge. Returns (paths, problems, facts for the summary).
    """
    n_panels = 3 if folds is not None and len(folds) else 2
    fig = plt.figure(figsize=(COL[width], mm(height_mm)))
    fig.set_layout_engine('none')
    axes = [fig.add_axes([0.08 + 0.32 * i, 0.3, 0.22, 0.5])
            for i in range(n_panels)]
    facts = {}

    # --- A: cumulative dS -------------------------------------------------
    ax = axes[0]
    for k, colour in zip(order, ramp):
        ds = np.sort(per_pair[k]['dS'].dropna().to_numpy())
        ax.plot(ds, np.arange(1, len(ds) + 1) / len(ds), color=colour,
                linewidth=LINEWIDTH['data'], zorder=3)
    for t in (1.0, 1.5):
        ax.axvline(t, color=NEUTRAL['rule'], linestyle='--',
                   linewidth=LINEWIDTH['thin'], zorder=2)
        ax.text(t, 1.02, f'{t:g}', ha='center', va='bottom',
                fontsize=FONT['annotation'], color=NEUTRAL['ink_light'])
    ax.set_xlim(0, ds_view)
    ax.set_ylim(0, 1)
    ax.set_xlabel('dS')
    ax.set_ylabel('Cumulative fraction of pairs')
    facts['A_beyond_view'] = {
        DISPLAY[k]: int((per_pair[k]['dS'] > ds_view).sum()) for k in order}

    # --- B: omega against dS ----------------------------------------------
    ax = axes[1]
    ax.fill_between(binned['ds_mid'], binned['omega_q1'], binned['omega_q3'],
                    color=NEUTRAL['grid'], zorder=2, linewidth=0)
    ax.plot(binned['ds_mid'], binned['omega_med'], color=NEUTRAL['ink'],
            linewidth=LINEWIDTH['data'], zorder=4)
    # Each comparison at its own median dS, so the ranking in the diagnostic
    # boxplot can be read straight off the trend.
    for k, colour in zip(order, ramp):
        d = per_pair[k].dropna(subset=['dS', 'dN_dS'])
        ax.plot([d['dS'].median()], [d['dN_dS'].median()], marker='o',
                markersize=3.4, markerfacecolor=colour,
                markeredgecolor=NEUTRAL['ink'], markeredgewidth=0.7, zorder=5)
    ax.axvline(1.5, color=NEUTRAL['rule'], linestyle='--',
               linewidth=LINEWIDTH['thin'], zorder=3)
    ax.set_xlim(0, ds_view)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('dS')
    ax.set_ylabel('dN/dS')
    rr = rho_range(within)
    lines = [f'{name.capitalize()} pairs: ρ = {range_text(lo, hi)}'
             for name, (lo, hi, _) in rr.items()]
    if lines:
        above_axes(ax, '\n'.join(lines))
    facts['B_rho_text'] = lines

    # --- C: fold change against cut-off -----------------------------------
    n_texts = []
    if n_panels == 3:
        ax = axes[2]
        x = np.arange(len(folds))
        ax.errorbar(x, folds['fold'],
                    yerr=[folds['fold'] - folds['lo'], folds['hi'] - folds['fold']],
                    fmt='o', markersize=3.4, color=NEUTRAL['ink'],
                    markerfacecolor=NEUTRAL['ink'],
                    markeredgecolor=NEUTRAL['ink'], markeredgewidth=0.7,
                    elinewidth=LINEWIDTH['thin'], capsize=1.8,
                    capthick=LINEWIDTH['thin'], zorder=4)
        ax.axhline(1.0, color=NEUTRAL['rule'], linestyle='--',
                   linewidth=LINEWIDTH['rule'], zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels(['All' if pd.isna(c) or c is None else f'≤{c:g}'
                            for c in folds['cutoff']])
        ax.set_xlim(-0.5, len(folds) - 0.5)
        ax.set_xlabel('dS cut-off')
        ax.set_ylabel('dN/dS fold change\n(eliminated / retained)')
        lo = min(0.95, folds['lo'].min() * 0.97)
        ax.set_ylim(lo, folds['hi'].max() * 1.16)
        for xi, r in zip(x, folds.itertuples()):
            n_texts.append(ax.text(
                xi, r.hi + (ax.get_ylim()[1] - lo) * 0.03, f'n = {r.n_elim:,}',
                ha='center', va='bottom', fontsize=FONT['annotation'],
                color=NEUTRAL['ink_light']))
        if annotated_label:
            above_axes(ax, annotated_label)
        facts['C_label'] = annotated_label or ''

    # --- legend: two rows, ascarid pairs then outgroup pairs ---------------
    handles = [Line2D([], [], color=c, linewidth=LINEWIDTH['data'])
               for c in ramp]
    names = [DISPLAY[k] for k in order]
    ncol = 3 if len(order) == 6 else len(order)
    if ncol == 3:            # row-major: matplotlib fills legends by column
        idx = [0, 3, 1, 4, 2, 5]
        handles, names = [handles[i] for i in idx], [names[i] for i in idx]
    leg = fig.legend(handles=handles, labels=names, loc='lower center',
                     bbox_to_anchor=(0.5, mm(0.8) / mm(height_mm)), ncol=ncol,
                     frameon=False, fontsize=FONT['legend'], handlelength=1.5,
                     columnspacing=1.6, handletextpad=0.5, borderaxespad=0.0,
                     borderpad=0.0)

    # --- layout -------------------------------------------------------------
    letter_texts = ([fig.text(0, 1, s, ha='left', va='top',
                              fontsize=FONT['panel_letter'], fontweight='bold',
                              color=NEUTRAL['ink'])
                     for s in 'ABC'[:n_panels]] if letters else [])
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    leg_h = leg.get_window_extent(r).height
    letter_h = max((t.get_window_extent(r).height for t in letter_texts),
                   default=0.0)
    place_row(fig, axes, shares=[1.0] * n_panels,
              gaps_mm=[8.0] * (n_panels - 1),
              min_top_px=letter_h + mm(0.8) * fig.dpi if letters else 0.0,
              bottom_extra_px=leg_h + mm(2.0) * fig.dpi)
    # n labels: "n = 333" over every interval when they fit; else "n =" on
    # the first only; else bare numbers (the legend notes then define them).
    # Each form is tried at the panel's own x-limits, widened just enough to
    # keep the first and last labels inside the spines.
    if n_texts:
        n_elim = [int(v) for v in folds['n_elim']]
        forms = (('"n = 333" over every interval',
                  [f'n = {v:,}' for v in n_elim]),
                 ('"n =" on the first interval only',
                  [f'n = {n_elim[0]:,}'] + [f'{v:,}' for v in n_elim[1:]]),
                 ('bare numbers over the intervals; "n =" did not fit',
                  [f'{v:,}' for v in n_elim]))
        for desc, texts in forms:
            for t, s_ in zip(n_texts, texts):
                t.set_text(s_)
            axes[2].set_xlim(-0.5, len(folds) - 0.5)
            widened = keep_inside_x(axes[2], n_texts)
            if n_labels_fit(axes[2], n_texts):
                break
        facts['C_n_labels'] = desc
        facts['C_xlim_widened'] = round(widened, 3)
        logging.info('  C: n labels %s; x-limits widened by %.3f a side', desc,
                     widened)
    if letters:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        W, H = fig.bbox.width, fig.bbox.height
        pad = mm(0.8) * fig.dpi
        for ax, t in zip(axes, letter_texts):
            t.set_position((max(ax.get_tightbbox(r).x0, pad) / W,
                            1.0 - pad / H))
    problems = list(audit(fig, verbose=False)) + text_problems(fig)
    for problem in problems:
        logging.warning(f'FIGURE AUDIT [{stem}]: {problem}')
    paths = save(fig, outdir, stem, formats=[e.lstrip('.') for e in EXTS])
    plt.close(fig)
    logging.info(f"Saved: {stem}{'/'.join(EXTS)}")
    return paths, problems, facts


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def fingerprint(path):
    st = os.stat(path)
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return {'path': os.path.abspath(path), 'bytes': st.st_size,
            'mtime_utc': datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
            'sha256_16': h.hexdigest()[:16]}


def git_commit():
    try:
        out = subprocess.run(['git', '-C', _HERE, 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            dirty = subprocess.run(['git', '-C', _HERE, 'status', '--porcelain'],
                                   capture_output=True, text=True, timeout=5)
            return out.stdout.strip() + ('+dirty' if dirty.stdout.strip() else '')
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='dS saturation diagnostic for the Toxocara comparisons '
                    '(supplementary fig. S11; S12 before renumbering).')
    for flag, helptext in (('--as-pu', 'Ascaris vs Parascaris'),
                           ('--as-bsc', 'Ascaris vs Baylisascaris (use _v2)'),
                           ('--pu-bsc', 'Parascaris vs Baylisascaris (use _v2)'),
                           ('--tox-as', 'Toxocara vs Ascaris'),
                           ('--tox-pu', 'Toxocara vs Parascaris'),
                           ('--tox-bsc', 'Toxocara vs Baylisascaris (use _v2)')):
        p.add_argument(flag, required=True, help=f'{helptext} rate table')
    p.add_argument('--annotated', default=None,
                   help='Annotated Toxocara table carrying region calls. Panel C '
                        'is skipped without it.')
    p.add_argument('--region-col', default='tox_region',
                   help="Region column in --annotated (default: tox_region -- "
                        "Toxocara's own calls, not projected from Ascaris)")
    p.add_argument('--cutoffs', type=float, nargs='+',
                   default=[2.0, 1.5, 1.2, 1.0, 0.8],
                   help='dS cut-offs for panel C, applied in addition to the '
                        'no-cut-off case')
    p.add_argument('--ds-view', type=float, default=3.0,
                   help='Upper dS limit of panels A and B (display only)')
    p.add_argument('--n-boot', type=int, default=2000)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--width', choices=['single', 'onehalf', 'double'],
                   default='double')
    p.add_argument('--height-mm', type=float, default=62.0)
    p.add_argument('--no-panel-letters', action='store_true',
                   help='Omit the A-C panel letters (drawn by default, since '
                        'the legend refers to the panels)')
    p.add_argument('--panel-letters', action='store_true',
                   help=argparse.SUPPRESS)   # old opt-in; letters are default
    p.add_argument('--annotated-label', default='T. canis × Ascaris',
                   help='Name of the --annotated comparison, shown above panel '
                        'C (default "T. canis × Ascaris"; empty string for none)')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'saturation_diagnostic_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    # mode='w': FileHandler appends by default, so rerunning into an
    # existing directory concatenates logs while the figures and
    # provenance.json are overwritten -- the log then describes more
    # runs than the directory contains.
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(os.path.join(args.outdir, 'saturation.log'),
                                     mode='w'),
                  logging.StreamHandler(sys.stdout)])
    for noisy in ('matplotlib.font_manager', 'fontTools', 'PIL'):
        logging.getLogger(noisy).setLevel(logging.ERROR)

    logging.info('=' * 60)
    logging.info('dS saturation diagnostic')
    logging.info('=' * 60)

    spec = [('As x Pu', args.as_pu), ('As x Bsc', args.as_bsc),
            ('Pu x Bsc', args.pu_bsc), ('Tox x As', args.tox_as),
            ('Tox x Pu', args.tox_pu), ('Tox x Bsc', args.tox_bsc)]
    for _, path in spec:
        if not os.path.isfile(path):
            sys.exit(f'ERROR: input not found: {path}')

    rng = np.random.default_rng(args.seed)
    per_pair, labels = {}, {}
    for name, path in spec:
        df = read_rate_table(path, require=('dS', 'dN_dS'))
        per_pair[name] = df
        labels[name] = name
        logging.info('%-9s %7s pairs, median dS %.4f', name, f'{len(df):,}',
                     df['dS'].median())

    # Ramp assigned in ascending median dS, matching Figure 1C exactly.
    order = sorted(per_pair, key=lambda k: per_pair[k]['dS'].median())
    ramp = pair_ramp(len(order))

    # --- A: exceedance -----------------------------------------------------
    logging.info('\n--- fraction of pairs above each dS threshold ---')
    exc_rows = []
    for k in order:
        e = exceedance(per_pair[k]['dS'])
        exc_rows.append({'comparison': k, 'n': len(per_pair[k]),
                         'median_dS': float(per_pair[k]['dS'].median()),
                         **{f'frac_gt_{t:g}': e[t] for t in THRESHOLDS}})
        logging.info('  %-9s median dS %.4f   ' + '   '.join(
            f'>{t:g}: %5.1f%%' for t in THRESHOLDS), k,
            per_pair[k]['dS'].median(), *[100 * e[t] for t in THRESHOLDS])
    exc = pd.DataFrame(exc_rows)
    exc.to_csv(os.path.join(args.outdir, 'saturation_exceedance.tsv'),
               sep='\t', index=False)

    # --- B: omega against dS ----------------------------------------------
    pooled = pd.concat([per_pair[k][['dS', 'dN_dS']] for k in order],
                       ignore_index=True)
    binned, rho, rho_p = omega_vs_ds(pooled)
    logging.info('\n--- dN/dS against dS, all comparisons pooled ---')
    logging.info('  Spearman rho = %.3f, P = %s, n = %s', rho, fmt_p(rho_p),
                 f'{len(pooled.dropna()):,}')
    logging.info('  omega at the lowest dS bin %.4f, at the highest %.4f',
                 binned['omega_med'].iloc[0], binned['omega_med'].iloc[-1])
    binned.to_csv(os.path.join(args.outdir, 'omega_vs_ds_bins.tsv'),
                  sep='\t', index=False)

    within = omega_ds_within(per_pair)
    logging.info('  within-comparison Spearman (pooled rho could otherwise be '
                 'two clusters at different heights):')
    for _, r in within.iterrows():
        logging.info('    %-9s rho = %+.3f  P = %-14s n = %s',
                     r['comparison'], r['rho'], fmt_p(r['p']), f"{int(r['n']):,}")
    within.to_csv(os.path.join(args.outdir, 'omega_vs_ds_within.tsv'),
                  sep='\t', index=False)

    # --- C: fold change against cut-off -----------------------------------
    folds, skipped_cuts = None, []
    if args.annotated:
        if not os.path.isfile(args.annotated):
            sys.exit(f'ERROR: annotated table not found: {args.annotated}')
        ann = read_rate_table(args.annotated, require=('dS', 'dN_dS'))
        if args.region_col not in ann.columns:
            sys.exit(f"ERROR: no '{args.region_col}' column in {args.annotated}. "
                     f'Columns: {ann.columns.tolist()}')
        cutoffs = [None] + list(args.cutoffs)
        folds, skipped_cuts = fold_by_cutoff(ann, args.region_col, cutoffs, rng,
                                             n_boot=args.n_boot)
        logging.info('\n--- eliminated/retained dN/dS fold change by dS cut-off ---')
        for r in folds.itertuples():
            logging.info('  cutoff %-5s fold %.3f [%.3f, %.3f]  '
                         'n = %s / %s  P = %s',
                         'all' if pd.isna(r.cutoff) else f'<={r.cutoff:g}',
                         r.fold, r.lo, r.hi, f'{r.n_ret:,}', f'{r.n_elim:,}',
                         fmt_p(r.p))
        folds.to_csv(os.path.join(args.outdir, 'fold_by_ds_cutoff.tsv'),
                     sep='\t', index=False)
    else:
        logging.warning('No --annotated table given; panel C (the one that '
                        'actually tests whether the region contrast depends on '
                        'saturated pairs) is skipped.')

    paths, problems, fig_facts = plot_saturation(
        per_pair, order, ramp, binned, within, folds, args.outdir,
        'saturation_diagnostic', width=args.width, height_mm=args.height_mm,
        letters=not args.no_panel_letters, ds_view=args.ds_view,
        annotated_label=args.annotated_label if folds is not None else None)

    # --- summary -----------------------------------------------------------
    lines = ['dS saturation diagnostic (supplementary fig. S11; S12 before '
             'renumbering)', '=' * 60, '',
             'A. Fraction of pairs above each dS threshold', '',
             f"  {'comparison':<11}{'n':>9}{'median dS':>11}" +
             ''.join(f'{f">{t:g}":>9}' for t in THRESHOLDS), '']
    for _, r in exc.iterrows():
        pcts = ''.join(f'{100 * r[f"frac_gt_{t:g}"]:>8.1f}%' for t in THRESHOLDS)
        lines.append(f'  {r["comparison"]:<11}{r["n"]:>9,}'
                     f'{r["median_dS"]:>11.4f}{pcts}')
    lines += ['', 'B. dN/dS against dS, all six comparisons pooled', '',
              f'  Spearman rho = {rho:.3f}, P = {fmt_p(rho_p)}, '
              f'n = {len(pooled.dropna()):,}',
              f"  median dN/dS in the lowest dS bin  {binned['omega_med'].iloc[0]:.4f} "
              f"(dS {binned['ds_mid'].iloc[0]:.3f})",
              f"  median dN/dS in the highest dS bin {binned['omega_med'].iloc[-1]:.4f} "
              f"(dS {binned['ds_mid'].iloc[-1]:.3f})", '',
              '  A negative rho means dN/dS falls as dS rises. That is the',
              '  saturation signature: inflated dS estimates near saturation',
              '  push omega down. It is also why absolute dN/dS must not be',
              '  compared between comparisons of different divergence depth,',
              '  and why the six-way dN/dS boxplot is a pipeline check rather',
              '  than a manuscript figure.', '',
              '  Within each comparison separately, which removes comparison',
              '  identity as an explanation for the pooled trend:', '',
              f"  {'comparison':<11}{'n':>9}{'median dS':>11}{'rho':>8}{'P':>15}", '']
    for _, r in within.iterrows():
        lines.append(f"  {r['comparison']:<11}{int(r['n']):>9,}"
                     f"{r['median_dS']:>11.4f}{r['rho']:>+8.3f}"
                     f"{fmt_p(r['p']):>15}")
    lines += ['',
              '  The comparisons cluster in dS, so a pooled correlation alone',
              '  could be two groups at different heights rather than any',
              '  effect of dS. The within-comparison values are what close',
              '  that off -- the Toxocara comparisons span the saturation',
              '  threshold internally.']
    lo_grp = within[within['median_dS'] < 1.0]
    hi_grp = within[within['median_dS'] >= 1.0]
    if len(lo_grp) and len(hi_grp):
        lo_sign = 'positive' if (lo_grp['rho'] > 0).all() else (
            'negative' if (lo_grp['rho'] < 0).all() else 'mixed')
        hi_sign = 'positive' if (hi_grp['rho'] > 0).all() else (
            'negative' if (hi_grp['rho'] < 0).all() else 'mixed')
        lines += ['',
                  f'  Observed: the correlation is {lo_sign} in all '
                  f'{len(lo_grp)} comparisons below',
                  f'  median dS 1.0, and {hi_sign} in all {len(hi_grp)} above it.']
        if lo_sign == 'positive' and hi_sign == 'negative':
            lines += ['  A sign that flips with divergence depth is what',
                      '  saturation predicts and what a biological gradient in',
                      '  constraint does not: past saturation the multiple-hit',
                      '  correction overshoots, dS is overestimated, and omega',
                      '  is mechanically depressed. Toxocara omega values are',
                      '  therefore underestimates, which is a further reason not',
                      '  to compare them with the ascarid ones directly.']
    if folds is not None and len(folds):
        lines += ['', 'C. Eliminated/retained dN/dS fold change by dS cut-off', '',
                  f"  {'cutoff':<9}{'fold':>8}{'95% CI':>20}{'n ret':>9}"
                  f"{'n elim':>9}{'P':>15}", '']
        for r in folds.itertuples():
            ci = f'[{r.lo:.3f}, {r.hi:.3f}]'
            tag = 'all' if pd.isna(r.cutoff) else f'<={r.cutoff:g}'
            lines.append(f"  {tag:<9}"
                         f'{r.fold:>8.3f}{ci:>20}{r.n_ret:>9,}{r.n_elim:>9,}'
                         f'{fmt_p(r.p):>15}')
        first, last = folds.iloc[0], folds.iloc[-1]
        lines += ['',
                  f'  The fold change moves from {first.fold:.3f} with every pair '
                  f'to {last.fold:.3f}',
                  f'  at the strictest cut-off, while n falls from '
                  f'{int(first.n_elim):,} to {int(last.n_elim):,}',
                  '  eliminated genes. Read the intervals as well as the points:',
                  '  a flat line whose intervals widen to uselessness is not',
                  '  evidence of stability, it is loss of power.']
        if skipped_cuts:
            lines += ['',
                      '  Cut-offs dropped for too few genes (not plotted, so the',
                      '  panel does not imply they were never tried):']
            for cut, nr, ne in skipped_cuts:
                lines.append(f'    <={cut:g}: {nr:,} retained / {ne:,} eliminated')

    # --- legend notes: what every mark on the figure is --------------------
    rr = rho_range(within)
    beyond = '; '.join(f'{k} {v:,}' for k, v in fig_facts['A_beyond_view'].items()
                       if v)
    lines += ['', 'Legend notes (what saturation_diagnostic draws)', '',
              '  Comparisons (legend, ascarid pairs then outgroup pairs): '
              + '; '.join(DISPLAY[k] for k in order)
              + '. Light to dark = increasing median dS, as Figure 1C.',
              '  A: cumulative fraction of pairs by dS; dashed lines at dS = 1.0 '
              f'and 1.5. x-axis to dS = {args.ds_view:g}; pairs beyond it: '
              + (beyond or 'none') + '.',
              '  B: line = median dN/dS per bin of pooled dS (quantile bins, '
              '>= 50 pairs each); grey band = interquartile range per bin; '
              "points = each comparison's median dS and median dN/dS; dashed "
              'line at dS = 1.5.',
              '     Above the axes: Spearman rho of dN/dS on dS within each '
              'comparison, as the range over the '
              + ' and over the '.join(
                  f'{name} pairs ({n} comparisons, {range_text(lo, hi)})'
                  for name, (lo, hi, n) in rr.items())
              + f'. The pooled rho ({rho:.3f}, P '
              + (fmt_p(rho_p) if fmt_p(rho_p).startswith('<')
                 else f'= {fmt_p(rho_p)}')
              + ') is not drawn: it cannot be read on its own (see B above).']
    if folds is not None and len(folds):
        lines.append(
            f"  C: {args.annotated_label or 'the --annotated comparison'}, "
            f'{args.region_col} regions; points = median eliminated / median '
            'retained dN/dS (All = no cut-off), bars = 95% percentile bootstrap '
            f'interval ({args.n_boot:,} replicates, regions resampled '
            'separately); n = eliminated genes at the cut-off '
            f"({fig_facts.get('C_n_labels', '')}); dashed line = no difference.")
    lines += ['  Panel letters: ' + ('drawn' if not args.no_panel_letters
                                     else 'omitted') + '.',
              '  AUDIT: ' + '; '.join(problems) if problems else '  audit: clean']
    summary = '\n'.join(lines) + '\n'
    with open(os.path.join(args.outdir, 'saturation_summary.txt'), 'w') as fh:
        fh.write(summary)

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(),
            'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv,
            'args': vars(args),
            'inputs': {n: fingerprint(p) for n, p in spec},
            'annotated': fingerprint(args.annotated) if args.annotated else None,
            'rows': {n: int(len(per_pair[n])) for n, _ in spec},
            'figure': [os.path.basename(x) for x in paths],
            'drawn': fig_facts,
            'audit': problems,
        }, fh, indent=2, default=str)

    print()
    print(summary)
    logging.info('Summary: %s',
                 os.path.join(args.outdir, 'saturation_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
