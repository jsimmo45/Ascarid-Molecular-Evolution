#!/usr/bin/env python3
"""
cbr_meta_profile.py - dN/dS, dN and dS as a function of signed distance to a CBR.

Figure 3A, plus the single-chromosome track of Figure 3C (Figure 4A and 4C
before the renumbering of 2026-10-27). Replaces `kaks_landscape.py`.

Writes:

    cbr_meta_profile.{svg,pdf,png}        Figure 3A, three stacked panels
    cbr_meta_profile_dnds.*               dN/dS alone, if 3A is cut to one panel
    chromosome_track.*                    single-chromosome track (Figure 3C), with
                                          a key, and off-scale genes drawn hollow
                                          on the top edge rather than clipped
    cbr_meta_profile_binned.tsv           every plotted number, per bin per metric
    cbr_meta_profile_summary.txt          step size at the boundary, tests
    provenance.json                       inputs, arguments and versions


Changes 2026-10-30 (laid-out Fig. 3, panels A and C)
----------------------------------------------------
* The fold labels of 3A are plain Arial text, stars first, as in Fig. 4:
  "**** 1.69×", with a Unicode times sign (fold_label()). They were mathtext
  ($\\times$), which matplotlib writes glyph by glyph in its math font, so the
  whole label (digits and six-pointed asterisks included) opened in Illustrator
  in DejaVu Sans. The hand edit listed in script_to_figure_map.md section 2
  ("before any rerun") is now in the code. The single-panel variant uses the
  same label.
* The region labels ("Retained", "Eliminated") are upright, as in 3D: italic is
  for species names.
* The track (3C) is drawn for chromosome 2 by the wrapper (TRACK_CHROMOSOME,
  default chr02,chr01; was chr01): both of its terminal eliminated tracts are
  in the assembly, it carries the internal break that splits it into 2a and
  2b, and it is the chromosome drawn in Fig. 1A (chromosome 1 has an
  eliminated tract at its right end only). --chromosome now takes a
  comma-separated list: the first is the panel (chromosome_track.*, the name
  in the manifest), the others are drawn the same way as
  chromosome_track_<name>.* for comparison, and the summary gives each
  one's genes by region, eliminated tracts and CBRs. chr02's tracts are
  shorter than chr01's internal one, so compare the two before choosing.
* Statistics and tables are unchanged.

It no longer computes its own distances
---------------------------------------
The old version reparsed the GFF, re-derived CBR positions and recomputed a
signed distance, all of which `annotate_regions.py` had already done and written
to `kaks_annotated.tsv`. Two derivations of the same quantity is how Figure 3A
and Figure 3D came to disagree about which side of a CBR was positive: this
script signed **retained** positive while `cbr_changepoint_relax.py` fitted its
slopes with **eliminated** positive, so the two panels of one figure ran in
opposite directions.

It now reads `signed_cbr_dist1` from the table and draws it. Eliminated is
positive and on the right, throughout.

Given `--eliminated-bed` and `--chrom-sizes` it also **recomputes** the distance
through `region_lib` and asserts agreement with the table before drawing
anything. That is not redundancy for its own sake — it is the check that the
table was built with the CBR definition this figure assumes, and it fails loudly
if the table is stale. The recomputation is also what supplies the per-bin count
of contributing CBRs.


Ribbons are bootstrap intervals on the median
---------------------------------------------
The old version drew `groupby.sem()` around a median. The standard error of the
mean is not an interval for the median, and it is badly behaved on the
right-skewed per-gene rate distributions here: a bin whose genes include one
large omega gets a wide ribbon for a reason that has nothing to do with the
uncertainty in its median. That is what produced the unstable spike with a huge
ribbon on the dS panel right at the boundary, which is the first place a referee
looks.

Ribbons are now percentile bootstrap intervals on the median, matching every
other comparison in the paper.


Colour
------
All three traces are drawn in `NEUTRAL["ink"]`. Drawing dN red and dS blue - as
the old version did - collides with the region convention that blue means
retained and red means eliminated, in the one panel where region is also encoded
by the background shading. The metric is named on each panel's y-axis; it does
not need a hue.


Usage
-----
    python cbr_meta_profile.py \\
      --annotated results/annotate_regions_20260817/kaks_annotated.tsv \\
      --eliminated-bed AG.v50.eliminated_strict.bed \\
      --chrom-sizes AG.v50.fa.fai \\
      --species Ascaris \\
      --outdir results/cbr_meta_profile_20260817

Dependencies: numpy, pandas, matplotlib, scipy, paper_colors, region_lib
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from paper_colors import (REGION, REGION_SHADE, NEUTRAL, FONT, LINEWIDTH,
                          apply_style, audit, save, figure, panel_letter, mm)
from region_lib import (parse_bed, merge_intervals, load_chrom_sizes,
                        derive_cbrs, signed_distance, stars,
                        benjamini_hochberg, fingerprint, git_commit,
                        setup_logging)

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Panels of Figure 3A, top to bottom. dN/dS leads because it is the claim; dN
#: and dS below it are the decomposition that shows the claim is driven by dN.
METRICS = [('Ka_Ks', 'dN/dS'), ('Ka', 'dN'), ('Ks', 'dS')]


# ---------------------------------------------------------------------------
# Profile
# ---------------------------------------------------------------------------

def bootstrap_median_ci(values, rng, n_boot=2000, alpha=0.05):
    """Percentile bootstrap interval for a median.

    Resampling genes within the bin, which is the sampling unit. Returns
    (lo, hi), or (nan, nan) for a bin too small to resample meaningfully.
    """
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    if v.size < 3:
        return np.nan, np.nan
    idx = rng.integers(0, v.size, size=(n_boot, v.size))
    meds = np.median(v[idx], axis=1)
    return (float(np.quantile(meds, alpha / 2)),
            float(np.quantile(meds, 1 - alpha / 2)))


def build_profile(df, metric, dist_col, rng, window_kb=50, max_dist_kb=1000,
                  min_genes=10, n_boot=2000, cbr_col=None):
    """Median of `metric` in signed-distance bins, with bootstrap intervals.

    A bin is kept only if it holds at least `min_genes` genes. Without that
    floor the tail bins are medians of two or three genes drawn with a wide
    ribbon, and they dominate the reader's impression of a panel whose message
    is in the bins near zero.
    """
    w = window_kb * 1000.0
    max_dist = max_dist_kb * 1000.0
    edges = np.arange(-max_dist, max_dist + w, w)
    centres = (edges[:-1] + edges[1:]) / 2.0

    sub = df.dropna(subset=[dist_col, metric])
    sub = sub[(sub[dist_col] >= -max_dist) & (sub[dist_col] <= max_dist)]
    b = np.digitize(sub[dist_col].to_numpy(float), edges) - 1

    rows = []
    for i, c in enumerate(centres):
        v = sub[metric].to_numpy(float)[b == i]
        lo, hi = bootstrap_median_ci(v, rng, n_boot) if v.size else (np.nan, np.nan)
        row = {'metric': metric, 'centre_kb': c / 1000.0,
               'side': 'eliminated' if c > 0 else 'retained',
               'n_genes': int(v.size),
               'median': float(np.median(v)) if v.size else np.nan,
               'lo95': lo, 'hi95': hi,
               'q1': float(np.quantile(v, 0.25)) if v.size else np.nan,
               'q3': float(np.quantile(v, 0.75)) if v.size else np.nan}
        if cbr_col is not None:
            ids = sub[cbr_col].to_numpy()[b == i]
            row['n_cbr'] = int(pd.unique(ids[ids != '']).size) if ids.size else 0
        rows.append(row)
    prof = pd.DataFrame(rows)
    prof['used'] = prof['n_genes'] >= min_genes
    return prof


def boundary_step(df, metric, dist_col, within_kb=150.0):
    """Median each side of the CBR within `within_kb`, the fold, and a test.

    This is the number the panel is drawn to show, so it is computed rather than
    read off the trace. Restricted to the near window: the far tail on the
    eliminated side comes from a handful of long tracts and is not what "the step
    at the boundary" means.
    """
    d = df.dropna(subset=[dist_col, metric])
    w = within_kb * 1000.0
    ret = d.loc[(d[dist_col] < 0) & (d[dist_col] >= -w), metric].to_numpy(float)
    el = d.loc[(d[dist_col] > 0) & (d[dist_col] <= w), metric].to_numpy(float)
    if ret.size < 5 or el.size < 5:
        return None
    u, p = mannwhitneyu(ret, el, alternative='two-sided')
    mr, me = float(np.median(ret)), float(np.median(el))
    return {'metric': metric, 'within_kb': within_kb,
            'n_retained': int(ret.size), 'n_eliminated': int(el.size),
            'median_retained': mr, 'median_eliminated': me,
            'fold': me / mr if mr else np.nan, 'U': float(u), 'p': float(p)}


# ---------------------------------------------------------------------------
# Cross-check against the table
# ---------------------------------------------------------------------------

def recompute_and_check(df, eliminated_bed, chrom_sizes_path, terminal_tol,
                        dist_col, tol_bp=10.0):
    """Recompute the signed distance and assert it matches the table.

    The point is not the number - the table's number is used either way - it is
    that the table was built with the CBR definition this figure assumes. A
    stale table is silent otherwise: the profile would draw, and its boundary
    would sit somewhere subtly wrong.

    `tol_bp` is 10 bp: enough to absorb the rounding in the table's text
    representation, far below any real disagreement. A different CBR set moves a
    distance by thousands of base pairs, not by ones.

    Returns (per-gene CBR identifier, diagnostic). The identifier is what gives
    the per-bin count of contributing CBRs.
    """
    elim = parse_bed(eliminated_bed)
    merged = merge_intervals(elim)
    sizes, inferred = load_chrom_sizes(chrom_sizes_path, elim)
    cbrs, cbr_diag = derive_cbrs(merged, sizes, terminal_tol)

    recomputed = np.full(len(df), np.nan)
    cbr_id = np.full(len(df), '', dtype=object)
    have = df['chrom1'].notna() & df['midpoint1'].notna()
    for chrom, block in df[have].groupby('chrom1', sort=False):
        pos = block['midpoint1'].to_numpy(float)
        sd, idx = signed_distance(chrom, pos, cbrs, merged)
        loc = df.index.get_indexer(block.index)
        recomputed[loc] = sd
        cbr_id[loc] = [f'{chrom}:{i}' if not np.isnan(s) else ''
                       for s, i in zip(sd, idx)]

    both = ~np.isnan(recomputed) & df[dist_col].notna().to_numpy()
    delta = np.abs(recomputed[both] - df[dist_col].to_numpy(float)[both])
    n_bad = int((delta > tol_bp).sum())
    diag = {'n_compared': int(both.sum()), 'n_disagreeing': n_bad,
            'max_abs_delta_bp': float(delta.max()) if delta.size else 0.0,
            'cbrs': cbr_diag, 'chrom_sizes_inferred': inferred}
    if n_bad:
        sys.exit(
            f'ERROR: the annotated table disagrees with a fresh computation of '
            f'{dist_col} for {n_bad:,} of {int(both.sum()):,} genes '
            f'(max |delta| {delta.max():,.0f} bp).\n'
            f'  The table was built with a different CBR definition -- most '
            f'likely by the old annotate_kaks_regions.py, or with a different '
            f'--terminal-tol or a different eliminated BED.\n'
            f'  Rerun annotate_regions.py and point ANNOTATED_TABLE at it.')
    logging.info('  cross-check: %s genes agree with a fresh computation '
                 '(max |delta| %.1f bp)', f"{int(both.sum()):,}",
                 diag['max_abs_delta_bp'])
    return cbr_id, diag


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _trim(profiles):
    """x-limits from the bins that survived the count floor, both sides.

    The axis previously ran to +/- 1000 kb with shading to match while the data
    stopped near 700 kb on one side, so a third of the panel was empty coloured
    background asserting a range the analysis does not cover.
    """
    lo = min(p.loc[p['used'], 'centre_kb'].min() for p in profiles)
    hi = max(p.loc[p['used'], 'centre_kb'].max() for p in profiles)
    return lo, hi


def _draw(ax, prof, xlo, xhi, ylabel, show_ribbon=True):
    d = prof[prof['used']]
    ax.axvspan(xlo, 0, color=REGION_SHADE['retained'], zorder=0, lw=0)
    ax.axvspan(0, xhi, color=REGION_SHADE['eliminated'], zorder=0, lw=0)
    ax.axvline(0, color=NEUTRAL['cbr'], lw=1.0, zorder=4)
    if show_ribbon:
        ax.fill_between(d['centre_kb'], d['lo95'], d['hi95'],
                        color=NEUTRAL['ink'], alpha=0.18, lw=0, zorder=2)
    ax.plot(d['centre_kb'], d['median'], color=NEUTRAL['ink'], lw=1.2, zorder=3)
    ax.set_xlim(xlo, xhi)
    ax.set_ylabel(ylabel)
    return d


def _region_labels(ax, xlo, xhi, y=0.93):
    """Name the two sides inside the axes, in the region colours. Upright, as
    in Fig. 3D: italic is for species names (2026-10-30)."""
    ax.text(xlo * 0.55, y, 'Retained', transform=ax.get_xaxis_transform(),
            ha='center', va='top', color=REGION['retained'])
    ax.text(xhi * 0.55, y, 'Eliminated', transform=ax.get_xaxis_transform(),
            ha='center', va='top', color=REGION['eliminated'])


def fold_label(fold, q):
    """The step label of a profile panel, "**** 1.69×": the stars of `q`
    first, then the fold with a Unicode times sign, as Fig. 4 prints its folds.
    Plain text, not mathtext, so it opens in Illustrator in the label typeface
    (2026-10-30; a mathtext $\\times$ put the whole label in DejaVu Sans)."""
    return f'{stars(q)} {fold:.2f}×'


def plot_stack(profiles, steps, outdir, stem, species='', width='single',
               height_mm=130.0, panel_letters=False, show_ribbon=True):
    """Figure 3A: dN/dS, dN and dS stacked on a shared distance axis.

    Explicit mm margins rather than constrained_layout, which treats the axes as
    the flexible element and will shrink a panel toward zero rather than fail.
    """
    n = len(profiles)
    fig, ax0 = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax0)
    w_in, h_in = fig.get_size_inches()
    left, right, bottom, top, gap = mm(15.0), mm(3.0), mm(11.0), mm(3.0), mm(3.0)
    ah = (h_in - bottom - top - gap * (n - 1)) / n

    xlo, xhi = _trim(profiles)
    axes = []
    for i, ((metric, label), prof) in enumerate(zip(METRICS, profiles)):
        b = bottom + (n - 1 - i) * (ah + gap)
        ax = fig.add_axes([left / w_in, b / h_in,
                           (w_in - left - right) / w_in, ah / h_in])
        d = _draw(ax, prof, xlo, xhi, label, show_ribbon)
        lo = min(0.0, float(np.nanmin(d['lo95'])))
        # The top panel carries the region labels as well as its fold annotation,
        # so it needs more headroom and the two are placed on separate lines.
        ax.set_ylim(lo, float(np.nanmax(d['hi95'])) * (1.24 if i == 0 else 1.12))
        if i == 0:
            _region_labels(ax, xlo, xhi, y=0.97)
        if i < n - 1:
            ax.set_xticklabels([])
        if panel_letters:
            panel_letter(ax, 'ABC'[i])
        # The step is the panel's claim; state it rather than leave it to be
        # read off the trace.
        st = next((s for s in steps if s and s['metric'] == metric), None)
        if st:
            ax.annotate(fold_label(st['fold'], st.get('q', st['p'])),
                        xy=(0.985, 0.80 if i == 0 else 0.93),
                        xycoords='axes fraction',
                        ha='right', va='top', color=NEUTRAL['ink'],
                        fontsize=plt.rcParams['xtick.labelsize'])
        axes.append(ax)

    axes[-1].set_xlabel('Distance from nearest CBR (kb)')
    if species:
        axes[0].set_title(species, fontstyle='italic', pad=2)
    for problem in audit(fig, verbose=False):
        logging.warning('  audit: %s', problem)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths


def plot_single(prof, step, outdir, stem, ylabel='dN/dS', species='',
                width='single', height_mm=55.0, show_ribbon=True):
    """dN/dS alone, in case Figure 3A is cut to one panel at layout."""
    fig, ax0 = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax0)
    w_in, h_in = fig.get_size_inches()
    left, right, bottom, top = mm(15.0), mm(3.0), mm(11.0), mm(4.0)
    ax = fig.add_axes([left / w_in, bottom / h_in,
                       (w_in - left - right) / w_in,
                       (h_in - bottom - top) / h_in])
    xlo, xhi = _trim([prof])
    d = _draw(ax, prof, xlo, xhi, ylabel, show_ribbon)
    # Headroom for two lines of annotation: the region labels sit at the top of
    # the axes and the fold annotation goes below them. At 1.15 they collided,
    # which is visible in any assembly that uses this single-panel variant.
    ax.set_ylim(min(0.0, float(np.nanmin(d['lo95']))),
                float(np.nanmax(d['hi95'])) * 1.28)
    _region_labels(ax, xlo, xhi, y=0.97)
    ax.set_xlabel('Distance from nearest CBR (kb)')
    if step:
        ax.annotate(fold_label(step['fold'], step.get('q', step['p'])),
                    xy=(0.985, 0.82), xycoords='axes fraction', ha='right',
                    va='top', color=NEUTRAL['ink'],
                    fontsize=plt.rcParams['xtick.labelsize'])
    if species:
        ax.set_title(species, fontstyle='italic', pad=2)
    for problem in audit(fig, verbose=False):
        logging.warning('  audit: %s', problem)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths


def _legend_above(fig, ax, handles, labels, handler_map=None, gap_mm=2.0,
                  top_pad_mm=0.8):
    """Put the key above the axes, in two rows, or three if two won't fit.

    The axes box is shortened from the top to make room, and the key is never
    pushed into the axes, so it cannot sit on data; audit() checks this.

    Two rows is the fixed layout rather than a fallback. One row of six entries
    is about 178 mm in DejaVu Sans, which is the font the compute nodes render
    with, against 167 mm of panel width. A layout that fits on a laptop with
    Arial but wraps on Isaac would give the same figure a different height
    depending on where it was built. Three rows are used only if two are wider
    than the panel, as can happen with four-digit gene counts.

    Returns (legend, number of columns used).
    """
    w_in, h_in = fig.get_size_inches()
    pos = ax.get_position()
    rend = fig.canvas.get_renderer()
    n = len(handles)
    panel_w = pos.width * w_in
    tried = sorted({int(np.ceil(n / 2)), int(np.ceil(n / 3))}, reverse=True)
    for ncol in tried:
        leg = ax.legend(handles, labels, ncol=ncol, loc='lower left',
                        bbox_to_anchor=(0.0, 1.0), handler_map=handler_map,
                        frameon=False, fontsize=FONT['legend'],
                        borderaxespad=0.0, borderpad=0.0)
        fig.canvas.draw()
        bb = leg.get_window_extent(rend)
        if bb.width / fig.dpi <= panel_w:
            break
    else:
        logging.warning('  chromosome-track key is %.0f mm wide even in %d '
                        'rows, against %.0f mm of panel; shorten its labels',
                        bb.width / fig.dpi * 25.4, int(np.ceil(n / ncol)),
                        panel_w * 25.4)
    top_in = h_in - mm(top_pad_mm) - bb.height / fig.dpi - mm(gap_mm)
    bottom_in = pos.y0 * h_in
    ax.set_position([pos.x0, pos.y0, pos.width, (top_in - bottom_in) / h_in])
    leg.set_bbox_to_anchor((0.0, 1.0 + mm(gap_mm) / (top_in - bottom_in)),
                           transform=ax.transAxes)
    logging.info('  key: %d rows x %d columns, %.0f of %.0f mm; plotting area '
                 '%.1f mm high', int(np.ceil(n / ncol)), ncol,
                 bb.width / fig.dpi * 25.4, panel_w * 25.4,
                 (top_in - bottom_in) * 25.4)
    return leg, ncol


def plot_chromosome(df, chrom, elim_merged, cbrs, outdir, stem, metric='Ka_Ks',
                    ylabel='dN/dS', width='double', height_mm=57.0,
                    window_genes=25, cap_pct=98.0, rng=None):
    """One chromosome, per-gene points under a sliding-window median.

    Figure 3C: it grounds the meta-profile in a real chromosome, where the
    reader can see the eliminated tracts and the CBR positions that the profile
    collapses onto one axis.

    The y-axis stops at the `cap_pct` percentile so the bulk of the genes stays
    readable. Genes above the cap are not dropped from the page. Each one is
    drawn as a hollow circle in its region colour, on the top edge at its true
    position, and the key and the summary give their number. Before 2026-10-08
    they were clipped out of the axes with nothing to show they existed: a
    silent cap.

    The key sits above the axes, so the chromosome name has moved from a title
    into the x-axis label. Height went from 52 to 57 mm so the plotting
    area stays at the 37 mm it was.

    Returns a dict of what was drawn (counts per region, the cap, the
    off-scale genes) for the summary and provenance, or None if the chromosome
    has too few genes to draw.
    """
    d = df[(df['chrom1'] == chrom)].dropna(subset=[metric, 'midpoint1'])
    d = d.sort_values('midpoint1')
    if len(d) < 20:
        logging.warning('  %s: only %d genes, skipping track', chrom, len(d))
        return None

    fig, ax0 = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax0)
    w_in, h_in = fig.get_size_inches()
    # The top margin is provisional: _legend_above() shortens the axes from the
    # top to fit the key once the key's size is known.
    left, right, bottom, top = mm(15.0), mm(3.0), mm(11.0), mm(4.0)
    ax = fig.add_axes([left / w_in, bottom / h_in,
                       (w_in - left - right) / w_in,
                       (h_in - bottom - top) / h_in])

    x = d['midpoint1'].to_numpy(float) / 1e6
    y = d[metric].to_numpy(float)
    # Anything that is not retained or eliminated (ambiguous, unannotated) is
    # drawn grey under one key entry, so no colour on the panel goes unkeyed.
    region = d['region1'].to_numpy(dtype=object)
    region = np.where(np.isin(region, ['retained', 'eliminated']), region,
                      'unassigned')
    cap = float(np.quantile(y, cap_pct / 100.0))
    over = y > cap

    starts, ends = elim_merged.get(chrom, (np.array([]), np.array([])))
    for s, e in zip(starts, ends):
        ax.axvspan(s / 1e6, e / 1e6, color=REGION_SHADE['eliminated'],
                   zorder=0, lw=0)
    chrom_cbrs = cbrs.get(chrom, np.array([]))
    # CBR lines are drawn as segments in data coordinates spanning the fixed
    # y-range, not with axvline. audit() maps a Line2D through the data
    # transform, which puts axvline's axes-fraction ends outside the panel, and
    # it would then report a false collision with the key above the axes.
    for c in chrom_cbrs:
        ax.plot([c / 1e6, c / 1e6], [0.0, cap], color=NEUTRAL['cbr'], lw=0.6,
                alpha=0.85, zorder=1, solid_capstyle='butt')

    present = set(region.tolist())
    order = [r for r in ('retained', 'eliminated', 'unassigned') if r in present]
    colour = {r: REGION.get(r, NEUTRAL['ink_light']) for r in order}
    n_by_region, n_over_by_region = {}, {}
    for r in order:
        m = region == r
        n_by_region[r] = int(m.sum())
        n_over_by_region[r] = int((m & over).sum())
        keep = m & ~over
        # Retained first and eliminated on top, so the minority class is never
        # hidden under the majority.
        ax.scatter(x[keep], y[keep], s=2.0, color=colour[r], alpha=0.45,
                   linewidths=0, zorder=2.5 if r == 'eliminated' else 2,
                   rasterized=True)
        hi = m & over
        if hi.any():
            # Off-scale genes: hollow, in the region colour, pinned to the top
            # edge at their true x. clip_on=False so the axes boundary they sit
            # on does not cut the marker in half.
            ax.scatter(x[hi], np.full(int(hi.sum()), cap), s=7.0,
                       facecolors='none', edgecolors=colour[r],
                       linewidths=LINEWIDTH['thin'], zorder=4, clip_on=False)

    roll = pd.Series(y).rolling(window_genes, center=True,
                                min_periods=max(5, window_genes // 3)).median()
    ax.plot(x, roll, color=NEUTRAL['ink'], lw=1.0, zorder=3)

    ax.set_ylim(0, cap)
    ax.set_xlim(0, x.max() * 1.005)
    ax.set_xlabel(f'Position on {chrom} (Mb)')
    ax.set_ylabel(ylabel)

    # --- key: every mark on the panel, in the panel's own colours ---------
    handles, labels = [], []
    for r in order:
        handles.append(Line2D([], [], ls='none', marker='o', markersize=3.0,
                              markerfacecolor=colour[r], markeredgecolor='none',
                              alpha=0.45))
        labels.append(f'{r.capitalize()} (n = {n_by_region[r]:,})')
    n_over = int(over.sum())
    if n_over:
        # One entry for both regions, as a paired handle in the regions that
        # actually have off-scale genes.
        hollow = tuple(Line2D([], [], ls='none', marker='o', markersize=3.6,
                              markerfacecolor='none', markeredgecolor=colour[r],
                              markeredgewidth=LINEWIDTH['thin'])
                       for r in order if n_over_by_region[r])
        handles.append(hollow if len(hollow) > 1 else hollow[0])
        labels.append(f'Off scale, > {cap:.2f} (n = {n_over:,})')
    handles.append(Line2D([], [], color=NEUTRAL['ink'], lw=1.0))
    labels.append(f'Running median, {window_genes} genes')
    if chrom_cbrs.size:
        handles.append(Line2D([], [], ls='none', marker='|', markersize=7.0,
                              markeredgewidth=0.8, color=NEUTRAL['cbr']))
        labels.append('CBR')
    if starts.size:
        handles.append(Patch(facecolor=REGION_SHADE['eliminated'],
                             edgecolor='none'))
        labels.append('Eliminated tract')
    leg, ncol = _legend_above(
        fig, ax, handles, labels,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})

    for problem in audit(fig, verbose=False):
        logging.warning('  audit: %s', problem)
    paths = save(fig, outdir, stem)
    plt.close(fig)

    split = ', '.join(f'{n_over_by_region[r]} {r}' for r in order
                      if n_over_by_region[r])
    logging.info('  wrote %s (%d genes; y-axis capped at the %.0fth percentile, '
                 '%s %.3f; %d off-scale gene(s) drawn hollow on the top edge%s)',
                 ', '.join(os.path.basename(p) for p in paths), len(d),
                 cap_pct, ylabel, cap, n_over, f': {split}' if split else '')
    return {'chromosome': chrom, 'metric': metric, 'n_genes': int(len(d)),
            'n_by_region': n_by_region, 'cap_pct': float(cap_pct),
            'cap_value': cap, 'n_off_scale': n_over,
            'n_off_scale_by_region': n_over_by_region,
            'window_genes': int(window_genes), 'n_cbrs': int(chrom_cbrs.size),
            'key_rows': int(np.ceil(len(handles) / ncol)),
            'files': [os.path.basename(p) for p in paths]}


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='dN/dS, dN and dS meta-profiles against signed CBR distance '
                    '(Figure 3A), plus a single-chromosome track.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', required=True,
                   help='kaks_annotated.tsv from annotate_regions.py. Must carry '
                        'signed_cbr_dist1 -- tables from the old '
                        'annotate_kaks_regions.py do not.')
    p.add_argument('--eliminated-bed', default=None,
                   help='Eliminated regions BED. Used to cross-check the table '
                        'and to draw the chromosome track. Strongly recommended: '
                        'without it a stale table is not detected.')
    p.add_argument('--chrom-sizes', default=None,
                   help='.fai or two-column TSV, for the cross-check')
    p.add_argument('--terminal-tol', type=int, default=1000,
                   help='Must match the value annotate_regions.py was run with')
    p.add_argument('--check-tol-bp', type=float, default=10.0,
                   help='Tolerance for the table cross-check (default: 10 bp, '
                        'enough to absorb text rounding and far below any real '
                        'disagreement)')
    p.add_argument('--dist-col', default='signed_cbr_dist1')
    p.add_argument('--window-kb', type=float, default=50.0)
    p.add_argument('--max-dist-kb', type=float, default=1000.0,
                   help='Bins are built out to here; the drawn axis is then '
                        'trimmed to the bins that pass --min-genes')
    p.add_argument('--min-genes', type=int, default=10,
                   help='Minimum genes for a bin to be drawn (default: 10)')
    p.add_argument('--step-window-kb', type=float, default=150.0,
                   help='Window each side of the CBR for the reported step')
    p.add_argument('--n-boot', type=int, default=2000)
    p.add_argument('--no-ribbon', dest='show_ribbon', action='store_false')
    p.add_argument('--chromosome', default=None,
                   help='Chromosome for the track panel; default is the one '
                        'with the most annotated genes. Several, comma-'
                        'separated, draw the first as chromosome_track.* and '
                        'each other as chromosome_track_<name>.*, for choosing '
                        'between them (chr02,chr01 in the wrapper)')
    p.add_argument('--species', default='')
    p.add_argument('--width', choices=['single', 'onehalf', 'double'],
                   default='single')
    p.add_argument('--height-mm', type=float, default=130.0)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--panel-letters', action='store_true')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'cbr_meta_profile_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'cbr_meta_profile')
    # matplotlib's PDF writer subsets fonts through fontTools, which logs every
    # table it prunes at INFO -- several hundred lines per figure, burying the
    # lines this script writes. Its warnings still come through.
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.info('cbr_meta_profile.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()
    rng = np.random.default_rng(args.seed)

    df = pd.read_csv(args.annotated, sep='\t')
    logging.info('Loaded %s rows from %s', f'{len(df):,}',
                 os.path.basename(args.annotated))
    if args.dist_col not in df.columns:
        sys.exit(f'ERROR: {args.dist_col!r} is not a column of {args.annotated}.\n'
                 f'  This table predates annotate_regions.py. Rerun it and point '
                 f'ANNOTATED_TABLE at the new run.\n'
                 f'  columns: {list(df.columns)}')
    missing = [m for m, _ in METRICS if m not in df.columns]
    if missing:
        sys.exit(f'ERROR: missing rate columns {missing} in {args.annotated}')

    # --- cross-check ------------------------------------------------------
    cbr_id, check_diag, elim_merged, cbrs = None, None, {}, {}
    if args.eliminated_bed:
        logging.info('Cross-checking %s against a fresh computation',
                     args.dist_col)
        elim = parse_bed(args.eliminated_bed)
        elim_merged = merge_intervals(elim)
        sizes, _ = load_chrom_sizes(args.chrom_sizes, elim)
        cbrs, _ = derive_cbrs(elim_merged, sizes, args.terminal_tol)
        cbr_id, check_diag = recompute_and_check(
            df, args.eliminated_bed, args.chrom_sizes, args.terminal_tol,
            args.dist_col, tol_bp=args.check_tol_bp)
        df['_cbr_id'] = cbr_id
    else:
        logging.warning('No --eliminated-bed: the table is taken on trust, and a '
                        'stale one will not be detected.')

    # --- profiles ---------------------------------------------------------
    profiles, steps = [], []
    for metric, label in METRICS:
        prof = build_profile(df, metric, args.dist_col, rng,
                             window_kb=args.window_kb,
                             max_dist_kb=args.max_dist_kb,
                             min_genes=args.min_genes, n_boot=args.n_boot,
                             cbr_col='_cbr_id' if cbr_id is not None else None)
        profiles.append(prof)
        st = boundary_step(df, metric, args.dist_col, args.step_window_kb)
        steps.append(st)
        drawn = prof[prof['used']]
        logging.info('  %-6s %d of %d bins pass the %d-gene floor, '
                     'x from %.0f to %.0f kb', label, len(drawn), len(prof),
                     args.min_genes, drawn['centre_kb'].min(),
                     drawn['centre_kb'].max())

    q = benjamini_hochberg([s['p'] if s else np.nan for s in steps])
    for s, qi in zip(steps, q):
        if s:
            s['q'] = float(qi)

    pd.concat(profiles).to_csv(
        os.path.join(args.outdir, 'cbr_meta_profile_binned.tsv'),
        sep='\t', index=False, float_format='%.6g')

    # --- figures ----------------------------------------------------------
    logging.info('Figure 3A')
    plot_stack(profiles, steps, args.outdir, 'cbr_meta_profile',
               species=args.species, width=args.width,
               height_mm=args.height_mm, panel_letters=args.panel_letters,
               show_ribbon=args.show_ribbon)
    plot_single(profiles[0], steps[0], args.outdir, 'cbr_meta_profile_dnds',
                ylabel=METRICS[0][1], species=args.species, width=args.width,
                show_ribbon=args.show_ribbon)

    track, comparison_tracks = None, []
    if elim_merged:
        chroms = ([c.strip() for c in args.chromosome.split(',') if c.strip()]
                  if args.chromosome else [])
        if not chroms:
            counts = df.dropna(subset=['Ka_Ks'])['chrom1'].value_counts()
            counts = counts[counts.index.isin(elim_merged)]
            chroms = [counts.index[0]] if len(counts) else []
        if len(set(chroms)) != len(chroms):
            raise SystemExit(f'--chromosome names a chromosome twice: {chroms}')
        for k, chrom in enumerate(chroms):
            # The first is the figure panel (chromosome_track.*, the name the
            # manifest uses); the others are drawn identically for comparison.
            stem = 'chromosome_track' if k == 0 else f'chromosome_track_{chrom}'
            logging.info('Chromosome track: %s (%s)', chrom, stem)
            drawn = plot_chromosome(df, chrom, elim_merged, cbrs, args.outdir,
                                    stem, rng=rng)
            if drawn:
                starts, ends = elim_merged.get(chrom, (np.array([]), np.array([])))
                drawn['n_tracts'] = int(len(starts))
                drawn['tract_kb'] = float(np.sum(np.asarray(ends) -
                                                 np.asarray(starts)) / 1e3)
            if k == 0:
                track = drawn
            elif drawn:
                comparison_tracks.append(drawn)
    else:
        logging.info('No --eliminated-bed: skipping the chromosome track')

    # --- summary ----------------------------------------------------------
    lines = ['dN/dS, dN and dS against signed distance to the nearest CBR', '',
             f'Table               {os.path.basename(args.annotated)}, '
             f'{len(df):,} rows',
             f'Distance column     {args.dist_col} (positive = eliminated)',
             f'Bin width           {args.window_kb:g} kb',
             f'Minimum per bin     {args.min_genes} genes',
             f'Bootstrap           {args.n_boot:,} resamples, percentile', '']
    if check_diag:
        lines += [f"Cross-check         {check_diag['n_compared']:,} genes "
                  f"recomputed, {check_diag['n_disagreeing']} disagree "
                  f"(max |delta| {check_diag['max_abs_delta_bp']:.1f} bp)",
                  f"CBRs                {check_diag['cbrs']['n_cbrs']:,} on "
                  f"{check_diag['cbrs']['n_sequences_with_cbrs']} sequences", '']

    lines += [f'Step at the boundary, within {args.step_window_kb:g} kb each side',
              '']
    for (metric, label), s in zip(METRICS, steps):
        if not s:
            lines.append(f'  {label:<6} not enough genes')
            continue
        lines.append(f"  {label:<6} retained {s['median_retained']:.4f} "
                     f"(n = {s['n_retained']:,})  ->  eliminated "
                     f"{s['median_eliminated']:.4f} (n = {s['n_eliminated']:,})  "
                     f"{s['fold']:.2f}x, P = {s['p']:.3g}, q = {s['q']:.3g} "
                     f"({stars(s['q'])})")
    def track_lines(t, title):
        nr = t['n_by_region']
        no = t['n_off_scale_by_region']
        split = ', '.join(f'{nr[r]:,} {r}' for r in nr)
        osplit = ', '.join(f'{no[r]} {r}' for r in no if no[r])
        return ['',
                f"{title:<20}{t['chromosome']}, "
                f"{t['n_genes']:,} genes ({split}); file {t['files'][0]}",
                f"                    eliminated tracts: {t['n_tracts']}, "
                f"{t['tract_kb']:,.0f} kb in all; {t['n_cbrs']} CBRs",
                f"                    y-axis capped at the "
                f"{t['cap_pct']:.0f}th percentile "
                f"({METRICS[0][1]} {t['cap_value']:.3f})",
                f"                    {t['n_off_scale']} "
                f"{'gene' if t['n_off_scale'] == 1 else 'genes'} above the "
                f"cap, drawn as hollow circles on",
                f"                    the top edge"
                + (f" ({osplit})" if osplit else ''),
                f"                    running median over "
                f"{t['window_genes']} genes"]

    if track:
        lines += track_lines(track, 'Chromosome track')
    for t in comparison_tracks:
        lines += track_lines(t, 'Comparison track')
    lines += ['',
              'Notes', '',
              '  Distances are read from the annotated table, not recomputed for',
              '  this figure. Two derivations of one quantity is how 3A and 3D',
              '  came to disagree about which side of a CBR is positive.', '',
              '  Ribbons are percentile bootstrap intervals on the median. The',
              '  earlier version drew the standard error of the mean around a',
              '  median, which is not an interval for it and which produced the',
              '  unstable wide ribbon on the dS panel at the boundary.', '',
              '  All three traces are drawn in the neutral ink. Region is carried',
              '  by the background shading, so a red dN trace and a blue dS trace',
              '  would collide with the convention in the one panel that uses',
              '  both.', '',
              f'  Bins holding fewer than {args.min_genes} genes are computed and',
              '  written to the TSV but not drawn, and the axis is trimmed to the',
              '  bins that remain rather than run to the nominal maximum over',
              '  empty shading.', '']
    summary = '\n'.join(lines)
    with open(os.path.join(args.outdir, 'cbr_meta_profile_summary.txt'), 'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(), 'python': platform.python_version(),
                   'packages': {'numpy': np.__version__,
                                'pandas': pd.__version__,
                                'matplotlib': matplotlib.__version__},
                   'argv': sys.argv, 'args': vars(args),
                   'inputs': {'annotated': fingerprint(args.annotated),
                              'eliminated_bed': fingerprint(args.eliminated_bed),
                              'chrom_sizes': fingerprint(args.chrom_sizes)},
                   'cross_check': check_diag,
                   'chromosome_track': track,
                   'comparison_tracks': comparison_tracks,
                   'steps': steps}, fh, indent=2, default=float)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
