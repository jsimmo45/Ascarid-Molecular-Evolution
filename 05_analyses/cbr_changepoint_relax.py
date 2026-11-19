#!/usr/bin/env python3
"""
cbr_changepoint_relax.py - how sharp the constraint boundary is, and how relaxed.

Figure 3D (per-arm slopes and the step at the CBR; Figure 4D before the
renumbering of 2026-10-27, so the dated notes of 2026-10-12 say 4D) and
supplementary fig. S3B (the relaxation exponent k; S8B before the 2026-10-12
renumbering). These were
one three-panel figure; they are now two, because they answer different
questions and belong in different places.

Writes:

    cbr_slope_step.{svg,pdf,png}     Figure 3D
    S8B_relaxation_k.*               supplementary fig. S3B. The file keeps its
                                     old stem so the assembly manifest
                                     (supplementary_figure_sources.tsv) and
                                     figure_locations stay valid; the assembly
                                     copies it to FigS03_B_relaxation_k.
    cbr_changepoint_stats.tsv        every fitted number, with bootstrap CIs
    cbr_changepoint_summary.txt      the same, formatted for reading, with the
                                     gene counts the two legends need
    provenance.json                  inputs, arguments and versions


Changes 2026-10-30 (Fig. 3D and S3B labels in one typeface)
-----------------------------------------------------------
* The axis label log10(dN/dS) no longer uses Unicode subscript digits (₁₀).
  Arial has none, so Illustrator set them in Myriad Pro in the laid-out Fig. 3D
  (the house rule of 2026-10-15). The label is set as plain text so the layout
  reserves its room, then redrawn as runs with the 10 as its own 6 pt Arial text
  lowered 0.15 of the label size (plain_text.subscript_runs): on the y axis of
  Fig. 3D with plain_text.replace_ylabel, on the x axis of S3B with
  replace_xlabel. Statistics, tables and every other mark are unchanged.
* 2026-10-31: fontTools' INFO logging is silenced, as in cbr_meta_profile.py;
  it added several hundred lines per PDF to the SLURM log.

Changes 2026-10-26 (Jianbin Wang's review)
------------------------------------------
* The step panel's x-axis reads "Distance from nearest CBR (kb)", as in the
  profile and gene-density panels beside it, instead of "Signed distance".
  The sign is still carried by the axis: retained to the left (negative),
  eliminated to the right, with each arm labelled above the axes.

Changes 2026-10-12 (figure review of supplementary fig. S3)
-----------------------------------------------------------
* Plain text, no mathtext: the axis label is log10(dN/dS) with upright dN/dS
  and Unicode subscript digits, the step reads "1.78x" with a Unicode times
  sign, and k is plain text. Matplotlib writes mathtext glyph by glyph in its
  math font, so these labels opened in Illustrator in a second typeface and
  with italic N and S, against the house rule of 2026-10-09.
* No species titles. --species is still accepted and goes into the summary;
  the figure specification puts the species in the caption.
* S3B: the legend and k sit above the axes, so no line crosses them (the legend
  used to sit in the upper right, crossed by both mean lines and the neutral
  rule). The "neutral" label sits at the top of the omega = 1 rule, above both
  densities, instead of at its foot over the filled curves; the x-range is
  extended past omega = 1 when the label would otherwise leave the axes.
* Figure 4D: the region labels sit above the axes, one per arm, instead of on
  the point cloud, and are upright (italic is for species names).
* The summary counts the genes each filter removes, per region, so the S3B
  legend can say where 9,002 / 320 come from (dN/dS = 0 has no logarithm;
  above --omega-max is dropped as saturated dS), and the genes behind the
  Figure 4D n (a CBR on the sequence, inside the slope window).
No statistic changes: every fit and bootstrap runs before any figure is drawn,
from the same random stream in the same order, and cbr_changepoint_stats.tsv
is byte-identical to the previous version's on the same input.


What Figure 3D tests
--------------------
Within each region *separately*, log10(omega) is regressed on signed distance
from the nearest CBR, restricted to |d| <= `--slope-window`:

    log10 omega(d) = alpha + beta * d

Two flat slopes plus a large difference between the two arms' boundary
intercepts (alpha_E - alpha_R, the **step**) is what "sharp boundary" means
operationally: constraint does not grade with distance inside either compartment,
and it changes discontinuously between them.

**The step now carries a bootstrap confidence interval.** It is the headline
number of the panel and the earlier version reported it as a point estimate while
giving CIs to the slopes, which are the supporting result.

**A free changepoint search is reported but is not evidence.** The sign of the
distance is defined by region membership, and region membership changes exactly
at the CBR, so a search for the best-fitting breakpoint will return zero whenever
the two groups differ at all. It is printed as a consistency check and labelled
as circular; do not quote it.


What S3B tests
--------------
The RELAX relationship omega_test = omega_ref^k, applied to the per-gene omega
distributions rather than to branches:

    k = mean(log omega_eliminated) / mean(log omega_retained)

k < 1 is relaxation - eliminated genes shifted toward neutrality. **This is not
HyPhy RELAX**, which labels branches; PDE is ancestral across all three lineages
here, so the contrast is between gene sets and not between branches. The
docstring said so before and the summary now says so too, because it is the first
thing a referee will ask.

Inference on k is the bootstrap interval, resampling each gene set independently.
A Welch test on log omega is reported beside it because it is the exact null for
k = 1 (equal means of log omega); Mann-Whitney is reported as the
distribution-free companion but tests stochastic dominance, which is a slightly
different statement.


It reads the distance, it does not derive it
--------------------------------------------
`signed_cbr_dist1` comes from `annotate_regions.py`. The old version took the
unsigned `min_cbr_dist1` and re-signed it by category, and had a second mode that
re-derived CBRs from a BED. Both are gone: one derivation, in one place, checked
here against a fresh `region_lib` computation before anything is fitted.

Usage
-----
    python cbr_changepoint_relax.py \\
      --annotated results/annotate_regions_20260817/kaks_annotated.tsv \\
      --eliminated-bed AG.v50.eliminated_strict.bed \\
      --chrom-sizes AG.v50.fa.fai \\
      --slope-window 150 --species Ascaris \\
      --outdir results/cbr_changepoint_20260817

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
from scipy.stats import (gaussian_kde, linregress, mannwhitneyu, norm,
                         ttest_ind)

from paper_colors import (REGION, REGION_SHADE, NEUTRAL, FONT, apply_style,
                          audit, save, figure, panel_letter, mm)
from plain_text import replace_xlabel, replace_ylabel, subscript_runs
from region_lib import (parse_bed, merge_intervals, load_chrom_sizes,
                        derive_cbrs, signed_distance, fmt_p, stars,
                        benjamini_hochberg, fingerprint, git_commit,
                        setup_logging)

_HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------

def fit_arms(x, y, window_kb, min_n=10):
    """Per-region regression of log10(omega) on signed distance, in kb.

    The intercept is at d = 0, i.e. the fitted boundary value for that arm, so
    the difference between the two intercepts is the step across the CBR.
    Restricted to |d| <= window_kb: beyond that the eliminated side is carried by
    a handful of long tracts and the fit starts describing those tracts rather
    than the boundary.
    """
    arms = {}
    for name, side in (('retained', x < 0), ('eliminated', x > 0)):
        m = side & (np.abs(x) <= window_kb)
        if m.sum() < min_n:
            arms[name] = None
            continue
        lr = linregress(x[m], y[m])
        arms[name] = {'slope_per_100kb': lr.slope * 100.0,
                      'intercept': lr.intercept, 'p': lr.pvalue,
                      'stderr_per_100kb': lr.stderr * 100.0,
                      'intercept_stderr': float(getattr(lr, 'intercept_stderr',
                                                        np.nan)),
                      'n': int(m.sum())}
    step = None
    if arms['retained'] and arms['eliminated']:
        step = arms['eliminated']['intercept'] - arms['retained']['intercept']
    return arms, step


def step_analytic_p(arms):
    """Wald test on the step, from the two intercept standard errors.

    The bootstrap gives the interval, but its p-value is floored at 1 / n_boot -
    it cannot resolve below its own resolution, and a headline number reported as
    'P < 0.0005' when the truth is 10^-20 reads as weak evidence. The two arms are
    fitted on disjoint genes, so their intercepts are independent and
    SE(step) = sqrt(SE_R^2 + SE_E^2).
    """
    a, b = arms.get('retained'), arms.get('eliminated')
    if not a or not b:
        return np.nan, np.nan
    se = np.sqrt(a['intercept_stderr'] ** 2 + b['intercept_stderr'] ** 2)
    if not np.isfinite(se) or se == 0:
        return np.nan, np.nan
    z = (b['intercept'] - a['intercept']) / se
    return float(se), float(2.0 * norm.sf(abs(z)))


def bootstrap_arms(x, y, window_kb, n_boot, rng, min_n=10):
    """Resample genes within the window and refit, for CIs on slopes and step.

    Resampling is stratified by arm: each side is resampled to its own size, so a
    draw cannot happen to empty one arm and the two sides keep the sample sizes
    the point estimate was computed at.
    """
    keep = np.abs(x) <= window_kb
    xs, ys = x[keep], y[keep]
    left, right = np.flatnonzero(xs < 0), np.flatnonzero(xs > 0)
    if left.size < min_n or right.size < min_n:
        return {}
    out = {'slope_retained': [], 'slope_eliminated': [], 'step': []}
    for _ in range(n_boot):
        li = left[rng.integers(0, left.size, left.size)]
        ri = right[rng.integers(0, right.size, right.size)]
        lr_l = linregress(xs[li], ys[li])
        lr_r = linregress(xs[ri], ys[ri])
        out['slope_retained'].append(lr_l.slope * 100.0)
        out['slope_eliminated'].append(lr_r.slope * 100.0)
        out['step'].append(lr_r.intercept - lr_l.intercept)
    return {k: np.asarray(v) for k, v in out.items()}


def ci95(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (np.nan, np.nan)
    return tuple(np.percentile(v, [2.5, 97.5]))


def boot_p(v, null=0.0):
    """Two-sided bootstrap p: twice the smaller tail mass beyond the null.

    Floored at 1/n_boot rather than reported as zero, because a bootstrap cannot
    resolve a p-value below its own resolution and printing 0 invites it to be
    quoted as one.
    """
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    frac = min((v <= null).mean(), (v >= null).mean())
    return max(2.0 * frac, 1.0 / v.size)


def free_changepoint(x, y, n_grid=200):
    """Best-fitting step location by grid search.

    **Circular by construction** and reported only as a consistency check: the
    sign of x is region membership, and region membership changes at the CBR, so
    any difference between the two groups puts the best breakpoint at zero. It is
    kept because a result *away* from zero would be informative - it would mean
    the fitted discontinuity is not where the CBR is.
    """
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    lo, hi = np.percentile(xs, [5, 95])
    best = None
    for tau in np.linspace(lo, hi, n_grid):
        l, r = ys[xs < tau], ys[xs >= tau]
        if l.size < 5 or r.size < 5:
            continue
        sse = ((l - l.mean()) ** 2).sum() + ((r - r.mean()) ** 2).sum()
        if best is None or sse < best['sse']:
            best = {'tau_kb': float(tau), 'level_left': float(l.mean()),
                    'level_right': float(r.mean()), 'sse': float(sse)}
    return best


def aic(sse, n, k):
    return n * np.log(sse / n) + 2 * k


def relax_k(logw_ret, logw_eli):
    """k = mean(log omega eliminated) / mean(log omega retained)."""
    mu_r, mu_e = float(np.mean(logw_ret)), float(np.mean(logw_eli))
    if abs(mu_r) < 1e-9:
        return np.nan, mu_r, mu_e
    return mu_e / mu_r, mu_r, mu_e


# ---------------------------------------------------------------------------
# Gene accounting
# ---------------------------------------------------------------------------

REGIONS = ('retained', 'eliminated')

#: Filter categories in the order main() applies the filters. Every gene of the
#: table falls in exactly one of them.
FILTER_COLS = ('not_finite', 'zero', 'below_min', 'above_max', 'usable')


def omega_filter_counts(d, omega_col, region_col, omega_min, omega_max):
    """Genes removed by each dN/dS filter, per region.

    not_finite  no estimate, or an infinite one (dS = 0)
    zero        dN/dS = 0: no non-synonymous change, and no logarithm
    below_min   0 < dN/dS < omega_min
    above_max   dN/dS > omega_max: saturated dS gives yn00 artefacts
    usable      the genes every fit and both figures use

    Mirrors the two filtering lines of main() exactly (finite and > 0, then
    inside [omega_min, omega_max]); main() checks that `usable` equals what
    survives them, so the S3B legend counts can be quoted from this table.
    Returns a DataFrame indexed by region: in_table, then FILTER_COLS.
    """
    w = pd.to_numeric(d[omega_col], errors='coerce').to_numpy(float)
    finite = np.isfinite(w)
    cat = np.full(w.shape, 'usable', dtype=object)
    cat[~finite] = 'not_finite'
    wf = np.where(finite, w, 0.0)                 # no comparisons on NaN
    cat[finite & (wf <= 0)] = 'zero'
    cat[finite & (wf > 0) & (wf < omega_min)] = 'below_min'
    cat[finite & (wf > 0) & (wf >= omega_min) & (wf > omega_max)] = 'above_max'
    reg = d[region_col].to_numpy()
    tab = pd.DataFrame(0, index=list(REGIONS), columns=list(FILTER_COLS),
                       dtype=int)
    for r in REGIONS:
        for c in FILTER_COLS:
            tab.loc[r, c] = int(((reg == r) & (cat == c)).sum())
    tab.insert(0, 'in_table', [int((reg == r).sum()) for r in REGIONS])
    parts = tab[list(FILTER_COLS)].sum(axis=1)
    if not (parts == tab['in_table']).all():
        raise RuntimeError(f'filter categories do not add up:\n{tab}')
    return tab


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

#: The words of the log10(dN/dS) label, set first as plain text so the layout
#: reserves its room; _set_log_omega_label() then redraws it with the 10 as a
#: drawn subscript. No Unicode subscript digits (Arial has none) and no
#: mathtext, so it prints in the label typeface with dN/dS upright.
LOG_OMEGA_LABEL = 'log10(dN/dS)'


def _log_omega_runs(fontsize):
    """log, a lowered 6 pt 10, then (dN/dS): runs for plain_text.draw_runs."""
    return subscript_runs('log', '10', '(dN/dS)', fontsize=fontsize)


def _set_log_omega_label(ax, axis):
    """Redraw the plain log10(dN/dS) label of `axis` ('x' or 'y') with a drawn
    subscript. Call once the layout is final (the label is measured)."""
    if axis == 'y':
        return replace_ylabel(ax, _log_omega_runs(ax.yaxis.label.get_fontsize()))
    return replace_xlabel(ax, _log_omega_runs(ax.xaxis.label.get_fontsize()))


def _axes(width, height_mm, left=15.0, right=3.0, bottom=13.0, top=4.0):
    """Explicit mm margins. constrained_layout shrinks the axes rather than fail."""
    fig, ax0 = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax0)
    w_in, h_in = fig.get_size_inches()
    ax = fig.add_axes([mm(left) / w_in, mm(bottom) / h_in,
                       (w_in - mm(left) - mm(right)) / w_in,
                       (h_in - mm(bottom) - mm(top)) / h_in])
    return fig, ax


def _finish(fig, outdir, stem):
    for problem in audit(fig, verbose=False):
        logging.warning('  audit: %s', problem)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths


def plot_slope_step(x, y, arms, step, boots, outdir, stem, window_kb,
                    n_bins=24, width='single', height_mm=65.0,
                    panel_letters=False, cap_pct=99.0):
    """Figure 3D: per-gene points, binned means, two fitted arms, the step.

    Binned summaries here are **mean +/- SE of log10(omega)**, not the median and
    bootstrap ribbon used in Figure 3A. That is deliberate and not an
    inconsistency: 3A describes a distribution, this panel displays a fitted
    linear model, and the quantity the model estimates is the conditional mean of
    log omega. Showing medians beside a least-squares line would draw the eye to a
    summary the line is not fitting.

    The region labels, with the n of each arm's fit, sit above the axes: inside,
    they sat on the point cloud. The 5 mm top margin holds them (65 mm keeps the
    45 mm axes of the 64 mm version); a panel letter needs 7.5 mm.
    """
    top = 7.5 if panel_letters else 5.0
    fig, ax = _axes(width, height_mm, bottom=15.0, top=top)
    xlo, xhi = -window_kb, window_kb
    ax.axvspan(xlo, 0, color=REGION_SHADE['retained'], zorder=0, lw=0)
    ax.axvspan(0, xhi, color=REGION_SHADE['eliminated'], zorder=0, lw=0)

    keep = np.abs(x) <= window_kb
    ax.scatter(x[keep], y[keep], s=1.6, linewidths=0, zorder=1, alpha=0.30,
               rasterized=True,
               c=[REGION['retained'] if v < 0 else REGION['eliminated']
                  for v in x[keep]])

    edges = np.linspace(xlo, xhi, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    mu = np.full(n_bins, np.nan)
    se = np.full(n_bins, np.nan)
    for j in range(n_bins):
        m = keep & (x >= edges[j]) & (x < edges[j + 1])
        if m.sum() >= 3:
            mu[j] = y[m].mean()
            se[j] = y[m].std(ddof=1) / np.sqrt(m.sum())
    ok = np.isfinite(mu)
    ax.errorbar(centres[ok], mu[ok], yerr=se[ok], fmt='o', ms=2.0,
                color=NEUTRAL['ink'], ecolor=NEUTRAL['ink'], elinewidth=0.6,
                capsize=1.2, lw=0, zorder=3)

    for name, xr in (('retained', (xlo, 0)), ('eliminated', (0, xhi))):
        a = arms[name]
        if a is None:
            continue
        xs = np.linspace(*xr, 60)
        ax.plot(xs, a['intercept'] + a['slope_per_100kb'] / 100.0 * xs,
                color=REGION[name], lw=1.4, zorder=5)

    ax.axvline(0, color=NEUTRAL['cbr'], lw=1.0, zorder=4)

    if step is not None:
        iR, iE = arms['retained']['intercept'], arms['eliminated']['intercept']
        ax.annotate('', xy=(0, iE), xytext=(0, iR), zorder=6,
                    arrowprops=dict(arrowstyle='<->', color=NEUTRAL['ink'],
                                    lw=0.9, shrinkA=0, shrinkB=0))
        lo, hi = ci95(boots.get('step', []))
        txt = f'step {10 ** step:.2f}×'           # Unicode times sign
        if np.isfinite(lo):
            txt += f' [{10 ** lo:.2f}, {10 ** hi:.2f}]'
        # Bottom left, not beside the arrow. The arrow sits between the two
        # fitted intercepts, which is exactly where the eliminated regression
        # line and its binned points are, so a label offset from it lands on the
        # data every time. There is only one arrow; the reader connects them.
        ax.annotate(txt, xy=(0.02, 0.04), xycoords='axes fraction',
                    ha='left', va='bottom', zorder=6,
                    fontsize=plt.rcParams['xtick.labelsize'])

    # Cap the y-axis on percentiles of the drawn points. Omega spans three
    # decades and the extreme low tail otherwise sets the range, pushing every
    # fitted line and the step into the top fifth of the panel. Every fit uses
    # the complete data; this bounds only what is displayed.
    yin = y[keep]
    ylo = float(np.percentile(yin, 100.0 - cap_pct))
    yhi = float(np.percentile(yin, cap_pct))
    pad = 0.16 * (yhi - ylo)
    ax.set_ylim(ylo - pad, yhi + pad)

    ax.set_xlim(xlo, xhi)
    ax.set_xlabel('Distance from nearest CBR (kb)')
    ax.set_ylabel(LOG_OMEGA_LABEL)
    # One label per arm, centred over it, above the top spine. Upright: italic
    # is for species names.
    for name, xc in (('retained', xlo / 2.0), ('eliminated', xhi / 2.0)):
        n = arms[name]['n'] if arms[name] else 0
        ax.annotate(f'{name.capitalize()} (n = {n:,})', xy=(xc, 1.0),
                    xycoords=('data', 'axes fraction'), xytext=(0, 2.0),
                    textcoords='offset points', ha='center', va='bottom',
                    color=REGION[name], fontsize=FONT['legend'],
                    annotation_clip=False)
    if panel_letters:
        panel_letter(ax, 'D')
    _set_log_omega_label(ax, 'y')             # the layout is final here
    return _finish(fig, outdir, stem)


def _place_neutral_label(fig, ax, means, dens_max, offset_pt=2.5):
    """Label the omega = 1 rule at its top, clear of every line and density.

    The label goes right of the rule. When it would run past the right edge (no
    gene sits that close to neutrality on that side), the x-range is extended
    rather than the label moved: the extension shows only empty density axis.
    It goes left of the rule only if a mean line would pass under it on the
    right. Finally the y-range is raised just enough that both densities stay
    below the label. Returns what was done, for the summary and provenance.
    """
    renderer = fig.canvas.get_renderer()
    pad_px = offset_pt * fig.dpi / 72.0
    info = {'side': 'right', 'x_extended_to': None, 'y_top_raised_to': None,
            'mean_under_label': False}

    def put(side):
        sign = 1.0 if side == 'right' else -1.0
        return ax.annotate('neutral', xy=(0.0, 0.98),
                           xycoords=('data', 'axes fraction'),
                           xytext=(sign * offset_pt, 0.0),
                           textcoords='offset points', va='top',
                           ha='left' if side == 'right' else 'right',
                           color=NEUTRAL['rule'], fontsize=FONT['annotation'],
                           zorder=6)

    def extent(label):
        """(window bbox, data x0, data x1) of a label, after a fresh draw."""
        fig.canvas.draw()
        bb = label.get_window_extent(renderer)
        inv = ax.transData.inverted()
        return bb, inv.transform((bb.x0, bb.y0))[0], inv.transform((bb.x1, bb.y0))[0]

    def blocked(x0, x1):
        return any(x0 <= m <= x1 for m in means)

    label = put('right')
    bb, x0, x1 = extent(label)
    axbb = ax.get_window_extent(renderer)
    if bb.x1 > axbb.x1 - pad_px:
        lo, hi = ax.get_xlim()
        x_rule = ax.transData.transform((0.0, 0.0))[0]
        need = (bb.x1 - x_rule) + pad_px            # px right of the rule
        if lo < 0.0 and need < axbb.width:
            # (0 - lo) / (hi' - lo) * W + need <= W, solved for hi'
            hi_new = lo + (0.0 - lo) * axbb.width / (axbb.width - need)
            ax.set_xlim(lo, hi_new)
            info['x_extended_to'] = float(hi_new)
            bb, x0, x1 = extent(label)
    if blocked(x0, x1):
        alt = put('left')
        bb_l, x0_l, x1_l = extent(alt)
        if not blocked(x0_l, x1_l) and bb_l.x0 >= axbb.x0:
            label.remove()
            label, bb, info['side'] = alt, bb_l, 'left'
        else:
            alt.remove()
            info['mean_under_label'] = True
            logging.warning('  S3B: a mean line passes under the "neutral" '
                            'label on both sides of the rule')

    # Raise the top until both densities end below the label. The label is
    # placed in axes fraction, so its own position does not move with ylim.
    fig.canvas.draw()
    bb = label.get_window_extent(renderer)
    bottom_frac = ax.transAxes.inverted().transform((bb.x0, bb.y0))[1]
    y0, y1 = ax.get_ylim()
    need_top = y0 + (dens_max - y0) / (bottom_frac - 0.03)
    if need_top > y1:
        ax.set_ylim(y0, need_top)
        info['y_top_raised_to'] = float(need_top)
    return info


def plot_relaxation_k(logw_ret, logw_eli, kstat, outdir, stem,
                      width='single', height_mm=63.0, panel_letters=False):
    """Supplementary fig. S3B: the two log-omega distributions and k.

    Inside the axes there is only data, the omega = 1 rule and its label. The
    legend (left) and k (right, one line beside each legend row) sit above the
    axes, so no line crosses them: the mean lines and the rule run the full
    height. The 9 mm top margin holds the two rows; 63 mm keeps the 41 mm axes
    of the 58 mm version, whose top margin was 4 mm.

    Returns (paths, layout), layout recording any extension of the axis ranges.
    """
    fig, ax = _axes(width, height_mm, top=9.0)
    grid = np.linspace(min(logw_ret.min(), logw_eli.min()),
                       max(logw_ret.max(), logw_eli.max()), 400)
    dens_max, means = 0.0, []
    for arr, name in ((logw_ret, 'retained'), (logw_eli, 'eliminated')):
        dens = gaussian_kde(arr)(grid)
        dens_max = max(dens_max, float(dens.max()))
        ax.fill_between(grid, dens, color=REGION[name], alpha=0.30, lw=0,
                        zorder=2)
        ax.plot(grid, dens, color=REGION[name], lw=1.2, zorder=3,
                label=f'{name.capitalize()} (n = {len(arr):,})')
        ax.axvline(arr.mean(), color=REGION[name], ls=(0, (3, 2)), lw=0.9,
                   zorder=4)
        means.append(float(arr.mean()))
    # omega = 1 is neutrality, i.e. log10 omega = 0.
    ax.axvline(0.0, color=NEUTRAL['rule'], lw=0.8, zorder=1)
    ax.set_xlabel(LOG_OMEGA_LABEL)
    ax.set_ylabel('Density')

    # Legend above the axes, left; k above the axes, right, each line level
    # with a legend row (measured after a draw, not assumed from font metrics).
    leg = ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=1,
                    frameon=False, borderpad=0.0, borderaxespad=0.3,
                    handlelength=1.4, handletextpad=0.5,
                    fontsize=FONT['legend'])
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_axes = ax.transAxes.inverted()
    k_lines = ((f"k = {kstat['k']:.3f} "
                f"[{kstat['ci'][0]:.3f}, {kstat['ci'][1]:.3f}]")
               .replace('-', '−'),           # true minus, if any
               'k < 1: relaxation')
    for row, line in zip(leg.get_texts(), k_lines):
        bb = row.get_window_extent(renderer)
        y_row = to_axes.transform((0.0, 0.5 * (bb.y0 + bb.y1)))[1]
        ax.text(1.0, y_row, line, transform=ax.transAxes, ha='right',
                va='center', fontsize=FONT['legend'], clip_on=False)

    layout = _place_neutral_label(fig, ax, means, dens_max)
    layout['xlim'] = [float(v) for v in ax.get_xlim()]
    layout['ylim'] = [float(v) for v in ax.get_ylim()]

    # The header must not run into itself or off the figure.
    fig.canvas.draw()
    texts = [t.get_window_extent(renderer) for t in ax.texts]
    leg_bb = leg.get_window_extent(renderer)
    k_left = min(bb.x0 for bb, t in zip(texts, ax.texts)
                 if t.get_text() in k_lines)
    if leg_bb.x1 > k_left:
        logging.warning('  S3B: the legend runs into the k annotation '
                        '(%.1f mm overlap)',
                        (leg_bb.x1 - k_left) / fig.dpi * 25.4)
    top_px = fig.bbox.y1
    if max([leg_bb.y1] + [bb.y1 for bb in texts]) > top_px:
        logging.warning('  S3B: the header runs off the top of the figure')
    if panel_letters:
        panel_letter(ax, 'B')
    _set_log_omega_label(ax, 'x')             # after every range extension
    return _finish(fig, outdir, stem), layout


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Per-arm slopes and the step at the CBR (Figure 3D), and the '
                    'relaxation exponent k (supplementary fig. S3B, written to '
                    'S8B_relaxation_k.*).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', required=True,
                   help='kaks_annotated.tsv from annotate_regions.py. Must carry '
                        'signed_cbr_dist1.')
    p.add_argument('--eliminated-bed', default=None,
                   help='Eliminated regions BED, for the cross-check against the '
                        'table. Strongly recommended.')
    p.add_argument('--chrom-sizes', default=None)
    p.add_argument('--terminal-tol', type=int, default=1000,
                   help='Must match annotate_regions.py')
    p.add_argument('--check-tol-bp', type=float, default=10.0)
    p.add_argument('--dist-col', default='signed_cbr_dist1')
    p.add_argument('--omega-col', default='Ka_Ks')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--slope-window', type=float, default=150.0,
                   help='Window each side of the CBR for the arm fits, kb. The '
                        'drawn axis stops here too: showing distances the model '
                        'was not fitted on invites them to be read as fitted.')
    p.add_argument('--omega-min', type=float, default=1e-4)
    p.add_argument('--omega-max', type=float, default=10.0,
                   help='Drop omega above this: saturated dS gives yn00 artifacts')
    p.add_argument('--cap-pct', type=float, default=99.0,
                   help='Percentile bounding the drawn y-range of Figure 3D. '
                        'Every fit uses the complete data.')
    p.add_argument('--n-bins', type=int, default=24)
    p.add_argument('--n-boot', type=int, default=2000)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--species', default='',
                   help='Recorded in the summary. No longer drawn as a title: '
                        'the species goes in the caption.')
    p.add_argument('--width', choices=['single', 'onehalf', 'double'],
                   default='single')
    p.add_argument('--panel-letters', action='store_true')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'cbr_changepoint_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'cbr_changepoint')
    # matplotlib's PDF writer subsets fonts through fontTools, which logs every
    # table it prunes at INFO -- several hundred lines per figure, burying the
    # lines this script writes (as cbr_meta_profile.py). Warnings still show.
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.info('cbr_changepoint_relax.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()
    rng = np.random.default_rng(args.seed)

    df = pd.read_csv(args.annotated, sep='\t')
    logging.info('Loaded %s rows from %s', f'{len(df):,}',
                 os.path.basename(args.annotated))
    for col in (args.dist_col, args.omega_col, args.region_col):
        if col not in df.columns:
            sys.exit(f'ERROR: {col!r} is not a column of {args.annotated}.\n'
                     f'  If this is signed_cbr_dist1, the table predates '
                     f'annotate_regions.py -- rerun it and repoint '
                     f'ANNOTATED_TABLE.\n  columns: {list(df.columns)}')

    # --- cross-check ------------------------------------------------------
    check_diag = None
    if args.eliminated_bed:
        logging.info('Cross-checking %s against a fresh computation', args.dist_col)
        elim = parse_bed(args.eliminated_bed)
        merged = merge_intervals(elim)
        sizes, inferred = load_chrom_sizes(args.chrom_sizes, elim)
        cbrs, cbr_diag = derive_cbrs(merged, sizes, args.terminal_tol)
        recomputed = np.full(len(df), np.nan)
        have = df['chrom1'].notna() & df['midpoint1'].notna()
        for chrom, block in df[have].groupby('chrom1', sort=False):
            sd, _ = signed_distance(chrom, block['midpoint1'].to_numpy(float),
                                    cbrs, merged)
            recomputed[df.index.get_indexer(block.index)] = sd
        both = ~np.isnan(recomputed) & df[args.dist_col].notna().to_numpy()
        delta = np.abs(recomputed[both] - df[args.dist_col].to_numpy(float)[both])
        n_bad = int((delta > args.check_tol_bp).sum())
        if n_bad:
            sys.exit(f'ERROR: the table disagrees with a fresh computation of '
                     f'{args.dist_col} for {n_bad:,} of {int(both.sum()):,} genes '
                     f'(max |delta| {delta.max():,.0f} bp). The table was built '
                     f'with a different CBR definition. Rerun annotate_regions.py.')
        check_diag = {'n_compared': int(both.sum()), 'n_disagreeing': 0,
                      'max_abs_delta_bp': float(delta.max()) if delta.size else 0.0,
                      'cbrs': cbr_diag, 'chrom_sizes_inferred': inferred}
        logging.info('  cross-check: %s genes agree (max |delta| %.1f bp)',
                     f"{int(both.sum()):,}", check_diag['max_abs_delta_bp'])
    else:
        logging.warning('No --eliminated-bed: the table is taken on trust.')

    # --- prepare ----------------------------------------------------------
    d = df[df[args.region_col].isin(['retained', 'eliminated'])].copy()
    d[args.omega_col] = pd.to_numeric(d[args.omega_col], errors='coerce')
    n0 = len(d)
    # Count what each filter removes before applying the filters, so the
    # legends can say where their n come from.
    filt = omega_filter_counts(d, args.omega_col, args.region_col,
                               args.omega_min, args.omega_max)
    d = d[np.isfinite(d[args.omega_col]) & (d[args.omega_col] > 0)]
    d = d[(d[args.omega_col] >= args.omega_min) &
          (d[args.omega_col] <= args.omega_max)]
    d['logw'] = np.log10(d[args.omega_col])
    logging.info('  %s of %s rows usable (omega in [%g, %g])',
                 f'{len(d):,}', f'{n0:,}', args.omega_min, args.omega_max)
    for r in REGIONS:
        f = filt.loc[r]
        logging.info('    %-10s %s in table: %s not finite, %s zero, %s below '
                     '%g, %s above %g, %s usable', r, f'{f.in_table:,}',
                     f'{f.not_finite:,}', f'{f.zero:,}', f'{f.below_min:,}',
                     args.omega_min, f'{f.above_max:,}', args.omega_max,
                     f'{f.usable:,}')

    logw_ret = d.loc[d[args.region_col] == 'retained', 'logw'].to_numpy()
    logw_eli = d.loc[d[args.region_col] == 'eliminated', 'logw'].to_numpy()
    logging.info('  retained %s, eliminated %s', f'{len(logw_ret):,}',
                 f'{len(logw_eli):,}')
    if (filt.loc['retained', 'usable'], filt.loc['eliminated', 'usable']) != \
            (len(logw_ret), len(logw_eli)):
        sys.exit('ERROR: the filter accounting disagrees with the filtered '
                 'table; omega_filter_counts() no longer mirrors main().')

    spatial = d.dropna(subset=[args.dist_col])
    x = spatial[args.dist_col].to_numpy(float) / 1000.0     # kb
    y = spatial['logw'].to_numpy(float)
    in_win = int((np.abs(x) <= args.slope_window).sum())
    logging.info('  %s genes within +/-%g kb of a CBR', f'{in_win:,}',
                 args.slope_window)
    if in_win < 30:
        sys.exit('ERROR: fewer than 30 genes inside the slope window.')

    # Figure 3D accounting: usable genes -> a CBR on their sequence -> inside
    # the window. The arms are the sign of the distance, which must agree with
    # the region label (positive = eliminated side, by annotate_regions.py).
    reg_all = d[args.region_col].to_numpy()
    reg_sp = spatial[args.region_col].to_numpy()
    win = np.abs(x) <= args.slope_window
    acc4d = {
        'usable': {r: int((reg_all == r).sum()) for r in REGIONS},
        'no_cbr_on_sequence': {r: int(((reg_all == r) &
                                       d[args.dist_col].isna().to_numpy()).sum())
                               for r in REGIONS},
        'with_distance': {r: int((reg_sp == r).sum()) for r in REGIONS},
        'in_window_arm': {'retained': int((win & (x < 0)).sum()),
                          'eliminated': int((win & (x > 0)).sum())},
        'in_window_at_zero': int((win & (x == 0)).sum()),
        'sign_disagrees_with_region': int((((reg_sp == 'retained') & (x > 0)) |
                                           ((reg_sp == 'eliminated') & (x < 0)))
                                          .sum()),
    }
    if acc4d['sign_disagrees_with_region']:
        logging.warning('  %s genes have a distance sign that disagrees with '
                        'their region label; the arms follow the sign',
                        f"{acc4d['sign_disagrees_with_region']:,}")

    # --- Figure 3D fits ---------------------------------------------------
    arms, step = fit_arms(x, y, args.slope_window)
    boots = bootstrap_arms(x, y, args.slope_window, args.n_boot, rng)
    free = free_changepoint(x[np.abs(x) <= args.slope_window],
                            y[np.abs(x) <= args.slope_window])

    rows = []
    for name in ('retained', 'eliminated'):
        a = arms[name]
        if not a:
            continue
        lo, hi = ci95(boots.get(f'slope_{name}', []))
        rows.append({'quantity': f'slope_{name}', 'unit': 'log10 omega / 100 kb',
                     'estimate': a['slope_per_100kb'], 'lo95': lo, 'hi95': hi,
                     'p': a['p'], 'n': a['n']})
    step_se, step_p = step_analytic_p(arms)
    step_pboot = boot_p(boots.get('step', []), 0.0)
    if step is not None:
        lo, hi = ci95(boots.get('step', []))
        rows.append({'quantity': 'step_at_CBR', 'unit': 'log10 omega',
                     'estimate': step, 'lo95': lo, 'hi95': hi,
                     'se': step_se, 'p': step_p, 'p_bootstrap': step_pboot,
                     'n': in_win})
        rows.append({'quantity': 'step_at_CBR_fold', 'unit': 'fold',
                     'estimate': 10 ** step, 'lo95': 10 ** lo, 'hi95': 10 ** hi,
                     'se': np.nan, 'p': step_p, 'p_bootstrap': step_pboot,
                     'n': in_win})

    # --- S3B --------------------------------------------------------------
    k, mu_r, mu_e = relax_k(logw_ret, logw_eli)
    kv = []
    for _ in range(args.n_boot):
        rr = logw_ret[rng.integers(0, len(logw_ret), len(logw_ret))]
        ee = logw_eli[rng.integers(0, len(logw_eli), len(logw_eli))]
        val = relax_k(rr, ee)[0]
        if np.isfinite(val):
            kv.append(val)
    kci = ci95(kv)
    p_welch = float(ttest_ind(logw_ret, logw_eli, equal_var=False).pvalue)
    p_mwu = float(mannwhitneyu(logw_ret, logw_eli,
                               alternative='two-sided').pvalue)
    kstat = {'k': k, 'ci': kci, 'p_boot': boot_p(kv, 1.0), 'p_welch': p_welch,
             'p_mwu': p_mwu, 'mu_r': mu_r, 'mu_e': mu_e}
    rows.append({'quantity': 'relaxation_k', 'unit': 'dimensionless',
                 'estimate': k, 'lo95': kci[0], 'hi95': kci[1],
                 'p': p_welch, 'n': len(logw_ret) + len(logw_eli)})

    stats = pd.DataFrame(rows)
    # BH over the two slopes and the step: three tests asking whether constraint
    # varies with position. k is a different question and is not in the family.
    fam = stats['quantity'].isin(['slope_retained', 'slope_eliminated',
                                  'step_at_CBR'])
    stats['q'] = np.nan
    stats.loc[fam, 'q'] = benjamini_hochberg(stats.loc[fam, 'p'].to_numpy())
    stats.to_csv(os.path.join(args.outdir, 'cbr_changepoint_stats.tsv'),
                 sep='\t', index=False, float_format='%.6g')

    # --- figures ----------------------------------------------------------
    # Drawn after every statistic, and they draw no random numbers, so the
    # figures cannot change a number in the table.
    logging.info('Figure 3D')
    plot_slope_step(x, y, arms, step, boots, args.outdir, 'cbr_slope_step',
                    args.slope_window, n_bins=args.n_bins, width=args.width,
                    panel_letters=args.panel_letters, cap_pct=args.cap_pct)
    logging.info('Supplementary fig. S3B (S8B_relaxation_k)')
    _, s3b_layout = plot_relaxation_k(logw_ret, logw_eli, kstat, args.outdir,
                                      'S8B_relaxation_k', width=args.width,
                                      panel_letters=args.panel_letters)

    # --- summary ----------------------------------------------------------
    L = ['Constraint boundary at CBRs, and the relaxation exponent', '',
         f'Table               {os.path.basename(args.annotated)}']
    if args.species:
        L.append(f'Species             {args.species}')
    L += [f'Usable genes        {len(d):,} '
          f'(omega in [{args.omega_min:g}, {args.omega_max:g}])',
          f'Slope window        +/- {args.slope_window:g} kb '
          f'({in_win:,} genes)',
          f'Bootstrap           {args.n_boot:,} resamples, stratified by arm', '']
    if check_diag:
        L += [f"Cross-check         {check_diag['n_compared']:,} genes agree "
              f"(max |delta| {check_diag['max_abs_delta_bp']:.1f} bp)",
              f"CBRs                {check_diag['cbrs']['n_cbrs']:,} on "
              f"{check_diag['cbrs']['n_sequences_with_cbrs']} sequences", '']

    # Gene accounting, in the order the filters are applied.
    head = ('in table', 'not finite', 'dN/dS = 0', f'< {args.omega_min:g}',
            f'> {args.omega_max:g}', 'usable')
    L += ['Genes per region, by the filter that removes them', '',
          '  ' + f'{"":<12}' + ''.join(f'{h:>12}' for h in head)]
    for r in REGIONS:
        L.append('  ' + f'{r:<12}' + ''.join(
            f'{int(v):>12,}' for v in filt.loc[r, ['in_table'] +
                                              list(FILTER_COLS)]))
    fr, fe = filt.loc['retained'], filt.loc['eliminated']
    L += ['',
          '  not finite: no estimate, or infinite (dS = 0). dN/dS = 0 has no',
          '  logarithm. Above the upper bound, saturated dS gives yn00 artefacts.',
          '',
          f'  For the S3B legend: n = {fr.usable:,} retained and {fe.usable:,} '
          f'eliminated genes with',
          f'  {args.omega_min:g} <= dN/dS <= {args.omega_max:g}; excluded: '
          f'dN/dS = 0 ({fr.zero:,} and {fe.zero:,}),',
          f'  dN/dS > {args.omega_max:g} ({fr.above_max:,} and '
          f'{fe.above_max:,}), no finite estimate ({fr.not_finite:,} and '
          f'{fe.not_finite:,})'
          + (f', 0 < dN/dS < {args.omega_min:g} ({fr.below_min:,} and '
             f'{fe.below_min:,})' if fr.below_min + fe.below_min else '')
          + '.', '']
    a4 = acc4d

    def pair(label, counts):
        return (f"    {label:<48}retained {counts['retained']:,}, "
                f"eliminated {counts['eliminated']:,}")

    L += ['  Figure 3D, from the usable genes:',
          pair('on a sequence without a CBR (no distance)',
               a4['no_cbr_on_sequence']),
          pair('with a distance', a4['with_distance']),
          pair(f'inside +/- {args.slope_window:g} kb (the arms; the panel n)',
               a4['in_window_arm'])
          + (f" ({a4['in_window_at_zero']:,} at exactly 0, in neither arm)"
             if a4['in_window_at_zero'] else ''),
          f"    {'distance sign disagreeing with region':<48}"
          f"{a4['sign_disagrees_with_region']:,}", '']

    L += ['Figure 3D -- per-arm slopes and the step', '']
    for name in ('retained', 'eliminated'):
        a = arms[name]
        if not a:
            L.append(f'  {name:<12} too few genes')
            continue
        lo, hi = ci95(boots.get(f'slope_{name}', []))
        q = float(stats.loc[stats['quantity'] == f'slope_{name}', 'q'].iloc[0])
        L.append(f'  {name:<12} slope {a["slope_per_100kb"]:+.4f} '
                 f'[{lo:+.4f}, {hi:+.4f}] log10 omega per 100 kb, '
                 f'P = {a["p"]:.3g}, q = {q:.3g} ({stars(q)}), n = {a["n"]:,}')
        L.append(f'  {"":<12} boundary log10 omega = {a["intercept"]:.4f} '
                 f'(omega = {10 ** a["intercept"]:.4f})')
    if step is not None:
        lo, hi = ci95(boots.get('step', []))
        pq = float(stats.loc[stats['quantity'] == 'step_at_CBR', 'q'].iloc[0])
        L += ['',
              f'  step          {step:+.4f} [{lo:+.4f}, {hi:+.4f}] log10 units',
              f'                {10 ** step:.3f}x [{10 ** lo:.3f}, {10 ** hi:.3f}] '
              f'fold',
              f'                Wald P = {step_p:.3g}, q = {pq:.3g} '
              f'({stars(pq)});  bootstrap P = {step_pboot:.3g} '
              f'(floored at 1/{args.n_boot:,})']
    if free:
        L += ['',
              f"  free changepoint at {free['tau_kb']:+.1f} kb "
              f"(levels {free['level_left']:.3f} -> {free['level_right']:.3f})",
              '  CIRCULAR -- the sign of the distance is region membership, which',
              '  changes at the CBR, so this lands at zero whenever the groups',
              '  differ at all. Consistency check only; do not quote it.']

    L += ['', 'Supplementary fig. S3B (file S8B_relaxation_k) -- relaxation '
              'exponent', '',
          f'  mean log10 omega   retained {mu_r:.4f}, eliminated {mu_e:.4f}',
          f'  k = {k:.4f}  [{kci[0]:.4f}, {kci[1]:.4f}]',
          f'  k {"<" if k < 1 else ">"} 1  ==>  '
          f'{"relaxation" if k < 1 else "intensification"}',
          f'  bootstrap P(k = 1) = {kstat["p_boot"]:.3g}',
          f'  Welch on log omega P = {p_welch:.3g}   '
          f'(the exact null for k = 1)',
          f'  Mann-Whitney       P = {p_mwu:.3g}   '
          f'(distribution-free, tests dominance)', '',
          'Notes', '',
          '  This is NOT HyPhy RELAX. Canonical RELAX contrasts *branches*; PDE is',
          '  ancestral across all three lineages here, so the contrast is between',
          '  gene sets. The k relationship is the same, the labelling is not, and',
          '  the Methods must say so.', '',
          '  The step is the difference between the two arms\' fitted boundary',
          '  intercepts, so it is an extrapolation to d = 0 from each side rather',
          '  than a comparison of the genes nearest the CBR. That is the point:',
          '  flat slopes plus a large intercept difference is what makes the',
          '  boundary sharp rather than gradual.', '',
          '  Points are drawn only within the slope window. Showing distances the',
          '  model was not fitted on invites them to be read as fitted.', '',
          '  Binned summaries are the mean and SE of log10 omega, not the median',
          '  and bootstrap interval used in Figure 3A. The panel displays a',
          '  least-squares fit, and the conditional mean is the quantity that fit',
          '  estimates.', '',
          '  S3B draws Gaussian kernel densities of log10 omega (scipy default',
          '  bandwidth, Scott\'s rule) with each group\'s mean as a dashed line;',
          '  k is the ratio of those two means.']
    if s3b_layout.get('x_extended_to') is not None:
        L.append(f"  The S3B x-axis runs to {s3b_layout['x_extended_to']:.2f}, "
                 f"past the largest value, so the neutral label fits.")
    if s3b_layout.get('side') == 'left':
        L.append('  The neutral label sits left of the rule: a mean line would '
                 'pass under it on the right.')
    L.append('')
    summary = '\n'.join(L)
    with open(os.path.join(args.outdir, 'cbr_changepoint_summary.txt'), 'w') as fh:
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
                   'arms': arms, 'step': step, 'step_se': step_se,
                   'step_p_wald': step_p, 'step_p_bootstrap': step_pboot,
                   'step_ci': list(ci95(boots.get('step', []))),
                   'free_changepoint': free,
                   'relaxation': kstat,
                   'gene_filter_counts': {r: {c: int(v) for c, v in
                                              filt.loc[r].items()}
                                          for r in REGIONS},
                   'figure_4d_genes': acc4d,
                   's3b_layout': s3b_layout}, fh, indent=2, default=float)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
