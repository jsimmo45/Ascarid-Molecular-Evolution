#!/usr/bin/env python3
"""
dnds_decomposition.py - decompose the dN/dS contrast into its dN and dS parts.

Elevated dN/dS in eliminated genes has several possible causes, and they are
distinguishable by looking at the numerator and denominator separately:

    elevated dN, unchanged dS   relaxed purifying selection
    elevated dN and dS          higher mutation rate plus relaxed selection
    elevated dN, reduced dS     relaxed selection plus synonymous constraint
    unchanged dN, reduced dS    apparent inflation from a low denominator

Reads the annotated rate table from the region-annotation step and writes:

    dnds_decomposition.{svg,pdf,png}   dN, dS and dN/dS by region
    dnds_consensus.{svg,pdf,png}       the same, restricted to genes whose
                                       region call agrees between both species
    dn_vs_ds_scatter.{svg,pdf,png}     dN against dS, coloured by region
    decomposition_summary.txt          medians, fold changes, tests
    provenance.json                    inputs, arguments and versions


Notes on method
---------------
**Two region definitions.** The primary figure uses the *Ascaris* annotation,
which anchors the result to the most completely annotated genome. The consensus
figure keeps only ortholog pairs whose retained/eliminated call agrees between
*Ascaris* and *Parascaris*, which controls for misannotation in either genome at
the cost of some genes. The consensus contrast is the larger of the two, so the
single-species version reported in the manuscript is the conservative one.

**Display capping does not affect the statistics.** Points above the 99th
percentile of the pooled values are left out of the jitter so the panels are
readable, but the boxes, the whiskers and every Mann-Whitney test use the
complete distribution. The drawing limit is raised to the higher upper whisker
end if that lies above the percentile, so no whisker runs past the points. The
summary counts the points not drawn, per panel and group, so nothing is cut
silently.

**The scatter is clipped, and says so.** Restricting the axes to the dense region
hides a small number of pairs; the panel labels carry the full n and the number
clipped, rather than quietly reporting the surviving count as though it were the
sample size.

Changes 2026-10-12 (figure review of Fig 2A and supplementary fig. S1)
---------------------------------------------------------------------
* Each bracket sits above the highest drawn point of its panel. It used to sit
  at 1.06 x the 98th percentile, inside the top of the cloud, and in the dN/dS
  panel it crossed the eliminated whisker.
* The summary gives, per panel and group, how many points above the display cap
  are not drawn (the log had the total only), so the legend can quote them.
* The scatter's off-axis count moves from inside the axes, where the reference
  lines ran through it, to a line under the legend. The figure is 4 mm taller
  (76 mm) to hold it; the axes box is unchanged.
No statistic changes: tests, medians and folds use every value as before.

Changes 2026-10-14 (boxes from every gene, as in Fig 2D)
--------------------------------------------------------
* Boxes and whiskers are computed from every gene (box_stats) and drawn with
  ax.bxp. ax.boxplot computed them from the drawn values, after the points
  above the 99th percentile of the pooled values had been removed. Most of
  those points are eliminated genes, so the eliminated quartiles, median line
  and whisker sat below their true values, and the drawn median differed from
  the median in the summary.
* The drawing limit is the 99th percentile or the higher upper whisker end,
  whichever is higher. Points above it are still not drawn, and still counted.
* The summary lists each group's quartiles, whisker ends and the number of
  genes beyond the whiskers, for the legend.
No statistic changes: tests, medians and folds use every value as before.

Usage
-----
    python dnds_decomposition.py \\
      --annotated kaks_regions/kaks_annotated.tsv \\
      --outdir results/dnds_decomposition_20260815

Input needs columns for the nonsynonymous rate, the synonymous rate and their
ratio, under either the dN/dS or the legacy Ka/Ks names, plus `region1` and
`consensus_region`.

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
from matplotlib.colors import to_rgba
from scipy.stats import mannwhitneyu

_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, 'lib'), os.path.join(os.path.dirname(_HERE), 'lib')):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, _cand)

from paper_colors import (COL, FONT, LINEWIDTH, NEUTRAL, REGION, apply_style,
                          audit, mm, panel_letter, save)
from rate_tables import read_rate_table

apply_style()

EXTS = ('.svg', '.pdf', '.png')

#: column, axis label. Order is the order of the sub-panels, left to right.
METRICS = [
    ('dN',    'dN (nonsynonymous)'),
    ('dS',    'dS (synonymous)'),
    ('dN_dS', 'dN/dS'),
]

#: Box order, left to right: retained blue/left, eliminated red/right.
REGIONS = ('retained', 'eliminated')

#: The box_stats keys ax.bxp() reads.
BXP_KEYS = ('med', 'q1', 'q3', 'whislo', 'whishi', 'fliers')


def box_stats(values, whis=1.5):
    """Box and whisker statistics for one group, in the form ``ax.bxp`` draws.

    Computed from every value passed, never from the values drawn. Median and
    quartiles are ``np.percentile`` with linear interpolation, which is what
    pandas ``median()`` gives in the summary, so the box on the page is the box
    in the table. Whiskers follow Tukey, matplotlib's default rule: each ends at
    the most extreme value within ``whis`` x IQR of the box, or at the box edge
    if no value lies beyond it. The same function as gene_structure.py (Figure
    2D) and codon_usage_bias.py (Figure 2B), without the log-axis option, which
    these linear axes do not need.
    """
    v = np.asarray(values, dtype=float)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    lo_lim = q1 - whis * (q3 - q1)
    hi_lim = q3 + whis * (q3 - q1)
    up = v[v <= hi_lim]
    whishi = float(up.max()) if len(up) and up.max() >= q3 else float(q3)
    down = v[v >= lo_lim]
    whislo = float(down.min()) if len(down) and down.min() <= q1 else float(q1)
    return {'med': float(med), 'q1': float(q1), 'q3': float(q3),
            'whislo': whislo, 'whishi': whishi,
            # Read unconditionally by bxp() in older matplotlib releases.
            'fliers': np.array([]),
            'beyond_whiskers': int(((v < whislo) | (v > whishi)).sum())}


def stars(p):
    """Significance marker. Figure 2 uses markers throughout; exact values go in
    the summary and the legend, so the panels stay uncluttered and one
    convention holds across the whole figure."""
    for cut, mark in ((1e-4, '****'), (1e-3, '***'), (1e-2, '**'), (5e-2, '*')):
        if p < cut:
            return mark
    return 'ns'


def bracket(ax, x1, x2, y, p, height=0.045):
    """Significance bracket spanning two box positions."""
    lo, hi = ax.get_ylim()
    tick = (hi - lo) * height * 0.35
    ax.plot([x1, x1, x2, x2], [y, y + tick, y + tick, y],
            color=NEUTRAL['ink'], linewidth=LINEWIDTH['thin'],
            clip_on=False, solid_joinstyle='miter')
    ax.text((x1 + x2) / 2, y + tick * 1.15, stars(p), ha='center', va='bottom',
            fontsize=FONT['annotation'], fontweight='bold', color=NEUTRAL['ink'])


def save_figure(fig, outdir, stem):
    for problem in audit(fig, verbose=False):
        logging.warning(f"FIGURE AUDIT [{stem}]: {problem}")
    # save() suppresses the per-run timestamp each format would otherwise
    # embed, so rerunning on unchanged input reproduces the file byte for byte.
    save(fig, outdir, stem, formats=[e.lstrip('.') for e in EXTS])
    plt.close(fig)
    logging.info(f"Saved: {stem}{'/'.join(EXTS)}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_decomposition(df, region_col, outdir, stem, rng, width='double',
                       height_mm=58.0, letter=None, cap_pct=99.0,
                       max_points=1500, box_alpha=0.40):
    """dN, dS and dN/dS side by side for retained versus eliminated genes.

    One function serves both the single-species and the consensus versions.
    They were separate copies of the same sixty lines, which is how two panels
    of one figure end up with different conventions.

    The test, the boxes and the whiskers use every gene. Points above the
    drawing limit (the ``cap_pct``-th percentile of the pooled values, raised
    to the higher upper whisker end if that is above it) are left out of the
    jitter only, and counted per group in the returned results.
    """
    fig, axes = plt.subplots(1, len(METRICS), figsize=(COL[width], mm(height_mm)))
    results = {}

    for ax, (col, ylabel) in zip(axes, METRICS):
        ret = df.loc[df[region_col] == 'retained', col].dropna()
        elim = df.loc[df[region_col] == 'eliminated', col].dropna()
        if len(ret) < 5 or len(elim) < 5:
            ax.text(0.5, 0.5, 'insufficient data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=FONT['annotation'])
            ax.set_axis_off()
            continue

        # Tests use every value; the cap only governs which points are drawn.
        _, p = mannwhitneyu(ret, elim, alternative='two-sided')
        results[col] = dict(n_ret=len(ret), n_elim=len(elim),
                            med_ret=float(ret.median()),
                            med_elim=float(elim.median()),
                            fold=float(elim.median() / ret.median()), p=float(p))

        # Box and whiskers from every gene (box_stats), never from the points
        # drawn. ax.boxplot() computed them from the capped values, and most of
        # the capped points are eliminated genes, so the eliminated box sat low.
        box = {key: box_stats(vals.to_numpy(dtype=float))
               for key, vals in zip(REGIONS, (ret, elim))}

        # Drawing limit: the cap_pct-th percentile of the pooled values, raised
        # to the higher upper whisker end if that lies above it, so that no
        # whisker runs past the drawn points. Rates are >= 0 and the axis starts
        # at 0, so there is no lower limit.
        cap = float(np.percentile(pd.concat([ret, elim]), cap_pct))
        limit = max(cap, *(b['whishi'] for b in box.values()))
        shown = [ret[ret <= limit], elim[elim <= limit]]
        results[col]['n_clipped_ret'] = int(len(ret) - len(shown[0]))
        results[col]['n_clipped_elim'] = int(len(elim) - len(shown[1]))
        results[col]['n_clipped'] = (results[col]['n_clipped_ret'] +
                                     results[col]['n_clipped_elim'])
        results[col]['cap'] = cap
        results[col]['limit'] = float(limit)
        results[col]['box'] = {key: {k: b[k] for k in ('q1', 'med', 'q3', 'whislo',
                                                       'whishi', 'beyond_whiskers')}
                               for key, b in box.items()}

        # bxp() draws the precomputed statistics; boxplot() would recompute them
        # from whatever it is handed, which is how the capped values set them.
        bp = ax.bxp([{k: box[key][k] for k in BXP_KEYS} for key in REGIONS],
                    positions=[0, 1], widths=0.5, patch_artist=True,
                    showfliers=False, zorder=2)
        for box_patch, key in zip(bp['boxes'], REGIONS):
            # Alpha in the face colour, not via set_alpha(): patch-level alpha
            # fades the edge too, so a black outline prints mid-grey. The fill is
            # kept light because the jitter is drawn on top of it -- an opaque box
            # hides every point in the interquartile range, which is most of them.
            box_patch.set_facecolor(to_rgba(REGION[key], box_alpha))
            box_patch.set_alpha(None)
            box_patch.set_edgecolor('black')
            box_patch.set_linewidth(LINEWIDTH['thin'])
        # Summary statistics back above the cloud.
        for part, z in (('whiskers', 4), ('caps', 4), ('medians', 5)):
            for artist in bp[part]:
                artist.set_color('black')
                artist.set_zorder(z)
                artist.set_linewidth(LINEWIDTH['data'] if part == 'medians'
                                     else LINEWIDTH['thin'])

        for i, (vals, key) in enumerate(zip(shown, REGIONS)):
            v = vals.to_numpy()
            if len(v) > max_points:
                v = v[rng.choice(len(v), max_points, replace=False)]
            ax.scatter(np.full(len(v), i) + rng.normal(0, 0.06, len(v)), v,
                       s=1.5, alpha=0.28, color=REGION[key], edgecolors='none',
                       zorder=3, rasterized=True)

        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'Retained\n(n = {len(ret):,})',
                            f'Eliminated\n(n = {len(elim):,})'])
        ax.set_ylabel(ylabel)
        # Bracket above the highest drawn point, so it never crosses a point or
        # a whisker (each whisker end is a gene's value, so it is among the
        # drawn points); the axis is extended to hold it and its label.
        top = float(max(s.max() for s in shown))
        y_br = top * 1.04
        ax.set_ylim(0, top * 1.19)
        bracket(ax, 0, 1, y_br, p)

    if letter:
        panel_letter(axes[0], letter)
    save_figure(fig, outdir, stem)
    return results


def plot_dn_vs_ds(df, outdir, stem, width='single', height_mm=76.0,
                  ds_max=1.5, dn_max=0.6, letter=None):
    """dN against dS, coloured by region.

    Axis limits and the point filter are the same numbers, so nothing is counted
    that cannot be seen. The off-axis counts are annotated in the corner rather
    than folded into the legend labels, which keeps the legend narrow enough to
    fit the column.

    Layout is set explicitly in millimetres instead of by constrained_layout.
    The reference lines are labelled where they leave the axes, so the legend
    carries two entries rather than five: at 89 mm, five entries of legend text
    are wider than the figure, and constrained_layout responds to that by
    shrinking the axes towards zero rather than by reporting a problem.
    """
    fig = plt.figure(figsize=(COL[width], mm(height_mm)))
    fig.set_layout_engine('none')

    # Margins in millimetres, converted to figure fractions. Fixed numbers mean
    # the axes box cannot be renegotiated by anything drawn later.
    w_mm, h_mm = COL[width] * 25.4, height_mm
    # The bottom margin holds tick labels, the axis label, the legend row and
    # one line for the off-axis count.
    left, right, bottom, top = 12.5, 3.0, height_mm - 52.0, 2.5
    ax = fig.add_axes([left / w_mm, bottom / h_mm,
                       (w_mm - left - right) / w_mm,
                       (h_mm - bottom - top) / h_mm])

    # Eliminated drawn last: 322 points under 9,085 are otherwise invisible.
    off_counts = []
    for key, label, z in (('retained', 'Retained', 2),
                          ('eliminated', 'Eliminated', 3)):
        sub = df.loc[df['region1'] == key].dropna(subset=['dN', 'dS'])
        vis = sub[(sub['dN'] <= dn_max) & (sub['dS'] <= ds_max)]
        off = len(sub) - len(vis)
        if off:
            off_counts.append(f'{label.lower()} {off}')
        ax.scatter(vis['dS'], vis['dN'], s=1.5, alpha=0.30, color=REGION[key],
                   edgecolors='none', zorder=z, rasterized=True,
                   label=f'{label} (n = {len(sub):,})')

    ax.set_xlim(0, ds_max)
    ax.set_ylim(0, dn_max)
    ax.set_xlabel('dS (synonymous)')
    ax.set_ylabel('dN (nonsynonymous)')

    # --- reference lines, labelled in place ---------------------------------
    # Rotation has to match the line's slope on the page, not in the data, so it
    # is computed from the axes box in inches. Both are known here because the
    # box was set above in millimetres.
    ax_w_in = (w_mm - left - right) / 25.4
    ax_h_in = (h_mm - bottom - top) / 25.4
    for omega, ls, lw, colour, text in (
            (1.0, '--', LINEWIDTH['rule'], NEUTRAL['rule'], 'dN = dS'),
            (0.5, ':', LINEWIDTH['thin'], NEUTRAL['ink_light'], 'dN/dS = 0.5'),
            (0.2, '-.', LINEWIDTH['thin'], NEUTRAL['ink_light'], 'dN/dS = 0.2')):
        ax.plot([0, ds_max], [0, ds_max * omega], linestyle=ls, color=colour,
                linewidth=lw, zorder=4)
        angle = np.degrees(np.arctan2(omega * ds_max * ax_h_in,
                                      dn_max * ax_w_in))
        # Where the line leaves the panel: through the top if it reaches dn_max
        # inside the x range, otherwise through the right edge.
        x_exit = min(ds_max, dn_max / omega)
        x_lab = x_exit * 0.62
        ax.text(x_lab, omega * x_lab, text + ' ', rotation=angle,
                rotation_mode='anchor', transform_rotates_text=True,
                ha='right', va='bottom', fontsize=FONT['annotation'],
                color=NEUTRAL['ink'], zorder=6,
                bbox=dict(facecolor='white', edgecolor='none',
                          alpha=0.75, pad=0.8))

    if off_counts:
        # Outside the axes, under the legend: inside, the reference lines and the
        # sparse tail of the cloud ran through it wherever it was put.
        fig.text(left / w_mm, 6.5 / h_mm,
                 'Not drawn (outside the axes): ' + ', '.join(off_counts),
                 ha='left', va='top', fontsize=FONT['annotation'],
                 color=NEUTRAL['ink'])

    # Two entries in one row: ~55 mm of the 73.5 mm axes width.
    leg = ax.legend(loc='upper left', bbox_to_anchor=(0.0, -0.22), ncol=2,
                    frameon=False, fontsize=FONT['legend'], handlelength=1.0,
                    columnspacing=1.0, handletextpad=0.5, borderaxespad=0.0,
                    markerscale=4)
    # The points are drawn at alpha 0.30 so the cloud reads; a swatch at 0.30 is
    # too pale to key against. matplotlib renamed this attribute in 3.7.
    for handle in (getattr(leg, 'legend_handles', None)
                   or getattr(leg, 'legendHandles', [])):
        handle.set_alpha(1.0)
    if letter:
        panel_letter(ax, letter)
    save_figure(fig, outdir, stem)


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
        description='Decompose the dN/dS contrast into dN and dS by genomic region.')
    p.add_argument('--annotated', required=True,
                   help='Annotated rate table from the region-annotation step')
    p.add_argument('--seed', type=int, default=0,
                   help='RNG seed for jitter and point subsampling. Fixed by '
                        'default so a rerun reproduces the same figure.')
    p.add_argument('--width', choices=['single', 'onehalf', 'double'],
                   default='double',
                   help='Column width of the decomposition strips (default: '
                        'double, i.e. a full-width row of three panels)')
    p.add_argument('--height-mm', type=float, default=58.0)
    p.add_argument('--scatter-width', choices=['single', 'onehalf', 'double'],
                   default='single')
    p.add_argument('--scatter-height-mm', type=float, default=76.0,
                   help='Height of the scatter; the axes box is 49.5 mm whatever '
                        'the height, the rest is margin (default: 76)')
    p.add_argument('--ds-max', type=float, default=1.5,
                   help='Upper dS limit of the scatter; also the point filter, '
                        'so nothing is counted that cannot be seen')
    p.add_argument('--dn-max', type=float, default=0.6)
    p.add_argument('--box-alpha', type=float, default=0.40,
                   help='Opacity of the box fill. The jitter is drawn over the '
                        'box, so a high value hides the points inside the '
                        'interquartile range (default: 0.40)')
    p.add_argument('--cap-pct', type=float, default=99.0,
                   help='Percentile of the pooled values above which points are '
                        'not drawn, for legibility; raised to the higher upper '
                        'whisker end if that is above it. Boxes, whiskers and '
                        'tests always use the complete distribution.')
    p.add_argument('--panel-letters', action='store_true',
                   help='Draw panel letters. Leave off if they are added during '
                        'figure assembly.')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'dnds_decomposition_' + datetime.now().strftime('%Y%m%d'))
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
        handlers=[logging.FileHandler(os.path.join(args.outdir, 'decomposition.log'),
                                     mode='w'),
                  logging.StreamHandler(sys.stdout)])
    for noisy in ('matplotlib.font_manager', 'fontTools', 'PIL'):
        logging.getLogger(noisy).setLevel(logging.ERROR)

    logging.info('=' * 60)
    logging.info('dN/dS decomposition')
    logging.info('=' * 60)

    if not os.path.isfile(args.annotated):
        sys.exit(f'ERROR: input not found: {args.annotated}')
    df = read_rate_table(args.annotated, require=('dN', 'dS', 'dN_dS'))
    for needed in ('region1',):
        if needed not in df.columns:
            sys.exit(f"ERROR: no '{needed}' column in {args.annotated}")
    logging.info(f'Loaded {len(df):,} annotated pairs')

    rng = np.random.default_rng(args.seed)

    schemes = [('region1', 'dnds_decomposition', 'Ascaris annotation', 'A')]
    if 'consensus_region' in df.columns:
        schemes.append(('consensus_region', 'dnds_consensus',
                        'two-species consensus', 'A'))
    else:
        logging.warning('no consensus_region column; consensus figure skipped')

    all_results = {}
    for col, stem, label, letter in schemes:
        logging.info(f'\n--- {label} ({col}) ---')
        res = plot_decomposition(
            df, col, args.outdir, stem, rng, width=args.width,
            height_mm=args.height_mm, cap_pct=args.cap_pct,
            box_alpha=args.box_alpha,
            letter=letter if args.panel_letters else None)
        all_results[col] = res
        for metric, r in res.items():
            logging.info(
                f"  {metric:<6} retained {r['med_ret']:.4f} (n = {r['n_ret']:,})  "
                f"eliminated {r['med_elim']:.4f} (n = {r['n_elim']:,})  "
                f"fold {r['fold']:.3f}  P = {r['p']:.3e}  [{stars(r['p'])}]")
            for key, b in r['box'].items():
                logging.info(
                    f"         {key:<10} box {b['q1']:.4f} / {b['med']:.4f} / "
                    f"{b['q3']:.4f}  whiskers {b['whislo']:.4f} to "
                    f"{b['whishi']:.4f}  ({b['beyond_whiskers']:,} beyond)")
            if r['limit'] > r['cap']:
                logging.info(f"         drawing limit raised from the "
                             f"{args.cap_pct:g}th percentile ({r['cap']:.4g}) to the "
                             f"upper whisker end ({r['limit']:.4g})")
            if r['n_clipped']:
                logging.info(f"         {r['n_clipped']} point(s) above {r['limit']:.4g} "
                             f"not drawn; the box, whiskers and test used all values")

    plot_dn_vs_ds(df, args.outdir, 'dn_vs_ds_scatter', width=args.scatter_width,
                  height_mm=args.scatter_height_mm, ds_max=args.ds_max,
                  dn_max=args.dn_max, letter='B' if args.panel_letters else None)

    # --- what the decomposition means -------------------------------------
    r = all_results['region1']
    dn, ds = r.get('dN'), r.get('dS')
    if dn and ds:
        dn_up = dn['fold'] > 1 and dn['p'] < 0.05
        ds_sig = ds['p'] < 0.05
        if dn_up and not ds_sig:
            verdict = ('Relaxed purifying selection: dN is elevated in eliminated '
                       'genes while dS is indistinguishable, so the dN/dS '
                       'difference cannot be attributed to a local difference in '
                       'mutation rate.')
        elif dn_up and ds_sig and ds['fold'] > 1:
            verdict = ('Elevated mutation rate together with relaxed selection: '
                       'both dN and dS are higher in eliminated genes.')
        elif dn_up and ds_sig and ds['fold'] < 1:
            verdict = ('Relaxed protein-level selection together with stronger '
                       'synonymous constraint: dN is elevated while dS is reduced.')
        elif not dn_up and ds_sig and ds['fold'] < 1:
            verdict = ('Apparent inflation: dN/dS is elevated through a reduced '
                       'denominator rather than an elevated numerator.')
        else:
            verdict = (f"Mixed: dN fold {dn['fold']:.2f} (P = {dn['p']:.2e}), "
                       f"dS fold {ds['fold']:.2f} (P = {ds['p']:.2e}).")
        logging.info(f'\n{verdict}')

    # --- summary -----------------------------------------------------------
    path = os.path.join(args.outdir, 'decomposition_summary.txt')
    with open(path, 'w') as fh:
        fh.write('dN/dS decomposition\n' + '=' * 60 + '\n\n')
        for col, stem, label, _ in schemes:
            fh.write(f'{label} ({col})\n')
            fh.write(f"  {'metric':<8}{'retained':>12}{'eliminated':>12}"
                     f"{'fold':>9}{'P':>12}{'':>6}\n")
            for metric, rr in all_results[col].items():
                fh.write(f"  {metric:<8}{rr['med_ret']:>12.4f}{rr['med_elim']:>12.4f}"
                         f"{rr['fold']:>9.3f}{rr['p']:>12.2e}{stars(rr['p']):>6}\n")
            n = next(iter(all_results[col].values()))
            fh.write(f"  n = {n['n_ret']:,} retained, {n['n_elim']:,} eliminated\n")
            fh.write('  boxes (every gene; whiskers to the most extreme gene within\n'
                     '  1.5 x IQR of the box):\n')
            fh.write(f"    {'metric':<8}{'group':<12}{'Q1':>9}{'median':>9}{'Q3':>9}"
                     f"{'whisker lo':>12}{'whisker hi':>12}{'beyond':>8}\n")
            for metric, rr in all_results[col].items():
                for key in REGIONS:
                    b = rr['box'][key]
                    fh.write(f"    {metric:<8}{key:<12}{b['q1']:>9.4f}{b['med']:>9.4f}"
                             f"{b['q3']:>9.4f}{b['whislo']:>12.4f}{b['whishi']:>12.4f}"
                             f"{b['beyond_whiskers']:>8,}\n")
            fh.write(f'  not drawn (above the {args.cap_pct:g}th percentile of the pooled '
                     f'values, or the\n  higher upper whisker end where that is higher; '
                     f'the boxes, whiskers and\n  tests use them):\n')
            for metric, rr in all_results[col].items():
                raised = (f" (whisker end; the {args.cap_pct:g}th percentile is "
                          f"{rr['cap']:.4g})" if rr['limit'] > rr['cap'] else '')
                fh.write(f"    {metric:<6} above {rr['limit']:.4g}{raised}: "
                         f"{rr['n_clipped_ret']} retained, {rr['n_clipped_elim']} "
                         f"eliminated\n")
            fh.write('\n')
        if dn and ds:
            fh.write(f'Interpretation\n  {verdict}\n\n')
        fh.write('Reporting notes\n')
        fh.write(f'  Points above the {args.cap_pct:g}th percentile (or above the\n')
        fh.write('  higher upper whisker end, where that is higher) are not drawn,\n')
        fh.write('  for legibility. The boxes, the whiskers and every test use the\n')
        fh.write('  complete distribution, so the drawn box is the box tabulated above.\n')
        fh.write('  The scatter is clipped to the axis limits; its legend gives the\n')
        fh.write('  full sample size and the line under it the number of pairs\n')
        fh.write('  falling outside.\n')
        fh.write('  Boxes: quartiles of every gene; centre line, the median.\n')
        fh.write('  Whiskers: Tukey, to the most extreme gene within 1.5 x IQR of\n')
        fh.write('  the box, as in Figure 2B and D.\n')
        fh.write('  Significance markers: **** P < 0.0001, *** P < 0.001,\n')
        fh.write('  ** P < 0.01, * P < 0.05, ns otherwise.\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(),
                   'run_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
                   'host': platform.node(), 'python': sys.version.split()[0],
                   'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                                'matplotlib': matplotlib.__version__},
                   'argv': sys.argv, 'args': vars(args),
                   'input': fingerprint(args.annotated),
                   'rows': int(len(df))}, fh, indent=2)

    logging.info(f'\nSummary: {path}')
    logging.info('Done.')


if __name__ == '__main__':
    main()
