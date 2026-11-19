#!/usr/bin/env python3
"""
tox_pairwise_decomposition.py -- eliminated vs retained dN/dS, dN and dS in the
pairwise comparisons of T. canis with each ascarid (supplementary fig. S11).

The pairwise counterpart of Figure 5C/5D. It reads the annotated tables that
annotate_tox_kaks.py already wrote (kaks_Tox_*_annotated/kaks_annotated.tsv)
-- no annotation is redone -- and puts the five comparisons the May script
drew as five separate files into one figure:

  T. canis x Ascaris        labelled by T. canis (tox_region)
  T. canis x Ascaris        labelled by Ascaris (partner_region)
  T. canis x P. univalens   labelled by T. canis
  T. canis x P. univalens   labelled by P. univalens
  T. canis x B. schroederi  labelled by T. canis (there is no B. schroederi
                            elimination annotation)

Statistics, per comparison and metric
-------------------------------------
  Medians, fold change (eliminated / retained median), Mann-Whitney P, and a
  95% percentile bootstrap interval for the fold (each region resampled
  separately, --n-boot replicates, seeded). BH q within each metric across the
  comparisons (the 'metric' family Figure 5C uses), so dN and dS never share a
  correction family. A comparison's genes are those labelled retained or
  eliminated in its region column, with a finite value of the metric.

Panels (one 185 mm strip, figS11_tox_pairwise)
----------------------------------------------
  A dN/dS, B dN, C dS: fold change with its 95% interval, one row per
  comparison, log axis, dashed line at 1 (no difference). Points are ink:
  each is an eliminated/retained ratio, not a region. Asterisks from BH q in
  a column at the right of each panel; who labelled the genes, and n, in
  columns right of C.

Absolute dN/dS in these comparisons is not comparable with the ascarid pairs
(dS near saturation, time dependence); only the within-comparison contrast is
read, which is what a fold shows.

Outputs (in --outdir)
---------------------
  figS11_tox_pairwise.{svg,pdf,png}
  tox_pairwise_stats.tsv
  tox_pairwise_summary.txt
  tox_pairwise_decomposition.log
  provenance.json

Usage
-----
  python tox_pairwise_decomposition.py \\
    --set 'T. canis × Ascaris' kaks_Tox_As_annotated/kaks_annotated.tsv tox_region 'T. canis' \\
    --set 'T. canis × Ascaris' kaks_Tox_As_annotated/kaks_annotated.tsv partner_region Ascaris \\
    ... --outdir results/tox_pairwise_<stamp>

Dependencies: numpy, pandas, scipy, matplotlib; paper_colors, region_lib and
rate_tables from lib/.
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
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from matplotlib.transforms import ScaledTranslation
from scipy.stats import mannwhitneyu
from scipy import __version__ as scipy_version

# A script that cannot import the shared style must stop, not fall back to
# matplotlib defaults.
from paper_colors import (NEUTRAL, FONT, LINEWIDTH, apply_style, audit, figure,
                          save, mm)
from region_lib import (benjamini_hochberg, fingerprint, git_commit,
                        setup_logging, stars)
from rate_tables import read_rate_table

_HERE = os.path.dirname(os.path.abspath(__file__))

METRICS = ('dN_dS', 'dN', 'dS')
METRIC_LABEL = {'dN_dS': 'dN/dS', 'dN': 'dN', 'dS': 'dS'}
REGIONS = ('retained', 'eliminated')


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def p_text(p):
    """Plain-text P for the summary."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p == 0:
        return '< 1e-300'
    return f'{p:.2g}' if p < 0.001 else f'{p:.3f}'


def fold_ci(ret, elim, n_boot, rng, level=0.95):
    """Percentile bootstrap interval for median(elim) / median(ret), each
    region resampled separately. Replicates with a zero retained median are
    dropped and counted."""
    r = rng.choice(ret, size=(n_boot, len(ret)), replace=True)
    e = rng.choice(elim, size=(n_boot, len(elim)), replace=True)
    mr, me = np.median(r, axis=1), np.median(e, axis=1)
    ok = mr > 0
    f = me[ok] / mr[ok]
    if not len(f):
        return np.nan, np.nan, int((~ok).sum())
    a = (1 - level) / 2
    return (float(np.quantile(f, a)), float(np.quantile(f, 1 - a)),
            int((~ok).sum()))


def compare(df, region_col, metric, n_boot, rng):
    """One comparison x metric: medians, fold, Mann-Whitney, bootstrap CI."""
    v = pd.to_numeric(df[metric], errors='coerce')
    ok = np.isfinite(v) & df[region_col].isin(REGIONS)
    r = v[ok & (df[region_col] == 'retained')].to_numpy(float)
    e = v[ok & (df[region_col] == 'eliminated')].to_numpy(float)
    out = {'metric': metric, 'n_retained': len(r), 'n_eliminated': len(e),
           'median_retained': np.median(r) if len(r) else np.nan,
           'median_eliminated': np.median(e) if len(e) else np.nan,
           'fold': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan, 'p': np.nan,
           'n_boot_dropped': 0}
    if len(r) >= 5 and len(e) >= 5:
        if out['median_retained'] > 0:
            out['fold'] = out['median_eliminated'] / out['median_retained']
        out['p'] = float(mannwhitneyu(e, r, alternative='two-sided')[1])
        out['ci_lo'], out['ci_hi'], out['n_boot_dropped'] = fold_ci(
            r, e, n_boot, rng)
    return out


# ---------------------------------------------------------------------------
# Figure helpers (shared pattern with the other converted scripts)
# ---------------------------------------------------------------------------

def blank_figure(width, height_mm):
    """A journal-width canvas with no axes and no layout engine; panels are
    placed explicitly, because constrained_layout shrinks a panel to nothing
    to make room for text it cannot fit and only warns."""
    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax)
    return fig


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
        leg = ax.get_legend()
        if leg is not None:
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


def place_row(fig, axes, shares, gaps_mm, pad_mm=0.8, iterations=4):
    """Lay out panels side by side from their measured text; the row fills
    the canvas whatever font the node renders with."""
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
        avail = W - 2 * pad - sum(lefts) - sum(rights) - sum(gaps)
        widths = avail * shares
        x = pad
        y0, hh = pad + bottom, H - 2 * pad - bottom - top
        for i, ax in enumerate(axes):
            x += lefts[i]
            ax.set_position([x / W, y0 / H, widths[i] / W, hh / H])
            x += widths[i] + rights[i] + (gaps[i] if i < len(gaps) else 0.0)


def finish(fig, outdir, stem):
    """Audit, check text, save in every format, close."""
    problems = list(audit(fig, verbose=False)) + text_problems(fig)
    for prob in problems:
        logging.warning('  audit %s: %s', stem, prob)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths, problems


def column_text(ax, y, s, dx_pt, **kw):
    """Text in a column right of the axes, `dx_pt` points from the spine."""
    tr = ax.get_yaxis_transform() + ScaledTranslation(dx_pt / 72, 0,
                                                      ax.figure.dpi_scale_trans)
    opts = dict(ha='left', va='center_baseline', fontsize=FONT['annotation'],
                color=NEUTRAL['ink'], clip_on=False)
    opts.update(kw)
    return ax.text(1.0, y, s, transform=tr, **opts)


def nice_log_ticks(lo, hi):
    """Doublings (0.25, 0.5, 1, 2, 4, 8) inside [lo, hi]: evenly spaced on a
    log axis, so the labels cannot crowd however narrow the panel is."""
    ticks = [t for t in (0.125, 0.25, 0.5, 1, 2, 4, 8, 16) if lo <= t <= hi]
    return ticks if len(ticks) >= 2 else [1, 2]


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def plot_s11(stats, sets, outdir, stem, args):
    """Three forest panels (dN/dS, dN, dS) sharing the comparison rows."""
    fig = blank_figure('double', args.height_mm)
    axes = [fig.add_axes([0.25 + 0.22 * k, 0.2, 0.18, 0.62]) for k in range(3)]
    ys = -np.arange(len(sets), dtype=float)
    finite = stats[np.isfinite(stats['ci_lo']) & np.isfinite(stats['ci_hi'])]
    lo = min(0.7, float(finite['ci_lo'].min()) * 0.9) if len(finite) else 0.5
    hi = max(2.1, float(finite['ci_hi'].max()) * 1.1) if len(finite) else 2.1
    for k, (ax, m) in enumerate(zip(axes, METRICS)):
        st = stats[stats['metric'] == m].set_index('set_index')
        ax.axvline(1.0, color=NEUTRAL['rule'], lw=LINEWIDTH['rule'],
                   ls=(0, (4, 2)), zorder=1)
        for i, y in enumerate(ys):
            if i not in st.index:
                continue
            r = st.loc[i]
            if np.isfinite(r['ci_lo']):
                ax.plot([r['ci_lo'], r['ci_hi']], [y, y], color=NEUTRAL['ink'],
                        lw=LINEWIDTH['thin'], solid_capstyle='butt', zorder=2)
            if np.isfinite(r['fold']):
                ax.plot([r['fold']], [y], ls='none', marker='o', ms=3.6,
                        mfc=NEUTRAL['ink'], mec=NEUTRAL['ink'], zorder=3)
            # Asterisks just right of the interval (or the point).
            x_end = r['ci_hi'] if np.isfinite(r['ci_hi']) else r['fold']
            if np.isfinite(x_end):
                ax.text(x_end, y, '  ' + stars(r['q']), ha='left',
                        va='center_baseline', fontsize=FONT['annotation'],
                        color=NEUTRAL['ink'])
        ax.set_xscale('log')
        ax.set_xlim(lo, hi * 1.45)         # room on the right for the asterisks
        ticks = nice_log_ticks(lo, hi)
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(ys[-1] - 0.6, ys[0] + 0.6)
        ax.yaxis.set_major_locator(FixedLocator(ys))
        ax.set_yticklabels([s['label'] for s in sets] if k == 0 else [])
        ax.tick_params(axis='y', length=0)
        ax.set_xlabel(f'{METRIC_LABEL[m]} fold change\n(eliminated / retained)')
    # who labelled the genes, and n, right of the dS panel
    last = axes[-1]
    st = stats[stats['metric'] == 'dN_dS'].set_index('set_index')
    column_text(last, ys[0] + 0.85, 'Labelled by', 6.0,
                color=NEUTRAL['ink_light'])
    column_text(last, ys[0] + 0.85, 'n retained / eliminated', 62.0,
                color=NEUTRAL['ink_light'])
    for i, (y, s) in enumerate(zip(ys, sets)):
        column_text(last, y, s['labelled_by'], 6.0)
        if i in st.index:
            column_text(last, y, f"{int(st.loc[i, 'n_retained']):,} / "
                                 f"{int(st.loc[i, 'n_eliminated']):,}", 62.0)
    place_row(fig, axes, shares=[1, 1, 1], gaps_mm=[5.0, 5.0])
    paths, problems = finish(fig, outdir, stem)
    return paths, problems, {'x_range': (float(lo), float(hi * 1.45))}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Eliminated vs retained dN/dS, dN and dS in the T. canis '
                    'pairwise comparisons (supplementary fig. S11).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--set', nargs=4, action='append', required=True,
                   metavar=('LABEL', 'TABLE', 'REGION_COL', 'LABELLED_BY'),
                   help='One comparison: its label, annotated table, the region '
                        'column to use and who that column\'s labels come from. '
                        'Repeat, top row first.')
    p.add_argument('--n-boot', type=int, default=2000)
    p.add_argument('--height-mm', type=float, default=52.0)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.n_boot < 200:
        p.error('--n-boot must be >= 200')
    if args.outdir is None:
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.outdir = os.path.join('results', f'tox_pairwise_{stamp}')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'tox_pairwise_decomposition')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('tox_pairwise_decomposition.py -- supplementary fig. S11')
    apply_style()

    sets, rows, cache = [], [], {}
    for i, (label, table, region_col, who) in enumerate(args.set):
        if not os.path.isfile(table):
            sys.exit(f'ERROR: not found: {table}  ({label})')
        if table not in cache:
            cache[table] = read_rate_table(table, require=METRICS)
        df = cache[table]
        if region_col not in df.columns:
            sys.exit(f'ERROR: column {region_col!r} not in {table}. Columns: '
                     f'{list(df.columns)}')
        counts = df[region_col].value_counts(dropna=False).to_dict()
        logging.info('  %-26s %-15s labelled by %-13s rows %s; %s', label,
                     region_col, who, f'{len(df):,}', counts)
        sets.append({'label': label, 'table': table, 'region_col': region_col,
                     'labelled_by': who, 'n_rows': len(df),
                     'region_counts': {str(k): int(v) for k, v in counts.items()}})
        rng = np.random.default_rng(args.seed + i)   # one stream per comparison
        for m in METRICS:
            r = compare(df, region_col, m, args.n_boot, rng)
            r.update(set_index=i, label=label, labelled_by=who,
                     region_col=region_col, table=table)
            rows.append(r)
    stats = pd.DataFrame(rows)
    stats['q'] = np.nan
    for m, idx in stats.groupby('metric').groups.items():
        stats.loc[idx, 'q'] = benjamini_hochberg(
            stats.loc[idx, 'p'].to_numpy(float))
    cols = ['set_index', 'label', 'labelled_by', 'region_col', 'metric',
            'n_retained', 'n_eliminated', 'median_retained', 'median_eliminated',
            'fold', 'ci_lo', 'ci_hi', 'p', 'q', 'n_boot_dropped', 'table']
    stats[cols].to_csv(os.path.join(args.outdir, 'tox_pairwise_stats.tsv'),
                       sep='\t', index=False, float_format='%.6g')
    for r in stats.itertuples():
        logging.info('  %-26s %-13s %-6s fold %.3f [%.3f, %.3f]  P %s  q %s',
                     r.label, r.labelled_by, METRIC_LABEL[r.metric], r.fold,
                     r.ci_lo, r.ci_hi, p_text(r.p), p_text(r.q))

    logging.info('S11: forest strip')
    paths, problems, fig_info = plot_s11(stats, sets, args.outdir,
                                         'figS11_tox_pairwise', args)

    # --- summary ----------------------------------------------------------------
    L = ['T. canis pairwise comparisons, eliminated vs retained '
         '(supplementary fig. S11)', '=' * 72, '']
    for s in sets:
        L.append(f"  {s['label']:<26} labelled by {s['labelled_by']:<13} "
                 f"({s['region_col']}); {s['n_rows']:,} rows in {s['table']}")
    L += ['', f"{'comparison':<26} {'by':<13} {'metric':<6} {'n ret':>6} "
          f"{'n elim':>6} {'med ret':>8} {'med elim':>8} {'fold':>6} "
          f"{'95% CI':>13} {'P':>9} {'q':>9}"]
    for r in stats.itertuples():
        L.append(f"{r.label:<26} {r.labelled_by:<13} {METRIC_LABEL[r.metric]:<6} "
                 f"{r.n_retained:>6,} {r.n_eliminated:>6,} "
                 f"{r.median_retained:>8.4f} {r.median_eliminated:>8.4f} "
                 f"{r.fold:>6.2f} {r.ci_lo:>6.2f}-{r.ci_hi:<6.2f} "
                 f"{p_text(r.p):>9} {p_text(r.q):>9}")
    L += ['', f'Fold = median eliminated / median retained; 95% percentile '
          f'bootstrap, {args.n_boot} replicates, each region resampled '
          f'separately (seed {args.seed} + comparison index).',
          'q = BH within each metric across the comparisons (dN and dS never '
          'share a family).', '']
    up = {m: stats[(stats['metric'] == m)] for m in METRICS}
    L += ['What the text can claim', '-' * 72]
    for m in ('dN_dS', 'dN'):
        t = up[m]
        allup = bool(((t['fold'] > 1) & (t['q'] < 0.05)).all())
        L.append(f"  {METRIC_LABEL[m]}: " + (
            f'higher in eliminated genes in every comparison (folds '
            f"{t['fold'].min():.2f}-{t['fold'].max():.2f}; largest q "
            f"{p_text(t['q'].max())})." if allup else
            'NOT higher in every comparison: ' + '; '.join(
                f"{r.label} by {r.labelled_by} fold {r.fold:.2f}, q {p_text(r.q)}"
                for r in t.itertuples() if not (r.fold > 1 and r.q < 0.05))))
    t = up['dS']
    sig = t[t['q'] < 0.05]
    L.append(f"  dS: folds {t['fold'].min():.2f}-{t['fold'].max():.2f}; "
             + (f'{len(sig)} of {len(t)} significant after BH ('
                + '; '.join(f'{r.label} by {r.labelled_by} {r.fold:.2f}'
                            for r in sig.itertuples()) + ')'
                if len(sig) else 'none significant after BH') + '.')
    L += ['', 'Drawn (figS11_tox_pairwise)',
          f"  x-axis {fig_info['x_range'][0]:.2f}-{fig_info['x_range'][1]:.2f} "
          '(log); every interval is inside it.',
          '  AUDIT: ' + '; '.join(problems) if problems else '  audit: clean', '']
    with open(os.path.join(args.outdir, 'tox_pairwise_summary.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(),
                   'python': platform.python_version(),
                   'versions': {'numpy': np.__version__,
                                'pandas': pd.__version__,
                                'scipy': scipy_version,
                                'matplotlib': matplotlib.__version__},
                   'argv': sys.argv,
                   'inputs': {t: fingerprint(t) for t in cache},
                   'sets': sets,
                   'figure': [os.path.basename(p) for p in paths],
                   'audit': problems},
                  fh, indent=2, default=str)
    logging.info('Summary: %s', os.path.join(args.outdir,
                                             'tox_pairwise_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
