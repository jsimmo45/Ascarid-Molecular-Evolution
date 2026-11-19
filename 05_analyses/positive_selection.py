#!/usr/bin/env python3
"""
positive_selection.py - genes with omega > 1, and whether they concentrate in
eliminated regions.

Supplementary fig. S3A and supplementary table S3 (S8A and T3 before the
2026-10-12 renumbering).

Writes:

    S8A_omega_above_one.{svg,pdf,png}   fig. S3A. The file keeps its old stem so
                                        the assembly manifest and
                                        figure_locations stay valid; the
                                        assembly copies it to FigS03_A.
    omega_above_one_genes.tsv           table S3A, sorted by omega
    positive_selection_stats.tsv        counts, proportions, tests, sensitivity
                                        (table S3B)
    positive_selection_summary.txt      the same, formatted for reading
    ds_diagnostic.*                     dS and alignment length of the omega > 1
                                        set against the background -- a check,
                                        not a manuscript panel
    provenance.json                     inputs, arguments and versions


"omega > 1" is not "positively selected"
----------------------------------------
The earlier version called these genes positively selected. A single pairwise
omega above 1, with no test attached, is not evidence of positive selection: it is
a point estimate of a ratio whose variance explodes as dS goes to zero, and with
~9,400 pairs a tail of estimation noise above 1 is expected even under strict
purifying selection everywhere. Nothing here fits a model in which omega > 1 is
tested against omega = 1, because pairwise yn00 does not provide one.

So the category is named for what it is - omega > 1 - throughout the outputs, and
two things are reported that let a reader judge how much of the tail is real:

  - the dS and alignment-length distributions of the omega > 1 set against the
    background. If those genes sit at systematically low dS or on short
    alignments, the tail is dominated by estimation noise.
  - the count and the enrichment test recomputed across a range of minimum-dS
    thresholds. A result that survives `--min-ds 0.05` is worth more than one
    that evaporates.

**The manuscript should say "genes with dN/dS > 1" and state the caveat**, or run
a codon-model test on that gene set and report that instead. Calling 46 genes
positively selected on this evidence is the sentence a referee will pick out.


What changed besides the naming
-------------------------------
- **The stacked bar chart is gone.** It plotted "omega <= 1" (about 99.5% of
  genes) stacked under "omega > 1", so the quantity of interest was a sliver at
  the top of a full-height bar. The panel is now the two proportions themselves,
  with Wilson intervals.
- **The Ka/Ks histogram is gone.** S3B (cbr_changepoint_relax.py) is the
  log-omega density for the same two gene sets; the two panels of one supplementary figure were showing the same
  distribution twice.
- **The odds ratio is reported eliminated-over-retained**, so OR > 1 means
  eliminated genes are enriched. The old orientation returned OR < 1 for
  enrichment and needed a paragraph of prose to interpret.
- **The odds ratio has a confidence interval.**
- `fmt_p` replaces `p = 1.00e+00`.


Changes 2026-10-12
------------------
- **The odds ratio is the conditional maximum-likelihood estimate, with the
  exact 95% interval that inverts Fisher's test** (scipy's
  `contingency.odds_ratio(kind='conditional')`, R's `fisher.test`). It was the
  Haldane-Anscombe estimate (0.5 added to every cell) with a Woolf interval,
  which disagreed with Methods and with every other odds ratio in the paper
  (`microsynteny.py`, `functional_enrichment.py`, the same code). The 0.5
  correction is for zero cells, and these tables have none; it pulled the
  estimate toward 1: 0.93 (0.18 to 4.73) against the conditional 0.63 (0.015
  to 3.69), while the Results sentence prints the counts, 1 of 322 against 45
  of 9,085, from which a reader computes about 0.63. Fisher's P is unchanged.
- `positive_selection_stats.tsv` carries an `or_method` column: 'scipy', or
  'fallback' (the same estimator, from the noncentral hypergeometric
  distribution, when scipy is older than 1.10).
- The S8A annotation prints a bound below 0.1 to two significant digits (an
  exact lower bound of 0.015 printed as 0.01 with two decimals).
- `make_supplementary_tables.py` checks that S3B holds conditional estimates.
- Figure text, from the review of fig. S3: no species title (--species goes
  into the summary; the figure specification puts it in the caption), and the
  y label and the diagnostic's tick labels are plain text with dN/dS upright
  and a Unicode less-than-or-equal sign, not mathtext. The statistics and
  tables are unchanged.

Usage
-----
    python positive_selection.py \\
      --annotated results/annotate_regions_20260817/kaks_annotated.tsv \\
      --species Ascaris \\
      --outdir results/positive_selection_20260817

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
from matplotlib.colors import to_rgba
from scipy.stats import fisher_exact, mannwhitneyu, norm
from scipy import __version__ as scipy_version

from paper_colors import (REGION, REGION_ORDER, NEUTRAL, apply_style, audit,
                          save, figure, panel_letter, mm)
from region_lib import (fmt_p, stars, benjamini_hochberg, fingerprint,
                        git_commit, setup_logging)

_HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def wilson(k, n, alpha=0.05):
    """Wilson score interval for a binomial proportion.

    Used rather than the normal approximation because the proportions here are
    well under 1% on a few hundred genes, where the Wald interval routinely
    includes negative values.
    """
    if n == 0:
        return (np.nan, np.nan)
    z = norm.ppf(1 - alpha / 2)
    p = k / n
    d = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def _conditional_or_fallback(a, b, c, d, level):
    """Conditional MLE odds ratio and exact CI for scipy < 1.10 (no
    scipy.stats.contingency.odds_ratio): `a` follows Fisher's noncentral
    hypergeometric distribution given the margins; the limits are where the
    observed count sits in either (1 - level)/2 tail, as in R's fisher.test.
    The same code as microsynteny.py and functional_enrichment.py."""
    from scipy.optimize import brentq
    from scipy.stats import nchypergeom_fisher
    M, n_col, N = a + b + c + d, a + c, a + b
    lo_sup, hi_sup = max(0, N - (M - n_col)), min(N, n_col)
    alpha = (1 - level) / 2

    def dist(log_or):
        return nchypergeom_fisher(M, n_col, N, np.exp(log_or))

    bound = 20.0
    mle = np.nan
    if a == lo_sup:
        mle, lo = 0.0, 0.0
    else:
        lo = np.exp(brentq(lambda x: dist(x).sf(a - 1) - alpha, -bound, bound))
    if a == hi_sup:
        mle, hi = np.inf, np.inf
    else:
        hi = np.exp(brentq(lambda x: dist(x).cdf(a) - alpha, -bound, bound))
    if lo_sup < a < hi_sup:
        mle = np.exp(brentq(lambda x: dist(x).mean() - a, -bound, bound))
    return mle, lo, hi


def odds_ratio_ci(a, b, c, d, alpha=0.05):
    """OR = odds(omega > 1 | eliminated) / odds(omega > 1 | retained), with CI.

    Oriented eliminated-over-retained so that OR > 1 means eliminated genes are
    enriched, matching the direction of every other comparison in the paper.

    The conditional maximum-likelihood estimate with the exact interval that
    inverts Fisher's test, the pairing that goes with the Fisher P reported
    beside it and the estimator of every other odds ratio in the paper. It is
    finite unless a cell is zero (then 0 or infinity, with a finite exact
    bound on the other side), and the counts here have no zero cell.

    a, b: eliminated with and without omega > 1
    c, d: retained with and without omega > 1

    Returns (odds ratio, lower, upper, method), method 'scipy' or 'fallback'.
    """
    level = 1 - alpha
    try:
        from scipy.stats.contingency import odds_ratio
    except ImportError:
        mle, lo, hi = _conditional_or_fallback(a, b, c, d, level)
        return float(mle), float(lo), float(hi), 'fallback'
    res = odds_ratio(np.array([[a, b], [c, d]]), kind='conditional')
    ci = res.confidence_interval(confidence_level=level)
    return float(res.statistic), float(ci.low), float(ci.high), 'scipy'


def fmt_or(x):
    """An odds ratio or bound for the panel: two significant digits below 0.1
    (an exact lower bound of 0.015 prints as 0.01 with two decimals)."""
    if x is None or np.isnan(x):
        return 'NA'
    if np.isinf(x):
        return '\u221e'
    if x == 0:
        return '0'
    return f'{x:.2g}' if x < 0.1 else f'{x:.2f}'


def enrichment(df, region_col, omega_col, threshold):
    """Counts, proportions and the Fisher test for one omega threshold."""
    hi = df[omega_col] > threshold
    n_el = int((df[region_col] == 'eliminated').sum())
    n_re = int((df[region_col] == 'retained').sum())
    k_el = int((hi & (df[region_col] == 'eliminated')).sum())
    k_re = int((hi & (df[region_col] == 'retained')).sum())
    if min(n_el, n_re) == 0:
        return None
    # fisher_exact on [[elim_hi, elim_lo], [ret_hi, ret_lo]] gives the
    # eliminated-over-retained orientation directly.
    _, p = fisher_exact([[k_el, n_el - k_el], [k_re, n_re - k_re]],
                        alternative='two-sided')
    or_, lo, hi_ci, how = odds_ratio_ci(k_el, n_el - k_el, k_re, n_re - k_re)
    return {'threshold': threshold,
            'n_eliminated': n_el, 'k_eliminated': k_el,
            'n_retained': n_re, 'k_retained': k_re,
            'prop_eliminated': k_el / n_el, 'prop_retained': k_re / n_re,
            'ci_eliminated': wilson(k_el, n_el),
            'ci_retained': wilson(k_re, n_re),
            'odds_ratio': or_, 'or_lo95': lo, 'or_hi95': hi_ci,
            'or_method': how, 'p': float(p)}


def ds_sensitivity(df, region_col, omega_col, ds_col, threshold, cuts):
    """Recompute the enrichment after dropping genes below each dS cut.

    omega is a ratio, and its sampling variance grows without bound as dS goes to
    zero, so the omega > 1 tail is where estimation noise concentrates. If the
    enrichment is an artefact of that noise it will not survive a dS floor.
    """
    rows = []
    for cut in cuts:
        sub = df[df[ds_col] >= cut]
        e = enrichment(sub, region_col, omega_col, threshold)
        if e is None:
            continue
        e['min_ds'] = cut
        e['n_total'] = int(len(sub))
        rows.append(e)
    return rows


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _axes(width, height_mm, left=16.0, right=3.0, bottom=14.0, top=4.0):
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


def plot_proportions(res, outdir, stem, threshold, width='single',
                     height_mm=58.0, bar_alpha=0.40, panel_letters=False):
    """Fig. S3A: the proportion of genes with omega > threshold, by region.

    Two bars and their Wilson intervals. The earlier version stacked
    "omega <= 1" beneath "omega > 1", which meant the quantity the panel exists
    to show was a two-pixel sliver on top of a full-height bar and the reader had
    to take the annotation's word for it.
    """
    fig, ax = _axes(width, height_mm)
    order = REGION_ORDER                       # retained left, eliminated right
    vals = [res[f'prop_{r}'] * 100 for r in order]
    cis = [tuple(v * 100 for v in res[f'ci_{r}']) for r in order]
    ks = [res[f'k_{r}'] for r in order]
    ns = [res[f'n_{r}'] for r in order]

    pos = np.arange(len(order))
    for i, r in enumerate(order):
        ax.bar(pos[i], vals[i], width=0.6, zorder=2,
               color=to_rgba(REGION[r], bar_alpha), edgecolor='black',
               linewidth=0.8)
        ax.plot([pos[i], pos[i]], [cis[i][0], cis[i][1]], color='black',
                lw=0.9, zorder=3, solid_capstyle='butt')
        for y in cis[i]:
            ax.plot([pos[i] - 0.10, pos[i] + 0.10], [y, y], color='black',
                    lw=0.9, zorder=3)
        ax.annotate(f'{ks[i]}/{ns[i]:,}', xy=(pos[i], cis[i][1]),
                    xytext=(0, 3), textcoords='offset points', ha='center',
                    va='bottom', fontsize=plt.rcParams['xtick.labelsize'],
                    zorder=4)

    ax.set_xticks(pos)
    ax.set_xticklabels([f'{r.capitalize()}\n(n = {n:,})'
                        for r, n in zip(order, ns)])
    ax.set_xlim(-0.65, len(order) - 0.35)
    ax.set_ylim(0, max(c[1] for c in cis) * 1.45)
    ax.set_ylabel(f'Genes with dN/dS > {threshold:g} (%)')

    ax.annotate(f"OR = {fmt_or(res['odds_ratio'])} "
                f"[{fmt_or(res['or_lo95'])}, {fmt_or(res['or_hi95'])}]\n"
                f"{fmt_p(res['p'])}  {stars(res['p'])}",
                xy=(0.03, 0.97), xycoords='axes fraction', ha='left', va='top',
                fontsize=plt.rcParams['legend.fontsize'])
    if panel_letters:
        panel_letter(ax, 'A')
    return _finish(fig, outdir, stem)


def plot_ds_diagnostic(df, hi_mask, ds_col, len_col, outdir, stem,
                       threshold, width='double', height_mm=58.0):
    """Diagnostic, not a manuscript panel.

    Where the omega > 1 genes sit in dS and in alignment length relative to
    everything else. Both are the axes along which a pairwise omega estimate
    becomes unreliable, so this is the panel that says how seriously to take the
    tail.
    """
    cols = [c for c in (ds_col, len_col) if c and c in df.columns]
    if not cols:
        logging.info('  no dS or alignment-length column; skipping diagnostic')
        return None
    fig, ax0 = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax0)
    w_in, h_in = fig.get_size_inches()
    left, right, bottom, top, gap = mm(16.0), mm(3.0), mm(14.0), mm(5.0), mm(16.0)
    aw = (w_in - left - right - gap * (len(cols) - 1)) / len(cols)

    for i, col in enumerate(cols):
        ax = fig.add_axes([(left + i * (aw + gap)) / w_in, bottom / h_in,
                           aw / w_in, (h_in - bottom - top) / h_in])
        bg = df.loc[~hi_mask, col].dropna().to_numpy(float)
        fg = df.loc[hi_mask, col].dropna().to_numpy(float)
        if bg.size < 5 or fg.size < 5:
            continue
        data = [bg, fg]
        bp = ax.boxplot(data, positions=[0, 1], widths=0.55, patch_artist=True,
                        showfliers=False, whis=(5, 95), zorder=2)
        for patch, colour in zip(bp['boxes'],
                                 [NEUTRAL['ink_light'], REGION['eliminated']]):
            patch.set_facecolor(to_rgba(colour, 0.40))
            patch.set_edgecolor('black')
            patch.set_linewidth(0.8)
            patch.set_alpha(None)
        for part in ('whiskers', 'caps'):
            for a in bp[part]:
                a.set_color('black')
                a.set_linewidth(0.8)
        for med in bp['medians']:
            med.set_color('black')
            med.set_linewidth(1.4)
        u, p = mannwhitneyu(bg, fg, alternative='two-sided')
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'dN/dS \u2264 {threshold:g}\n(n = {bg.size:,})',
                            f'> {threshold:g}\n(n = {fg.size:,})'])
        ax.set_ylabel(col)
        ax.set_xlim(-0.6, 1.6)
        ax.annotate(f'{fmt_p(p)}', xy=(0.5, 0.97), xycoords='axes fraction',
                    ha='center', va='top',
                    fontsize=plt.rcParams['xtick.labelsize'])
    return _finish(fig, outdir, stem)


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Genes with dN/dS > 1 and their distribution between '
                    'regions (supplementary fig. S3A, written to '
                    'S8A_omega_above_one.*; table S3).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', required=True,
                   help='kaks_annotated.tsv from annotate_regions.py')
    p.add_argument('--threshold', type=float, default=1.0,
                   help='dN/dS above which a gene is counted (default: 1.0). '
                        'This is a threshold, not a test.')
    p.add_argument('--omega-col', default='Ka_Ks')
    p.add_argument('--ds-col', default='Ks')
    p.add_argument('--dn-col', default='Ka')
    p.add_argument('--len-col', default='align_len',
                   help='Alignment length column, for the diagnostic. Ignored if '
                        'absent.')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--min-ds', type=float, default=0.0,
                   help='Drop genes below this dS before the primary analysis. '
                        'Default 0: the primary result uses everything and the '
                        'sensitivity table shows what a floor would do.')
    p.add_argument('--ds-cuts', type=float, nargs='*',
                   default=[0.0, 0.01, 0.02, 0.05, 0.10],
                   help='dS floors for the sensitivity table')
    p.add_argument('--omega-max', type=float, default=None,
                   help='Optional ceiling on dN/dS. Off by default: capping the '
                        'tail is exactly what this analysis is about.')
    p.add_argument('--species', default='',
                   help='Recorded in the summary. No longer drawn as a title: '
                        'the species goes in the caption.')
    p.add_argument('--bar-alpha', type=float, default=0.40)
    p.add_argument('--width', choices=['single', 'onehalf', 'double'],
                   default='single')
    p.add_argument('--panel-letters', action='store_true')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'positive_selection_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'positive_selection')
    logging.info('positive_selection.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()

    df = pd.read_csv(args.annotated, sep='\t')
    logging.info('Loaded %s rows from %s', f'{len(df):,}',
                 os.path.basename(args.annotated))
    for col in (args.omega_col, args.region_col):
        if col not in df.columns:
            sys.exit(f'ERROR: {col!r} is not a column of {args.annotated}\n'
                     f'  columns: {list(df.columns)}')

    d = df[df[args.region_col].isin(['retained', 'eliminated'])].copy()
    d[args.omega_col] = pd.to_numeric(d[args.omega_col], errors='coerce')
    d = d[np.isfinite(d[args.omega_col]) & (d[args.omega_col] >= 0)]
    if args.omega_max is not None:
        n_before = len(d)
        d = d[d[args.omega_col] <= args.omega_max]
        logging.info('  --omega-max %g dropped %s genes', args.omega_max,
                     f'{n_before - len(d):,}')
    has_ds = args.ds_col in d.columns
    if args.min_ds > 0 and has_ds:
        n_before = len(d)
        d = d[pd.to_numeric(d[args.ds_col], errors='coerce') >= args.min_ds]
        logging.info('  --min-ds %g dropped %s genes', args.min_ds,
                     f'{n_before - len(d):,}')
    logging.info('  %s genes with a region call', f'{len(d):,}')

    res = enrichment(d, args.region_col, args.omega_col, args.threshold)
    if res is None:
        sys.exit('ERROR: one region has no genes.')
    k_tot = res['k_retained'] + res['k_eliminated']
    logging.info('  dN/dS > %g: %s of %s genes (%.2f%%)', args.threshold,
                 f'{k_tot:,}', f'{len(d):,}', 100.0 * k_tot / len(d))
    logging.info('  retained %d/%s (%.2f%%), eliminated %d/%s (%.2f%%)',
                 res['k_retained'], f"{res['n_retained']:,}",
                 100 * res['prop_retained'], res['k_eliminated'],
                 f"{res['n_eliminated']:,}", 100 * res['prop_eliminated'])
    logging.info('  OR (eliminated / retained) = %.3g [%.3g, %.3g], P = %.3g '
                 '(conditional MLE, exact 95%% CI; %s)', res['odds_ratio'],
                 res['or_lo95'], res['or_hi95'], res['p'], res['or_method'])

    # --- sensitivity -------------------------------------------------------
    sens = []
    if has_ds and args.ds_cuts:
        d2 = d.copy()
        d2[args.ds_col] = pd.to_numeric(d2[args.ds_col], errors='coerce')
        d2 = d2.dropna(subset=[args.ds_col])
        sens = ds_sensitivity(d2, args.region_col, args.omega_col, args.ds_col,
                              args.threshold, args.ds_cuts)
        if sens:
            q = benjamini_hochberg([s['p'] for s in sens])
            for s, qi in zip(sens, q):
                s['q'] = float(qi)
            logging.info('  dS sensitivity:')
            for s in sens:
                logging.info('    min dS %.2f: %d + %d genes above threshold, '
                             'OR %.3g [%.3g, %.3g], P = %.3g', s['min_ds'],
                             s['k_retained'], s['k_eliminated'],
                             s['odds_ratio'], s['or_lo95'], s['or_hi95'], s['p'])

    # --- outputs -----------------------------------------------------------
    hi = d[args.omega_col] > args.threshold
    keep = [c for c in (['gene1', 'gene2', args.dn_col, args.ds_col,
                         args.omega_col, args.region_col, 'consensus_region',
                         'elim_type1', 'chrom1', 'start1', 'end1',
                         'signed_cbr_dist1', 'cbr_proximal1', args.len_col])
            if c in d.columns]
    table = d.loc[hi, keep].sort_values(args.omega_col, ascending=False)
    table.to_csv(os.path.join(args.outdir, 'omega_above_one_genes.tsv'),
                 sep='\t', index=False, float_format='%.10g')
    logging.info('  wrote omega_above_one_genes.tsv (%s genes)', f'{len(table):,}')

    stats_rows = [{'analysis': 'primary', 'min_ds': args.min_ds, **{
        k: v for k, v in res.items() if not k.startswith('ci_')}}]
    for s in sens:
        stats_rows.append({'analysis': 'ds_sensitivity',
                           **{k: v for k, v in s.items()
                              if not k.startswith('ci_')}})
    pd.DataFrame(stats_rows).to_csv(
        os.path.join(args.outdir, 'positive_selection_stats.tsv'),
        sep='\t', index=False, float_format='%.6g')

    logging.info('Supplementary fig. S3A (S8A_omega_above_one)')
    plot_proportions(res, args.outdir, 'S8A_omega_above_one', args.threshold,
                     width=args.width,
                     bar_alpha=args.bar_alpha,
                     panel_letters=args.panel_letters)
    logging.info('Diagnostic')
    plot_ds_diagnostic(d, hi, args.ds_col if has_ds else None, args.len_col,
                       args.outdir, 'ds_diagnostic', args.threshold)

    # --- summary -----------------------------------------------------------
    L = [f'Genes with dN/dS > {args.threshold:g}', '',
         f'Table            {os.path.basename(args.annotated)}']
    if args.species:
        L.append(f'Species          {args.species}')
    L += [f'Genes with a region call  {len(d):,}',
          f'Above threshold  {k_tot:,} ({100.0 * k_tot / len(d):.2f}%)', '',
          'By region', '']
    for r in REGION_ORDER:
        lo, hi_ci = res[f'ci_{r}']
        L.append(f"  {r:<12} {res[f'k_{r}']:>4} / {res[f'n_{r}']:>6,}  "
                 f"= {100 * res[f'prop_{r}']:5.2f}%  "
                 f"[{100 * lo:.2f}, {100 * hi_ci:.2f}] (Wilson)")
    L += ['',
          f"  OR (eliminated / retained) = {res['odds_ratio']:.3f} "
          f"[{res['or_lo95']:.3f}, {res['or_hi95']:.3f}]",
          '  (conditional maximum-likelihood estimate, exact 95% CI)',
          f"  Fisher exact P = {res['p']:.3g} ({stars(res['p'])})",
          "  OR > 1 means eliminated genes are enriched.", '']

    if sens:
        L += ['dS sensitivity  (the omega > 1 tail is where estimation noise',
              'concentrates, because Var(omega) grows as dS -> 0)', '',
              f"  {'min dS':>8} {'n genes':>9} {'ret':>6} {'elim':>6} "
              f"{'OR':>7} {'95% CI':>17} {'P':>10} {'q':>10}"]
        for s in sens:
            ci = f"[{s['or_lo95']:.3g}, {s['or_hi95']:.3g}]"
            L.append(f"  {s['min_ds']:>8.2f} {s['n_total']:>9,} "
                     f"{s['k_retained']:>6} {s['k_eliminated']:>6} "
                     f"{s['odds_ratio']:>7.3g} {ci:>17} {s['p']:>10.3g} "
                     f"{s['q']:>10.3g}")
        L.append('')

    L += ['Notes', '',
          '  These are genes whose pairwise dN/dS point estimate exceeds',
          f'  {args.threshold:g}. That is a threshold, not a test. Pairwise yn00',
          '  provides no likelihood-ratio test of omega = 1, the variance of a',
          '  ratio grows without bound as its denominator goes to zero, and with',
          f'  {len(d):,} pairs a tail above 1 is expected under purifying',
          '  selection everywhere.', '',
          '  Write "genes with dN/dS > 1" in the manuscript, not "positively',
          '  selected genes", unless a codon-model test is run on this set and',
          '  reported instead. The dS sensitivity table above is the evidence for',
          '  how much of the tail survives a plausible noise floor.', '',
          '  The odds ratio is oriented eliminated over retained, so OR > 1 is',
          '  enrichment in eliminated regions. The earlier version reported the',
          '  reciprocal and needed prose to say which way it pointed.', '',
          '  Odds ratios are conditional maximum-likelihood estimates with exact',
          '  95% intervals (the inversion of Fisher\'s test), as for every other',
          '  odds ratio in the paper. Before 2026-10-12 this script reported the',
          '  Haldane-Anscombe estimate with a Woolf interval.', '',
          '  ds_diagnostic.* is a check, not a manuscript panel: it shows where',
          '  the above-threshold genes sit in dS and alignment length relative to',
          '  the rest.', '']
    summary = '\n'.join(L)
    with open(os.path.join(args.outdir,
                           'positive_selection_summary.txt'), 'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(), 'python': platform.python_version(),
                   'packages': {'numpy': np.__version__,
                                'pandas': pd.__version__,
                                'matplotlib': matplotlib.__version__,
                                'scipy': scipy_version},
                   'argv': sys.argv, 'args': vars(args),
                   'inputs': {'annotated': fingerprint(args.annotated)},
                   'primary': res, 'ds_sensitivity': sens,
                   'n_genes': int(len(d)), 'n_above_threshold': int(k_tot)},
                  fh, indent=2, default=float)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
