#!/usr/bin/env python3
"""
microsynteny.py -- conservation of local gene order between Ascaris and
P. univalens, for retained, eliminated and CBR-spanning gene pairs
(supplementary fig. S7, table T8).

Definition (unchanged from the May script, so Results' counts stand)
--------------------------------------------------------------------
Genes are ranked by midpoint along each sequence, in each species, over
every mRNA model in the miniprot GFF. For each pair of genes ADJACENT in
Ascaris whose two genes both have a P. univalens ortholog in the rate table
(and both orthologs are placed in the P. univalens GFF), the pair is
CONSERVED if the two orthologs lie on one P. univalens sequence within
--max-distance ranks of each other (default 3, which tolerates a single-gene
insertion or a small tandem duplication). Pairs are classed by the Ascaris
region of their two genes: retained (both retained), eliminated (both
eliminated) or CBR-spanning (one of each); a pair with an unannotated or
ambiguous gene is counted and left out.

Tests
-----
  * Conserved fraction, retained vs eliminated: two-sided Fisher's exact test,
    conditional MLE odds ratio with its exact 95% CI (R's fisher.test). The
    May script computed the test but recorded only the counts, so the P in
    Results (9.7e-51) was calculated by hand; it is now written out.
  * The two comparisons with the 12 or so CBR-spanning pairs, for
    completeness; BH q over the three pairwise tests.
  * Threshold sweep: the same test at rank distances 1, 2, 3, 5, 10, 20 and 50,
    so the choice of 3 can be seen not to drive the result. Panel B draws the
    whole sweep.

Panels (one 185 mm strip, figS7_microsynteny)
---------------------------------------------
  A  % of adjacent pairs conserved, retained / CBR-spanning / eliminated, with
     Wilson 95% intervals; the Fisher test drawn as a bracket.
  B  Cumulative % of pairs whose orthologs lie within rank distance d, against
     d on a log axis, retained and eliminated; the dashed line marks the
     threshold, so each curve's height there is panel A's bar. Pairs whose
     orthologs are on different sequences never enter the curve; the share is
     given in the key. CBR-spanning pairs are too few for a curve.

What changed from the May script
--------------------------------
  * The Fisher P and an odds ratio with CI are reported (summary, tests TSV).
  * Panel B replaces a histogram cut at rank distance 50, which dropped every
    pair beyond it without saying so, and did not show different-sequence
    pairs at all.
  * Ties in midpoint (overlapping models) are broken by start, end and ID;
    the May script left them to an unstable sort, so the gene order could
    depend on the numpy build. Counts can differ from May's by a pair or two
    for that reason only -- the summary prints both.
  * House figure pipeline; n per filtering step; nothing capped silently.

Outputs (in --outdir)
---------------------
  figS7_microsynteny.{svg,pdf,png}
  microsynteny_pairs.tsv            every evaluated adjacent pair
  microsynteny_summary.tsv          T8: per class n, conserved, %, Wilson CI
  microsynteny_threshold_sweep.tsv  per threshold and class, with the Fisher test
  microsynteny_tests.tsv            the pairwise Fisher tests at --max-distance
  microsynteny_summary.txt
  microsynteny.log
  provenance.json

Usage
-----
  python microsynteny.py --annotated results/annotate_regions_20260817/kaks_annotated.tsv \\
    --gff1 As_miniprot.gff --gff2 Pu_miniprot.gff --panel-letters \\
    --outdir results/microsynteny_<stamp>

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
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from scipy.stats import fisher_exact
from scipy import __version__ as scipy_version

# A script that cannot import the shared style must stop, not fall back to
# matplotlib defaults.
from paper_colors import (REGION, NEUTRAL, FONT, LINEWIDTH, apply_style, audit,
                          figure, save, mm)
from region_lib import (benjamini_hochberg, fingerprint, fmt_p, git_commit,
                        parse_gff_models, setup_logging, stars)
from rate_tables import read_rate_table

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Pair classes in drawing order: retained left, eliminated right.
CLASSES = ('retained', 'boundary', 'eliminated')
CLASS_LABEL = {'retained': 'Retained', 'boundary': 'CBR-spanning',
               'eliminated': 'Eliminated'}
CLASS_COLOR = {'retained': REGION['retained'], 'boundary': NEUTRAL['ink_light'],
               'eliminated': REGION['eliminated']}
Z95 = 1.959963984540054


# ---------------------------------------------------------------------------
# Gene order and pairs
# ---------------------------------------------------------------------------

def gene_order(models, label):
    """{gene_id: (sequence, rank)} and {sequence: [gene_id, ...]}.

    Rank is the 0-based position by midpoint along the sequence, over every
    mRNA model in the GFF. Ties (overlapping models with one midpoint) are
    broken by start, end and ID so that the order is the same on any machine.
    """
    m = models.copy()
    n_dup = int(m['gene_id'].duplicated().sum())
    if n_dup:
        logging.warning('  %s: %d repeated model IDs; the first record of each is '
                        'kept, as before', label, n_dup)
        m = m.drop_duplicates('gene_id', keep='first')
    m['midpoint'] = (m['start'] + m['end']) // 2
    m = m.sort_values(['chrom', 'midpoint', 'start', 'end', 'gene_id'],
                      kind='mergesort')
    m['rank'] = m.groupby('chrom', sort=False).cumcount()
    rank_of = dict(zip(m['gene_id'], zip(m['chrom'], m['rank'])))
    order = {c: g['gene_id'].tolist() for c, g in m.groupby('chrom', sort=False)}
    n_tied = int(m.duplicated(['chrom', 'midpoint'], keep=False).sum())
    logging.info('  %s: %s models on %d sequences (%s share a midpoint with '
                 'another model)', label, f'{len(m):,}', len(order),
                 f'{n_tied:,}')
    return rank_of, order, n_tied


def pair_class(ra, rb):
    """retained / eliminated / boundary / unknown for two region labels."""
    if ra == 'retained' and rb == 'retained':
        return 'retained'
    if ra == 'eliminated' and rb == 'eliminated':
        return 'eliminated'
    if {ra, rb} == {'retained', 'eliminated'}:
        return 'boundary'
    return 'unknown'


def adjacent_pairs(order1, ortholog, region, rank2):
    """Every Ascaris-adjacent pair, with the reason it is or is not evaluated.

    Returns (pairs DataFrame of evaluated pairs, step counts). The rank
    DISTANCE in species 2 is stored rather than the conserved flag, so any
    threshold can be applied afterwards.
    """
    rows = []
    steps = {'adjacent model pairs in Ascaris': 0,
             'both genes have an ortholog': 0,
             'both orthologs placed in the P. univalens GFF': 0}
    for chrom, genes in order1.items():
        for ga, gb in zip(genes[:-1], genes[1:]):
            steps['adjacent model pairs in Ascaris'] += 1
            oa, ob = ortholog.get(ga), ortholog.get(gb)
            if oa is None or ob is None:
                continue
            steps['both genes have an ortholog'] += 1
            if oa not in rank2 or ob not in rank2:
                continue
            steps['both orthologs placed in the P. univalens GFF'] += 1
            (ca, ra), (cb, rb) = rank2[oa], rank2[ob]
            same = ca == cb
            rows.append({'geneA': ga, 'geneB': gb, 'chrom1': chrom,
                         'regionA': region.get(ga, 'unknown'),
                         'regionB': region.get(gb, 'unknown'),
                         'orthA': oa, 'orthB': ob,
                         'chrom2A': ca, 'chrom2B': cb, 'same_chrom2': same,
                         'rank_dist2': abs(ra - rb) if same else np.nan})
    df = pd.DataFrame(rows)
    if len(df):
        df['pair_class'] = [pair_class(a, b) for a, b in
                            zip(df['regionA'], df['regionB'])]
    return df, steps


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def wilson(k, n, z=Z95):
    """Wilson score interval for k of n, as proportions."""
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def _conditional_or_fallback(a, b, c, d, level):
    """Conditional MLE odds ratio and exact CI for scipy < 1.10 (no
    scipy.stats.contingency.odds_ratio): a follows Fisher's noncentral
    hypergeometric distribution given the margins; the limits are where the
    observed count sits in either (1 - level)/2 tail, as in R's fisher.test."""
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


def conditional_or(a, b, c, d, level=0.95):
    """Conditional MLE odds ratio and exact CI for [[a, b], [c, d]]."""
    try:
        from scipy.stats.contingency import odds_ratio
    except ImportError:
        mle, lo, hi = _conditional_or_fallback(a, b, c, d, level)
        return float(mle), float(lo), float(hi), 'fallback'
    res = odds_ratio(np.array([[a, b], [c, d]]), kind='conditional')
    ci = res.confidence_interval(confidence_level=level)
    return float(res.statistic), float(ci.low), float(ci.high), 'scipy'


def class_counts(pairs, max_d):
    """Per class: n, conserved, %, Wilson CI, and where the rest went."""
    rows = []
    for cls in CLASSES:
        sub = pairs[pairs['pair_class'] == cls]
        n = len(sub)
        cons = int((sub['same_chrom2'] & (sub['rank_dist2'] <= max_d)).sum())
        other = int((~sub['same_chrom2']).sum())
        lo, hi = wilson(cons, n)
        rows.append({'pair_class': cls, 'max_distance': max_d, 'n_pairs': n,
                     'n_conserved': cons,
                     'pct_conserved': 100 * cons / n if n else np.nan,
                     'pct_lo': 100 * lo, 'pct_hi': 100 * hi,
                     'n_same_seq_beyond': n - cons - other,
                     'n_different_seq': other,
                     'pct_different_seq': 100 * other / n if n else np.nan})
    return pd.DataFrame(rows).set_index('pair_class')


def fisher_pair(tab, g1, g2):
    """Fisher's exact test of conserved fraction, g1 vs g2.

    Table [[g1 conserved, g1 not], [g2 conserved, g2 not]], so the odds ratio
    is the odds of conservation in g1 over the odds in g2.
    """
    a, n1 = int(tab.loc[g1, 'n_conserved']), int(tab.loc[g1, 'n_pairs'])
    c, n2 = int(tab.loc[g2, 'n_conserved']), int(tab.loc[g2, 'n_pairs'])
    b, d = n1 - a, n2 - c
    out = {'group1': g1, 'group2': g2, 'conserved1': a, 'n1': n1,
           'conserved2': c, 'n2': n2, 'odds_ratio_sample': np.nan,
           'odds_ratio': np.nan, 'or_lo': np.nan, 'or_hi': np.nan,
           'or_method': '', 'p': np.nan}
    if min(n1, n2) == 0:
        return out
    or_s, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
    or_c, lo, hi, how = conditional_or(a, b, c, d)
    out.update(odds_ratio_sample=float(or_s), odds_ratio=or_c, or_lo=lo,
               or_hi=hi, or_method=how, p=float(p))
    return out


def num_word(n):
    """'three' for 3, as the text spells small numbers."""
    words = ('zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven',
             'eight', 'nine', 'ten')
    return words[n] if 0 <= n < len(words) else f'{n:,}'


def p_text(p):
    """Plain-text P for the summary (fmt_p writes mathtext for the figure)."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p == 0:
        return '< 1e-300'
    return f'{p:.2g}' if p < 0.001 else f'{p:.3f}'


def p_manuscript(p):
    """P in the manuscript's notation, e.g. 9.7 × 10⁻⁵¹."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p == 0:
        return '< 10⁻³⁰⁰'
    if p >= 0.001:
        return f'{p:.2g}'
    e = int(np.floor(np.log10(p)))
    m = p / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    sup = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    return f'{m:.1f} × 10{str(e).translate(sup)}'


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
    for leg in fig.legends:
        out += list(leg.get_texts())
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


def place_row(fig, axes, shares, gaps_mm, pad_mm=0.8, top_extra_mm=0.0,
              iterations=4):
    """Lay out panels side by side from their measured text; the row fills
    the canvas whatever font the node renders with. `top_extra_mm` reserves a
    band above everything (here, for the panel letters)."""
    W, H = fig.bbox.width, fig.bbox.height
    pad = mm(pad_mm) * fig.dpi
    extra = mm(top_extra_mm) * fig.dpi
    gaps = [mm(g) * fig.dpi for g in gaps_mm]
    shares = np.asarray(shares, float) / np.sum(shares)
    for _ in range(iterations):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        tight = [ax.get_tightbbox(r) for ax in axes]
        boxes = [ax.get_window_extent(r) for ax in axes]
        lefts = [max(0.0, b.x0 - t.x0) for t, b in zip(tight, boxes)]
        rights = [max(0.0, t.x1 - b.x1) for t, b in zip(tight, boxes)]
        top = max(max(0.0, t.y1 - b.y1) for t, b in zip(tight, boxes)) + extra
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


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

def draw_conserved(ax, tab, test, max_d, bar_alpha=0.40):
    """A: % of adjacent pairs conserved per class, Wilson 95% CI, Fisher
    bracket between retained and eliminated."""
    x = np.arange(len(CLASSES))
    for i, cls in enumerate(CLASSES):
        r = tab.loc[cls]
        if not r['n_pairs']:
            continue
        ax.bar(i, r['pct_conserved'], width=0.62,
               color=to_rgba(CLASS_COLOR[cls], bar_alpha), edgecolor='black',
               linewidth=LINEWIDTH['thin'], zorder=2)
        ax.plot([i, i], [r['pct_lo'], r['pct_hi']], color='black',
                lw=LINEWIDTH['thin'], solid_capstyle='butt', zorder=3)
        for y in (r['pct_lo'], r['pct_hi']):
            ax.plot([i - 0.08, i + 0.08], [y, y], color='black',
                    lw=LINEWIDTH['thin'], zorder=3)
        ax.text(i, r['pct_hi'] + 2.5, f"{r['pct_conserved']:.1f}%",
                ha='center', va='bottom', fontsize=FONT['annotation'])
    top = 116.0
    ax.plot([0, 0, 2, 2], [top, top + 3, top + 3, top], color='black',
            lw=LINEWIDTH['thin'], zorder=4)
    ax.text(1, top + 4.5, f"{stars(test['q'])}  {fmt_p(test['p'])}",
            ha='center', va='bottom', fontsize=FONT['annotation'])
    ax.set_ylim(0, 132)
    ax.yaxis.set_major_locator(FixedLocator([0, 25, 50, 75, 100]))
    ax.set_xlim(-0.6, len(CLASSES) - 0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{CLASS_LABEL[c]}\n(n = {int(tab.loc[c, 'n_pairs']):,})"
                        for c in CLASSES])
    ax.tick_params(axis='x', length=0)
    ax.set_ylabel(f'Adjacent pairs conserved (%)\n(orthologs within '
                  f'{max_d} genes)')


def draw_cumulative(ax, pairs, tab, max_d):
    """B: cumulative % of pairs with orthologs within rank distance d.

    The denominator is every pair in the class, so a curve levels off at the
    share whose orthologs share a sequence; the rest are on different
    sequences and are named in the key. Each curve is carried flat to the
    right edge, past its last pair.
    """
    curves, handles, labels = {}, [], []
    d_max = 1.0
    for cls in ('retained', 'eliminated'):
        sub = pairs[pairs['pair_class'] == cls]
        n = len(sub)
        if not n:
            continue
        d = np.sort(sub.loc[sub['same_chrom2'], 'rank_dist2'].to_numpy(float))
        u, c = np.unique(d, return_counts=True)
        y = 100 * np.cumsum(c) / n
        if len(u):
            d_max = max(d_max, float(u.max()))
            # Start at d = 1 with the share already at distance <= 1.
            y1 = y[u <= 1][-1] if (u <= 1).any() else 0.0
            xs = np.concatenate([[1.0], u[u > 1]])
            ys = np.concatenate([[y1], y[u > 1]])
        else:
            xs, ys = np.array([1.0]), np.array([0.0])
        curves[cls] = (xs, ys)
        handles.append(Line2D([], [], color=REGION[cls], lw=LINEWIDTH['data']))
        labels.append(f"{CLASS_LABEL[cls]} (n = {n:,}; "
                      f"{tab.loc[cls, 'pct_different_seq']:.1f}% on another "
                      'sequence)')
    x_right = 10 ** np.ceil(np.log10(max(d_max * 1.0001, 10.0)))
    for cls, (xs, ys) in curves.items():
        ax.step(np.append(xs, x_right), np.append(ys, ys[-1]), where='post',
                color=REGION[cls], lw=LINEWIDTH['data'], zorder=3)
    ax.axvline(max_d, color=NEUTRAL['rule'], lw=LINEWIDTH['rule'],
               ls=(0, (4, 2)), zorder=1)
    ax.text(max_d, 101.5, f' threshold ({max_d})', ha='left', va='bottom',
            fontsize=FONT['annotation'], color=NEUTRAL['rule'])
    ax.set_xscale('log')
    ax.set_xlim(1, x_right)
    ax.set_ylim(0, 108)
    ax.yaxis.set_major_locator(FixedLocator([0, 25, 50, 75, 100]))
    # Plain integers, not 10^k: mathtext would set the exponents in the math
    # font, a second typeface in the SVG.
    decades = 10.0 ** np.arange(0, int(round(np.log10(x_right))) + 1)
    ax.xaxis.set_major_locator(FixedLocator(decades))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:,.0f}'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('Rank distance between the two orthologs (genes)')
    ax.set_ylabel('Pairs within this distance (%)')
    ax.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.0, 1.0),
              ncol=1, frameon=False, fontsize=FONT['legend'],
              borderaxespad=0.3, handlelength=1.6, labelspacing=0.3)
    return {'x_right': float(x_right), 'max_rank_distance': float(d_max)}


def plot_s7(pairs, tab, test, outdir, stem, args):
    """Supplementary fig. S7 as one 185 mm strip."""
    fig = blank_figure('double', args.height_mm)
    h_in = fig.get_size_inches()[1]
    ax_a = fig.add_axes([0.08, 0.2, 0.25, 0.6])
    ax_b = fig.add_axes([0.45, 0.2, 0.50, 0.6])
    draw_conserved(ax_a, tab, test, args.max_distance)
    info = draw_cumulative(ax_b, pairs, tab, args.max_distance)
    top_extra = 4.5 if args.panel_letters else 0.0
    place_row(fig, [ax_a, ax_b], shares=[0.36, 0.64], gaps_mm=[9.0],
              top_extra_mm=top_extra)
    if args.panel_letters:
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        y = 1.0 - mm(0.8) / h_in
        for ax, letter in ((ax_a, 'A'), (ax_b, 'B')):
            x = ax.get_tightbbox(rend).x0 / fig.bbox.width
            fig.text(x, y, letter, ha='left', va='top',
                     fontsize=FONT['panel_letter'], fontweight='bold')
    paths, problems = finish(fig, outdir, stem)
    return paths, problems, info


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def write_summary(path, ctx):
    """Human-readable results, ending with the sentence they support."""
    a = ctx['args']
    tab, tests, sweep = ctx['tab'], ctx['tests'], ctx['sweep']
    main = tests[(tests['group1'] == 'retained') &
                 (tests['group2'] == 'eliminated')].iloc[0]
    L = ['Micro-synteny, Ascaris vs P. univalens (supplementary fig. S7, '
         'table T8)', '=' * 72, '',
         f'Rate table   {a.annotated}',
         f'GFF 1        {a.gff1}',
         f'GFF 2        {a.gff2}',
         f'Conserved    orthologs on one sequence within {a.max_distance} '
         'ranks', '',
         'Pairs at each step']
    for k, v in ctx['steps'].items():
        L.append(f'  {k:<48} {v:>8,}')
    for cls in CLASSES + ('unknown',):
        n = ctx['class_n'].get(cls, 0)
        L.append(f'    {CLASS_LABEL.get(cls, "unannotated or ambiguous gene"):<46} '
                 f'{n:>8,}')
    if ctx['n_dup_orth']:
        L.append(f"  {ctx['n_dup_orth']} P. univalens genes are the ortholog of "
                 'more than one Ascaris gene (see log)')
    L.append('')

    L += ['-' * 72, f'Conserved fraction (rank distance <= {a.max_distance}); '
          'Wilson 95% CI', '-' * 72,
          f"{'class':<14} {'conserved':>10} {'n':>7} {'%':>7} {'95% CI':>15} "
          f"{'other seq %':>12}"]
    for cls in CLASSES:
        r = tab.loc[cls]
        L.append(f"{CLASS_LABEL[cls]:<14} {int(r['n_conserved']):>10,} "
                 f"{int(r['n_pairs']):>7,} {r['pct_conserved']:>7.1f} "
                 f"{r['pct_lo']:>7.1f}-{r['pct_hi']:<7.1f} "
                 f"{r['pct_different_seq']:>12.1f}")
    L += ['', "Fisher's exact tests (two-sided; OR = odds of conservation, "
          'group 1 over group 2;', 'conditional MLE with exact 95% CI; BH q over '
          'the three)']
    for t in tests.itertuples():
        L.append(f'  {CLASS_LABEL[t.group1]:<13} vs {CLASS_LABEL[t.group2]:<13}'
                 f'  OR {t.odds_ratio:7.2f} [{t.or_lo:.2f}, {t.or_hi:.2f}]  '
                 f'P = {p_text(t.p):<9} q = {p_text(t.q)}')
    L.append(f"  (sample odds ratio, retained vs eliminated: "
             f"{main['odds_ratio_sample']:.2f})")
    L.append('')

    L += ['-' * 72, 'Threshold sweep (retained vs eliminated)', '-' * 72,
          f"{'max d':>5} {'retained %':>11} {'CBR-span %':>11} "
          f"{'eliminated %':>13} {'OR':>8} {'P':>10}"]
    for d, block in sweep.groupby('max_distance'):
        b = block.set_index('pair_class')
        L.append(f"{int(d):>5} {b.loc['retained', 'pct_conserved']:>11.1f} "
                 f"{b.loc['boundary', 'pct_conserved']:>11.1f} "
                 f"{b.loc['eliminated', 'pct_conserved']:>13.1f} "
                 f"{b.loc['retained', 'odds_ratio']:>8.2f} "
                 f"{p_text(b.loc['retained', 'p']):>10}")
    worst = sweep[sweep['pair_class'] == 'retained']['p'].max()
    wide = sweep.pivot(index='max_distance', columns='pair_class',
                       values='pct_conserved')
    p_by_d = sweep[sweep['pair_class'] == 'retained'].set_index('max_distance')['p']
    fails = [int(d) for d in wide.index
             if not (wide.loc[d, 'retained'] > wide.loc[d, 'eliminated']
                     and p_by_d.loc[d] < 0.05)]
    L += [f'  Largest P across the sweep: {p_text(worst)}.' +
          (f'  Retained not significantly higher at: {fails}.' if fails else
           '  Retained significantly higher at every threshold.'), '']

    may = ctx['may']
    L += ['-' * 72, 'Against the May run (Results quotes these)', '-' * 72]
    same = all(int(tab.loc[c, 'n_conserved']) == may[c][0] and
               int(tab.loc[c, 'n_pairs']) == may[c][1] for c in CLASSES)
    for c in CLASSES:
        L.append(f"  {CLASS_LABEL[c]:<13} May {may[c][0]:,}/{may[c][1]:,}   now "
                 f"{int(tab.loc[c, 'n_conserved']):,}/{int(tab.loc[c, 'n_pairs']):,}")
    L.append('  identical' if same else
             '  DIFFERENT: if only by a pair or two, the cause is the midpoint '
             'tie-break (see the docstring); otherwise check the inputs.')
    L.append('')

    fi = ctx['figure']
    L += ['-' * 72, 'Drawn (figS7_microsynteny)', '-' * 72,
          '  A: bars = % conserved, whiskers = Wilson 95% CI; asterisks from BH '
          'q, exact P printed.',
          f"  B: x-axis 1 to {fi['x_right']:g} (largest rank distance between "
          f"same-sequence orthologs {fi['max_rank_distance']:g}); every pair "
          'is inside the axis.',
          f"     CBR-spanning pairs (n = {int(tab.loc['boundary', 'n_pairs'])}) "
          'are not drawn in B.']
    L.append('  AUDIT: ' + '; '.join(ctx['problems']) if ctx['problems']
             else '  audit: clean')
    L.append('')

    r, e, bd = tab.loc['retained'], tab.loc['eliminated'], tab.loc['boundary']
    lo_, hi_ = sorted((r['pct_conserved'], e['pct_conserved']))
    between = ('intermediate' if lo_ < bd['pct_conserved'] < hi_ else
               'NOT intermediate (rephrase)')
    L += ['-' * 72, 'What the text can claim', '-' * 72,
          f"  Of adjacent Ascaris gene pairs in which both genes have a "
          f"P. univalens ortholog, {int(r['n_conserved']):,} of "
          f"{int(r['n_pairs']):,} retained pairs ({r['pct_conserved']:.1f}%) "
          f"stay within {num_word(a.max_distance)} genes of each other in "
          f"P. univalens, against {int(e['n_conserved']):,} "
          f"of {int(e['n_pairs']):,} eliminated pairs ({e['pct_conserved']:.1f}%; "
          f"Fisher's exact P = {p_manuscript(main['p'])}), with pairs spanning "
          f"a breakage region {between} at {int(bd['n_conserved'])} of "
          f"{int(bd['n_pairs'])} ({bd['pct_conserved']:.1f}%).",
          f"  Odds ratio {main['odds_ratio']:.1f} (95% CI {main['or_lo']:.1f}-"
          f"{main['or_hi']:.1f})." + (
              f" The direction holds at every threshold from "
              f"{int(sweep['max_distance'].min())} to "
              f"{int(sweep['max_distance'].max())} genes (largest P "
              f"{p_manuscript(worst)})." if not fails else
              f" NOT robust to the threshold: fails at {fails}; say so."), '']
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Micro-synteny conservation, retained vs eliminated '
                    '(supplementary fig. S7, table T8).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', '--kaks-table', dest='annotated',
                   required=True,
                   help='kaks_annotated.tsv from annotate_regions.py '
                        '(ANNOTATED_TABLE in config.sh)')
    p.add_argument('--gff1', required=True, help='Ascaris miniprot GFF')
    p.add_argument('--gff2', required=True, help='P. univalens miniprot GFF')
    p.add_argument('--gene-col1', default='gene1')
    p.add_argument('--gene-col2', default='gene2')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--max-distance', type=int, default=3,
                   help='Rank distance in P. univalens at or below which a pair '
                        'is conserved (default 3)')
    p.add_argument('--sweep', type=int, nargs='+',
                   default=[1, 2, 3, 5, 10, 20, 50],
                   help='Thresholds for the sensitivity sweep')
    p.add_argument('--may-counts', default='4644/4877,8/12,58/125',
                   help='Counts from the May run, retained,CBR-spanning,'
                        'eliminated as conserved/n, for the comparison in the '
                        'summary (default: the values Results quotes)')
    p.add_argument('--height-mm', type=float, default=60.0)
    p.add_argument('--panel-letters', action='store_true',
                   help='Draw A and B (off by default, per the house rule; S7 '
                        'is emitted whole, so its wrapper turns them on)')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.max_distance < 1 or min(args.sweep) < 1:
        p.error('--max-distance and every --sweep value must be >= 1')
    try:
        vals = [tuple(int(v) for v in s.split('/'))
                for s in args.may_counts.split(',')]
        assert len(vals) == 3 and all(len(v) == 2 for v in vals)
        args.may = dict(zip(CLASSES, vals))
    except (ValueError, AssertionError):
        p.error('--may-counts must look like 4644/4877,8/12,58/125')
    if args.outdir is None:
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.outdir = os.path.join('results', f'microsynteny_{stamp}')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'microsynteny')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('microsynteny.py -- supplementary fig. S7, table T8')
    for k in ('annotated', 'gff1', 'gff2'):
        f = getattr(args, k)
        logging.info('  %-9s %s', k, os.path.abspath(f))
        if not os.path.isfile(f):
            sys.exit(f'ERROR: not found: {f}')
    logging.info('  outdir    %s', os.path.abspath(args.outdir))
    apply_style()

    # --- orthologs and regions ------------------------------------------------
    tab_in = read_rate_table(args.annotated, require=())
    for col in (args.gene_col1, args.gene_col2, args.region_col):
        if col not in tab_in.columns:
            sys.exit(f'ERROR: column {col!r} not in {args.annotated}. '
                     f'Columns: {list(tab_in.columns)}')
    tab_in = tab_in.dropna(subset=[args.gene_col1, args.gene_col2])
    g1 = tab_in[args.gene_col1].astype(str)
    g2 = tab_in[args.gene_col2].astype(str)
    n_dup1 = int(g1.duplicated().sum())
    if n_dup1:
        logging.warning('  %d Ascaris genes appear in more than one row; the '
                        'last row is used, as in the May script', n_dup1)
    ortholog = dict(zip(g1, g2))
    region = dict(zip(g1, tab_in[args.region_col].astype(str)))
    dup2 = g2[g2.duplicated(keep=False)]
    n_dup_orth = int(dup2.nunique())
    if n_dup_orth:
        logging.warning('  %d P. univalens genes are listed as the ortholog of '
                        'more than one Ascaris gene, e.g. %s', n_dup_orth,
                        sorted(dup2.unique())[:3])
    logging.info('Ortholog pairs: %s (regions: %s)', f'{len(ortholog):,}',
                 pd.Series(list(region.values())).value_counts().to_dict())

    # --- gene order -------------------------------------------------------------
    rank1, order1, ties1 = gene_order(parse_gff_models(args.gff1), 'Ascaris')
    rank2, _, ties2 = gene_order(parse_gff_models(args.gff2), 'P. univalens')
    miss1 = sum(1 for g in ortholog if g not in rank1)
    miss2 = sum(1 for g in ortholog.values() if g not in rank2)
    logging.info('  orthologs not in the GFFs: Ascaris %s, P. univalens %s',
                 f'{miss1:,}', f'{miss2:,}')
    if miss1 > 0.5 * len(ortholog) or miss2 > 0.5 * len(ortholog):
        sys.exit('ERROR: most rate-table genes are not model IDs in the GFFs. '
                 'Are the GFFs the ones the table was built from?\n'
                 f'  table: {list(ortholog.items())[:2]}\n'
                 f'  GFF 1: {list(rank1)[:2]}\n  GFF 2: {list(rank2)[:2]}')

    # --- pairs --------------------------------------------------------------------
    pairs, steps = adjacent_pairs(order1, ortholog, region, rank2)
    if pairs.empty:
        sys.exit('ERROR: no adjacent pair has both genes with a placed ortholog.')
    class_n = pairs['pair_class'].value_counts().to_dict()
    for k, v in steps.items():
        logging.info('  %-48s %s', k, f'{v:,}')
    logging.info('  classes: %s', class_n)
    pairs['conserved'] = pairs['same_chrom2'] & (pairs['rank_dist2'] <=
                                                 args.max_distance)
    pairs.to_csv(os.path.join(args.outdir, 'microsynteny_pairs.tsv'), sep='\t',
                 index=False)

    tab = class_counts(pairs, args.max_distance)
    tab.reset_index().to_csv(os.path.join(args.outdir,
                                          'microsynteny_summary.tsv'),
                             sep='\t', index=False, float_format='%.6g')
    tests = pd.DataFrame([fisher_pair(tab, 'retained', 'eliminated'),
                          fisher_pair(tab, 'retained', 'boundary'),
                          fisher_pair(tab, 'boundary', 'eliminated')])
    tests['q'] = benjamini_hochberg(tests['p'].to_numpy(float))
    tests.insert(0, 'max_distance', args.max_distance)
    tests.to_csv(os.path.join(args.outdir, 'microsynteny_tests.tsv'), sep='\t',
                 index=False, float_format='%.6g')
    for t in tests.itertuples():
        logging.info('  %s vs %s: OR %.2f [%.2f, %.2f], P %s, q %s', t.group1,
                     t.group2, t.odds_ratio, t.or_lo, t.or_hi, p_text(t.p),
                     p_text(t.q))

    sweep_rows = []
    for d in sorted(set(args.sweep) | {args.max_distance}):
        td = class_counts(pairs, d)
        f = fisher_pair(td, 'retained', 'eliminated')
        for cls in CLASSES:
            row = td.loc[cls].to_dict()
            row['pair_class'] = cls
            if cls == 'retained':
                row.update(odds_ratio=f['odds_ratio'], or_lo=f['or_lo'],
                           or_hi=f['or_hi'], p=f['p'])
            sweep_rows.append(row)
    sweep = pd.DataFrame(sweep_rows)
    sweep.to_csv(os.path.join(args.outdir, 'microsynteny_threshold_sweep.tsv'),
                 sep='\t', index=False, float_format='%.6g')

    # --- figure -------------------------------------------------------------------
    logging.info('S7: micro-synteny strip')
    main_test = tests.iloc[0].to_dict()
    paths, problems, fig_info = plot_s7(pairs, tab, main_test, args.outdir,
                                        'figS7_microsynteny', args)

    write_summary(os.path.join(args.outdir, 'microsynteny_summary.txt'), {
        'args': args, 'steps': steps, 'class_n': class_n, 'tab': tab,
        'tests': tests, 'sweep': sweep, 'may': args.may,
        'n_dup_orth': n_dup_orth, 'figure': fig_info, 'problems': problems})

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
                   'inputs': {k: fingerprint(getattr(args, k))
                              for k in ('annotated', 'gff1', 'gff2')},
                   'steps': steps, 'classes': class_n,
                   'tied_midpoints': {'ascaris': ties1, 'parascaris': ties2},
                   'figure': [os.path.basename(p) for p in paths],
                   'audit': problems},
                  fh, indent=2, default=str)
    logging.info('Summary: %s', os.path.join(args.outdir,
                                             'microsynteny_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
