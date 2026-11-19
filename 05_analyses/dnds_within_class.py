#!/usr/bin/env python3
"""
dnds_within_class.py - does functional class account for the dN/dS excess of
eliminated genes?

Reviewer control for Figure 6C (2026-10-19, results draft note 70). Two of the
73 functional classes tested in 6C are enriched among eliminated genes after
Bonferroni correction (Unknown and Proteasome), so the text can no longer say
that eliminated genes are a functionally random set. If the enriched classes
evolve faster for reasons of their own, the eliminated-gene excess could be
partly a matter of which classes the genes belong to rather than of
elimination. This script asks whether the excess holds within classes, and
draws the answer as Supplementary Figure S12D (2026-10-19), beside S12C on the
bottom row of S12.

What it does
------------
1. Joins the annotated rate table to the expression workbook through the
   position bridge (`--bridge workbook_model_bridge.tsv`), with
   functional_enrichment.py's workbook loader and class column, so every gene
   carries the class Figure 6C counts. The analysis set is the genes with a
   class, the denominators of 6C.

2. Re-runs the 6C class tests (`functional_enrichment.category_tests`, same
   minimum class size) and lists the classes whose frequency differs between
   eliminated and retained genes after Bonferroni correction. The summary
   prints the number tested, which must equal 6C's.

3. Eliminated against retained dN/dS, as a ratio of medians with a percentile
   bootstrap interval (each group resampled separately) and a two-sided
   Mann-Whitney P (`germline_restricted_dnds.ratio_contrast`, the S12C
   estimator), on: every analysed gene (the 1.77-fold headline); the genes
   with a class; those without the Bonferroni classes; those without the
   unknown-type classes (any class whose name matches `--unknown-pattern`,
   default 'unknown|uncharacteri[sz]ed', case-insensitive: "Unknown",
   "Unknown / membrane", "Secreted / Unknown", "Uncharacterized" and the like,
   which the workbook keeps as separate classes); and without either.

4. Class-matched contrast: each eliminated gene is paired, without
   replacement and in random order, with a retained gene of its own class
   (`germline_restricted_dnds.match_once`, as the S12C matched rows). The
   point estimate is the median fold over `--n-match` matchings of the
   observed genes; the interval comes from `--n-match-boot` replicates that
   resample both groups before matching, so it carries sampling as well as
   matching variability. The median of the matched retained genes over the
   median of every retained gene with a class measures how far class
   composition alone moves the retained baseline.

5. van Elteren's stratified rank-sum test (van Elteren 1960; the
   implementation of codon_usage_bias.py, where it stratifies by GC3s):
   genes are ranked within each class and the per-class rank sums combined
   with weights 1/(N + 1). A class contributes when it holds at least
   `--min-per-stratum` genes of each region. Positive Z: eliminated genes
   rank higher.

6. Regression, the Figure 6E and S12 model form
   (`multivariate_insulation.fit_model`: standardized OLS, elimination coded
   0/1): dN/dS, raw and ranked, on elimination status alone, and with class
   fixed effects (one indicator per class with at least `--min-class-genes`
   genes, the largest class as reference, smaller classes pooled into one
   indicator). The change in the elimination coefficient is the share of the
   effect that class membership absorbs.

7. Per class, for the classes with at least `--min-elim-class` eliminated
   genes (the 6C forest set): the eliminated-to-retained fold with interval
   and P, Benjamini-Hochberg across these classes.

8. Retained genes only, per tested class: median dN/dS in the class against
   every other retained gene with a class (fold and Mann-Whitney P, BH across
   classes). This says whether an enriched class evolves faster on its own,
   which is how its enrichment could raise the eliminated median.

Supplementary Figure S12D
-------------------------
A forest of eliminated-to-retained dN/dS folds in the style of S12C, 89 mm wide
(`--d-width`) and as tall as S12C (`--d-height-mm`, which
run_dnds_within_class.sh reads from the newest S12C_germline_contrasts.svg),
so C | D make the bottom row. Rows, top to bottom: the gene sets of step 3
(every gene with a class; without the Bonferroni classes; without the
unknown-type classes; without both) and the class-matched fold; then the
classes of the 6C forest, in 6C's order (odds ratio, largest first) and with
6C's labels. Filled circles, 95% intervals, log2 axis, dashed rule at 1; the
n column (eliminated / retained genes, or pairs) is a table axes, as in 6C.
Row labels are drawn one text per line, right-aligned: matplotlib writes a
multi-line tick label to SVG line by line, left-anchored where the node's font
put each line, so in Illustrator's Arial a wrapped label ends short of the
others (6C's "Transcription machinery", 1.7 mm). The stratified test and the
models are single numbers; they go in the legend, not on the panel.

Writes (into --outdir)
----------------------
    S12D_function_classes.{svg,pdf,png}   the panel above
    class_controls.tsv            rows of steps 3-6: the contrast, its gene set,
                                  n, medians, fold, interval, P
    per_class_contrasts.tsv       step 7
    retained_by_class.tsv         step 8, every tested class, with the 6C test
    class_tests.tsv               the step-2 class tests (as 6C's table)
    class_models.tsv              every coefficient of every model in step 6,
                                  class indicators named
    matched_replicates.tsv        one row per matching replicate
    genes_with_class.tsv          the analysis set: gene, region, dN/dS, class
    dnds_within_class_summary.txt counts at every step, every result, outcome
    provenance.json
    dnds_within_class.log

Usage
-----
    python dnds_within_class.py \\
      --kaks results/annotate_regions_20260817/kaks_annotated.tsv \\
      --annotation Ascaris_genes_15714_annotation_expression.xlsx \\
      --bridge results/workbook_model_bridge_20261017_103009/workbook_model_bridge.tsv \\
      --outdir results/dnds_within_class_20261019_150000

Dependencies: numpy, pandas (+ openpyxl), scipy, matplotlib, paper_colors,
region_lib, rate_tables, functional_enrichment, germline_restricted_dnds,
multivariate_insulation, workbook_model_bridge
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import platform
import re
import sys
import textwrap
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator
from matplotlib.transforms import ScaledTranslation
from scipy.stats import mannwhitneyu, norm, rankdata
from scipy import __version__ as scipy_version

# A script that cannot import the shared style must stop, not fall back to
# defaults: a figure drawn without it would look plausible and be wrong.
from paper_colors import NEUTRAL, FONT, LINEWIDTH, apply_style, audit, save
from region_lib import (benjamini_hochberg, fingerprint, git_commit,
                        setup_logging)
from rate_tables import read_rate_table
# One join, one class column and one set of class tests with Figure 6C, and
# 6C's layout helpers and label style for S12D.
from functional_enrichment import (load_workbook, category_tests, CLASS_COL,
                                   blank_figure, add_axes_mm, side_text,
                                   place_with_table, sentence_case, wrap_class,
                                   offcanvas_text)
# The S12C estimators: ratio of medians with bootstrap interval, and the
# matching pass.
from germline_restricted_dnds import ratio_contrast, match_once
# The Figure 6E / S12 model: standardized OLS, t-based 95% intervals.
from multivariate_insulation import fit_model
from workbook_model_bridge import read_bridge, shared_gene_counts

_HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def van_elteren(values, is_elim, strata, min_per_group=1):
    """Stratified Wilcoxon rank-sum test (van Elteren 1960).

    Ranks within each stratum (here, each functional class) and combines the
    per-stratum rank sums of eliminated genes with the locally best weights
    1/(N + 1), with the tie-corrected variance. The formula is that of
    codon_usage_bias.py's van_elteren; with one stratum it reduces to the
    normal approximation of the Mann-Whitney test without continuity
    correction. A stratum with fewer than `min_per_group` genes of either
    region, or with every value tied, contributes nothing.

    Returns a dict (z, p, n_strata, n_retained, n_eliminated) or None.
    Positive z: eliminated genes rank higher within their strata.
    """
    v = np.asarray(values, float)
    e = np.asarray(is_elim, bool)
    codes, uniq = pd.factorize(pd.Series(strata), sort=True)
    num = den = 0.0
    used = n_r = n_e = 0
    for k in range(len(uniq)):
        m = codes == k
        ret, eli = v[m & ~e], v[m & e]
        n1, n2 = len(ret), len(eli)
        if n1 < max(1, min_per_group) or n2 < max(1, min_per_group):
            continue
        N = n1 + n2
        allv = np.concatenate([ret, eli])
        ranks = rankdata(allv)
        W = ranks[n1:].sum()
        expected = n2 * (N + 1) / 2.0
        _, counts = np.unique(allv, return_counts=True)
        tie = float((counts.astype(float) ** 3 - counts).sum())
        var = (n1 * n2 / 12.0) * ((N + 1) - tie / (N * (N - 1.0)))
        if var <= 0:
            continue
        w = 1.0 / (N + 1.0)
        num += w * (W - expected)
        den += w * w * var
        used += 1
        n_r += n1
        n_e += n2
    if used == 0 or den <= 0:
        return None
    z = num / np.sqrt(den)
    return {'z': float(z), 'p': float(2 * norm.sf(abs(z))), 'n_strata': used,
            'n_retained': n_r, 'n_eliminated': n_e}


def candidates_by_class(e_codes, r_codes):
    """For each eliminated gene, the positions of the retained genes of its
    class (an empty array if the class has none)."""
    order = np.argsort(r_codes, kind='stable')
    sorted_codes = r_codes[order]
    lo = np.searchsorted(sorted_codes, e_codes, side='left')
    hi = np.searchsorted(sorted_codes, e_codes, side='right')
    return [order[a:b] for a, b in zip(lo, hi)]


def class_matched_contrast(e_cls, e_v, r_cls, r_v, n_rep, n_boot, rng,
                           level=0.95):
    """Eliminated against retained dN/dS with every eliminated gene paired to a
    retained gene of its own class, without replacement.

    Observed replicates (`n_rep`) give the point estimate, the spread due to
    matching alone and the composition effect (median of the matched retained
    genes over the median of every retained gene); resample-then-match
    replicates (`n_boot`) give the interval. Returns (summary dict,
    per-replicate frame).
    """
    e_cls, r_cls = np.asarray(e_cls, object), np.asarray(r_cls, object)
    e_v, r_v = np.asarray(e_v, float), np.asarray(r_v, float)
    codes, _ = pd.factorize(pd.Series(np.concatenate([e_cls, r_cls])))
    e_c, r_c = codes[:len(e_cls)], codes[len(e_cls):]
    cands = candidates_by_class(e_c, r_c)
    n_none = int(sum(len(c) == 0 for c in cands))
    ret_median = float(np.median(r_v))
    rows = []
    for rep in range(n_rep):
        pick = match_once(cands, len(r_v), rng)
        ok = pick >= 0
        if ok.sum() < 5:
            continue
        me, mr = float(np.median(e_v[ok])), float(np.median(r_v[pick[ok]]))
        rows.append({'kind': 'observed', 'replicate': rep,
                     'n_matched': int(ok.sum()), 'median_eliminated': me,
                     'median_matched_retained': mr,
                     'fold': me / mr if mr > 0 else np.nan,
                     'composition': mr / ret_median if ret_median > 0 else np.nan,
                     'p': float(mannwhitneyu(e_v[ok], r_v[pick[ok]],
                                             alternative='two-sided')[1])})
    for rep in range(n_boot):
        ie = rng.integers(0, len(e_v), len(e_v))
        ir = rng.integers(0, len(r_v), len(r_v))
        pick = match_once(candidates_by_class(e_c[ie], r_c[ir]), len(ir), rng)
        ok = pick >= 0
        if ok.sum() < 5:
            continue
        me = float(np.median(e_v[ie][ok]))
        mr = float(np.median(r_v[ir][pick[ok]]))
        rows.append({'kind': 'bootstrap', 'replicate': rep,
                     'n_matched': int(ok.sum()), 'median_eliminated': me,
                     'median_matched_retained': mr,
                     'fold': me / mr if mr > 0 else np.nan,
                     'composition': np.nan, 'p': np.nan})
    reps = pd.DataFrame(rows, columns=['kind', 'replicate', 'n_matched',
                                       'median_eliminated',
                                       'median_matched_retained', 'fold',
                                       'composition', 'p'])
    obs = reps[reps['kind'] == 'observed']
    bts = reps[reps['kind'] == 'bootstrap']
    alpha = (1 - level) / 2
    out = {'n_eliminated': len(e_v), 'n_retained_pool': len(r_v),
           'n_no_candidate': n_none, 'retained_median_all': ret_median,
           'n_rep_observed': len(obs), 'n_rep_bootstrap': len(bts),
           'fold': np.nan, 'match_lo': np.nan, 'match_hi': np.nan,
           'ci_lo': np.nan, 'ci_hi': np.nan, 'median_n_matched': np.nan,
           'median_eliminated': np.nan, 'median_matched_retained': np.nan,
           'composition': np.nan, 'composition_lo': np.nan,
           'composition_hi': np.nan, 'pct_rep_p_lt_0.05': np.nan,
           'median_p': np.nan}
    if len(obs):
        f = obs['fold'].dropna()
        cmp_ = obs['composition'].dropna()
        out.update(fold=float(f.median()), match_lo=float(f.quantile(alpha)),
                   match_hi=float(f.quantile(1 - alpha)),
                   median_n_matched=float(obs['n_matched'].median()),
                   median_eliminated=float(obs['median_eliminated'].median()),
                   median_matched_retained=float(
                       obs['median_matched_retained'].median()),
                   composition=float(cmp_.median()),
                   composition_lo=float(cmp_.quantile(alpha)),
                   composition_hi=float(cmp_.quantile(1 - alpha)),
                   median_p=float(obs['p'].median()),
                   **{'pct_rep_p_lt_0.05': float(100 * (obs['p'] < 0.05).mean())})
    if len(bts):
        f = bts['fold'].dropna()
        if len(f):
            out.update(ci_lo=float(f.quantile(alpha)),
                       ci_hi=float(f.quantile(1 - alpha)))
    return out, reps


def class_design(frame, min_class_genes):
    """Class indicator columns for the fixed-effects model.

    One 0/1 column per class with at least `min_class_genes` genes in `frame`,
    except the largest, which is the reference; the remaining small classes
    share one column. Returns (design frame with dN_dS, dN_dS_rank,
    is_eliminated and the indicators; {column: class label}; reference class).
    """
    sizes = frame['cls'].value_counts()
    big = sorted(sizes[sizes >= min_class_genes].index)
    if not big:
        return None, {}, None
    # The largest class is the reference; ties broken by name, so the choice
    # is reproducible. The reference does not change the elimination term.
    top = sizes[sizes == sizes.max()].index
    ref = sorted(top)[0]
    mf = pd.DataFrame({'dN_dS': frame['dN_dS'].to_numpy(float),
                       'is_eliminated': frame['is_eliminated'].to_numpy(int)},
                      index=frame.index)
    mf['dN_dS_rank'] = rankdata(mf['dN_dS'])
    label = {}
    for i, c in enumerate(big):
        if c == ref:
            continue
        col = f'class_{i:03d}'
        mf[col] = (frame['cls'] == c).astype(int).to_numpy()
        label[col] = c
    small = ~frame['cls'].isin(big)
    if small.any():
        mf['class_small'] = small.astype(int).to_numpy()
        label['class_small'] = (f'classes with fewer than {min_class_genes} '
                                f'genes, pooled')
    return mf, label, ref


def fit_quietly(mf, response, binary, label):
    """fit_model with its per-coefficient log lines suppressed (a class model
    has some 70 terms); warnings still reach the log."""
    root = logging.getLogger()
    old = root.level
    root.setLevel(logging.WARNING)
    try:
        return fit_model(mf, response, [], binary, label)
    finally:
        root.setLevel(old)


def median_fold_p(a, b):
    """(median a, median b, fold a/b, Mann-Whitney P) without an interval."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ma = float(np.median(a)) if len(a) else np.nan
    mb = float(np.median(b)) if len(b) else np.nan
    fold = ma / mb if (len(a) and len(b) and mb > 0) else np.nan
    p = (float(mannwhitneyu(a, b, alternative='two-sided')[1])
         if len(a) and len(b) else np.nan)
    return ma, mb, fold, p


# ---------------------------------------------------------------------------
# Figure S12D
# ---------------------------------------------------------------------------

#: Width of a journal column in mm, by the names paper_colors.figure() takes.
WIDTH_MM = {'single': 89.0, 'onehalf': 136.0, 'double': 185.0}
#: Row-label wrap width: 6C's, so a class reads the same in both panels.
LABEL_WRAP = 21


def row_layout(rows, block_gap=0.45):
    """y of each row, top row at 0, and the pitch of each: a wrapped label
    takes more height, as in 6C, and a gap separates the blocks."""
    pitch = [max(1.0, 0.62 * (r['label'].count('\n') + 1) + 0.38) for r in rows]
    ys, y = [], 0.0
    for i, r in enumerate(rows):
        if i:
            y -= (pitch[i - 1] + pitch[i]) / 2
            if r['block'] != rows[i - 1]['block']:
                y -= block_gap
        ys.append(y)
    return np.array(ys), pitch


def function_height_mm(rows):
    """Natural canvas height: 4.2 mm per row unit, as 6C and S12C, plus the key
    line, the x-axis and its label."""
    ys, pitch = row_layout(rows)
    span = ys[0] - ys[-1] + (pitch[0] + pitch[-1]) / 2 + 0.5
    return float(min(150.0, 20.6 + 4.2 * span))


def draw_row_labels(ax, ys, labels, fontsize, pad_pt=3.0):
    """Row labels left of the axes, one right-aligned text per line.

    A multi-line tick label is written to SVG line by line, each line
    left-anchored where the node's font (DejaVu Sans on Isaac) put it, so in
    Illustrator's Arial it ends short of the single-line labels. A single-line
    text with ha='right' is written with text-anchor end and stays aligned in
    any font.
    """
    fig = ax.figure
    color = matplotlib.rcParams['ytick.color']
    line_pt = 1.15 * fontsize
    for y, lab in zip(ys, labels):
        lines = lab.split('\n')
        for j, ln in enumerate(lines):
            dy_pt = ((len(lines) - 1) / 2.0 - j) * line_pt
            tr = ax.get_yaxis_transform() + ScaledTranslation(
                -pad_pt / 72.0, dy_pt / 72.0, fig.dpi_scale_trans)
            ax.text(0.0, y, ln, transform=tr, ha='right', va='center_baseline',
                    fontsize=fontsize, color=color, clip_on=False)


def draw_function_forest(ax, ax_t, rows):
    """S12D: one eliminated-to-retained fold per row, with its 95% interval.

    `rows`: dicts with label, block, fold, lo, hi and n_text. Log2 axis with
    S12C's ticks, dashed rule at 1, arrowheads where an interval runs off the
    axis. The n column is written in the table axes `ax_t`.
    """
    ys, pitch = row_layout(rows)
    ends = [v for r in rows for v in (r['lo'], r['hi'])
            if v is not None and np.isfinite(v) and v > 0]
    nice = [0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
    lo, hi = min(ends + [0.9]), max(ends + [1.1])
    x_lo = max([t for t in nice if t <= lo] or [nice[0]])
    x_hi = min([t for t in nice if t >= hi] or [nice[-1]])
    n_clipped = 0
    for r, yy in zip(rows, ys):
        f, a, b = r['fold'], r['lo'], r['hi']
        if f is None or not np.isfinite(f):
            side_text(ax_t, yy, r['n_text'], dx_pt=0.0, x=0.0)
            continue
        if a is not None and b is not None and np.isfinite(a) and np.isfinite(b):
            ax.plot([max(a, x_lo), min(b, x_hi)], [yy, yy], color=NEUTRAL['ink'],
                    lw=LINEWIDTH['thin'], solid_capstyle='butt', zorder=2)
            for end, lim, mk in ((a, x_lo, '<'), (b, x_hi, '>')):
                if (mk == '<' and end < lim) or (mk == '>' and end > lim):
                    ax.plot([lim], [yy], marker=mk, ms=2.8, color=NEUTRAL['ink'],
                            zorder=3, clip_on=False)
                    n_clipped += 1
        ax.plot([min(max(f, x_lo), x_hi)], [yy], ls='none', marker='o', ms=3.4,
                mfc=NEUTRAL['ink'], mec=NEUTRAL['ink'], mew=0.7, zorder=4)
        side_text(ax_t, yy, r['n_text'], dx_pt=0.0, x=0.0)
    y_hi = ys[0] + pitch[0] / 2 + 0.25
    y_lo = ys[-1] - pitch[-1] / 2 - 0.25
    side_text(ax_t, ys[0] + pitch[0] / 2 + 0.12, 'n', dx_pt=0.0, x=0.0,
              va='bottom')
    # The rule at 1 in data coordinates (S12C): axvline's y runs in axes units,
    # which audit() reads as data units.
    ax.plot([1.0, 1.0], [y_lo, y_hi], color=NEUTRAL['rule'],
            lw=LINEWIDTH['rule'], ls=(0, (4, 2)), zorder=1)
    ax.set_xscale('log', base=2)
    ax.set_xlim(x_lo, x_hi)
    ticks = [t for t in nice if x_lo <= t <= x_hi]
    # Over more than about two and a half doublings the in-between ticks (0.75,
    # 1.5, 3) crowd their neighbours at 89 mm; keep powers of two only.
    if np.log2(x_hi / x_lo) > 2.6:
        ticks = [t for t in ticks if float(np.log2(t)).is_integer()]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f'{t:g}' for t in ticks])
    ax.set_yticks(ys)
    ax.set_yticklabels([''] * len(ys))
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(y_lo, y_hi)
    ax_t.set_ylim(y_lo, y_hi)
    draw_row_labels(ax, ys, [r['label'] for r in rows], FONT['tick'])
    ax.set_xlabel('dN/dS fold change (ratio of medians)')
    h = [Line2D([], [], ls='none', marker='o', ms=3.4, mfc=NEUTRAL['ink'],
                mec=NEUTRAL['ink'], mew=0.7)]
    ax.legend(h, ['Eliminated / retained'], loc='lower left',
              bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=FONT['legend'],
              borderaxespad=0.4, handletextpad=0.3)
    return {'x_range': (x_lo, x_hi), 'n_interval_ends_clipped': n_clipped}


def plot_function_panel(rows, outdir, stem, args):
    """S12D at --d-width, --d-height-mm tall (default: from the rows).

    The boxes are starting values: place_with_table sets the plot and the n
    column from the drawn text, so the panel fills the column in any font.
    Returns (paths, problems, info).
    """
    height = args.d_height_mm or function_height_mm(rows)
    width_mm = WIDTH_MM[args.d_width]
    fig = blank_figure(args.d_width, height)
    ax = add_axes_mm(fig, left=26.0, bottom=11.0, right=17.0, top=7.6)
    ax_t = add_axes_mm(fig, left=width_mm - 16.0, bottom=11.0, right=1.0, top=7.6)
    ax_t.axis('off')
    info = draw_function_forest(ax, ax_t, rows)
    place_with_table(fig, ax, ax_t)
    problems = list(audit(fig, verbose=False)) + offcanvas_text(fig)
    for prob in problems:
        logging.warning('  audit %s: %s', stem, prob)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    info.update(height_mm=float(height), width_mm=width_mm,
                height_source='--d-height-mm' if args.d_height_mm else 'rows')
    return paths, problems, info


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Does functional class account for the dN/dS excess of '
                    'eliminated genes? Reviewer control for Figure 6C.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--kaks', required=True,
                   help='Annotated rate table (ANNOTATED_TABLE, from '
                        'annotate_regions.py)')
    p.add_argument('--annotation', required=True,
                   help='Ascaris_genes_15714_annotation_expression.xlsx')
    p.add_argument('--bridge', required=True,
                   help='workbook_model_bridge.tsv (workbook_model_bridge.py): '
                        'model ID -> workbook gene, by position, as Figure 6')
    p.add_argument('--gene-col', default='gene1')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--min-class-genes', type=int, default=5,
                   help='A class is tested, and gets its own indicator in the '
                        'model, if it holds at least this many genes with a '
                        'class (default 5, as Figure 6C)')
    p.add_argument('--bonferroni-alpha', type=float, default=0.05,
                   help='Classes with a Bonferroni-adjusted P below this are '
                        'the ones excluded in the first exclusion row (default '
                        '0.05)')
    p.add_argument('--unknown-pattern', default=r'unknown|uncharacteri[sz]ed',
                   help='Regex, case-insensitive, marking the unknown-type '
                        "classes (default 'unknown|uncharacteri[sz]ed')")
    p.add_argument('--min-elim-class', type=int, default=10,
                   help='Per-class contrasts for classes with at least this '
                        'many eliminated genes (default 10, the 6C forest set)')
    p.add_argument('--min-group', type=int, default=5,
                   help='A contrast is computed only if both groups hold at '
                        'least this many genes (default 5)')
    p.add_argument('--min-per-stratum', type=int, default=1,
                   help='van Elteren: a class contributes if it holds at least '
                        'this many genes of each region (default 1)')
    p.add_argument('--n-boot', type=int, default=2000,
                   help='Bootstrap replicates per contrast (default 2000)')
    p.add_argument('--n-match', type=int, default=500,
                   help='Matchings of the observed genes (default 500)')
    p.add_argument('--n-match-boot', type=int, default=1000,
                   help='Resample-then-match replicates for the matched '
                        'interval (default 1000)')
    p.add_argument('--d-width', choices=['single', 'onehalf', 'double'],
                   default='single',
                   help='Width of S12D: 89, 136 or 185 mm (default single, '
                        'beside S12C on the bottom row)')
    p.add_argument('--d-height-mm', type=float, default=None,
                   help='Height of S12D (default: from the rows; the wrapper '
                        "passes S12C's height)")
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None,
                   help='Default: results/dnds_within_class_<YYYYmmdd_HHMMSS>')
    args = p.parse_args()
    try:
        args.unknown_re = re.compile(args.unknown_pattern, re.I)
    except re.error as exc:
        sys.exit(f'ERROR: --unknown-pattern {args.unknown_pattern!r}: {exc}')
    for name in ('min_class_genes', 'min_elim_class', 'min_group',
                 'min_per_stratum'):
        if getattr(args, name) < 1:
            sys.exit(f"ERROR: --{name.replace('_', '-')} must be at least 1")
    for name in ('n_boot', 'n_match', 'n_match_boot'):
        if getattr(args, name) < 100:
            sys.exit(f"ERROR: --{name.replace('_', '-')} must be at least 100")
    if not 0 < args.bonferroni_alpha < 1:
        sys.exit('ERROR: --bonferroni-alpha must lie between 0 and 1')
    if args.d_height_mm is not None and not 30.0 <= args.d_height_mm <= 150.0:
        sys.exit(f'ERROR: --d-height-mm {args.d_height_mm:g} is outside 30-150 mm')
    if args.outdir is None:
        # Seconds in the stamp: two runs on one day must not share a folder.
        args.outdir = os.path.join(
            'results', 'dnds_within_class_'
            + datetime.now().strftime('%Y%m%d_%H%M%S'))
    return args


def fmt_fold(r):
    """'1.80 [1.52, 2.11]' or 'NA (reason)'."""
    if r is None or not np.isfinite(r.get('fold', np.nan)):
        return 'NA' + (f" ({r['note']})" if r is not None and r.get('note') else '')
    ci = (f" [{r['ci_lo']:.2f}, {r['ci_hi']:.2f}]"
          if np.isfinite(r.get('ci_lo', np.nan)) else ' [no interval]')
    return f"{r['fold']:.2f}{ci}"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'dnds_within_class')
    for noisy in ('fontTools', 'matplotlib.font_manager'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('dnds_within_class.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()
    # Separate streams, so the matching does not depend on how many bootstrap
    # replicates the contrasts drew, and the reverse.
    rng_boot, rng_match = (np.random.default_rng(s) for s in
                           np.random.SeedSequence(args.seed).spawn(2))

    # --- rate table ---------------------------------------------------------
    rates = read_rate_table(args.kaks, require=('dN_dS',))
    for col in (args.gene_col, args.region_col):
        if col not in rates.columns:
            sys.exit(f'ERROR: column {col!r} not in {args.kaks}\n'
                     f'  columns: {list(rates.columns)}')
    rates = rates.rename(columns={args.gene_col: 'gene1',
                                  args.region_col: 'region1'})
    rates['gene1'] = rates['gene1'].astype(str)
    rates['dN_dS'] = pd.to_numeric(rates['dN_dS'], errors='coerce')
    n_rows = len(rates)
    df = rates[rates['region1'].isin(['retained', 'eliminated'])].copy()
    logging.info('Rate table: %s rows, %s retained or eliminated',
                 f'{n_rows:,}', f'{len(df):,}')
    if df['region1'].nunique() < 2:
        sys.exit('ERROR: the rate table needs both retained and eliminated rows.')

    # --- join, as Figure 6 -----------------------------------------------------
    name_of, bridge_info = read_bridge(args.bridge)
    df['annot_name'] = df['gene1'].map(name_of)
    wb = load_workbook(args.annotation)
    if CLASS_COL not in wb.columns:
        sys.exit(f'ERROR: workbook has no {CLASS_COL!r} column; there are no '
                 f'functional classes to test.\n  columns: {list(wb.columns)}')
    df = df.merge(wb[['name', CLASS_COL]], left_on='annot_name',
                  right_on='name', how='left')
    df['cls'] = df[CLASS_COL].where(df[CLASS_COL].notna())
    df.loc[df['cls'].notna(), 'cls'] = (df.loc[df['cls'].notna(), 'cls']
                                        .astype(str).str.strip())
    df['is_eliminated'] = df['region1'].eq('eliminated').astype(int)

    def counts(mask):
        return (int((mask & df['region1'].eq('retained')).sum()),
                int((mask & df['region1'].eq('eliminated')).sum()))

    n_steps = {
        'analysed (retained or eliminated)': counts(df['gene1'].notna()),
        'given a workbook gene (bridge)': counts(df['annot_name'].notna()),
        'matched to the workbook': counts(df['name'].notna()),
        'with a functional class': counts(df['cls'].notna()),
        'with a class and dN/dS': counts(df['cls'].notna() & df['dN_dS'].notna()),
    }
    for k, (r, e) in n_steps.items():
        logging.info('  %-34s retained %s, eliminated %s', k, f'{r:,}', f'{e:,}')
    shared = shared_gene_counts(df['annot_name'], df['region1'])
    if n_steps['with a functional class'] == (0, 0):
        sys.exit('ERROR: no analysed gene has a functional class; check the '
                 'bridge against the rate table and the workbook.')

    # Class tests on every gene with a class, as Figure 6C counts them; the dN/dS
    # analyses on those that also have dN/dS (every analysed gene does).
    has_class = df[df['cls'].notna()].copy()
    hc = has_class[has_class['dN_dS'].notna()].copy()
    if hc['region1'].nunique() < 2:
        sys.exit('ERROR: the genes with a class are all of one region.')
    hc[['gene1', 'region1', 'dN_dS', 'cls']].rename(columns={'cls': 'class'}).to_csv(
        os.path.join(args.outdir, 'genes_with_class.tsv'), sep='\t',
        index=False, float_format='%.6g')

    # --- the 6C class tests -----------------------------------------------------
    # category_tests reads its category column by name; give it the cleaned one.
    sizes = has_class['cls'].value_counts()
    tested = sorted(sizes[sizes >= args.min_class_genes].index)
    tests = category_tests(has_class, 'cls', tested, 'class')
    tests = tests.sort_values('p').reset_index(drop=True)
    tests.to_csv(os.path.join(args.outdir, 'class_tests.tsv'), sep='\t',
                 index=False, float_format='%.6g')
    differ = tests[tests['p_bonferroni'] < args.bonferroni_alpha]
    sig_classes = list(differ['category'])
    unknown_classes = sorted(c for c in sizes.index if args.unknown_re.search(c))
    logging.info('  %d classes tested; Bonferroni P < %g: %s', len(tests),
                 args.bonferroni_alpha, sig_classes or 'none')
    logging.info('  unknown-type classes (pattern %r): %d', args.unknown_pattern,
                 len(unknown_classes))

    # --- step 3: contrasts on gene sets ---------------------------------------
    def contrast(frame, key, label):
        a = frame.loc[frame['region1'] == 'eliminated', 'dN_dS'].to_numpy(float)
        b = frame.loc[frame['region1'] == 'retained', 'dN_dS'].to_numpy(float)
        res = ratio_contrast(a, b, rng_boot, args.n_boot, args.min_group)
        res.update(key=key, label=label, method='ratio of medians, bootstrap, '
                   'Mann-Whitney')
        logging.info('  %-24s %s  P = %.3g  (n %d / %d)', key, fmt_fold(res),
                     res['p'], res['n_num'], res['n_den'])
        return res

    all_genes = df[df['dN_dS'].notna()]
    in_sig = hc['cls'].isin(sig_classes)
    in_unk = hc['cls'].isin(unknown_classes)
    rows = [contrast(all_genes, 'all_analysed', 'every analysed gene'),
            contrast(hc, 'with_class', 'genes with a functional class'),
            contrast(hc[~in_sig], 'without_bonferroni',
                     'without the classes that differ after Bonferroni ('
                     + ('; '.join(sig_classes) if sig_classes else 'none') + ')'),
            contrast(hc[~in_unk], 'without_unknown_type',
                     f'without the {len(unknown_classes)} unknown-type classes'),
            contrast(hc[~(in_sig | in_unk)], 'without_either',
                     'without both sets')]

    # --- step 4: class-matched contrast ---------------------------------------
    e_sel, r_sel = hc['region1'] == 'eliminated', hc['region1'] == 'retained'
    matched, reps = class_matched_contrast(
        hc.loc[e_sel, 'cls'], hc.loc[e_sel, 'dN_dS'], hc.loc[r_sel, 'cls'],
        hc.loc[r_sel, 'dN_dS'], args.n_match, args.n_match_boot, rng_match)
    reps.to_csv(os.path.join(args.outdir, 'matched_replicates.tsv'), sep='\t',
                index=False, float_format='%.6g')
    logging.info('  class-matched            %s  (%s of %d eliminated genes '
                 'matched, median)', fmt_fold(matched),
                 f"{matched['median_n_matched']:.0f}", matched['n_eliminated'])
    rows.append({'key': 'class_matched',
                 'label': 'each eliminated gene against a retained gene of its '
                          'own class',
                 'method': f"median fold over {matched['n_rep_observed']} "
                           f"matchings; interval from "
                           f"{matched['n_rep_bootstrap']} resample-then-match "
                           f"replicates",
                 'n_num': matched['n_eliminated'],
                 'n_den': matched['n_retained_pool'],
                 'median_num': matched['median_eliminated'],
                 'median_den': matched['median_matched_retained'],
                 'fold': matched['fold'], 'ci_lo': matched['ci_lo'],
                 'ci_hi': matched['ci_hi'], 'p': matched['median_p'],
                 'note': f"median P over matchings; P < 0.05 in "
                         f"{matched['pct_rep_p_lt_0.05']:.0f}% of them"})

    # --- step 5: van Elteren ----------------------------------------------------
    ve = van_elteren(hc['dN_dS'], hc['is_eliminated'], hc['cls'],
                     args.min_per_stratum)
    if ve is None:
        logging.warning('  van Elteren: no class holds genes of both regions')
    else:
        logging.info('  van Elteren              Z = %+.2f, P = %.3g (%d classes, '
                     '%d retained / %d eliminated genes)', ve['z'], ve['p'],
                     ve['n_strata'], ve['n_retained'], ve['n_eliminated'])
        rows.append({'key': 'van_elteren',
                     'label': 'ranks within each class, combined across classes',
                     'method': f"van Elteren stratified rank-sum, classes with "
                               f">= {args.min_per_stratum} gene(s) of each region",
                     'n_num': ve['n_eliminated'], 'n_den': ve['n_retained'],
                     'z': ve['z'], 'p': ve['p'],
                     'note': f"{ve['n_strata']} classes contribute"})

    # --- step 6: regression with class fixed effects --------------------------
    mf, design_label, ref_class = class_design(hc, args.min_class_genes)
    model_rows, model_summ = [], {}
    if mf is None:
        logging.warning('  no class holds %d genes; no fixed-effects model',
                        args.min_class_genes)
    else:
        ind = list(design_label)
        for response in ('dN_dS', 'dN_dS_rank'):
            for name, binary in (('base', ['is_eliminated']),
                                 ('class', ['is_eliminated'] + ind)):
                m = fit_quietly(mf, response, binary, f'{response}:{name}')
                if m is None:
                    continue
                for nm, b, se, lo, hi, pv in zip(m['names'], m['beta'], m['se'],
                                                 m['ci_lo'], m['ci_hi'], m['p']):
                    model_rows.append({'response': response, 'model': name,
                                       'term': nm,
                                       'class': design_label.get(nm, ''),
                                       'beta': b, 'se': se, 'ci_lo': lo,
                                       'ci_hi': hi, 'p': pv, 'n': m['n'],
                                       'n_retained': m['n_retained'],
                                       'n_eliminated': m['n_eliminated'],
                                       'r2': m['r2'],
                                       'n_class_terms': len(m['names']) - 2})
                i = m['names'].index('is_eliminated')
                model_summ[(response, name)] = {
                    'beta': float(m['beta'][i]), 'ci_lo': float(m['ci_lo'][i]),
                    'ci_hi': float(m['ci_hi'][i]), 'p': float(m['p'][i]),
                    'r2': float(m['r2']), 'n': int(m['n']),
                    'n_class_terms': len(m['names']) - 2}
                logging.info('  model %-10s %-5s elimination beta = %+.3f '
                             '[%+.3f, %+.3f], P = %.3g (%d class terms, R2 %.4f)',
                             response, name, m['beta'][i], m['ci_lo'][i],
                             m['ci_hi'][i], m['p'][i], len(m['names']) - 2,
                             m['r2'])
    pd.DataFrame(model_rows).to_csv(os.path.join(args.outdir, 'class_models.tsv'),
                                    sep='\t', index=False, float_format='%.6g')
    for response in ('dN_dS', 'dN_dS_rank'):
        b0, b1 = model_summ.get((response, 'base')), model_summ.get((response, 'class'))
        if b0 and b1:
            rows.append({'key': f'model_{response}',
                         'label': 'elimination coefficient, alone -> with class '
                                  'fixed effects',
                         'method': 'standardized OLS (Figure 6E form); '
                                   f'{b1["n_class_terms"]} class indicators, '
                                   f'reference {ref_class!r}',
                         'n_num': int(mf['is_eliminated'].sum()),
                         'n_den': int((mf['is_eliminated'] == 0).sum()),
                         'beta_base': b0['beta'], 'beta_class': b1['beta'],
                         'ci_lo': b1['ci_lo'], 'ci_hi': b1['ci_hi'],
                         'p': b1['p'],
                         'note': f"share absorbed by class "
                                 f"{100 * (1 - b1['beta'] / b0['beta']):.0f}%"
                                 if b0['beta'] else ''})

    controls = pd.DataFrame(rows)
    order = ['key', 'label', 'method', 'n_num', 'n_den', 'median_num',
             'median_den', 'fold', 'ci_lo', 'ci_hi', 'U', 'p', 'z', 'beta_base',
             'beta_class', 'n_boot_dropped', 'note']
    controls = controls.reindex(columns=order).rename(
        columns={'n_num': 'n_eliminated', 'n_den': 'n_retained',
                 'median_num': 'median_eliminated',
                 'median_den': 'median_retained'})
    controls.to_csv(os.path.join(args.outdir, 'class_controls.tsv'), sep='\t',
                    index=False, float_format='%.6g')

    # --- step 7: per class ------------------------------------------------------
    n_elim_cls = hc.loc[e_sel, 'cls'].value_counts()
    per_classes = sorted(n_elim_cls[n_elim_cls >= args.min_elim_class].index,
                         key=lambda c: (-n_elim_cls[c], c))
    per_rows = []
    for c in per_classes:
        sub = hc[hc['cls'] == c]
        a = sub.loc[sub['region1'] == 'eliminated', 'dN_dS'].to_numpy(float)
        b = sub.loc[sub['region1'] == 'retained', 'dN_dS'].to_numpy(float)
        res = ratio_contrast(a, b, rng_boot, args.n_boot, args.min_group)
        res.update({'class': c,
                    'bonferroni_class': c in sig_classes,
                    'unknown_type': c in unknown_classes})
        per_rows.append(res)
    per = pd.DataFrame(per_rows)
    if len(per):
        ok = per['p'].notna()
        per['q_bh'] = np.nan
        if ok.any():
            per.loc[ok, 'q_bh'] = benjamini_hochberg(per.loc[ok, 'p'].to_numpy(float))
        per = per.rename(columns={'n_num': 'n_eliminated', 'n_den': 'n_retained',
                                  'median_num': 'median_eliminated',
                                  'median_den': 'median_retained'})
        cols = ['class', 'n_eliminated', 'n_retained', 'median_eliminated',
                'median_retained', 'fold', 'ci_lo', 'ci_hi', 'U', 'p', 'q_bh',
                'bonferroni_class', 'unknown_type', 'n_boot_dropped', 'note']
        per = per[cols]
    per.to_csv(os.path.join(args.outdir, 'per_class_contrasts.tsv'), sep='\t',
               index=False, float_format='%.6g')

    # --- step 8: retained genes, each tested class against the rest -------------
    ret = hc[r_sel]
    rb_rows = []
    for c in tested:
        inn = ret.loc[ret['cls'] == c, 'dN_dS'].to_numpy(float)
        out = ret.loc[ret['cls'] != c, 'dN_dS'].to_numpy(float)
        mi, mo, fold, pv = median_fold_p(inn, out)
        rb_rows.append({'class': c, 'n_retained_in': len(inn),
                        'n_retained_other': len(out), 'median_in': mi,
                        'median_other': mo, 'fold': fold, 'p': pv})
    rb = pd.DataFrame(rb_rows)
    if len(rb):
        ok = rb['p'].notna()
        rb['q_bh'] = np.nan
        if ok.any():
            rb.loc[ok, 'q_bh'] = benjamini_hochberg(rb.loc[ok, 'p'].to_numpy(float))
        rb = rb.merge(tests[['category', 'n_eliminated_in', 'odds_ratio',
                             'p_bonferroni']].rename(
                          columns={'category': 'class',
                                   'n_eliminated_in': 'n_eliminated_in_class',
                                   'odds_ratio': 'odds_ratio_6C',
                                   'p_bonferroni': 'p_bonferroni_6C'}),
                      on='class', how='left')
        rb['unknown_type'] = rb['class'].isin(unknown_classes)
        rb = rb.sort_values('p').reset_index(drop=True)
    rb.to_csv(os.path.join(args.outdir, 'retained_by_class.tsv'), sep='\t',
              index=False, float_format='%.6g')

    # --- S12D ---------------------------------------------------------------------
    row_of = {r['key']: r for r in rows}
    d_rows = []
    for key, label in (('with_class', 'All with a class'),
                       ('without_bonferroni', 'Without enriched'),
                       ('without_unknown_type', 'Without unknown-type'),
                       ('without_either', 'Without both')):
        r = row_of[key]
        d_rows.append({'key': key, 'label': label, 'block': 'sets',
                       'fold': r['fold'], 'lo': r['ci_lo'], 'hi': r['ci_hi'],
                       'n_text': f"{int(r['n_num']):,} / {int(r['n_den']):,}"})
    d_rows.append({'key': 'class_matched', 'label': 'Class-matched',
                   'block': 'sets', 'fold': matched['fold'],
                   'lo': matched['ci_lo'], 'hi': matched['ci_hi'],
                   'n_text': (f"{int(round(matched['median_n_matched'])):,} pairs"
                              if np.isfinite(matched['median_n_matched']) else '')})
    # The 6C forest set in 6C's order (plot_class_forest: odds ratio, largest
    # first), with 6C's labels.
    forest = (tests[tests['n_eliminated_in'] >= args.min_elim_class]
              .sort_values('odds_ratio', ascending=False))
    per_of = ({r['class']: r for r in per.to_dict('records')} if len(per) else {})
    for c in forest['category']:
        r = per_of.get(c)
        d_rows.append({'key': f'class:{c}',
                       'label': wrap_class(sentence_case(c), LABEL_WRAP),
                       'block': 'classes',
                       'fold': r['fold'] if r else np.nan,
                       'lo': r['ci_lo'] if r else np.nan,
                       'hi': r['ci_hi'] if r else np.nan,
                       'n_text': (f"{int(r['n_eliminated']):,} / "
                                  f"{int(r['n_retained']):,}") if r else ''})
    logging.info('S12D: dN/dS folds by gene set and by class')
    _, d_problems, d_info = plot_function_panel(d_rows, args.outdir,
                                                'S12D_function_classes', args)

    # --- outcome ------------------------------------------------------------------
    excl_keys = ['without_bonferroni', 'without_unknown_type', 'without_either']
    excl_above = [np.isfinite(row_of[k].get('ci_lo', np.nan))
                  and row_of[k]['ci_lo'] > 1 for k in excl_keys]
    matched_above = np.isfinite(matched['ci_lo']) and matched['ci_lo'] > 1
    ve_sig = ve is not None and ve['p'] < 0.05 and ve['z'] > 0
    # Composition: does pairing within class move the retained baseline? The
    # spread is over matchings only, so this is descriptive, not a test.
    comp_part = (np.isfinite(matched['composition_lo'])
                 and matched['composition_lo'] > 1)
    if all(excl_above) and matched_above and ve_sig and not comp_part:
        reading = ('the excess holds within functional classes and without the '
                   'enriched and unknown-type classes; class composition does '
                   'not account for it.')
    elif all(excl_above) and matched_above and ve_sig:
        reading = (f"the excess holds within functional classes and without the "
                   f"enriched and unknown-type classes; class composition "
                   f"accounts for part of it (class-matched fold "
                   f"{matched['fold']:.2f} against {row_of['with_class']['fold']:.2f} "
                   f"unmatched; composition alone x{matched['composition']:.2f}).")
    elif matched_above or ve_sig:
        reading = ('the excess holds within classes by at least one stratified '
                   'test, but not by every control: see the rows above.')
    else:
        reading = ('within classes the excess is not resolved; class '
                   'composition may account for part or all of it.')

    # --- summary --------------------------------------------------------------------
    def fold_line(r):
        return (f"{fmt_fold(r):<22} P = {r['p']:<10.3g} "
                f"n {int(r['n_num']):,} / {int(r['n_den']):,}   "
                f"medians {r['median_num']:.4f} / {r['median_den']:.4f}")

    L = ['dN/dS within functional classes, eliminated vs retained', '',
         f'Rate table   {os.path.basename(args.kaks)}, {n_rows:,} rows',
         f'Workbook     {os.path.basename(args.annotation)}, {len(wb):,} genes',
         f"Join         position bridge {os.path.basename(args.bridge)}: "
         f"{bridge_info['assigned']:,} of {bridge_info['models']:,} models given "
         f"a workbook gene (as Figure 6)",
         f'Class column {CLASS_COL!r}', '',
         'Genes at each step                    retained  eliminated']
    for k, (r, e) in n_steps.items():
        L.append(f'  {k:<34} {r:>8,}  {e:>10,}')
    L += [f'  Sharing their workbook gene with another analysed gene: retained '
          f'{shared[0]:,}, eliminated {shared[1]:,} (each keeps its class, as '
          f'in Figure 6).', '',
          f'Class tests (as Figure 6C): {len(tests)} classes with >= '
          f'{args.min_class_genes} genes; Bonferroni P < {args.bonferroni_alpha:g}: '
          f'{len(differ)}']
    for r in differ.itertuples():
        L.append(f'    {r.category}: {r.n_eliminated_in} / {r.n_retained_in:,} '
                 f'genes, OR {r.odds_ratio:.2f} [{r.or_ci_lo:.2f}, '
                 f'{r.or_ci_hi:.2f}], Bonferroni P = {r.p_bonferroni:.3g} '
                 f"({'enriched' if r.odds_ratio > 1 else 'depleted'})")
    unk_in_set = hc['cls'].isin(unknown_classes)
    L += [f'Unknown-type classes (pattern {args.unknown_pattern!r}, any case): '
          f'{len(unknown_classes)} present, holding '
          f"{int((unk_in_set & r_sel).sum()):,} retained and "
          f"{int((unk_in_set & e_sel).sum()):,} eliminated genes:",
          *(textwrap.wrap('; '.join(unknown_classes), width=96,
                          initial_indent='    ', subsequent_indent='    ')
            if unknown_classes else ['    none']),
          '',
          'Eliminated / retained dN/dS (fold = ratio of medians [95% bootstrap '
          f'CI, {args.n_boot} resamples]; two-sided Mann-Whitney P)']
    for k in ['all_analysed', 'with_class'] + excl_keys:
        r = row_of[k]
        L.append(f"  {k:<22} {fold_line(r)}")
        L.append(f"  {'':<22} {r['label']}")
    L += ['', 'Class-matched (each eliminated gene against a retained gene of its '
              'own class, without replacement)',
          f"  fold {matched['fold']:.3f}; interval {matched['ci_lo']:.3f}-"
          f"{matched['ci_hi']:.3f} ({matched['n_rep_bootstrap']} resample-then-"
          f"match replicates); matching alone {matched['match_lo']:.3f}-"
          f"{matched['match_hi']:.3f} ({matched['n_rep_observed']} matchings)",
          f"  {matched['median_n_matched']:.0f} of {matched['n_eliminated']} "
          f"eliminated genes matched (median); {matched['n_no_candidate']} have no "
          f"retained gene in their class",
          f"  medians {matched['median_eliminated']:.4f} (eliminated) / "
          f"{matched['median_matched_retained']:.4f} (matched retained); "
          f"Mann-Whitney P < 0.05 in {matched['pct_rep_p_lt_0.05']:.0f}% of "
          f"matchings (median P {matched['median_p']:.3g})",
          f"  composition alone: matched retained median / every retained gene "
          f"with a class ({matched['retained_median_all']:.4f}) = "
          f"{matched['composition']:.3f} (matching spread "
          f"{matched['composition_lo']:.3f}-{matched['composition_hi']:.3f})", '']
    if ve is not None:
        L += [f"van Elteren stratified rank-sum across classes: Z = {ve['z']:+.2f}, "
              f"P = {ve['p']:.3g}; {ve['n_strata']} classes with >= "
              f"{args.min_per_stratum} gene(s) of each region, holding "
              f"{ve['n_retained']:,} retained and {ve['n_eliminated']:,} "
              f"eliminated genes (positive Z: eliminated rank higher)", '']
    else:
        L += ['van Elteren: not computed (no class holds genes of both regions)', '']
    if model_summ:
        L += ['Models (standardized OLS, as Figure 6E; elimination coded 0/1, so '
              'its coefficient is', '  the eliminated-minus-retained difference in '
              f'SD of the response). Reference class {ref_class!r}.',
              '  response    model   beta_elim [95% CI]          P          '
              'class terms  R2']
        for (response, name), s in model_summ.items():
            L.append(f"  {response:<11} {name:<7} {s['beta']:+.3f} "
                     f"[{s['ci_lo']:+.3f}, {s['ci_hi']:+.3f}]   {s['p']:<10.3g} "
                     f"{s['n_class_terms']:>11}  {s['r2']:.4f}")
        for response in ('dN_dS', 'dN_dS_rank'):
            b0 = model_summ.get((response, 'base'))
            b1 = model_summ.get((response, 'class'))
            if b0 and b1 and b0['beta']:
                L.append(f"  {response}: class fixed effects change the elimination "
                         f"coefficient by {100 * (b1['beta'] / b0['beta'] - 1):+.1f}%")
        L.append('')
    if len(per):
        L += [f'Per class, classes with >= {args.min_elim_class} eliminated genes '
              '(the 6C forest set); BH q across them',
              '  n elim / ret    fold [95% CI]           P          q          class']
        for r in per.to_dict('records'):
            ci = (f"[{r['ci_lo']:.2f}, {r['ci_hi']:.2f}]"
                  if np.isfinite(r['ci_lo']) else '[no interval]')
            flag = (' (Bonferroni class)' if r['bonferroni_class'] else '') + \
                   (' (unknown-type)' if r['unknown_type'] else '')
            L.append(f"  {int(r['n_eliminated']):5d} / {int(r['n_retained']):<7,} "
                     f"{r['fold']:5.2f} {ci:<16} {r['p']:<10.3g} "
                     f"{r['q_bh']:<10.3g} {r['class']}{flag}")
        L.append('')
    if len(rb):
        L += ['Retained genes only: each tested class against every other retained '
              'gene with a class',
              '  (does the class evolve faster on its own?) Bonferroni and '
              'unknown-type classes, then the',
              '  five smallest P:']
        focus = rb[rb['class'].isin(sig_classes) | rb['unknown_type']]
        rest = rb[~rb.index.isin(focus.index)].head(5)
        for r in pd.concat([focus, rest]).to_dict('records'):
            L.append(f"  {int(r['n_retained_in']):6,} genes  median "
                     f"{r['median_in']:.4f} vs {r['median_other']:.4f}  fold "
                     f"{r['fold']:5.2f}  P = {r['p']:<9.3g} q = {r['q_bh']:<9.3g} "
                     f"{r['class']}")
        L.append('')
    L += ['Figure audit (paper_colors.audit, plus text past the canvas edge)',
          '  S12D: ' + ('clean' if not d_problems else '; '.join(d_problems)), '',
          'Legend notes (S12D; S12A is functional_enrichment.py\'s S6, S12B and C '
          'germline_restricted_dnds.py\'s)',
          f"  Panel: S12D {d_info['width_mm']:.0f} x {d_info['height_mm']:.1f} mm "
          f"(height from {d_info['height_source']}), bottom row beside S12C. "
          f"Panel letters are added at assembly.",
          f"  x-axis {d_info['x_range'][0]:g} to {d_info['x_range'][1]:g} (log2); "
          f"{d_info['n_interval_ends_clipped']} interval ends run off it "
          f"(arrowheads). Folds: ratio of medians, 95% bootstrap intervals "
          f"({args.n_boot:,} resamples per group); class-matched: median over "
          f"{args.n_match} matchings, interval from {args.n_match_boot:,} "
          f"resample-then-match replicates. n: eliminated / retained genes "
          f"(pairs for the matched row).",
          '  Rows as drawn (label, fold [95% CI], n):']
    for r in d_rows:
        ci = (f"[{r['lo']:.2f}, {r['hi']:.2f}]"
              if r['lo'] is not None and np.isfinite(r['lo']) else '[no interval]')
        fold = f"{r['fold']:.2f}" if np.isfinite(r['fold']) else 'NA'
        L.append(f"    {r['label'].replace(chr(10), ' '):<26} {fold:>5} {ci:<14} "
                 f"{r['n_text']}")
    L += ['', f'Outcome: {reading}', '']
    summary = '\n'.join(L)
    with open(os.path.join(args.outdir, 'dnds_within_class_summary.txt'), 'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(), 'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'scipy': scipy_version},
            'argv': sys.argv,
            'args': {k: v for k, v in vars(args).items() if k != 'unknown_re'},
            'inputs': {k: fingerprint(v) for k, v in
                       (('kaks', args.kaks), ('annotation', args.annotation),
                        ('bridge', args.bridge))},
            'join': {'mode': 'position bridge', 'bridge': bridge_info,
                     'shared_workbook_genes': {'retained': shared[0],
                                               'eliminated': shared[1]}},
            'gene_counts': {k: {'retained': r, 'eliminated': e}
                            for k, (r, e) in n_steps.items()},
            'classes_tested': len(tests),
            'bonferroni_classes': sig_classes,
            'unknown_type_classes': unknown_classes,
            'reference_class': ref_class,
            'matched': matched,
            'van_elteren': ve,
            'outcome': reading,
            'figure_S12D': {'info': d_info, 'audit': d_problems,
                            'rows': [{k: r[k] for k in ('key', 'label', 'fold',
                                                         'lo', 'hi', 'n_text')}
                                     for r in d_rows]},
        }, fh, indent=2, default=str)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
