#!/usr/bin/env python3
"""
functional_enrichment.py - gene function and expression stage, eliminated vs retained.

Figure 6A, 6B and 6C, and Supplementary Figure S12A (still written under its
earlier name, S6_dnds_vs_expression; S12B and C come from
germline_restricted_dnds.py). Replaces the May version of the same name in the
miniprot directory.

Writes:

    6AB_expression_stage.{svg,pdf,png}   6A and 6B as one 185 mm strip on shared rows
    6C_functional_class_forest.*         odds ratio per functional class, 89 mm
    S6_dnds_vs_expression.*              S12A: dN/dS against mean expression, 89 mm
    genes_with_annotation.tsv            one row per analysed gene, workbook columns joined
    expression_category_tests.tsv        germline / somatic / low: Fisher, BH, Wilson CIs
    expression_stage_tests.tsv           the ten stages: Fisher, Bonferroni, BH
    functional_class_tests.tsv           every class tested: counts, fold, odds ratio, exact CI
    enriched_set_tests.tsv               the workbook's stage-enriched gene sets
    dnds_by_stage.tsv                    n, median, quartiles of dN/dS per row of 6B
    functional_enrichment_summary.txt    counts at every filter, tests, figure audit
    provenance.json                      inputs, arguments, versions


What changed from the May script, and why
-----------------------------------------
**Denominators are the genes the workbook can classify.** The old Fisher tables
divided by every gene in the region, so a gene with no workbook match counted as
"not in this class". When the two regions match the workbook at different rates,
that difference enters every eliminated proportion relative to the retained one
(about 9% with the name join of 2026-10-09). Expression tests now use genes
with a stage of peak expression; class tests use genes with a functional class.
The summary reports the count at every step.

**Early embryo is germline.** Pre-elimination embryos still carry the germline
genome. The old category panel already counted it as germline while the stage
panel's colour key did not, so two panels of one figure disagreed. Germline is
testis, ovary, zygote and early embryo; somatic is late embryo, larvae, muscle,
intestine and carcass.

**6A and 6B share their rows.** Both panels are "something, by stage of peak
expression", so they are drawn as one strip: one row per stage, grouped under a
bold row for its category, labelled once. 6A is a paired dot plot of the share of
genes in each row, retained beside eliminated with Wilson 95% intervals, replacing
two stacked bars that needed ten colours and could only be compared at the
bottom of the stack; its P column gives the Fisher test for each row
(uncorrected; Bonferroni and BH are in the tables). 6B gives dN/dS for the same
rows, the category rows pooling their stages.

**6B no longer uses red and blue.** Germline was drawn red and somatic blue,
which in this paper means eliminated and retained. Boxes are neutral ink and the
grouping is written on the panel.

**6C is a forest plot of odds ratios with exact 95% intervals**, restricted to
classes with at least `--forest-min-elim` eliminated genes. The bar chart of fold
enrichment it replaces reached sevenfold on a single gene and showed no interval.
The interval is the conditional (exact) one that inverts Fisher's test, the same
one quoted for the top classes on 2026-08-15, so the plotted point is the
conditional maximum-likelihood odds ratio rather than the fold enrichment. For
classes this small the two agree closely; fold enrichment is in the table. The
count column is eliminated genes per class; retained counts are in the table.

**Expression is joined by position (2026-10-17).** The rate table is keyed on
miniprot model IDs; the workbook on `ascaris_gene_NNNNN`. Until 2026-10-17 the
join was the model's GFF `Target=` name (`target_bridge`), but the query
proteins come from a later annotation release than the workbook, and each
release numbers its genes serially along its own assembly, so equal numbers are
different genes: the workbook's elimination call agreed with region1 at chance
and gene lengths were uncorrelated (results draft note 60). `--bridge` reads
`workbook_model_bridge.tsv`, which `workbook_model_bridge.py` builds by lifting
each workbook gene from v37 to v50 by sequence and giving each model the gene it
overlaps; on the analysed genes the two elimination calls then agree for 95% of
eliminated and 99.9% of retained genes. `--name-join` keeps the old join only to
reproduce the superseded figures; one of the two is required.

**Figure text (hand edits of 2026-10-15, now in the code).** 6C class names in
sentence case; S6's axis "log2", not a Unicode subscript; S6's P and 6C's
Bonferroni P with any exponent raised as its own 6 pt text (`plain_text.draw_p`)
instead of mathtext.

**2026-10-19 (results draft note 67).** S12A draws the genes with zero mean
expression at log2(0.1) instead of leaving them out, so its n is 6A's less any
gene without dN/dS or stage TPMs (both counted in the summary). The 2 of its
"log2" is a subscript: its own 6 pt text, lowered 1.2 pt, drawn where the plain
label was once the layout is final (`plain_text.replace_xlabel`). 6C gains a
column of Bonferroni-adjusted P per class, right of n; the note above the panel
said only the smallest, which hid a second class below 0.05.


Tests
-----
Two-sided Fisher's exact test on each 2x2 table (eliminated vs retained, in vs
out of the category). BH within each family; Bonferroni as well, because Results
quotes Bonferroni for stages and classes. The three expression categories also
get a 2x3 chi-square test of homogeneity. dN/dS across stages: Kruskal-Wallis,
and germline vs somatic by Mann-Whitney, reported in the summary only.


Usage
-----
    python functional_enrichment.py \\
      --kaks results/annotate_regions_20260817/kaks_annotated.tsv \\
      --annotation Ascaris_genes_15714_annotation_expression.xlsx \\
      --bridge results/workbook_model_bridge_20261017_103009/workbook_model_bridge.tsv \\
      --outdir results/functional_enrichment_20261017_120000

Dependencies: numpy, pandas (+ openpyxl), scipy, matplotlib, paper_colors,
plain_text, region_lib, rate_tables, workbook_model_bridge
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import platform
import sys
import textwrap
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator
from matplotlib.transforms import ScaledTranslation
from scipy.stats import (fisher_exact, chi2_contingency, kruskal, mannwhitneyu,
                         spearmanr)
from scipy import __version__ as scipy_version

# A script that cannot import the shared style must stop, not fall back to
# defaults: a figure drawn without it would look plausible and be wrong.
from paper_colors import (REGION, NEUTRAL, FONT, LINEWIDTH, MIN_AXES_MM,
                          apply_style, audit, save, figure, mm)
from plain_text import draw_p, replace_xlabel, subscript_runs, width_pt
from region_lib import (parse_gff_models, select_primary_alignment,
                        benjamini_hochberg, fingerprint, git_commit,
                        setup_logging)
from rate_tables import read_rate_table
from workbook_model_bridge import read_bridge, shared_gene_counts

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Stage columns of the workbook, in the order the paper groups them. The
#: workbook also carries `muscle.1`, `intestine.1` and `carcass.1`, which are
#: exact duplicates of the unsuffixed columns and are ignored.
STAGES = ['testis', 'ovary', 'zygote', 'early_embryo',
          'late_embryo', 'larvae', 'muscle', 'intestine', 'carcass']

#: Pre-elimination embryos still carry the germline genome, so early embryo is
#: germline. Changing either set changes 6A, 6B and every category test.
GERMLINE = ['testis', 'ovary', 'zygote', 'early_embryo']
SOMATIC = ['late_embryo', 'larvae', 'muscle', 'intestine', 'carcass']
LOW = 'low_expression'

CATEGORY_OF = {**{s: 'germline' for s in GERMLINE},
               **{s: 'somatic' for s in SOMATIC}, LOW: 'low'}
CATEGORY_LABEL = {'germline': 'Germline', 'somatic': 'Somatic',
                  'low': 'Low expression'}
STAGE_LABEL = {'testis': 'Testis', 'ovary': 'Ovary', 'zygote': 'Zygote',
               'early_embryo': 'Early embryo', 'late_embryo': 'Late embryo',
               'larvae': 'Larvae', 'muscle': 'Muscle', 'intestine': 'Intestine',
               'carcass': 'Carcass', LOW: 'Low expression'}

PEAK_COL = 'expression'            # the workbook's stage-of-peak-expression label
CLASS_COL = "AL's functional class"
ENRICHED_COL = 'enriched_genes'

Z95 = 1.959963984540054


# ---------------------------------------------------------------------------
# Joining the rate table to the workbook
# ---------------------------------------------------------------------------

def target_bridge(gff_path):
    """THE SUPERSEDED NAME JOIN: model ID -> workbook name by query-protein name.

    The GFF `Target=` holds the query protein (`Ascaris_gene_00002_t01`);
    stripping the transcript suffix and lowercasing gives a workbook-style name
    (`ascaris_gene_00002`). This is wrong: the query proteins come from a later
    annotation release than the workbook, which numbers its genes differently,
    so equal numbers are different genes (results draft note 60). Kept for
    `--name-join`, which reproduces the superseded figures, and for
    workbook_model_bridge.py, which scores it beside the position bridge. Use
    `--bridge` for anything new.

    Returns ({model_id: name}, {model_id: raw target}, alignment diagnostics).
    """
    models = parse_gff_models(gff_path)
    best, diag = select_primary_alignment(models)
    best = best.dropna(subset=['target'])
    names = (best['target'].astype(str)
             .str.replace(r'_t\d+$', '', regex=True).str.lower())
    return (dict(zip(best['gene_id'].astype(str), names)),
            dict(zip(best['gene_id'].astype(str), best['target'].astype(str))),
            diag)


def load_workbook(path):
    """The annotation/expression workbook, with the columns this script needs."""
    wb = pd.read_excel(path)
    need = ['name', PEAK_COL] + STAGES
    missing = [c for c in need if c not in wb.columns]
    if missing:
        sys.exit(f'ERROR: workbook {path} lacks column(s) {missing}.\n'
                 f'  columns: {list(wb.columns)}')
    wb = wb.copy()
    wb['name'] = wb['name'].astype(str).str.strip().str.lower()
    for s in STAGES:
        wb[s] = pd.to_numeric(wb[s], errors='coerce')
    wb['mean_expression'] = wb[STAGES].mean(axis=1)
    if 'average' in wb.columns:
        # The workbook's own `average` is this mean rounded to 0.1. A larger
        # disagreement would mean the stage columns are not the ones it used.
        dev = (pd.to_numeric(wb['average'], errors='coerce')
               - wb['mean_expression']).abs()
        logging.info('  workbook `average` vs mean of the nine stage columns: '
                     'max |difference| %.3f', float(np.nanmax(dev.to_numpy())))
    dup = int(wb['name'].duplicated().sum())
    if dup:
        logging.warning('  %d duplicated gene names in the workbook; the first '
                        'row of each is used', dup)
        wb = wb.drop_duplicates('name', keep='first')
    logging.info('Workbook: %s genes, %s columns', f'{len(wb):,}', wb.shape[1])
    return wb


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def wilson(k, n, z=Z95):
    """Wilson score interval for k successes in n trials, as proportions."""
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def _conditional_or_fallback(a, b, c, d, level):
    """Conditional MLE odds ratio and exact CI, for scipy without
    `scipy.stats.contingency.odds_ratio` (added in 1.10).

    Conditional on all margins, the count `a` follows Fisher's noncentral
    hypergeometric distribution with the odds ratio as its parameter. The MLE
    solves E[a] = observed; the exact limits are the odds ratios at which the
    observed count sits in the upper or lower `(1 - level) / 2` tail -- the
    construction R's fisher.test uses.
    """
    from scipy.optimize import brentq
    from scipy.stats import nchypergeom_fisher
    M, n_elim, N = a + b + c + d, a + c, a + b
    lo_sup, hi_sup = max(0, N - (M - n_elim)), min(N, n_elim)
    alpha = (1 - level) / 2

    def dist(log_or):
        return nchypergeom_fisher(M, n_elim, N, np.exp(log_or))

    bound = 20.0                      # odds ratios from 2e-9 to 5e8
    # At the edge of the support the MLE and one limit are 0 or infinite.
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
    """Conditional MLE odds ratio and exact CI for [[a, b], [c, d]].

    Rows: in / not in the category. Columns: eliminated / retained. The ratio is
    the odds of being in the category for an eliminated gene over the odds for
    a retained gene, so > 1 is enrichment among eliminated genes.
    """
    try:
        from scipy.stats.contingency import odds_ratio
    except ImportError:
        mle, lo, hi = _conditional_or_fallback(a, b, c, d, level)
        return float(mle), float(lo), float(hi), 'fallback'
    res = odds_ratio(np.array([[a, b], [c, d]]), kind='conditional')
    ci = res.confidence_interval(confidence_level=level)
    return float(res.statistic), float(ci.low), float(ci.high), 'scipy'


def category_tests(df, cat_col, categories, family):
    """One Fisher test per category, eliminated vs retained. Tidy rows.

    `df` must already be restricted to genes for which `cat_col` is known;
    that restriction is what makes the denominators right.
    """
    is_elim = (df['region1'] == 'eliminated').to_numpy()
    n_e, n_r = int(is_elim.sum()), int((~is_elim).sum())
    rows = []
    for cat in categories:
        in_cat = (df[cat_col] == cat).to_numpy()
        a = int((in_cat & is_elim).sum())
        b = int((in_cat & ~is_elim).sum())
        c, d = n_e - a, n_r - b
        p_fisher = fisher_exact([[a, b], [c, d]], alternative='two-sided')[1]
        or_, or_lo, or_hi, how = conditional_or(a, b, c, d)
        pe, pr = (a / n_e if n_e else np.nan), (b / n_r if n_r else np.nan)
        e_lo, e_hi = wilson(a, n_e)
        r_lo, r_hi = wilson(b, n_r)
        rows.append({'family': family, 'category': cat,
                     'n_eliminated_in': a, 'n_retained_in': b,
                     'n_eliminated_total': n_e, 'n_retained_total': n_r,
                     'pct_eliminated': 100 * pe, 'pct_eliminated_lo': 100 * e_lo,
                     'pct_eliminated_hi': 100 * e_hi,
                     'pct_retained': 100 * pr, 'pct_retained_lo': 100 * r_lo,
                     'pct_retained_hi': 100 * r_hi,
                     'fold_enrichment': (pe / pr) if pr else np.inf,
                     'odds_ratio': or_, 'or_ci_lo': or_lo, 'or_ci_hi': or_hi,
                     'or_method': how, 'p': p_fisher})
    out = pd.DataFrame(rows)
    if len(out):
        out['p_bonferroni'] = np.minimum(out['p'] * len(out), 1.0)
        out['q_bh'] = benjamini_hochberg(out['p'].to_numpy())
    return out


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

def blank_figure(width, height_mm):
    """A journal-width canvas with no axes and no layout engine.

    Panels are placed explicitly. constrained_layout treats the axes as the
    flexible element and will shrink a panel to nothing to make room for a
    label or legend, emitting only a warning; fixed geometry fails visibly.
    """
    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax)
    return fig


def add_axes_mm(fig, left, bottom, right, top):
    """Axes box given by its margins from the canvas edges, in mm."""
    w_in, h_in = fig.get_size_inches()
    return fig.add_axes([mm(left) / w_in, mm(bottom) / h_in,
                         (w_in - mm(left) - mm(right)) / w_in,
                         (h_in - mm(bottom) - mm(top)) / h_in])


def hboxplot(ax, values, y, widths=0.62):
    """One horizontal box at height `y`, whiskers at the 5th-95th percentile.

    matplotlib 3.10 replaced `vert=False` with `orientation='horizontal'` and
    3.11 warns on the old form; the cluster's matplotlib may predate the new
    keyword, so try it first and fall back.
    """
    kw = dict(positions=[y], widths=widths, patch_artist=True,
              showfliers=False, whis=(5, 95), zorder=2)
    try:
        return ax.boxplot(values, orientation='horizontal', **kw)
    except TypeError:
        return ax.boxplot(values, vert=False, **kw)


def style_box(bp, face, box_alpha):
    """House box style: translucent face, opaque black edge, median on top.

    Alpha goes into the face colour, never `set_alpha()`, which would fade the
    outline too. Draw order: box 2, jitter 3, whiskers and caps 4, median 5.
    """
    for patch in bp['boxes']:
        patch.set_facecolor(to_rgba(face, box_alpha))
        patch.set_edgecolor('black')
        patch.set_linewidth(LINEWIDTH['thin'])
        patch.set_zorder(2)
    for part in ('whiskers', 'caps'):
        for art in bp[part]:
            art.set_color('black')
            art.set_linewidth(LINEWIDTH['thin'])
            art.set_zorder(4)
    for med in bp['medians']:
        med.set_color('black')
        med.set_linewidth(LINEWIDTH['data'])
        med.set_zorder(5)


def side_text(ax, y, s, dx_pt=3.0, x=1.0, **kw):
    """Text in a column just right of the axes, at data height `y`.

    The offset is in points, so the column sits the same distance from the
    spine whatever width the layout gives the axes. `x` is the anchor in axes
    fractions: 1.0, the right spine, for a column hung off a panel; 0.0 for a
    column in a table axes (6C).
    """
    fig = ax.figure
    tr = ax.get_yaxis_transform() + ScaledTranslation(dx_pt / 72, 0,
                                                      fig.dpi_scale_trans)
    # center_baseline, not center: a column of digits has no descenders, so
    # centring the text box would sit every number visibly above its row.
    opts = dict(ha='left', va='center_baseline', fontsize=FONT['annotation'],
                color=NEUTRAL['ink'], clip_on=False)
    opts.update(kw)
    return ax.text(x, y, s, transform=tr, **opts)


def fit_to_canvas(fig, ax, pad_mm=0.8, iterations=4, min_frac=0.5):
    """Shrink the axes box until everything attached to it lies on the canvas.

    Tick labels, the column written beside the panel, the legend above it and
    the axis labels all hang outside the axes box, and their widths depend on
    the font the node renders with (DejaVu Sans on Isaac is wider than Arial).
    Margins fixed by hand either waste the column or cut a label off; this
    measures the drawn text and moves only the edges that overflow. Text placed
    in axes coordinates moves with the box, hence the iteration. It refuses to
    shrink the box below `min_frac` of its starting size -- text that needs
    more than that has to be shortened, and audit() will say so.
    """
    w_px, h_px = fig.bbox.width, fig.bbox.height
    pad = mm(pad_mm) * fig.dpi
    start = ax.get_position()
    for _ in range(iterations):
        fig.canvas.draw()
        bb = ax.get_tightbbox(fig.canvas.get_renderer())
        pos = ax.get_position()
        dx0 = max(0.0, pad - bb.x0) / w_px
        dx1 = max(0.0, bb.x1 - (w_px - pad)) / w_px
        dy0 = max(0.0, pad - bb.y0) / h_px
        dy1 = max(0.0, bb.y1 - (h_px - pad)) / h_px
        if not (dx0 or dx1 or dy0 or dy1):
            return
        new_w = pos.width - dx0 - dx1
        new_h = pos.height - dy0 - dy1
        if new_w < min_frac * start.width or new_h < min_frac * start.height:
            logging.warning('  text around the panel needs more room than the '
                            'canvas allows; not shrinking further')
            return
        ax.set_position([pos.x0 + dx0, pos.y0 + dy0, new_w, new_h])


def drawn_texts(fig):
    """Every text object that is actually drawn.

    `fig.findobj(Text)` also returns labels of ticks the locator placed just
    outside the view limits, which are never drawn and would be reported as
    running off the canvas.
    """
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
        out += [t for t in (ax.title, getattr(ax, '_left_title', None),
                            getattr(ax, '_right_title', None)) if t is not None]
        leg = ax.get_legend()
        if leg is not None:
            out += list(leg.get_texts())
    return [t for t in out if t.get_visible() and (t.get_text() or '').strip()]


def offcanvas_text(fig, tol_px=1.0):
    """Drawn text that runs past the canvas edge, which `audit()` does not
    check. A label cut off at the edge saves without complaint and is only
    noticed once the figure is assembled."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    out = []
    for t in drawn_texts(fig):
        bb = t.get_window_extent(rend)
        if (bb.x0 < -tol_px or bb.y0 < -tol_px or bb.x1 > W + tol_px
                or bb.y1 > H + tol_px):
            out.append(f'text {t.get_text()[:28]!r} runs past the canvas edge')
    return out


def finish(fig, outdir, stem, fit=True):
    """Fit to the canvas, audit, save in every format, close.

    Returns (paths, problems); problems go into the summary as well as the log.
    """
    if fit:
        for ax in fig.axes:
            fit_to_canvas(fig, ax)
    problems = list(audit(fig, verbose=False)) + offcanvas_text(fig)
    for prob in problems:
        logging.warning('  audit %s: %s', stem, prob)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths, problems


def signed(x, nd):
    """Fixed-point number with a true minus sign (U+2212), as on the tick
    labels; an f-string writes a hyphen, which is shorter and sits lower. A
    value that rounds to zero prints unsigned."""
    s = f'{x:.{nd}f}'
    if float(s) == 0:
        s = s.lstrip('-')
    return s.replace('-', '\u2212')


def short_p(p):
    """P for a narrow panel column: two significant figures, trailing zeros
    kept (0.20, not 0.2), floored at 0.001. Exact values are in the tables."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p < 0.001:
        return '< 0.001'
    if p >= 0.995:
        return '1.0'
    decimals = max(0, 1 - int(np.floor(np.log10(p))))
    out = f'{p:.{decimals}f}'
    # Rounding can carry into a new digit (0.0996 -> 0.100); trim back to two.
    if float(out) >= 10 ** (np.floor(np.log10(p)) + 1):
        out = f'{float(out):.{max(0, decimals - 1)}f}'
    return out


def region_handles(n_ret, n_elim, marker='o', ms=3.6):
    return ([Line2D([], [], ls='none', marker=marker, ms=ms,
                    mfc=REGION['retained'], mec=REGION['retained']),
             Line2D([], [], ls='none', marker=marker, ms=ms,
                    mfc=REGION['eliminated'], mec=REGION['eliminated'])],
            [f'Retained (n = {n_ret:,})', f'Eliminated (n = {n_elim:,})'])


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def stage_rows():
    """Row structure shared by 6A and 6B, top to bottom: (kind, key, label).

    Low expression is both a category and its only stage, so it is one row.
    """
    rows = [('category', 'germline', CATEGORY_LABEL['germline'])]
    rows += [('stage', s, STAGE_LABEL[s]) for s in GERMLINE]
    rows += [('category', 'somatic', CATEGORY_LABEL['somatic'])]
    rows += [('stage', s, STAGE_LABEL[s]) for s in SOMATIC]
    rows += [('category', 'low', CATEGORY_LABEL['low'])]
    return rows


def row_positions(rows, gap=0.45):
    """y for each row, top row at 0, with a gap above every category row."""
    ys, y = [], 0.0
    for i, (kind, _, _) in enumerate(rows):
        if kind == 'category' and i:
            y -= gap
        ys.append(y)
        y -= 1.0
    return np.array(ys)


def place_strip(fig, ax_a, ax_b, share_a=0.45, gap_mm=5.0, pad_mm=0.8,
                iterations=4):
    """Lay out two panels side by side on one canvas from their measured text.

    Left of A: the row labels. Right of A: its P column. Right of B: its n
    column. Above both: the legend and the reference label. Everything is
    measured from the drawn text and the remaining width split `share_a` to A,
    so the strip fills 185 mm exactly whatever font the node has.
    """
    W, H = fig.bbox.width, fig.bbox.height
    pad, gap = mm(pad_mm) * fig.dpi, mm(gap_mm) * fig.dpi
    for _ in range(iterations):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        ta, tb = ax_a.get_tightbbox(r), ax_b.get_tightbbox(r)
        pa, pb = ax_a.get_window_extent(r), ax_b.get_window_extent(r)
        left = max(0.0, pa.x0 - ta.x0)
        a_right = max(0.0, ta.x1 - pa.x1)
        b_left = max(0.0, pb.x0 - tb.x0)
        b_right = max(0.0, tb.x1 - pb.x1)
        top = max(ta.y1 - pa.y1, tb.y1 - pb.y1, 0.0)
        bottom = max(pa.y0 - ta.y0, pb.y0 - tb.y0, 0.0)
        avail = W - 2 * pad - left - a_right - gap - b_left - b_right
        wa = avail * share_a
        wb = avail - wa
        x_a = pad + left
        x_b = x_a + wa + a_right + gap + b_left
        y0 = pad + bottom
        hh = H - 2 * pad - bottom - top
        ax_a.set_position([x_a / W, y0 / H, wa / W, hh / H])
        ax_b.set_position([x_b / W, y0 / H, wb / W, hh / H])


def place_with_table(fig, ax, ax_t, gap_pt=3.0, pad_mm=0.8, min_table_mm=None,
                     iterations=4):
    """Lay out a panel and the table of columns to its right, from drawn text.

    6C's n and P columns are a table axes of their own (axis off, the panel's
    rows), not text hung off the panel's edge: a column of values aligned with
    the rows is part of the panel, and audit() measures the span of a figure's
    axes to catch a legend taking the page. As text on the panel's edge, the
    two columns left the plot 46% of the 89 mm width and failed that check.
    The table is as wide as its widest text, but at least `min_table_mm`
    (default: audit()'s floor for an axes, plus 0.5 mm), and starts `gap_pt` right
    of the panel; the panel takes what is left of the canvas after its own
    measured text (row labels, tick labels, axis label).
    """
    if min_table_mm is None:
        min_table_mm = MIN_AXES_MM + 0.5
    W, H = fig.bbox.width, fig.bbox.height
    pad = mm(pad_mm) * fig.dpi
    gap = gap_pt * fig.dpi / 72.0
    t_min = mm(min_table_mm) * fig.dpi
    for _ in range(iterations):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        pp, tp = ax.get_window_extent(r), ax.get_tightbbox(r)
        tb = ax_t.get_window_extent(r)
        ext = [t.get_window_extent(r) for t in ax_t.texts
               if t.get_visible() and (t.get_text() or '').strip()]
        t_w = (max(b.x1 for b in ext) - tb.x0) if ext else 0.0
        t_top = max((b.y1 for b in ext), default=pp.y1)
        left = max(0.0, pp.x0 - tp.x0)
        bottom = max(0.0, pp.y0 - tp.y0)
        top = max(0.0, tp.y1 - pp.y1, t_top - pp.y1)
        right = max(gap + max(t_w, t_min), tp.x1 - pp.x1)
        w, h = W - 2 * pad - left - right, H - 2 * pad - bottom - top
        x0, y0 = pad + left, pad + bottom
        ax.set_position([x0 / W, y0 / H, w / W, h / H])
        ax_t.set_position([(x0 + w + gap) / W, y0 / H, (right - gap) / W, h / H])


def plot_stage_strip(cat_tests, stage_tests, has_peak, outdir, stem, rng,
                     metric_label='dN/dS', cap_pct=99.0, jitter_max=600,
                     box_alpha=0.40, height_mm=66.0):
    """6A + 6B on shared rows, 185 mm.

    6A: share of genes with peak expression in each row, retained just above
    eliminated so the two Wilson intervals never overlap. P column: Fisher's
    exact test for the row, uncorrected (BH for categories and Bonferroni for
    stages are in the tables and summary).

    6B: dN/dS for the genes in each row, both regions pooled. Category rows pool
    their stages and are drawn as boxes only, since their genes are the stage
    rows' genes; stage rows carry jitter on top of a translucent box. Boxes are
    quartiles, whiskers the 5th-95th percentile, no fliers. The x-axis stops at
    the `cap_pct` percentile of all values, which bounds only what is drawn;
    jittered points beyond it sit on the edge and are counted.
    """
    rows = stage_rows()
    ys = row_positions(rows)
    tests = {('category', r.category): r for r in cat_tests.itertuples()}
    tests.update({('stage', r.category): r for r in stage_tests.itertuples()})

    fig = blank_figure('double', height_mm)
    ax_a = add_axes_mm(fig, left=24, bottom=11, right=105, top=8)
    ax_b = add_axes_mm(fig, left=100, bottom=11, right=12, top=8)

    # --- 6A ---------------------------------------------------------------
    dodge, xmax = 0.17, 0.0
    for (kind, key, _), yy in zip(rows, ys):
        t = tests.get((kind, key)) or tests.get(('stage', key))
        if t is None:
            continue
        for reg, off in (('retained', dodge), ('eliminated', -dodge)):
            v = getattr(t, f'pct_{reg}')
            lo, hi = getattr(t, f'pct_{reg}_lo'), getattr(t, f'pct_{reg}_hi')
            ax_a.plot([lo, hi], [yy + off] * 2, color=REGION[reg],
                      lw=LINEWIDTH['thin'], solid_capstyle='butt', zorder=2)
            ax_a.plot([v], [yy + off], ls='none', marker='o', ms=3.4,
                      mfc=REGION[reg], mec=REGION[reg], zorder=3)
            xmax = max(xmax, hi)
        side_text(ax_a, yy, short_p(t.p),
                  color=NEUTRAL['ink'] if kind == 'category' else NEUTRAL['ink_light'],
                  fontweight='bold' if kind == 'category' else 'normal')
    side_text(ax_a, ys[0] + 0.62, 'P', va='bottom')
    ax_a.set_xlim(0, xmax * 1.05)
    ax_a.set_xlabel('Genes by stage of peak expression (%)')
    n_ret = int(cat_tests['n_retained_total'].iloc[0])
    n_elim = int(cat_tests['n_eliminated_total'].iloc[0])
    h, l = region_handles(n_ret, n_elim, ms=3.4)
    ax_a.legend(h, l, loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2,
                frameon=False, fontsize=FONT['legend'], borderaxespad=0.4,
                handletextpad=0.3, columnspacing=1.2)

    # --- 6B ---------------------------------------------------------------
    vals_all = has_peak['dN_dS'].dropna().to_numpy()
    xcap = float(np.percentile(vals_all, cap_pct))
    pooled_median = float(np.median(vals_all))
    rows_out = []
    for (kind, key, _), yy in zip(rows, ys):
        if kind == 'category':
            sel = has_peak['category'] == key
        else:
            sel = has_peak[PEAK_COL] == key
        v = has_peak.loc[sel, 'dN_dS'].dropna().to_numpy()
        n_drawn = 0
        if len(v) >= 5:
            bp = hboxplot(ax_b, v, yy, widths=0.62)
            style_box(bp, NEUTRAL['ink_light'], box_alpha)
            if kind == 'stage':
                show = v if len(v) <= jitter_max else rng.choice(
                    v, jitter_max, replace=False)
                n_drawn = len(show)
                ax_b.scatter(np.minimum(show, xcap * 1.02),
                             yy + rng.uniform(-0.22, 0.22, len(show)),
                             s=1.1, linewidths=0, color=NEUTRAL['ink_light'],
                             alpha=0.45, zorder=3, rasterized=True)
        side_text(ax_b, yy, f'{len(v):,}',
                  fontweight='bold' if kind == 'category' else 'normal')
        rows_out.append({'row': kind, 'key': key, 'n': len(v),
                         'median': float(np.median(v)) if len(v) else np.nan,
                         'q1': float(np.percentile(v, 25)) if len(v) else np.nan,
                         'q3': float(np.percentile(v, 75)) if len(v) else np.nan,
                         'n_beyond_axis': int((v > xcap).sum()),
                         'n_jitter_drawn': n_drawn})
    side_text(ax_b, ys[0] + 0.62, 'n', va='bottom')
    ax_b.axvline(pooled_median, color=NEUTRAL['rule'], lw=LINEWIDTH['hairline'],
                 ls=(0, (2, 2)), zorder=1)
    ax_b.text(pooled_median, 1.0, f'all genes, median {pooled_median:.2f}',
              transform=ax_b.get_xaxis_transform(), ha='left', va='bottom',
              fontsize=FONT['annotation'], color=NEUTRAL['ink'])
    ax_b.set_xlim(0, xcap * 1.04)
    ax_b.set_xlabel(metric_label)

    # --- shared rows ------------------------------------------------------
    for ax in (ax_a, ax_b):
        ax.set_yticks(ys)
        ax.set_ylim(ys[-1] - 0.7, ys[0] + 0.7)
        ax.tick_params(axis='y', length=0)
        ax.spines['left'].set_visible(False)
    ax_a.set_yticklabels([lab for _, _, lab in rows])
    for t, (kind, _, _) in zip(ax_a.get_yticklabels(), rows):
        t.set_fontweight('bold' if kind == 'category' else 'normal')
    ax_b.set_yticklabels([])

    place_strip(fig, ax_a, ax_b)
    paths, problems = finish(fig, outdir, stem, fit=False)
    stage_tab = pd.DataFrame(rows_out)
    jitter_dropped = int((stage_tab.loc[stage_tab['row'] == 'stage', 'n']
                          - stage_tab.loc[stage_tab['row'] == 'stage',
                                          'n_jitter_drawn']).sum())
    return paths, problems, stage_tab, {
        'x_axis_cap': xcap, 'cap_pct': cap_pct,
        'n_beyond_axis': int(stage_tab.loc[stage_tab['row'] == 'stage',
                                           'n_beyond_axis'].sum()),
        'n_jitter_not_drawn': jitter_dropped, 'jitter_max': jitter_max,
        'pooled_median': pooled_median}


def sentence_case(name):
    """A workbook class name in sentence case, for a figure label.

    The workbook writes classes in title case ("Signal Transduction",
    "Secreted / Unknown"). Only the first letter of the label stays capital; a
    word keeps its own case when any hyphen-separated part of it is an acronym
    or a single capital ("GPI", "C-type"), carries a digit, or has a capital
    after its first letter ("membrAnchor"). Tables keep the workbook's names.
    """
    def keep(word):
        core = word.strip('()[],;:')
        parts = [p for p in core.split('-') if p]
        return any((p.isupper() and p.isalpha()) or any(c.isdigit() for c in p)
                   or any(c.isupper() for c in p[1:]) for p in parts)

    words = str(name).split(' ')
    out = [w if keep(w) else w.lower() for w in words]
    for i, w in enumerate(out):
        if w and w[0].isalpha():
            if not keep(words[i]):
                out[i] = w[0].upper() + w[1:]
            break
    return ' '.join(out)


def wrap_class(name, width=21):
    """Wrap a class name at its ' / ' separators, then by words if needed."""
    parts = [p.strip() for p in str(name).split('/')]
    lines, cur = [], ''
    for p in parts:
        cand = f'{cur} / {p}' if cur else p
        if len(cand) <= width or not cur:
            cur = cand
        else:
            lines.append(cur + ' /')
            cur = p
    lines.append(cur)
    out = []
    for ln in lines:
        out += textwrap.wrap(ln, width) or ['']
    return '\n'.join(out)


def plot_class_forest(tests, outdir, stem, min_elim, n_tested):
    """6C - odds ratio and exact 95% CI for each class with >= min_elim
    eliminated genes, on a log axis with a rule at 1.

    Markers are neutral ink: an odds ratio is not a region, and drawing
    significant rows red would read as "eliminated" to anyone who learned the
    convention from Figure 2. Rows are ordered by odds ratio. Right of the
    panel: the number of eliminated genes in the class, then its
    Bonferroni-adjusted P over all `n_tested` classes. Until 2026-10-19 the
    panel named only the smallest Bonferroni P, so a second class that also
    passed (Proteasome, 0.017, beside Unknown, 0.002) was invisible (results
    draft note 67).
    """
    sub = tests[tests['n_eliminated_in'] >= min_elim].copy()
    if sub.empty:
        logging.warning('  6C: no class has >= %d eliminated genes', min_elim)
        return None, [], sub
    sub = sub.sort_values('odds_ratio', ascending=False).reset_index(drop=True)
    labels = [wrap_class(sentence_case(c)) for c in sub['category']]
    n_lines = [lab.count('\n') + 1 for lab in labels]
    # Row pitch follows the number of label lines so wrapped names never touch.
    pitch = [max(1.0, 0.62 * k + 0.38) for k in n_lines]
    ys, y = [], 0.0
    for i, p in enumerate(pitch):
        if i:
            y -= (pitch[i - 1] + p) / 2
        ys.append(y)
    ys = np.array(ys)
    span_rows = ys[0] - ys[-1] + (pitch[0] + pitch[-1]) / 2 + 0.5
    height_mm = float(min(200.0, 22.0 + 4.2 * span_rows))

    fig = blank_figure('single', height_mm)
    # Starting boxes only: place_with_table sets both from the drawn text.
    ax = add_axes_mm(fig, left=30.0, bottom=10.5, right=20.0, top=9.0)
    ax_t = add_axes_mm(fig, left=70.0, bottom=10.5, right=1.0, top=9.0)
    ax_t.axis('off')
    # The P column starts 6 pt past the widest n, as Arial sets them.
    fs = FONT['annotation']
    n_txt = [f"{int(v)}" for v in sub['n_eliminated_in']]
    p_txt = [short_p(v) for v in sub['p_bonferroni']]
    w_n = max(width_pt(fig, t, fs) for t in n_txt + ['n'])
    w_p = max(width_pt(fig, t, fs) for t in p_txt + ['P'])
    p_dx = w_n + 6.0

    finite_lo = sub['or_ci_lo'][sub['or_ci_lo'] > 0]
    finite_hi = sub['or_ci_hi'][np.isfinite(sub['or_ci_hi'])]
    lo = min(finite_lo.min() if len(finite_lo) else 0.5, 0.8)
    hi = max(finite_hi.max() if len(finite_hi) else 2.0, 1.25)
    e_lo = max(-6, int(np.floor(np.log2(lo))))
    e_hi = min(6, int(np.ceil(np.log2(hi))))
    x_lo, x_hi = 2.0 ** e_lo, 2.0 ** e_hi
    n_clipped = 0
    for i, ((_, r), yy) in enumerate(zip(sub.iterrows(), ys)):
        a, b = r['or_ci_lo'], r['or_ci_hi']
        a_c = max(a, x_lo)
        b_c = min(b, x_hi) if np.isfinite(b) else x_hi
        ax.plot([a_c, b_c], [yy, yy], color=NEUTRAL['ink'],
                lw=LINEWIDTH['thin'], solid_capstyle='butt', zorder=2)
        # An interval running off the axis gets an arrowhead, so a clipped
        # bound is never mistaken for a bound that ends there.
        if a < x_lo:
            ax.plot([x_lo], [yy], marker='<', ms=3, color=NEUTRAL['ink'],
                    zorder=3, clip_on=False)
            n_clipped += 1
        if not np.isfinite(b) or b > x_hi:
            ax.plot([x_hi], [yy], marker='>', ms=3, color=NEUTRAL['ink'],
                    zorder=3, clip_on=False)
            n_clipped += 1
        est = min(max(r['odds_ratio'], x_lo), x_hi)
        ax.plot([est], [yy], ls='none', marker='s', ms=3.2,
                mfc=NEUTRAL['ink'], mec=NEUTRAL['ink'], zorder=4)
        side_text(ax_t, yy, n_txt[i], dx_pt=0.0, x=0.0)
        side_text(ax_t, yy, p_txt[i], dx_pt=p_dx, x=0.0)
    head_y = ys[0] + pitch[0] / 2 + 0.2
    side_text(ax_t, head_y, 'n', dx_pt=0.0, x=0.0, va='bottom')
    side_text(ax_t, head_y, 'P', dx_pt=p_dx, x=0.0, va='bottom')
    ax.axvline(1.0, color=NEUTRAL['rule'], lw=LINEWIDTH['rule'],
               ls=(0, (4, 2)), zorder=1)

    ax.set_xscale('log', base=2)
    ax.set_xlim(x_lo, x_hi)
    step = max(1, int(np.ceil((e_hi - e_lo) / 6)))
    ticks = [2.0 ** e for e in range(e_lo, e_hi + 1, step)]
    if 1.0 not in ticks:
        ticks = sorted(set(ticks) | {1.0})
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f'{t:g}' for t in ticks])
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, linespacing=1.0)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(ys[-1] - pitch[-1] / 2 - 0.25, ys[0] + pitch[0] / 2 + 0.25)
    ax_t.set_ylim(ax.get_ylim())
    ax.set_xlabel('Odds ratio, eliminated vs retained (95% CI)')
    # What the P column is, on a line of its own above the column heads and
    # right-aligned with the P column (on the heads' line it ran into "n").
    # The column is in short_p form (< 0.001 at the floor); exact values are
    # in functional_class_tests.tsv.
    ax_t.text(0.0, 1.0, f'P, Bonferroni-adjusted for {n_tested} classes',
              transform=ax_t.transAxes + ScaledTranslation(
                  (p_dx + w_p) / 72, 10.0 / 72, fig.dpi_scale_trans),
              ha='right', va='baseline', fontsize=fs, color=NEUTRAL['ink'])
    place_with_table(fig, ax, ax_t)
    paths, problems = finish(fig, outdir, stem, fit=False)
    sub['n_ci_ends_clipped'] = n_clipped
    return paths, problems, sub


def plot_dnds_vs_expression(df, outdir, stem, rng, metric_label='dN/dS',
                            cap_pct=99.5, height_mm=72.0):
    """S12A (file S6_dnds_vs_expression) - per-gene dN/dS against mean
    expression over the nine stages.

    Spearman on every plotted gene. A gene with zero mean expression is drawn
    at log2(0.1), the left edge: the 0.1 TPM pseudocount is there to give it a
    value. Until 2026-10-19 these genes were left out ("no log value"), so S12
    had fewer genes than 6A (results draft note 67). A gene
    without dN/dS or without stage TPMs cannot be drawn; both counts are
    returned. Points above the `cap_pct` percentile are drawn as open triangles
    on the top edge and counted on the panel, so a clipped point is never
    mistaken for an absent one.
    """
    n_no_dnds = int(df['dN_dS'].isna().sum())
    n_no_expr = int((df['dN_dS'].notna() & df['mean_expression'].isna()).sum())
    sub = df.dropna(subset=['dN_dS', 'mean_expression']).copy()
    n_zero = int((sub['mean_expression'] <= 0).sum())
    n_neg = int((sub['mean_expression'] < 0).sum())
    if n_neg:
        logging.warning('  S12A: %d genes with a negative mean TPM, drawn as zero',
                        n_neg)
    sub['log2_expr'] = np.log2(sub['mean_expression'].clip(lower=0.0) + 0.1)
    rho, p = spearmanr(sub['log2_expr'], sub['dN_dS'])
    ycap = float(np.percentile(sub['dN_dS'], cap_pct))
    # A margin at each end of the y-axis. With the axis at exactly (0, cap), a
    # gene at dN/dS = 0 sat on the x-axis and was cut in half by it, and a gene
    # just under the cap was cut by the top of the axes (both in the 2026-10-09
    # build). The triangles for genes above the cap go on the new top edge, so
    # they no longer share a row with genes inside the range.
    y_lo, y_hi = -0.02 * ycap, 1.03 * ycap

    fig = blank_figure('single', height_mm)
    ax = add_axes_mm(fig, left=13.0, bottom=10.5, right=3.0, top=12.0)
    n_above = {}
    for reg in ('retained', 'eliminated'):        # eliminated drawn on top
        r = sub[sub['region1'] == reg]
        inside = r[r['dN_dS'] <= ycap]
        above = r[r['dN_dS'] > ycap]
        n_above[reg] = len(above)
        ax.scatter(inside['log2_expr'], inside['dN_dS'],
                   s=2.0 if reg == 'retained' else 3.0, linewidths=0,
                   color=to_rgba(REGION[reg], 0.35 if reg == 'retained' else 0.8),
                   zorder=2 if reg == 'retained' else 3, rasterized=True)
        if len(above):
            ax.scatter(above['log2_expr'], np.full(len(above), y_hi),
                       marker='^', s=6, facecolors='none',
                       edgecolors=REGION[reg], linewidths=0.4, zorder=4,
                       clip_on=False)
    ax.set_ylim(y_lo, y_hi)
    # The plain label reserves the room; once the layout is final the 2 of
    # log2 is redrawn as a subscript, its own 6 pt text lowered 1.2 pt
    # (plain_text.replace_xlabel). Not a Unicode subscript: Arial has none, so
    # Illustrator set "log\u2082" in Myriad Pro (hand-fixed 2026-10-15).
    ax.set_xlabel('Mean expression, log2(TPM + 0.1)')
    ax.set_ylabel(metric_label)
    n_r = int((sub['region1'] == 'retained').sum())
    n_e = int((sub['region1'] == 'eliminated').sum())
    tot_above = sum(n_above.values())
    # Key above the axes in two rows: the two regions, then the open-triangle
    # marker for points above the axis with their count -- a blue and a red
    # triangle side by side, as drawn (one black triangle matched neither).
    # The Spearman result sits right-aligned on the second row.
    h, l = region_handles(n_r, n_e, ms=3.2)
    tri = tuple(Line2D([], [], ls='none', marker='^', ms=3.2, mfc='none',
                       mec=REGION[reg], mew=0.5)
                for reg in ('retained', 'eliminated'))
    blank = Line2D([], [], ls='none')
    ax.legend([h[0], tri, h[1], blank],
              [l[0], f'Above the axis ({tot_above})', l[1], ''],
              loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2,
              frameon=False, fontsize=FONT['legend'], borderaxespad=0.3,
              handletextpad=0.3, columnspacing=1.2, labelspacing=0.3,
              handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
    # Plain text with any exponent raised as its own 6 pt text; fmt_p wrote
    # the exponent in mathtext, which Illustrator opens in another typeface.
    draw_p(ax, 1.0, 1.0, p, prefix=f'Spearman ρ = {signed(rho, 3)}, ',
           transform=ax.transAxes, ha='right', va='bottom', dy_pt=3.0,
           fontsize=FONT['annotation'], color=NEUTRAL['ink'])
    fit_to_canvas(fig, ax)
    replace_xlabel(ax, subscript_runs('Mean expression, log', '2', '(TPM + 0.1)',
                                      fontsize=ax.xaxis.label.get_fontsize()))
    paths, problems = finish(fig, outdir, stem, fit=False)
    return paths, problems, {'rho': float(rho), 'p': float(p), 'n': len(sub),
                             'n_retained': n_r, 'n_eliminated': n_e,
                             'n_zero_expression_drawn': n_zero,
                             'zero_expression_x': float(np.log2(0.1)),
                             'n_no_dnds': n_no_dnds,
                             'n_no_expression': n_no_expr,
                             'n_input': len(df),
                             'y_cap': ycap, 'cap_pct': cap_pct,
                             'n_above_cap': n_above}


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Functional class and expression stage of eliminated vs '
                    'retained genes; Figure 6A-C and S6.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--kaks', required=True,
                   help='Annotated rate table (ANNOTATED_TABLE, from '
                        'annotate_regions.py)')
    p.add_argument('--annotation', required=True,
                   help='Ascaris_genes_15714_annotation_expression.xlsx')
    p.add_argument('--bridge', default=None,
                   help='workbook_model_bridge.tsv (workbook_model_bridge.py): '
                        'model ID -> workbook gene, by position. The join to use')
    p.add_argument('--name-join', action='store_true',
                   help='Join by the GFF Target= name instead, as before '
                        '2026-10-17. Wrong (results draft note 60); only to '
                        'reproduce the superseded figures. Needs --gff1')
    p.add_argument('--gff1', default=None,
                   help='Ascaris miniprot GFF; used only by --name-join')
    p.add_argument('--gene-col', default='gene1')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--min-class-genes', type=int, default=5,
                   help='A class is tested if it holds at least this many '
                        'analysed genes, both regions together (default 5, as '
                        'before)')
    p.add_argument('--forest-min-elim', type=int, default=10,
                   help='6C draws classes with at least this many eliminated '
                        'genes (default 10)')
    p.add_argument('--dnds-cap-pct', type=float, default=99.0,
                   help='6B x-axis stops at this percentile; drawing only')
    p.add_argument('--jitter-max', type=int, default=600,
                   help='6B draws at most this many jittered points per stage; '
                        'boxes and statistics use every gene (default 600)')
    p.add_argument('--s6-cap-pct', type=float, default=99.5,
                   help='S6 y-axis stops at this percentile; points above are '
                        'drawn on the edge and counted')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if bool(args.bridge) == bool(args.name_join):
        sys.exit('ERROR: give --bridge workbook_model_bridge.tsv (the position '
                 'join), or --name-join to reproduce the superseded name join; '
                 'exactly one of the two.')
    if args.name_join and not args.gff1:
        sys.exit('ERROR: --name-join needs --gff1.')
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'functional_enrichment_'
            + datetime.now().strftime('%Y%m%d_%H%M%S'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'functional_enrichment')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('functional_enrichment.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()
    rng = np.random.default_rng(args.seed)

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
    if df.empty:
        sys.exit('ERROR: no row is labelled retained or eliminated.')

    # --- join -----------------------------------------------------------------
    align_diag, bridge_info = None, None
    if args.bridge:
        name_of, bridge_info = read_bridge(args.bridge)
        df['annot_name'] = df['gene1'].map(name_of)
        join_step, join_col = 'given a workbook gene (bridge)', 'annot_name'
    else:
        logging.warning('--name-join: the superseded Target-name join, wrong '
                        'for most genes (results draft note 60)')
        name_of, target_of, align_diag = target_bridge(args.gff1)
        df['target'] = df['gene1'].map(target_of)
        df['annot_name'] = df['gene1'].map(name_of)
        join_step, join_col = 'with a GFF Target', 'target'
    wb = load_workbook(args.annotation)
    keep = ['name', PEAK_COL, 'mean_expression'] + STAGES + \
        [c for c in (CLASS_COL, ENRICHED_COL) if c in wb.columns]
    df = df.merge(wb[keep], left_on='annot_name', right_on='name', how='left')
    if CLASS_COL not in df.columns:
        logging.warning("Workbook has no %r column; 6C is skipped", CLASS_COL)

    def counts(mask):
        return (int((mask & (df['region1'] == 'retained')).sum()),
                int((mask & (df['region1'] == 'eliminated')).sum()))

    n_steps = {
        'analysed (retained or eliminated)': counts(df['gene1'].notna()),
        join_step: counts(df[join_col].notna()),
        'matched to the workbook': counts(df['name'].notna()),
        'with a stage of peak expression': counts(df[PEAK_COL].notna()),
    }
    if CLASS_COL in df.columns:
        n_steps['with a functional class'] = counts(df[CLASS_COL].notna())
    for k, (r, e) in n_steps.items():
        logging.info('  %-34s retained %s, eliminated %s', k, f'{r:,}', f'{e:,}')
    if n_steps['matched to the workbook'] == (0, 0):
        sys.exit('ERROR: no rate-table gene matched a workbook name.\n'
                 f"  bridged names: {df['annot_name'].dropna().head(3).tolist()}\n"
                 f"  workbook:      {wb['name'].head(3).tolist()}")
    # Two analysed models can take one workbook gene (by position: a gene the
    # later release splits in two; by name: two models of one query protein).
    # Each keeps the shared gene's values; the count goes in the summary.
    shared = shared_gene_counts(df['annot_name'], df['region1'])
    if sum(shared):
        logging.info('  %s analysed genes share their workbook gene with another '
                     'analysed gene (retained %s, eliminated %s)',
                     f'{sum(shared):,}', f'{shared[0]:,}', f'{shared[1]:,}')

    df['category'] = df[PEAK_COL].map(CATEGORY_OF)
    other = df[PEAK_COL].notna() & df['category'].isna()
    n_other = int(other.sum())
    if n_other:
        logging.warning('  %d genes have an unrecognised peak-stage label and are '
                        'left out of the expression tests: %s', n_other,
                        sorted(df.loc[other, PEAK_COL].astype(str).unique())[:5])

    out_cols = ['gene1'] + [c for c in ('gene2',) if c in df.columns] + \
        ['region1'] + [c for c in ('dN', 'dS', 'dN_dS') if c in df.columns] + \
        [c for c in ('target',) if c in df.columns] + \
        ['annot_name', PEAK_COL, 'category'] + \
        [c for c in (ENRICHED_COL, CLASS_COL) if c in df.columns] + \
        STAGES + ['mean_expression']
    df[out_cols].to_csv(os.path.join(args.outdir, 'genes_with_annotation.tsv'),
                        sep='\t', index=False, float_format='%.6g')

    # --- expression tests -----------------------------------------------------
    has_peak = df[df['category'].notna()].copy()
    cat_tests = category_tests(has_peak, 'category',
                               ['germline', 'somatic', 'low'], 'category')
    stage_tests = category_tests(has_peak, PEAK_COL, STAGES + [LOW], 'stage')
    tab = np.array([[int(((has_peak['region1'] == reg) &
                          (has_peak['category'] == c)).sum())
                     for c in ('germline', 'somatic', 'low')]
                    for reg in ('eliminated', 'retained')])
    chi2, chi_p, chi_df, _ = chi2_contingency(tab, correction=False)
    cat_tests['chi2_homogeneity'] = chi2
    cat_tests['chi2_df'] = chi_df
    cat_tests['chi2_p'] = chi_p
    cat_tests.to_csv(os.path.join(args.outdir, 'expression_category_tests.tsv'),
                     sep='\t', index=False, float_format='%.6g')
    stage_tests.to_csv(os.path.join(args.outdir, 'expression_stage_tests.tsv'),
                       sep='\t', index=False, float_format='%.6g')

    if ENRICHED_COL in has_peak.columns:
        with_set = has_peak[has_peak[ENRICHED_COL].notna()]
        sets = sorted(with_set[ENRICHED_COL].astype(str).unique())
        category_tests(with_set, ENRICHED_COL, sets, 'enriched_set').to_csv(
            os.path.join(args.outdir, 'enriched_set_tests.tsv'), sep='\t',
            index=False, float_format='%.6g')

    # --- functional classes ---------------------------------------------------
    class_tests = pd.DataFrame()
    if CLASS_COL in df.columns:
        has_class = df[df[CLASS_COL].notna()].copy()
        has_class[CLASS_COL] = has_class[CLASS_COL].astype(str).str.strip()
        sizes = has_class[CLASS_COL].value_counts()
        tested = sorted(sizes[sizes >= args.min_class_genes].index)
        class_tests = category_tests(has_class, CLASS_COL, tested, 'class')
        class_tests = class_tests.sort_values('p').reset_index(drop=True)
        class_tests.to_csv(os.path.join(args.outdir, 'functional_class_tests.tsv'),
                           sep='\t', index=False, float_format='%.6g')

    # --- figures --------------------------------------------------------------
    audits = {}
    logging.info('6A + 6B: expression stage strip')
    _, audits['6AB'], by_stage, b_info = plot_stage_strip(
        cat_tests, stage_tests, has_peak, args.outdir, '6AB_expression_stage',
        rng, cap_pct=args.dnds_cap_pct, jitter_max=args.jitter_max)
    by_stage.to_csv(os.path.join(args.outdir, 'dnds_by_stage.tsv'), sep='\t',
                    index=False, float_format='%.6g')
    stage_groups = [has_peak.loc[has_peak[PEAK_COL] == s, 'dN_dS'].dropna()
                    for s in STAGES + [LOW]]
    stage_groups = [g for g in stage_groups if len(g) >= 5]
    kw_h, kw_p = (kruskal(*stage_groups) if len(stage_groups) > 1
                  else (np.nan, np.nan))
    g_v = has_peak.loc[has_peak['category'] == 'germline', 'dN_dS'].dropna()
    s_v = has_peak.loc[has_peak['category'] == 'somatic', 'dN_dS'].dropna()
    gs_p = (mannwhitneyu(g_v, s_v, alternative='two-sided')[1]
            if len(g_v) >= 5 and len(s_v) >= 5 else np.nan)

    forest = pd.DataFrame()
    if len(class_tests):
        logging.info('6C: functional class forest')
        _, audits['6C'], forest = plot_class_forest(
            class_tests, args.outdir, '6C_functional_class_forest',
            args.forest_min_elim, len(class_tests))

    logging.info('S12A (file S6_dnds_vs_expression): dN/dS vs expression level')
    _, audits['S12A'], s6 = plot_dnds_vs_expression(
        df[df['name'].notna()], args.outdir, 'S6_dnds_vs_expression', rng,
        cap_pct=args.s6_cap_pct)

    # --- summary --------------------------------------------------------------
    if args.bridge:
        join_line = (f"Join         position bridge {os.path.basename(args.bridge)}"
                     f": {bridge_info['assigned']:,} of {bridge_info['models']:,} "
                     f"models given a workbook gene")
    else:
        join_line = (f"Join         SUPERSEDED name join (results draft note 60): "
                     f"GFF Target= of {os.path.basename(args.gff1)}, "
                     f"{align_diag['n_ids']:,} models")
    L = ['Functional class and expression stage, eliminated vs retained', '',
         f'Rate table   {os.path.basename(args.kaks)}, {n_rows:,} rows',
         f'Workbook     {os.path.basename(args.annotation)}, {len(wb):,} genes',
         join_line, '',
         'Genes at each step                    retained  eliminated  (% elim.)']
    for k, (r, e) in n_steps.items():
        L.append(f'  {k:<34} {r:>8,}  {e:>10,}  ({100 * e / max(r + e, 1):.2f}%)')
    an_r, an_e = n_steps['analysed (retained or eliminated)']
    m_r, m_e = n_steps['matched to the workbook']
    L += ['', f'  Workbook match rate: retained {100 * m_r / max(an_r, 1):.1f}%, '
              f'eliminated {100 * m_e / max(an_e, 1):.1f}%. Every test uses '
              f'matched genes only.',
          f'  Sharing their workbook gene with another analysed gene: retained '
          f'{shared[0]:,}, eliminated {shared[1]:,} (each keeps its values).',
          '  (The May script divided by every gene in the region, which let the '
          'difference in match',
          '  rate into every eliminated proportion.)']
    if n_other:
        L.append(f'  {n_other} genes with an unrecognised peak-stage label are '
                 f'left out of the expression tests.')
    L.append('')

    L += ['Expression category (6A bold rows). Germline = testis, ovary, zygote, '
          'early embryo; somatic =',
          'late embryo, larvae, muscle, intestine, carcass. Denominator = genes '
          'with a stage of peak expression.',
          '  category   retained % [95% CI]     eliminated % [95% CI]   '
          'OR [95% CI]          P          BH q']
    for r in cat_tests.itertuples():
        L.append(f'  {r.category:<9} {r.pct_retained:5.1f} [{r.pct_retained_lo:4.1f}, '
                 f'{r.pct_retained_hi:4.1f}]      {r.pct_eliminated:5.1f} '
                 f'[{r.pct_eliminated_lo:4.1f}, {r.pct_eliminated_hi:4.1f}]      '
                 f'{r.odds_ratio:4.2f} [{r.or_ci_lo:4.2f}, {r.or_ci_hi:4.2f}]   '
                 f'{r.p:<10.3g} {r.q_bh:.3g}')
    L += [f'  Homogeneity, 2 x 3: chi-square = {chi2:.2f}, df = {chi_df}, '
          f'P = {chi_p:.3g}', '',
          'Expression stage (6A rows; the panel shows uncorrected P)',
          '  stage            n elim  n ret   fold    P          Bonferroni  BH q']
    for r in stage_tests.itertuples():
        L.append(f'  {r.category:<15} {r.n_eliminated_in:6d}  {r.n_retained_in:6,}  '
                 f'{r.fold_enrichment:5.2f}   {r.p:<10.3g} {r.p_bonferroni:<11.3g} '
                 f'{r.q_bh:.3g}')
    L += ['', 'dN/dS by stage of peak expression (6B), both regions pooled',
          '  row       key              n       median   Q1       Q3']
    for r in by_stage.itertuples():
        L.append(f'  {r.row:<9} {r.key:<15} {r.n:6,}  {r.median:.4f}   '
                 f'{r.q1:.4f}   {r.q3:.4f}')
    L += [f"  Pooled median {b_info['pooled_median']:.4f}. Kruskal-Wallis across "
          f"the stages: H = {kw_h:.2f}, P = {kw_p:.3g}. Germline vs somatic, "
          f"Mann-Whitney: P = {gs_p:.3g}.",
          '  Neither test is drawn on the panel.',
          f"  Drawn: x-axis stops at the {b_info['cap_pct']:g}th percentile "
          f"({b_info['x_axis_cap']:.3f}); {b_info['n_beyond_axis']:,} jittered "
          f"genes lie beyond it, on the edge.",
          f"  Jitter: {b_info['n_jitter_not_drawn']:,} points not drawn (at most "
          f"{b_info['jitter_max']} per stage); boxes and statistics use all genes.",
          '']

    if len(class_tests):
        n_has_r, n_has_e = n_steps['with a functional class']
        L += [f'Functional class (6C). Denominator = genes with a class '
              f'(retained {n_has_r:,}, eliminated {n_has_e:,}).',
              f'  {len(class_tests)} classes tested (>= {args.min_class_genes} '
              f'genes); smallest Bonferroni-adjusted P = '
              f"{class_tests['p_bonferroni'].min():.3g}; "
              f"smallest BH q = {class_tests['q_bh'].min():.3g}.",
              '  OR = conditional MLE; CI = exact, inverting Fisher\'s test '
              f"(computed by {class_tests['or_method'].iloc[0]}).",
              '  Ten smallest P:',
              '  n elim / ret    fold     OR    [95% CI]            P          '
              'Bonf.   class']
        for r in class_tests.head(10).itertuples():
            L.append(f'  {r.n_eliminated_in:5d} / {r.n_retained_in:<7,} '
                     f'{r.fold_enrichment:6.2f}  {r.odds_ratio:6.2f}  '
                     f'[{r.or_ci_lo:6.2f}, {r.or_ci_hi:6.2f}]   {r.p:<10.3g} '
                     f'{r.p_bonferroni:<7.3g} {r.category}')
        if len(forest):
            L += [f'  Drawn in 6C: {len(forest)} classes with >= '
                  f'{args.forest_min_elim} eliminated genes; '
                  f"{int(forest['n_ci_ends_clipped'].iloc[0])} interval ends "
                  f'run off the axis (arrowheads).',
                  f'  Not drawn: {len(class_tests) - len(forest)} tested '
                  f'classes with fewer eliminated genes.']
            excl_1 = forest[(forest['or_ci_lo'] > 1) | (forest['or_ci_hi'] < 1)]
            L.append(f'  Drawn classes whose 95% CI excludes 1: {len(excl_1)}'
                     + (': ' + '; '.join(excl_1['category']) if len(excl_1) else ''))
            L.append(f'  P column of 6C, Bonferroni over {len(class_tests)} '
                     'classes (top to bottom; short form as drawn, exact value):')
            for r in forest.itertuples():
                L.append(f'    {short_p(r.p_bonferroni):<8} {r.p_bonferroni:<10.3g} '
                         f'{r.category}')
            n_pass = int((class_tests['p_bonferroni'] < 0.05).sum())
            n_pass_drawn = int((forest['p_bonferroni'] < 0.05).sum())
            L.append(f'  Bonferroni P < 0.05: {n_pass} of {len(class_tests)} tested '
                     f'classes, {n_pass_drawn} of them drawn.')
        top_fold = class_tests.sort_values('fold_enrichment',
                                           ascending=False).head(5)
        L += ['  Five largest fold enrichments, all tested classes (the claim '
              'in Results):']
        for r in top_fold.itertuples():
            L.append(f'    {r.fold_enrichment:5.2f}x on {r.n_eliminated_in} '
                     f'eliminated genes, OR CI [{r.or_ci_lo:.2f}, '
                     f'{r.or_ci_hi:.2f}]  {r.category}')
        L.append('')

    pk_r, pk_e = n_steps['with a stage of peak expression']
    L += [f"dN/dS vs mean expression (S12A, file S6_dnds_vs_expression): "
          f"Spearman rho = {s6['rho']:+.4f}, "
          f"P = {s6['p']:.3g}, n = {s6['n']:,} "
          f"(retained {s6['n_retained']:,}, eliminated {s6['n_eliminated']:,}).",
          f"  {s6['n_zero_expression_drawn']:,} genes with zero mean expression "
          f"drawn at log2(0.1) = {s6['zero_expression_x']:.2f}, the left edge "
          f"(left out before 2026-10-19).",
          f"  Not drawn: {s6['n_no_dnds']:,} of the {s6['n_input']:,} matched genes "
          f"have no dN/dS; {s6['n_no_expression']:,} have dN/dS but no stage TPMs.",
          f"  6A counts {pk_r:,} retained and {pk_e:,} eliminated genes (a peak "
          f"stage; dN/dS not needed); S12A {s6['n_retained']:,} and "
          f"{s6['n_eliminated']:,}.",
          f"  y-axis stops at the {s6['cap_pct']:g}th percentile "
          f"({s6['y_cap']:.3f}), plus a 3% margin (2% below 0) so no point is "
          f"cut by an axis; above the cap, drawn as open triangles on the top "
          f"edge: "
          f"retained {s6['n_above_cap']['retained']}, eliminated "
          f"{s6['n_above_cap']['eliminated']}.", '']
    L += ['Figure audit (paper_colors.audit, plus text past the canvas edge)']
    for k, probs in audits.items():
        L.append(f'  {k}: ' + ('clean' if not probs else '; '.join(probs)))
    L += ['', 'Notes', '',
          '  6A P column: Fisher exact, two-sided, uncorrected; BH for the three '
          'categories and',
          '  Bonferroni for the ten stages are in the tables above. 6C intervals '
          'are exact conditional',
          '  95% CIs for the odds ratio. Proportions and odds ratios use matched '
          'genes only.', '']
    summary = '\n'.join(L)
    with open(os.path.join(args.outdir, 'functional_enrichment_summary.txt'),
              'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(), 'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'scipy': scipy_version,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv, 'args': vars(args),
            'inputs': {k: fingerprint(v) for k, v in
                       (('kaks', args.kaks), ('annotation', args.annotation),
                        ('bridge', args.bridge), ('gff1', args.gff1))},
            'join': {'mode': 'position bridge' if args.bridge
                     else 'name (superseded)', 'bridge': bridge_info,
                     'alignment': align_diag,
                     'shared_workbook_genes': {'retained': shared[0],
                                               'eliminated': shared[1]}},
            'gene_counts': {k: {'retained': r, 'eliminated': e}
                            for k, (r, e) in n_steps.items()},
            'germline_stages': GERMLINE, 'somatic_stages': SOMATIC,
            'figure_audit': audits,
        }, fh, indent=2, default=str)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
