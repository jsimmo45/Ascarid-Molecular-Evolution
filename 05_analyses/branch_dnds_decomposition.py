#!/usr/bin/env python3
"""
branch_dnds_decomposition.py -- Decompose branch dN/dS into dN and dS.

`branch_dnds_four.py` established that eliminated genes have elevated dN/dS on
all four terminal branches. This script asks WHICH COMPONENT drives it, per
lineage: is the elevation in the numerator (dN) with the denominator (dS)
unchanged? That is the signature of relaxed purifying selection rather than a
mutation-rate difference, and it is the mechanistic core of the paper.

No new computation. dN and dS per branch are already columns in
branch_dnds_four_annotated.tsv; this only re-analyses and plots them.

Produces Figure 5 and supplementary fig. S8 WHOLE, as one grid each (since
2026-10-22; see below): rows dN/dS, dN and dS, columns the four terminal
branches. `--figure legacy` still writes the old separate panels (5C, 5D and
their S14 counterparts; S14 is S8's frozen-map number, kept in file names).


MASKING CONVENTION -- read this before interpreting output
----------------------------------------------------------
The per-branch omega cap removes a gene from a branch when that branch's dN/dS
exceeds it. Because dN/dS = dN/dS, that cap preferentially removes LOW-dS
genes, at different rates in the retained and eliminated classes. Filtering on
dN/dS and then testing dS is therefore circular, and biases the one comparison
whose NULL result carries the mechanistic argument.

  dN, dN/dS  ->  capped subset          (matches the n in Figure 5A)
  dS         ->  ALL genes, no cap      (default; --mask-ds restores the old
                                         circular behaviour for comparison)

This applies to the stats table, both figures and the summary. The bias runs in
the expected direction: on the Toxocara branch the cap drops ~17.9% of
eliminated but only ~14.9% of retained genes, all of them low-dS, which
inflates the eliminated median and overstates the dS fold change by roughly one
star.


MULTIPLE TESTING
----------------
BH q-values are computed across the four branches within each
(grouping, metric, masking) family and reported as q_value. `--fdr-family`
chooses the family the FIGURE ASTERISKS are drawn from; the summary reports the
dS verdict under every family so a reviewer question is already answered in the
file rather than needing a rerun.


CONVERSION NOTES, 2026-08-18
----------------------------
Style pass plus the shared loader, as agreed. The analysis is unchanged: the
masking convention, the FDR families and the sensitivity table were already
correct and are preserved verbatim in behaviour.

  * `usable_mask`, `compare`, `split_and_compare`, the species order and the
    BH implementation now come from `branch_lib` / `region_lib`. This file and
    `branch_dnds_four.py` previously held four pairs of near-duplicates, so the
    two halves of Figure 5 could disagree about which genes are usable -- and
    about whether these tests need FDR correction at all.
  * Built at journal column width with sizes in mm. The old 20 x 10 inch canvas
    with bbox_inches='tight' printed its 15 pt titles at about 5.6 pt once
    placed at 185 mm, below the 6 pt floor.
  * Box alpha moved into the face colour. `patch.set_alpha(0.75)` also fades
    the black edge to light grey, and the box loses its boundary at print size.
  * Species names render italic via `paper_colors.SPECIES_LABEL`.
  * dN / dS / dN/dS in every label; no bare "omega" a reader can see.
  * provenance.json, and `make_outdir()` when --outdir is not given.


CHANGES, 2026-10-22 (figure only; no number changes)
----------------------------------------------------
Figure 5 was rebuilt as three rows of four panels (results draft notes 76-78):
5B (branch medians) duplicated 5A and 5D (fold bars) duplicated the folds
printed over 5C, so both were dropped, and 5C's two rows moved up under 5A.
Assembled by hand, that needed 5C scaled to 5A's panel width, which lengthened
its ticks, shifted its first column and left its tick labels ragged. This
script now draws the whole figure, and S8 the same way:

  * `plot_branch_grid()`: A = branch dN/dS (5A's content, which
    branch_dnds_four.py drew), B = branch dN, C = branch dS, at 185 mm with
    every axes placed in millimetres. Columns line up by construction.
  * One log decade has the same height in every row (`--decade-mm`, default
    7 mm), so a fold change is the same vertical distance in A, B and C. Row
    heights follow from each row's decades; the zero band's share is measured
    (`_zero_band_fraction`), not assumed.
  * Species names once, as column headers over row A; y tick labels on the
    first column only; the dN/dS = 1 rule across row A, labelled on the last
    panel; panel letters A-C at 10 pt (on by default; `--no-panel-letters`).
  * Row A's asterisks come from the same family as B and C: within the metric,
    across the four branches and both labellings (n = 8), which is what the
    legend and Methods state. branch_dnds_four.py's 5A used the four branches
    of one labelling; at these q values the asterisks are the same.
  * `check_grid_spacing()` measures the drawn text and stops the run if a row's
    tick labels would touch the next row's fold line, or anything leaves the
    page. The measured gaps and decade heights go to the summary.
  * Legend notes for the grid in the summary: per row and branch, n, zeros,
    values above the axis, fold and q, and which lower whiskers reach the band.
  * Statistics, `dnds_decomposition_stats.tsv` and the summary's tables are
    unchanged; `--figure legacy` output is identical to before (its jitter has
    its own generator, seeded as before).

Later the same day (figure only; no number changes):
  * The grid's jittered points are a stratified sample (`_jitter_sample`): at
    most 1,500 per box as before, but the genes in the band at 0 and the rest
    are sampled separately at the same rate, and a band that holds any gene
    gets at least one point. A plain random 1,500 of the ~5,800 retained genes
    drew a lone zero about one time in four, so the first real run printed a
    count of 1 beside three empty bands (Figure 5B and S8B, Toxocara dN; S8C,
    P. univalens dS). The run now stops if a band with a count is drawn
    without a point, checked on the drawn points themselves. Boxes, whiskers,
    counts and tests always used every gene and are unchanged; only which
    points are drawn changes. The legacy panels are untouched.


CHANGES, 2026-10-14 (figure only; no number changes)
----------------------------------------------------
  * 5C's log axes are labelled 0.001 ... 10 instead of 10^-3 ... 10^1
    (`branch_lib.decade_label`). The power form was mathtext, which reached
    the assembled figure in DejaVu Sans with 4.9 pt exponents.
  * The counts beside the zero bands are 7 pt (FONT['annotation']); they were
    6 pt.


CHANGES, 2026-10-10
-------------------
  * Upstream: `branch_dnds_four.py` now gives codeml the UNROOTED tree. The
    table it wrote before that date split the Toxocara branch in two at the
    root, so every Toxocara value in 5C, 5D and S14 changes on the rerun.
    `branch_lib.load_annotated` refuses a table made with the rooted tree.
  * The capped rows (dN/dS, and dN) also drop a branch with dN = dS = 0
    (`branch_lib.usable_mask`): its dN/dS was wherever codeml's randomly
    started optimiser stopped, so whether it passed the cap was chance. The
    uncapped dS rows keep these genes, as dS = 0 is a real estimate.
  * Zeros. dN or dS printed as 0 cannot sit on a log axis; the panels said
    "N at 0" once per panel while drawing those genes' whiskers -- and on the
    Toxocara dN panel the whole retained box -- down to the axis floor. They
    are now drawn in a band at 0 below a scale break
    (`branch_lib.zero_band_axis`, shared with 5A), with the count beside each
    box. Whiskers are the 5th and 95th percentiles; 1.5 IQR whiskers on a log
    axis reach the minimum of any right-skewed sample.
  * The stars and fold sit above each panel. Inside it they ran into the
    whisker caps on the Toxocara dN panel and into the jitter elsewhere.
  * The footer is gone: the cap, the dS masking, the asterisk family, the zero
    counts and the clipped counts are written to the summary as legend notes.
  * 5D: the dN bars were region red, which read as "eliminated" on a panel of
    eliminated/retained ratios. Both metrics are now drawn in ink, as
    paper_colors requires for every metric.

Later the same day (cosmetic only; no number changes):
  * 5C and 5D: the bottom margin is measured (`branch_lib.fit_bottom_margin`,
    1 mm under the lowest label) instead of a fixed fraction. 5D had 5-6 mm of
    empty figure under its labels, 5C about 2.
  * 5D: the dashed rule at 1 was drawn over the bars (same zorder, added
    later). It now sits beneath them, and the bars are opaque -- the same
    tones as before, composited onto white -- so the rule does not show
    through the light dS bars.
  * Summary: the closing "claim to make in the text" was a fixed sentence
    saying no dS difference survives correction. Since the unrooted-tree
    rerun two Toxocara dS tests do, and the sentence contradicted the table
    above it. It is now computed from the figure's correction family.

Usage
-----
  python branch_dnds_decomposition.py \
    --annotated results/branch_dnds_four_<stamp>/branch_dnds_four_annotated.tsv \
    --omega-cap 5 --outdir results/dnds_decomposition_<stamp>

  Writes fig5_branch_grid_as_region.{svg,pdf,png} (Figure 5) and
  figS14_branch_grid_tox_region.{svg,pdf,png} (supplementary fig. S8).
  Add --figure both for the old 5C/5D files as well.

Dependencies: pandas, numpy, scipy, matplotlib
"""

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

from paper_colors import (REGION, NEUTRAL, COL, COL_MM, MAX_HEIGHT_MM,
                          LINEWIDTH, apply_style, audit, save, make_outdir, mm,
                          SPECIES_LABEL, FONT, panel_letter)
from region_lib import (benjamini_hochberg, fmt_p, stars, fingerprint,
                        git_commit, setup_logging)
from branch_lib import (SPECIES_ORDER, SPECIES_KEY, METRICS, METRIC_LABEL,
                        usable_mask, split_and_compare, load_annotated,
                        at_zero, as_plotted, zero_band_limits, zero_band_axis,
                        mark_axis_break, fit_bottom_margin)

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Figure 5D bar colours. dN and dS are metrics, not regions, so neither takes
#: a region colour -- paper_colors draws every metric in ink. Until 2026-10-10
#: the dN bars were REGION['eliminated'] red, which read as "eliminated" on a
#: panel whose bars are each an eliminated/retained ratio.
COL_DN_BAR = NEUTRAL['ink']
COL_DS_BAR = NEUTRAL['ink_light']

#: Strength of each bar colour. Applied by compositing onto white rather than
#: as alpha, so the bars are opaque and the reference rule beneath them does
#: not show through (see on_white()).
DN_BAR_STRENGTH = 0.85
DS_BAR_STRENGTH = 0.45


def on_white(color, strength):
    """`color` at `strength` composited onto a white page, as an opaque RGB
    tuple. Looks the same as facecolor alpha on white paper, but hides what is
    drawn underneath -- here the dashed rule at 1, which a translucent light
    dS bar let through."""
    if not 0.0 <= strength <= 1.0:
        raise ValueError(f'strength must be in [0, 1], got {strength}')
    r, g, b, _ = to_rgba(color)
    return tuple(strength * c + (1.0 - strength) for c in (r, g, b))

STATS_COLS = ['ret_n', 'elim_n', 'ret_median', 'elim_median', 'fold_change',
              'p_value', 'branch', 'metric', 'grouping', 'masked', 'q_value',
              'q_star']

# Multiple-testing families available for the figure asterisks. Keys are
# --fdr-family values; the tuple names the stats columns defining one family.
#
# WHY THIS MATTERS MORE THAN IT LOOKS: BH computes q = p * m / rank. Adding the
# eight astronomically significant dN tests to the dS family raises every dS
# test's RANK by eight while only doubling m, so m/rank FALLS and the dS tests
# come out MORE significant, not less. On the v2 data, pooling dN with dS turns
# three dS comparisons significant (q = 0.017, 0.045, 0.045) that are null in
# every dS-only family. dN and dS answer different questions and must not share
# a correction family. 'mixed' exists only so that claim can be demonstrated.
FDR_FAMILIES = {
    'metric':          ('metric', 'masked'),               # default; dS -> n=8
    'metric-grouping': ('grouping', 'metric', 'masked'),   # dS -> n=4
    'mixed':           ('grouping', 'masked'),             # NOT RECOMMENDED
}


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------

def add_star_qvalues(stats, family):
    """Add `q_star`, the statistic the figure asterisks are drawn from."""
    if family == 'none':
        stats['q_star'] = stats['p_value'].astype(float)
        return stats
    if family not in FDR_FAMILIES:
        raise ValueError(f'unknown --fdr-family {family!r}; choose from '
                         f"{sorted(FDR_FAMILIES)} or 'none'")
    keys = list(FDR_FAMILIES[family])
    stats['q_star'] = np.nan
    for _, idx in stats.groupby(keys).groups.items():
        stats.loc[idx, 'q_star'] = benjamini_hochberg(
            stats.loc[idx, 'p_value'].to_numpy(dtype=float))
    return stats


def star_column(star_stat):
    """Stats column the asterisks read. Single source of truth."""
    return 'p_value' if star_stat == 'p' else 'q_star'


def star_caption(star_stat, family, n_per_family=None):
    """One-line description of the asterisks, for the figure footer."""
    if star_stat == 'p' or family == 'none':
        return 'asterisks: uncorrected P (Mann-Whitney U)'
    scope = {'metric': 'within each metric, across branches and labellings',
             'metric-grouping': 'within each metric and region labelling',
             'mixed': 'across dN and dS together [NOT RECOMMENDED]'}[family]
    extra = f', n = {n_per_family}' if n_per_family else ''
    return f'asterisks: BH-FDR q, corrected {scope}{extra}'


def ds_row_masking(mask_ds):
    """Which `masked` value the dS rows and panels use. Single source."""
    return bool(mask_ds)


def build_stats(df, region_col, omega_cap, ds_max=np.inf):
    """Per-branch retained-vs-eliminated stats for dN/dS, dN and dS.

    Two families of rows:
      masked=True   capped subset  (omega, dN, dS)
      masked=False  no cap         (dN, dS)

    The masked=False dS row is the one to quote for any dS claim. The
    masked=False dN row is included so a reviewer can confirm the dN result is
    not an artifact of the cap either.
    """
    rows = []
    for sp in SPECIES_ORDER:
        fold_masked = {}
        for metric in METRICS:
            keep = usable_mask(df, sp, metric, omega_cap, True, ds_max)
            res, _, _ = split_and_compare(df, keep, region_col, f'{sp}_{metric}')
            res.update({'branch': sp, 'metric': metric,
                        'grouping': region_col, 'masked': True})
            rows.append(res)
            fold_masked[metric] = res['fold_change']

        f_dn, f_ds = fold_masked.get('dN'), fold_masked.get('dS')
        if np.isfinite(f_dn) and np.isfinite(f_ds) and f_ds > 0:
            logging.info('  %-14s fold(dN/dS)=%.3f vs fold(dN)/fold(dS)=%.3f '
                         '[capped]', sp, fold_masked['omega'], f_dn / f_ds)

        for metric in ('dN', 'dS'):
            keep = usable_mask(df, sp, metric, omega_cap, False, ds_max)
            res, _, _ = split_and_compare(df, keep, region_col, f'{sp}_{metric}')
            res.update({'branch': sp, 'metric': metric,
                        'grouping': region_col, 'masked': False})
            rows.append(res)

        # How much of the branch the cap removes, per class. An asymmetry here
        # is exactly what biases a post-cap dS comparison.
        capped = usable_mask(df, sp, 'dS', omega_cap, True, ds_max)
        allg = usable_mask(df, sp, 'dS', omega_cap, False, ds_max)
        for cls in ('retained', 'eliminated'):
            in_cls = df[region_col] == cls
            n_all = int((allg & in_cls).sum())
            n_cap = int((capped & in_cls).sum())
            if n_all:
                logging.info('    %-14s %-10s cap drops %d/%d (%.1f%%)',
                             sp, cls, n_all - n_cap, n_all,
                             100 * (n_all - n_cap) / n_all)

    stats = pd.DataFrame(rows)
    stats['q_value'] = np.nan
    for _, idx in stats.groupby(['grouping', 'metric', 'masked']).groups.items():
        stats.loc[idx, 'q_value'] = benjamini_hochberg(
            stats.loc[idx, 'p_value'].to_numpy(dtype=float))
    # Filled by add_star_qvalues() only AFTER every grouping has been built,
    # because the default FDR family pools across region labellings.
    stats['q_star'] = np.nan
    return stats[STATS_COLS]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_dn_ds_panels(df, region_col, outdir, stem, omega_cap, mask_ds, ds_max,
                      rng, stats, star_stat='q', fdr_family='metric',
                      fam_n=None, yscale='log', width='double',
                      height_mm=92.0, panel_letters=False):
    """Figure 5C -- 2 x 4 grid, dN on top and dS below, per terminal branch.

    dN uses the capped subset; dS uses all genes unless mask_ds. The two rows
    therefore carry DIFFERENT n by design -- that is the point of the masking
    fix, and the n in each tick label states which set was used.

    Significance is LOOKED UP from `stats` rather than recomputed here, so the
    asterisks on this figure and on the fold-change bars can never disagree.

    Returns {'zeros': {(metric, branch, region): n}, 'above': {(metric,
    branch, region): n}, 'limits': {metric: (lo, hi)}} for the legend notes.
    """
    apply_style()
    scol = star_column(star_stat)

    def lookup(sp, metric, masked):
        r = stats[(stats['branch'] == sp) & (stats['metric'] == metric)
                  & (stats['masked'] == masked)]
        if not len(r):
            return np.nan, np.nan
        return float(r[scol].iloc[0]), float(r['fold_change'].iloc[0])

    fig, axes = plt.subplots(2, 4, figsize=(COL[width], mm(height_mm)),
                             layout='none')

    info = {'zeros': {}, 'above': {}, 'limits': {}}
    log = yscale == 'log'
    for row, metric in enumerate(('dN', 'dS')):
        apply_cap = True if metric == 'dN' else ds_row_masking(mask_ds)

        # ONE SCALE PER ROW, AND IT HAS TO BE LOG.
        #
        # The old panel shared a LINEAR limit across the row, set by the
        # largest branch. Toxocara is the deepest comparison in the tree, so
        # its absolute dN and dS run roughly an order of magnitude above the
        # ascarid branches -- and the shared linear limit pressed the other
        # three panels into the bottom tenth of the axis, where a 1.8x
        # difference is invisible. Absolute dN and dS are not comparable
        # across branches of very different length in the first place; that is
        # the same trap recorded for the Toxocara panels of S12 and for the
        # genome traces.
        #
        # A log axis keeps the comparability the shared limit was reaching for
        # -- a fold change is a constant vertical offset at any branch length
        # -- while leaving every panel readable. Zeros get the band below the
        # scale break (branch_lib.zero_band_axis, the same as 5A). The lower
        # limit is the decade of the smallest non-zero value, so nothing
        # non-zero is clipped below; the upper is the decade above the 99.5th
        # percentile, and what lies above it is counted for the legend.
        vals = []
        for sp in SPECIES_ORDER:
            keep = usable_mask(df, sp, metric, omega_cap, apply_cap, ds_max)
            v = pd.to_numeric(df[keep][f'{sp}_{metric}'],
                              errors='coerce').dropna().to_numpy()
            if len(v):
                vals.append(v)
        if log:
            ylo, yhi, _ = zero_band_limits(vals, metric)
        else:
            ylo = 0.0
            yhi = (max(np.percentile(v, 97) for v in vals) * 1.30
                   if vals else 1.0)
        info['limits'][metric] = (ylo, yhi)

        for col, sp in enumerate(SPECIES_ORDER):
            ax = axes[row, col]
            keep = usable_mask(df, sp, metric, omega_cap, apply_cap, ds_max)
            _, ret, elim = split_and_compare(df, keep, region_col,
                                             f'{sp}_{metric}')

            if len(ret) < 5 or len(elim) < 5:
                ax.text(0.5, 0.5, 'insufficient data', ha='center', va='center',
                        transform=ax.transAxes, fontsize=FONT['annotation'])
                ax.set_xticks([])
                if row == 0:
                    ax.set_title(SPECIES_LABEL[SPECIES_KEY[sp]],
                                 fontsize=FONT['panel_title'])
                continue

            data = [as_plotted(ret, metric), as_plotted(elim, metric)]
            bp = ax.boxplot(data, positions=[0, 1], widths=0.55,
                            patch_artist=True, showfliers=False, whis=(5, 95),
                            zorder=2)
            for patch, key in zip(bp['boxes'], ('retained', 'eliminated')):
                patch.set_facecolor(to_rgba(REGION[key], 0.40))
                patch.set_edgecolor('black')
                patch.set_alpha(None)
            for med in bp['medians']:
                med.set_color('black')
                med.set_zorder(5)
            for art in bp['whiskers'] + bp['caps']:
                art.set_zorder(4)

            for i, (v, key) in enumerate(zip(data, ('retained', 'eliminated'))):
                info['zeros'][(metric, sp, key)] = int(at_zero(v, metric).sum())
                info['above'][(metric, sp, key)] = int((v > yhi).sum())
                if len(v) > 1500:
                    v = v[rng.choice(len(v), 1500, replace=False)]
                ax.scatter(np.full(len(v), i) + rng.normal(0, 0.07, len(v)), v,
                           s=1.6, alpha=0.18, color=REGION[key], linewidths=0,
                           zorder=3, rasterized=True)

            sig, fold = lookup(sp, metric, apply_cap)

            ax.set_xlim(-0.6, 1.6)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f'Retained\n(n = {len(ret):,})',
                                f'Eliminated\n(n = {len(elim):,})'])
            if log:
                mark_axis_break(ax, zero_band_axis(ax, ylo, yhi))
                for i, key in enumerate(('retained', 'eliminated')):
                    z = info['zeros'][(metric, sp, key)]
                    if z:
                        # House annotation size (7 pt); it was 1 pt smaller.
                        ax.text(i + 0.30, 0.0, f'{z:,}', ha='left',
                                va='center', zorder=6,
                                fontsize=FONT['annotation'],
                                color=NEUTRAL['ink_light'])
            else:
                ax.set_ylim(ylo, yhi)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            if col == 0:
                ax.set_ylabel(METRIC_LABEL[metric])
            if row == 0:
                # Padded to sit above the stars-and-fold line.
                ax.set_title(SPECIES_LABEL[SPECIES_KEY[sp]],
                             fontsize=FONT['panel_title'],
                             pad=FONT['annotation'] + 5)
            # Above the panel, not inside it: inside, it ran into both whisker
            # caps on the Toxocara dN panel and into the jitter elsewhere.
            ax.text(0.5, 1.03, f'{stars(sig)}  {fold:.2f}×',
                    transform=ax.transAxes, ha='center', va='bottom',
                    fontsize=FONT['annotation'])

    if panel_letters:
        panel_letter(axes[0, 0], 'C')

    # No footer: the cap, the dS masking, the asterisk family and the counts
    # of zeros and clipped values go to the summary as legend notes. bottom is
    # a starting value only; fit_bottom_margin() sets it from the drawn labels.
    fig.subplots_adjust(left=0.085, right=0.995, top=0.915, bottom=0.105,
                        wspace=0.28, hspace=0.42)
    logging.debug('  %s: bottom margin moved %.1f mm', stem,
                  fit_bottom_margin(fig))
    audit(fig)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    for metric in ('dN', 'dS'):
        n_above = sum(n for (m, _, _), n in info['above'].items() if m == metric)
        if n_above:
            logging.info('  %s axis stops at %g: %s values above it are '
                         'clipped (legend notes in the summary)', metric,
                         info['limits'][metric][1], f'{n_above:,}')
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return info


def plot_fold_changes(stats, region_col, outdir, stem, mask_ds,
                      star_stat='q', fdr_family='metric', fam_n=None,
                      width='single', height_mm=56.0, panel_letters=False):
    """Figure 5D -- fold change (eliminated / retained) for dN and dS per branch.

    dN comes from the capped rows, dS from the unmasked rows (unless mask_ds),
    so the bars agree with the panels and with the quotable stats. A dS bar at
    1 beside a tall dN bar is the relaxed-selection signature, and the dashed
    rule at 1 is the reference that makes it readable.
    """
    apply_style()
    s = stats[stats['grouping'] == region_col]
    ds_masked = ds_row_masking(mask_ds)
    scol = star_column(star_stat)

    def pick(sp, metric, masked):
        r = s[(s['branch'] == sp) & (s['metric'] == metric)
              & (s['masked'] == masked)]
        if not len(r):
            return np.nan, np.nan
        return float(r['fold_change'].iloc[0]), float(r[scol].iloc[0])

    dn = [pick(sp, 'dN', True) for sp in SPECIES_ORDER]
    ds = [pick(sp, 'dS', ds_masked) for sp in SPECIES_ORDER]
    dn_f, dn_p = [a for a, _ in dn], [b for _, b in dn]
    ds_f, ds_p = [a for a, _ in ds], [b for _, b in ds]

    finite = [v for v in dn_f + ds_f if np.isfinite(v)]
    ymax = (max(finite) * 1.34) if finite else 1.0

    x = np.arange(len(SPECIES_ORDER))
    w = 0.36
    fig, ax = plt.subplots(figsize=(COL[width], mm(height_mm)),
                           layout='none')
    # Rule beneath the bars: zorder 1 against the bars' 2. Before 2026-10-10
    # both were at 1 and the rule, added later, was drawn across every bar.
    # The bars are opaque (on_white), so the rule does not show through the
    # light dS bars either; it reads between the bars, where it is needed.
    ax.bar(x - w / 2, dn_f, w, label='dN', zorder=2,
           facecolor=on_white(COL_DN_BAR, DN_BAR_STRENGTH), edgecolor='black',
           linewidth=0.6)
    ax.bar(x + w / 2, ds_f, w, label='dS', zorder=2,
           facecolor=on_white(COL_DS_BAR, DS_BAR_STRENGTH), edgecolor='black',
           linewidth=0.6)
    ax.axhline(1.0, color=NEUTRAL['rule'], ls='--', lw=0.9, zorder=1)

    for xi in range(len(SPECIES_ORDER)):
        for off, val, p in ((-w / 2, dn_f[xi], dn_p[xi]),
                            (w / 2, ds_f[xi], ds_p[xi])):
            if not np.isfinite(val):
                continue
            ax.text(xi + off, val + 0.02 * ymax, f'{val:.2f}', ha='center',
                    va='bottom', fontsize=FONT['annotation'])
            ax.text(xi + off, val + 0.11 * ymax, stars(p), ha='center',
                    va='bottom', fontsize=FONT['annotation'])

    ax.set_xticks(x)
    ax.set_xticklabels([SPECIES_LABEL[SPECIES_KEY[sp]] for sp in SPECIES_ORDER])
    ax.set_ylabel('Fold change\n(eliminated / retained)')
    ax.set_ylim(0, ymax)
    # Reserve the space, then place: loc='upper left' put the key over the
    # tallest dN bar, which audit() flags. ncol=2 keeps it one short row.
    ax.legend(frameon=False, fontsize=FONT['legend'], ncol=2,
              handlelength=1.2, loc='lower left', bbox_to_anchor=(0.0, 1.005),
              columnspacing=1.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if panel_letters:
        panel_letter(ax, 'D')

    # bottom is a starting value only; fit_bottom_margin() sets it from the
    # drawn labels. The fixed 0.20 left 5-6 mm of empty figure under them.
    fig.subplots_adjust(left=0.20, right=0.99, top=0.90, bottom=0.20)
    logging.debug('  %s: bottom margin moved %.1f mm', stem,
                  fit_bottom_margin(fig))
    audit(fig)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))


# ---------------------------------------------------------------------------
# Figure 5 / S8 as one grid (2026-10-22)
# ---------------------------------------------------------------------------

#: Rows of the grid, top to bottom: (metric, panel letter).
GRID_ROWS = (('omega', 'A'), ('dN', 'B'), ('dS', 'C'))

#: Grid geometry, in mm at final size. Columns: `left` holds the y label and
#: the y tick labels of the first column (13.9 mm, the spine position of the
#: hand-built Figure 5 v4), `col_gap` separates the columns and `right` is left
#: at the right edge. Rows, top to bottom: `header` holds the species names,
#: `fold` the stars-and-fold line over each row, `xlab` the two-line tick labels
#: under it and `row_gap` the space between one row's tick labels and the next
#: row's fold line; `bottom` is what remains under the last row. `ylabel_x` is
#: where the right edge of each row's (rotated) y label sits.
GRID_MM = {'left': 13.9, 'right': 0.3, 'col_gap': 4.0,
           'header': 4.5, 'fold': 4.0, 'xlab': 8.0, 'row_gap': 2.5,
           'bottom': 0.5, 'ylabel_x': 3.9}

#: Jittered genes drawn per box at most; boxes always use every gene. Larger
#: boxes are sampled by `_jitter_sample`, which never leaves a band empty.
GRID_JITTER_MAX = 1500

#: Points between the top of an axes and the baseline region of its fold line.
GRID_FOLD_OFFSET_PT = 2.0

#: Smallest gap check_grid_spacing() accepts between neighbouring text blocks.
GRID_MIN_GAP_MM = 0.8


def _zero_band_fraction(lo, hi):
    """Axes fraction at which zero_band_axis() puts `lo` -- the share of an
    axes' height taken by the band for 0 and the margin under it.

    Measured on a scratch axes rather than derived from zero_band_axis()'s
    constants (linscale 0.3, a lower limit of -0.6 lo), so a change to the band
    cannot leave the grid's decades unequal unnoticed.
    """
    scratch = plt.figure(figsize=(1.0, 1.0), layout='none')
    try:
        return float(zero_band_axis(scratch.add_axes([0.0, 0.0, 1.0, 1.0]),
                                    lo, hi))
    finally:
        plt.close(scratch)


def _grid_row_data(df, region_col, metric, omega_cap, mask_ds, ds_max):
    """Values, masking and axis limits for one row of the grid.

    dN/dS and dN use the capped subset; dS uses every quartet unless --mask-ds
    (see the module docstring). The limits are pooled over the four branches,
    as the separate panels pooled them: over the two region classes for dN/dS
    (Figure 5A, with the cap as the ceiling) and over every usable value for
    dN and dS (Figure 5C), so the axes match the figures they replace.
    """
    masked = True if metric in ('omega', 'dN') else ds_row_masking(mask_ds)
    per_sp, pooled = {}, []
    for sp in SPECIES_ORDER:
        keep = usable_mask(df, sp, metric, omega_cap, masked, ds_max)
        _, ret, elim = split_and_compare(df, keep, region_col, f'{sp}_{metric}')
        per_sp[sp] = (ret.to_numpy(dtype=float), elim.to_numpy(dtype=float))
        if metric == 'omega':
            pooled.extend(per_sp[sp])
        else:
            v = pd.to_numeric(df[keep][f'{sp}_{metric}'],
                              errors='coerce').dropna().to_numpy(dtype=float)
            if len(v):
                pooled.append(v)
    hi_cap = omega_cap if metric == 'omega' else np.inf
    lo, hi, _ = zero_band_limits(pooled, metric, hi_cap=hi_cap)
    return {'metric': metric, 'masked': masked, 'per_sp': per_sp,
            'lo': lo, 'hi': hi, 'n_dec': int(round(np.log10(hi / lo)))}


def _jitter_sample(v, metric, n_max, rng):
    """The values of one box to draw as jittered points, and how many of them
    lie in the band at 0.

    A box of at most `n_max` genes is drawn whole. A larger one is sampled to
    `n_max`, stratified: the genes in the band (`at_zero`) and the rest are
    sampled separately at the same rate, n_max / n, so the band keeps its
    share of the points (to within rounding), and a band that holds any gene
    gets at least one point.

    Why: the box statistics and the count beside the band use every gene, but
    a plain random sample of 1,500 from ~5,800 retained genes includes a lone
    zero only about one time in four, and the first real run printed "1"
    beside three empty bands. The minimum of one costs the rest at most one
    point.

    `v` is a box's values as plotted (`as_plotted`: zeros exactly 0). Returns
    (sampled values in their original order, number of them in the band).
    """
    v = np.asarray(v, dtype=float)
    band = at_zero(v, metric)
    n, n_band = len(v), int(band.sum())
    if n <= n_max:
        return v, n_band
    rate = n_max / n
    k_band = min(n_band, max(1, int(round(n_band * rate)))) if n_band else 0
    k_rest = min(n_max - k_band, n - n_band)
    parts = []
    if k_band:
        parts.append(rng.choice(np.flatnonzero(band), k_band, replace=False))
    if k_rest:
        parts.append(rng.choice(np.flatnonzero(~band), k_rest, replace=False))
    idx = np.sort(np.concatenate(parts)) if parts else np.array([], dtype=int)
    return v[idx], k_band


def _extent_mm(artist, renderer, dpi):
    """Drawn extent of an artist as (x0, y0, x1, y1) in mm from the figure's
    lower-left corner; None for an artist with nothing drawn."""
    bb = artist.get_window_extent(renderer)
    if bb.width <= 0 and bb.height <= 0:
        return None
    k = 25.4 / dpi
    return (float(bb.x0 * k), float(bb.y0 * k), float(bb.x1 * k),
            float(bb.y1 * k))


def _touches_data(ax, text, renderer, max_points=3):
    """True if `text` covers a box, a whisker, cap or median, or more than
    `max_points` jittered points in `ax`. A few faint points (alpha 0.2,
    1.3 pt) under a label are tolerated; the structure of a box plot is not.
    Dashed lines (the reference rule the label belongs to) do not count."""
    tb = text.get_window_extent(renderer)
    n_points = 0
    for coll in ax.collections:
        pts = np.asarray(coll.get_offsets(), dtype=float)
        if not len(pts):
            continue
        pts = pts[np.isfinite(pts).all(axis=1)]
        d = ax.transData.transform(pts)
        n_points += int(((d[:, 0] >= tb.x0) & (d[:, 0] <= tb.x1)
                         & (d[:, 1] >= tb.y0) & (d[:, 1] <= tb.y1)).sum())
    if n_points > max_points:
        return True
    for patch in ax.patches:
        if patch.get_window_extent(renderer).overlaps(tb):
            return True
    for line in ax.lines:
        if line.get_linestyle() not in ('-', 'solid'):
            continue
        xy = np.asarray(line.get_xydata(), dtype=float)
        xy = xy[np.isfinite(xy).all(axis=1)]
        if len(xy) < 1:
            continue
        d = ax.transData.transform(xy)
        if len(d) > 1:                     # densify: a segment can cross a box
            t = np.linspace(0.0, 1.0, 50)[:, None, None]
            d = (d[:-1][None] + t * (d[1:] - d[:-1])[None]).reshape(-1, 2)
        if ((d[:, 0] >= tb.x0) & (d[:, 0] <= tb.x1)
                & (d[:, 1] >= tb.y0) & (d[:, 1] <= tb.y1)).any():
            return True
    return False


def place_rule_label(fig, axes, label='dN/dS = 1', y=1.0):
    """Direct-label the dashed rule at `y` at the right end of a panel.

    The rule crosses every panel of the row, so the label can sit in any of
    them. It goes on the rightmost panel where it covers no data -- normally
    the last, as on the hand-built Figure 5 -- and moves left only when a box,
    whisker or point lies under it. Returns (panel index, clear); clear is
    False, with a warning, when no panel has room, and the index is None when
    no panel has data (nothing drawn).
    """
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    kw = dict(ha='right', va='bottom', fontsize=FONT['annotation'],
              color=NEUTRAL['rule'], zorder=6)
    drawn = [k for k, ax in enumerate(axes) if ax.collections]   # has data
    for k in reversed(drawn):
        t = axes[k].text(1.58, y, label, **kw)
        if not _touches_data(axes[k], t, ren):
            return k, True
        t.remove()
    if not drawn:
        logging.warning('  no panel of the row has data; "%s" not drawn.', label)
        return None, False
    axes[drawn[-1]].text(1.58, y, label, **kw)
    logging.warning('  the "%s" label covers data on every panel; left on the '
                    'rightmost one with data. Check it on the figure.', label)
    return drawn[-1], False


def check_grid_spacing(fig, rows, headers, letters, min_gap_mm=GRID_MIN_GAP_MM):
    """Measure the drawn grid and stop if text blocks touch or leave the page.

    The rows are placed from fixed allowances (GRID_MM), which are estimates of
    how much room two-line tick labels and a fold line take. This checks the
    estimate on the drawn figure instead of trusting it: the species headers
    against row A's fold line, each row's tick labels against the next row's
    fold line and panel letter, the last row's labels against the bottom edge,
    and every text against the sides. It also measures one decade on every row,
    which must be the same height throughout.

    Returns a dict of the measurements for the summary; raises SystemExit, with
    the numbers, on any failure.
    """
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    dpi = fig.dpi
    w_mm = fig.get_figwidth() * 25.4
    h_mm = fig.get_figheight() * 25.4
    problems, out = [], {'gaps_mm': {}, 'decade_mm': {}}

    def lowest(texts):
        ys = [e[1] for e in (_extent_mm(t, ren, dpi) for t in texts) if e]
        return min(ys) if ys else None

    def highest(texts):
        ys = [e[3] for e in (_extent_mm(t, ren, dpi) for t in texts) if e]
        return max(ys) if ys else None

    for i, r in enumerate(rows):
        tick_texts = [t for ax in r['axes'] for t in ax.get_xticklabels()
                      if t.get_text()]
        r['_ticks_low'] = lowest(tick_texts)
        r['_fold_high'] = highest(r['fold_texts'])
        if r['_fold_high'] is None:        # every panel 'insufficient data'
            r['_fold_high'] = max(_extent_mm(ax.patch, ren, dpi)[3]
                                  for ax in r['axes'])
        if r['_ticks_low'] is None:
            r['_ticks_low'] = min(_extent_mm(ax.patch, ren, dpi)[1]
                                  for ax in r['axes'])
        # One decade, measured between the two lowest labelled decades.
        ax = r['axes'][0]
        y0 = ax.transData.transform((0, r['lo']))[1]
        y1 = ax.transData.transform((0, r['lo'] * 10.0))[1]
        out['decade_mm'][r['metric']] = round(float((y1 - y0) * 25.4 / dpi), 3)

    top_texts = list(headers) + ([letters[0]] if letters else [])
    gap = lowest(headers) - rows[0]['_fold_high']
    out['gaps_mm']['header/A'] = round(gap, 2)
    if gap < min_gap_mm:
        problems.append(f'species headers sit {gap:.2f} mm above row A\'s fold '
                        f'line (minimum {min_gap_mm} mm)')
    for upper, lower in zip(rows[:-1], rows[1:]):
        nxt = [*lower['fold_texts']] + ([lower['letter_text']]
                                        if lower.get('letter_text') else [])
        gap = upper['_ticks_low'] - highest(nxt)
        key = f"{upper['letter']}/{lower['letter']}"
        out['gaps_mm'][key] = round(gap, 2)
        if gap < min_gap_mm:
            problems.append(f'row {upper["letter"]} tick labels sit {gap:.2f} mm '
                            f'above row {lower["letter"]}\'s fold line '
                            f'(minimum {min_gap_mm} mm)')
    bottom = rows[-1]['_ticks_low']
    out['gaps_mm']['bottom'] = round(bottom, 2)
    if bottom < 0.2:
        problems.append(f'the last row\'s tick labels end {bottom:.2f} mm from '
                        'the bottom edge')
    top = highest(top_texts)
    out['gaps_mm']['top'] = round(float(h_mm - top), 2)
    if top > h_mm + 1e-6:
        problems.append(f'text runs {top - h_mm:.2f} mm above the top edge')
    tb = fig.get_tightbbox(ren)            # inches
    if tb.x0 * 25.4 < -1e-6 or tb.x1 * 25.4 > w_mm + 1e-6:
        problems.append(f'drawn extent {tb.x0 * 25.4:.2f}-{tb.x1 * 25.4:.2f} mm '
                        f'leaves the {w_mm:.0f} mm page')

    decs = list(out['decade_mm'].values())
    if max(decs) - min(decs) > 0.05:
        problems.append(f'decade heights differ between rows: {out["decade_mm"]}')

    if problems:
        raise SystemExit('ERROR: grid layout check failed:\n  '
                         + '\n  '.join(problems)
                         + '\nAdjust GRID_MM (xlab, row_gap, header) and rerun.')
    logging.info('  grid spacing (mm): %s; decade (mm): %s',
                 out['gaps_mm'], out['decade_mm'])
    return out


def plot_branch_grid(df, region_col, outdir, stem, omega_cap, mask_ds, ds_max,
                     rng, stats, star_stat='q', decade_mm=7.0,
                     panel_letters=True):
    """Figure 5 (as_region) or supplementary fig. S8 (tox_region), whole.

    Three rows of four panels at 185 mm: A branch dN/dS, B branch dN, C branch
    dS, one column per terminal branch in SPECIES_ORDER. Every axes is placed
    in millimetres (GRID_MM), so the columns of all three rows share one left
    edge and one width, and nothing is scaled after drawing.

    Each row is a log axis with the band for 0 below a scale break
    (zero_band_axis). A decade is `decade_mm` tall in every row, so a fold
    change is the same vertical distance in A, B and C; each row's height
    follows from its number of decades and the band's measured share.

    Significance is looked up from `stats` (the q_star family) rather than
    recomputed, as in the separate panels. Row A therefore takes its asterisks
    from the dN/dS family across both labellings (n = 8), the family the legend
    and Methods state.

    Returns the figure facts and legend counts for the summary.
    """
    apply_style()
    scol = star_column(star_stat)
    G = GRID_MM

    def lookup(sp, metric, masked):
        r = stats[(stats['branch'] == sp) & (stats['metric'] == metric)
                  & (stats['masked'] == masked)]
        if not len(r):
            return np.nan, np.nan
        return float(r[scol].iloc[0]), float(r['fold_change'].iloc[0])

    # Rows: values, limits and heights, before anything is drawn.
    rows = []
    for metric, letter in GRID_ROWS:
        r = _grid_row_data(df, region_col, metric, omega_cap, mask_ds, ds_max)
        r['letter'] = letter
        r['f_lo'] = _zero_band_fraction(r['lo'], r['hi'])
        if not 0.0 < r['f_lo'] < 0.5:
            raise SystemExit(f"ERROR: zero band takes {r['f_lo']:.2f} of row "
                             f"{letter}'s height; zero_band_axis() changed?")
        r['h_mm'] = decade_mm * r['n_dec'] / (1.0 - r['f_lo'])
        rows.append(r)

    w_mm = COL_MM['double']
    panel_w = (w_mm - G['left'] - G['right'] - 3 * G['col_gap']) / 4.0
    h_mm = (G['header'] + sum(G['fold'] + r['h_mm'] + G['xlab'] for r in rows)
            + (len(rows) - 1) * G['row_gap'] + G['bottom'])
    if h_mm > MAX_HEIGHT_MM:
        raise SystemExit(f'ERROR: the grid is {h_mm:.1f} mm tall, over the '
                         f'{MAX_HEIGHT_MM:g} mm type area; lower --decade-mm.')
    fig = plt.figure(figsize=(mm(w_mm), mm(h_mm)), layout='none')

    info = {'stem': stem, 'width_mm': w_mm, 'height_mm': round(h_mm, 2),
            'decade_mm': decade_mm, 'panel_width_mm': round(panel_w, 2),
            'rows': {}}
    band_top = h_mm - G['header']                   # mm from the bottom edge
    for r in rows:
        r['band_top'] = band_top
        axes_top = band_top - G['fold']
        axes_bottom = axes_top - r['h_mm']
        r['axes'] = [fig.add_axes([(G['left'] + c * (panel_w + G['col_gap']))
                                   / w_mm, axes_bottom / h_mm,
                                   panel_w / w_mm, r['h_mm'] / h_mm])
                     for c in range(len(SPECIES_ORDER))]
        band_top = axes_bottom - G['xlab'] - G['row_gap']

    for r in rows:
        metric, lo, hi = r['metric'], r['lo'], r['hi']
        rinfo = {'masked': r['masked'], 'lo': lo, 'hi': hi,
                 'axes_height_mm': round(r['h_mm'], 2), 'n': {}, 'zeros': {},
                 'band_points': {}, 'above': {}, 'fold': {}, 'q': {}}
        r['fold_texts'] = []
        for c, (ax, sp) in enumerate(zip(r['axes'], SPECIES_ORDER)):
            ret, elim = r['per_sp'][sp]
            rinfo['n'][sp] = (len(ret), len(elim))
            if len(ret) < 5 or len(elim) < 5:
                # Keep the row's scales, so the empty panel still lines up.
                mark_axis_break(ax, zero_band_axis(ax, lo, hi))
                ax.set_xlim(-0.6, 1.6)
                if c:
                    ax.tick_params(axis='y', labelleft=False)
                ax.text(0.5, 0.5, 'insufficient data', ha='center',
                        va='center', transform=ax.transAxes,
                        fontsize=FONT['annotation'])
                ax.set_xticks([])
                continue

            data = [as_plotted(ret, metric), as_plotted(elim, metric)]
            # Box statistics from every gene (all values are passed); only the
            # jitter is subsampled.
            bp = ax.boxplot(data, positions=[0, 1], widths=0.55,
                            patch_artist=True, showfliers=False, whis=(5, 95),
                            zorder=2)
            for patch, key in zip(bp['boxes'], ('retained', 'eliminated')):
                # Alpha in the face colour, never set_alpha() on the patch.
                patch.set_facecolor(to_rgba(REGION[key], 0.40))
                patch.set_edgecolor('black')
                patch.set_alpha(None)
            for med in bp['medians']:
                med.set_color('black')
                med.set_zorder(5)
            for art in bp['whiskers'] + bp['caps']:
                art.set_zorder(4)

            for i, (v, key) in enumerate(zip(data, ('retained', 'eliminated'))):
                n_zero = int(at_zero(v, metric).sum())
                rinfo['zeros'][(sp, key)] = n_zero
                rinfo['above'][(sp, key)] = int((v > hi).sum())
                # At most GRID_JITTER_MAX points, sampled so that a band with
                # a count beside it always shows at least one of them.
                v, _ = _jitter_sample(v, metric, GRID_JITTER_MAX, rng)
                sc = ax.scatter(np.full(len(v), i)
                                + rng.normal(0, 0.07, len(v)), v,
                                s=1.6, alpha=0.20, color=REGION[key],
                                linewidths=0, zorder=3, rasterized=True)
                # Checked on the points as drawn, not on the sampler's word.
                drawn0 = int(np.count_nonzero(
                    np.asarray(sc.get_offsets())[:, 1] == 0.0))
                rinfo['band_points'][(sp, key)] = drawn0
                if n_zero and not drawn0:
                    raise SystemExit(
                        f"ERROR: {stem} {r['letter']} {sp} {key}: {n_zero:,} "
                        "at 0 but no point drawn in the band.")

            mark_axis_break(ax, zero_band_axis(ax, lo, hi))
            for i, key in enumerate(('retained', 'eliminated')):
                z = rinfo['zeros'][(sp, key)]
                if z:
                    ax.text(i + 0.30, 0.0, f'{z:,}', ha='left', va='center',
                            zorder=6, fontsize=FONT['annotation'],
                            color=NEUTRAL['ink_light'])
            if metric == 'omega':
                ax.axhline(1.0, color=NEUTRAL['rule'], ls='--',
                           lw=LINEWIDTH['rule'], zorder=1)

            ax.set_xlim(-0.6, 1.6)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f'Retained\n(n = {len(ret):,})',
                                f'Eliminated\n(n = {len(elim):,})'])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            if c:
                # One scale per row: tick marks on every panel, numbers on
                # the first only.
                ax.tick_params(axis='y', labelleft=False)

            sig, fold = lookup(sp, metric, r['masked'])
            rinfo['fold'][sp], rinfo['q'][sp] = fold, sig
            if np.isfinite(fold):
                # Above the axes at a fixed offset, so every row's fold line
                # sits the same distance over its panels.
                r['fold_texts'].append(ax.annotate(
                    f'{stars(sig)}  {fold:.2f}×', xy=(0.5, 1.0),
                    xycoords='axes fraction',
                    xytext=(0.0, GRID_FOLD_OFFSET_PT),
                    textcoords='offset points', ha='center', va='bottom',
                    fontsize=FONT['annotation'], annotation_clip=False))

        first = r['axes'][0]
        first.set_ylabel(METRIC_LABEL[metric])
        # Right edge of every row's y label at the same distance from the
        # figure edge: all first-column axes share one left edge and width.
        first.yaxis.set_label_coords((G['ylabel_x'] - G['left']) / panel_w,
                                     0.5)
        info['rows'][metric] = rinfo

    # Direct label on the rule instead of a key (house rule), placed once
    # every panel of row A is drawn, so it can avoid their data.
    k, clear = place_rule_label(fig, rows[0]['axes'])
    info['rule_label'] = {'panel': None if k is None else SPECIES_ORDER[k],
                          'clear': clear}

    # Species once, as column headers over row A.
    headers = []
    for c, sp in enumerate(SPECIES_ORDER):
        x_mid = G['left'] + c * (panel_w + G['col_gap']) + panel_w / 2.0
        headers.append(fig.text(x_mid / w_mm, (h_mm - 0.4) / h_mm,
                                SPECIES_LABEL[SPECIES_KEY[sp]], ha='center',
                                va='top', fontsize=FONT['panel_title'],
                                fontstyle='italic'))

    letters = []
    if panel_letters:
        for k, r in enumerate(rows):
            # A at the top-left corner, beside the headers; B and C level with
            # the top of their row's fold line. House size, 10 pt bold.
            y = (h_mm - 0.2) if k == 0 else r['band_top']
            r['letter_text'] = fig.text(0.0, y / h_mm, r['letter'],
                                        fontsize=FONT['panel_letter'],
                                        fontweight='bold', ha='left', va='top')
            letters.append(r['letter_text'])

    info['spacing'] = check_grid_spacing(fig, rows, headers, letters)
    audit(fig)
    # Saved with constrained layout off. savefig restores a figure's layout
    # engine after each format, and a figure created without one gets the
    # rcParams default back -- constrained, under apply_style() -- so the PDF
    # and PNG ran a constrained-layout pass (matplotlib 3.11). It cannot move
    # axes placed with add_axes(), but it warns on every save; with the
    # default off, nothing runs and the three files share one geometry.
    with plt.rc_context({'figure.constrained_layout.use': False}):
        paths = save(fig, outdir, stem)
    plt.close(fig)
    for r in rows:
        n_above = sum(info['rows'][r['metric']]['above'].values())
        if n_above:
            logging.info('  %s row %s stops at %g: %s value(s) above it are not '
                         'drawn (legend notes in the summary)', stem,
                         r['letter'], r['hi'], f'{n_above:,}')
    logging.info('  wrote %s (%.1f x %.1f mm)',
                 ', '.join(os.path.basename(p) for p in paths), w_mm, h_mm)
    return info


def grid_legend_notes(args, grid_info, fam_sizes):
    """Lines for the Figure 5 and S8 legends, from the grid as drawn."""
    ds_masked = ds_row_masking(args.mask_ds)
    fam = ', '.join(f'{lab} n = {fam_sizes[m]}' for m, lab in
                    (('omega', 'dN/dS'), ('dN', 'dN'), ('dS', 'dS'))
                    if fam_sizes.get(m))
    L = ['-' * 55,
         'Figure legend notes: the grid (Figure 5 = as_region, S8 = tox_region)',
         '-' * 55,
         'Rows: A branch dN/dS, B branch dN, C branch dS; columns Ascaris, '
         'P. univalens, B. schroederi, T. canis.',
         f'A and B: branches masked where dN/dS > {args.omega_cap:g} or '
         'dN = dS = 0 (the same set, so the same n).',
         'C: ' + ('masked as A and B [CIRCULAR -- --mask-ds was set].'
                  if ds_masked else 'every quartet, unmasked (the cap is a '
                  'filter on dS).'),
         star_caption(args.star_stat, args.fdr_family) + (f' ({fam})'
                                                         if fam else '') + '.',
         'Boxes: quartiles of every gene; whiskers: 5th-95th percentile; '
         f'points: up to {GRID_JITTER_MAX:,} genes per box, sampled at the '
         'same rate in the band and above it, with at least one point in any '
         'band that holds a gene.',
         'Log axes; one decade is the same height in every row. Below the '
         'scale break, the band at 0.',
         'A: the band holds branches at codeml\'s dN/dS floor of 0.0001 (no '
         'nonsynonymous change detected).',
         'B, C: the bands hold values codeml printed as 0.0000 (below 0.00005).',
         'The grey number beside a box is the count in its band.', '']
    for g, info in grid_info.items():
        L.append(f'{g}  ({info["stem"]}, {info["width_mm"]:g} x '
                 f'{info["height_mm"]:g} mm, decade {info["decade_mm"]:g} mm)')
        rl = info.get('rule_label', {})
        if rl and rl['panel'] is None:
            L.append('  "dN/dS = 1" label not drawn: no panel of A has data')
        elif rl:
            L.append(f'  "dN/dS = 1" label on the {rl["panel"]} panel of A'
                     + ('' if rl['clear'] else ' -- it covers data there; '
                        'check the figure'))
        for metric, lab in (('omega', 'A dN/dS'), ('dN', 'B dN'),
                            ('dS', 'C dS')):
            ri = info['rows'].get(metric)
            if not ri:
                continue
            L.append(f'  {lab}: axis {ri["lo"]:g} to {ri["hi"]:g}')
            for sp in SPECIES_ORDER:
                n_r, n_e = ri['n'].get(sp, (0, 0))
                z_r = ri['zeros'].get((sp, 'retained'), 0)
                z_e = ri['zeros'].get((sp, 'eliminated'), 0)
                a_r = ri['above'].get((sp, 'retained'), 0)
                a_e = ri['above'].get((sp, 'eliminated'), 0)
                b_r = ri['band_points'].get((sp, 'retained'), 0)
                b_e = ri['band_points'].get((sp, 'eliminated'), 0)
                fold, q = ri['fold'].get(sp, np.nan), ri['q'].get(sp, np.nan)
                L.append(f'    {sp:<14} n {n_r:,} / {n_e:,}   at 0 {z_r:,} / '
                         f'{z_e:,} (points drawn {b_r:,} / {b_e:,})   above '
                         f'axis {a_r} / {a_e}   fold {fold:.2f}   q {q:.2g} '
                         f'{stars(q)}')
                for key, n, z in (('retained', n_r, z_r),
                                  ('eliminated', n_e, z_e)):
                    if n and z / n > 0.05:
                        L.append(f'      {key} lower whisker (5th percentile) '
                                 f'reaches the band: {z:,} of {n:,} '
                                 f'({100 * z / n:.1f}%) at 0')
        L.append('')
    return L


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def legend_notes(args, fig_info, fam_n):
    """Lines for the Figure 5C / S14C legend: everything the footer used to say,
    plus the counts the panels draw at 0 and the counts clipped above."""
    L = ['-' * 55, 'Figure legend notes (5C and S14C; stars also apply to 5D)',
         '-' * 55,
         f'dN: capped at branch dN/dS <= {args.omega_cap:g}, on that branch '
         'only.',
         'dS: ' + ('capped as dN [CIRCULAR -- --mask-ds was set].'
                   if ds_row_masking(args.mask_ds) else
                   'uncapped (the cap is a filter on dS).'),
         star_caption(args.star_stat, args.fdr_family, fam_n) + '.',
         'Boxes: quartiles; whiskers: 5th-95th percentile; jitter: up to '
         '1,500 genes per box.']
    if args.yscale == 'log':
        # "Printed as 0.0000", not "dN/dS at the floor": the 5A band counts
        # the floor, and on a long branch the two differ (branch_lib).
        L.append('Genes whose dN or dS codeml printed as 0.0000 are drawn in '
                 'the band below')
        L.append('the scale break; the grey number beside each box is that '
                 'count (retained / eliminated):')
    for g, info in fig_info.items():
        for metric in ('dN', 'dS'):
            lo, hi = info['limits'][metric]
            z = info['zeros']
            a = info['above']
            L.append(f'  {g:<11} {metric}  ' + '   '.join(
                f"{sp} {z.get((metric, sp, 'retained'), 0):,} / "
                f"{z.get((metric, sp, 'eliminated'), 0):,}"
                for sp in SPECIES_ORDER))
            n_above = sum(n for (m, _, _), n in a.items() if m == metric)
            L.append(f'  {"":<11} {metric}  axis {lo:g} to {hi:g}; '
                     f'{n_above:,} value(s) above it clipped'
                     + (' (' + ', '.join(
                         f"{sp} {a.get((metric, sp, 'retained'), 0)}/"
                         f"{a.get((metric, sp, 'eliminated'), 0)}"
                         for sp in SPECIES_ORDER) + ')' if n_above else ''))
    return L + ['']


def write_summary(path, stats_all, args, n_genes, groupings, fig_info=None,
                  fam_n=None, grid_info=None, fam_sizes=None):
    """Human-readable summary, with the multiple-testing block for dS.

    `grid_info` (from plot_branch_grid) adds the legend notes of the grid;
    `fig_info` (from plot_dn_ds_panels) those of the old separate panels.
    """
    ds_masked = ds_row_masking(args.mask_ds)

    with open(path, 'w') as fh:
        fh.write('Per-branch dN / dS decomposition\n')
        fh.write('=' * 55 + '\n\n')
        fh.write(f'Input: {args.annotated}\n')
        fh.write(f'Genes: {n_genes:,}\n')
        fh.write(f'Per-branch dN/dS cap: {args.omega_cap:g} '
                 '(applied to that branch only)\n')
        fh.write(f'dS rows used for plots and claims: masked={ds_masked}'
                 f"{'  [CIRCULAR -- --mask-ds was set]' if ds_masked else ''}\n")
        if np.isfinite(args.ds_max):
            fh.write(f'dS upper bound: {args.ds_max:g}\n')
        fh.write('\nInterpretation: dN fold > 1 with dS fold ~ 1 indicates '
                 'relaxed purifying\nselection rather than a mutation-rate '
                 'difference.\n\n')
        fh.write('MASKING: dN and dN/dS use the capped subset (matches Figure '
                 '5A). dS uses\nALL genes -- the cap is a filter on dS, so '
                 'testing dS after applying it is\ncircular. Quote the '
                 'masked=False rows for any dS claim.\n\n')
        fh.write('CAUTION: the Toxocara terminal branch is the longest in the '
                 'tree and its dS\napproaches saturation, so a null dS result '
                 "there is weaker evidence of 'no\ndifference' than the same "
                 'result on the Ascaris or Parascaris branches. The\nconverse '
                 'also holds: a positive dS result there is weak evidence of a '
                 'real\nsynonymous-rate difference, because alignment and '
                 'gene-model noise inflate dN\nand dS together.\n\n')
        if grid_info:
            fh.write('\n'.join(grid_legend_notes(args, grid_info,
                                                 fam_sizes or {})) + '\n')
            fh.write('Grid layout, measured on the drawn figure (mm): '
                     + '; '.join(f'{g}: gaps {gi["spacing"]["gaps_mm"]}, '
                                 f'decade {gi["spacing"]["decade_mm"]}'
                                 for g, gi in grid_info.items()) + '\n\n')
        if fig_info:
            fh.write('\n'.join(legend_notes(args, fig_info, fam_n)) + '\n')

        fh.write('-' * 55 + '\n')
        fh.write('fold(dN) / fold(dS) -- dN/dS elevation not attributable to '
                 'dS\n')
        fh.write('-' * 55 + '\n')
        fh.write(f"{'grouping':<12} {'branch':<15} {'fold dN':>8} "
                 f"{'fold dS':>8} {'ratio':>8}\n")
        for g in groupings:
            for sp in SPECIES_ORDER:
                sel = stats_all[(stats_all['grouping'] == g)
                                & (stats_all['branch'] == sp)]
                dn = sel[(sel['metric'] == 'dN') & (sel['masked'])]
                ds = sel[(sel['metric'] == 'dS')
                         & (sel['masked'] == ds_masked)]
                if not len(dn) or not len(ds):
                    continue
                a = float(dn['fold_change'].iloc[0])
                b = float(ds['fold_change'].iloc[0])
                ratio = (a / b if np.isfinite(a) and np.isfinite(b) and b > 0
                         else np.nan)
                fh.write(f'{g:<12} {sp:<15} {a:>8.3f} {b:>8.3f} '
                         f'{ratio:>8.3f}\n')
        fh.write('\n')

        pooled = stats_all[(stats_all['metric'] == 'dS')
                           & (~stats_all['masked'])].copy()
        if len(pooled):
            pooled['q_pooled'] = benjamini_hochberg(
                pooled['p_value'].to_numpy(dtype=float))
            n_sig = int((pooled['q_pooled'] < 0.05).sum())
            fh.write('-' * 55 + '\n')
            fh.write('Multiple testing across ALL unmasked dS tests '
                     f'(n = {len(pooled)}, both groupings)\n')
            fh.write('This is the correct family if both region labellings are '
                     'reported.\n')
            fh.write('-' * 55 + '\n')
            fh.write(f"{'grouping':<12} {'branch':<15} {'fold':>7} "
                     f"{'p':>11} {'q_pooled':>11}\n")
            for _, r in pooled.iterrows():
                fh.write(f"{r['grouping']:<12} {r['branch']:<15} "
                         f"{r['fold_change']:>7.3f} {r['p_value']:>11.3g} "
                         f"{r['q_pooled']:>11.3g}\n")
            fh.write(f'\nSignificant at pooled BH q < 0.05: {n_sig} of '
                     f"{len(pooled)} (min q = {pooled['q_pooled'].min():.3g})\n")
            if n_sig == 0:
                fh.write('=> No synonymous-rate difference survives correction '
                         'for multiple testing\n   on any branch, while every '
                         'dN difference is significant by many orders\n   of '
                         'magnitude.\n')
            fh.write('\n')

        ds_rows = stats_all[(stats_all['metric'] == 'dS')
                            & (stats_all['masked'] == ds_masked)]
        if len(ds_rows):
            fh.write('-' * 55 + '\n')
            fh.write('SENSITIVITY: dS verdict under each multiple-testing '
                     'family\n')
            fh.write('-' * 55 + '\n')
            fh.write(f'Figures were drawn with --star-stat {args.star_stat} '
                     f'--fdr-family {args.fdr_family}\n\n')
            fh.write(f"{'family':<24} {'n':>4} {'min q':>9} {'# q<0.05':>9}"
                     f"   {'significant branches':<40}\n")
            for fam in list(FDR_FAMILIES) + ['none']:
                tmp = add_star_qvalues(stats_all.copy(), fam)
                sub = tmp[(tmp['metric'] == 'dS')
                          & (tmp['masked'] == ds_masked)]
                if fam == 'none':
                    n_txt = '-'
                else:
                    sizes = sorted(set(
                        tmp.groupby(list(FDR_FAMILIES[fam])).size().tolist()))
                    n_txt = str(sizes[0]) if len(sizes) == 1 else '*'
                hits = sub[sub['q_star'] < 0.05]
                who = ', '.join(f'{r.grouping}/{r.branch}'
                                for _, r in hits.iterrows()) or '(none)'
                fh.write(f"{('none (raw P)' if fam == 'none' else fam):<24} "
                         f"{n_txt:>4} {sub['q_star'].min():>9.4f} "
                         f"{len(hits):>9}   {who:<40}"
                         f"{'  <- NOT VALID' if fam == 'mixed' else ''}\n")
            fh.write("\nWHY 'mixed' IS WRONG: BH sets q = p * m / rank. "
                     'Pooling the eight dN tests\ninto the dS family raises '
                     "every dS test's rank by eight while m only doubles,\nso "
                     'm/rank falls and dS comes out MORE significant. dN and '
                     'dS answer different\nquestions and do not belong in one '
                     'correction family.\n\n')
            # Computed, not fixed text. Until 2026-10-10 this line always said
            # that no dS difference survives correction -- true when written,
            # false after the unrooted-tree rerun, where two Toxocara dS tests
            # do, while the table directly above it said so.
            if args.star_stat == 'q' and args.fdr_family in FDR_FAMILIES:
                tmp = add_star_qvalues(stats_all.copy(), args.fdr_family)
                sub = tmp[(tmp['metric'] == 'dS')
                          & (tmp['masked'] == ds_masked)]
                hits = sub[sub['q_star'] < 0.05]
                if not len(hits):
                    fh.write('What the text can claim (family '
                             f'{args.fdr_family}): no synonymous-rate '
                             'difference survives\ncorrection, on either '
                             'region labelling.\n\n')
                else:
                    fh.write('What the text can claim (family '
                             f'{args.fdr_family}): {len(hits)} of {len(sub)} '
                             'dS tests survive correction:\n')
                    for _, r in hits.iterrows():
                        fh.write(f'  {r.grouping}/{r.branch}: fold '
                                 f'{r.fold_change:.2f}, q = {r.q_star:.2g}\n')
                    fh.write('Do not claim unchanged dS on these branches; '
                             'every other dS test is null.\n\n')

        fh.write('-' * 55 + '\n')
        fh.write('Full stats table. q_value = BH across the four branches '
                 'within each\ngrouping x metric x masking family. q_star = '
                 'the value the figure asterisks\nwere drawn from, under '
                 f'--fdr-family {args.fdr_family}.\n')
        fh.write('-' * 55 + '\n')
        fh.write(stats_all.to_string(index=False))
        fh.write('\n')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Decompose per-branch dN/dS into dN and dS.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', required=True,
                   help='branch_dnds_four_annotated.tsv from branch_dnds_four.py')
    p.add_argument('--omega-cap', type=float, default=5.0,
                   help='Per-branch dN/dS cap; must match branch_dnds_four.py '
                        '(default: 5). Applied to dN/dS and dN only.')
    p.add_argument('--mask-ds', action='store_true',
                   help='Apply the cap to dS as well. CIRCULAR -- the cap is a '
                        'filter on dS. Provided only for side-by-side checking.')
    p.add_argument('--ds-max', type=float, default=float('inf'),
                   help='Drop dS above this as an estimation failure. Default '
                        'inf. Try 2 as a saturation check on Toxocara.')
    p.add_argument('--groupings', default=None,
                   help='Comma-separated region columns. Default: auto-detect '
                        'as_region and tox_region.')
    p.add_argument('--star-stat', choices=['q', 'p'], default='q',
                   help="Statistic the asterisks are drawn from. 'q' (default) "
                        "uses BH-FDR; 'p' uses raw P and puts a '*' on the "
                        'Baylisascaris and Toxocara dS bars that no correction '
                        'supports -- contradicting the manuscript text.')
    p.add_argument('--fdr-family', default='metric',
                   choices=sorted(FDR_FAMILIES) + ['none'],
                   help="Family for --star-stat q. 'metric' (default) corrects "
                        'within each metric across branches and labellings. '
                        "'mixed' pools dN and dS, which is wrong here and makes "
                        'dS look MORE significant; it exists only for the '
                        'sensitivity table.')
    p.add_argument('--yscale', choices=['log', 'linear'], default='log',
                   help='y-axis of Figure 5C. A shared LINEAR limit across a '
                        'row is set by Toxocara, whose branch is an order of '
                        'magnitude longer, and flattens the other three '
                        'panels (default: log)')
    p.add_argument('--figure', choices=['grid', 'legacy', 'both'],
                   default='grid',
                   help="'grid' (default): Figure 5 and S8 whole, three rows "
                        "of four panels (plot_branch_grid). 'legacy': the old "
                        "separate 5C/5D and S14C/S14D panels. 'both': all.")
    p.add_argument('--decade-mm', type=float, default=7.0,
                   help='Height of one log decade, in mm, in every row of the '
                        'grid (default: 7). Sets the figure height: about '
                        '150 mm at 7.')
    p.add_argument('--panel-letters', action=argparse.BooleanOptionalAction,
                   default=None,
                   help='Panel letters: A-C on the grid (default on), C/D on '
                        'the legacy panels (default off).')
    p.add_argument('--seed', type=int, default=0,
                   help='RNG seed for jitter and subsampling (default: 0)')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if not 4.0 <= args.decade_mm <= 12.0:
        # Below 4 mm the dN row (about 4.5 decades) is under the 15 mm axes
        # floor audit() enforces; above 12 the figure leaves the type area.
        p.error(f'--decade-mm {args.decade_mm:g} is outside 4-12 mm')
    if args.figure in ('grid', 'both') and args.yscale != 'log':
        p.error('the grid is drawn on log axes only; use --figure legacy '
                'for --yscale linear')
    if args.outdir is None:
        args.outdir = make_outdir('dnds_decomposition', base='results')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'dnds_decomposition')
    apply_style()

    # Arial is absent on Isaac compute nodes by design (svg.fonttype='none'
    # keeps the SVG referencing it by name for Illustrator). Silence the
    # thousands of findfont warnings so the log stays readable.
    for noisy in ('matplotlib.font_manager', 'matplotlib.texmanager',
                  'fontTools', 'fontTools.subset', 'fontTools.ttLib', 'PIL'):
        logging.getLogger(noisy).setLevel(logging.ERROR)

    logging.info('Per-branch dN / dS decomposition')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))

    if args.ds_max <= 0:
        sys.exit('ERROR: --ds-max must be positive.')
    if args.omega_cap <= 0:
        sys.exit('ERROR: --omega-cap must be positive.')

    groupings = ([g.strip() for g in args.groupings.split(',')]
                 if args.groupings else None)
    df, groupings = load_annotated(args.annotated, groupings)

    ds_masked = ds_row_masking(args.mask_ds)
    logging.info('dS masking: cap %s -> plots and claims use masked=%s rows',
                 'APPLIED (circular)' if ds_masked else 'NOT applied', ds_masked)

    # One generator per figure type, each seeded with --seed: the legacy panels
    # draw exactly the jitter they drew before the grid existed, whichever
    # figures are requested.
    rng = np.random.default_rng(args.seed)
    rng_grid = np.random.default_rng(args.seed)

    # Stats for EVERY grouping are built before any figure is drawn: the
    # default FDR family pools across region labellings, so computing q_star
    # one grouping at a time would silently give the 'metric-grouping' family.
    all_stats = []
    for g in groupings:
        logging.info('Grouping %s: %s', g, df[g].value_counts().to_dict())
        all_stats.append(build_stats(df, g, args.omega_cap, ds_max=args.ds_max))
    stats_all = pd.concat(all_stats, ignore_index=True)
    add_star_qvalues(stats_all, args.fdr_family)

    fam_n = None
    if args.star_stat == 'q' and args.fdr_family != 'none':
        ds_fam = stats_all[(stats_all['metric'] == 'dS')
                           & (stats_all['masked'] == ds_masked)]
        sizes = ds_fam.groupby(list(FDR_FAMILIES[args.fdr_family])).size()
        fam_n = int(sizes.iloc[0]) if len(set(sizes.tolist())) == 1 else None
        logging.info('Asterisks: %s (family=%s); dS family size(s) %s',
                     args.star_stat, args.fdr_family, sorted(set(sizes.tolist())))

    # Family size behind each grid row's asterisks, for its legend line.
    fam_sizes = {}
    if args.star_stat == 'q' and args.fdr_family != 'none':
        for metric, mk in (('omega', True), ('dN', True), ('dS', ds_masked)):
            sel = stats_all[(stats_all['metric'] == metric)
                            & (stats_all['masked'] == mk)]
            sizes = set(sel.groupby(list(FDR_FAMILIES[args.fdr_family]))
                        .size().tolist())
            if len(sizes) == 1:
                fam_sizes[metric] = sizes.pop()

    letters_grid = True if args.panel_letters is None else args.panel_letters
    letters_legacy = bool(args.panel_letters)
    fig_info, grid_info = {}, {}
    for g in groupings:
        stats = stats_all[stats_all['grouping'] == g]
        for sp in SPECIES_ORDER:
            for metric, mk in (('dN', True), ('dS', ds_masked)):
                sel = stats[(stats.branch == sp) & (stats.metric == metric)
                            & (stats.masked == mk)]
                if not len(sel):
                    continue
                r = sel.iloc[0]
                logging.info('  %s %-14s %s (masked=%s): ret %.4f (n=%s) '
                             'elim %.4f (n=%s) fold %.2f %s q_star=%.3g -> %s',
                             g, sp, metric, mk, r.ret_median, f'{int(r.ret_n):,}',
                             r.elim_median, f'{int(r.elim_n):,}', r.fold_change,
                             fmt_p(r.p_value).replace('$', ''), r.q_star,
                             stars(r[star_column(args.star_stat)]))

        # File names keep the frozen-map numbers (S14 is now S8); the
        # assembly manifests map them to the published numbers.
        tag = 'fig5' if g == 'as_region' else 'figS14'
        if args.figure in ('grid', 'both'):
            grid_info[g] = plot_branch_grid(
                df, g, args.outdir, f'{tag}_branch_grid_{g}', args.omega_cap,
                args.mask_ds, args.ds_max, rng_grid, stats, args.star_stat,
                decade_mm=args.decade_mm, panel_letters=letters_grid)
        if args.figure in ('legacy', 'both'):
            fig_info[g] = plot_dn_ds_panels(
                df, g, args.outdir, f'{tag}C_dn_ds_{g}', args.omega_cap,
                args.mask_ds, args.ds_max, rng, stats, args.star_stat,
                args.fdr_family, fam_n, yscale=args.yscale,
                panel_letters=letters_legacy)
            plot_fold_changes(stats_all, g, args.outdir,
                              f'{tag}D_fold_change_{g}', args.mask_ds,
                              args.star_stat, args.fdr_family, fam_n,
                              panel_letters=letters_legacy)

    stats_path = os.path.join(args.outdir, 'dnds_decomposition_stats.tsv')
    stats_all.to_csv(stats_path, sep='\t', index=False, float_format='%.6g')
    logging.info('Stats: %s', stats_path)

    write_summary(os.path.join(args.outdir, 'decomposition_summary.txt'),
                  stats_all, args, len(df), groupings, fig_info=fig_info,
                  fam_n=fam_n, grid_info=grid_info, fam_sizes=fam_sizes)

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(),
                   'python': platform.python_version(),
                   'argv': sys.argv, 'n_genes': int(len(df)),
                   'groupings': groupings,
                   'figure': args.figure,
                   'grid': {g: {k: v for k, v in gi.items() if k != 'rows'}
                            for g, gi in grid_info.items()},
                   'inputs': {'annotated': fingerprint(args.annotated)}},
                  fh, indent=2)
    logging.info('Done.')


if __name__ == '__main__':
    main()
