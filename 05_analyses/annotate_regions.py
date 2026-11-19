#!/usr/bin/env python3
"""
annotate_regions.py - region, tract type and CBR distance for every ortholog pair.

Upstream. Writes `kaks_annotated.tsv`, which almost every downstream script reads,
plus supplementary fig. S7 (S5 before the 2026-10-12 renumbering) and the table
behind table S1. Replaces `annotate_kaks_regions.py`.

Writes:

    kaks_annotated.tsv              the annotated rate table -- the real product
    S5A_elim_subtypes.{svg,pdf,png} fig. S7A, terminal vs internal eliminated
    S5B_region_x_proximity.*        fig. S7B, region x CBR proximity, four
                                    classes. Both keep their old S5 stems so the
                                    assembly manifest and figure_locations stay
                                    valid; the assembly copies them to FigS07_*.
    region_boxplots.*               retained vs eliminated, diagnostic
    consensus_boxplots.*            pairs where both species agree, diagnostic
    annotate_regions_summary.txt    counts, medians, tests
    annotate_regions_stats.tsv      every test, with BH q-values
    provenance.json                 inputs, arguments and versions


Tests
-----
Two-sided Mann-Whitney, BH-corrected within each family (per metric):

    region              retained vs eliminated (species-1 calls)
    consensus           the same, on pairs both species call alike
    elim_subtype        terminal vs internal eliminated tracts (S7A)
    proximity_pooled    CBR-proximal vs interior, all genes pooled
    region_x_proximity  the four region x proximity classes (S7B)

`proximity_pooled` is the comparison that `region_x_proximity` takes apart:
CBR-proximal genes are higher pooled, but not within either region. Results
quotes its P, so it is computed here rather than by a one-off command. It is a
family of one test per metric, so its q equals its P and it moves no other q.


Sign convention
---------------
`signed_cbr_dist` is **positive on the eliminated side**, negative on the
retained side. This is new: the old script wrote only an unsigned `min_cbr_dist`,
so every downstream script that wanted a signed axis re-derived the sign itself,
and they did not all choose the same direction — Figure 4B was drawn
eliminated-left while Figure 4D was fitted eliminated-right. Writing the signed
distance once, here, is what stops that recurring.

Fitted slopes and boundary steps in `cbr_changepoint_relax.py` do not change: they
were already computed in the eliminated-positive frame.


What counts as a CBR
--------------------
Delegated to `region_lib.derive_cbrs`, which is also what `gene_density_cbr.py`
uses. Panels 4A, 4B and 4D therefore share one CBR set. Two things it does that
the old code did not:

  - Eliminated tracts are **merged** before boundaries are taken, so the junction
    between two abutting tracts is not a CBR sitting inside eliminated sequence.
  - Boundaries at a sequence end are **dropped** (`--terminal-tol`): a terminal
    tract runs to the telomere, and its distal edge is where the assembly stops.

On *Ascaris* this takes the CBR count from every BED edge down to the boundaries
that are actually breakage regions, and it changes which genes are called
CBR-proximal.


One locus per gene
------------------
A miniprot query can align more than once. The old code took the first record in
the file; this takes the best (lowest `Rank`, then highest `Identity`, then
longest span). A secondary alignment can sit on a different chromosome from the
primary, so for those genes the choice decides the region call outright. The
summary reports how many calls this changes.


Changes 2026-10-12 (figure review of supplementary fig. S7)
-----------------------------------------------------------
- Brackets sit above the highest point drawn under them. They were placed at
  1.06 x the 99th percentile while points up to 1.28 x were drawn, so they
  crossed the upper tail (in S7B, the one eliminated gene above 1).
- Only points up to the cap (the --cap-pct percentile of the highest group) are
  drawn, and the points above it are counted per panel and group in the
  summary. The summary used to say that points outside the central 99% were
  omitted while drawing points up to 1.28 x the 99th percentile and clipping
  the rest uncounted.
- Hatching means CBR-proximal and nothing else. S7A no longer hatches the
  internal tracts: both boxes are eliminated, in the eliminated hue, and the
  tick labels name the subtype.
- S7B has no title. The definition of CBR-proximal goes in the caption, and
  the summary prints it.
No statistic changes: every test and kaks_annotated.tsv are computed before any
figure is drawn. The random stream only places jitter.

Usage
-----
    python annotate_regions.py \\
      --kaks kaks_results/kaks_filtered.tsv \\
      --gff1 As_miniprot.gff --gff2 Pu_miniprot.gff \\
      --retained1 AG.v50.retained_loose.bed \\
      --eliminated1 AG.v50.eliminated_strict.bed \\
      --retained2 germ_to_soma_mapping.bed \\
      --eliminated2 pu_v1.13_eliminated.bed \\
      --chrom-sizes1 AG.v50.fa.fai \\
      --species1 Ascaris --species2 Parascaris \\
      --outdir results/annotate_regions_20260817

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
from scipy.stats import mannwhitneyu

from paper_colors import (REGION, NEUTRAL, PROXIMITY_HATCH, apply_style, audit,
                          save, figure, panel_letter, mm)
from region_lib import (parse_gff_models, select_primary_alignment, parse_bed,
                        merge_intervals, load_chrom_sizes, derive_cbrs,
                        annotate_positions, label_arrays, fmt_p, stars,
                        benjamini_hochberg, fingerprint, git_commit,
                        setup_logging)

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Drawn left-to-right. Region sets the hue, proximity sets the texture, so the
#: four classes of S7B read as two regions x two proximities rather than as four
#: unrelated categories.
PROXIMITY_ORDER = [('retained', 'interior'), ('retained', 'cbr_proximal'),
                   ('eliminated', 'interior'), ('eliminated', 'cbr_proximal')]

#: Terminal and internal eliminated tracts are *both* eliminated, so both take the
#: eliminated hue. Drawing one blue would say "retained" to any reader who has
#: learned the convention from Figure 2. The subtype is carried by the tick
#: label alone: hatching means CBR-proximal in S7B, and one texture must mean
#: one thing across the figure.
SUBTYPE_LABEL = {'end': 'Terminal', 'internal': 'Internal'}


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

def annotate_species(gff, retained_bed, eliminated_bed, chrom_sizes_path,
                     cbr_window, terminal_tol, label):
    """Region calls and CBR distances for one species. Returns (frame, diagnostics)."""
    logging.info('--- %s ---', label)
    models = parse_gff_models(gff)
    genes, align_diag = select_primary_alignment(models)
    genes = genes.assign(midpoint=(genes['start'] + genes['end']) // 2)

    retained = parse_bed(retained_bed)
    eliminated = parse_bed(eliminated_bed, keep_label=True)
    ret_merged = merge_intervals(retained)
    elim_merged = merge_intervals(eliminated)
    elim_labels = label_arrays(eliminated, elim_merged)

    sizes, inferred = load_chrom_sizes(chrom_sizes_path, genes, eliminated)
    cbrs, cbr_diag = derive_cbrs(elim_merged, sizes, terminal_tol)
    cbr_diag['chrom_sizes_inferred'] = inferred

    annotated, counts = annotate_positions(genes, ret_merged, elim_merged,
                                           elim_labels, cbrs, cbr_window)
    diag = {'alignment': align_diag, 'cbrs': cbr_diag, 'region_counts': counts,
            'has_subtype_labels': bool(elim_labels),
            'n_cbr_proximal': int(annotated['cbr_proximal'].sum())}
    return annotated, diag


def merge_with_rates(kaks, annot1, annot2, gene1_col, gene2_col):
    """Join both species' annotations onto the rate table.

    `consensus_region` is the region only where the two species agree; where they
    disagree it is `discordant`, which is a real category and not a missing value
    — a gene retained in one species and eliminated in the other is informative,
    and collapsing it into NaN loses that.
    """
    def side(annot, gene_col, n):
        cols = ['gene_id', 'chrom', 'start', 'end', 'midpoint', 'region',
                'elim_type', 'cbr_proximal', 'min_cbr_dist', 'signed_cbr_dist']
        out = annot[cols].copy()
        out.columns = [gene_col] + [f'{c}{n}' for c in cols[1:]]
        return out

    merged = kaks.merge(side(annot1, gene1_col, 1), on=gene1_col, how='left')
    merged = merged.merge(side(annot2, gene2_col, 2), on=gene2_col, how='left')

    r1, r2 = merged['region1'], merged['region2']
    known = r1.isin(['retained', 'eliminated']) & r2.isin(['retained', 'eliminated'])
    merged['consensus_region'] = np.where(
        ~known, 'unannotated', np.where(r1 == r2, r1, 'discordant'))

    merged['cbr_proximal_either'] = (
        merged['cbr_proximal1'].fillna(False).astype(bool) |
        merged['cbr_proximal2'].fillna(False).astype(bool))

    merged['proximity1'] = np.where(
        merged['cbr_proximal1'].fillna(False).astype(bool),
        'cbr_proximal', 'interior')
    return merged


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def pairwise_tests(df, group_col, groups, value_col, family):
    """Mann-Whitney across every pair of groups, as tidy rows.

    Two-sided throughout: the eliminated-higher direction is the finding, not an
    assumption, and a one-sided test would build it into the method.
    """
    rows = []
    for i, g1 in enumerate(groups):
        for g2 in groups[i + 1:]:
            v1 = df.loc[df[group_col] == g1, value_col].dropna()
            v2 = df.loc[df[group_col] == g2, value_col].dropna()
            if len(v1) < 5 or len(v2) < 5:
                continue
            u, p = mannwhitneyu(v1, v2, alternative='two-sided')
            rows.append({'family': family, 'metric': value_col,
                         'group1': str(g1), 'group2': str(g2),
                         'n1': len(v1), 'n2': len(v2),
                         'median1': v1.median(), 'median2': v2.median(),
                         'fold': (v2.median() / v1.median()
                                  if v1.median() else np.nan),
                         'U': u, 'p': p})
    return rows


def group_summary(df, group_col, groups, value_col):
    rows = []
    for g in groups:
        v = df.loc[df[group_col] == g, value_col].dropna()
        rows.append({'group': str(g), 'n': len(v), 'median': v.median(),
                     'q1': v.quantile(0.25), 'q3': v.quantile(0.75),
                     'mean': v.mean()})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def box_panel(ax, data, colours, hatches, labels, rng, ylabel,
              box_alpha=0.40, cap_pct=99.0, jitter_max=1500,
              neutrality=None):
    """Box over jitter, house style.

    Ordering and opacity are the two things that go wrong here, so both are
    explicit. The box fill is translucent (`box_alpha`) but the *edge is not*:
    `patch.set_alpha()` dims the outline along with the fill and turns a black box
    grey, so the opacity goes into the face colour and `set_alpha` is left unset.
    Draw order is fill (2) under jitter (3) under whiskers and caps (4) under the
    median (5), so the median stays readable over the densest part of the cloud.

    `cap_pct` bounds only what is *drawn*. Every test uses the complete
    distribution. Without it a handful of high-omega outliers set the y-range and
    squash every box into the bottom tenth of the panel. The cap is the cap_pct
    percentile of the highest group; points above it are not drawn, and are
    counted rather than clipped silently by the axes.

    Returns {'cap', 'lo', 'top' (highest drawn point per group), 'hidden'
    (points above the cap per group), 'n' (per group)}. The y-range is set by
    set_box_ylim(), after the brackets are placed.
    """
    n = len(data)
    pos = np.arange(n)
    bp = ax.boxplot(data, positions=pos, widths=0.55, patch_artist=True,
                    showfliers=False, whis=(5, 95), zorder=2)

    for patch, colour, hatch in zip(bp['boxes'], colours, hatches):
        patch.set_facecolor(to_rgba(colour, box_alpha))
        patch.set_edgecolor('black')
        patch.set_linewidth(0.8)
        patch.set_hatch(hatch)
        patch.set_alpha(None)
        patch.set_zorder(2)
    for part in ('whiskers', 'caps'):
        for artist in bp[part]:
            artist.set_color('black')
            artist.set_linewidth(0.8)
            artist.set_zorder(4)
    for med in bp['medians']:
        med.set_color('black')
        med.set_linewidth(1.4)
        med.set_zorder(5)

    cap = max(np.percentile(d, cap_pct) for d in data if len(d))
    lo = min(0.0, min(d.min() for d in data if len(d)))
    shown = [d[d <= cap] for d in data]
    for i, (d, colour) in enumerate(zip(shown, colours)):
        show = d if len(d) <= jitter_max else rng.choice(d, jitter_max,
                                                         replace=False)
        ax.scatter(np.full(len(show), i) + rng.normal(0, 0.055, len(show)),
                   show, s=1.6, linewidths=0, color=colour, alpha=0.35,
                   zorder=3, rasterized=True)

    # The omega = 1 rule is drawn only when it falls inside the range the data
    # already occupy. Forcing the axis out to 1.0 for a rule when the boxes sit
    # near 0.1 is what made the old panels unreadable.
    if neutrality is not None and lo < neutrality <= cap:
        ax.axhline(neutrality, color=NEUTRAL['rule'], lw=0.8, ls=(0, (4, 2)),
                   zorder=1)

    ax.set_xticks(pos)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.6, n - 0.4)
    # Whisker caps are drawn too; the 95th percentile never exceeds the
    # highest drawn point, but take the maximum anyway.
    whisk = [float(np.percentile(d, 95)) if len(d) else lo for d in data]
    return {'cap': float(cap), 'lo': float(lo),
            'top': [max(float(s.max()) if len(s) else lo, w)
                    for s, w in zip(shown, whisk)],
            'hidden': [int((d > cap).sum()) for d in data],
            'n': [int(len(d)) for d in data]}


def set_box_ylim(ax, info, brackets):
    """Place brackets above the highest point drawn under each, then set ylim.

    `brackets` is [(i, j, text)]. Each bracket sits 4% of the drawn range above
    the highest point (or whisker) of the groups it spans, so it never crosses
    the data; brackets over disjoint groups need not stack. The upper limit
    leaves room for the tallest bracket and its label.
    """
    lo, cap = info['lo'], info['cap']
    span = cap - lo if cap > lo else 1.0
    ys = [max(info['top'][min(i, j):max(i, j) + 1]) + 0.04 * span
          for i, j, _ in brackets]
    top = max([y + 0.13 * span for y in ys] + [cap + 0.05 * span])
    ax.set_ylim(lo, top)
    for (i, j, text), y in zip(brackets, ys):
        add_bracket(ax, i, j, y, text)
    return {'ylim': [float(lo), float(top)], 'bracket_y': [float(y) for y in ys]}


def add_bracket(ax, i, j, y, text):
    """Significance bracket. `ns` is drawn like any other result — the panel has
    to state the negative, not leave the reader to infer it from a missing mark."""
    h = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.022
    ax.plot([i, i, j, j], [y, y + h, y + h, y], color='black', lw=0.7,
            clip_on=False, zorder=6)
    ax.text((i + j) / 2, y + h * 1.15, text, ha='center', va='bottom',
            fontsize=plt.rcParams['xtick.labelsize'], zorder=6)


def plot_subtypes(df, metric, ylabel, outdir, stem, rng, tests,
                  width='single', height_mm=62.0, cap_pct=99.0,
                  box_alpha=0.40, panel_letters=False):
    """Fig. S7A - terminal versus internal eliminated tracts.

    Both groups are eliminated, so both take the eliminated hue and neither is
    hatched; the tick labels name the subtype. The previous version drew one blue
    and one red, which reads as retained versus eliminated to anyone who learned
    the convention earlier in the paper, and then hatched the internal tracts,
    while S7B uses the same hatch for CBR-proximal.

    Returns {'groups', 'paths', 'cap', 'hidden', 'n', ...} for the summary.
    """
    order = [g for g in ('end', 'internal') if (df['elim_type1'] == g).sum() >= 5]
    if len(order) < 2:
        logging.warning('  S7A: fewer than two eliminated subtypes with n >= 5')
        return None
    data = [df.loc[df['elim_type1'] == g, metric].dropna().to_numpy() for g in order]

    fig, ax = _axes(width, height_mm)
    info = box_panel(ax, data,
                     [REGION['eliminated']] * len(order),
                     [''] * len(order),
                     [f'{SUBTYPE_LABEL[g]}\n(n = {len(d):,})'
                      for g, d in zip(order, data)],
                     rng, ylabel, box_alpha=box_alpha, cap_pct=cap_pct,
                     neutrality=1.0)
    t = next((r for r in tests if r['family'] == 'elim_subtype'
              and r['metric'] == metric), None)
    brackets = [(0, 1, f"{stars(t['q'])}  {fmt_p(t['p'])}")] if t else []
    info.update(set_box_ylim(ax, info, brackets))
    if panel_letters:
        panel_letter(ax, 'A')
    info['groups'] = [SUBTYPE_LABEL[g] for g in order]
    info['paths'] = _finish(fig, ax, outdir, stem)
    return info


def plot_region_x_proximity(df, metric, ylabel, outdir, stem, rng, tests,
                            cbr_window_kb, width='double', height_mm=62.0,
                            cap_pct=99.0, box_alpha=0.40, panel_letters=False):
    """Fig. S7B - region crossed with CBR proximity, four classes.

    Hue is region and texture is proximity, so the panel decomposes into two
    comparisons a reader can hold at once. Four unrelated categorical colours,
    which is what this was, forced the reader to consult the legend for every box
    and broke the region convention at the same time.

    Tick labels are two lines with the count on the second, replacing the
    underscored single-line labels that physically overlapped at journal width.
    The definition of CBR-proximal (`cbr_window_kb`) goes in the caption, not in
    a title; it is returned for the summary.

    Returns {'groups', 'paths', 'cap', 'hidden', 'n', ...} for the summary.
    """
    groups, data, colours, hatches, labels = [], [], [], [], []
    for region, prox in PROXIMITY_ORDER:
        sel = (df['region1'] == region) & (df['proximity1'] == prox)
        v = df.loc[sel, metric].dropna().to_numpy()
        if len(v) < 5:
            continue
        key = f'{region}|{prox}'
        df.loc[sel, '_rxp'] = key
        groups.append(key)
        data.append(v)
        colours.append(REGION[region])
        hatches.append(PROXIMITY_HATCH[prox])
        labels.append(f"{region.capitalize()}, "
                      f"{'CBR-proximal' if prox == 'cbr_proximal' else 'interior'}\n"
                      f"(n = {len(v):,})")
    if len(groups) < 2:
        logging.warning('  S7B: fewer than two classes with n >= 5')
        return None

    fig, ax = _axes(width, height_mm)
    info = box_panel(ax, data, colours, hatches, labels, rng, ylabel,
                     box_alpha=box_alpha, cap_pct=cap_pct, neutrality=1.0)

    # Only the two within-region proximity contrasts are bracketed. The
    # between-region contrast is Figure 2A's result and does not need restating
    # here; bracketing all six pairs would bury the two comparisons the panel is
    # actually for.
    # The two pairs span disjoint groups, so their brackets need not stack.
    brackets = []
    for region in ('retained', 'eliminated'):
        a, b = f'{region}|interior', f'{region}|cbr_proximal'
        if a in groups and b in groups:
            t = next((r for r in tests
                      if r['family'] == 'region_x_proximity'
                      and r['metric'] == metric
                      and {r['group1'], r['group2']} == {a, b}), None)
            if t:
                brackets.append((groups.index(a), groups.index(b),
                                 f"{stars(t['q'])}  {fmt_p(t['p'])}"))
    info.update(set_box_ylim(ax, info, brackets))
    if panel_letters:
        panel_letter(ax, 'B')
    info['groups'] = [lab.split('\n')[0] for lab in labels]
    info['cbr_window_kb'] = float(cbr_window_kb)
    info['paths'] = _finish(fig, ax, outdir, stem)
    return info


def plot_simple(df, group_col, order, metric, ylabel, outdir, stem, rng, tests,
                family, width='single', height_mm=62.0, cap_pct=99.0,
                box_alpha=0.40):
    """Two-group diagnostic panel: retained vs eliminated, or the consensus set."""
    order = [g for g in order if (df[group_col] == g).sum() >= 5]
    if len(order) < 2:
        return None
    data = [df.loc[df[group_col] == g, metric].dropna().to_numpy() for g in order]
    fig, ax = _axes(width, height_mm)
    info = box_panel(ax, data, [REGION[g] for g in order], [''] * len(order),
                     [f'{g.capitalize()}\n(n = {len(d):,})'
                      for g, d in zip(order, data)],
                     rng, ylabel, box_alpha=box_alpha, cap_pct=cap_pct,
                     neutrality=1.0)
    t = next((r for r in tests if r['family'] == family
              and r['metric'] == metric), None)
    brackets = [(0, 1, f"{stars(t['q'])}  {fmt_p(t['p'])}")] if t else []
    info.update(set_box_ylim(ax, info, brackets))
    info['groups'] = [g.capitalize() for g in order]
    info['paths'] = _finish(fig, ax, outdir, stem)
    return info


def _axes(width, height_mm):
    """Explicit margins in mm.

    `constrained_layout` treats the axes as the flexible element: give it
    anything it cannot fit and it shrinks the panel toward zero and emits only a
    warning. Fixed margins fail loudly instead, and `audit()` checks the result.
    """
    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax)
    w_in, h_in = fig.get_size_inches()
    left, right, bottom, top = mm(14.0), mm(3.0), mm(14.0), mm(5.0)
    ax = fig.add_axes([left / w_in, bottom / h_in,
                       (w_in - left - right) / w_in,
                       (h_in - bottom - top) / h_in])
    return fig, ax


def _finish(fig, ax, outdir, stem):
    for problem in audit(fig, verbose=False):
        logging.warning('  audit: %s', problem)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Annotate a pairwise rate table by region, tract type and '
                    'CBR distance; draw supplementary fig. S7 (written to '
                    'S5A_* and S5B_*).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--kaks', required=True, help='Filtered pairwise rate TSV')
    p.add_argument('--gff1', required=True)
    p.add_argument('--gff2', required=True)
    p.add_argument('--retained1', required=True)
    p.add_argument('--eliminated1', required=True,
                   help='4th column, if present, is the tract type (end/internal)')
    p.add_argument('--retained2', required=True)
    p.add_argument('--eliminated2', required=True)
    p.add_argument('--chrom-sizes1', default=None,
                   help='.fai or two-column TSV. Needed to tell a terminal tract '
                        'boundary from a breakage region; without it a telomeric '
                        'edge is counted as a CBR.')
    p.add_argument('--chrom-sizes2', default=None)
    p.add_argument('--cbr-window', type=float, default=100000,
                   help='Distance within which a gene is CBR-proximal (bp)')
    p.add_argument('--terminal-tol', type=int, default=1000,
                   help='Tract boundaries within this many bp of a sequence end '
                        'are telomeric, not breakage regions (default: 1000). '
                        'Must match gene_density_cbr.py or Figure 4 panels will '
                        'use different CBR sets.')
    p.add_argument('--gene1-col', default='gene1')
    p.add_argument('--gene2-col', default='gene2')
    p.add_argument('--metric', default='Ka_Ks',
                   help='Primary metric for the panels (default: Ka_Ks)')
    p.add_argument('--metric-label', default='dN/dS')
    p.add_argument('--extra-metrics', nargs='*', default=['Ka', 'Ks'],
                   help='Also tested and tabulated, not drawn (default: Ka Ks)')
    p.add_argument('--species1', default='Ascaris')
    p.add_argument('--species2', default='Parascaris')
    p.add_argument('--cap-pct', type=float, default=99.0,
                   help='Percentile bounding what is drawn; every test uses the '
                        'complete distribution')
    p.add_argument('--box-alpha', type=float, default=0.40)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--panel-letters', action='store_true',
                   help='Off by default: panel letters are added at assembly')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'annotate_regions_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'annotate_regions')
    logging.info('annotate_regions.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()
    rng = np.random.default_rng(args.seed)

    kaks = pd.read_csv(args.kaks, sep='\t')
    logging.info('Loaded %s rate rows from %s', f'{len(kaks):,}',
                 os.path.basename(args.kaks))
    for col in (args.gene1_col, args.gene2_col, args.metric):
        if col not in kaks.columns:
            sys.exit(f'ERROR: column {col!r} not in {args.kaks}\n'
                     f'  columns: {list(kaks.columns)}')

    annot1, diag1 = annotate_species(args.gff1, args.retained1, args.eliminated1,
                                     args.chrom_sizes1, args.cbr_window,
                                     args.terminal_tol, args.species1)
    annot2, diag2 = annotate_species(args.gff2, args.retained2, args.eliminated2,
                                     args.chrom_sizes2, args.cbr_window,
                                     args.terminal_tol, args.species2)

    merged = merge_with_rates(kaks, annot1, annot2, args.gene1_col, args.gene2_col)

    n1 = int(merged['region1'].notna().sum())
    n2 = int(merged['region2'].notna().sum())
    logging.info('Join: %s/%s rows carry a %s annotation, %s a %s annotation',
                 f'{n1:,}', f'{len(merged):,}', args.species1, f'{n2:,}',
                 args.species2)
    # A gene-identifier namespace mismatch between the rate table and the GFF is
    # the failure this script has actually hit before, and it produces an empty
    # table rather than an error. Fail here, with both namespaces on screen.
    for n, sp, gff, col in ((n1, args.species1, args.gff1, args.gene1_col),
                            (n2, args.species2, args.gff2, args.gene2_col)):
        if n == 0:
            sys.exit(f'ERROR: no {col} value matched an identifier in {gff}.\n'
                     f'  rate table: {list(kaks[col].head(3))}\n'
                     f'  GFF:        {list(annot1["gene_id"].head(3))}')

    out_table = os.path.join(args.outdir, 'kaks_annotated.tsv')
    # %.10g, not %.6g. This table carries genomic coordinates and CBR distances
    # up to ~2.4e8, and six significant figures rounds a distance by several base
    # pairs -- which then propagates into the changepoint fit and into
    # Supplementary Table T1. Ten covers the largest coordinate exactly and still
    # trims the rate columns to a sane width.
    merged.to_csv(out_table, sep='\t', index=False, float_format='%.10g')
    logging.info('Wrote %s (%s rows, %d columns)', out_table, f'{len(merged):,}',
                 merged.shape[1])

    # --- tests ------------------------------------------------------------
    metrics = [args.metric] + [m for m in (args.extra_metrics or [])
                               if m in merged.columns]
    tests = []
    for metric in metrics:
        tests += pairwise_tests(merged, 'region1', ['retained', 'eliminated'],
                                metric, 'region')
        tests += pairwise_tests(merged, 'consensus_region',
                                ['retained', 'eliminated'], metric, 'consensus')
        tests += pairwise_tests(merged, 'elim_type1', ['end', 'internal'],
                                metric, 'elim_subtype')
        # Pooled across regions: the contrast the next family takes apart.
        # Its own family, so adding it changes no existing q-value.
        tests += pairwise_tests(merged, 'proximity1',
                                ['interior', 'cbr_proximal'], metric,
                                'proximity_pooled')
        merged['_rxp'] = merged['region1'].astype(str) + '|' + \
            merged['proximity1'].astype(str)
        tests += pairwise_tests(
            merged, '_rxp',
            [f'{r}|{p}' for r, p in PROXIMITY_ORDER], metric,
            'region_x_proximity')

    stats = pd.DataFrame(tests)
    if stats.empty:
        sys.exit('ERROR: no comparison had two groups of at least 5 genes.')
    # BH within each family of tests. Correcting across families would penalise
    # the S7 comparisons for the number of diagnostics run beside them.
    stats['q'] = np.nan
    for fam, block in stats.groupby(['family', 'metric']):
        stats.loc[block.index, 'q'] = benjamini_hochberg(block['p'].to_numpy())
    stats.to_csv(os.path.join(args.outdir, 'annotate_regions_stats.tsv'),
                 sep='\t', index=False, float_format='%.6g')
    tests = stats.to_dict('records')

    # --- figures ----------------------------------------------------------
    # Drawn after every test; the random stream only places jitter.
    drawn = {}
    logging.info('S7A (S5A_elim_subtypes): eliminated subtypes')
    drawn['S7A (S5A_elim_subtypes)'] = plot_subtypes(
        merged, args.metric, args.metric_label, args.outdir,
        'S5A_elim_subtypes', rng, tests, cap_pct=args.cap_pct,
        box_alpha=args.box_alpha, panel_letters=args.panel_letters)

    logging.info('S7B (S5B_region_x_proximity): region x CBR proximity')
    drawn['S7B (S5B_region_x_proximity)'] = plot_region_x_proximity(
        merged, args.metric, args.metric_label, args.outdir,
        'S5B_region_x_proximity', rng, tests, args.cbr_window / 1000,
        cap_pct=args.cap_pct, box_alpha=args.box_alpha,
        panel_letters=args.panel_letters)

    logging.info('Diagnostics')
    drawn['region_boxplots (diagnostic)'] = plot_simple(
        merged, 'region1', ['retained', 'eliminated'], args.metric,
        args.metric_label, args.outdir, 'region_boxplots', rng, tests,
        'region', cap_pct=args.cap_pct, box_alpha=args.box_alpha)
    drawn['consensus_boxplots (diagnostic)'] = plot_simple(
        merged, 'consensus_region', ['retained', 'eliminated'],
        args.metric, args.metric_label, args.outdir,
        'consensus_boxplots', rng, tests, 'consensus',
        cap_pct=args.cap_pct, box_alpha=args.box_alpha)

    # --- summary ----------------------------------------------------------
    lines = ['Regional annotation of pairwise rates', '',
             f'Rate table          {os.path.basename(args.kaks)}, '
             f'{len(kaks):,} rows',
             f'Annotated           {len(merged):,} rows, {merged.shape[1]} columns',
             f'CBR window          {args.cbr_window/1000:g} kb',
             f'Terminal tolerance  {args.terminal_tol:,} bp', '']

    for sp, diag in ((args.species1, diag1), (args.species2, diag2)):
        c, a, cb = diag['region_counts'], diag['alignment'], diag['cbrs']
        lines += [f'{sp}',
                  f"  gene models       {a['n_records']:,} -> {a['n_ids']:,} "
                  f"identifiers",
                  f"  multiply aligned  {a['n_ids_multiply_aligned']:,}; "
                  f"{a['n_changed_vs_first_in_file']:,} pick a different locus "
                  f"than first-in-file",
                  f"                    ({a['n_moved_to_another_sequence']:,} of "
                  f"those move to another sequence, changing the region call)",
                  f"  CBRs              {cb['n_cbrs']:,} on "
                  f"{cb['n_sequences_with_cbrs']} sequences "
                  f"({cb['n_terminal_dropped']:,} of {cb['n_boundaries']:,} tract "
                  f"boundaries dropped as terminal)",
                  f"  retained          {c['retained']:,}",
                  f"  eliminated        {c['eliminated']:,}",
                  f"  ambiguous         {c['ambiguous']:,}",
                  f"  unannotated       {c['unannotated']:,}",
                  f"  CBR-proximal      {diag['n_cbr_proximal']:,}"]
        if cb['chrom_sizes_inferred']:
            lines.append('  NOTE: no sequence lengths given; a telomeric tract '
                         'edge may be counted as a CBR')
        lines.append('')

    for col, name in (('region1', f'{args.species1} region'),
                      ('consensus_region', 'consensus region'),
                      ('elim_type1', 'eliminated subtype'),
                      ('proximity1', 'CBR proximity, all regions pooled'),
                      ('_rxp', 'region x CBR proximity')):
        order = sorted(merged[col].dropna().unique())
        s = group_summary(merged, col, order, args.metric)
        lines += [f'{args.metric_label} by {name}',
                  s.to_string(index=False, float_format=lambda v: f'{v:.4f}'), '']

    lines += [f'Tests ({args.metric_label}, Mann-Whitney, BH within family; '
              f'proximity_pooled is a family of one, so its q equals its P)']
    for r in stats[stats['metric'] == args.metric].itertuples():
        lines.append(f'  {r.family:<20} {r.group1} vs {r.group2}: '
                     f'n {r.n1:,}/{r.n2:,}, medians {r.median1:.4f}/'
                     f'{r.median2:.4f}, P = {r.p:.3g}, q = {r.q:.3g} '
                     f'({stars(r.q)})')
    lines += ['',
              'Notes', '',
              '  signed_cbr_dist is POSITIVE on the eliminated side. Downstream',
              '  scripts should use it rather than re-deriving a sign from',
              '  min_cbr_dist, which is how Figures 4B and 4D came to disagree.',
              '',
              '  CBRs are boundaries of MERGED eliminated tracts, excluding those',
              '  within the terminal tolerance of a sequence end. This matches',
              '  gene_density_cbr.py exactly, so every panel of Figure 4 uses one',
              '  CBR set.',
              '',
              '  One locus per gene is the BEST alignment (Rank, then Identity,',
              '  then span), not the first record in the file. See the per-species',
              '  counts above for how many calls that changes.',
              '',
              f'  CBR-proximal (S7B, for the caption): within '
              f'{args.cbr_window / 1000:g} kb of a CBR.',
              '',
              f'  Panels draw points up to a cap, the {args.cap_pct:g}th percentile',
              '  of the highest group; points above it are not drawn and are',
              '  counted below. Every test uses the complete distribution. Boxes',
              '  are the interquartile range with the median; whiskers run from',
              '  the 5th to the 95th percentile.',
              '']
    for name, info in drawn.items():
        if not info:
            lines.append(f'  {name}: not drawn (too few genes)')
            continue
        lines.append(f"  {name}: cap {info['cap']:.4g}; points not drawn")
        w = max(len(g) for g in info['groups'])
        lines += [f'    {g:<{w}}  {h:>5,} of {n:,}' for g, h, n in
                  zip(info['groups'], info['hidden'], info['n'])]
    lines.append('')
    summary = '\n'.join(lines)
    with open(os.path.join(args.outdir, 'annotate_regions_summary.txt'), 'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(), 'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv, 'args': vars(args),
            'inputs': {k: fingerprint(v) for k, v in
                       (('kaks', args.kaks), ('gff1', args.gff1),
                        ('gff2', args.gff2), ('retained1', args.retained1),
                        ('eliminated1', args.eliminated1),
                        ('retained2', args.retained2),
                        ('eliminated2', args.eliminated2),
                        ('chrom_sizes1', args.chrom_sizes1),
                        ('chrom_sizes2', args.chrom_sizes2))},
            'species1': diag1, 'species2': diag2,
            'rows_annotated': int(len(merged)),
            'join': {'species1_matched': n1, 'species2_matched': n2},
            'figures': {name: ({k: v for k, v in info.items() if k != 'paths'}
                               if info else None)
                        for name, info in drawn.items()},
        }, fh, indent=2, default=float)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
