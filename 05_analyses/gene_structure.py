#!/usr/bin/env python3
"""
gene_structure.py - gene architecture in retained versus eliminated regions.

Figure 2D. Compares three per-gene measurements between regions:

    genomic span     end - start + 1 of the gene model
    intron fraction  1 - (coding length / genomic span)
    CDS length       length of the extracted coding sequence

and computes a fourth that no earlier analysis measured:

    intergenic distance   bp to the nearest neighbouring gene model

Writes:

    gene_structure_panel.{svg,pdf,png}   Figure 2D strip, 185 x 58 mm
    span_boxplot.{svg,pdf,png}           }
    intron_fraction_boxplot.*            } single metrics, for supplements
    cds_length_boxplot.*                 }
    intergenic_boxplot.*                 }
    gene_structure_per_gene.tsv          per-gene table, including intergenic
    gene_structure_summary.tsv           medians, IQRs, tests
    gene_structure_summary.txt           the same, formatted for reading
    provenance.json                      inputs, arguments and versions

Notes on method
---------------
**Lengths are log-scaled in the panels.** Span and CDS length are strongly
right-skewed; a linear axis set by the longest gene squashes both boxes into the
bottom few percent of the panel and makes a difference at P = 6 x 10^-13 look
like no difference at all. Intron fraction is bounded on [0, 1] and stays linear.

**GFF coordinates are 1-based and inclusive**, so a span is end - start + 1. The
earlier version used end - start, which is one base short on every gene. The
error is negligible for a span but it also fed the intron fraction, where it
biases the ratio in a direction that varies with gene size.

**Intergenic distance is computed from every gene model in the GFF, not only
from the genes being compared.** Restricting to the ortholog set would skip over
whatever lies between two orthologs and report the gap to the next *ortholog*
rather than to the next gene, inflating every distance by an amount that depends
on how densely that region is covered by orthologs. Gaps are measured between
consecutive non-overlapping models per chromosome; overlapping models give a
negative gap, which is set to zero and counted in the summary.

**Intergenic distance is a mediator, not a clean covariate.** Eliminated regions
are gene-poor as a property of what was eliminated, so spacing is downstream of
the same process that defines the groups. It is written to the per-gene table for
use as a covariate in the multivariate model, where that caveat can be stated;
it is deliberately *not* a panel of Figure 2D, where it would read as an
independent finding and would duplicate the gene-density profile of Figure 4B.

**Restricted to *Ascaris*.** *Parascaris* models are miniprot projections and
*Toxocara* models are TransDecoder ORFs on assembled transcripts, so intron
structure and genomic span are not measured comparably in either.

Changes 2026-10-13 (panel D checked against the conventions of Fig 2A-C)
------------------------------------------------------------------------
* **Boxes and whiskers come from every gene.** They were computed from the
  genes left after trimming the pooled values to their 1st-99th percentiles,
  so the box drawn was not exactly the box in the summary, and on the two log
  panels the lower whisker stopped at the trim rather than at a Tukey limit.
* **Whiskers are Tukey's 1.5 x IQR**, as in Fig 2A-C. On the log panels the
  IQR is taken on the log scale: in kilobases the lower limit of a
  right-skewed length is negative, so the lower whisker would always run to
  the shortest gene, however far below the box it sits.
* **The axis range is the 1st-99th percentile of the pooled values, widened
  to take in every whisker end.** Genes outside it are not drawn and are
  counted per metric and region in the summary. They used to be dropped
  without a count, and the summary called the trim "the central 99%" when
  it was the central 98%.
* **Significance markers come from the two-sided Mann-Whitney P**, the basis
  of the markers in Fig 2A-C and of the P values Results quotes. They came
  from BH q across the three metrics, beside a summary line saying panel A
  did the same, which it does not. q is still computed and tabulated, and
  the summary says whether it would change any marker.
* The jitter's sampling step (at most 1,200 genes per region, now
  --max-points) is reported, and --cap-pct reaches the panels: they ignored it
  and always used 99.
* Log-axis tick labels are plain numbers without a forced decimal ("10", not
  "10.0").

No statistic changes: gene_structure_per_gene.tsv and gene_structure_summary.tsv
are byte-identical to the previous version's on the same inputs.

Usage
-----
    python gene_structure.py \\
      --cds-fasta As_cds.fa \\
      --gff As_miniprot.gff \\
      --rates-table kaks_regions/kaks_annotated.tsv \\
      --gene-col gene1 --region-col region1 \\
      --outdir results/gene_structure_20260817

Requires numpy, pandas, scipy and matplotlib, plus paper_colors.py, which is
located automatically in lib/.
"""

import argparse
import hashlib
import json
import logging
import os
import platform
import subprocess
import sys
import textwrap
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.ticker import FuncFormatter, LogLocator
from scipy import stats

_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, 'lib'), os.path.join(os.path.dirname(_HERE), 'lib')):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, _cand)

from paper_colors import (COL, FONT, LINEWIDTH, NEUTRAL, REGION, apply_style,
                          audit, mm, panel_letter, save)

apply_style()

EXTS = ('.svg', '.pdf', '.png')

#: column, axis label, log scale. Order is the order of the sub-panels, left to
#: right, and matches the 05/22 legend: span, intron fraction, CDS length.
PANEL_METRICS = [
    ('genomic_span_kb', 'Genomic span (kb)', True),
    ('intron_fraction', 'Intron fraction\n(1 − CDS / span)', False),
    ('cds_length_kb', 'CDS length (kb)', True),
]


def stars(p):
    """Significance marker. Figure 2 uses markers throughout, so panel D matches
    panel A rather than printing e-notation beside it. Same cut-offs as
    dnds_decomposition.py and codon_usage_bias.py, applied, as there, to the
    two-sided Mann-Whitney P."""
    for cut, mark in ((1e-4, '****'), (1e-3, '***'), (1e-2, '**'), (5e-2, '*')):
        if p < cut:
            return mark
    return 'ns'


#: Marker cut-offs as legend text, from stars() above.
MARKER_KEY = '* P < 0.05, ** P < 0.01, *** P < 0.001, **** P < 10^-4'


def ordinal(x):
    """1 -> '1st', 2 -> '2nd', 99 -> '99th', 0.5 -> '0.5th', 99.5 -> '99.5th'."""
    if float(x) != int(x):
        return f'{x:g}th'
    n = int(x)
    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')
    if n % 100 in (11, 12, 13):
        suffix = 'th'
    return f'{n}{suffix}'


def box_stats(values, log=False, whis=1.5):
    """Box and whisker statistics for one group, in the form ``ax.bxp`` draws.

    Computed from every value passed, never from the values drawn. Median and
    quartiles are ``np.percentile`` with linear interpolation, which is what
    pandas ``median()`` and ``quantile()`` give in the summary table, so the
    box on the page is the box in the table.

    Whiskers follow Tukey, as matplotlib's default does in Fig 2A-C: each ends
    at the most extreme value within ``whis`` x IQR of the box, or at the box
    edge if no value lies beyond it. With ``log`` the IQR is taken on log10
    values. That is what the rule means on a log axis; in linear units the
    lower limit of a right-skewed length is negative, so the lower whisker
    would always run to the smallest value.
    """
    v = np.asarray(values, dtype=float)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    if log:
        l1, l3 = np.log10(q1), np.log10(q3)
        lo_lim = 10 ** (l1 - whis * (l3 - l1))
        hi_lim = 10 ** (l3 + whis * (l3 - l1))
    else:
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


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return 'n/a'
    if p <= 0:
        return '< 10^-300'
    if p >= 1e-3:
        return f'{p:.3f}'
    exp = int(np.floor(np.log10(p)))
    return f'{p / 10 ** exp:.2f} x 10^{exp}'


def benjamini_hochberg(pvals):
    """BH-adjusted q-values for one family of tests.

    The family is the three metrics drawn in panel D. Pooling them with tests
    from other panels would change every q on this panel, so the family is
    defined by what the reader sees side by side.
    """
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = n - rank
        prev = min(prev, p[idx] * n / i)
        q[idx] = prev
    return q


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_fasta_lengths(path):
    """{id: sequence length in bp}, keyed on the first token of the header."""
    lengths, gid, total = {}, None, 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if gid is not None:
                    lengths[gid] = total
                gid, total = line[1:].split()[0], 0
            else:
                total += len(line)
    if gid is not None:
        lengths[gid] = total
    logging.info('Parsed %s CDS lengths from %s', f'{len(lengths):,}', path)
    return lengths


def parse_gff_models(path):
    """Every mRNA feature as a row: id, chrom, start, end, strand, span.

    Span is end - start + 1: GFF coordinates are 1-based and inclusive.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'mRNA':
                continue
            gid = None
            for field in parts[8].split(';'):
                if field.startswith('ID='):
                    gid = field[3:].strip()
                    break
            if not gid:
                continue
            target = None
            for field in parts[8].split(';'):
                if field.startswith('Target='):
                    target = field[7:].split()[0].strip()
                    break
            start, end = int(parts[3]), int(parts[4])
            rows.append({'gene_id': gid, 'chrom': parts[0],
                         'start': start, 'end': end, 'strand': parts[6],
                         'target': target,
                         'genomic_span': end - start + 1})
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit(f'ERROR: no mRNA features found in {path}')
    logging.info('Parsed %s gene models from %s across %d sequences',
                 f'{len(df):,}', path, df['chrom'].nunique())
    return df


def compute_intergenic(models):
    """Distance to the nearest neighbouring *locus*, per gene.

    Overlapping models are merged into loci first. On a miniprot GFF this is not
    a nicety: a query protein can align more than once (Rank=1, Rank=2, ...) and
    the secondary alignments frequently sit on or beside the primary one, so
    treating every model as a separate gene puts a large spike of zero-length
    gaps into the distribution. A spike at exactly zero ties the rank test in
    both groups and destroys the power to detect any real difference, so an
    unmerged null would be uninterpretable.

    Uses every model in the GFF, not only the genes under comparison: the gap to
    the next *ortholog* is not the gap to the next *gene*, and the difference
    depends on ortholog density, which differs between the two regions being
    compared.

    Returns the frame with per-gene distances added and a diagnostic dict.
    """
    merged_rows, gene_locus, n_models_overlapping = [], {}, 0
    same_target = diff_target = 0

    for chrom, block in models.groupby('chrom', sort=False):
        b = block.sort_values(['start', 'end'])
        locus_start = locus_end = None
        members, targets = [], []
        loci = []
        for r in b.itertuples():
            if locus_start is None:
                locus_start, locus_end = r.start, r.end
                members, targets = [r.gene_id], [r.target]
                continue
            if r.start <= locus_end:            # overlaps the open locus
                n_models_overlapping += 1
                if r.target is not None and r.target in targets:
                    same_target += 1
                else:
                    diff_target += 1
                locus_end = max(locus_end, r.end)
                members.append(r.gene_id)
                targets.append(r.target)
            else:
                loci.append((locus_start, locus_end, members))
                locus_start, locus_end = r.start, r.end
                members, targets = [r.gene_id], [r.target]
        if locus_start is not None:
            loci.append((locus_start, locus_end, members))

        n = len(loci)
        starts = np.array([l[0] for l in loci], dtype=float)
        ends = np.array([l[1] for l in loci], dtype=float)
        up = np.full(n, np.nan)
        down = np.full(n, np.nan)
        if n > 1:
            up[1:] = starts[1:] - ends[:-1] - 1
            down[:-1] = starts[1:] - ends[:-1] - 1
        for i, (_, _, members) in enumerate(loci):
            for gid in members:
                gene_locus[gid] = (up[i], down[i], len(members))
        merged_rows.append(n)

    df = models.copy()
    df['intergenic_up'] = df['gene_id'].map(
        lambda g: gene_locus.get(g, (np.nan, np.nan, 1))[0])
    df['intergenic_down'] = df['gene_id'].map(
        lambda g: gene_locus.get(g, (np.nan, np.nan, 1))[1])
    df['locus_models'] = df['gene_id'].map(
        lambda g: gene_locus.get(g, (np.nan, np.nan, 1))[2])
    df['intergenic_nearest'] = df[['intergenic_up', 'intergenic_down']].min(
        axis=1, skipna=True)

    diag = {'n_models': int(len(models)), 'n_loci': int(sum(merged_rows)),
            'n_models_overlapping': int(n_models_overlapping),
            'overlap_same_target': int(same_target),
            'overlap_diff_target': int(diff_target)}
    return df, diag


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def save_figure(fig, outdir, stem):
    for problem in audit(fig, verbose=False):
        logging.warning('FIGURE AUDIT [%s]: %s', stem, problem)
    save(fig, outdir, stem, formats=[e.lstrip('.') for e in EXTS])
    plt.close(fig)
    logging.info("Saved: %s%s", stem, '/'.join(EXTS))


def _boxplot_axes(ax, df, metric, ylabel, rng, log=False, box_alpha=0.40,
                  max_points=1200, cap_pct=99.0, marker=None):
    """Retained versus eliminated for one metric. Returns the test result and
    the accounting of what is drawn.

    The test, the box and the whiskers use every gene. ``cap_pct`` and the
    whisker ends set the axis range: the (100 - cap_pct)th to cap_pct-th
    percentile of the pooled values, widened to take in every whisker end.
    Genes outside the range are not drawn, and are counted per region in
    ``res['not_drawn']`` for the summary.
    """
    regions = ('retained', 'eliminated')
    ret = df.loc[df['region'] == 'retained', metric].dropna()
    eli = df.loc[df['region'] == 'eliminated', metric].dropna()
    n_nonpositive = {'retained': 0, 'eliminated': 0}
    if log:
        # A log axis cannot show 0; the test in main() drops the same genes.
        n_nonpositive = {'retained': int((ret <= 0).sum()),
                         'eliminated': int((eli <= 0).sum())}
        ret, eli = ret[ret > 0], eli[eli > 0]
    if len(ret) < 5 or len(eli) < 5:
        ax.text(0.5, 0.5, 'insufficient data', transform=ax.transAxes,
                ha='center', va='center', fontsize=FONT['annotation'])
        ax.set_axis_off()
        return None

    _, p = stats.mannwhitneyu(ret, eli, alternative='two-sided')
    res = dict(n_ret=len(ret), n_elim=len(eli), med_ret=float(ret.median()),
               med_elim=float(eli.median()),
               fold=float(eli.median() / ret.median()) if ret.median() else np.nan,
               p=float(p), n_nonpositive=n_nonpositive)

    # Box and whiskers from every gene (box_stats), never from the genes drawn.
    box = {key: box_stats(vals.to_numpy(dtype=float), log=log)
           for key, vals in zip(regions, (ret, eli))}

    # Axis range: the (100 - cap_pct)th to cap_pct-th percentile of the
    # pooled values, widened to every whisker end so that no whisker is cut by
    # the axis. Points are drawn only inside it; the counts outside go to the
    # summary.
    pooled = pd.concat([ret, eli])
    lo, hi = np.percentile(pooled, [100 - cap_pct, cap_pct])
    lo = min(lo, *(b['whislo'] for b in box.values()))
    hi = max(hi, *(b['whishi'] for b in box.values()))
    shown = [ret[(ret >= lo) & (ret <= hi)], eli[(eli >= lo) & (eli <= hi)]]
    res['n_clipped'] = int(sum(len(a) - len(b) for a, b in
                               zip((ret, eli), shown)))
    res['range'] = (float(lo), float(hi))
    res['not_drawn'] = {key: {'below': int((vals < lo).sum()),
                              'above': int((vals > hi).sum())}
                        for key, vals in zip(regions, (ret, eli))}
    res['box'] = {key: {k: b[k] for k in ('q1', 'med', 'q3', 'whislo',
                                          'whishi', 'beyond_whiskers')}
                  for key, b in box.items()}

    if log:
        ax.set_yscale('log')
        # These distributions span barely more than one decade, so the default
        # log locator offers two labelled ticks and the reader cannot put a
        # number on a box. Label the 2x and 5x minors as plain numbers.
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=(2.0, 5.0),
                                              numticks=12))
        # '{:g}', not ScalarFormatter: with 0.5 on the axis ScalarFormatter
        # gives every label one decimal place ("10.0", "200.0").
        for fmt_axis in (ax.yaxis.set_major_formatter,
                         ax.yaxis.set_minor_formatter):
            fmt_axis(FuncFormatter(lambda v, _pos: f'{v:g}'))

    # bxp() draws precomputed statistics; boxplot() would recompute them from
    # whatever it is handed, which is how the trimmed values used to set them.
    bp = ax.bxp([{k: box[key][k] for k in ('med', 'q1', 'q3', 'whislo',
                                           'whishi', 'fliers')}
                 for key in regions],
                positions=[0, 1], widths=0.5, patch_artist=True,
                showfliers=False, zorder=2)
    for patch, key in zip(bp['boxes'], regions):
        patch.set_facecolor(to_rgba(REGION[key], box_alpha))
        patch.set_alpha(None)
        patch.set_edgecolor('black')
        patch.set_linewidth(LINEWIDTH['thin'])
    for part, z in (('whiskers', 4), ('caps', 4), ('medians', 5)):
        for artist in bp[part]:
            artist.set_color('black')
            artist.set_zorder(z)
            artist.set_linewidth(LINEWIDTH['data'] if part == 'medians'
                                 else LINEWIDTH['thin'])

    res['jittered'] = {}
    for i, (vals, key) in enumerate(zip(shown, regions)):
        v = vals.to_numpy()
        if len(v) > max_points:
            v = v[rng.choice(len(v), max_points, replace=False)]
        res['jittered'][key] = {'drawn': int(len(v)), 'in_range': int(len(vals))}
        ax.scatter(np.full(len(v), i) + rng.normal(0, 0.06, len(v)), v,
                   s=1.2, alpha=0.28, color=REGION[key], edgecolors='none',
                   zorder=3, rasterized=True)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Retained\n(n = {len(ret):,})',
                        f'Eliminated\n(n = {len(eli):,})'])
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.6, 1.6)

    # Bracket above every value inside the axis range, which includes both
    # whisker ends (each is a gene's value). On a log axis the arithmetic has
    # to happen in log space or the bracket lands somewhere unrelated to the
    # boxes.
    top = max(s.max() for s in shown)
    bottom = min(s.min() for s in shown)
    if log:
        lo_l, hi_l = np.log10(bottom), np.log10(top)
        span = hi_l - lo_l
        y = 10 ** (hi_l + span * 0.10)
        ax.set_ylim(10 ** (lo_l - span * 0.06), 10 ** (hi_l + span * 0.30))
        tick = 10 ** (hi_l + span * 0.14)
        text_y = 10 ** (hi_l + span * 0.16)
    else:
        span = top - bottom
        y = top + span * 0.08
        ax.set_ylim(bottom - span * 0.05, top + span * 0.30)
        tick = y + span * 0.05
        text_y = y + span * 0.06
    ax.plot([0, 0, 1, 1], [y, tick, tick, y], color=NEUTRAL['ink'],
            linewidth=LINEWIDTH['thin'], solid_joinstyle='miter', zorder=6)
    ax.text(0.5, text_y, marker if marker is not None else stars(p),
            ha='center', va='bottom', fontsize=FONT['annotation'],
            fontweight='bold', color=NEUTRAL['ink'], zorder=6)
    return res


def plot_structure_strip(df, outdir, stem, rng, markers, width='double',
                         height_mm=58.0, letter=None, box_alpha=0.40,
                         cap_pct=99.0, max_points=1200):
    """Figure 2D: span, intron fraction and CDS length side by side. Returns
    each metric's test result and drawing accounting (_boxplot_axes)."""
    w_mm, h_mm = COL[width] * 25.4, height_mm
    fig = plt.figure(figsize=(COL[width], mm(height_mm)))
    fig.set_layout_engine('none')
    left, right, bottom, top, gutter = 15.0, 3.0, 15.0, 4.5, 17.0
    each = (w_mm - left - right - gutter * 2) / 3
    y0, hgt = bottom / h_mm, (h_mm - bottom - top) / h_mm

    results = {}
    for i, (metric, ylabel, log) in enumerate(PANEL_METRICS):
        ax = fig.add_axes([(left + i * (each + gutter)) / w_mm, y0,
                           each / w_mm, hgt])
        results[metric] = _boxplot_axes(ax, df, metric, ylabel, rng, log=log,
                                        box_alpha=box_alpha,
                                        max_points=max_points,
                                        cap_pct=cap_pct,
                                        marker=markers.get(metric))
        if letter and i == 0:
            panel_letter(ax, letter)
    save_figure(fig, outdir, stem)
    return results


def plot_single(df, metric, ylabel, outdir, stem, rng, log=False,
                width='single', height_mm=62.0, box_alpha=0.40,
                cap_pct=99.0, max_points=1200, marker=None):
    w_mm, h_mm = COL[width] * 25.4, height_mm
    fig = plt.figure(figsize=(COL[width], mm(height_mm)))
    fig.set_layout_engine('none')
    left, right, bottom, top = 17.0, 3.0, 16.0, 4.0
    ax = fig.add_axes([left / w_mm, bottom / h_mm,
                       (w_mm - left - right) / w_mm,
                       (h_mm - bottom - top) / h_mm])
    res = _boxplot_axes(ax, df, metric, ylabel, rng, log=log,
                        box_alpha=box_alpha, max_points=max_points,
                        cap_pct=cap_pct, marker=marker)
    save_figure(fig, outdir, stem)
    return res


# ---------------------------------------------------------------------------
# Legend notes
# ---------------------------------------------------------------------------

#: Short names for the notes, with units.
NOTE_NAMES = {'genomic_span_kb': 'genomic span (kb)',
              'intron_fraction': 'intron fraction',
              'cds_length_kb': 'CDS length (kb)',
              'intergenic_kb': 'nearest gene (kb)'}


def _para(text, width=76):
    """One prose note, wrapped with the summary's two-space indent."""
    return textwrap.wrap(text, width=width, initial_indent='  ',
                         subsequent_indent='  ', break_on_hyphens=False)


def legend_notes(panel, single, markers, markers_q, args, missing_span,
                 missing_cds):
    """Everything the Figure 2D legend has to say about what is drawn, from
    this run's numbers: box and whisker rule, marker basis, axis range and the
    genes outside it per region, whisker ends, the jitter's sampling, and the
    genes a metric could not use. Returns summary lines (two-space indent)."""
    lo_pct, hi_pct = ordinal(100 - args.cap_pct), ordinal(args.cap_pct)
    drawn = [m for m, _, _ in PANEL_METRICS if panel.get(m)]
    logs = [NOTE_NAMES[m].split(' (')[0] for m, _, log in PANEL_METRICS
            if log and m in drawn]
    out = ['Legend notes (Figure 2D, from this run)']
    out += _para('Boxes: quartiles of every gene; centre line, the median (the '
                 'medians in the table above).')
    whisk = ('Whiskers: Tukey, to the most extreme gene within 1.5 x IQR of '
             'the box, as in panels A-C.')
    if logs:
        whisk += (f' On the log axes ({" and ".join(logs)}) the IQR is taken '
                  'on the log scale: in kilobases the lower limit of a '
                  'right-skewed length is negative, so the lower whisker would '
                  'always run to the shortest gene.')
    out += _para(whisk)
    differ = [m for m in markers if markers[m] != markers_q.get(m)]
    mark = 'Markers: two-sided Mann-Whitney P, as in panels A-C.'
    if differ:
        mark += (' BH q across the three metrics would change: ' +
                 '; '.join(f'{NOTE_NAMES[m]} {markers[m]} -> {markers_q[m]}'
                           for m in differ) + '.')
    else:
        mark += ' BH q across the three metrics gives the same markers.'
    out += _para(mark) + [f'    {MARKER_KEY}']
    out += _para(f'Axis range: the {lo_pct} to {hi_pct} percentile of the '
                 'pooled values, widened to take in every whisker end. Genes '
                 'outside it are not drawn; the boxes, whiskers and tests use '
                 'them.')
    out += [f"    {'':<20}{'':>18}{'not drawn (below / above)':>34}",
            f"    {'metric':<20}{'axis range':>18}{'retained':>17}{'eliminated':>17}"]
    for m in drawn:
        r = panel[m]
        nd = r['not_drawn']
        rng_txt = f"{r['range'][0]:.4g} to {r['range'][1]:.4g}"
        ret_txt = f"{nd['retained']['below']:,} / {nd['retained']['above']:,}"
        eli_txt = f"{nd['eliminated']['below']:,} / {nd['eliminated']['above']:,}"
        out.append(f"    {NOTE_NAMES[m]:<20}{rng_txt:>18}{ret_txt:>17}{eli_txt:>17}")
    out.append(f"    {'whisker ends':<20}{'':>18}{'retained':>17}{'eliminated':>17}")
    for m in drawn:
        b = panel[m]['box']
        out.append(f"    {NOTE_NAMES[m]:<20}{'lower':>18}"
                   f"{b['retained']['whislo']:>17.4g}{b['eliminated']['whislo']:>17.4g}")
        out.append(f"    {'':<20}{'upper':>18}"
                   f"{b['retained']['whishi']:>17.4g}{b['eliminated']['whishi']:>17.4g}")
    jit = []
    for key in ('retained', 'eliminated'):
        pairs = [(panel[m]['jittered'][key]['drawn'],
                  panel[m]['jittered'][key]['in_range']) for m in drawn]
        if all(d == n for d, n in pairs):
            jit.append(f'{key}, every gene in range')
        else:
            n_lo, n_hi = min(n for _, n in pairs), max(n for _, n in pairs)
            span_txt = f'{n_lo:,}' if n_lo == n_hi else f'{n_lo:,} to {n_hi:,}'
            jit.append(f'{key}, {args.max_points:,} of the {span_txt} in range')
    out += _para(f'Points: at most {args.max_points:,} genes per region are '
                 f'drawn in each panel, chosen at random (--seed {args.seed}): '
                 + '; '.join(jit) + '. Boxes, whiskers and tests use every '
                 'gene.')
    if missing_span or missing_cds:
        out += _para(f'Genes left out of the metrics that need them: '
                     f'{missing_span:,} without a GFF model, {missing_cds:,} '
                     'without a CDS entry; n per panel is in its tick labels.')
    for m, r in single.items():
        if r and sum(r['n_nonpositive'].values()):
            c = r['n_nonpositive']
            out += _para(f"{NOTE_NAMES[m].capitalize()}: {c['retained']:,} "
                         f"retained and {c['eliminated']:,} eliminated genes "
                         'at 0 are left out of its log axis and its test.')
    out += _para('The single-metric figures (*_boxplot) follow the same rules; '
                 'intergenic_boxplot is not a panel of Figure 2D.')
    return out


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
        description='Gene architecture in retained versus eliminated regions '
                    '(Figure 2D), plus per-gene intergenic distance.')
    p.add_argument('--cds-fasta', required=True, help='CDS FASTA, nucleotide')
    p.add_argument('--gff', required=True,
                   help='Gene models. Spans and intergenic distances come from '
                        'here; intergenic uses every model, not only the genes '
                        'being compared.')
    p.add_argument('--rates-table', required=True,
                   help='Table carrying gene IDs and per-gene region calls')
    p.add_argument('--gene-col', default='gene1')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--box-alpha', type=float, default=0.40,
                   help='Opacity of the box fill; the jitter is drawn over it '
                        '(default: 0.40)')
    p.add_argument('--cap-pct', type=float, default=99.0,
                   help='Axis range: the (100 - x)th to x-th percentile of the '
                        'pooled values, widened to take in every whisker end. '
                        'Genes outside it are not drawn and are counted in the '
                        'summary; boxes, whiskers and tests use every gene. '
                        'Before 2026-10-13 the panels ignored this option and '
                        'always used 99 (default: 99)')
    p.add_argument('--max-points', type=int, default=1200,
                   help='Most genes per region drawn as jittered points, chosen '
                        'at random with --seed; boxes, whiskers and tests use '
                        'every gene (default: 1200)')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--width', choices=['single', 'onehalf', 'double'],
                   default='double')
    p.add_argument('--height-mm', type=float, default=58.0)
    p.add_argument('--panel-letters', action='store_true')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'gene_structure_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    # mode='w': FileHandler appends by default, so rerunning into an existing
    # directory concatenates logs while the figures are overwritten.
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(os.path.join(args.outdir, 'gene_structure.log'),
                                      mode='w'),
                  logging.StreamHandler(sys.stdout)])
    for noisy in ('matplotlib.font_manager', 'fontTools', 'PIL'):
        logging.getLogger(noisy).setLevel(logging.ERROR)

    logging.info('=' * 60)
    logging.info('Gene structure')
    logging.info('=' * 60)

    for path in (args.cds_fasta, args.gff, args.rates_table):
        if not os.path.isfile(path):
            sys.exit(f'ERROR: input not found: {path}')

    rng = np.random.default_rng(args.seed)
    cds_lengths = parse_fasta_lengths(args.cds_fasta)
    models = parse_gff_models(args.gff)
    models, diag = compute_intergenic(models)
    logging.info('Intergenic: %s models merged into %s loci; %s models '
                 'overlapped a neighbour (%s sharing a query protein, %s not)',
                 f"{diag['n_models']:,}", f"{diag['n_loci']:,}",
                 f"{diag['n_models_overlapping']:,}",
                 f"{diag['overlap_same_target']:,}",
                 f"{diag['overlap_diff_target']:,}")
    if diag['n_models_overlapping'] > 0.05 * diag['n_models']:
        logging.warning('%.0f%% of models overlap a neighbour. Distances are '
                        'measured between merged loci, so this does not put a '
                        'spike of zeros into the distribution -- but check the '
                        'same-query share above: a high value means the GFF '
                        'carries redundant secondary alignments rather than '
                        'genuinely overlapping genes.',
                        100 * diag['n_models_overlapping'] / diag['n_models'])

    table = pd.read_csv(args.rates_table, sep='\t')
    for col in (args.gene_col, args.region_col):
        if col not in table.columns:
            sys.exit(f"ERROR: no '{col}' column in {args.rates_table}. "
                     f'Columns: {table.columns.tolist()}')

    regions = table[[args.gene_col, args.region_col]].rename(
        columns={args.gene_col: 'gene_id', args.region_col: 'region'})
    regions = regions[regions['region'].isin(('retained', 'eliminated'))]

    df = regions.merge(models, on='gene_id', how='left')
    df['cds_length'] = df['gene_id'].map(cds_lengths)
    missing_span = int(df['genomic_span'].isna().sum())
    missing_cds = int(df['cds_length'].isna().sum())
    if missing_span or missing_cds:
        logging.warning('%s genes have no GFF model and %s have no CDS entry; '
                        'they are dropped per metric, not globally',
                        f'{missing_span:,}', f'{missing_cds:,}')

    df['genomic_span_kb'] = df['genomic_span'] / 1000.0
    df['cds_length_kb'] = df['cds_length'] / 1000.0
    valid = df['genomic_span'] > 0
    df['intron_fraction'] = np.where(
        valid, 1.0 - (df['cds_length'] / df['genomic_span']), np.nan)
    out_of_range = int(((df['intron_fraction'] < 0) |
                        (df['intron_fraction'] > 1)).sum())
    if out_of_range:
        logging.warning('%s genes have an intron fraction outside [0, 1] and '
                        'are clipped. A CDS longer than its own span means the '
                        'FASTA and the GFF disagree for that gene, so check the '
                        'count rather than ignoring it.', f'{out_of_range:,}')
    df['intron_fraction'] = df['intron_fraction'].clip(0, 1)
    df['intergenic_kb'] = df['intergenic_nearest'] / 1000.0

    logging.info('Merged %s genes (%s retained, %s eliminated)',
                 f'{len(df):,}', f'{(df.region == "retained").sum():,}',
                 f'{(df.region == "eliminated").sum():,}')

    # --- tests, with BH across the three panel metrics ----------------------
    METRICS = PANEL_METRICS + [('intergenic_kb',
                                'Distance to nearest gene (kb)', True)]
    raw = {}
    for metric, _, log in METRICS:
        ret = df.loc[df.region == 'retained', metric].dropna()
        eli = df.loc[df.region == 'eliminated', metric].dropna()
        if log:
            ret, eli = ret[ret > 0], eli[eli > 0]
        if len(ret) < 5 or len(eli) < 5:
            continue
        _, p = stats.mannwhitneyu(ret, eli, alternative='two-sided')
        raw[metric] = {'p': float(p), 'med_ret': float(ret.median()),
                       'med_elim': float(eli.median()),
                       'n_ret': len(ret), 'n_elim': len(eli)}

    panel_keys = [m for m, _, _ in PANEL_METRICS if m in raw]
    qvals = benjamini_hochberg([raw[m]['p'] for m in panel_keys])
    for m, q in zip(panel_keys, qvals):
        raw[m]['q'] = float(q)
    # Markers from P, as in panels A-C and as Results quotes; q stays in the
    # tables, and the summary says whether it would change any marker.
    markers = {m: stars(raw[m]['p']) for m in panel_keys}
    markers_q = {m: stars(raw[m]['q']) for m in panel_keys}

    for metric, label, _ in METRICS:
        r = raw.get(metric)
        if not r:
            continue
        q_txt = f", q = {fmt_p(r['q'])}" if 'q' in r else ''
        logging.info('  %-18s retained %10.4f  eliminated %10.4f  P = %s%s',
                     metric, r['med_ret'], r['med_elim'], fmt_p(r['p']), q_txt)

    # --- figures ------------------------------------------------------------
    panel = plot_structure_strip(df, args.outdir, 'gene_structure_panel', rng,
                                 markers, width=args.width,
                                 height_mm=args.height_mm,
                                 box_alpha=args.box_alpha,
                                 cap_pct=args.cap_pct,
                                 max_points=args.max_points,
                                 letter='D' if args.panel_letters else None)
    single = {}
    for metric, label, log in METRICS:
        single[metric] = plot_single(df, metric, label, args.outdir,
                                     metric.replace('_kb', '') + '_boxplot',
                                     rng, log=log, box_alpha=args.box_alpha,
                                     cap_pct=args.cap_pct,
                                     max_points=args.max_points,
                                     marker=markers.get(metric))
    notes = legend_notes(panel, single, markers, markers_q, args,
                         missing_span, missing_cds)
    logging.info('LEGEND NOTES (Figure 2D)')
    for line in notes:
        logging.info('%s', line)

    # --- outputs ------------------------------------------------------------
    keep = ['gene_id', 'region', 'chrom', 'start', 'end', 'strand',
            'genomic_span', 'cds_length', 'intron_fraction',
            'intergenic_up', 'intergenic_down', 'intergenic_nearest',
            'locus_models']
    df[keep].to_csv(os.path.join(args.outdir, 'gene_structure_per_gene.tsv'),
                    sep='\t', index=False)

    summary_rows = []
    for metric, label, _ in METRICS:
        for region in ('retained', 'eliminated'):
            vals = df.loc[df.region == region, metric].dropna()
            summary_rows.append({'metric': metric, 'region': region,
                                 'n': len(vals), 'median': vals.median(),
                                 'mean': vals.mean(),
                                 'Q1': vals.quantile(0.25),
                                 'Q3': vals.quantile(0.75)})
        if metric in raw:
            summary_rows.append({'metric': metric, 'region': 'MWU',
                                 'n': None, 'median': raw[metric]['p'],
                                 'mean': raw[metric].get('q'),
                                 'Q1': None, 'Q3': None})
    pd.DataFrame(summary_rows).to_csv(
        os.path.join(args.outdir, 'gene_structure_summary.tsv'),
        sep='\t', index=False)

    lines = ['Gene structure (Figure 2D)', '=' * 60, '',
             f'CDS      {args.cds_fasta}',
             f'Models   {args.gff}',
             f'Regions  {args.region_col} in {args.rates_table}', '',
             f"{'metric':<20}{'retained':>12}{'eliminated':>12}{'P':>16}{'q':>16}", '']
    for metric, label, _ in METRICS:
        r = raw.get(metric)
        if not r:
            continue
        lines.append(f"  {metric:<18}{r['med_ret']:>12.4f}{r['med_elim']:>12.4f}"
                     f"{fmt_p(r['p']):>16}"
                     f"{fmt_p(r['q']) if 'q' in r else '-':>16}")
    lines += ['',
              '  q is Benjamini-Hochberg across the three metrics drawn in',
              '  panel D. Intergenic distance is not part of that family and',
              '  carries no q: it is not a panel of Figure 2D. The markers on',
              '  the panels come from P, not q (legend notes below).', '']
    lines += notes
    lines += ['',
              'Reporting notes',
              '  Span and CDS length are drawn on a log axis. Both are strongly',
              '  right-skewed, and a linear axis set by the longest gene puts',
              '  both boxes in the bottom few percent of the panel -- which is',
              '  how a difference at P = 6 x 10^-13 came to look like no',
              '  difference in the previous version.', '',
              'Intergenic distance',
              '  Computed from every gene model in the GFF, not only from the',
              '  genes being compared: the gap to the next ortholog is not the',
              '  gap to the next gene, and the difference depends on ortholog',
              '  density, which is exactly what differs between these regions.',
              f"  {diag['n_models']:,} models merge into {diag['n_loci']:,} loci; "
              f"{diag['n_models_overlapping']:,} models overlap a",
              f"  neighbour ({diag['overlap_same_target']:,} share a query protein, "
              f"{diag['overlap_diff_target']:,} do not). Distances",
              '  are measured between loci, not between models: a query protein',
              '  can align more than once and the secondary alignments sit on or',
              '  beside the primary, so per-model gaps would put a large spike of',
              '  zeros into the distribution, tie the rank test in both groups',
              '  and make any null uninterpretable.',
              '  It is a mediator, not a clean covariate: eliminated regions are',
              '  gene-poor as a property of what was eliminated, so spacing sits',
              '  downstream of the same process that defines the groups.',
              '  Conditioning on it can remove real signal. It is written to the',
              '  per-gene table for the multivariate model, where that can be',
              '  stated, and is deliberately not a panel here.']
    if out_of_range:
        lines += ['', f'  {out_of_range:,} genes had an intron fraction outside '
                      '[0, 1] and were clipped.']
    summary = '\n'.join(lines) + '\n'
    with open(os.path.join(args.outdir, 'gene_structure_summary.txt'), 'w') as fh:
        fh.write(summary)

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(),
            'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv,
            'args': vars(args),
            'inputs': {k: fingerprint(v) for k, v in
                       (('cds_fasta', args.cds_fasta), ('gff', args.gff),
                        ('rates_table', args.rates_table))},
            'genes_analysed': int(len(df)),
            'intergenic': diag,
            'intron_fraction_clipped': out_of_range,
            'markers_from': 'two-sided Mann-Whitney P',
            'markers': markers,
            'markers_if_q': markers_q,
            'panel_d': panel,
            'single_metric_figures': single,
            'legend_notes': notes,
        }, fh, indent=2)

    print()
    print(summary)
    logging.info('Summary: %s',
                 os.path.join(args.outdir, 'gene_structure_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
