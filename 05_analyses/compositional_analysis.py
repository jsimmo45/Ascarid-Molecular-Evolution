#!/usr/bin/env python3
"""
compositional_analysis.py -- base composition of eliminated and retained
Ascaris genes, and whether it could produce the dN/dS contrast
(Supplementary Figure S2; the file keeps its name from before the renumbering,
figS1_composition, and "S1" in the notes below is this figure).

2026-10-19 (results draft note 67): the C key's "Above the axis" marker is a
blue and a red open triangle side by side, as drawn (one black triangle matched
neither; fixed in S12A by note 40, not here), and every P is plain text with
its exponent raised as its own 6 pt text (plain_text.draw_p), not mathtext.

Composition is measured on the Ascaris coding sequence of every analysed
Ascaris x P. univalens ortholog pair (the rate table annotated by
annotate_regions.py, ANNOTATED_TABLE in config.sh) and set against the
pairwise dN/dS of the same pairs.

Measures
--------
  GC3s  G or C at synonymous third positions: the third base of every sense
        codon except ATG (Met) and TGG (Trp), which have no synonymous third
        position, counted over unambiguous codons. Stop codons, internal or
        terminal, are not counted. This is the definition codon_usage_bias.py
        uses for Figure 2 (and codonW's), so S1 and Figure 2 quote one number.
  GC    G + C over all unambiguous bases of the coding sequence, terminal stop
        codon excluded.

Gene set. With --codon-per-gene (codon_usage_bias.py's per-gene table, the
wrapper's default), S1 analyses exactly the genes that table holds, so the GC3s
test here and in Results' codon paragraph are one test on one gene set.
Without it, genes with fewer than --min-codons informative codons (sense codons
other than ATG and TGG; default 100) are excluded. Both counts are reported.
(2026-10-12: the first real run used S1's own rule, kept 9,064 retained genes
against Fig 2's 9,085, and so printed P = 1.1e-45 beside Fig 2's 1.4e-45.)

Tests
-----
  1. GC3s and GC, eliminated vs retained: Mann-Whitney, BH q over the two.
  2. dN/dS vs GC3s: Spearman, across all genes (the figure Results quotes) and
     within each region. The pooled coefficient mixes the region effect into
     the composition effect; the within-retained one does not.
  3. GC3s-matched contrast. Results argues that because eliminated genes are
     GC3s-richer and GC3s correlates negatively with dN/dS, composition biases
     the dN/dS contrast DOWNWARD. That is an inference from two numbers; this
     tests it directly. Each eliminated gene is paired with a retained gene
     within --caliper of its GC3s (drawn without replacement), the median
     dN/dS fold and Mann-Whitney P are recomputed, and this is repeated
     --n-matched times -- the design of the GC3s-matched codon-usage control.
     A matched fold at or above the unmatched one confirms the sentence.

Panels (one 185 mm strip, figS1_composition)
--------------------------------------------
  A  GC3s by region     box = quartiles, whiskers = 5th-95th percentiles,
  B  GC by region       jitter on top; asterisks from BH q, exact P printed
  C  dN/dS against GC3s, eliminated drawn last so its 322 points are not
     buried under 9,085 retained ones; y-axis stops at the --dnds-cap-pct
     percentile and genes above it are open triangles on the top edge,
     counted in the key.

What changed from the May script
--------------------------------
  * GC3 counted ATG and TGG, whose third base is always G, so it was inflated
    by a per-gene amount that depends on Met/Trp content. Now GC3s, matching
    codon_usage_bias.py. The May numbers (0.487 vs 0.443) are superseded.
  * The gene-density meta-profile is gone: Figure 4B (gene_density_cbr.py)
    replaced it, and the May version counted every BED edge as a CBR.
  * House figure pipeline (paper_colors); e-notation; retained left / blue,
    eliminated right / red; n in the tick labels; nothing capped silently.
  * New: the GC3s-matched contrast, within-region Spearman, and an optional
    per-gene cross-check against codon_usage_bias.py output (--codon-per-gene).

Outputs (in --outdir)
---------------------
  figS1_composition.{svg,pdf,png}
  composition_per_gene.tsv      every analysed pair, with GC3s, GC and QC counts
  composition_tests.tsv         every test, with BH q where a family applies
  matched_replicates.tsv        one row per GC3s-matched replicate
  composition_summary.txt       results, and the Results sentence they support
  compositional_analysis.log
  provenance.json

Usage
-----
  python compositional_analysis.py \\
    --annotated results/annotate_regions_20260817/kaks_annotated.tsv \\
    --cds As_cds.fa --panel-letters \\
    --outdir results/compositional_analysis_<stamp>

Dependencies: numpy, pandas, scipy, matplotlib; paper_colors, region_lib and
rate_tables from lib/.
"""

from __future__ import annotations

import argparse
import gzip
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
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation
from scipy.stats import mannwhitneyu, spearmanr
from scipy import __version__ as scipy_version

# A script that cannot import the shared style must stop, not fall back to
# matplotlib defaults.
from paper_colors import (REGION, NEUTRAL, FONT, LINEWIDTH, apply_style, audit,
                          figure, save, mm)
from region_lib import (benjamini_hochberg, fingerprint, git_commit,
                        setup_logging, stars)
from rate_tables import read_rate_table
from plain_text import draw_p

_HERE = os.path.dirname(os.path.abspath(__file__))

REGIONS = ('retained', 'eliminated')
STOPS = frozenset({'TAA', 'TAG', 'TGA'})
#: Sense codons with no synonymous third position.
SINGLE_CODON = frozenset({'ATG', 'TGG'})
_BASES = 'ACGT'
#: Codons whose third base is a synonymous site: every unambiguous sense codon
#: except ATG and TGG (standard code).
INFORMATIVE = frozenset(a + b + c for a in _BASES for b in _BASES
                        for c in _BASES) - STOPS - SINGLE_CODON


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------

def read_fasta(path, wanted=None):
    """{id: sequence} for the records in `wanted` (all records if None).

    The ID is the first whitespace-delimited token of the header. A repeated
    ID keeps its first record; repeats are returned as {id: identical?} so the
    caller can stop when a repeated ID is one the analysis uses and the copies
    disagree -- which copy is "right" is not a question to answer silently.
    """
    seqs, repeats = {}, {}
    n_records = 0
    opener = gzip.open if path.endswith('.gz') else open
    cur, chunks = None, []

    def store():
        if cur is None or (wanted is not None and cur not in wanted):
            return
        s = ''.join(chunks).upper()
        if cur in seqs:
            repeats[cur] = repeats.get(cur, True) and seqs[cur] == s
        else:
            seqs[cur] = s

    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                store()
                n_records += 1
                toks = line[1:].split()
                if not toks:
                    sys.exit(f'ERROR: {path}: a FASTA header with no ID '
                             f'(record {n_records:,}).')
                cur, chunks = toks[0], []
            else:
                chunks.append(line.strip())
        store()
    if n_records == 0:
        sys.exit(f'ERROR: no FASTA records in {path}')
    return seqs, repeats, n_records


def composition(seq):
    """GC3s, GC and gene-model QC counts for one coding sequence.

    Codons are read from the first base in frame; bases beyond the last whole
    codon are counted as `trailing_bases` and ignored. A final stop codon is
    excluded from everything; a stop codon elsewhere is counted as
    `internal_stops` and excluded from GC3s, but its bases stay in GC because
    they are part of the sequence as annotated.
    """
    s = seq.replace('U', 'T')
    n_full = len(s) // 3
    codons = [s[3 * i:3 * i + 3] for i in range(n_full)]
    terminal_stop = bool(codons) and codons[-1] in STOPS
    body = codons[:-1] if terminal_stop else codons

    n_inf = gc3s = internal = ambiguous = n_single = 0
    for c in body:
        if c in INFORMATIVE:
            n_inf += 1
            if c[2] in 'GC':
                gc3s += 1
        elif c in STOPS:
            internal += 1
        elif c in SINGLE_CODON:
            n_single += 1
        else:
            ambiguous += 1
    body_str = s[:3 * len(body)]
    gc = body_str.count('G') + body_str.count('C')
    acgt = gc + body_str.count('A') + body_str.count('T')
    return {'n_codons': len(body), 'n_informative': n_inf,
            'n_met_trp': n_single, 'internal_stops': internal,
            'ambiguous_codons': ambiguous, 'terminal_stop': terminal_stop,
            'trailing_bases': len(s) % 3,
            'GC3s': gc3s / n_inf if n_inf else np.nan,
            'GC': gc / acgt if acgt else np.nan}


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def p_text(p):
    """Plain-text P for the summary (the figure uses plain_text.draw_p)."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p == 0:
        return '< 1e-300'
    return f'{p:.2g}' if p < 0.001 else f'{p:.3f}'


def p_manuscript(p):
    """P in the manuscript's notation, e.g. 1.4 × 10⁻⁴⁵."""
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


def signed(x, nd):
    """x to `nd` decimals with a true minus sign (U+2212), as typeset. A value
    that rounds to zero prints unsigned."""
    s = f'{x:.{nd}f}'
    if float(s) == 0:
        s = s.lstrip('-')
    return s.replace('-', '\u2212')


def region_test(df, col):
    """Eliminated vs retained on `col`: n, medians, quartiles, Mann-Whitney."""
    r = df.loc[df['region'] == 'retained', col].dropna().to_numpy(float)
    e = df.loc[df['region'] == 'eliminated', col].dropna().to_numpy(float)
    out = {'metric': col, 'n_retained': len(r), 'n_eliminated': len(e),
           'median_retained': np.median(r) if len(r) else np.nan,
           'median_eliminated': np.median(e) if len(e) else np.nan,
           'q1_retained': np.percentile(r, 25) if len(r) else np.nan,
           'q3_retained': np.percentile(r, 75) if len(r) else np.nan,
           'q1_eliminated': np.percentile(e, 25) if len(e) else np.nan,
           'q3_eliminated': np.percentile(e, 75) if len(e) else np.nan,
           'U': np.nan, 'p': np.nan}
    out['difference'] = out['median_eliminated'] - out['median_retained']
    out['fold'] = (out['median_eliminated'] / out['median_retained']
                   if out['median_retained'] else np.nan)
    if len(r) >= 5 and len(e) >= 5:
        out['U'], out['p'] = mannwhitneyu(e, r, alternative='two-sided')
    return out


def spearman_block(df, x, y):
    """Spearman of y on x: pooled over both regions, then within each."""
    rows = []
    for name, sub in (('all genes', df),
                      ('retained', df[df['region'] == 'retained']),
                      ('eliminated', df[df['region'] == 'eliminated'])):
        s = sub[[x, y]].dropna()
        rho, p = (spearmanr(s[x], s[y]) if len(s) >= 10 else (np.nan, np.nan))
        rows.append({'subset': name, 'n': len(s), 'rho': float(rho),
                     'p': float(p)})
    return rows


def matched_contrast(ret_gc, ret_val, elim_gc, elim_val, caliper, n_rep, seed):
    """GC3s-matched eliminated vs retained comparison of `val` (dN/dS).

    Each replicate visits the eliminated genes in random order and gives each
    one a retained partner drawn at random from the retained genes within
    `caliper` of its GC3s that no earlier eliminated gene has taken. An
    eliminated gene with no free partner is dropped from that replicate and
    counted. Returns (per-replicate DataFrame, n eliminated with no retained
    gene inside the caliper at all).
    """
    order = np.argsort(ret_gc, kind='stable')
    rg, rv = ret_gc[order], ret_val[order]
    lo = np.searchsorted(rg, elim_gc - caliper, side='left')
    hi = np.searchsorted(rg, elim_gc + caliper, side='right')
    n_no_window = int((hi <= lo).sum())
    rng = np.random.default_rng(seed)
    rows = []
    for rep in range(n_rep):
        used = np.zeros(len(rg), dtype=bool)
        pick = np.full(len(elim_gc), -1)
        for i in rng.permutation(len(elim_gc)):
            if hi[i] <= lo[i]:
                continue
            cand = np.arange(lo[i], hi[i])
            cand = cand[~used[cand]]
            if not len(cand):
                continue
            j = cand[rng.integers(len(cand))]
            used[j] = True
            pick[i] = j
        ok = pick >= 0
        e, m = elim_val[ok], rv[pick[ok]]
        med_e, med_m = np.median(e), np.median(m)
        p = mannwhitneyu(e, m, alternative='two-sided')[1]
        p_gc = mannwhitneyu(elim_gc[ok], rg[pick[ok]], alternative='two-sided')[1]
        rows.append({'replicate': rep, 'n_matched': int(ok.sum()),
                     'median_eliminated': med_e, 'median_matched_retained': med_m,
                     'fold': med_e / med_m if med_m else np.nan, 'p': p,
                     'gc3s_balance_p': p_gc})
    return pd.DataFrame(rows), n_no_window


def codon_gene_set(path):
    """Genes with a finite GC3s in codon_usage_bias.py's per-gene table.

    That table is Figure 2's gene set. S1 adopts it, when given, so that the
    GC3s test in S1 and the one in Results' codon paragraph are the same test
    on the same genes, whatever minimum-codon rule the codon script applied.
    Returns (set, gene column, value column) or None if the table has no
    usable columns.
    """
    tab = pd.read_csv(path, sep='\t')
    lower = {c.lower(): c for c in tab.columns}
    gcol = next((lower[k] for k in ('gene_id', 'gene1', 'gene') if k in lower),
                None)
    vcol = next((lower[k] for k in ('gc3s', 'gc3') if k in lower), None)
    if gcol is None or vcol is None:
        logging.warning('  %s has no gene/GC3s column (columns: %s); S1 keeps '
                        'its own codon filter', path, list(tab.columns))
        return None
    ok = pd.to_numeric(tab[vcol], errors='coerce').notna()
    genes = set(tab.loc[ok, gcol].astype(str))
    logging.info('  codon_usage_bias.py gene set: %s genes with %s (%s rows)',
                 f'{len(genes):,}', vcol, f'{len(tab):,}')
    return genes, gcol, vcol


def cross_check(per_gene, path):
    """Compare per-gene GC3s with codon_usage_bias.py output.

    Accepts a gene column named gene_id, gene1 or gene and a GC3s column named
    GC3s or, failing that, GC3 (the converted codon script computes GC3s; the
    column may keep the old name). Returns a dict for the summary, or None.
    """
    tab = pd.read_csv(path, sep='\t')
    lower = {c.lower(): c for c in tab.columns}
    gcol = next((lower[k] for k in ('gene_id', 'gene1', 'gene') if k in lower),
                None)
    vcol = next((lower[k] for k in ('gc3s', 'gc3') if k in lower), None)
    if gcol is None or vcol is None:
        logging.warning('  cross-check skipped: %s has no gene/GC3s column '
                        '(columns: %s)', path, list(tab.columns))
        return None
    ref = tab.set_index(tab[gcol].astype(str))[vcol]
    ref = ref[~ref.index.duplicated(keep='first')]
    mine = per_gene.set_index('gene1')['GC3s']
    shared = mine.index.intersection(ref.index)
    if not len(shared):
        logging.warning('  cross-check: no gene in common with %s', path)
        return {'path': path, 'column': vcol, 'n_shared': 0}
    diff = (mine.loc[shared].astype(float)
            - pd.to_numeric(ref.loc[shared], errors='coerce')).abs()
    n_bad = int((diff > 1e-6).sum())
    out = {'path': path, 'column': vcol, 'n_shared': int(len(shared)),
           'n_differ': n_bad, 'max_abs_diff': float(diff.max()),
           'examples': list(diff[diff > 1e-6].sort_values(ascending=False)
                            .index[:5])}
    logging.info('  cross-check vs %s (%s): %s genes shared, %s differ by '
                 '> 1e-6 (max |diff| %.4g)', os.path.basename(path), vcol,
                 f'{len(shared):,}', f'{n_bad:,}', out['max_abs_diff'])
    return out


# ---------------------------------------------------------------------------
# Figure helpers (shared pattern with functional_enrichment.py and
# multivariate_insulation.py)
# ---------------------------------------------------------------------------

def blank_figure(width, height_mm):
    """A journal-width canvas with no axes and no layout engine; panels are
    placed explicitly, because constrained_layout shrinks a panel to nothing
    to make room for text it cannot fit and only warns."""
    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax)
    return fig


def style_box(bp, faces, box_alpha):
    """Translucent face (alpha in the colour, not set_alpha, so the black edge
    stays black); box 2, jitter 3, whiskers and caps 4, median 5."""
    for patch, face in zip(bp['boxes'], faces):
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
    """Lay out panels side by side from their measured text.

    Each panel's own labels (left), anything right of it and the legend or
    text above are measured, the gaps between neighbours are kept, and the
    remaining width is split in `shares`. The row fills the canvas exactly
    whatever font the node renders with. `top_extra_mm` reserves a band above
    everything (here, for the panel letters).
    """
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
    """Audit, check text, save in every format, close. Returns
    (paths, problems); problems go into the summary as well as the log."""
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

def draw_box(ax, df, col, ylabel, test, rng, jitter_max, range_pct,
             box_alpha=0.40):
    """A or B: `col` by region, retained left, eliminated right.

    Box = quartiles, whiskers = 5th-95th percentiles, jitter on top. The
    y-range covers the central `range_pct` of the pooled values; jitter
    beyond it is drawn on the edge and counted. Drawing only -- the test uses
    every gene. Returns what was not drawn as it is.
    """
    data = [df.loc[df['region'] == g, col].dropna().to_numpy(float)
            for g in REGIONS]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.55, patch_artist=True,
                    showfliers=False, whis=(5, 95), zorder=2)
    style_box(bp, [REGION[g] for g in REGIONS], box_alpha)
    pooled = np.concatenate(data)
    lo, hi = np.percentile(pooled, range_pct)
    n_edge, n_skip = {}, {}
    for i, (g, d) in enumerate(zip(REGIONS, data)):
        show = d if len(d) <= jitter_max else rng.choice(d, jitter_max,
                                                         replace=False)
        n_skip[g] = len(d) - len(show)
        n_edge[g] = int(((show < lo) | (show > hi)).sum())
        ax.scatter(i + rng.normal(0, 0.06, len(show)), np.clip(show, lo, hi),
                   s=1.4, linewidths=0, color=REGION[g], alpha=0.35, zorder=3,
                   rasterized=True)
    span = hi - lo
    top = hi + 0.06 * span
    ax.plot([0, 0, 1, 1], [top, top + 0.03 * span, top + 0.03 * span, top],
            color='black', lw=LINEWIDTH['thin'], zorder=6)
    ax.text(0.5, top + 0.045 * span, stars(test['q']), ha='center',
            va='bottom', fontsize=FONT['annotation'])
    draw_p(ax, 0.5, top + 0.16 * span, test['p'], ha='center', va='bottom',
           fontsize=FONT['annotation'])
    ax.set_ylim(lo - 0.04 * span, top + 0.30 * span)
    ax.set_xlim(-0.55, 1.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Retained\n(n = {len(data[0]):,})',
                        f'Eliminated\n(n = {len(data[1]):,})'])
    ax.tick_params(axis='x', length=0)
    ax.set_ylabel(ylabel)
    return {'y_range': (float(lo), float(hi)), 'range_pct': list(range_pct),
            'jitter_not_drawn': n_skip, 'jitter_on_edge': n_edge,
            'jitter_max': jitter_max}


def draw_scatter(ax, df, rho, p, cap_pct):
    """C: dN/dS against GC3s, both regions, eliminated drawn last.

    The y-axis stops at the `cap_pct` percentile of the plotted dN/dS; genes
    above are drawn as open triangles on the top edge and counted in the key,
    so a clipped point is never mistaken for an absent one.
    """
    sub = df.dropna(subset=['GC3s', 'dN_dS'])
    ycap = float(np.percentile(sub['dN_dS'], cap_pct))
    # A margin at each end of the y-axis, as in S6. With the axis at exactly
    # (0, cap), a gene at dN/dS = 0 sits on the x-axis and is cut in half by
    # it, and a gene just under the cap is cut by the top of the axes. The
    # triangles for genes above the cap go on the new top edge.
    y_lo, y_hi = -0.02 * ycap, 1.03 * ycap
    n_above, n_by = {}, {}
    for reg in REGIONS:                       # eliminated drawn on top
        r = sub[sub['region'] == reg]
        n_by[reg] = len(r)
        inside, above = r[r['dN_dS'] <= ycap], r[r['dN_dS'] > ycap]
        n_above[reg] = len(above)
        elim = reg == 'eliminated'
        ax.scatter(inside['GC3s'], inside['dN_dS'], s=3.0 if elim else 2.0,
                   linewidths=0,
                   color=to_rgba(REGION[reg], 0.85 if elim else 0.35),
                   zorder=3 if elim else 2, rasterized=True)
        if len(above):
            ax.scatter(above['GC3s'], np.full(len(above), y_hi), marker='^',
                       s=6, facecolors='none', edgecolors=REGION[reg],
                       linewidths=0.4, zorder=4, clip_on=False)
    x0, x1 = sub['GC3s'].min(), sub['GC3s'].max()
    pad = 0.03 * (x1 - x0 if x1 > x0 else 1.0)
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel('GC3s')
    ax.set_ylabel('dN/dS')

    # Key above the axes, two rows: the two regions, then the open-triangle
    # marker with its count -- a blue and a red triangle side by side, as
    # drawn (one black triangle matched neither). The Spearman result sits
    # right-aligned on the second row; place_spearman() moves it up a row if
    # it would collide.
    handles = [Line2D([], [], ls='none', marker='o', ms=3.2,
                      mfc=REGION[g], mec=REGION[g]) for g in REGIONS]
    tri = tuple(Line2D([], [], ls='none', marker='^', ms=3.2, mfc='none',
                       mec=REGION[g], mew=0.5) for g in REGIONS)
    blank = Line2D([], [], ls='none')
    labels = [f'Retained (n = {n_by["retained"]:,})',
              f'Above the axis ({sum(n_above.values()):,})',
              f'Eliminated (n = {n_by["eliminated"]:,})', '']
    leg = ax.legend([handles[0], tri, handles[1], blank], labels,
                    loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2,
                    frameon=False, fontsize=FONT['legend'], borderaxespad=0.3,
                    handletextpad=0.3, columnspacing=1.2, labelspacing=0.3,
                    handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
    # Plain text, the exponent raised as its own 6 pt text: a list of texts.
    txts = draw_p(ax, 1.0, 1.0, p, prefix=f'Spearman ρ = {signed(rho, 3)}, ',
                  transform=ax.transAxes, ha='right', va='bottom', dy_pt=3.0,
                  fontsize=FONT['annotation'], color=NEUTRAL['ink'])
    return leg, txts, {'y_cap': ycap, 'cap_pct': cap_pct, 'n_above': n_above,
                       'n_plotted': n_by}


def place_spearman(fig, ax, leg, txts):
    """Raise the Spearman texts above the key if they overlap a key label.

    `txts` is draw_p's list (the text and its raised exponent); they move
    together, by the same offset added to each one's own transform, so the
    junction between them is kept. Call once.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(rend) for t in txts]
    x0, x1 = min(b.x0 for b in boxes), max(b.x1 for b in boxes)
    y0, y1 = min(b.y0 for b in boxes), max(b.y1 for b in boxes)
    for t in leg.get_texts():
        if not t.get_text():
            continue
        b = t.get_window_extent(rend)
        if min(x1, b.x1) > max(x0, b.x0) and min(y1, b.y1) > max(y0, b.y0):
            lb = leg.get_window_extent(rend)
            dy_pt = (lb.y1 - y0) / fig.dpi * 72 + 1.5
            for tx in txts:
                tx.set_transform(tx.get_transform() + ScaledTranslation(
                    0, dy_pt / 72, fig.dpi_scale_trans))
            return True
    return False


def plot_s1(use, fin, tests, rho_all, outdir, stem, rng, args):
    """Supplementary fig. S1 as one 185 mm strip: A GC3s, B GC, C dN/dS vs
    GC3s.

    `use` = every gene passing the codon filter (A, B and their tests);
    `fin` = those with a finite dN/dS (C and the Spearman test).
    """
    fig = blank_figure('double', args.height_mm)
    w_in, h_in = fig.get_size_inches()
    ax_a = fig.add_axes([0.06, 0.2, 0.15, 0.6])
    ax_b = fig.add_axes([0.30, 0.2, 0.15, 0.6])
    ax_c = fig.add_axes([0.55, 0.2, 0.40, 0.6])
    info = {}
    info['A'] = draw_box(ax_a, use, 'GC3s', 'GC3s', tests['GC3s'], rng,
                         args.jitter_max, args.box_range_pct)
    info['B'] = draw_box(ax_b, use, 'GC', 'GC content (coding sequence)',
                         tests['GC'], rng, args.jitter_max, args.box_range_pct)
    leg, txts, info['C'] = draw_scatter(ax_c, fin, rho_all['rho'], rho_all['p'],
                                        args.dnds_cap_pct)
    top_extra = 4.5 if args.panel_letters else 0.0
    place_row(fig, [ax_a, ax_b, ax_c], shares=[0.235, 0.235, 0.53],
              gaps_mm=[5.0, 7.0], top_extra_mm=top_extra)
    if place_spearman(fig, ax_c, leg, txts):
        logging.info('  Spearman text raised above the key to clear it')
        place_row(fig, [ax_a, ax_b, ax_c], shares=[0.235, 0.235, 0.53],
                  gaps_mm=[5.0, 7.0], top_extra_mm=top_extra)
    if args.panel_letters:
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        y = 1.0 - mm(0.8) / h_in
        for ax, letter in ((ax_a, 'A'), (ax_b, 'B'), (ax_c, 'C')):
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
    t3, tgc = ctx['tests']['GC3s'], ctx['tests']['GC']
    sp = {r['subset']: r for r in ctx['spearman']}
    m, um = ctx['matched'], ctx['unmatched']
    L = ['Base composition of eliminated and retained genes (supplementary '
         'fig. S2)', '=' * 72, '',
         f'Rate table   {a.annotated}',
         f'CDS FASTA    {a.cds}  ({ctx["n_fasta"]:,} records)',
         f"Gene set     {ctx['gene_set']['source']}", '',
         'Genes at each step                               retained  eliminated']
    for k, (r, e) in ctx['steps'].items():
        L.append(f'  {k:<46} {r:>8,}  {e:>10,}')
    gs = ctx['gene_set']
    if gs['codon_only'] is not None:
        L.append(f"  The Fig 2 set includes {gs['codon_only'][0]:,} retained and "
                 f"{gs['codon_only'][1]:,} eliminated genes with fewer than "
                 f'{a.min_codons} informative codons (sense codons other than '
                 'ATG and TGG),')
        L.append(f"  and omits {gs['own_only'][0]:,} / {gs['own_only'][1]:,} "
                 'that pass that rule. Methods describes the codon filter as '
                 f'>= {a.min_codons} informative codons excluding Met, Trp and '
                 'stops;')
        if gs['codon_only'][0] or gs['codon_only'][1]:
            L.append('  that does not match the set the codon script kept. Check '
                     'its filter before Methods is final.')
        else:
            L.append('  the two agree.')
    L.append('')

    L += ['-' * 72, '1. Composition, eliminated vs retained (Mann-Whitney; BH q '
          'over the two)', '-' * 72,
          f"{'metric':<6} {'median ret':>10} {'median elim':>11} {'diff':>8} "
          f"{'P':>9} {'q':>9}"]
    for t in (t3, tgc):
        L.append(f"{t['metric']:<6} {t['median_retained']:>10.4f} "
                 f"{t['median_eliminated']:>11.4f} {t['difference']:>+8.4f} "
                 f"{p_text(t['p']):>9} {p_text(t['q']):>9}")
    L += [f"  IQR, GC3s: retained {t3['q1_retained']:.4f}-{t3['q3_retained']:.4f},"
          f" eliminated {t3['q1_eliminated']:.4f}-{t3['q3_eliminated']:.4f}", '']

    L += ['-' * 72, '2. dN/dS vs GC3s (Spearman)', '-' * 72]
    for name in ('all genes', 'retained', 'eliminated'):
        r = sp[name]
        L.append(f"  {name:<11} rho = {r['rho']:+.4f}   P = {p_text(r['p']):<9} "
                 f"n = {r['n']:,}")
    L += ['  The pooled coefficient includes the region effect (eliminated '
          'genes are both', '  GC3s-richer and faster-evolving), which pulls '
          'it toward zero; the within-retained', '  coefficient is free of it.',
          '']

    L += ['-' * 72, f'3. GC3s-matched dN/dS contrast (caliper {a.caliper}, '
          f'{a.n_matched} replicates, seed {a.seed})', '-' * 72,
          f"  Unmatched, same genes: eliminated/retained median dN/dS "
          f"{um['fold']:.3f}-fold (P = {p_text(um['p'])}; n = "
          f"{um['n_eliminated']:,} / {um['n_retained']:,})"]
    if m is not None and len(m):
        L += [f"  Matched: {m['fold'].median():.3f}-fold (2.5-97.5% of "
              f"replicates {m['fold'].quantile(0.025):.3f}-"
              f"{m['fold'].quantile(0.975):.3f}); median P = "
              f"{p_text(m['p'].median())}; P < 0.05 in "
              f"{100 * (m['p'] < 0.05).mean():.0f}% of replicates",
              f"  Eliminated genes matched per replicate: median "
              f"{int(m['n_matched'].median()):,} of {um['n_eliminated']:,} "
              f"(without any retained gene inside the caliper: "
              f"{ctx['n_no_window']})",
              f"  GC3s balance after matching: median Mann-Whitney P = "
              f"{p_text(m['gc3s_balance_p'].median())} (a working match "
              'leaves no GC3s difference)', '']
    else:
        L += ['  Not run (--n-matched 0, or fewer than 5 eliminated genes).', '']

    cc = ctx['cross_check']
    L += ['-' * 72, '4. Cross-check against codon_usage_bias.py per-gene GC3s',
          '-' * 72]
    if cc is None:
        L.append('  Not run (no --codon-per-gene table, or it had no usable '
                 'columns).')
    elif cc['n_shared'] == 0:
        L.append(f"  {cc['path']}: no gene in common.")
    else:
        verdict = ('identical' if cc['n_differ'] == 0 else
                   f"DIFFERENT for {cc['n_differ']:,} genes (max |diff| "
                   f"{cc['max_abs_diff']:.4g}; e.g. {', '.join(cc['examples'])})")
        L.append(f"  {cc['path']} (column {cc['column']}): {cc['n_shared']:,} "
                 f'genes shared; GC3s {verdict}.')
        if cc['n_differ']:
            L.append('  S1 and Figure 2 would quote different GC3s for these '
                     'genes: find out why before either is used.')
    L.append('')

    q = ctx['qc']
    L += ['-' * 72, '5. Gene-model QC (analysed genes with a CDS record)', '-' * 72,
          f"  internal stop codon(s)      {q['internal_stops']:,}",
          f"  ambiguous codon(s)          {q['ambiguous']:,}",
          f"  length not a multiple of 3  {q['trailing']:,}",
          f"  no terminal stop codon      {q['no_terminal_stop']:,}",
          '  Internal stops and ambiguous codons are left out of GC3s; such genes '
          'are kept.', '']

    fi = ctx['figure']
    L += ['-' * 72, 'Drawn (figS1_composition)', '-' * 72]
    for key, name in (('A', 'GC3s'), ('B', 'GC')):
        d = fi[key]
        L.append(f"  {key} ({name}): y-range = central {d['range_pct'][0]:g}-"
                 f"{d['range_pct'][1]:g} percentiles; jitter on the edge: "
                 f"retained {d['jitter_on_edge']['retained']:,}, eliminated "
                 f"{d['jitter_on_edge']['eliminated']:,}; jitter not drawn: "
                 f"retained {d['jitter_not_drawn']['retained']:,}, eliminated "
                 f"{d['jitter_not_drawn']['eliminated']:,} (at most "
                 f"{d['jitter_max']:,} per region; boxes and tests use every "
                 'gene)')
    c = fi['C']
    L += [f"  C: y-axis stops at the {c['cap_pct']:g}th percentile of dN/dS "
          f"({c['y_cap']:.3f}), plus a 3% margin (2% below 0) so no point is cut "
          f"by an axis; above the cap, drawn as open triangles on the top edge: "
          f"retained {c['n_above']['retained']:,}, eliminated "
          f"{c['n_above']['eliminated']:,}.",
          '  Box = quartiles, whiskers = 5th-95th percentiles. Asterisks from '
          'BH q over panels A and B; the exact P is printed.']
    if ctx['problems']:
        L.append('  AUDIT: ' + '; '.join(ctx['problems']))
    else:
        L.append('  audit: clean')
    L.append('')

    # --- what the text can claim ------------------------------------------
    L += ['-' * 72, 'What the text can claim', '-' * 72]
    higher = t3['difference'] > 0
    neg = sp['all genes']['rho'] < 0 and sp['all genes']['p'] < 0.05
    if higher and t3['p'] < 0.05 and neg:
        L.append('  Replacement for the second control in "Eliminated genes '
                 'evolve under relaxed purifying selection":')
        L.append(f"    Second, eliminated genes have higher GC content at "
                 f"synonymous third positions (median GC3s "
                 f"{t3['median_eliminated']:.3f} versus "
                 f"{t3['median_retained']:.3f}; P = "
                 f"{p_manuscript(t3['p'])}), but GC3s correlates negatively "
                 f"with dN/dS across all genes (Spearman ρ = "
                 f"{signed(sp['all genes']['rho'], 3)}; P = "
                 f"{p_manuscript(sp['all genes']['p'])}), so the compositional "
                 'difference biases the observed contrast downward '
                 '(supplementary fig. S2, Supplementary Material online).')
        if m is not None and len(m):
            if m['fold'].median() >= um['fold']:
                L.append(f"  Direct test agrees: matched on GC3s, the contrast is "
                         f"{m['fold'].median():.2f}-fold against "
                         f"{um['fold']:.2f}-fold unmatched. Optional clause: "
                         f"\"; matched on GC3s, eliminated genes still have "
                         f"{m['fold'].median():.2f}-fold higher dN/dS\".")
            else:
                L.append(f"  CAUTION: matched on GC3s the contrast is SMALLER "
                         f"({m['fold'].median():.2f}- vs {um['fold']:.2f}-fold). "
                         'Do not claim the composition effect is downward; '
                         'report the matched fold instead.')
    else:
        L.append('  The premise does not hold on these numbers (GC3s not higher '
                 'in eliminated genes, or no negative correlation): rewrite '
                 'the sentence from the tables above.')
    L.append(f"  GC3s here must equal Results' codon paragraph (GC3s "
             f"{t3['median_eliminated']:.3f} versus "
             f"{t3['median_retained']:.3f}; P = {p_manuscript(t3['p'])}); see "
             'section 4.')
    L.append('')
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='GC3s and GC of eliminated vs retained Ascaris genes, and '
                    'dN/dS against GC3s (supplementary fig. S2).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', required=True,
                   help='kaks_annotated.tsv from annotate_regions.py '
                        '(ANNOTATED_TABLE in config.sh)')
    p.add_argument('--cds', required=True,
                   help='Ascaris CDS FASTA keyed like the gene column '
                        '(As_cds.fa; .gz accepted)')
    p.add_argument('--codon-per-gene', default=None,
                   help='Optional codon_usage_per_gene.tsv from '
                        'codon_usage_bias.py, for a per-gene GC3s cross-check')
    p.add_argument('--own-filter', action='store_true',
                   help='Use S1\'s own --min-codons rule even when '
                        '--codon-per-gene is given (default: adopt the codon '
                        'table\'s gene set, so S1 and Fig 2 test the same genes)')
    p.add_argument('--gene-col', default='gene1',
                   help='Ascaris gene column (default: gene1)')
    p.add_argument('--region-col', default='region1',
                   help='Ascaris region column (default: region1)')
    p.add_argument('--min-codons', type=int, default=100,
                   help='Minimum informative codons per gene (sense codons '
                        'other than ATG and TGG; default 100, as in '
                        'codon_usage_bias.py)')
    p.add_argument('--caliper', type=float, default=0.01,
                   help='GC3s caliper for the matched contrast (default 0.01, '
                        'as in the codon-usage control)')
    p.add_argument('--n-matched', type=int, default=500,
                   help='Matched replicates (default 500)')
    p.add_argument('--dnds-cap-pct', type=float, default=99.5,
                   help='Panel C y-axis stops at this percentile; points above '
                        'are drawn on the edge and counted (default 99.5)')
    p.add_argument('--box-range-pct', type=float, nargs=2, default=(0.25, 99.75),
                   metavar=('LO', 'HI'),
                   help='Panels A and B y-range, as pooled percentiles '
                        '(default 0.25 99.75); drawing only')
    p.add_argument('--jitter-max', type=int, default=1500,
                   help='At most this many jittered points per region in A and '
                        'B; boxes and tests use every gene (default 1500)')
    p.add_argument('--height-mm', type=float, default=62.0,
                   help='Figure height (width is the 185 mm double column)')
    p.add_argument('--panel-letters', action='store_true',
                   help='Draw A, B, C (off by default, per the house rule; S1 '
                        'is emitted whole, so its wrapper turns them on)')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.min_codons < 1 or args.n_matched < 0 or args.caliper <= 0:
        p.error('--min-codons must be >= 1, --n-matched >= 0, --caliper > 0')
    if not (0 < args.dnds_cap_pct <= 100):
        p.error('--dnds-cap-pct must be in (0, 100]')
    lo, hi = args.box_range_pct
    if not (0 <= lo < hi <= 100):
        p.error('--box-range-pct needs 0 <= LO < HI <= 100')
    if args.outdir is None:
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.outdir = os.path.join('results', f'compositional_analysis_{stamp}')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'compositional_analysis')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('compositional_analysis.py -- supplementary fig. S2')
    logging.info('  rate table: %s', os.path.abspath(args.annotated))
    logging.info('  CDS:        %s', os.path.abspath(args.cds))
    logging.info('  outdir:     %s', os.path.abspath(args.outdir))
    for f in (args.annotated, args.cds):
        if not os.path.isfile(f):
            sys.exit(f'ERROR: not found: {f}')
    if args.codon_per_gene and not os.path.isfile(args.codon_per_gene):
        sys.exit(f'ERROR: --codon-per-gene not found: {args.codon_per_gene}')
    apply_style()
    rng = np.random.default_rng(args.seed)

    # --- rate table ---------------------------------------------------------
    rates = read_rate_table(args.annotated, require=('dN_dS',))
    for col in (args.gene_col, args.region_col):
        if col not in rates.columns:
            sys.exit(f'ERROR: column {col!r} not in {args.annotated}. '
                     f'Columns: {list(rates.columns)}')
    rates = rates.rename(columns={args.gene_col: 'gene1',
                                  args.region_col: 'region'})
    rates['gene1'] = rates['gene1'].astype(str)
    for c in ('dN', 'dS', 'dN_dS'):
        if c in rates.columns:
            rates[c] = pd.to_numeric(rates[c], errors='coerce')
    n_rows = len(rates)
    df = rates[rates['region'].isin(REGIONS)].copy()
    logging.info('Rate table: %s rows, %s retained or eliminated',
                 f'{n_rows:,}', f'{len(df):,}')
    if df.empty:
        sys.exit('ERROR: no row is labelled retained or eliminated in '
                 f'{args.region_col!r}. Values: '
                 f"{rates['region'].value_counts().head().to_dict()}")
    n_dup = int(df['gene1'].duplicated().sum())
    if n_dup:
        logging.warning('  %d rows repeat an Ascaris gene already in the table; '
                        'each row is analysed as given', n_dup)

    def counts(mask):
        return (int((mask & (df['region'] == 'retained')).sum()),
                int((mask & (df['region'] == 'eliminated')).sum()))

    steps = {'analysed (retained or eliminated)': counts(df['gene1'].notna())}

    # --- sequences ------------------------------------------------------------
    wanted = set(df['gene1'])
    seqs, repeats, n_fasta = read_fasta(args.cds, wanted)
    logging.info('CDS FASTA: %s records; %s of %s analysed genes found',
                 f'{n_fasta:,}', f'{len(seqs):,}', f'{len(wanted):,}')
    bad = sorted(k for k, same in repeats.items() if not same)
    if bad:
        sys.exit(f'ERROR: {len(bad)} analysed genes have more than one CDS '
                 f'record with different sequences in {args.cds}, e.g. '
                 f'{bad[:5]}. Resolve the duplicates first.')
    if repeats:
        logging.warning('  %d analysed genes appear twice in the FASTA with '
                        'identical sequence; one copy used', len(repeats))
    if len(seqs) < 0.5 * len(wanted):
        sys.exit('ERROR: fewer than half the analysed genes have a CDS record. '
                 'Is the FASTA keyed like the gene column?\n'
                 f"  table IDs: {sorted(wanted)[:3]}\n"
                 f"  FASTA IDs found: {sorted(seqs)[:3]}")

    comp = pd.DataFrame.from_dict({g: composition(s) for g, s in seqs.items()},
                                  orient='index')
    comp.index.name = 'gene1'
    df = df.merge(comp, left_on='gene1', right_index=True, how='left')
    has_cds = df['n_codons'].notna()
    steps['with a CDS record'] = counts(has_cds)
    passes = has_cds & (df['n_informative'] >= args.min_codons)
    steps[f'>= {args.min_codons} informative codons (S1 rule)'] = counts(passes)
    # Gene set: Figure 2's, when the codon table is given (see codon_gene_set),
    # else S1's own minimum-codon rule.
    analysed = passes
    gene_set = {'source': f'S1 rule (>= {args.min_codons} informative codons)',
                'own_only': None, 'codon_only': None}
    cg = (codon_gene_set(args.codon_per_gene)
          if args.codon_per_gene and not args.own_filter else None)
    if cg is not None:
        in_codon = df['gene1'].isin(cg[0])
        analysed = has_cds & in_codon
        steps['in the codon_usage_bias.py gene set (Fig 2)'] = counts(analysed)
        gene_set = {'source': 'codon_usage_bias.py per-gene table (Fig 2 gene '
                              'set)',
                    'own_only': counts(passes & ~in_codon),
                    'codon_only': counts(analysed & ~passes)}
        logging.info('  gene set = Fig 2: %s retained, %s eliminated (S1 rule '
                     'alone would drop %s / %s of them and add %s / %s)',
                     *map(lambda v: f'{v:,}', steps[
                         'in the codon_usage_bias.py gene set (Fig 2)']),
                     *map(lambda v: f'{v:,}', gene_set['codon_only']),
                     *map(lambda v: f'{v:,}', gene_set['own_only']))
    qc = {'internal_stops': int((df.loc[has_cds, 'internal_stops'] > 0).sum()),
          'ambiguous': int((df.loc[has_cds, 'ambiguous_codons'] > 0).sum()),
          'trailing': int((df.loc[has_cds, 'trailing_bases'] > 0).sum()),
          'no_terminal_stop': int((~df.loc[has_cds, 'terminal_stop']
                                   .astype(bool)).sum())}
    logging.info('  QC: internal stops %s, ambiguous codons %s, length not a '
                 'multiple of 3 %s, no terminal stop %s', *map(str, qc.values()))

    out_cols = (['gene1'] + [c for c in ('gene2',) if c in df.columns] +
                ['region'] + [c for c in ('chrom1',) if c in df.columns] +
                [c for c in ('dN', 'dS', 'dN_dS') if c in df.columns] +
                ['GC3s', 'GC', 'n_codons', 'n_informative', 'n_met_trp',
                 'internal_stops', 'ambiguous_codons', 'terminal_stop',
                 'trailing_bases'])
    per_gene = df[out_cols].assign(passes_min_codons=passes.to_numpy(),
                                   analysed=analysed.to_numpy())
    per_gene.rename(columns={'region': 'region1'}).to_csv(
        os.path.join(args.outdir, 'composition_per_gene.tsv'), sep='\t',
        index=False, float_format='%.6g')

    use = df[analysed].copy()
    steps['analysed, finite dN/dS (panel C, tests 2-3)'] = counts(
        analysed & df['dN_dS'].notna() & np.isfinite(df['dN_dS']))
    for k, (r, e) in steps.items():
        logging.info('  %-46s retained %s, eliminated %s', k, f'{r:,}', f'{e:,}')
    if min(counts(analysed)) < 5:
        sys.exit('ERROR: fewer than 5 analysed genes in a region.')

    # --- tests ----------------------------------------------------------------
    tests = {c: region_test(use, c) for c in ('GC3s', 'GC')}
    qs = benjamini_hochberg([tests['GC3s']['p'], tests['GC']['p']])
    tests['GC3s']['q'], tests['GC']['q'] = float(qs[0]), float(qs[1])
    for c in ('GC3s', 'GC'):
        t = tests[c]
        logging.info('  %-4s retained %.4f, eliminated %.4f, P %s, q %s', c,
                     t['median_retained'], t['median_eliminated'],
                     p_text(t['p']), p_text(t['q']))

    fin = use[np.isfinite(use['dN_dS'])].copy()
    spear = spearman_block(fin, 'GC3s', 'dN_dS')
    for r in spear:
        logging.info('  dN/dS vs GC3s, %-10s rho %+.4f, P %s, n %s', r['subset'],
                     r['rho'], p_text(r['p']), f"{r['n']:,}")

    unmatched = region_test(fin, 'dN_dS')
    r_mask = fin['region'] == 'retained'
    matched, n_no_window = None, 0
    if args.n_matched and unmatched['n_eliminated'] >= 5:
        matched, n_no_window = matched_contrast(
            fin.loc[r_mask, 'GC3s'].to_numpy(float),
            fin.loc[r_mask, 'dN_dS'].to_numpy(float),
            fin.loc[~r_mask, 'GC3s'].to_numpy(float),
            fin.loc[~r_mask, 'dN_dS'].to_numpy(float),
            args.caliper, args.n_matched, args.seed + 1)
        matched.to_csv(os.path.join(args.outdir, 'matched_replicates.tsv'),
                       sep='\t', index=False, float_format='%.6g')
        logging.info('  GC3s-matched dN/dS: %.3f-fold (unmatched %.3f), median '
                     'P %s, balance P %s', matched['fold'].median(),
                     unmatched['fold'], p_text(matched['p'].median()),
                     p_text(matched['gc3s_balance_p'].median()))

    rows = []
    for c in ('GC3s', 'GC'):
        t = tests[c]
        rows.append({'test': f'{c} eliminated vs retained', 'family':
                     'composition', 'n_retained': t['n_retained'],
                     'n_eliminated': t['n_eliminated'],
                     'median_retained': t['median_retained'],
                     'median_eliminated': t['median_eliminated'],
                     'effect': t['difference'], 'effect_type': 'difference',
                     'p': t['p'], 'q': t['q'], 'method': 'Mann-Whitney'})
    for r in spear:
        rows.append({'test': f"dN/dS vs GC3s, {r['subset']}", 'family':
                     'spearman', 'n_retained': np.nan, 'n_eliminated': np.nan,
                     'median_retained': np.nan, 'median_eliminated': np.nan,
                     'effect': r['rho'], 'effect_type': 'rho', 'p': r['p'],
                     'q': np.nan, 'method': f"Spearman, n = {r['n']}"})
    rows.append({'test': 'dN/dS eliminated vs retained, unmatched',
                 'family': 'dnds', 'n_retained': unmatched['n_retained'],
                 'n_eliminated': unmatched['n_eliminated'],
                 'median_retained': unmatched['median_retained'],
                 'median_eliminated': unmatched['median_eliminated'],
                 'effect': unmatched['fold'], 'effect_type': 'fold',
                 'p': unmatched['p'], 'q': np.nan, 'method': 'Mann-Whitney'})
    if matched is not None and len(matched):
        rows.append({'test': 'dN/dS eliminated vs retained, GC3s-matched',
                     'family': 'dnds',
                     'n_retained': int(matched['n_matched'].median()),
                     'n_eliminated': int(matched['n_matched'].median()),
                     'median_retained': matched['median_matched_retained'].median(),
                     'median_eliminated': matched['median_eliminated'].median(),
                     'effect': matched['fold'].median(), 'effect_type':
                     'fold (median over replicates)',
                     'p': matched['p'].median(), 'q': np.nan,
                     'method': f'Mann-Whitney per replicate, caliper '
                               f'{args.caliper}, {args.n_matched} replicates'})
    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, 'composition_tests.tsv'),
                              sep='\t', index=False, float_format='%.6g')

    cc = (cross_check(per_gene[analysed.to_numpy()], args.codon_per_gene)
          if args.codon_per_gene else None)

    # --- figure ---------------------------------------------------------------
    logging.info('S2: composition strip')
    rho_all = next(r for r in spear if r['subset'] == 'all genes')
    if len(fin) < 10:
        sys.exit('ERROR: fewer than 10 genes with a finite dN/dS.')
    paths, problems, fig_info = plot_s1(use, fin, tests, rho_all, args.outdir,
                                        'figS1_composition', rng, args)

    write_summary(os.path.join(args.outdir, 'composition_summary.txt'), {
        'args': args, 'n_fasta': n_fasta, 'steps': steps, 'tests': tests,
        'spearman': spear, 'matched': matched, 'unmatched': unmatched,
        'n_no_window': n_no_window, 'cross_check': cc, 'qc': qc,
        'gene_set': gene_set,
        'figure': fig_info, 'problems': problems})

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
                   'inputs': {'annotated': fingerprint(args.annotated),
                              'cds': fingerprint(args.cds),
                              'codon_per_gene': fingerprint(args.codon_per_gene)},
                   'n_rows': int(n_rows),
                   'steps': {k: list(v) for k, v in steps.items()},
                   'figure': [os.path.basename(p) for p in paths],
                   'audit': problems},
                  fh, indent=2, default=str)
    logging.info('Summary: %s', os.path.join(args.outdir,
                                             'composition_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
