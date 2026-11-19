#!/usr/bin/env python3
"""
kaks_compartments.py -- dN, dS and dN/dS of retained Ascaris genes by somatic
A/B compartment (Supplementary Figure S13A; the file keeps its name from before
the renumbering, figS9A_compartments).

A post-elimination (somatic) Hi-C eigenvector resolves compartments far better
than any pre-elimination stage, but it is defined only on retained DNA. This
script therefore asks the question that eigenvector can answer -- among
RETAINED genes, does constraint differ between A and B compartments? -- and
leaves the eliminated-vs-retained contrast to the pre-elimination eigenvector
controls in multivariate_insulation.py (table S9).

Method
------
  1. Positions are the germline gene midpoints in the annotated table
     (midpoint1, from each model's best alignment).
  2. Each retained gene is moved to somatic coordinates with the germline-to-
     somatic map (--coord-map): the retained segment containing its midpoint
     gives the somatic chromosome and the offset.
  3. EV1 is the value of the somatic 100 kb bin containing that position.
  4. EV1 is ORIENTED PER SOMATIC CHROMOSOME so that it correlates positively
     with local GC (G+C in --gc-window around each gene, from the germline
     assembly). FAN-C's sign can invert between chromosomes -- four of 21 were
     inverted in the testis eigenvector -- and the May script took the raw sign,
     so A and B could be swapped on some chromosomes and the A-vs-B comparison
     diluted towards nothing. A = EV1 > 0 after orientation; EV1 = 0 is
     unassigned.
  5. Tests on retained genes: Mann-Whitney A vs B for dN/dS, dN and dS (BH q
     over the three), and Spearman of each against EV1.

Checks written to the summary
-----------------------------
  * Eliminated genes that received a somatic position: should be 0; any is a
    sign that the coordinate map covers eliminated sequence.
  * Coordinate-map format and consistency: segments per somatic chromosome,
    and each segment's length against the eigenvector track's extent. A
    somatic chromosome built from more than one germline segment cannot be
    placed from a 4-column map, and is refused.
  * Per-chromosome Spearman(EV1, GC) before and after orientation. A
    chromosome still near zero after orientation is listed: a mirrored
    coordinate conversion (somatic chromosome reverse-complemented relative to
    germline, with no strand column) looks like that.
  * With --soma-genome, the strand and offset of every somatic chromosome are
    read from the two assemblies (unique 60-mer probes) and the map is
    corrected where it disagrees. This settles what the correlation heuristic
    can only suggest: the first real run (2026-10-12) flagged seven somatic
    chromosomes whose EV1 tracked GC better read backwards.

What changed from the May script
--------------------------------
  * EV1 orientation by GC (above); the May A/B labels used the raw sign.
  * Positions from the annotated table, not the first GFF record per model.
  * Midpoints are 1-based and the bins 0-based half-open; the lookup converts.
  * Removed: the eliminated-gene compartment enrichment and the compartment x
    region panels -- with a somatic eigenvector no eliminated gene has a
    compartment, so both were empty comparisons.
  * Colour: compartments are a subclass within retained genes, so both boxes
    take the retained blue, and the compartment is named on the axis. The May
    figure drew A red and B blue, which reads as eliminated vs retained
    everywhere else in the paper. B was hatched until 2026-10-19, when the
    hatch was removed by hand in the compiled S13A; now in the code.
  * House figure pipeline; nothing capped silently (the May script dropped
    values above the 99th percentile before drawing).
  * 2026-10-19 (results draft note 67): genes above the y-axis limit are open
    triangles at the limit, not dots. As dots they made a flat row at the cap
    (dN/dS 0.99, dN 0.48, dS 1.4), which sits mid-axis because the axis
    extends above it for the bracket, and read as real values: the 6D defect
    of note 31. The axis starts 2% below 0 so a gene at 0 is not cut by it,
    and P is plain text with its exponent raised as its own 6 pt text
    (plain_text.draw_p), not mathtext.

Outputs (in --outdir)
---------------------
  figS9A_compartments.{svg,pdf,png}  S13A: three box panels, 185 mm strip
  compartment_per_gene.tsv           every table row, with somatic position,
                                     EV1 (raw and oriented) and compartment
  compartment_tests.tsv
  ev1_orientation.tsv                per somatic chromosome, before and after
  coord_map_check.tsv                per somatic chromosome
  compartment_summary.txt
  kaks_compartments.log
  provenance.json

Usage
-----
  python kaks_compartments.py \\
    --annotated results/annotate_regions_20260817/kaks_annotated.tsv \\
    --eigenvector .../postpde/eigenvectors/as_10day_iced_100kb.matrix.eigenvector \\
    --coord-map .../v50_postPDE/prepde_to_postpde_coordinates.bed \\
    --genome AG.v50.fa --ev-label 'somatic (10 day)' \\
    --outdir results/kaks_compartments_<stamp>

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
import re
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from scipy.stats import mannwhitneyu, spearmanr
from scipy import __version__ as scipy_version

# A script that cannot import the shared style must stop, not fall back to
# matplotlib defaults.
from paper_colors import (REGION, FONT, LINEWIDTH, apply_style, audit, figure,
                          save, mm)
from region_lib import (benjamini_hochberg, fingerprint, git_commit,
                        setup_logging, stars)
from rate_tables import read_rate_table
from plain_text import draw_p

_HERE = os.path.dirname(os.path.abspath(__file__))

METRICS = ('dN_dS', 'dN', 'dS')
METRIC_LABEL = {'dN_dS': 'dN/dS', 'dN': 'dN', 'dS': 'dS'}
COMPARTMENTS = ('A', 'B')


# ---------------------------------------------------------------------------
# Inputs (track and genome helpers as in multivariate_insulation.py)
# ---------------------------------------------------------------------------

def norm_chrom(c):
    """Strip zero-padding from numeric runs: chr06 -> chr6, chrX1 stays chrX1.

    The table, the tracks, the map and the FASTA come from different pipelines
    and do not agree on padding; every source passes through this.
    """
    c = str(c).strip()
    body = c[3:] if c.startswith('chr') else c
    return 'chr' + re.sub(r'(\d+)', lambda m: str(int(m.group(1))), body)


def load_track(path, value_name):
    """A 4-column bedGraph-style track: chrom, start, end, value.

    UCSC `track`/`browser`/`#` lines are skipped; NA, NaN and '.' values are
    dropped. Coordinates stay 0-based half-open, as in the file.
    """
    rows, skipped = [], 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith(('track', 'browser', '#')):
                skipped += 1
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 4:
                skipped += 1
                continue
            rows.append((norm_chrom(p[0]), p[1], p[2], p[3]))
    if not rows:
        sys.exit(f'ERROR: no data rows in {path}')
    df = pd.DataFrame(rows, columns=['chrom', 'start', 'end', value_name])
    for c in ('start', 'end', value_name):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    n_all = len(df)
    df = df.dropna(subset=['start', 'end', value_name]).copy()
    df[['start', 'end']] = df[['start', 'end']].astype(np.int64)
    logging.info('  %s: %s bins with a value (%s without; %d header lines)',
                 os.path.basename(path), f'{len(df):,}', f'{n_all - len(df):,}',
                 skipped)
    return df, n_all


def lookup_track(chrom, pos0, track, value_name):
    """Value of the bin containing each 0-based position, or NaN."""
    by_chrom = {c: (g['start'].to_numpy(), g['end'].to_numpy(),
                    g[value_name].to_numpy())
                for c, g in track.sort_values(['chrom', 'start'])
                .groupby('chrom', sort=False)}
    out = np.full(len(chrom), np.nan)
    n_nochrom = n_gap = 0
    for i, (c, p) in enumerate(zip(chrom, pos0)):
        if not isinstance(c, str) or not np.isfinite(p):
            continue
        if c not in by_chrom:
            n_nochrom += 1
            continue
        s, e, v = by_chrom[c]
        j = int(np.searchsorted(s, p, side='right')) - 1
        if 0 <= j < len(s) and s[j] <= p < e[j]:
            out[i] = v[j]
        else:
            n_gap += 1
    return out, {'sequence_absent': n_nochrom, 'no_bin': n_gap}


def load_genome(path, wanted):
    """FASTA sequences for the wanted (normalised) sequence names only
    (.gz accepted)."""
    seqs, name, buf = {}, None, []
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None and name in wanted:
                    seqs[name] = ''.join(buf).upper()
                name = norm_chrom(line[1:].split()[0])
                buf = []
            elif name in wanted:
                buf.append(line.strip())
    if name is not None and name in wanted:
        seqs[name] = ''.join(buf).upper()
    missing = sorted(set(wanted) - set(seqs))
    logging.info('  genome: %d of %d needed sequences found%s', len(seqs),
                 len(wanted), f' (missing: {missing[:5]})' if missing else '')
    return seqs


def assign_gc(chrom, mid1, genome, window):
    """G+C fraction in a `window`-bp window centred on each 1-based midpoint;
    N and other non-ACGT characters are left out of the denominator."""
    half = window // 2
    out = np.full(len(chrom), np.nan)
    for i, (c, m) in enumerate(zip(chrom, mid1)):
        seq = genome.get(c)
        if seq is None or not np.isfinite(m):
            continue
        c0 = int(m) - 1
        sub = seq[max(0, c0 - half):min(len(seq), c0 + half)]
        gc = sub.count('G') + sub.count('C')
        at = sub.count('A') + sub.count('T')
        if gc + at:
            out[i] = gc / (gc + at)
    return out


# ---------------------------------------------------------------------------
# Germline -> somatic coordinates
# ---------------------------------------------------------------------------

def load_coord_map(path):
    """The germline-to-somatic map, in one of three layouts.

      4 columns  germ_chrom germ_start germ_end soma_chrom
                 (the May layout: the somatic chromosome starts at germ_start,
                 same orientation)
      5 columns  ... soma_chrom strand          (+ or -)
      6+ columns ... soma_chrom soma_start soma_end [strand]

    germ_start / soma_start are 0-based, BED style. Returns (DataFrame with
    germ_chrom, germ_start, germ_end, soma_chrom, soma_start, strand; layout
    name).
    """
    raw = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)
    if raw.shape[1] < 4:
        sys.exit(f'ERROR: {path} has {raw.shape[1]} columns; need at least 4 '
                 '(germ_chrom germ_start germ_end soma_chrom).')
    df = pd.DataFrame({'germ_chrom': raw[0].map(norm_chrom),
                       'germ_start': pd.to_numeric(raw[1], errors='coerce'),
                       'germ_end': pd.to_numeric(raw[2], errors='coerce'),
                       'soma_chrom': raw[3].map(norm_chrom)})
    layout = '4-column (offset = germ_start, + strand)'
    df['soma_start'] = 0
    df['strand'] = '+'
    if raw.shape[1] >= 6 and pd.to_numeric(raw[4], errors='coerce').notna().all() \
            and pd.to_numeric(raw[5], errors='coerce').notna().all():
        df['soma_start'] = pd.to_numeric(raw[4])
        layout = '6-column (explicit somatic start)'
        if raw.shape[1] >= 7 and raw[6].isin(['+', '-']).all():
            df['strand'] = raw[6]
            layout += ' with strand'
    elif raw.shape[1] >= 5 and raw[4].isin(['+', '-']).all():
        df['strand'] = raw[4]
        layout = '5-column (strand)'
    elif raw.shape[1] > 4:
        logging.warning('  coordinate map: columns beyond the 4th are neither a '
                        'strand nor somatic coordinates; ignored (first row: %s)',
                        raw.iloc[0].tolist())
    bad = df[['germ_start', 'germ_end']].isna().any(axis=1) | \
        (df['germ_end'] <= df['germ_start'])
    if bad.any():
        sys.exit(f'ERROR: {int(bad.sum())} rows of {path} have unusable germline '
                 f'coordinates, e.g. {raw[bad].iloc[0].tolist()}')
    df[['germ_start', 'germ_end', 'soma_start']] = \
        df[['germ_start', 'germ_end', 'soma_start']].astype(np.int64)
    per_soma = df.groupby('soma_chrom').size()
    multi = per_soma[per_soma > 1]
    if len(multi) and not layout.startswith('6'):
        sys.exit(f'ERROR: {len(multi)} somatic chromosomes are built from more '
                 f'than one germline segment ({", ".join(multi.index[:5])}), and '
                 'this map has no somatic start column to place them. Supply a '
                 '6-column map.')
    logging.info('  coordinate map: %d segments, %d germline -> %d somatic '
                 'sequences; layout %s', len(df), df['germ_chrom'].nunique(),
                 df['soma_chrom'].nunique(), layout)
    return df, layout


def to_somatic(chrom, mid1, cmap):
    """Somatic chromosome and 0-based position for 1-based germline midpoints.

    Positions outside every mapped segment get (None, NaN). For a '-' segment
    the offset is counted from the segment's end.
    """
    seg = {c: g.sort_values('germ_start') for c, g in cmap.groupby('germ_chrom')}
    out_c = np.full(len(chrom), None, dtype=object)
    out_p = np.full(len(chrom), np.nan)
    for i, (c, m) in enumerate(zip(chrom, mid1)):
        g = seg.get(c)
        if g is None or not np.isfinite(m):
            continue
        p0 = int(m) - 1
        s = g['germ_start'].to_numpy()
        j = int(np.searchsorted(s, p0, side='right')) - 1
        if j < 0 or p0 >= g['germ_end'].iloc[j]:
            continue
        r = g.iloc[j]
        off = (p0 - r['germ_start'] if r['strand'] == '+'
               else r['germ_end'] - 1 - p0)
        out_c[i] = r['soma_chrom']
        out_p[i] = r['soma_start'] + off
    return out_c, out_p


_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def orient_from_sequence(cmap, germ, soma, k=60, n_probes=9, pad=200_000,
                         min_agree=3):
    """Strand and offset of every somatic chromosome, read from the sequences.

    A 4-column map says only which germline segment became which somatic
    chromosome; it assumes the somatic sequence starts at germ_start and runs
    the same way. Here `n_probes` k-mers spread along each somatic chromosome
    are searched (exact, unique matches only) in its germline segment +/- `pad`,
    forwards and reverse-complemented. A strand is accepted when at least
    `min_agree` probes give the same constant, and the map row is rewritten
    from it:
        +  soma = germ - c          (soma_start = germ_start - c)
        -  soma = c - germ          (soma_start = c - (germ_end - 1))
    all 0-based. Returns (corrected map, per-chromosome table).
    """
    out = cmap.copy()
    rows = []
    for idx, r in cmap.iterrows():
        g, sq = germ.get(r['germ_chrom']), soma.get(r['soma_chrom'])
        row = {'soma_chrom': r['soma_chrom'], 'germ_chrom': r['germ_chrom'],
               'segment_bp': int(r['germ_end'] - r['germ_start']),
               'soma_bp': len(sq) if sq is not None else np.nan,
               'probes': 0, 'fwd_hits': 0, 'rev_hits': 0, 'strand': 'unresolved',
               'agrees_with_map': None, 'offset_shift_bp': np.nan}
        if g is None or sq is None or len(sq) < k * 2:
            row['strand'] = 'no sequence'
            rows.append(row)
            continue
        lo = max(0, int(r['germ_start']) - pad)
        hi = min(len(g), int(r['germ_end']) + pad)
        region = g[lo:hi]
        rc = revcomp(region)
        L = len(region)
        fwd, rev = [], []
        for f in np.linspace(0.05, 0.95, n_probes):
            i = int(f * (len(sq) - k))
            w = sq[i:i + k]
            if 'N' in w or len(w) < k:
                continue
            row['probes'] += 1
            j = region.find(w)
            if j >= 0 and region.find(w, j + 1) < 0:
                fwd.append(lo + j - i)                    # germ = soma + c
            j = rc.find(w)
            if j >= 0 and rc.find(w, j + 1) < 0:
                rev.append(lo + L - 1 - j + i)            # germ = c - soma
        row['fwd_hits'], row['rev_hits'] = len(fwd), len(rev)
        best = None
        for strand, consts in (('+', fwd), ('-', rev)):
            if consts:
                vals, cnt = np.unique(consts, return_counts=True)
                c, n = int(vals[cnt.argmax()]), int(cnt.max())
                if n >= min_agree and (best is None or n > best[2]):
                    best = (strand, c, n)
        if best is not None:
            strand, c, _ = best
            row['strand'] = strand
            if strand == '+':
                new_start = int(r['germ_start']) - c
                row['agrees_with_map'] = (r['strand'] == '+' and
                                          new_start == int(r['soma_start']))
                row['offset_shift_bp'] = new_start - int(r['soma_start'])
            else:
                new_start = c - (int(r['germ_end']) - 1)
                row['agrees_with_map'] = r['strand'] == '-' and \
                    new_start == int(r['soma_start'])
            out.loc[idx, 'strand'] = strand
            out.loc[idx, 'soma_start'] = new_start
        rows.append(row)
    tab = pd.DataFrame(rows)
    logging.info('  sequence check: %d +, %d -, %d unresolved, %d without '
                 'sequence; %d disagree with the map',
                 int((tab['strand'] == '+').sum()),
                 int((tab['strand'] == '-').sum()),
                 int((tab['strand'] == 'unresolved').sum()),
                 int((tab['strand'] == 'no sequence').sum()),
                 int((tab['agrees_with_map'] == False).sum()))  # noqa: E712
    return out, tab


def coord_map_check(cmap, track):
    """Per somatic chromosome: mapped segment length against the eigenvector
    track's extent. Equal within a few kb (telomere addition) is expected."""
    extent = track.groupby('chrom')['end'].max()
    rows = []
    for sc, g in cmap.groupby('soma_chrom'):
        span = int((g['soma_start'] + g['germ_end'] - g['germ_start']).max())
        ext = extent.get(sc, np.nan)
        rows.append({'soma_chrom': sc, 'n_segments': len(g),
                     'germ_chrom': ','.join(sorted(g['germ_chrom'].unique())),
                     'mapped_extent_bp': span,
                     'track_extent_bp': ext,
                     'difference_bp': (ext - span) if np.isfinite(ext) else np.nan})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# EV1 orientation
# ---------------------------------------------------------------------------

def orientation(df, flip, mirror_margin=0.15):
    """Per somatic chromosome Spearman(EV1, GC); with `flip`, negate EV1 on
    every chromosome where it is negative. Returns (oriented EV1, table).

    The table also carries |Spearman| for EV1 read at the MIRRORED position
    (`ev1_mirror`, the same offset counted from the other end of the somatic
    chromosome). A sign flip cannot repair a map that runs the wrong way along
    a chromosome; if reading it backwards correlates with GC clearly better
    (by `mirror_margin`) than reading it forwards, the chromosome is flagged.
    """
    ev = df['ev1_raw'].to_numpy(float).copy()
    rows = []
    for c, g in df.dropna(subset=['ev1_raw', 'gc']).groupby('soma_chrom'):
        row = {'soma_chrom': c, 'n_genes': len(g), 'rho_before': np.nan,
               'flipped': False, 'rho_after': np.nan,
               'abs_rho_mirrored': np.nan, 'possible_mirror': False}
        if len(g) < 5 or g['ev1_raw'].std(ddof=0) == 0:
            rows.append(row)
            continue
        r = float(spearmanr(g['ev1_raw'], g['gc'])[0])
        do = flip and r < 0
        if do:
            ev[df.index.get_indexer(g.index)] *= -1
        gm = g.dropna(subset=['ev1_mirror'])
        if len(gm) >= 5 and gm['ev1_mirror'].std(ddof=0) > 0:
            rm = abs(float(spearmanr(gm['ev1_mirror'], gm['gc'])[0]))
            row.update(abs_rho_mirrored=rm,
                       possible_mirror=bool(rm > abs(r) + mirror_margin))
        row.update(rho_before=r, flipped=do, rho_after=-r if do else r)
        rows.append(row)
    return ev, pd.DataFrame(rows)


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
    """P in the manuscript's notation: 0.49, 2.1 × 10⁻⁵."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p >= 0.001:
        return f'{p:.2g}'
    e = int(np.floor(np.log10(p)))
    m = p / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    sup = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    return f'{m:.1f} × 10{str(e).translate(sup)}'


def signed(x, nd):
    """x to `nd` decimals with a true minus sign (U+2212), as typeset."""
    return f'{x:.{nd}f}'.replace('-', '−')


def ab_tests(df):
    """A vs B per metric (Mann-Whitney; BH q over the metrics) and Spearman
    of each metric against oriented EV1."""
    rows = []
    for m in METRICS:
        a = df.loc[df['compartment'] == 'A', m].dropna().to_numpy(float)
        b = df.loc[df['compartment'] == 'B', m].dropna().to_numpy(float)
        p = (mannwhitneyu(a, b, alternative='two-sided')[1]
             if len(a) >= 5 and len(b) >= 5 else np.nan)
        s = df[['ev1', m]].dropna()
        rho, rp = (spearmanr(s['ev1'], s[m]) if len(s) >= 10 else (np.nan, np.nan))
        rows.append({'metric': m, 'n_A': len(a), 'n_B': len(b),
                     'median_A': np.median(a) if len(a) else np.nan,
                     'median_B': np.median(b) if len(b) else np.nan,
                     'p': p, 'n_spearman': len(s), 'rho_ev1': float(rho),
                     'p_rho': float(rp)})
    out = pd.DataFrame(rows)
    out['q'] = benjamini_hochberg(out['p'].to_numpy(float))
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


def draw_metric(ax, df, metric, test, rng, jitter_max, top_pct, box_alpha=0.40):
    """One box panel: `metric` in A vs B compartment genes, all retained.

    Box = quartiles, whiskers = 5th-95th percentiles, jitter on top; both boxes
    retained blue. The y-axis runs from 2% below 0, so a gene at 0 is not cut
    by the axis, to a limit: the `top_pct` percentile of the pooled values,
    raised to the highest whisker end if that is higher, plus 3% so a triangle
    sits clear of a whisker cap at the limit. Genes above the limit are drawn
    at it as open triangles (the 6D convention), so a clamped value is never
    read as a real one. Counts are of all genes, then of the triangles in the
    jittered sample. Drawing only.
    """
    data = [df.loc[df['compartment'] == c, metric].dropna().to_numpy(float)
            for c in COMPARTMENTS]
    kw = dict(positions=[0, 1], widths=0.55, patch_artist=True,
              showfliers=False, whis=(5, 95), zorder=2)
    bp = ax.boxplot(data, **kw)
    for patch in bp['boxes']:
        patch.set_facecolor(to_rgba(REGION['retained'], box_alpha))
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
    pct_value = float(np.percentile(np.concatenate(data), top_pct))
    hi = max([pct_value] + [float(np.percentile(d, 95)) for d in data if len(d)])
    hi *= 1.03
    n_above, n_tri, n_skip = {}, {}, {}
    for i, (c, d) in enumerate(zip(COMPARTMENTS, data)):
        show = d if len(d) <= jitter_max else rng.choice(d, jitter_max,
                                                         replace=False)
        n_skip[c] = len(d) - len(show)
        n_above[c] = int((d > hi).sum())
        x = i + rng.normal(0, 0.06, len(show))
        above = show > hi
        n_tri[c] = int(above.sum())
        ax.scatter(x[~above], show[~above], s=1.4, linewidths=0,
                   color=REGION['retained'], alpha=0.35, zorder=3,
                   rasterized=True)
        if n_tri[c]:
            ax.scatter(x[above], np.full(n_tri[c], hi), marker='^', s=4.0,
                       facecolors='none', edgecolors=REGION['retained'],
                       linewidths=0.4, zorder=4, clip_on=False)
    top = hi * 1.06
    ax.plot([0, 0, 1, 1], [top, top * 1.03, top * 1.03, top], color='black',
            lw=LINEWIDTH['thin'], zorder=6)
    ax.text(0.5, top * 1.045, stars(test['q']), ha='center', va='bottom',
            fontsize=FONT['annotation'])
    draw_p(ax, 0.5, top * 1.15, test['p'], ha='center', va='bottom',
           fontsize=FONT['annotation'])
    ax.set_ylim(-0.02 * hi, top * 1.28)
    ax.set_xlim(-0.55, 1.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'{c} compartment\n(n = {len(d):,})'
                        for c, d in zip(COMPARTMENTS, data)])
    ax.tick_params(axis='x', length=0)
    ax.set_ylabel(METRIC_LABEL[metric])
    return {'y_pct_value': pct_value, 'y_limit': hi, 'top_pct': top_pct,
            'above_limit': n_above, 'triangles_drawn': n_tri,
            'jitter_not_drawn': n_skip, 'jitter_max': jitter_max}


def plot_s9a(df, tests, outdir, stem, rng, args):
    """S13A: dN/dS, dN, dS by compartment, one 185 mm strip."""
    fig = blank_figure('double', args.height_mm)
    axes = [fig.add_axes([0.07 + 0.33 * k, 0.2, 0.25, 0.6]) for k in range(3)]
    info = {}
    for ax, m in zip(axes, METRICS):
        t = tests.set_index('metric').loc[m]
        info[m] = draw_metric(ax, df, m, t, rng, args.jitter_max, args.top_pct)
    place_row(fig, axes, shares=[1, 1, 1], gaps_mm=[8.0, 8.0])
    paths, problems = finish(fig, outdir, stem)
    return paths, problems, info


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def write_summary(path, ctx):
    a = ctx['args']
    tests, orient, check = ctx['tests'], ctx['orient'], ctx['check']
    L = ['dN, dS and dN/dS of retained genes by somatic compartment '
         '(supplementary fig. S13A)', '=' * 72, '',
         f'Rate table    {a.annotated}',
         f'Eigenvector   {a.eigenvector}' +
         (f'  ({a.ev_label})' if a.ev_label else ''),
         f'Coord map     {a.coord_map}  [{ctx["layout"]}]',
         f'Genome        {a.genome}  (GC in {a.gc_window:,} bp windows)',
         f'Orientation   EV1 sign per somatic chromosome so that '
         f'Spearman(EV1, GC) > 0: {"ON" if a.orient_ev_by_gc else "OFF"}', '',
         'Genes at each step                               retained  eliminated']
    for k, (r, e) in ctx['steps'].items():
        L.append(f'  {k:<46} {r:>8,}  {e:>10,}')
    if ctx['steps']['with a somatic position'][1]:
        L.append('  WARNING: eliminated genes received somatic positions -- the '
                 'coordinate map covers eliminated sequence. Check it.')
    L += [f"  EV1 lookup misses: somatic sequence absent from the track "
          f"{ctx['lookup']['sequence_absent']:,}, position in no bin "
          f"{ctx['lookup']['no_bin']:,}; EV1 exactly 0 (unassigned) "
          f"{ctx['n_zero']:,}", '']

    L += ['-' * 72, 'A vs B, retained genes (Mann-Whitney; BH q over the three)',
          '-' * 72,
          f"{'metric':<7} {'n A':>6} {'n B':>6} {'median A':>9} {'median B':>9} "
          f"{'P':>9} {'q':>9}   Spearman vs EV1"]
    for t in tests.itertuples():
        L.append(f"{METRIC_LABEL[t.metric]:<7} {t.n_A:>6,} {t.n_B:>6,} "
                 f"{t.median_A:>9.4f} {t.median_B:>9.4f} {p_text(t.p):>9} "
                 f"{p_text(t.q):>9}   rho = {t.rho_ev1:+.4f}, P = "
                 f"{p_text(t.p_rho)} (n = {t.n_spearman:,})")
    if ctx['unoriented'] is not None:
        u = ctx['unoriented']
        L += ['', 'Without orientation (raw FAN-C sign; the May analysis):']
        for t in u.itertuples():
            L.append(f"  {METRIC_LABEL[t.metric]:<6} median A {t.median_A:.4f}, "
                     f"B {t.median_B:.4f}, P = {p_text(t.p)}; rho = "
                     f"{t.rho_ev1:+.4f}, P = {p_text(t.p_rho)}")
    L.append('')

    L += ['-' * 72, 'EV1 orientation per somatic chromosome (ev1_orientation.tsv)',
          '-' * 72]
    fl = orient[orient['flipped']]
    L.append(f"  {len(orient)} chromosomes; flipped: {len(fl)}"
             + (f" ({', '.join(fl['soma_chrom'])})" if len(fl) else ''))
    mir = orient[orient['possible_mirror']]
    if len(mir):
        L.append(f"  READ BACKWARDS, EV1 tracks GC better on {len(mir)}: " +
                 ', '.join(f"{r.soma_chrom} (|rho| {abs(r.rho_before):.2f} "
                           f"forwards, {r.abs_rho_mirrored:.2f} backwards)"
                           for r in mir.itertuples()))
        if ctx['seq_check'] is None:
            L.append('  The coordinate map may run the wrong way along these '
                     'somatic chromosomes. Rerun with --soma-genome to settle it '
                     'from sequence.')
        else:
            L.append('  Heuristic only (GC and EV1 are both smooth along a '
                     'chromosome, so a mirrored read can correlate by chance); '
                     'the sequence check above decides orientation.')
    weak = orient[orient['rho_after'] < a.weak_rho]
    if len(weak):
        L.append(f"  Spearman(EV1, GC) still below {a.weak_rho} after orientation "
                 f"on {len(weak)}: " + ', '.join(
                     f"{r.soma_chrom} ({r.rho_after:+.2f}, n = {r.n_genes})"
                     for r in weak.itertuples()))
        L.append('  A weak correlation is expected on small or compartment-poor '
                 'chromosomes; on a large one it can mean the coordinate map '
                 'mirrors that chromosome (no strand column). Inspect them.')
    else:
        L.append(f'  Every chromosome has Spearman(EV1, GC) >= {a.weak_rho} '
                 'after orientation.')
    L.append(f"  median rho after orientation {orient['rho_after'].median():+.3f}")
    L.append('')

    sc = ctx['seq_check']
    L += ['-' * 72, 'Strand and offset from sequence (coord_map_sequence_check.tsv)',
          '-' * 72]
    if sc is None:
        L.append('  Not run (no --soma-genome). The map is taken as given: '
                 'every somatic chromosome assumed to start at germ_start and '
                 'run the same way.')
    else:
        n_plus = int((sc['strand'] == '+').sum())
        n_minus = int((sc['strand'] == '-').sum())
        bad = sc[sc['agrees_with_map'] == False]  # noqa: E712
        unres = sc[~sc['strand'].isin(['+', '-'])]
        L.append(f"  {len(sc)} somatic chromosomes: {n_plus} run the same way as "
                 f"the germline, {n_minus} reverse-complemented; "
                 f"{len(unres)} unresolved.")
        if len(bad):
            L.append(f"  {len(bad)} differ from the map's assumption: " + ', '.join(
                f"{r.soma_chrom} ({r.strand}"
                + (f', offset {r.offset_shift_bp:+,.0f} bp' if r.strand == '+'
                   else '') + ')' for r in bad.itertuples()))
            L.append('  ' + ('Used as given (--trust-map).' if a.trust_map else
                             'Positions were converted with the corrected strand '
                             'and offset.'))
        else:
            L.append("  Every resolved chromosome agrees with the map.")
        if len(unres):
            L.append('  Unresolved (map used as given): ' + ', '.join(
                f"{r.soma_chrom} ({r.strand}; {r.fwd_hits} fwd / {r.rev_hits} "
                f"rev unique hits of {r.probes})" for r in unres.itertuples()))
    L.append('')

    L += ['-' * 72, 'Coordinate-map check (coord_map_check.tsv)', '-' * 72]
    d = check['difference_bp'].abs()
    L.append(f"  {len(check)} somatic chromosomes; track extent minus mapped "
             f"extent: median {np.nanmedian(check['difference_bp']):,.0f} bp, "
             f"largest |difference| {np.nanmax(d) if d.notna().any() else np.nan:,.0f} bp")
    big = check[d > a.extent_tol]
    if len(big):
        L.append(f"  Differ by more than {a.extent_tol:,} bp: " +
                 ', '.join(f"{r.soma_chrom} ({r.difference_bp:+,.0f})"
                           for r in big.head(8).itertuples()))
    nolen = check['track_extent_bp'].isna().sum()
    if nolen:
        L.append(f'  {nolen} mapped somatic chromosomes are absent from the '
                 'eigenvector track.')
    L.append('')

    fi, probs = ctx['figure'], ctx['problems']
    L += ['-' * 72, 'Drawn (figS9A_compartments)', '-' * 72]
    for m in METRICS:
        f = fi[m]
        L.append(f"  {METRIC_LABEL[m]}: y-axis limit {f['y_limit']:.4f} (the "
                 f"{f['top_pct']:g}th percentile, {f['y_pct_value']:.4f}, or the "
                 f"highest whisker end, plus 3%); genes above it A "
                 f"{f['above_limit']['A']:,}, B {f['above_limit']['B']:,}, drawn "
                 f"as open triangles at the limit A "
                 f"{f['triangles_drawn']['A']:,}, B {f['triangles_drawn']['B']:,}; "
                 f"jitter not drawn A {f['jitter_not_drawn']['A']:,}, B "
                 f"{f['jitter_not_drawn']['B']:,} (at most {f['jitter_max']:,} "
                 'per box)')
    L.append('  Box = quartiles, whiskers = 5th-95th percentiles; both boxes '
             'retained blue; the axis starts 2% below 0.')
    L.append('  Asterisks from BH q; exact P printed.')
    L.append('  AUDIT: ' + '; '.join(probs) if probs else '  audit: clean')
    L.append('')

    t = tests.set_index('metric').loc['dN_dS']
    L += ['-' * 72, 'What the text can claim', '-' * 72,
          f"  Among retained genes, for which a somatic compartment identity is "
          f"defined, dN/dS is {'uncorrelated' if t['p_rho'] >= 0.05 else 'CORRELATED'}"
          f" with the compartment eigenvector (Spearman ρ = "
          f"{signed(t['rho_ev1'], 3)}; P = {p_manuscript(t['p_rho'])}) and "
          f"{'indistinguishable' if t['p'] >= 0.05 else 'DIFFERS'} between A- and "
          f"B-compartment genes (median {t['median_A']:.3f} versus "
          f"{t['median_B']:.3f}; P = {p_manuscript(t['p'])}; supplementary fig. "
          'S13A, Supplementary Material online).']
    if t['p_rho'] < 0.05 or t['p'] < 0.05:
        L.append('  The capitalised word means the null statement no longer holds: '
                 'rewrite the sentence.')
    L.append('')
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='dN, dS and dN/dS of retained genes by somatic A/B '
                    'compartment (supplementary fig. S13A).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', '--kaks', dest='annotated', required=True,
                   help='kaks_annotated.tsv from annotate_regions.py')
    p.add_argument('--eigenvector', required=True,
                   help='Somatic (post-elimination) EV1 track, 4 columns')
    p.add_argument('--coord-map', required=True,
                   help='Germline-to-somatic coordinate map (see load_coord_map)')
    p.add_argument('--genome', required=True,
                   help='Germline genome FASTA, for the GC used to orient EV1')
    p.add_argument('--soma-genome', default=None,
                   help='Somatic (post-elimination) assembly FASTA. When given, '
                        'the strand and offset of every somatic chromosome are '
                        'read from sequence and the map corrected where it '
                        'disagrees (.gz accepted)')
    p.add_argument('--trust-map', action='store_true',
                   help='With --soma-genome, report the sequence check but '
                        'convert with the map as given')
    p.add_argument('--gc-window', type=int, default=20000)
    p.add_argument('--orient-ev-by-gc', action=argparse.BooleanOptionalAction,
                   default=True,
                   help='Flip EV1 per somatic chromosome to correlate positively '
                        'with GC (default on; --no-orient-ev-by-gc reproduces '
                        'the May labels)')
    p.add_argument('--ev-label', default=None,
                   help='Eigenvector stage for the summary, e.g. "somatic (10 day)"')
    p.add_argument('--gene-col', default='gene1')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--chrom-col', default='chrom1')
    p.add_argument('--pos-col', default='midpoint1',
                   help='1-based germline midpoint (default midpoint1; from '
                        'start1/end1 if absent)')
    p.add_argument('--top-pct', type=float, default=99.5,
                   help='y-axis top as a pooled percentile; drawing only')
    p.add_argument('--jitter-max', type=int, default=1500)
    p.add_argument('--weak-rho', type=float, default=0.1,
                   help='List chromosomes whose Spearman(EV1, GC) stays below '
                        'this after orientation (default 0.1)')
    p.add_argument('--extent-tol', type=int, default=50000,
                   help='Flag somatic chromosomes whose mapped and track extents '
                        'differ by more than this (bp; default 50,000)')
    p.add_argument('--height-mm', type=float, default=58.0)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.outdir = os.path.join('results', f'kaks_compartments_{stamp}')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'kaks_compartments')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('kaks_compartments.py -- supplementary fig. S13A')
    for k in ('annotated', 'eigenvector', 'coord_map', 'genome',
              'soma_genome'):
        f = getattr(args, k)
        if f is None:
            continue
        logging.info('  %-11s %s', k, os.path.abspath(f))
        if not os.path.isfile(f):
            sys.exit(f'ERROR: not found: {f}')
    logging.info('  outdir      %s', os.path.abspath(args.outdir))
    apply_style()
    rng = np.random.default_rng(args.seed)

    # --- table ----------------------------------------------------------------
    tab = read_rate_table(args.annotated, require=METRICS)
    for col in (args.gene_col, args.region_col, args.chrom_col):
        if col not in tab.columns:
            sys.exit(f'ERROR: column {col!r} not in {args.annotated}. '
                     f'Columns: {list(tab.columns)}')
    if args.pos_col in tab.columns:
        mid = pd.to_numeric(tab[args.pos_col], errors='coerce')
    elif {'start1', 'end1'} <= set(tab.columns):
        mid = (pd.to_numeric(tab['start1'], errors='coerce')
               + pd.to_numeric(tab['end1'], errors='coerce')) // 2
        logging.info('  %s absent; midpoint computed from start1/end1',
                     args.pos_col)
    else:
        sys.exit(f'ERROR: neither {args.pos_col!r} nor start1/end1 in '
                 f'{args.annotated}')
    df = pd.DataFrame({
        'gene1': tab[args.gene_col].astype(str),
        'region1': tab[args.region_col],
        'chrom1': tab[args.chrom_col].map(lambda c: norm_chrom(c)
                                          if pd.notna(c) else np.nan),
        'midpoint1': mid.to_numpy(float)})
    for m in METRICS:
        df[m] = pd.to_numeric(tab[m], errors='coerce')
    steps = {}

    def count(label, mask):
        steps[label] = (int((mask & (df['region1'] == 'retained')).sum()),
                        int((mask & (df['region1'] == 'eliminated')).sum()))
        logging.info('  %-44s retained %s, eliminated %s', label,
                     f'{steps[label][0]:,}', f'{steps[label][1]:,}')

    logging.info('Table: %s rows', f'{len(df):,}')
    count('retained or eliminated', df['region1'].isin(['retained',
                                                        'eliminated']))
    count('with a germline position', df['chrom1'].notna()
          & df['midpoint1'].notna())

    # --- somatic coordinates and EV1 ------------------------------------------
    logging.info('Coordinate map')
    cmap, layout = load_coord_map(args.coord_map)
    logging.info('Germline assembly')
    genome = load_genome(args.genome, set(df['chrom1'].dropna())
                         | set(cmap['germ_chrom']))
    seq_check = None
    if args.soma_genome:
        logging.info('Somatic assembly: strand and offset of each somatic '
                     'chromosome from sequence')
        soma = load_genome(args.soma_genome, set(cmap['soma_chrom']))
        cmap_seq, seq_check = orient_from_sequence(cmap, genome, soma)
        del soma
        seq_check.to_csv(os.path.join(args.outdir, 'coord_map_sequence_check.tsv'),
                         sep='\t', index=False)
        if args.trust_map:
            logging.info('  --trust-map: map used as given')
        else:
            cmap = cmap_seq
            layout += '; strand and offset corrected from sequence'
    df['soma_chrom'], df['soma_pos0'] = to_somatic(
        df['chrom1'].to_numpy(), df['midpoint1'].to_numpy(), cmap)
    count('with a somatic position', df['soma_chrom'].notna())

    logging.info('Eigenvector')
    track, _ = load_track(args.eigenvector, 'ev1')
    df['ev1_raw'], lookup = lookup_track(df['soma_chrom'].to_numpy(),
                                         df['soma_pos0'].to_numpy(), track, 'ev1')
    count('with somatic EV1', df['ev1_raw'].notna())
    logging.info('  EV1 lookup misses: sequence absent %s, no bin %s',
                 f"{lookup['sequence_absent']:,}", f"{lookup['no_bin']:,}")
    # EV1 at the mirrored offset, for the orientation check only.
    extent = track.groupby('chrom')['end'].max().to_dict()
    mirror_pos = np.array([extent[c] - 1 - p if isinstance(c, str) and c in extent
                           and np.isfinite(p) else np.nan
                           for c, p in zip(df['soma_chrom'], df['soma_pos0'])])
    df['ev1_mirror'], _ = lookup_track(df['soma_chrom'].to_numpy(), mirror_pos,
                                       track, 'ev1')
    check = coord_map_check(cmap, track)
    check.to_csv(os.path.join(args.outdir, 'coord_map_check.tsv'), sep='\t',
                 index=False)

    logging.info('GC (%s bp windows, germline assembly)', f'{args.gc_window:,}')
    need = df['ev1_raw'].notna()
    df['gc'] = np.nan
    df.loc[need, 'gc'] = assign_gc(df.loc[need, 'chrom1'].to_numpy(),
                                   df.loc[need, 'midpoint1'].to_numpy(),
                                   genome, args.gc_window)
    del genome

    df['ev1'], orient = orientation(df, args.orient_ev_by_gc)
    orient.to_csv(os.path.join(args.outdir, 'ev1_orientation.tsv'), sep='\t',
                  index=False, float_format='%.4f')
    if args.orient_ev_by_gc:
        fl = orient[orient['flipped']]
        logging.info('  EV1 sign flipped on %d of %d somatic chromosomes%s',
                     len(fl), len(orient),
                     (': ' + ', '.join(fl['soma_chrom'])) if len(fl) else '')
    df['compartment'] = np.where(df['ev1'] > 0, 'A',
                                 np.where(df['ev1'] < 0, 'B', None))
    n_zero = int((df['ev1'] == 0).sum())
    df.to_csv(os.path.join(args.outdir, 'compartment_per_gene.tsv'), sep='\t',
              index=False, float_format='%.6g')

    ret = df[(df['region1'] == 'retained') & df['compartment'].notna()].copy()
    count('retained with a compartment (analysed)',
          (df['region1'] == 'retained') & df['compartment'].notna())
    if ret['compartment'].nunique() < 2 or len(ret) < 20:
        sys.exit('ERROR: fewer than 20 retained genes with a compartment, or only '
                 'one compartment. Check chromosome names across the table, the '
                 'map and the eigenvector (see the counts above).')

    tests = ab_tests(ret)
    unoriented = None
    if args.orient_ev_by_gc and orient['flipped'].any():
        raw = ret.copy()
        raw['ev1'] = raw['ev1_raw']
        raw['compartment'] = np.where(raw['ev1'] > 0, 'A',
                                      np.where(raw['ev1'] < 0, 'B', None))
        unoriented = ab_tests(raw[raw['compartment'].notna()])
    out = tests.assign(orientation='by GC' if args.orient_ev_by_gc else 'raw')
    if unoriented is not None:
        out = pd.concat([out, unoriented.assign(orientation='raw')])
    out.to_csv(os.path.join(args.outdir, 'compartment_tests.tsv'), sep='\t',
               index=False, float_format='%.6g')
    for t in tests.itertuples():
        logging.info('  %-6s A %.4f (n %s), B %.4f (n %s), P %s; rho vs EV1 %+.4f '
                     '(P %s)', METRIC_LABEL[t.metric], t.median_A, f'{t.n_A:,}',
                     t.median_B, f'{t.n_B:,}', p_text(t.p), t.rho_ev1,
                     p_text(t.p_rho))

    logging.info('S13A: compartment strip')
    paths, problems, fig_info = plot_s9a(ret, tests, args.outdir,
                                         'figS9A_compartments', rng, args)

    write_summary(os.path.join(args.outdir, 'compartment_summary.txt'), {
        'args': args, 'layout': layout, 'steps': steps, 'lookup': lookup,
        'n_zero': n_zero, 'tests': tests, 'unoriented': unoriented,
        'seq_check': seq_check,
        'orient': orient, 'check': check, 'figure': fig_info,
        'problems': problems})

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
                   'inputs': {k: fingerprint(getattr(args, k)) for k in
                              ('annotated', 'eigenvector', 'coord_map')},
                   'genome': os.path.abspath(args.genome),
                   'soma_genome': (os.path.abspath(args.soma_genome)
                                   if args.soma_genome else None),
                   'sequence_check': (seq_check.to_dict('records')
                                      if seq_check is not None else None),
                   'coord_map_layout': layout,
                   'steps': {k: list(v) for k, v in steps.items()},
                   'flipped': orient.loc[orient['flipped'], 'soma_chrom'].tolist(),
                   'figure': [os.path.basename(p) for p in paths],
                   'audit': problems},
                  fh, indent=2, default=str)
    logging.info('Summary: %s', os.path.join(args.outdir,
                                             'compartment_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
