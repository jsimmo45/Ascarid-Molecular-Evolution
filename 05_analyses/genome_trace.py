#!/usr/bin/env python3
"""
genome_trace.py - dN/dS, dN or dS along whole chromosomes, drawn to scale.

Figure 3 (dN/dS, one run per species) and supplementary fig. S6 (the dN and dS
runs, and the other two species; S4 before the 2026-10-12 renumbering). Renamed
from `plot_kaks_genome_trace.py`.

Writes, per run:

    <metric>_genome_trace_<species>.{svg,pdf,png}
    <metric>_genome_trace_windows.tsv   every sliding-window point, per chromosome;
                                        `isolated` = 1 for a window drawn as a dot
    <metric>_genome_trace_summary.txt   panel layout, per-chromosome counts,
                                        y-scale, truncation and label checks
    <metric>_provenance.json           argv, input hashes, row counts
        (all carry the metric prefix: one outdir holds a whole
         Ka_Ks/Ka/Ks loop without any file overwriting another)


The size problem, which was the real defect
-------------------------------------------
The earlier version built a 20 x N inch canvas with font sizes in points and
saved with ``bbox_inches='tight'``. Placed at the journal's 185 mm double-column
width that is a 2.7x reduction, so its 13 pt chromosome labels print at under
5 pt - below the 6 pt floor - and ``bbox_inches='tight'`` meant the final width
was whatever the content happened to be rather than a column width at all.

This version builds at exactly `--width` (185 mm by default) with heights in mm,
so a 7 pt label prints at 7 pt. Panels are still proportional to sequence length
and rows still hold one Mb-per-inch scale; there is simply less of it.

That constrains what fits. At 185 mm a row holding 75 Mb gives about 2.5 mm per
Mb, which is legible for a 500 kb sliding window and not for individual genes -
so the per-gene scatter is drawn small and rasterized, and the window line
carries the panel. If you need the per-gene detail, that is what the
single-chromosome track in `cbr_meta_profile.py` is for.


Sliding windows
---------------
Computed once per chromosome over all its genes and then clipped to each panel's
x-limits. Computing per slice instead would put an edge effect at every wrap
point, which would look like structure.

`--min-genes` defaults to 10, matching the CBR meta-profile's per-bin floor.
A 500 kb window holds roughly 30 genes in
*Ascaris*, so the floor only bites at chromosome ends and inside eliminated
tracts - which is exactly where a "median" of one gene would otherwise be drawn
as if it were comparable to the rest of the line.

Dropped windows leave gaps, and the trace is broken there rather than joined
across (`break_gaps`). A window whose neighbours on BOTH sides were dropped is
then a lone point between two breaks, and a line needs two points, so it used to
be drawn as nothing at all - 13 of 2,114 Baylisascaris windows in S4, several
of them among the highest. Such windows are now drawn as a dot of the trace
colour, flagged `isolated` in the windows table, and counted in the summary.

Markers are restricted to the panel that owns them. The line is the whole
chromosome clipped to each panel's x-limits, but the y-cap carets are drawn
with clip_on=False, so on a wrapped chromosome every slice used to draw (and
count) the carets of every other slice.


Changes 2026-10-12 (figure review of supplementary fig. S6)
-----------------------------------------------------------
- The y label is plain text, dN/dS, dN or dS, upright: METRIC_LABELS was
  mathtext (d$N$/d$S$), which printed italic N and S in the math font.
- The position axes have a unit: one "Position (Mb)" title under the grid,
  since every tick label is a position in Mb. The bottom margin grows by
  2.5 mm (11 to 13.5) to hold it.
- Every sequence the filters remove is listed in the summary with its length,
  the genes on it and the filter that removed it (no plotted gene,
  --include-chroms / --exclude-chroms, --min-chrom-mb, --max-chroms), so a
  legend can say why a scaffold is missing; chromosomes with too few genes for
  a trace are listed too.
- The legend n, "Genes plotted" and the default y-limit now count only the
  genes on the sequences drawn. They used to include genes on sequences the
  filters removed (unplaced fragments, scaffolds below --min-chrom-mb), so the
  legend counted genes that are not in the figure. The windows table and the
  per-chromosome medians are unchanged.

Species notes
-------------
  Ascaris        24 germline chromosomes; the grid works as intended.
  Parascaris     the germline is one 240 Mb chromosome. `--wrap-mb` cuts it
                 across rows; the continuous version is unreadable at column
                 width and is no longer produced. The assembly names the
                 chromosome chrX; `--label-sub '^chrX$' 'Germline chromosome'`
                 labels the slices 'Germline chromosome 0-60' and so on.
  Baylisascaris  chromosome-scale but with no eliminated-region calls. Run with
                 `--region-col none` so no gene is coloured by another species'
                 status. Its sequence names are 15-character GWH accessions
                 (GWHBAVG00000005), wider than many of the panels they sit
                 under; `--label-sub '^GWHBAVG0*' 'Scaffold '` labels them
                 'Scaffold 5' instead. The summary lists the mapping and
                 reports any adjacent panel labels that still collide.
  Toxocara       scaffold-level; a genome-wide trace is not meaningful.

Usage
-----
    python genome_trace.py \\
      --kaks-table results/annotate_regions_20260817/kaks_annotated.tsv \\
      --gff As_miniprot.gff \\
      --elim-bed AG.v50.eliminated_strict.bed \\
      --chrom-sizes AG.v50.chrom.sizes \\
      --species Ascaris --metric Ka_Ks \\
      --row-mb 75 --break-before '^chrX' \\
      --outdir results/genome_trace_20260817_ascaris

Dependencies: numpy, pandas, matplotlib, paper_colors, region_lib
"""

from __future__ import annotations

import argparse
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
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgba

from paper_colors import (REGION, REGION_SHADE, NEUTRAL, COL, apply_style,
                          audit, save, mm, MAX_HEIGHT_MM, MIN_AXES_MM)
from region_lib import fingerprint, git_commit, setup_logging

_HERE = os.path.dirname(os.path.abspath(__file__))

# Minimum on-page width, in mm, of the rug marker under an eliminated
# tract. 0.6 mm survives the 6% reduction OUP may apply and is still
# clearly a mark rather than a line.
MIN_TRACT_MM = 0.6

# Consecutive window centres more than this many steps apart are not joined by
# the trace. Shared by break_gaps() and isolated_windows(), so a window is drawn
# as a dot exactly when the line has nothing to join it to.
GAP_FACTOR = 1.5

# Diameter, in points, of the dot drawn for an isolated window: about twice the
# 0.9 pt trace, and solid ink, so it reads as a point on the line rather than as
# one of the pale, translucent per-gene points.
ISO_MS = 2.0

# Adjacent panel labels closer than this on the page are reported. Labels are
# centred under their panels, so this is the space left between two of them.
LABEL_GAP_MM = 0.5

# Plain text, dN/dS upright: mathtext printed italic N and S in the math font.
METRIC_LABELS = {'Ka_Ks': 'dN/dS', 'Ka': 'dN', 'Ks': 'dS'}

# The x-axis title under the whole grid. Every tick label is a position in Mb.
POSITION_LABEL = 'Position (Mb)'


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def natural_key(chrom):
    """chr1, chr2, ... chr10 rather than chr1, chr10, chr2; sex chromosomes last."""
    is_sex = 1 if chrom.upper().startswith(('CHRX', 'CHRY', 'X', 'Y')) else 0
    m = re.search(r'\d+', chrom)
    return (is_sex, int(m.group()) if m else 10 ** 9, chrom)


def load_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    sizes[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return sizes


def parse_elim_bed(path):
    """{chrom: [(start, end), ...]}, kept in BED coordinates - these are only
    ever used to draw a shaded span, never to classify a gene."""
    regions = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(('#', 'track')):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    regions.setdefault(parts[0], []).append(
                        (int(parts[1]), int(parts[2])))
                except ValueError:
                    continue
    return regions


def _id_candidates(raw):
    """Every plausible reading of a GFF ID attribute.

    Annotation pipelines wrap identifiers differently and the rate table may key
    on any level: a Cufflinks-style record carries 'ID=GENE^g88.t1,TRANS^g88.t1'
    while the ortholog table uses bare 'g88'. Emitting every reading and keeping
    whichever set intersects the table beats guessing, and it is why this script
    works on both miniprot and Augustus annotations.
    """
    cands = set()
    if not raw:
        return cands
    cands.add(raw)
    for piece in re.split(r'[,;]', raw):
        piece = piece.strip()
        if not piece:
            continue
        cands.add(piece)
        tail = piece.split('^')[-1] if '^' in piece else piece
        cands.add(tail)
        cands.add(re.sub(r'\.t\d+$', '', tail))
        cands.add(re.sub(r'\.\d+$', '', tail))
    cands.discard('')
    return cands


def get_gene_positions(gff_path, feature_types=('mRNA', 'transcript', 'gene',
                                                'match', 'CDS', 'exon')):
    """{feature_type: {gene_id: (chrom, midpoint)}}, plus every feature seen.

    Coordinates are aggregated per identifier (min start, max end) because
    formats such as 'match' emit one line per exon block; taking the first line
    would place a gene at its first exon rather than its span midpoint.
    """
    spans = {ft: {} for ft in feature_types}
    seen = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            ft = parts[2]
            seen[ft] = seen.get(ft, 0) + 1
            if ft not in spans:
                continue
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            attrs, chrom = parts[8], parts[0]
            raw_id = None
            for field in attrs.split(';'):
                field = field.strip()
                if field.startswith('ID='):
                    raw_id = field.split('=', 1)[1]
                    break
            if raw_id is None:
                m = re.search(r'(?:transcript_id|gene_id)[= ]"?([^";]+)', attrs)
                if m:
                    raw_id = m.group(1)
            if not raw_id:
                continue
            d = spans[ft]
            for cand in _id_candidates(raw_id):
                rec = d.get(cand)
                if rec is None:
                    d[cand] = [chrom, start, end]
                elif rec[0] == chrom:
                    rec[1] = min(rec[1], start)
                    rec[2] = max(rec[2], end)
    by_feature = {ft: {k: (v[0], (v[1] + v[2]) // 2) for k, v in d.items()}
                  for ft, d in spans.items()}
    return by_feature, seen


def pick_gene_positions(by_feature, table_ids, forced=None):
    """The GFF feature type whose identifiers best match the rate table."""
    if forced:
        d = by_feature.get(forced, {})
        return d, forced, len(table_ids & set(d))
    best = (None, {}, -1)
    for ft, d in by_feature.items():
        if not d:
            continue
        n = len(table_ids & set(d))
        logging.info("  GFF feature %-11s %8s ids, %8s match", f"'{ft}':",
                     f'{len(d):,}', f'{n:,}')
        if n > best[2]:
            best = (ft, d, n)
    return best[1], best[0], best[2]


def sliding_window(pos, vals, window_bp, step_bp, min_genes=10, stat='median'):
    """Windowed summary over genomic position.

    Returns (centres, values, counts, n_dropped). ``n_dropped`` is the number of
    windows that held at least one gene but fewer than ``min_genes``, so the
    caller can report the truncation instead of leaving it silent -- those
    windows are overwhelmingly at chromosome ends, which is exactly where the
    trace is most eye-catching and least trustworthy.
    """
    pos = np.asarray(pos, dtype=float)
    vals = np.asarray(vals, dtype=float)
    if pos.size == 0:
        return np.array([]), np.array([]), np.array([]), 0
    order = np.argsort(pos)
    pos, vals = pos[order], vals[order]
    func = np.median if stat == 'median' else np.mean
    centres, out, counts, dropped = [], [], [], 0
    for win_start in range(int(pos.min()), int(pos.max()) + 1, int(step_bp)):
        win_end = win_start + window_bp
        m = (pos >= win_start) & (pos < win_end)
        n = int(m.sum())
        if n >= min_genes:
            centres.append((win_start + win_end) / 2.0)
            out.append(float(func(vals[m])))
            counts.append(n)
        elif n:
            dropped += 1
    return (np.array(centres), np.array(out), np.array(counts), dropped)


def break_gaps(cen, val, step_bp):
    """Insert NaN wherever consecutive window centres are more than one step apart.

    Windows below ``--min-genes`` are skipped, not zero-filled, so the surviving
    centres can jump. ``ax.plot`` joins whatever it is given, which draws a
    straight segment across a region that was excluded precisely because it had
    too little data to summarise -- a line where there is no measurement. NaN
    breaks the path instead, so a gap in the data reads as a gap on the page.
    """
    cen = np.asarray(cen, dtype=float)
    val = np.asarray(val, dtype=float)
    if cen.size < 2:
        return cen, val
    gap = np.flatnonzero(np.diff(cen) > step_bp * GAP_FACTOR) + 1
    if gap.size == 0:
        return cen, val
    return (np.insert(cen, gap, np.nan), np.insert(val, gap, np.nan))


def isolated_windows(cen, step_bp):
    """Boolean mask of windows that the broken trace joins to nothing.

    After ``break_gaps`` a window is joined to the one before it, and to the one
    after it, only when their centres are at most GAP_FACTOR steps apart. A
    window joined on neither side is a single finite point between two NaNs,
    and ``ax.plot`` draws a line segment only between two points, so without a
    marker it does not appear on the page at all. The first and last windows
    count as unjoined on their outer side; a chromosome with one window has
    only an isolated one.
    """
    cen = np.asarray(cen, dtype=float)
    if cen.size == 0:
        return np.zeros(0, dtype=bool)
    if cen.size == 1:
        return np.ones(1, dtype=bool)
    joined = np.diff(cen) <= step_bp * GAP_FACTOR      # window i joined to i+1
    to_left = np.concatenate(([False], joined))
    to_right = np.concatenate((joined, [False]))
    return ~(to_left | to_right)


def panel_window_mask(cen, seg_s, seg_e, first, last):
    """Windows that belong to the panel spanning [seg_s, seg_e).

    Every window is owned by exactly one panel of its chromosome, so markers
    drawn per panel are drawn - and counted - once. Window centres can fall
    outside the chromosome (the last window starts within one step of the last
    gene, so its centre can lie up to half a window past it), so the first slice
    also owns everything before it and the last slice everything after it. An
    unwrapped chromosome is one panel that is both, and owns every window.
    """
    cen = np.asarray(cen, dtype=float)
    mine = np.ones(cen.size, dtype=bool)
    if not first:
        mine &= cen >= seg_s
    if not last:
        mine &= cen < seg_e
    return mine


def relabel(chroms, pattern, repl):
    """{chrom: panel label} after ``re.sub(pattern, repl, chrom, count=1)``.

    Exits on a bad pattern or replacement, and on two chromosomes ending up with
    the same label, since that would make two panels indistinguishable. A label
    the substitution empties is kept at its full name, with a warning.
    """
    try:
        pat = re.compile(pattern)
    except re.error as exc:
        sys.exit(f'ERROR: --label-sub pattern {pattern!r} is not a valid regex: '
                 f'{exc}')
    out = {}
    for ch in chroms:
        try:
            new = pat.sub(repl, ch, count=1).strip()
        except (re.error, IndexError) as exc:
            sys.exit(f'ERROR: --label-sub replacement {repl!r} failed on '
                     f'{ch!r}: {exc}')
        if not new:
            logging.warning("  --label-sub empties '%s'; keeping the full name",
                            ch)
            new = ch
        out[ch] = new
    seen = {}
    for ch, lab in out.items():
        seen.setdefault(lab, []).append(ch)
    clash = {lab: chs for lab, chs in seen.items() if len(chs) > 1}
    if clash:
        sys.exit('ERROR: --label-sub gives two or more chromosomes the same '
                 'label: ' + '; '.join(f"'{lab}' <- {', '.join(chs)}"
                                       for lab, chs in clash.items()))
    return out


def x_label_collisions(fig, axes_rows, min_gap_mm=LABEL_GAP_MM):
    """Adjacent panel labels in a row that overlap, or nearly do, on the page.

    audit() checks type size and legend placement but not text against text,
    which is how 15-character scaffold IDs under 14 mm panels reached S4. Draws
    the canvas first, because axis-label positions are only set during a draw.
    Returns [(row, left label, right label, gap in mm)] for every pair closer
    than ``min_gap_mm``; a negative gap is an overlap.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    px_per_mm = fig.dpi / 25.4
    hits = []
    for r, axes in enumerate(axes_rows, 1):
        boxes = [(ax.xaxis.label.get_text(),
                  ax.xaxis.label.get_window_extent(rend))
                 for ax in axes if ax.xaxis.label.get_text().strip()]
        for (t1, b1), (t2, b2) in zip(boxes, boxes[1:]):
            gap = (b2.x0 - b1.x1) / px_per_mm
            if gap < min_gap_mm:
                hits.append((r, t1, t2, gap))
    return hits


def pack_rows(panels, row_mb=None, ncols=5, break_before=None):
    """Group panels into rows.

    With `row_mb`, greedy first-fit on cumulative sequence length, so rows come
    out nearly equal width - panel width is proportional to sequence length, and
    a fixed panel count gives ragged rows whenever chromosome lengths vary.
    Without it, `ncols` panels per row.

    `break_before` is a regex; the FIRST matching chromosome starts a new row,
    so '^chrX' keeps the sex chromosomes together instead of stranding chrX1 at
    the end of an autosome row.
    """
    if not panels:
        return []
    if not row_mb or row_mb <= 0:
        n = max(1, int(ncols))
        return [panels[i:i + n] for i in range(0, len(panels), n)]
    try:
        break_pat = re.compile(break_before) if break_before else None
    except re.error as exc:
        raise ValueError(f'--break-before is not a valid regex: {exc}') from exc

    row_bp = row_mb * 1e6
    rows, cur, cur_bp, broken = [], [], 0.0, False
    for pan in panels:
        chrom, seg_s, seg_e = pan
        width = seg_e - seg_s
        forced = bool(break_pat and break_pat.search(chrom) and not broken)
        if forced:
            broken = True
        if cur and (forced or cur_bp + width > row_bp):
            rows.append(cur)
            cur, cur_bp = [], 0.0
        cur.append(pan)
        cur_bp += width
    if cur:
        rows.append(cur)
    assert sum(len(r) for r in rows) == len(panels), \
        'pack_rows lost or duplicated panels'
    return rows


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Whole-genome rate trace in a to-scale chromosome grid '
                    '(Figure 3; supplementary fig. S6 for dN, dS and the other '
                    'species).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--kaks-table', required=True)
    p.add_argument('--gff', required=True)
    p.add_argument('--elim-bed', default=None,
                   help='Eliminated regions BED. Omit for a species with no '
                        'region calls; nothing is shaded.')
    p.add_argument('--chrom-sizes', default=None,
                   help='Two-column chrom.sizes or .fai. Strongly recommended: '
                        'without it panels are sized from gene extent and the '
                        'chromosomes are not to scale.')
    p.add_argument('--gff-feature', default=None,
                   help='Force the GFF feature supplying positions (mRNA for '
                        'miniprot, transcript for Augustus). Default: whichever '
                        'matches the most table identifiers.')
    p.add_argument('--species', default='')
    p.add_argument('--metric', default='Ka_Ks', choices=['Ka_Ks', 'Ka', 'Ks'])
    p.add_argument('--gene-col', default='gene1')
    p.add_argument('--region-col', default='region1',
                   help="Column carrying retained/eliminated. Pass 'none' to "
                        "draw every gene in one neutral colour -- use when the "
                        "only labels available come from another species' "
                        "orthologs and would assert a call this genome does not "
                        "support.")
    p.add_argument('--chrom-col', default=None)
    p.add_argument('--window', type=int, default=500, help='Window, kb')
    p.add_argument('--step', type=int, default=100, help='Step, kb')
    p.add_argument('--min-genes', type=int, default=10,
                   help='Minimum genes in a window for it to be drawn '
                        '(default: 10). At 1 a single gene is drawn as if it '
                        'were a median comparable to the rest of the line.')
    p.add_argument('--stat', default='median', choices=['median', 'mean'])
    p.add_argument('--ncols', type=int, default=5,
                   help='Panels per row; ignored when --row-mb is given')
    p.add_argument('--row-mb', type=float, default=None,
                   help='Pack rows to about this many Mb instead of a fixed '
                        'panel count. Takes precedence over --ncols. Try 75 for '
                        'Ascaris; the summary prints the actual Mb per row.')
    p.add_argument('--break-before', default=None,
                   help="Regex; force a new row at its first match, e.g. "
                        "'^chrX'. Only meaningful with --row-mb.")
    p.add_argument('--wrap-mb', type=float, nargs='?', const=100.0, default=None,
                   help='Cut any chromosome longer than this across rows. '
                        'Needed for Parascaris (one 240 Mb germline '
                        'chromosome). X axes stay in absolute coordinates.')
    p.add_argument('--min-trailing-frac', type=float, default=0.25,
                   help='A final wrap slice shorter than this fraction of a full '
                        'slice is folded into the panel before it rather than '
                        'taking a whole row as a hairline.')
    p.add_argument('--ymax', type=float, default=None,
                   help='Fix the y-limit. Set it explicitly when comparing runs: '
                        'the default is a per-run quantile, so two species drawn '
                        'separately get different scales and look comparable when '
                        'they are not.')
    p.add_argument('--ymax-pct', type=float, default=99.0,
                   help='Percentile setting the default y-limit (default: 99)')
    p.add_argument('--min-chrom-mb', type=float, default=0.0)
    p.add_argument('--exclude-chroms', default=None)
    p.add_argument('--include-chroms', default=None)
    p.add_argument('--max-chroms', type=int, default=None)
    p.add_argument('--no-sex-highlight', action='store_true',
                   help="Do not accent chrX/chrY. Needed for Parascaris, whose "
                        "single germline chromosome is called chrX but is the "
                        "whole genome.")
    p.add_argument('--label-sub', nargs=2, metavar=('PATTERN', 'REPL'),
                   default=None,
                   help="Panel labels become re.sub(PATTERN, REPL, name, "
                        "count=1). Only the label changes; every table keeps "
                        "the full sequence name. For Baylisascaris: "
                        "--label-sub '^GWHBAVG0*' 'Scaffold ' turns "
                        "GWHBAVG00000005 into 'Scaffold 5'; for Parascaris, "
                        "--label-sub '^chrX$' 'Germline chromosome'.")
    p.add_argument('--width', choices=['single', 'onehalf', 'double'],
                   default='double')
    p.add_argument('--row-height-mm', type=float, default=26.0)
    p.add_argument('--no-legend', action='store_true')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        tag = f'_{args.species.lower()}' if args.species else ''
        args.outdir = os.path.join(
            'results',
            f'genome_trace_{datetime.now().strftime("%Y%m%d_%H%M%S")}{tag}')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, f'{args.metric}_genome_trace')
    # fontTools logs ~60 INFO lines per subsetted table while writing SVG/PDF and
    # buries everything else; matplotlib warns once per glyph run when Arial is
    # absent from the compute node. Neither affects the output.
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'fontTools.subset.timer', 'matplotlib.font_manager',
                  'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('genome_trace.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()

    # --- table ------------------------------------------------------------
    df = pd.read_csv(args.kaks_table, sep='\t')
    logging.info('Loaded %s rows from %s', f'{len(df):,}',
                 os.path.basename(args.kaks_table))
    need = [args.gene_col, args.metric]
    if args.region_col.lower() != 'none':
        need.append(args.region_col)
    for col in need:
        if col not in df.columns:
            sys.exit(f"ERROR: column '{col}' not in {args.kaks_table}\n"
                     f'  columns: {list(df.columns)}')

    elim_regions = {}
    if args.elim_bed:
        elim_regions = parse_elim_bed(args.elim_bed)
        logging.info('Eliminated regions on %d sequences', len(elim_regions))
    else:
        logging.warning('No --elim-bed: nothing will be shaded.')

    logging.info('Parsing gene positions from %s', os.path.basename(args.gff))
    by_feature, feat_counts = get_gene_positions(args.gff)
    table_ids = set(df[args.gene_col].astype(str))
    gene_pos, used_ft, n_match = pick_gene_positions(by_feature, table_ids,
                                                     forced=args.gff_feature)
    logging.info("Using feature '%s': %s positions, %s of %s table genes matched",
                 used_ft, f'{len(gene_pos):,}', f'{n_match:,}',
                 f'{len(table_ids):,}')
    if n_match == 0:
        ex_gff = next((sorted(d)[:5] for d in by_feature.values() if d), [])
        sys.exit(
            f"ERROR: no identifier in '{args.gene_col}' matches any GFF ID.\n"
            f'  table:    {sorted(table_ids)[:5]}\n'
            f'  GFF:      {ex_gff or "(none parsed)"}\n'
            f'  features: {sorted(feat_counts.items(), key=lambda kv: -kv[1])[:6]}\n'
            f'  Pass --gff-feature if the type you need is not in that list, or '
            f'check that --gff matches the annotation the proteins came from and '
            f'that --gene-col names the right species column.')

    df['_mid'] = df[args.gene_col].map(lambda g: gene_pos.get(g, (None, None))[1])
    if args.chrom_col and args.chrom_col in df.columns:
        df['_chrom'] = df[args.chrom_col]
    else:
        df['_chrom'] = df[args.gene_col].map(
            lambda g: gene_pos.get(g, (None, None))[0])
    df = df.dropna(subset=['_mid', '_chrom', args.metric])
    df['_mid'] = df['_mid'].astype(int)

    use_regions = (args.region_col.lower() != 'none'
                   and args.region_col in df.columns)
    if args.region_col.lower() != 'none' and args.region_col not in df.columns:
        logging.warning("Region column '%s' absent; drawing one colour",
                        args.region_col)
    if use_regions:
        df = df[df[args.region_col].isin(['retained', 'eliminated'])]
        n_ret = int((df[args.region_col] == 'retained').sum())
        n_elim = int((df[args.region_col] == 'eliminated').sum())
        logging.info('  %s genes (%s retained, %s eliminated)', f'{len(df):,}',
                     f'{n_ret:,}', f'{n_elim:,}')
    else:
        n_ret, n_elim = len(df), 0
        logging.info('  %s genes, not split by region', f'{len(df):,}')
    if df.empty:
        sys.exit('ERROR: no genes left after mapping positions.')

    # --- chromosomes ------------------------------------------------------
    # Every sequence a filter removes is recorded with the filter that removed
    # it, so the summary can say why a scaffold is missing from the figure.
    dropped = []                              # (filter, sequence, length bp)
    genes_by_chrom = df['_chrom'].value_counts()
    if args.chrom_sizes:
        sizes = load_chrom_sizes(args.chrom_sizes)
        logging.info('Loaded %d sequence lengths', len(sizes))
        have_genes = set(df['_chrom'])
        chroms = [c for c in sizes if c in have_genes]
        dropped += [('no plotted gene', c, sizes[c]) for c in sizes
                    if c not in have_genes]
        not_sized = sorted(have_genes - set(sizes))
        if not_sized:
            logging.warning('  %d sequence(s) carrying %s genes are absent from '
                            '--chrom-sizes and are not drawn (e.g. %s)',
                            len(not_sized),
                            f'{int(genes_by_chrom[not_sized].sum()):,}',
                            ', '.join(not_sized[:3]))
            dropped += [('absent from --chrom-sizes', c, np.nan)
                        for c in not_sized]
        to_scale = True
    else:
        logging.warning('No --chrom-sizes: panels are sized from gene extent and '
                        'the chromosomes are NOT to scale.')
        sizes = {c: int(g['_mid'].max()) for c, g in df.groupby('_chrom')}
        chroms = list(sizes)
        to_scale = False

    def drop(flag, gone):
        """Record the sequences one filter removed, and log what they held."""
        dropped.extend((flag, c, sizes[c]) for c in gone)
        if gone:
            logging.info('  %s removes %d sequence(s), %.2f Mb, %s genes', flag,
                         len(gone), sum(sizes[c] for c in gone) / 1e6,
                         f'{int(genes_by_chrom.reindex(gone).fillna(0).sum()):,}')

    for flag, pattern, keep in (('--include-chroms', args.include_chroms, True),
                                ('--exclude-chroms', args.exclude_chroms, False)):
        if pattern:
            pat = re.compile(pattern)
            before = len(chroms)
            gone = [c for c in chroms if bool(pat.search(c)) is not keep]
            chroms = [c for c in chroms if bool(pat.search(c)) is keep]
            logging.info("%s '%s': %d -> %d", flag, pattern, before, len(chroms))
            drop(f"{flag} '{pattern}'", gone)
    if args.min_chrom_mb > 0:
        gone = [c for c in chroms if sizes[c] < args.min_chrom_mb * 1e6]
        chroms = [c for c in chroms if sizes[c] >= args.min_chrom_mb * 1e6]
        drop(f'--min-chrom-mb {args.min_chrom_mb:g}', gone)
    if args.max_chroms:
        ranked = sorted(chroms, key=lambda c: -sizes[c])
        chroms, gone = ranked[:args.max_chroms], ranked[args.max_chroms:]
        drop(f'--max-chroms {args.max_chroms}', gone)
    chroms = sorted(chroms, key=natural_key)
    if not chroms:
        sys.exit('ERROR: no chromosomes left after filtering.')
    logging.info('%d chromosomes: %s', len(chroms), ', '.join(chroms))

    # Only genes on the sequences drawn count from here on: the legend n, the
    # summary and the default y-limit. They used to include genes on the
    # sequences removed above, which the figure does not show.
    n_genes_mapped = len(df)
    df = df[df['_chrom'].isin(chroms)]
    n_genes_off = n_genes_mapped - len(df)
    if use_regions:
        n_ret = int((df[args.region_col] == 'retained').sum())
        n_elim = int((df[args.region_col] == 'eliminated').sum())
    else:
        n_ret, n_elim = len(df), 0
    logging.info('  %s genes on the sequences drawn (%s on sequences not drawn)',
                 f'{len(df):,}', f'{n_genes_off:,}')

    # Panel labels. Only the drawn label changes; the windows table, summary
    # and provenance keep the full sequence names.
    label_of = {ch: ch for ch in chroms}
    if args.label_sub:
        label_of = relabel(chroms, *args.label_sub)
        n_changed = sum(label_of[c] != c for c in chroms)
        if n_changed == 0:
            logging.warning("--label-sub %r matched no chromosome name; labels "
                            "are unchanged", args.label_sub[0])
        else:
            ex = next(c for c in chroms if label_of[c] != c)
            logging.info("--label-sub: %d of %d labels shortened (e.g. %s -> "
                         "'%s')", n_changed, len(chroms), ex, label_of[ex])

    # --- panels -----------------------------------------------------------
    wrap_bp = args.wrap_mb * 1e6 if args.wrap_mb else None
    panels, folded = [], []
    for ch in chroms:
        L = float(sizes[ch])
        if wrap_bp and L > wrap_bp:
            bounds, seg = [], 0.0
            while seg < L:
                bounds.append((seg, min(seg + wrap_bp, L)))
                seg += wrap_bp
            # 240.073 Mb wrapped at 60 Mb ends with a 73 kb slice, which renders
            # as an empty hairline and takes a whole row.
            if len(bounds) > 1:
                s_last, e_last = bounds[-1]
                if (e_last - s_last) < args.min_trailing_frac * wrap_bp:
                    bounds[-2] = (bounds[-2][0], e_last)
                    bounds.pop()
                    folded.append((ch, (e_last - s_last) / 1e6))
            panels.extend((ch, s, e) for s, e in bounds)
        else:
            panels.append((ch, 0.0, L))
    for ch, mbs in folded:
        logging.info('  %s: trailing %.3f Mb folded into the previous panel',
                     ch, mbs)
    n_panels_of, slices_of = {}, {}
    for ch, s, e in panels:
        n_panels_of[ch] = n_panels_of.get(ch, 0) + 1
        slices_of.setdefault(ch, []).append((s, e))

    # --- windows ----------------------------------------------------------
    window_bp, step_bp = args.window * 1000, args.step * 1000
    window_lines, win_rows, dropped_of, iso_of = {}, [], {}, {}
    no_trace = []                 # drawn as points only: too few genes
    for ch in chroms:
        sub = df[df['_chrom'] == ch]
        if len(sub) <= 5:
            no_trace.append((ch, len(sub)))
            continue
        cen, val, cnt, ndrop = sliding_window(sub['_mid'].to_numpy(),
                                              sub[args.metric].to_numpy(),
                                              window_bp, step_bp,
                                              args.min_genes, args.stat)
        if ndrop:
            dropped_of[ch] = ndrop
        if cen.size:
            window_lines[ch] = (cen, val, cnt)
            # Computed once per chromosome, before the panels exist, so a wrap
            # point cannot make a window look isolated.
            iso_of[ch] = isolated_windows(cen, step_bp)
            win_rows.append(pd.DataFrame({'chrom': ch, 'centre_bp': cen,
                                          'value': val, 'n_genes': cnt,
                                          'metric': args.metric,
                                          'isolated': iso_of[ch].astype(int)}))
    n_iso = sum(int(m.sum()) for m in iso_of.values())
    if n_iso:
        logging.info('  %d isolated window(s), both neighbours below '
                     '--min-genes: drawn as dots', n_iso)
    if win_rows:
        pd.concat(win_rows).to_csv(
            os.path.join(args.outdir,
                         f'{args.metric}_genome_trace_windows.tsv'),
            sep='\t', index=False, float_format='%.6g')

    # --- layout -----------------------------------------------------------
    rows = pack_rows(panels, args.row_mb, args.ncols, args.break_before)
    row_totals = [sum(e - s for _, s, e in rw) for rw in rows]
    for r, (rw, tot) in enumerate(zip(rows, row_totals), 1):
        logging.info('  row %d: %6.1f Mb  [%s]', r, tot / 1e6,
                     ', '.join(c for c, _, _ in rw))

    fig_w_in = COL[args.width]
    legend_mm = 0.0 if args.no_legend else 7.0
    # Bottom: tick labels, the panel labels, and the position title under them.
    top_mm, bottom_mm, left_mm, right_mm = 2.0 + legend_mm, 13.5, 14.0, 2.0
    gap_mm = 9.0
    total_h_mm = (top_mm + bottom_mm + args.row_height_mm * len(rows)
                  + gap_mm * (len(rows) - 1))
    if total_h_mm > MAX_HEIGHT_MM:
        logging.warning('Figure is %.0f mm tall, over the %.0f mm type area. '
                        'Reduce --row-height-mm or pack fewer rows.',
                        total_h_mm, MAX_HEIGHT_MM)
    fig = plt.figure(figsize=(fig_w_in, mm(total_h_mm)))
    fig.set_layout_engine('none')
    h_in = mm(total_h_mm)

    ymax = (args.ymax if args.ymax else
            float(np.nanpercentile(df[args.metric], args.ymax_pct) * 1.10))
    logging.info('  y-limit %.4f (%s)', ymax,
                 'given' if args.ymax else f'{args.ymax_pct:g}th pct x 1.10')

    ylabel = METRIC_LABELS.get(args.metric, args.metric)
    max_total = max(row_totals)
    avail_w = fig_w_in - mm(left_mm) - mm(right_mm)
    panel_gap = mm(2.6)

    n_clipped, drew_sex = 0, False
    iso_dot_of, iso_caret_of = {}, {}   # isolated windows: as dots / as carets
    axes_rows = []                      # for the label collision check
    for r, row_panels in enumerate(rows):
        axes_rows.append([])
        # One Mb-per-inch scale across every row, so a short final row occupies
        # proportionally less width instead of being stretched to fill.
        row_w = avail_w * (row_totals[r] / max_total)
        inner = row_w - panel_gap * (len(row_panels) - 1)
        top = h_in - mm(top_mm) - r * (mm(args.row_height_mm) + mm(gap_mm))
        x = mm(left_mm)
        for c, (chrom, seg_s, seg_e) in enumerate(row_panels):
            w = inner * ((seg_e - seg_s) / row_totals[r])
            ax = fig.add_axes([x / fig_w_in,
                               (top - mm(args.row_height_mm)) / h_in,
                               w / fig_w_in, mm(args.row_height_mm) / h_in])
            axes_rows[-1].append(ax)
            x += w + panel_gap

            is_sex = (not args.no_sex_highlight
                      and chrom.upper().startswith(('CHRX', 'CHRY')))
            sub = df[(df['_chrom'] == chrom) & (df['_mid'] >= seg_s)
                     & (df['_mid'] < seg_e)]

            # An eliminated tract is drawn twice. The span is to scale and
            # carries the extent; the rug bar under it is a MARKER with a
            # minimum on-page width, because a 40 kb tract on a 30 mm panel is
            # 0.01 mm wide and disappears at print size no matter what colour
            # it is. Widening the span itself would misreport the extent, so
            # the guarantee goes on the marker and the caption says so.
            mb_per_in = (seg_e - seg_s) / 1e6 / max(w, 1e-6)
            min_mb = mb_per_in * (MIN_TRACT_MM / 25.4)
            for es, ee in elim_regions.get(chrom, []):
                if ee >= seg_s and es <= seg_e:
                    lo, hi = max(es, seg_s) / 1e6, min(ee, seg_e) / 1e6
                    ax.axvspan(lo, hi, zorder=0, lw=0,
                               color=to_rgba(REGION['eliminated'], 0.16))
                    span = max(hi - lo, min_mb)
                    mid = (lo + hi) / 2.0
                    ax.add_patch(Rectangle(
                        (mid - span / 2.0, 0.0), span, 0.055,
                        transform=ax.get_xaxis_transform(), clip_on=False,
                        facecolor=REGION['eliminated'], edgecolor='none',
                        zorder=6))

            if use_regions:
                for name, size in (('retained', 1.0), ('eliminated', 1.8)):
                    s = sub[sub[args.region_col] == name]
                    ax.scatter(s['_mid'] / 1e6, s[args.metric], s=size,
                               color=REGION[name], alpha=0.35, linewidths=0,
                               zorder=1 if name == 'retained' else 2,
                               rasterized=True)
            else:
                ax.scatter(sub['_mid'] / 1e6, sub[args.metric], s=1.0,
                           color=NEUTRAL['ink_light'], alpha=0.35, linewidths=0,
                           zorder=1, rasterized=True)

            if chrom in window_lines:
                cen, val, _cnt = window_lines[chrom]
                trace_col = NEUTRAL['sex_chrom'] if is_sex else NEUTRAL['ink']
                if is_sex:
                    drew_sex = True
                gx, gy = break_gaps(cen, val, step_bp)
                ax.plot(gx / 1e6, gy, '-', lw=0.9, zorder=3, color=trace_col)
                # The line above is the whole chromosome, clipped by the
                # x-limits. The markers below are not clipped (carets) or are
                # counted (both), so they take only this panel's windows -- or
                # a wrapped chromosome draws and counts each one per slice.
                slices = slices_of[chrom]
                mine = panel_window_mask(cen, seg_s, seg_e,
                                         first=(seg_s, seg_e) == slices[0],
                                         last=(seg_s, seg_e) == slices[-1])
                # A window above the y-cap silently vanishes, which reads as
                # missing data rather than as clipping. Caret it on the spine.
                over = mine & (val > ymax)
                if over.any():
                    ax.plot(cen[over] / 1e6, np.full(int(over.sum()), ymax),
                            marker='^', ls='', ms=2.2, mew=0, clip_on=False,
                            zorder=7, color=trace_col)
                    n_clipped += int(over.sum())
                # A window whose neighbours on both sides fell below
                # --min-genes is a lone point between two breaks, which the
                # line cannot draw. Mark it with a dot; one above the y-cap
                # already has its caret.
                iso = mine & iso_of[chrom]
                dots = iso & ~over
                if dots.any():
                    ax.plot(cen[dots] / 1e6, val[dots], marker='o', ls='',
                            ms=ISO_MS, mew=0, zorder=3, color=trace_col)
                iso_dot_of[chrom] = iso_dot_of.get(chrom, 0) + int(dots.sum())
                iso_caret_of[chrom] = (iso_caret_of.get(chrom, 0)
                                       + int((iso & over).sum()))

            ax.set_xlim(seg_s / 1e6, seg_e / 1e6)
            ax.set_ylim(0, ymax)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
            ax.tick_params(axis='x', pad=1.5)
            label = (label_of[chrom] if n_panels_of[chrom] == 1
                     else f'{label_of[chrom]} {seg_s / 1e6:g}–{seg_e / 1e6:g}')
            ax.set_xlabel(label, labelpad=1.5,
                          color=NEUTRAL['sex_chrom'] if is_sex else 'black')
            if c == 0:
                ax.set_ylabel(ylabel, labelpad=2)
            else:
                ax.tick_params(axis='y', labelleft=False, length=0)

    # One x-axis title for the whole grid, centred on the full row width.
    pos_title = fig.text((mm(left_mm) + avail_w / 2.0) / fig_w_in,
                         mm(1.0) / h_in, POSITION_LABEL, ha='center',
                         va='bottom', fontsize=plt.rcParams['axes.labelsize'])

    if not args.no_legend:
        if use_regions:
            handles = [
                Line2D([0], [0], marker='o', ls='', color=REGION['retained'],
                       label=f'Retained (n = {n_ret:,})', markersize=3),
                Line2D([0], [0], marker='o', ls='', color=REGION['eliminated'],
                       label=f'Eliminated (n = {n_elim:,})', markersize=3),
                Line2D([0], [0], color=NEUTRAL['ink'], lw=1.0,
                       label=f'{args.stat.capitalize()}, {args.window} kb window')]
            handles.append(Line2D([0], [0], marker='s', ls='', markersize=3,
                                  color=REGION['eliminated'],
                                  label='Eliminated tract'))
        else:
            handles = [
                Line2D([0], [0], marker='o', ls='', color=NEUTRAL['ink_light'],
                       label=f'Genes (n = {n_ret:,})', markersize=3),
                Line2D([0], [0], color=NEUTRAL['ink'], lw=1.0,
                       label=f'{args.stat.capitalize()}, {args.window} kb window')]
        # Two line colours are drawn whenever a sex chromosome is present, and
        # until now only one of them had a key. The gold tick labels direct-label
        # it, but a reader still has to infer that the two encode the same
        # quantity.
        if drew_sex:
            handles.append(Line2D([0], [0], color=NEUTRAL['sex_chrom'], lw=1.0,
                                  label='Sex chromosome'))
        # No suptitle: the caption names the species and the metric, and a title
        # inside the artwork is one more thing to keep consistent with it.
        fig.legend(handles=handles, loc='upper center',
                   bbox_to_anchor=(0.5, 1.0 - mm(1.0) / h_in),
                   ncol=len(handles), frameon=False, handletextpad=0.4,
                   columnspacing=1.4, borderaxespad=0.0)

    # A genome grid has many small panels by construction, so audit()'s minimum
    # axes size fires once per panel. Collapse those into one honest line rather
    # than 24 warnings that bury anything else the audit finds - but do not
    # silence it: a panel under the floor is still a legibility risk and the
    # count and the smallest dimension are reported.
    problems = audit(fig, verbose=False)
    small = [p for p in problems if 'floor is' in p]
    for problem in problems:
        if problem not in small:
            logging.warning('  audit: %s', problem)
    if small:
        logging.warning('  audit: %d of %d panels are below the %g mm minimum '
                        '(row height %g mm). Deliberate for a dense grid, but '
                        'check the proof at final size.',
                        len(small), len(panels), MIN_AXES_MM,
                        args.row_height_mm)

    # Text against text, which audit() does not check.
    collisions = x_label_collisions(fig, axes_rows)
    # The position title must sit below the last row's panel labels.
    rend = fig.canvas.get_renderer()
    pos_box = pos_title.get_window_extent(rend)
    lowest = min((ax.xaxis.label.get_window_extent(rend).y0
                  for ax in axes_rows[-1] if ax.xaxis.label.get_text()),
                 default=pos_box.y1 + 1)
    if pos_box.y1 > lowest:
        logging.warning("  labels: '%s' overlaps the last row's panel labels "
                        "by %.1f mm", POSITION_LABEL,
                        (pos_box.y1 - lowest) / fig.dpi * 25.4)
    for r, t1, t2, gap in collisions:
        logging.warning("  labels: row %d, '%s' and '%s' %s", r, t1, t2,
                        f'overlap by {-gap:.1f} mm' if gap < 0
                        else f'are only {gap:.1f} mm apart')
    if collisions:
        logging.warning('  labels: %d adjacent pair(s) collide; shorten them '
                        'with --label-sub or pack fewer panels per row',
                        len(collisions))
    else:
        logging.info('  labels: no adjacent panel labels closer than %g mm',
                     LABEL_GAP_MM)

    # Every isolated window must have been drawn once, as a dot or a caret. A
    # shortfall means a window owned by no panel, which should not happen.
    n_iso_dot = sum(iso_dot_of.values())
    n_iso_caret = sum(iso_caret_of.values())
    if n_iso_dot + n_iso_caret != n_iso:
        logging.warning('  isolated windows: %d found but %d drawn (%d dots, '
                        '%d carets)', n_iso, n_iso_dot + n_iso_caret,
                        n_iso_dot, n_iso_caret)
    stem = (f'{args.metric}_genome_trace'
            f'{"_" + args.species if args.species else ""}')
    paths = save(fig, args.outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))

    # --- summary ----------------------------------------------------------
    per_chrom = (df.groupby('_chrom')[args.metric]
                   .agg(['count', 'median']).reindex(chroms))
    L = [f'{args.metric} along {args.species or "chromosomes"}', '',
         f'Table            {os.path.basename(args.kaks_table)}',
         f'Genes plotted    {len(df):,}'
         + (f' ({n_ret:,} retained, {n_elim:,} eliminated)' if use_regions
            else ' (not split by region)')
         + (f'; {n_genes_off:,} more on sequences not drawn' if n_genes_off
            else ''),
         f'GFF feature      {used_ft} ({n_match:,}/{len(table_ids):,} matched)',
         f'Chromosomes      {len(chroms)}  ({len(panels)} panels)',
         f'To scale         {"yes" if to_scale else "NO -- no --chrom-sizes"}',
         f'Window           {args.window} kb, step {args.step} kb, '
         f'min {args.min_genes} genes, {args.stat}',
         f'y-limit          {ymax:.4f}'
         + ('  (given)' if args.ymax else
            f'  ({args.ymax_pct:g}th percentile x 1.10 of the genes drawn)'),
         f'Canvas           {args.width} column, {total_h_mm:.0f} mm tall', '',
         'Rows', '']
    for r, (rw, tot) in enumerate(zip(rows, row_totals), 1):
        L.append(f'  {r}  {tot / 1e6:6.1f} Mb   '
                 f'{", ".join(c for c, _, _ in rw)}')
    cw = max(10, max(len(c) for c in chroms))
    L += ['', 'Per chromosome', '',
          f'  {"chrom":<{cw}} {"genes":>7} {"median":>9} {"length Mb":>10}']
    for ch in chroms:
        row = per_chrom.loc[ch]
        L.append(f'  {ch:<{cw}} {int(row["count"]):>7,} {row["median"]:>9.4f} '
                 f'{sizes[ch] / 1e6:>10.2f}')
    if no_trace:
        L += ['', '  No window trace (5 genes or fewer; points only): '
              + ', '.join(f'{c} ({n})' for c, n in no_trace)]

    # Every sequence not drawn, grouped by the filter that removed it. A
    # legend has to be able to say why a scaffold is missing.
    L += ['', 'Sequences not drawn', '']
    if not dropped:
        L.append('  none')
    by_flag = {}
    for flag, ch, length in dropped:
        by_flag.setdefault(flag, []).append((ch, length))
    sub_pat = None
    if args.label_sub:
        sub_pat = re.compile(args.label_sub[0])
    for flag, items in by_flag.items():
        lens = np.array([v for _, v in items], dtype=float)
        n_g = int(genes_by_chrom.reindex([c for c, _ in items])
                  .fillna(0).sum())
        L.append(f'  {flag}: {len(items):,} sequence(s), '
                 f'{np.nansum(lens) / 1e6:.2f} Mb, {n_g:,} genes')
        items = sorted(items, key=lambda t: -(t[1] if np.isfinite(t[1]) else 0))
        show = items if len(items) <= 10 else items[:5]
        for ch, length in show:
            lab = sub_pat.sub(args.label_sub[1], ch, count=1).strip() \
                if sub_pat else ''
            L.append(f'    {ch:<{cw}} '
                     + (f'{length / 1e6:>8.2f} Mb' if np.isfinite(length)
                        else f'{"?":>8} Mb')
                     + f'  {int(genes_by_chrom.get(ch, 0)):>6,} genes'
                     + (f'  (label would be {lab!r})' if lab and lab != ch
                        else ''))
        if len(items) > len(show):
            L.append(f'    ... and {len(items) - len(show):,} smaller')
    if n_genes_off:
        L.append(f'  The legend n counts only the {len(df):,} genes on the '
                 f'sequences drawn, not these {n_genes_off:,}.')
    # No silent caps: both the floor and the y-cap remove things from the page,
    # and neither is visible in the artwork as a number.
    n_drop = sum(dropped_of.values())
    L += ['', 'Truncation', '',
          f'  Windows dropped below --min-genes {args.min_genes}: {n_drop}'
          + (f'  ({", ".join(f"{c} {n}" for c, n in sorted(dropped_of.items()))})'
             if dropped_of else '')]
    if n_drop:
        L += ['  These are almost all partial windows at chromosome ends. The',
              '  trace is BROKEN with a gap there rather than joined across, so',
              '  a straight segment is never drawn over an excluded region.']
    L += [f'  Windows above the y-cap {ymax:g}: {n_clipped}'
          + ('  (marked with a caret on the top spine)' if n_clipped else '')]
    if n_clipped:
        L += ['  State the cap and the count in the caption -- a clipped window',
              '  is otherwise indistinguishable from missing data.']
    iso_list = ', '.join(f'{c} {int(m.sum())}' for c, m in iso_of.items()
                         if m.any())
    L += [f'  Isolated windows, both neighbours below the floor: {n_iso}'
          + (f'  ({iso_list})' if n_iso else '')]
    if n_iso:
        L += [f'  Drawn as dots: {n_iso_dot}; above the y-cap, shown by their '
              f'caret: {n_iso_caret}.',
              '  A line needs two points, so without the dot such a window is not',
              '  drawn at all. Flagged isolated=1 in the windows table. Caption:',
              '  "a window whose neighbours on both sides fall below the floor is',
              '  drawn as a dot".']
    L += ['', 'Panel labels', '',
          f'  Adjacent labels closer than {LABEL_GAP_MM:g} mm: {len(collisions)}']
    for r, t1, t2, gap in collisions:
        L.append(f"    row {r}: '{t1}' / '{t2}'  "
                 + (f'overlap {-gap:.1f} mm' if gap < 0 else f'{gap:.1f} mm apart'))
    if args.label_sub:
        L += [f'  --label-sub {args.label_sub[0]!r} {args.label_sub[1]!r}; '
              'tables keep the full names:']
        L += [f'    {label_of[ch]:<{max(len(v) for v in label_of.values())}}  {ch}'
              for ch in chroms]
    L += ['', 'Notes', '',
          '  Built at exactly the journal column width with sizes in mm, so the',
          '  type prints at the size it is set. The earlier version built a 20',
          '  inch canvas and saved with bbox_inches=tight, which at 185 mm is a',
          '  2.7x reduction -- 13 pt labels printing under 5 pt.', '',
          '  Sliding windows are computed once per chromosome and then clipped',
          '  to each panel, so a wrap point does not create an edge effect that',
          '  looks like structure.', '',
          f'  --min-genes {args.min_genes}: a window holding fewer genes is not',
          '  drawn. The default is 10, matching the per-bin floor in the CBR',
          '  meta-profile so the two figures agree on what counts as a summary.',
          '  It was 3, and on the real Ascaris table the six largest terminal',
          '  windows were medians over 3-7 genes -- including the single window',
          '  that ran off the top of the axis. The most eye-catching feature of',
          '  the figure was its least reliable one.', '']
    if not args.ymax:
        L += ['  The y-limit is a percentile of THIS run. Two species drawn in',
              '  separate runs therefore get different scales. Pass --ymax to',
              '  compare them.', '']
    summary = '\n'.join(L)
    with open(os.path.join(
            args.outdir, f'{args.metric}_genome_trace_summary.txt'), 'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(
            args.outdir, f'{args.metric}_provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(), 'python': platform.python_version(),
                   'packages': {'numpy': np.__version__,
                                'pandas': pd.__version__,
                                'matplotlib': matplotlib.__version__},
                   'argv': sys.argv, 'args': vars(args),
                   'inputs': {'kaks_table': fingerprint(args.kaks_table),
                              'gff': fingerprint(args.gff),
                              'elim_bed': fingerprint(args.elim_bed),
                              'chrom_sizes': fingerprint(args.chrom_sizes)},
                   'gff_feature': used_ft, 'n_matched': int(n_match),
                   'n_genes': int(len(df)),
                   'n_genes_on_sequences_not_drawn': int(n_genes_off),
                   'sequences_not_drawn': [
                       {'filter': f, 'sequence': c,
                        'length_bp': (None if not np.isfinite(v) else int(v)),
                        'genes': int(genes_by_chrom.get(c, 0))}
                       for f, c, v in dropped],
                   'no_trace': [{'sequence': c, 'genes': int(n)}
                                for c, n in no_trace],
                   'chromosomes': chroms,
                   'n_panels': len(panels), 'to_scale': to_scale,
                   'ymax': ymax,
                   'rows_mb': [t / 1e6 for t in row_totals],
                   'n_windows_clipped': int(n_clipped),
                   'n_windows_isolated': int(n_iso),
                   'isolated_as_dots': int(n_iso_dot),
                   'isolated_as_carets': int(n_iso_caret),
                   'panel_labels': (label_of if args.label_sub else None),
                   'label_collisions': [
                       {'row': r, 'left': t1, 'right': t2,
                        'gap_mm': round(float(g), 2)}
                       for r, t1, t2, g in collisions]},
                  fh, indent=2, default=float)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
