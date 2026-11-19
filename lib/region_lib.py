#!/usr/bin/env python3
"""
region_lib.py - shared primitives for everything that talks about regions and CBRs.

Imported by `annotate_regions.py`, `gene_density_cbr.py`, `kaks_landscape.py` and
`cbr_changepoint_relax.py`.

The reason this module exists
----------------------------
Figure 4 has four panels drawn by four scripts, and every one of them needs an
answer to the same three questions: where are the CBRs, which side of one is a
given position on, and how far away is it. When each script answered those
questions in its own code the panels disagreed with each other — one signed
distance eliminated-positive and another eliminated-negative, one counted a
telomere as a breakage region and another did not, one treated two abutting
eliminated tracts as separated by a CBR. Those are not rounding differences; they
change which genes are called CBR-proximal and where the boundary sits on the
x-axis.

So the definitions live here, once. A change to what counts as a CBR propagates to
every panel of Figure 4 in one edit, and the panels cannot drift apart again.

Conventions fixed here
----------------------
**Coordinates are 1-based inclusive everywhere in this module.** GFF is already
1-based inclusive; BED is 0-based half-open and is converted on the way in, once,
in `parse_bed`. Mixing the two silently shifts every region edge by one base.

**Signed distance is positive on the eliminated side.** Eliminated is drawn right
and retained left throughout the paper, and `cbr_changepoint_relax.py` already
fits its slopes in the eliminated-positive frame.

**A CBR is a boundary of a merged eliminated tract that is not a sequence end.**
Two corrections are folded into that sentence — see `derive_cbrs`.
"""

from __future__ import annotations

import hashlib
import logging
import os
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Gene models
# ---------------------------------------------------------------------------

def parse_gff_models(path, feature='mRNA'):
    """Every `feature` record as a row: gene_id, chrom, start, end, strand, target,
    rank, identity.

    GFF coordinates are 1-based and inclusive, and are kept that way.

    `rank` and `identity` come from miniprot's `Rank=` and `Identity=` attributes
    and are what `select_primary_alignment` uses. They are None for GFFs that do
    not carry them, in which case selection falls back to the longest span.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != feature:
                continue
            gid = target = None
            rank = identity = None
            for field in parts[8].split(';'):
                if field.startswith('ID='):
                    gid = field[3:].strip()
                elif field.startswith('Target='):
                    target = field[7:].split()[0].strip()
                elif field.startswith('Rank='):
                    try:
                        rank = int(field[5:])
                    except ValueError:
                        pass
                elif field.startswith('Identity='):
                    try:
                        identity = float(field[9:])
                    except ValueError:
                        pass
            if gid is None and target is not None:
                gid = target
            if not gid:
                continue
            start, end = int(parts[3]), int(parts[4])
            rows.append({'gene_id': gid, 'chrom': parts[0], 'start': start,
                         'end': end, 'strand': parts[6], 'target': target,
                         'rank': rank, 'identity': identity,
                         'genomic_span': end - start + 1})
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit(f'ERROR: no {feature} features found in {path}')
    logging.info('Parsed %s %s records from %s across %d sequences',
                 f'{len(df):,}', feature, os.path.basename(path),
                 df['chrom'].nunique())
    return df


def select_primary_alignment(models, id_col='gene_id'):
    """One row per identifier: the best alignment, not the first one in the file.

    A miniprot query protein can align in more than one place. The earlier code
    took `drop_duplicates(keep='first')`, which is whichever record the file
    happened to list first — usually but not reliably the primary hit, and not
    stable if the GFF is ever re-sorted. A secondary alignment can sit on a
    different chromosome from the primary, so the choice decides that gene's
    region call outright.

    Preference order: lowest `Rank`, then highest `Identity`, then longest span.

    Returns (frame, diagnostic dict). The diagnostic records how many identifiers
    were multiply aligned and how many of those would have been called
    differently by the old first-in-file rule — that count is the size of the
    silent error being fixed.
    """
    n_before = len(models)
    dup_mask = models.duplicated(id_col, keep=False)
    n_multi = int(models.loc[dup_mask, id_col].nunique())

    first = models.drop_duplicates(id_col, keep='first').set_index(id_col)

    m = models.copy()
    m['_rank'] = m['rank'].fillna(np.inf)
    m['_ident'] = -m['identity'].fillna(-np.inf)
    m['_span'] = -m['genomic_span']
    best = (m.sort_values([id_col, '_rank', '_ident', '_span'])
              .drop_duplicates(id_col, keep='first')
              .drop(columns=['_rank', '_ident', '_span'])
              .reset_index(drop=True))

    b = best.set_index(id_col)
    shared = b.index.intersection(first.index)
    changed = int(((b.loc[shared, 'chrom'] != first.loc[shared, 'chrom']) |
                   (b.loc[shared, 'start'] != first.loc[shared, 'start'])).sum())
    moved_chrom = int((b.loc[shared, 'chrom'] != first.loc[shared, 'chrom']).sum())

    diag = {'n_records': n_before, 'n_ids': int(len(best)),
            'n_ids_multiply_aligned': n_multi,
            'n_changed_vs_first_in_file': changed,
            'n_moved_to_another_sequence': moved_chrom,
            'has_rank': bool(models['rank'].notna().any())}
    logging.info('  %s records -> %s identifiers; %s aligned more than once, '
                 '%s pick a different locus than first-in-file (%s a different '
                 'sequence)', f'{n_before:,}', f'{len(best):,}', f'{n_multi:,}',
                 f'{changed:,}', f'{moved_chrom:,}')
    return best, diag


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------

def parse_bed(path, keep_label=False):
    """BED intervals as 1-based inclusive coordinates.

    BED is 0-based half-open [start, end); the 1-based inclusive equivalent is
    [start + 1, end]. Every other coordinate in this codebase comes from a GFF
    and is already 1-based, so the conversion happens once, here.

    With `keep_label`, a 4th column is read as `label` (Ascaris uses it for
    end/internal). Files without one simply come back without the column.
    """
    probe = pd.read_csv(path, sep='\t', header=None, nrows=1, comment='#')
    has_label = keep_label and probe.shape[1] >= 4
    cols = ['chrom', 'start', 'end'] + (['label'] if has_label else [])
    usecols = [0, 1, 2] + ([3] if has_label else [])
    df = pd.read_csv(path, sep='\t', header=None, names=cols, usecols=usecols,
                     comment='#', dtype={0: str})
    df['start'] = df['start'].astype(int) + 1
    df['end'] = df['end'].astype(int)
    bad = df['end'] < df['start']
    if bad.any():
        logging.warning('  %d intervals in %s have end < start and were dropped',
                        int(bad.sum()), os.path.basename(path))
        df = df[~bad]
    logging.info('Parsed %s regions from %s across %d sequences%s',
                 f'{len(df):,}', os.path.basename(path), df['chrom'].nunique(),
                 ' (with labels)' if has_label else '')
    return df.reset_index(drop=True)


def merge_intervals(df):
    """Merge overlapping or bookended intervals per sequence.

    Returns {chrom: (starts, ends)}, sorted float arrays, 1-based inclusive.

    Bookended intervals (end + 1 == next start) merge too. Two abutting eliminated
    tracts are one tract; the join between them is a bookkeeping artefact of how
    the calls were made, not a breakage region, and counting it as one puts a
    spurious CBR in the middle of eliminated sequence.
    """
    out = {}
    for chrom, block in df.groupby('chrom', sort=False):
        b = block.sort_values('start')
        starts, ends = [], []
        for s, e in zip(b['start'].to_numpy(), b['end'].to_numpy()):
            if starts and s <= ends[-1] + 1:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)
        out[chrom] = (np.asarray(starts, float), np.asarray(ends, float))
    n_in = len(df)
    n_out = sum(len(v[0]) for v in out.values())
    if n_out < n_in:
        logging.info('  %s intervals merge into %s (overlapping or bookended)',
                     f'{n_in:,}', f'{n_out:,}')
    return out


def load_chrom_sizes(path, *frames):
    """Sequence lengths from a .fai or two-column TSV.

    Falls back to the largest coordinate present in the supplied frames, which
    underestimates every length. Returns (dict, inferred_flag) so the caller can
    put the caveat in its own summary rather than only in a log nobody reads.
    """
    if path:
        sizes = {}
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.split('\t') if '\t' in line else line.split()
                sizes[parts[0]] = int(parts[1])
        logging.info('Loaded %d sequence lengths from %s', len(sizes),
                     os.path.basename(path))
        return sizes, False

    sizes = {}
    for f in frames:
        if f is None or len(f) == 0:
            continue
        for chrom, e in f.groupby('chrom')['end'].max().items():
            sizes[chrom] = max(sizes.get(chrom, 0), int(e))
    logging.warning('No chromosome sizes given: lengths inferred from the last '
                    'annotated feature on each sequence. Every length is an '
                    'underestimate, so a tract boundary near a real sequence end '
                    'may not be recognised as terminal.')
    return sizes, True


# ---------------------------------------------------------------------------
# CBRs
# ---------------------------------------------------------------------------

def derive_cbrs(elim_merged, chrom_sizes, terminal_tol=1000):
    """CBR coordinates: the boundaries of merged eliminated tracts.

    Two corrections over taking every BED edge:

    1. The tracts are **merged first** (see `merge_intervals`), so the junction
       between two abutting eliminated tracts is not counted.
    2. A boundary within `terminal_tol` of a sequence end is **dropped**. A
       terminal eliminated tract runs to the telomere; its distal edge is where
       the assembly stops, not a breakage region. Keeping it inflates the CBR
       count and creates a phantom retained flank beyond the end of the sequence.

    Boundaries are placed at the half-base junction (start - 0.5, end + 0.5) so
    the first eliminated base has a small positive signed distance and the last
    retained base a small negative one, with nothing landing exactly on zero.

    Returns ({chrom: sorted array}, diagnostic dict).
    """
    out, n_raw, n_terminal, no_length = {}, 0, 0, []
    for chrom, (starts, ends) in elim_merged.items():
        L = chrom_sizes.get(chrom)
        cand = np.concatenate([starts - 0.5, ends + 0.5])
        n_raw += cand.size
        keep = cand > terminal_tol
        if L is not None:
            keep &= cand < (L - terminal_tol)
        else:
            no_length.append(chrom)
        n_terminal += int((~keep).sum())
        cand = np.unique(cand[keep])
        if cand.size:
            out[chrom] = cand
    if no_length:
        logging.warning('  no length for %d sequence(s) (%s...); their distal '
                        'tract boundaries cannot be tested for being terminal',
                        len(no_length), ', '.join(no_length[:3]))
    diag = {'n_boundaries': int(n_raw), 'n_terminal_dropped': int(n_terminal),
            'n_cbrs': int(sum(v.size for v in out.values())),
            'n_sequences_with_cbrs': len(out),
            'terminal_tol': int(terminal_tol)}
    logging.info('  %s CBRs on %d sequences (%s terminal boundaries dropped of %s)',
                 f"{diag['n_cbrs']:,}", diag['n_sequences_with_cbrs'],
                 f"{diag['n_terminal_dropped']:,}", f"{diag['n_boundaries']:,}")
    return out, diag


def nearest_cbr(pos, cbrs):
    """(index, absolute distance) of the nearest CBR. `cbrs` sorted ascending."""
    j = np.searchsorted(cbrs, pos)
    lo = np.clip(j - 1, 0, cbrs.size - 1)
    hi = np.clip(j, 0, cbrs.size - 1)
    d_lo = np.abs(pos - cbrs[lo])
    d_hi = np.abs(cbrs[hi] - pos)
    take_hi = d_hi < d_lo
    return np.where(take_hi, hi, lo), np.where(take_hi, d_hi, d_lo)


def in_intervals(pos, starts, ends):
    """Boolean per position: inside a sorted, non-overlapping interval set?"""
    pos = np.asarray(pos, dtype=float)
    if starts.size == 0:
        return np.zeros(pos.shape, dtype=bool)
    j = np.searchsorted(starts, pos, side='right') - 1
    ok = j >= 0
    j = np.clip(j, 0, starts.size - 1)
    return ok & (pos <= ends[j])


def which_interval(pos, starts, ends):
    """Index of the containing interval per position, or -1. Same rules as
    `in_intervals`; used to look up a per-interval label such as end/internal."""
    pos = np.asarray(pos, dtype=float)
    if starts.size == 0:
        return np.full(pos.shape, -1, dtype=int)
    j = np.searchsorted(starts, pos, side='right') - 1
    ok = (j >= 0)
    jc = np.clip(j, 0, starts.size - 1)
    ok &= (pos <= ends[jc])
    return np.where(ok, jc, -1)


def signed_distance(chrom, pos, cbrs_by_chrom, elim_merged):
    """Signed distance to the nearest CBR, in bp. **Positive = eliminated side.**

    Returns (signed distance, nearest-CBR index). NaN / -1 on sequences carrying
    no CBR; callers should exclude those sequences entirely rather than treating
    the distance as missing at random.
    """
    pos = np.asarray(pos, dtype=float)
    cbrs = cbrs_by_chrom.get(chrom)
    if cbrs is None or cbrs.size == 0:
        return np.full(pos.shape, np.nan), np.full(pos.shape, -1, dtype=int)
    idx, dist = nearest_cbr(pos, cbrs)
    starts, ends = elim_merged.get(chrom, (np.array([]), np.array([])))
    inside = in_intervals(pos, starts, ends)
    return np.where(inside, dist, -dist), idx


def annotate_positions(genes, retained_merged, elim_merged, elim_labels,
                       cbrs_by_chrom, cbr_window=100000):
    """Region call, tract label, CBR distance and signed distance, per gene.

    Vectorised per sequence. The earlier implementation looped every gene against
    every interval with `iterrows`, which is O(genes x regions) and takes minutes
    on *Parascaris*; this is two `searchsorted` calls per sequence.

    Adds:
      region          retained | eliminated | ambiguous | unannotated
      elim_type       the eliminated tract's label, where the BED carries one
      min_cbr_dist    absolute distance to the nearest CBR, bp
      signed_cbr_dist the same, **positive on the eliminated side**
      cbr_proximal    |distance| <= cbr_window

    `ambiguous` means the midpoint falls inside both a retained and an eliminated
    interval, which is a contradiction in the input BEDs rather than a property of
    the gene; it is reported rather than silently resolved.
    """
    g = genes.copy()
    n = len(g)
    region = np.full(n, 'unannotated', dtype=object)
    elim_type = np.full(n, None, dtype=object)
    min_dist = np.full(n, np.nan)
    signed = np.full(n, np.nan)

    for chrom, block in g.groupby('chrom', sort=False):
        idx = block.index.to_numpy()
        pos = block['midpoint'].to_numpy(float)

        rs, re_ = retained_merged.get(chrom, (np.array([]), np.array([])))
        es, ee = elim_merged.get(chrom, (np.array([]), np.array([])))
        in_ret = in_intervals(pos, rs, re_)
        which_el = which_interval(pos, es, ee)
        in_el = which_el >= 0

        r = np.where(in_ret & ~in_el, 'retained',
                     np.where(in_el & ~in_ret, 'eliminated',
                              np.where(in_ret & in_el, 'ambiguous',
                                       'unannotated')))
        region[g.index.get_indexer(idx)] = r

        labels = elim_labels.get(chrom)
        if labels is not None:
            lab = np.where(in_el, labels[np.clip(which_el, 0, None)], None)
            elim_type[g.index.get_indexer(idx)] = lab

        sd, _ = signed_distance(chrom, pos, cbrs_by_chrom, elim_merged)
        signed[g.index.get_indexer(idx)] = sd
        min_dist[g.index.get_indexer(idx)] = np.abs(sd)

    g['region'] = region
    g['elim_type'] = elim_type
    g['min_cbr_dist'] = min_dist
    g['signed_cbr_dist'] = signed
    g['cbr_proximal'] = np.where(np.isnan(min_dist), False,
                                 min_dist <= cbr_window)

    counts = {r: int((g['region'] == r).sum())
              for r in ('retained', 'eliminated', 'ambiguous', 'unannotated')}
    logging.info('  region calls: %s', counts)
    logging.info('  CBR-proximal (<= %.0f kb): %s of %s with a distance',
                 cbr_window / 1000, f"{int(g['cbr_proximal'].sum()):,}",
                 f"{int(g['min_cbr_dist'].notna().sum()):,}")
    if counts['ambiguous']:
        logging.warning('  %d genes fall inside both a retained and an '
                        'eliminated interval — the two BEDs overlap',
                        counts['ambiguous'])
    return g, counts


def label_arrays(bed, merged):
    """Per-sequence arrays of the 4th-column label, aligned to `merged` intervals.

    Merging can combine intervals with different labels; the merged interval takes
    the label of the interval containing its start, and disagreements are counted
    and logged rather than hidden.
    """
    if 'label' not in bed.columns:
        return {}
    out, n_conflict = {}, 0
    for chrom, (starts, ends) in merged.items():
        block = bed[bed['chrom'] == chrom].sort_values('start')
        bs = block['start'].to_numpy(float)
        bl = block['label'].to_numpy(dtype=object)
        j = np.searchsorted(bs, starts, side='left')
        j = np.clip(j, 0, max(len(bs) - 1, 0))
        lab = bl[j] if len(bs) else np.array([], dtype=object)
        for k, (s, e) in enumerate(zip(starts, ends)):
            inside = (bs >= s) & (bs <= e)
            if inside.sum() > 1 and len(set(bl[inside])) > 1:
                n_conflict += 1
        out[chrom] = lab
    if n_conflict:
        logging.warning('  %d merged tracts combine intervals with different '
                        'labels; the label of the first is used', n_conflict)
    return out


# ---------------------------------------------------------------------------
# Statistics helpers shared across scripts
# ---------------------------------------------------------------------------

def fmt_p(p):
    """Format a p-value for a panel or a table.

    Underflow to exactly 0.0 is real here — Mann-Whitney on ~10^4 genes returns
    it — and `int(log10(0))` raises OverflowError, so it is handled explicitly
    rather than left as a latent crash on the largest comparison.
    """
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return 'NA'
    if p == 0:
        return 'P < 10$^{-300}$'
    if p >= 0.001:
        return f'P = {p:.3g}'
    e = int(np.floor(np.log10(p)))
    m = p / 10 ** e
    return f'P = {m:.1f} \\times 10$^{{{e}}}$'.replace('\\times', r'$\times$')


def stars(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return 'NA'
    for thresh, mark in ((1e-4, '****'), (1e-3, '***'), (1e-2, '**'), (0.05, '*')):
        if p < thresh:
            return mark
    return 'ns'


def benjamini_hochberg(pvals):
    """BH-adjusted q-values, order preserved. NaNs pass through as NaN."""
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    if not ok.any():
        return q
    sub = p[ok]
    order = np.argsort(sub)
    ranked = sub[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def fingerprint(path):
    if not path or not os.path.exists(path):
        return None
    st = os.stat(path)
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return {'path': os.path.abspath(path), 'bytes': st.st_size,
            'mtime_utc': datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
            'sha256_16': h.hexdigest()[:16]}


def git_commit(here=None):
    here = here or os.path.dirname(os.path.abspath(__file__))
    try:
        out = subprocess.run(['git', '-C', here, 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            dirty = subprocess.run(['git', '-C', here, 'status', '--porcelain'],
                                   capture_output=True, text=True, timeout=5)
            return out.stdout.strip() + ('+dirty' if dirty.stdout.strip() else '')
    except Exception:
        pass
    return None


def setup_logging(outdir, stem):
    """File + stdout logging. `mode='w'`: FileHandler appends by default, so a
    rerun into an existing directory concatenates logs while the figures beside
    them are overwritten, and the log stops describing the files it sits next to."""
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(os.path.join(outdir, f'{stem}.log'),
                                      mode='w'),
                  logging.StreamHandler(sys.stdout)])
