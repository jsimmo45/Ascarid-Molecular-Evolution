#!/usr/bin/env python3
"""
hic_x_dosage.py -- Do Ascaris males carry one copy of each X? Read depth of
testis versus ovary Hi-C, chromosome by chromosome.

WHY
---
The Discussion's male-biased-mutation sentence (results draft, note 24) rests on
an X spending two-thirds of its generations in females, which requires males to
carry one copy of the X and females two. Supplementary fig. S16B supports that
from X-to-autosome read depth in a published population sample, but the sex of
those worms had to be re-called from the data, because the published sexing was
wrong for about a third of them. The lab's own Hi-C libraries were made from
dissected testes and ovaries, so the sex of each is known from anatomy and the
test needs nobody else's metadata.

HOW
---
Every read end of a valid pair is a draw of DNA from where it maps, wherever its
partner maps. Read-end coverage per bin is therefore proportional to copy number
times a capture efficiency set by the bin's sequence (restriction-site density,
mappability, GC). Both libraries were made with the same Arima kit and mapped to
the same assembly by the same pipeline, so that efficiency is shared, and it
cancels in the per-bin ratio

    r_b = (male_b / male autosomal total) / (female_b / female autosomal total)

which is 1 where males and females carry the same number of copies and 0.5
where males carry one and females two.

Each chromosome is summarized by the median r_b over its usable bins and called
by the nearest expectation on a log scale: one copy in males (log2 r = -1, r =
0.5) or two (log2 r = 0, r = 1), with the boundary between them at r = 2^-0.5 =
0.707. A median within 0.25 log2 units of that boundary (0.59 < r < 0.84) is
called intermediate rather than forced into a class, and one more than 0.5 log2
units beyond either expectation is called below one copy (r < 0.35) or above two
copies (r >= 1.41). Calling EACH X is the point: S16B pools the five X
chromosomes, so it shows that the X complement is single-copy in males, not that
every X is.

Usable bins exclude
  * bins overlapping eliminated DNA (--eliminated-bed). Gonads are germline and
    keep it, but somatic carryover would deplete it in one library and not the
    other, and that would read as a dosage difference;
  * bins whose read-end density in either library is below --min-frac of that
    library's median autosomal density (assembly gaps, unmappable sequence) or
    above --max-fold times it (collapsed repeats such as rDNA, whose coverage
    does not scale with copy number);
  * sequences that are neither X nor autosome (unplaced, mitochondrial).

No downsampling. Depth enters only as one scalar per library, and a 100 kb bin
holds thousands of read ends, so Poisson noise in r_b is a percent or two.
Downsampling before ICE matters for contact maps compared pixel by pixel, not
for coverage ratios.

WHAT A TESTIS SHOULD GIVE
-------------------------
Spermatogonia and testis somatic cells carry one X set against two autosome
sets (X:A = 0.5). Spermatocytes are the same at 4C. Spermatids and sperm are
haploid, half with the X set and half without, so their pool is again 0.5. Any
mix of these gives 0.5, so the expectation does not depend on how much of the
testis is post-meiotic. If the unpaired male X is compacted by meiotic silencing
and that changes ligation efficiency, r moves off 0.5 a little; it cannot make
or remove a factor of two.

TESTS
-----
  per chromosome   median r_b with IQR, the pooled ratio over usable bins, the
                   share of deviant bins (>= 0.5 log2 units from the class
                   expectation: 0.5 for X, 1 for autosomes), and the call
  X vs autosomes   exact permutation over every way of labelling k of the N
                   chromosomes X, on the difference in mean log2 median
                   (C(24, 5) = 42,504 for Ascaris; floor 1 / C(N, k))
  segments         runs of >= --min-run consecutive deviant usable bins: on an
                   X, a stretch present twice in males (r >= 0.707) or not at
                   all (r <= 0.354); on an autosome, a stretch present once
                   (r <= 0.707) or amplified in males (r >= 1.414). Consecutive
                   means consecutive among usable bins; the reported span
                   includes any masked bins inside it
  per library      each male library against the pooled females, and, when
                   there are two or more, each female library against the other
                   females (leave-one-out), as chromosome-level pooled ratios

INPUTS
------
HiC-Pro allValidPairs (read, chr1, pos1, strand1, chr2, pos2, strand2, size,
frag1, frag2, mapq1, mapq2, ...) or 4DN .pairs (columns from its '#columns:'
header), plain or gzipped. The format is detected per file. HiC-Pro sorts and
deduplicates allValidPairs, so the whole file is read; there is no option to
read only its start, which would be one chromosome.

Chromosome lengths come from --chrom-sizes (.fai or two-column). Names must
match the pairs files; a library whose first lines carry no known name stops
the run before it is read in full.

OUTPUTS (in --outdir)
---------------------
  hic_x_dosage_summary.txt  per-chromosome calls, the tests, n at every filter,
                            and what the text can claim
  per_chromosome.tsv        the pooled comparison per chromosome, plus the
                            per-library ratios
  per_bin.tsv               read ends per library, filters and r_b per bin
  segments.tsv              discordant segments (header only if none)
  per_library.tsv           read accounting per library
  hic_x_dosage.log, provenance.json

USAGE
-----
  python hic_x_dosage.py \\
    --male   /path/to/teste.allValidPairs \\
    --female /path/to/ovary.allValidPairs \\
    --chrom-sizes $DB_DIR/ascaris/v50/AG.v50.fa.fai \\
    --eliminated-bed $DB_DIR/ascaris/v50/AG.v50.eliminated_strict.bed \\
    --workers 2 --outdir results/hic_x_dosage_<stamp>

Dependencies: numpy, pandas; region_lib from lib/.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import logging
import math
import os
import platform
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from region_lib import (fingerprint, git_commit, load_chrom_sizes,
                        merge_intervals, parse_bed, setup_logging)

_HERE = os.path.dirname(os.path.abspath(__file__))

#: log2 of the male/female dosage expected for one and for two copies in males.
LOG2_ONE, LOG2_TWO = -1.0, 0.0
#: Class boundary: the log-scale midpoint of the two, 2 ** -0.5 = 0.7071.
LOG2_MID = (LOG2_ONE + LOG2_TWO) / 2.0
BOUNDARY = 2.0 ** LOG2_MID
#: A median within this many log2 units of the boundary is 'intermediate'.
INTERMEDIATE_HALF_WIDTH = 0.25
#: A bin (or median) this many log2 units or more from its expectation is
#: deviant: half the distance between one and two copies.
DEVIANT_LOG2 = 0.5
#: Expected log2 male/female ratio per chromosome class.
EXPECTED_LOG2 = {'X': LOG2_ONE, 'A': LOG2_TWO}

#: Above this many label assignments the permutation becomes Monte Carlo.
MAX_EXACT = 2_000_000

#: Log progress roughly this often while a library is being read.
PROGRESS_EVERY = 25_000_000

PAIRS_SUFFIXES = ('.gz', '.allValidPairs', '.validPairs', '.pairs')


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def p_text(p):
    """Plain-text P for the summary."""
    if p is None or not np.isfinite(p):
        return 'NA'
    return f'{p:.2g}' if p < 0.001 else f'{p:.3f}'


def open_text(path):
    """Open plain or gzipped text."""
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def library_label(path):
    """Short label from a file name: 'AG.v50.teste.allValidPairs' -> 'AG.v50.teste'."""
    name = os.path.basename(path)
    changed = True
    while changed:
        changed = False
        for suf in PAIRS_SUFFIXES:
            if name.endswith(suf) and len(name) > len(suf):
                name = name[:-len(suf)]
                changed = True
    return name


def quick_fingerprint(path, block=4 << 20):
    """Size, mtime and a SHA-256 over the first and last 4 MiB.

    region_lib.fingerprint hashes the whole file, which for a multi-gigabyte
    allValidPairs would read it a second time. Head, tail, size and mtime
    identify the file without doubling the job's I/O.
    """
    st = os.stat(path)
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        h.update(fh.read(block))
        if st.st_size > 2 * block:
            fh.seek(-block, os.SEEK_END)
            h.update(fh.read(block))
    return {'path': os.path.abspath(path), 'bytes': st.st_size,
            'mtime_utc': datetime.fromtimestamp(st.st_mtime,
                                                timezone.utc).isoformat(),
            'sha256_head_tail_16': h.hexdigest()[:16]}


def classify_ratio(ratio):
    """Copy-number call in males (females carrying two) for a male/female ratio:
    'one copy', 'two copies', 'intermediate', 'below one copy',
    'above two copies' or 'NA'."""
    if ratio is None or not np.isfinite(ratio) or ratio <= 0:
        return 'NA'
    lr = math.log2(ratio)
    if abs(lr - LOG2_MID) < INTERMEDIATE_HALF_WIDTH:
        return 'intermediate'
    if lr < LOG2_MID:
        return 'below one copy' if lr <= LOG2_ONE - DEVIANT_LOG2 else 'one copy'
    return 'above two copies' if lr >= LOG2_TWO + DEVIANT_LOG2 else 'two copies'


def deviant(ratios, cls):
    """True where a ratio is >= DEVIANT_LOG2 log2 units from its class
    expectation (0.5 for X, 1 for autosomes)."""
    with np.errstate(divide='ignore', invalid='ignore'):
        lr = np.log2(np.asarray(ratios, dtype=float))
    return np.abs(lr - EXPECTED_LOG2[cls]) >= DEVIANT_LOG2


def chrom_sort_key(name):
    """Autosomes by number, then X chromosomes by number, then anything else."""
    m = re.match(r'^chr(X?)(\d+)$', str(name))
    if m:
        return (1 if m.group(1) else 0, int(m.group(2)), str(name))
    return (2, 0, str(name))


# ---------------------------------------------------------------------------
# Reading pairs files
# ---------------------------------------------------------------------------

def sniff_format(path, n_check=2000):
    """
    Work out the column layout of a pairs file from its first lines.

    4DN .pairs files announce their columns in a '#columns:' header line; HiC-Pro
    allValidPairs have no header and a fixed layout (chr1, pos1 in fields 2-3,
    chr2, pos2 in fields 5-6, mapq1, mapq2 in fields 11-12 when present).

    Returns
    -------
    dict with keys format, cols {chr1, pos1, chr2, pos2[, mapq1, mapq2]} as
    0-based field indices, n_header (leading '#' lines), n_fields, and
    sample_ends ({chromosome name: read ends in the sampled lines}).
    """
    n_header, columns, first = 0, None, []
    with open_text(path) as fh:
        for line in fh:
            if not first and line.startswith('#'):
                n_header += 1
                if line.startswith('#columns:'):
                    columns = line[len('#columns:'):].split()
                continue
            if not line.strip():
                continue
            first.append(line.rstrip('\n').split('\t'))
            if len(first) >= n_check:
                break
    if not first:
        raise ValueError(f'{path}: no data lines')

    n_fields = len(first[0])
    if columns:
        need = ('chrom1', 'pos1', 'chrom2', 'pos2')
        missing = [c for c in need if c not in columns]
        if missing:
            raise ValueError(f'{path}: .pairs header lacks {missing}; '
                             f'columns are {columns}')
        cols = {'chr1': columns.index('chrom1'), 'pos1': columns.index('pos1'),
                'chr2': columns.index('chrom2'), 'pos2': columns.index('pos2')}
        if 'mapq1' in columns and 'mapq2' in columns:
            cols['mapq1'] = columns.index('mapq1')
            cols['mapq2'] = columns.index('mapq2')
        fmt = 'pairs'
    else:
        if n_fields < 7:
            raise ValueError(f'{path}: {n_fields} tab-separated fields and no '
                             "'#columns:' header, so neither HiC-Pro "
                             'allValidPairs (>= 7 fields) nor 4DN .pairs')
        cols = {'chr1': 1, 'pos1': 2, 'chr2': 4, 'pos2': 5}
        if n_fields >= 12:
            cols['mapq1'], cols['mapq2'] = 10, 11
        fmt = 'allValidPairs'

    # Validate on the sampled lines: integer positions and MAPQs, and chromosome
    # fields that are not numbers (a shifted layout puts positions there).
    bad = 0
    for f in first:
        try:
            for key in ('pos1', 'pos2', 'mapq1', 'mapq2'):
                if key in cols:
                    int(f[cols[key]])
            if f[cols['chr1']].lstrip('-').isdigit():
                raise ValueError
        except (ValueError, IndexError):
            bad += 1
    if bad > 0.01 * len(first):
        raise ValueError(f'{path}: {bad} of the first {len(first)} lines do not '
                         f'parse as {fmt} (fields used: {cols})')
    end_counts = {}
    for f in first:
        for key in ('chr1', 'chr2'):
            end_counts[f[cols[key]]] = end_counts.get(f[cols[key]], 0) + 1
    return {'format': fmt, 'cols': cols, 'n_header': n_header,
            'n_fields': n_fields, 'sample_ends': end_counts}


def _init_worker():
    """Logging in a worker process. Under fork the parent's handlers are
    inherited and nothing is done. Under spawn or forkserver (the default from
    Python 3.14) a worker starts unconfigured, and it logs to stdout only: the
    parent's log file is open without O_APPEND, so lines a worker appended would
    be overwritten by the parent's next write. Everything the log file must hold
    is therefore logged by the parent from the returned stats."""
    if logging.getLogger().handlers:
        return
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)])


def count_library(path, label, chrom_names, lengths, binsize, min_mapq,
                  chunk_rows):
    """
    Read-end counts per bin for one pairs file, in one pass.

    Parameters
    ----------
    chrom_names : list of str, the sequence order that defines the bin index
    lengths     : list of int, the matching sequence lengths
    min_mapq    : an end is counted only if its own MAPQ is >= this; ignored,
                  with a note, when the file carries no MAPQ

    Returns
    -------
    (label, counts, stats) -- counts is an int64 array over all bins.
    """
    t0 = time.time()
    fmt = sniff_format(path)
    cols = fmt['cols']
    name_to_idx = {c: i for i, c in enumerate(chrom_names)}
    # Fail before reading the whole file if most read ends in its first lines
    # land on names the sizes file does not have. Requiring a majority, not just
    # one shared name, catches two assemblies that share only their unplaced
    # scaffold names.
    ends = fmt['sample_ends']
    n_known = sum(n for c, n in ends.items() if c in name_to_idx)
    if n_known < 0.5 * sum(ends.values()):
        unknown = sorted((c for c in ends if c not in name_to_idx),
                         key=lambda c: -ends[c])
        raise ValueError(
            f'{label}: only {n_known} of {sum(ends.values())} read ends in the '
            f'first lines fall on sequences named in --chrom-sizes. Unknown '
            f'names: {", ".join(unknown[:6])}; sizes file has e.g. '
            f'{", ".join(chrom_names[:4])}. Wrong assembly or sizes file?')

    lengths = np.asarray(lengths, dtype=np.int64)
    nbins = -(-lengths // binsize)                        # ceil
    offsets = np.concatenate([[0], np.cumsum(nbins)[:-1]]).astype(np.int64)
    total_bins = int(nbins.sum())
    counts = np.zeros(total_bins, dtype=np.int64)

    has_mapq = 'mapq1' in cols
    use = sorted(set(cols.values()))
    by_index = {v: k for k, v in cols.items()}
    dtype = {i: ('category' if by_index[i].startswith('chr') else 'int64')
             for i in use}
    stats = {'label': label, 'path': os.path.abspath(path),
             'format': fmt['format'], 'has_mapq': has_mapq, 'pairs': 0,
             'ends_total': 0, 'ends_unknown_chrom': 0, 'ends_low_mapq': 0,
             'ends_beyond_length': 0, 'ends_counted': 0}
    unknown_names = {}

    logging.info('[%s] %s format, MAPQ %s; reading %s', label, fmt['format'],
                 f'>= {min_mapq} per end' if has_mapq else 'not in file',
                 os.path.abspath(path))
    reader = pd.read_csv(path, sep='\t', header=None, usecols=use, dtype=dtype,
                         skiprows=fmt['n_header'], chunksize=chunk_rows,
                         engine='c', na_filter=False, compression='infer')
    next_report = PROGRESS_EVERY
    for chunk in reader:
        stats['pairs'] += len(chunk)
        for end in ('1', '2'):
            chrom_col = chunk[cols['chr' + end]]
            cats = chrom_col.cat.categories
            lut = np.fromiter((name_to_idx.get(str(c), -1) for c in cats),
                              dtype=np.int64, count=len(cats))
            codes = chrom_col.cat.codes.to_numpy()
            ci = np.where(codes >= 0, lut[np.maximum(codes, 0)], -1)
            ok = ci >= 0
            stats['ends_total'] += len(ci)
            n_unknown = int((~ok).sum())
            if n_unknown:
                stats['ends_unknown_chrom'] += n_unknown
                bad_codes = np.unique(codes[~ok])
                for bc in bad_codes[:20]:
                    nm = 'NA' if bc < 0 else str(cats[bc])
                    unknown_names[nm] = unknown_names.get(nm, 0) + int(
                        (codes == bc).sum())
            if has_mapq:
                low = chunk[cols['mapq' + end]].to_numpy() < min_mapq
                stats['ends_low_mapq'] += int((ok & low).sum())
                ok &= ~low
            pos = chunk[cols['pos' + end]].to_numpy(np.int64)[ok]
            cio = ci[ok]
            stats['ends_beyond_length'] += int((pos > lengths[cio]).sum())
            # HiC-Pro and .pairs positions are 1-based; clip 0 and anything past
            # the sequence end into the first and last bin.
            b = np.minimum(np.maximum(pos - 1, 0) // binsize, nbins[cio] - 1)
            counts += np.bincount(offsets[cio] + b, minlength=total_bins)
            stats['ends_counted'] += int(ok.sum())
        if stats['pairs'] >= next_report:
            logging.info('[%s] %s pairs read (%.0f s)', label,
                         f"{stats['pairs']:,}", time.time() - t0)
            next_report += PROGRESS_EVERY

    stats['seconds'] = round(time.time() - t0, 1)
    stats['unknown_names'] = dict(sorted(unknown_names.items(),
                                         key=lambda kv: -kv[1])[:10])
    if stats['ends_counted'] == 0:
        raise ValueError(f'{label}: no read end was counted')
    if stats['ends_unknown_chrom'] > 0.5 * stats['ends_total']:
        raise ValueError(
            f"{label}: {stats['ends_unknown_chrom']:,} of {stats['ends_total']:,} "
            f'read ends fall on sequences missing from --chrom-sizes '
            f"(most common: {stats['unknown_names']}). Wrong assembly?")
    return label, counts, stats


# ---------------------------------------------------------------------------
# Bins, filters and ratios
# ---------------------------------------------------------------------------

def build_bins(chrom_names, lengths, binsize, x_pattern, autosome_pattern):
    """One row per bin, in the order count_library uses: chrom, 0-based start,
    end, length, class ('X', 'A' or 'other')."""
    xr, ar = re.compile(x_pattern), re.compile(autosome_pattern)
    rows = []
    for c, L in zip(chrom_names, lengths):
        cls = 'X' if xr.search(c) else 'A' if ar.search(c) else 'other'
        n = -(-int(L) // binsize)
        starts = np.arange(n, dtype=np.int64) * binsize
        ends = np.minimum(starts + binsize, int(L))
        rows.append(pd.DataFrame({'chrom': c, 'start': starts, 'end': ends,
                                  'length': ends - starts, 'class': cls}))
    return pd.concat(rows, ignore_index=True)


def flag_eliminated(bins, elim_merged):
    """True where a bin overlaps an eliminated interval.

    elim_merged is region_lib.merge_intervals output: {chrom: (starts, ends)},
    1-based inclusive. A bin [start, end) 0-based is [start + 1, end] 1-based.
    """
    flag = np.zeros(len(bins), dtype=bool)
    for chrom, (s, e) in elim_merged.items():
        sel = np.flatnonzero(bins['chrom'].to_numpy() == chrom)
        if not len(sel) or not len(s):
            continue
        bs = bins['start'].to_numpy()[sel] + 1
        be = bins['end'].to_numpy()[sel]
        k = np.searchsorted(e, bs, side='left')      # first interval ending >= bs
        hit = k < len(s)
        hit[hit] = s[k[hit]] <= be[hit]
        flag[sel] = hit
    return flag


def depth_flags(counts, bins, base_ok, min_frac, max_fold):
    """Low- and high-density flags for one library, against its own median
    autosomal density over `base_ok` bins. Density is per bp, so a short final
    bin of a chromosome is judged fairly."""
    dens = counts / bins['length'].to_numpy()
    ref = base_ok & (bins['class'].to_numpy() == 'A')
    med = float(np.median(dens[ref])) if ref.any() else np.nan
    if not np.isfinite(med) or med <= 0:
        raise ValueError('median autosomal read-end density is zero; is the '
                         'library empty, or are the autosome names wrong?')
    return dens < min_frac * med, dens > max_fold * med, med


def ratio_per_bin(male, female, usable, is_a):
    """r_b over usable bins (NaN elsewhere), each library scaled by its own
    usable autosomal total. Returns (ratio, male_total, female_total)."""
    m_tot = float(male[usable & is_a].sum())
    f_tot = float(female[usable & is_a].sum())
    r = np.full(len(male), np.nan)
    ok = usable & (female > 0)
    r[ok] = (male[ok] / m_tot) / (female[ok] / f_tot)
    return r, m_tot, f_tot


def pooled_ratio(male, female, sel, usable_a):
    """Chromosome-level ratio of summed read ends, each library scaled by its
    usable autosomal total."""
    fa = female[usable_a].sum()
    ma = male[usable_a].sum()
    fs = female[sel].sum()
    if fa == 0 or ma == 0 or fs == 0:
        return np.nan
    return float((male[sel].sum() / ma) / (fs / fa))


def chromosome_table(bins, r, usable, male, female, masks):
    """Per-chromosome summary of the pooled comparison."""
    is_a = bins['class'].to_numpy() == 'A'
    usable_a = usable & is_a
    rows = []
    chroms = [c for c in bins['chrom'].unique()
              if bins.loc[bins['chrom'] == c, 'class'].iloc[0] in ('A', 'X')]
    for chrom in sorted(chroms, key=chrom_sort_key):
        on = bins['chrom'].to_numpy() == chrom
        cls = bins.loc[on, 'class'].iloc[0]
        sel = on & usable
        rv = r[sel]
        rv = rv[np.isfinite(rv)]
        med = float(np.median(rv)) if len(rv) else np.nan
        rows.append({
            'chrom': chrom, 'class': cls,
            'length_mb': round(float(bins.loc[on, 'length'].sum()) / 1e6, 3),
            'n_bins': int(on.sum()), 'n_usable': int(sel.sum()),
            'n_eliminated': int((on & masks['eliminated']).sum()),
            'n_low': int((on & masks['low']).sum()),
            'n_high': int((on & masks['high']).sum()),
            'median_ratio': med,
            'q25': float(np.percentile(rv, 25)) if len(rv) else np.nan,
            'q75': float(np.percentile(rv, 75)) if len(rv) else np.nan,
            'pooled_ratio': pooled_ratio(male, female, sel, usable_a),
            'frac_deviant': float(deviant(rv, cls).mean()) if len(rv) else np.nan,
            'call': classify_ratio(med),
            'male_copies': 2.0 * med if np.isfinite(med) else np.nan})
    return pd.DataFrame(rows)


def chromosome_permutation(values, is_x, seed=0, n_mc=200_000):
    """Exact (or Monte Carlo) permutation test on per-chromosome values.

    Statistic: mean of the chromosomes labelled X minus mean of the rest. Every
    way of choosing which k of the N chromosomes are X is enumerated, so P is
    exact. Returns dict: stat, p_two, p_x_lower, n_assignments, exact.
    """
    m = np.asarray(values, dtype=float)
    is_x = np.asarray(is_x, dtype=bool)
    n, k = len(m), int(is_x.sum())
    res = {'stat': np.nan, 'p_two': np.nan, 'p_x_lower': np.nan,
           'n_assignments': 0, 'exact': True, 'n_x': k, 'n_a': n - k}
    if k < 1 or n - k < 1:
        return res
    obs = m[is_x].mean() - m[~is_x].mean()
    n_assign = math.comb(n, k)
    if n_assign <= MAX_EXACT:
        combos = np.fromiter(itertools.chain.from_iterable(
            itertools.combinations(range(n), k)), dtype=np.int64,
            count=n_assign * k).reshape(n_assign, k)
        sel = m[combos].sum(axis=1)
    else:
        rng = np.random.default_rng(seed)
        sel = np.array([m[rng.choice(n, k, replace=False)].sum()
                        for _ in range(n_mc)])
        res['exact'] = False
        n_assign = n_mc
    stats = sel / k - (m.sum() - sel) / (n - k)
    tol = 1e-12 * max(1.0, abs(obs))
    res.update(stat=float(obs), n_assignments=int(n_assign),
               p_two=float(np.mean(np.abs(stats) >= abs(obs) - tol)),
               p_x_lower=float(np.mean(stats <= obs + tol)))
    return res


def discordant_segments(bins, r, usable, min_run):
    """Runs of >= min_run consecutive usable bins that are deviant for their
    chromosome's class (see `deviant`)."""
    out = []
    cols = ['chrom', 'class', 'start', 'end', 'n_bins', 'n_masked_inside',
            'median_ratio']
    use_idx = np.flatnonzero(usable & np.isfinite(r))
    if not len(use_idx):
        return pd.DataFrame(columns=cols)
    sub = bins.iloc[use_idx]
    for chrom, block in sub.groupby('chrom', sort=False):
        cls = block['class'].iloc[0]
        idx = block.index.to_numpy()
        rv = r[idx]
        wrong = deviant(rv, cls)
        padded = np.concatenate([[False], wrong, [False]]).astype(np.int8)
        d = np.diff(padded)
        for s, e in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)):
            if e - s < min_run:
                continue
            first, last = idx[s], idx[e - 1]
            out.append({'chrom': chrom, 'class': cls,
                        'start': int(bins.at[first, 'start']),
                        'end': int(bins.at[last, 'end']),
                        'n_bins': int(e - s),
                        'n_masked_inside': int((last - first + 1) - (e - s)),
                        'median_ratio': float(np.median(rv[s:e]))})
    return pd.DataFrame(out, columns=cols)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def class_line(tab, cls, label):
    """One summary line for X or autosomes."""
    t = tab[tab['class'] == cls]
    if not len(t):
        return f'{label:<15} none'
    calls = t['call'].value_counts().to_dict()
    want = 'one copy' if cls == 'X' else 'two copies'
    return (f"{label:<15} per-chromosome median ratio "
            f"{t['median_ratio'].min():.3f}-{t['median_ratio'].max():.3f}; "
            f"{calls.get(want, 0)} of {len(t)} called {want}"
            + (f" ({', '.join(f'{v} {k}' for k, v in calls.items() if k != want)})"
               if len(calls) > 1 or want not in calls else ''))


def claim_text(tab, perm):
    """The sentence the data support, computed rather than assumed."""
    x = tab[tab['class'] == 'X']
    a = tab[tab['class'] == 'A']
    if not len(x) or not len(a):
        return 'No claim: X chromosomes or autosomes missing from the comparison.'
    x_ok = bool((x['call'] == 'one copy').all())
    a_ok = bool((a['call'] == 'two copies').all())
    # Three decimals: at two, five X medians of 0.447-0.454 printed as
    # '0.45-0.45' (first real run, 2026-10-11).
    xr = f"{x['median_ratio'].min():.3f}-{x['median_ratio'].max():.3f}"
    ar = f"{a['median_ratio'].min():.3f}-{a['median_ratio'].max():.3f}"
    if x_ok and a_ok:
        return (f'Testis Hi-C carries each of the {len(x)} X chromosomes at '
                f'{xr} times its ovary read depth, against {ar} for the '
                f'{len(a)} autosomes (exact permutation over chromosomes, '
                f'P = {p_text(perm["p_two"])}). Males carry a single copy of '
                'every X.')
    bad_x = x.loc[x['call'] != 'one copy', ['chrom', 'median_ratio', 'call']]
    bad_a = a.loc[a['call'] != 'two copies', ['chrom', 'median_ratio', 'call']]
    parts = []
    if len(bad_x):
        parts.append('X chromosomes not called one copy: ' + '; '.join(
            f"{r.chrom} {r.median_ratio:.2f} ({r.call})"
            for r in bad_x.itertuples()))
    if len(bad_a):
        parts.append('autosomes not called two copies: ' + '; '.join(
            f"{r.chrom} {r.median_ratio:.2f} ({r.call})"
            for r in bad_a.itertuples()))
    return ('Do NOT state that males carry one copy of every X from these data. '
            + ' | '.join(parts) + '. Inspect per_bin.tsv and segments.tsv '
            'before interpreting.')


def write_summary(path, args, lib_stats, males, females, bins, masks, usable,
                  tab, perm, segs, r, lib_tab, medians):
    """Human-readable results, ending with what the text can claim."""
    is_x = bins['class'].to_numpy() == 'X'
    is_a = bins['class'].to_numpy() == 'A'
    L = []
    L.append('Ascaris X dosage from Hi-C read depth: testis (male) vs ovary (female)')
    L.append('=' * 74)
    L.append('')
    for sex, labs in (('Male', males), ('Female', females)):
        for lab in labs:
            s = lib_stats[lab]
            L.append(f"{sex + ' library':<16} {lab}: {s['path']}")
    L.append(f"{'Bins':<16} {args.binsize:,} bp; ratio = (male / male usable "
             'autosomal total) / (female / female usable autosomal total)')
    L.append(f"{'MAPQ':<16} >= {args.min_mapq} per read end"
             + ('' if all(lib_stats[k]['has_mapq'] for k in lib_stats)
                else ' (not applied to libraries without MAPQ columns)'))
    L.append('')

    L.append('Read accounting (n at every stage)')
    L.append('-' * 34)
    hdr = (f"  {'library':<24} {'sex':<6} {'pairs':>14} {'ends':>14} "
           f"{'unknown chr':>12} {'low MAPQ':>12} {'counted':>14} "
           f"{'usable A':>13} {'usable X':>12} {'raw X:A':>8}")
    L.append(hdr)
    for rec in lib_tab.itertuples():
        L.append(f"  {rec.library:<24} {rec.sex:<6} {rec.pairs:>14,} "
                 f"{rec.ends_total:>14,} {rec.ends_unknown_chrom:>12,} "
                 f"{rec.ends_low_mapq:>12,} {rec.ends_counted:>14,} "
                 f"{rec.ends_usable_A:>13,} {rec.ends_usable_X:>12,} "
                 f"{rec.raw_x_to_a:>8.3f}")
    L.append('  raw X:A = read ends per bp on usable X bins / on usable autosomal '
             'bins; it includes capture bias, which the male/female ratio cancels.')
    L.append('')

    n_cls = {c: int((bins['class'] == c).sum()) for c in ('A', 'X', 'other')}
    L.append('Bins (n at every filter)')
    L.append('-' * 24)
    L.append(f"  total {len(bins):,}: autosomal {n_cls['A']:,}, X {n_cls['X']:,}, "
             f"other sequences {n_cls['other']:,} (excluded)")
    ax = is_a | is_x
    L.append(f"  overlapping eliminated DNA: {int((ax & masks['eliminated']).sum()):,}"
             + ('' if args.eliminated_bed else ' (no --eliminated-bed given)'))
    L.append(f"  below {args.min_frac:g} x median density in either library: "
             f"{int((ax & masks['low'] & ~masks['eliminated']).sum()):,}")
    L.append(f"  above {args.max_fold:g} x median density in either library: "
             f"{int((ax & masks['high'] & ~masks['eliminated'] & ~masks['low']).sum()):,}")
    L.append(f"  usable: autosomal {int((usable & is_a).sum()):,}, "
             f"X {int((usable & is_x).sum()):,}")
    L.append('  median autosomal read-end density (per kb): ' + ', '.join(
        f'{k} {v * 1000:.1f}' for k, v in medians.items()))
    L.append('')

    L.append('Per chromosome: pooled male / pooled female, median of per-bin ratios')
    L.append('-' * 70)
    L.append(f"  {'chrom':<8} {'class':<5} {'Mb':>7} {'bins':>5} {'usable':>6} "
             f"{'median':>7} {'IQR':>13} {'pooled':>7} {'deviant':>8}  call")
    for d in tab.to_dict('records'):
        L.append(f"  {d['chrom']:<8} {d['class']:<5} {d['length_mb']:>7.2f} "
                 f"{d['n_bins']:>5} {d['n_usable']:>6} {d['median_ratio']:>7.3f} "
                 f"{d['q25']:>6.3f}-{d['q75']:<6.3f} {d['pooled_ratio']:>7.3f} "
                 f"{100 * d['frac_deviant']:>7.1f}%  {d['call']}")
    L.append(f"  Calls: nearest of one copy (ratio 0.5) and two copies (1) on a "
             f"log scale, boundary {BOUNDARY:.3f}; 'intermediate' when the median "
             f"is within {INTERMEDIATE_HALF_WIDTH:g} log2 units of it "
             f"({2 ** (LOG2_MID - INTERMEDIATE_HALF_WIDTH):.3f}-"
             f"{2 ** (LOG2_MID + INTERMEDIATE_HALF_WIDTH):.3f}); 'below one copy' "
             f"< {2 ** (LOG2_ONE - DEVIANT_LOG2):.3f}, 'above two copies' >= "
             f"{2 ** (LOG2_TWO + DEVIANT_LOG2):.3f}. 'deviant' = share of usable "
             f"bins >= {DEVIANT_LOG2:g} log2 units from their class expectation.")
    L.append('')

    L.append(class_line(tab, 'X', 'X chromosomes'))
    L.append(class_line(tab, 'A', 'Autosomes'))
    rx = r[usable & is_x]
    ra = r[usable & is_a]
    rx, ra = rx[np.isfinite(rx)], ra[np.isfinite(ra)]
    if len(rx) and len(ra):
        L.append(f"{'Per-bin ratio':<15} X bins median {np.median(rx):.3f} "
                 f"(5th-95th {np.percentile(rx, 5):.3f}-{np.percentile(rx, 95):.3f}, "
                 f"n = {len(rx):,}); autosomal bins median {np.median(ra):.3f} "
                 f"(5th-95th {np.percentile(ra, 5):.3f}-{np.percentile(ra, 95):.3f}, "
                 f"n = {len(ra):,}). The autosomal spread is the noise floor.")
    L.append(f"{'X vs autosomes':<15} exact permutation on per-chromosome log2 "
             f"median ratios ({perm['n_assignments']:,} "
             f"{'exact' if perm['exact'] else 'Monte Carlo'} assignments of "
             f"{perm['n_x']} X among {perm['n_x'] + perm['n_a']}): "
             f"P = {p_text(perm['p_two'])} two-sided, "
             f"{p_text(perm['p_x_lower'])} for X lower")
    L.append('')

    L.append(f'Discordant segments (>= {args.min_run} consecutive usable bins '
             f'>= {DEVIANT_LOG2:g} log2 units from their class expectation)')
    L.append('-' * 70)
    if len(segs):
        for s in segs.to_dict('records'):
            L.append(f"  {s['chrom']}:{s['start']:,}-{s['end']:,} ({s['class']}) "
                     f"{s['n_bins']} bins, median ratio {s['median_ratio']:.3f}"
                     + (f", {s['n_masked_inside']} masked bins inside"
                        if s['n_masked_inside'] else ''))
    else:
        L.append('  none')
    L.append('')

    lib_cols = [c for c in tab.columns if c.startswith('ratio_')]
    if lib_cols:
        L.append('Per library, chromosome-level pooled ratio against the female '
                 'pool (leave-one-out for female libraries)')
        L.append('-' * 70)
        width = max(10, max(len(c) - 6 for c in lib_cols) + 2)
        L.append('  ' + f"{'chrom':<8}" + ''.join(f'{c[6:]:>{width}}'
                                                  for c in lib_cols))
        for d in tab.to_dict('records'):
            L.append('  ' + f"{d['chrom']:<8}" + ''.join(
                f'{d[c]:>{width}.3f}' for c in lib_cols))
        L.append('')

    L.append('What this can claim')
    L.append('-' * 19)
    L.append('  ' + claim_text(tab, perm))
    L.append('')
    L.append('Notes for Methods')
    L.append('-' * 17)
    L.append(f'  Read-end coverage of HiC-Pro valid pairs (MAPQ >= {args.min_mapq} '
             f'per end) in {args.binsize // 1000} kb bins, testis against ovary. '
             'Bins overlapping eliminated DNA, bins below '
             f'{args.min_frac:g} or above {args.max_fold:g} times the median '
             'autosomal density in either library, and sequences other than '
             'X chromosomes and autosomes were excluded. Each library was scaled '
             'by its usable autosomal total; the per-bin ratio cancels capture '
             'bias shared by the two libraries (same Arima kit, assembly and '
             'pipeline). Chromosomes were called one copy or two copies in males '
             'by the nearer expectation (0.5 or 1) on a log scale.')
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Male/female X dosage from testis vs ovary Hi-C read depth.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument('--male', nargs='+', required=True,
                   help='Pairs file(s) from testis (HiC-Pro allValidPairs or '
                        '.pairs, optionally .gz). Several are pooled.')
    p.add_argument('--female', nargs='+', required=True,
                   help='Pairs file(s) from ovary. Several are pooled.')
    p.add_argument('--chrom-sizes', required=True,
                   help='.fai or two-column sequence lengths of the assembly the '
                        'pairs were mapped to')
    p.add_argument('--eliminated-bed', default=None,
                   help='Eliminated intervals (BED); overlapping bins are '
                        'excluded. Strongly recommended.')
    p.add_argument('--binsize', type=int, default=100_000,
                   help='Bin size in bp (default: 100000)')
    p.add_argument('--min-mapq', type=int, default=30,
                   help='Per-end MAPQ floor (default: 30)')
    p.add_argument('--min-frac', type=float, default=0.2,
                   help='Drop bins below this fraction of the median autosomal '
                        'density in either library (default: 0.2)')
    p.add_argument('--max-fold', type=float, default=5.0,
                   help='Drop bins above this multiple of the median autosomal '
                        'density in either library (default: 5)')
    p.add_argument('--min-run', type=int, default=5,
                   help='Consecutive usable bins for a discordant segment '
                        '(default: 5)')
    p.add_argument('--x-pattern', default=r'^chrX',
                   help=r"Regex for X chromosomes (default: '^chrX')")
    p.add_argument('--autosome-pattern', default=r'^chr\d+$',
                   help=r"Regex for autosomes (default: '^chr\d+$')")
    p.add_argument('--workers', type=int, default=1,
                   help='Parallel readers, capped at the number of libraries '
                        '(default: 1)')
    p.add_argument('--chunk-rows', type=int, default=2_000_000,
                   help='Rows per read chunk (default: 2000000)')
    p.add_argument('--seed', type=int, default=0,
                   help='Seed for the Monte Carlo permutation, used only when '
                        'exact enumeration is too large (default: 0)')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()

    if args.binsize < 1000:
        p.error('--binsize must be >= 1000')
    if not 0 <= args.min_frac < 1 or args.max_fold <= 1:
        p.error('--min-frac must be in [0, 1) and --max-fold > 1')
    if args.min_run < 1 or args.workers < 1 or args.chunk_rows < 1000:
        p.error('--min-run and --workers must be >= 1, --chunk-rows >= 1000')
    paths = [os.path.abspath(x) for x in args.male + args.female]
    if len(set(paths)) != len(paths):
        p.error('the same file is given more than once')
    for x in args.male + args.female + [args.chrom_sizes]:
        if not os.path.isfile(x):
            p.error(f'not found: {x}')
    if args.eliminated_bed and not os.path.isfile(args.eliminated_bed):
        p.error(f'not found: {args.eliminated_bed}')
    if args.outdir is None:
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.outdir = os.path.join('results', f'hic_x_dosage_{stamp}')
    return args


def unique_labels(paths):
    """Labels from file names, made unique if two files share a stem."""
    labs, seen = [], {}
    for path in paths:
        lab = library_label(path)
        if lab in seen:
            seen[lab] += 1
            lab = f'{lab}_{seen[lab]}'
        else:
            seen[lab] = 1
        labs.append(lab)
    return labs


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'hic_x_dosage')
    logging.info('X dosage from testis vs ovary Hi-C read depth')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))

    # --- assembly and bins ---------------------------------------------------
    sizes, _ = load_chrom_sizes(args.chrom_sizes)
    chrom_names = list(sizes)
    lengths = [sizes[c] for c in chrom_names]
    bins = build_bins(chrom_names, lengths, args.binsize, args.x_pattern,
                      args.autosome_pattern)
    x_names = sorted({c for c in chrom_names if re.search(args.x_pattern, c)},
                     key=chrom_sort_key)
    a_names = [c for c in chrom_names if re.search(args.autosome_pattern, c)
               and c not in x_names]
    logging.info('  %d sequences: %d autosomes, %d X (%s), %d other; %s bins',
                 len(chrom_names), len(a_names), len(x_names),
                 ', '.join(x_names) or 'none',
                 len(chrom_names) - len(a_names) - len(x_names), f'{len(bins):,}')
    if not x_names or not a_names:
        sys.exit(f'ERROR: need X chromosomes and autosomes; --x-pattern '
                 f'{args.x_pattern!r} matched {len(x_names)}, --autosome-pattern '
                 f'{args.autosome_pattern!r} matched {len(a_names)}. Names: '
                 f'{chrom_names[:8]}')

    masks = {'eliminated': np.zeros(len(bins), dtype=bool)}
    if args.eliminated_bed:
        bed = parse_bed(args.eliminated_bed)
        shared = set(bed['chrom']) & set(chrom_names)
        logging.info('  eliminated BED: %d of %d sequences named in the sizes file',
                     len(shared), bed['chrom'].nunique())
        if not shared:
            sys.exit('ERROR: no sequence in --eliminated-bed is in --chrom-sizes; '
                     'they are from different assemblies.')
        masks['eliminated'] = flag_eliminated(bins, merge_intervals(bed))
        logging.info('  bins overlapping eliminated DNA: %s',
                     f"{int(masks['eliminated'].sum()):,}")
    else:
        logging.warning('  no --eliminated-bed: eliminated bins stay in, so '
                        'somatic carryover in either gonad could shift ratios.')

    # --- read the libraries --------------------------------------------------
    male_labels = unique_labels(args.male)
    female_labels = unique_labels(args.female)
    clash = set(male_labels) & set(female_labels)
    if clash:
        female_labels = [f'{x}_female' if x in clash else x for x in female_labels]
    jobs = ([(p, l, 'male') for p, l in zip(args.male, male_labels)]
            + [(p, l, 'female') for p, l in zip(args.female, female_labels)])
    workers = min(args.workers, len(jobs))
    logging.info('Reading %d libraries with %d worker(s)', len(jobs), workers)

    counts, lib_stats = {}, {}

    def collect(lab, c, s):
        counts[lab], lib_stats[lab] = c, s
        logging.info('[%s] done: %s pairs, %s read ends counted (%s on unknown '
                     'sequences, %s below MAPQ) in %.0f s', lab,
                     f"{s['pairs']:,}", f"{s['ends_counted']:,}",
                     f"{s['ends_unknown_chrom']:,}", f"{s['ends_low_mapq']:,}",
                     s['seconds'])

    if workers == 1:
        for path, lab, _ in jobs:
            collect(*count_library(path, lab, chrom_names, lengths,
                                   args.binsize, args.min_mapq, args.chunk_rows))
    else:
        with ProcessPoolExecutor(max_workers=workers,
                                 initializer=_init_worker) as ex:
            futs = [ex.submit(count_library, path, lab, chrom_names, lengths,
                              args.binsize, args.min_mapq, args.chunk_rows)
                    for path, lab, _ in jobs]
            for fut in as_completed(futs):
                collect(*fut.result())
    sex_of = {lab: sex for _, lab, sex in jobs}

    # --- filters ---------------------------------------------------------------
    male = np.sum([counts[l] for l in male_labels], axis=0).astype(np.int64)
    female = np.sum([counts[l] for l in female_labels], axis=0).astype(np.int64)
    is_a = bins['class'].to_numpy() == 'A'
    is_x = bins['class'].to_numpy() == 'X'
    base_ok = (is_a | is_x) & ~masks['eliminated']
    m_low, m_high, m_med = depth_flags(male, bins, base_ok, args.min_frac,
                                       args.max_fold)
    f_low, f_high, f_med = depth_flags(female, bins, base_ok, args.min_frac,
                                       args.max_fold)
    masks['low'] = m_low | f_low
    masks['high'] = m_high | f_high
    usable = base_ok & ~masks['low'] & ~masks['high']
    medians = {'male pool': m_med, 'female pool': f_med}
    logging.info('Usable bins: %s autosomal, %s X (low %s, high %s)',
                 f'{int((usable & is_a).sum()):,}', f'{int((usable & is_x).sum()):,}',
                 f"{int((base_ok & masks['low']).sum()):,}",
                 f"{int((base_ok & masks['high'] & ~masks['low']).sum()):,}")
    if not (usable & is_x).any():
        sys.exit('ERROR: no usable X bin. Check names and filters in the log.')

    # --- the comparison --------------------------------------------------------
    r, m_tot, f_tot = ratio_per_bin(male, female, usable, is_a)
    tab = chromosome_table(bins, r, usable, male, female, masks)

    # Per-library chromosome-level ratios against the female pool (male
    # libraries) or the other females (female libraries, when there are two+).
    usable_a = usable & is_a
    for lab in male_labels + female_labels:
        if sex_of[lab] == 'male':
            ref = female
        else:
            others = [l for l in female_labels if l != lab]
            if not others:
                continue
            ref = np.sum([counts[l] for l in others], axis=0)
        vals = []
        for chrom in tab['chrom']:
            sel = (bins['chrom'].to_numpy() == chrom) & usable
            vals.append(pooled_ratio(counts[lab], ref, sel, usable_a))
        tab[f'ratio_{lab}'] = vals

    # Chromosomes with no usable bin have no median and drop out of the test.
    med = tab['median_ratio'].to_numpy(float)
    has = np.isfinite(med) & (med > 0)
    if (~has).any():
        logging.warning('  no usable bins, left out of the permutation: %s',
                        ', '.join(tab.loc[~has, 'chrom']))
    perm = chromosome_permutation(np.log2(med[has]),
                                  (tab['class'] == 'X').to_numpy()[has],
                                  seed=args.seed)
    segs = discordant_segments(bins, r, usable, args.min_run)
    for d in tab.to_dict('records'):
        logging.info('  %-7s %s median %.3f (IQR %.3f-%.3f), pooled %.3f -> %s',
                     d['chrom'], d['class'], d['median_ratio'], d['q25'],
                     d['q75'], d['pooled_ratio'], d['call'])
    logging.info('X vs autosomes: exact permutation P = %s; %d discordant '
                 'segment(s)', p_text(perm['p_two']), len(segs))

    # --- tables ------------------------------------------------------------------
    lib_rows = []
    for lab in male_labels + female_labels:
        s, c = lib_stats[lab], counts[lab]
        ua, ux = int(c[usable & is_a].sum()), int(c[usable & is_x].sum())
        bp_a = float(bins.loc[usable & is_a, 'length'].sum())
        bp_x = float(bins.loc[usable & is_x, 'length'].sum())
        lib_rows.append({'library': lab, 'sex': sex_of[lab], 'path': s['path'],
                         'format': s['format'], 'has_mapq': s['has_mapq'],
                         'pairs': s['pairs'], 'ends_total': s['ends_total'],
                         'ends_unknown_chrom': s['ends_unknown_chrom'],
                         'ends_low_mapq': s['ends_low_mapq'],
                         'ends_beyond_length': s['ends_beyond_length'],
                         'ends_counted': s['ends_counted'],
                         'ends_usable_A': ua, 'ends_usable_X': ux,
                         'raw_x_to_a': (ux / bp_x) / (ua / bp_a)
                         if ua and bp_x and bp_a else np.nan,
                         'unknown_names': json.dumps(s['unknown_names']),
                         'seconds': s['seconds']})
    lib_tab = pd.DataFrame(lib_rows)

    per_bin = bins.copy()
    for lab in male_labels + female_labels:
        per_bin[f'ends_{lab}'] = counts[lab]
    per_bin['male_pool'] = male
    per_bin['female_pool'] = female
    per_bin['eliminated'] = masks['eliminated']
    per_bin['low'] = masks['low']
    per_bin['high'] = masks['high']
    per_bin['usable'] = usable
    per_bin['ratio'] = r

    out = args.outdir
    per_bin.to_csv(os.path.join(out, 'per_bin.tsv'), sep='\t', index=False,
                   float_format='%.6g')
    tab.to_csv(os.path.join(out, 'per_chromosome.tsv'), sep='\t', index=False,
               float_format='%.6g')
    segs.to_csv(os.path.join(out, 'segments.tsv'), sep='\t', index=False,
                float_format='%.6g')
    lib_tab.to_csv(os.path.join(out, 'per_library.tsv'), sep='\t', index=False,
                   float_format='%.6g')
    summary = os.path.join(out, 'hic_x_dosage_summary.txt')
    write_summary(summary, args, lib_stats, male_labels, female_labels, bins,
                  masks, usable, tab, perm, segs, r, lib_tab, medians)

    with open(os.path.join(out, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(),
                   'python': platform.python_version(),
                   'pandas': pd.__version__, 'numpy': np.__version__,
                   'argv': sys.argv,
                   'libraries': {lab: {'sex': sex_of[lab],
                                       **quick_fingerprint(lib_stats[lab]['path'])}
                                 for lab in male_labels + female_labels},
                   'inputs': {'chrom_sizes': fingerprint(args.chrom_sizes),
                              'eliminated_bed': fingerprint(args.eliminated_bed)},
                   'n_bins': int(len(bins)),
                   'n_usable': {'A': int((usable & is_a).sum()),
                                'X': int((usable & is_x).sum())}},
                  fh, indent=2)
    logging.info('Summary: %s', summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
