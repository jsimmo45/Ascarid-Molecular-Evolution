#!/usr/bin/env python3
"""
cbr_edge_check.py - does the retained bin next to the CBR in fig. 3D differ
from the rest of the retained arm?

Jianbin Wang's question (legend review 2026-10-26, meeting 2026-10-28): in
fig. 3D (Figure 4D before the 2026-10-26 renumbering) the retained bin closest
to the CBR has a mean log10(dN/dS) near the eliminated level. Does it mean
anything?

The figure's own tests do not answer that. The retained slope is fitted over
the whole window (+/- 150 kb), and supplementary fig. S7B compares genes within
100 kb of a CBR with interior genes, so an effect confined to the first
12.5 kb would barely move either. This script tests the narrow edge directly,
and lists the genes in it.

Same genes as the figure
------------------------
The table, filters and distance are those of cbr_changepoint_relax.py:
log10(dN/dS) for genes with --omega-min <= dN/dS <= --omega-max (its
omega_filter_counts() is imported and must account for every gene), the signed
distance to the nearest CBR in kb (negative = retained side), |d| <=
--slope-window. The bins are the figure's: --n-bins over +/- --slope-window,
each gene in the bin [lo, hi) holding its midpoint; the figure draws a bin's
mean when it holds at least --min-bin genes. The arm fits are the figure's
(fit_arms, imported). With --stats the fits are compared with the figure run's
cbr_changepoint_stats.tsv and the script stops if they differ, so the bins
tested here are the bins drawn there.

Tests
-----
For each arm separately (retained d < 0, eliminated d > 0) and each edge width
w in --edge-widths (the drawn bin width is always included), the genes within
w of the CBR are compared with the rest of that arm inside the window:

  mean log10(dN/dS) of edge and rest, their difference, Welch's t-test P
  (the figure plots means of log10 dN/dS, so this is the test of what is drawn)
  median dN/dS of edge and rest, their ratio, with a 95% percentile bootstrap
  interval (--n-boot resamples, each group resampled separately)
  Mann-Whitney U P (two-sided; the same on dN/dS and on its log)
  the mean residual of the edge genes from the arm's own fitted line, with a
  one-sample t-test P (is the drawn point off the line?)

Every test is run on all genes, and again on the genes whose span stays on
their own side of the boundary ("clear"): a gene is called by its midpoint, so
a retained gene can still reach across the CBR into eliminated sequence, and
part of its coding sequence is then eliminated. Genes overlapping any
--flag-bed are also left out of "clear". BH q across every row of the table.

The genes in each arm's edge bin are listed with their coordinates and rates,
the nearest CBR, the bp of their span on the far side of it, the gap from
their nearer end to it, and their overlap with each --flag-bed.

Genes beyond the fig. 3D y-axis (added 2026-10-30)
---------------------------------------------------
plot_slope_step() caps the y-axis at the (100 - cap)th and cap-th percentiles
of the drawn genes' log10(dN/dS), each widened by 16% of the range between
them, and matplotlib clips any point beyond: it is not drawn, although every
fit, bin mean and test uses it. The figure run does not count these genes, so
this script does: the same rule on the same genes (figure_ylim() below; cap
from --cap-pct, else the figure run's provenance.json beside --stats, else
the figure script's default of 99). They are listed in
cbr_edge_offaxis_genes.tsv, and "What to say" gives the sentence for the 3D
legend.

Writes, in --outdir
-------------------
  cbr_edge_bins.tsv        every bin: arm, edges, n, mean and SE of log10
                           dN/dS, median dN/dS, whether the figure draws it
  cbr_edge_tests.tsv       every test above, one row per arm x width x gene set
  cbr_edge_genes.tsv       the genes of the two edge bins
  cbr_edge_offaxis_genes.tsv  the genes beyond the fig. 3D y-axis (header only
                           when there are none)
  cbr_edge_summary.txt     the same, for reading, ending with what to say
  cbr_edge_check.log       the run log
  provenance.json          inputs (with checksums), arguments, versions

Usage
-----
    python cbr_edge_check.py \\
      --annotated kaks_regions/kaks_annotated.tsv \\
      --eliminated-bed AG.v50.eliminated_strict.bed \\
      --chrom-sizes AG.v50.fa.fai \\
      --stats results/cbr_changepoint_20261027/cbr_changepoint_stats.tsv \\
      --outdir results/cbr_edge_check_20261029

Dependencies: numpy, pandas, scipy, region_lib, cbr_changepoint_relax (both
beside this script or on PYTHONPATH; the latter needs paper_colors from lib/).
No figure.
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
import scipy
from scipy.stats import mannwhitneyu, ttest_1samp, ttest_ind

from region_lib import (parse_bed, merge_intervals, load_chrom_sizes,
                        derive_cbrs, signed_distance, nearest_cbr,
                        benjamini_hochberg, fingerprint, git_commit,
                        setup_logging)
# The figure's own definitions: one fit, one filter, so this script cannot test
# a different set of genes from the one fig. 3D draws.
from cbr_changepoint_relax import fit_arms, omega_filter_counts, REGIONS

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Gene sets every test is run on.
GENE_SETS = ('all', 'clear')


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def cross_check(df, args):
    """Recompute the signed distance from the BED, as cbr_changepoint_relax.py
    does before fitting, and stop if the table disagrees.

    Returns (cbrs_by_chrom, elim_merged, diagnostics). The CBRs and merged
    tracts are reused to place each listed gene against its boundary.
    """
    elim = parse_bed(args.eliminated_bed)
    merged = merge_intervals(elim)
    sizes, inferred = load_chrom_sizes(args.chrom_sizes, elim)
    cbrs, cbr_diag = derive_cbrs(merged, sizes, args.terminal_tol)
    recomputed = np.full(len(df), np.nan)
    have = df['chrom1'].notna() & df['midpoint1'].notna()
    for chrom, block in df[have].groupby('chrom1', sort=False):
        sd, _ = signed_distance(chrom, block['midpoint1'].to_numpy(float),
                                cbrs, merged)
        recomputed[df.index.get_indexer(block.index)] = sd
    both = ~np.isnan(recomputed) & df[args.dist_col].notna().to_numpy()
    delta = np.abs(recomputed[both] - df[args.dist_col].to_numpy(float)[both])
    n_bad = int((delta > args.check_tol_bp).sum())
    if n_bad:
        sys.exit(f'ERROR: the table disagrees with a fresh computation of '
                 f'{args.dist_col} for {n_bad:,} of {int(both.sum()):,} genes '
                 f'(max |delta| {delta.max():,.0f} bp): it was built with a '
                 f'different CBR definition or --terminal-tol. Use the BED and '
                 f'tolerance of annotate_regions.py.')
    diag = {'n_compared': int(both.sum()), 'n_disagreeing': 0,
            'max_abs_delta_bp': float(delta.max()) if delta.size else 0.0,
            'cbrs': cbr_diag, 'chrom_sizes_inferred': inferred}
    logging.info('  cross-check: %s genes agree (max |delta| %.1f bp)',
                 f'{diag["n_compared"]:,}', diag['max_abs_delta_bp'])
    return cbrs, merged, diag


def load_genes(args):
    """The genes and log10(dN/dS) values fig. 3D uses, with x in kb.

    Mirrors the 'prepare' block of cbr_changepoint_relax.main() line for line:
    the two region classes, dN/dS finite and > 0, inside [omega_min,
    omega_max], then the genes with a distance (a CBR on their sequence).
    """
    df = pd.read_csv(args.annotated, sep='\t')
    logging.info('Loaded %s rows from %s', f'{len(df):,}',
                 os.path.basename(args.annotated))
    for col in (args.dist_col, args.omega_col, args.region_col, 'chrom1',
                'midpoint1'):
        if col not in df.columns:
            sys.exit(f'ERROR: {col!r} is not a column of {args.annotated}.')
    d = df[df[args.region_col].isin(REGIONS)].copy()
    d[args.omega_col] = pd.to_numeric(d[args.omega_col], errors='coerce')
    filt = omega_filter_counts(d, args.omega_col, args.region_col,
                               args.omega_min, args.omega_max)
    d = d[np.isfinite(d[args.omega_col]) & (d[args.omega_col] > 0)]
    d = d[(d[args.omega_col] >= args.omega_min) &
          (d[args.omega_col] <= args.omega_max)].copy()
    for r in REGIONS:
        if int((d[args.region_col] == r).sum()) != int(filt.loc[r, 'usable']):
            sys.exit('ERROR: the filter accounting disagrees with the filtered '
                     'table; this script no longer mirrors '
                     'cbr_changepoint_relax.py.')
    d['logw'] = np.log10(d[args.omega_col])
    sp = d.dropna(subset=[args.dist_col]).copy()
    sp['x_kb'] = sp[args.dist_col].to_numpy(float) / 1000.0
    logging.info('  %s usable genes, %s with a CBR on their sequence',
                 f'{len(d):,}', f'{len(sp):,}')
    return df, sp, filt


def check_against_stats(arms, stats_path):
    """Compare the arm fits with the figure run's cbr_changepoint_stats.tsv.

    The table is written with float_format '%.6g', so each slope is compared
    at that precision (formatted the same way), and n exactly. Returns a list
    of lines for the summary; exits if a slope or n differs, because then this
    script is not looking at the genes the figure drew.
    """
    st = pd.read_csv(stats_path, sep='\t').set_index('quantity')
    lines, bad = [], []
    for name in REGIONS:
        q = f'slope_{name}'
        a = arms.get(name)
        if q not in st.index or a is None:
            bad.append(f'{q}: missing here or in {os.path.basename(stats_path)}')
            continue
        est, n = float(st.loc[q, 'estimate']), int(st.loc[q, 'n'])
        ok = (float(f'{a["slope_per_100kb"]:.6g}') == est and a['n'] == n)
        lines.append(f'  {name:<10} slope {a["slope_per_100kb"]:+.4f} (figure '
                     f'{est:+.4f}), n {a["n"]:,} (figure {n:,})  '
                     f'{"same" if ok else "DIFFERENT"}')
        if not ok:
            bad.append(lines[-1].strip())
    if bad:
        sys.exit('ERROR: the arm fits differ from the figure run '
                 f'({stats_path}):\n  ' + '\n  '.join(bad) +
                 '\n  Use the table, BED and options of that run.')
    return lines


# ---------------------------------------------------------------------------
# Bins and tests
# ---------------------------------------------------------------------------

def figure_bins(x, y, w, window_kb, n_bins, min_bin):
    """The bins of plot_slope_step() in cbr_changepoint_relax.py, with counts.

    Same edges (np.linspace over +/- window, n_bins + 1 edges), same membership
    (|x| <= window and lo <= x < hi), same drawing rule (n >= min_bin). The
    figure computes no mean for a bin it does not draw; here every bin is
    reported and `drawn` says which ones the figure shows.
    """
    keep = np.abs(x) <= window_kb
    edges = np.linspace(-window_kb, window_kb, n_bins + 1)
    rows = []
    for j in range(n_bins):
        lo, hi = edges[j], edges[j + 1]
        m = keep & (x >= lo) & (x < hi)
        n = int(m.sum())
        arm = ('retained' if hi <= 0 else 'eliminated' if lo >= 0
               else 'straddles_cbr')
        rows.append({
            'bin': j, 'arm': arm, 'lo_kb': lo, 'hi_kb': hi,
            'centre_kb': 0.5 * (lo + hi), 'n': n,
            'mean_log10_dnds': float(y[m].mean()) if n else np.nan,
            'se_log10_dnds': (float(y[m].std(ddof=1) / np.sqrt(n))
                              if n >= 2 else np.nan),
            'median_dnds': float(np.median(w[m])) if n else np.nan,
            'drawn': n >= min_bin})
    return pd.DataFrame(rows), edges


def arm_masks(x, arm, width, window_kb):
    """Edge and rest of one arm, with the figure's bin convention [lo, hi).

    retained:   edge -width <= x < 0,  rest -window <= x < -width
    eliminated: edge  0 <= x < width,  rest  width <= x <= window
    """
    if arm == 'retained':
        edge = (x >= -width) & (x < 0)
        rest = (x >= -window_kb) & (x < -width)
    else:
        edge = (x >= 0) & (x < width)
        rest = (x >= width) & (x <= window_kb)
    return edge, rest


#: plot_slope_step() in cbr_changepoint_relax.py widens its percentile limits by
#: this share of the range between them (`pad = 0.16 * (yhi - ylo)`).
YLIM_PAD = 0.16
#: Its default --cap-pct.
FIGURE_CAP_PCT = 99.0


def figure_ylim(y_drawn, cap_pct):
    """The fig. 3D y-limits, as plot_slope_step() sets them: the (100 - cap)th
    and cap-th percentiles of the drawn genes' log10(dN/dS) (numpy's default
    linear interpolation, as there), each widened by YLIM_PAD of their range."""
    ylo = float(np.percentile(y_drawn, 100.0 - cap_pct))
    yhi = float(np.percentile(y_drawn, cap_pct))
    pad = YLIM_PAD * (yhi - ylo)
    return ylo - pad, yhi + pad


def resolve_cap_pct(args):
    """--cap-pct if given; else the cap_pct of the figure run, from the
    provenance.json beside --stats; else the figure script's default."""
    if args.cap_pct is not None:
        return float(args.cap_pct), '--cap-pct'
    if args.stats:
        prov = os.path.join(os.path.dirname(os.path.abspath(args.stats)),
                            'provenance.json')
        if os.path.exists(prov):
            try:
                with open(prov) as fh:
                    return float(json.load(fh)['args']['cap_pct']), \
                        f'the figure run ({prov})'
            except (OSError, KeyError, TypeError, ValueError) as exc:
                logging.warning('  could not read cap_pct from %s (%s)', prov, exc)
    return FIGURE_CAP_PCT, ('the figure script default (no figure-run '
                            'provenance found)')


def off_axis(sp, x, y, window_kb, cap_pct):
    """Genes of the fig. 3D window whose log10(dN/dS) lies beyond its y-axis.

    Returns (info dict, frame of those genes). A point exactly at a limit is
    drawn, so only values strictly beyond count."""
    keep = np.abs(x) <= window_kb
    lo, hi = figure_ylim(y[keep], cap_pct)
    below, above = keep & (y < lo), keep & (y > hi)
    arm = np.where(x < 0, 'retained', 'eliminated')
    info = {'cap_pct': cap_pct, 'ylim_log10': [lo, hi],
            'n_drawn_window': int(keep.sum())}
    for name, m in (('below', below), ('above', above)):
        for r in REGIONS:
            info[f'n_{name}_{r}'] = int((m & (arm == r)).sum())
        info[f'n_{name}'] = int(m.sum())
    info['n_off'] = info['n_below'] + info['n_above']
    cols = [c for c in ('gene1', 'chrom1', 'start1', 'end1', 'midpoint1')
            if c in sp.columns]
    g = sp.loc[below | above, cols].copy()
    g.insert(0, 'arm', arm[below | above])
    g.insert(1, 'side', np.where(below[below | above], 'below', 'above'))
    g['x_kb'] = np.round(x[below | above], 3)
    g['log10_dnds'] = np.round(y[below | above], 4)
    g['dnds'] = np.round(10.0 ** y[below | above], 6)
    return info, g.sort_values(['side', 'log10_dnds'])


def boot_median_ratio(a, b, rng, n_boot, level=0.95):
    """Percentile interval of median(a) / median(b), each group resampled on
    its own. Resampled in chunks so memory stays small at any n."""
    out = []
    chunk = 250
    for start in range(0, n_boot, chunk):
        k = min(chunk, n_boot - start)
        ma = np.median(a[rng.integers(0, a.size, (k, a.size))], axis=1)
        mb = np.median(b[rng.integers(0, b.size, (k, b.size))], axis=1)
        good = mb > 0
        out.append(ma[good] / mb[good])
    v = np.concatenate(out) if out else np.array([])
    if v.size == 0:
        return np.nan, np.nan
    tail = 100 * (1 - level) / 2
    return tuple(float(t) for t in np.percentile(v, [tail, 100 - tail]))


def compare(sp, edge, rest, arm, arms, rng, n_boot, row, min_n=3):
    """Statistics of one edge (or crossing) set against the rest of its arm.

    Fills `row` with the means of log10 dN/dS and Welch's P, the medians of
    dN/dS with their ratio and bootstrap interval, Mann-Whitney's P, and the
    mean residual of the edge genes from the arm's fitted line with a
    one-sample t-test P.
    """
    x = sp['x_kb'].to_numpy(float)
    y = sp['logw'].to_numpy(float)
    w = 10.0 ** y
    ne, nr = int(edge.sum()), int(rest.sum())
    row.update({'n_edge': ne, 'n_rest': nr})
    if ne < min_n or nr < min_n:
        row['note'] = f'fewer than {min_n} genes in a group; not tested'
        return row
    ye, yr = y[edge], y[rest]
    a = arms.get(arm)
    resid = (ye - (a['intercept'] + a['slope_per_100kb'] / 100.0 * x[edge])
             if a else np.array([]))
    lo, hi = boot_median_ratio(w[edge], w[rest], rng, n_boot)
    row.update({
        'mean_log10_edge': float(ye.mean()), 'mean_log10_rest': float(yr.mean()),
        'diff_log10': float(ye.mean() - yr.mean()),
        'p_welch': float(ttest_ind(ye, yr, equal_var=False).pvalue),
        'median_dnds_edge': float(np.median(w[edge])),
        'median_dnds_rest': float(np.median(w[rest])),
        'median_ratio': float(np.median(w[edge]) / np.median(w[rest])),
        'ratio_lo95': lo, 'ratio_hi95': hi,
        'p_mwu': float(mannwhitneyu(ye, yr, alternative='two-sided').pvalue),
        'resid_mean_log10': float(resid.mean()) if resid.size else np.nan,
        'p_resid': (float(ttest_1samp(resid, 0.0).pvalue)
                    if resid.size >= min_n else np.nan),
        'note': ''})
    return row


def edge_test(sp, arm, width, window_kb, arms, clear_mask, gene_set, rng,
              n_boot):
    """One row of cbr_edge_tests.tsv: the edge of one arm against its rest."""
    x = sp['x_kb'].to_numpy(float)
    edge, rest = arm_masks(x, arm, width, window_kb)
    if gene_set == 'clear':
        edge &= clear_mask
        rest &= clear_mask
    row = {'arm': arm, 'width_kb': width, 'gene_set': gene_set}
    return compare(sp, edge, rest, arm, arms, rng, n_boot, row)


def crossing_test(sp, arm, window_kb, arms, clear_mask, rng, n_boot):
    """The genes of one arm whose span reaches across their CBR, against the
    arm's clear genes (window as the figure's). For a retained gene, crossing
    means part of its coding sequence lies in eliminated DNA."""
    x = sp['x_kb'].to_numpy(float)
    side = (x < 0) if arm == 'retained' else (x > 0)
    inside = side & (np.abs(x) <= window_kb)
    cross = sp['span_crosses_cbr'].astype(bool).to_numpy()
    row = {'arm': arm, 'width_kb': np.nan, 'gene_set': 'crossing'}
    return compare(sp, inside & cross, inside & clear_mask, arm, arms, rng,
                   n_boot, row)


# ---------------------------------------------------------------------------
# Genes against the boundary
# ---------------------------------------------------------------------------

def overlap_bp(start, end, starts, ends):
    """bp of [start, end] (1-based inclusive) inside a merged interval set."""
    if starts.size == 0:
        return 0
    lo = np.searchsorted(ends, start)                 # first interval ending >= start
    hi = np.searchsorted(starts, end, side='right')   # past the last starting <= end
    tot = 0.0
    for s, e in zip(starts[lo:hi], ends[lo:hi]):
        tot += max(0.0, min(end, e) - max(start, s) + 1)
    return int(tot)


def place_genes(sp, cbrs, merged, flag_sets, args):
    """Per gene: nearest CBR, span across it, gap to it, overlap with flags.

    Adds columns cbr_pos, span_crosses_cbr, span_far_side_bp, gap_to_cbr_kb and
    one <flag>_bp column per --flag-bed. Genes without start1/end1 get NaN, and
    then 'clear' cannot be formed (the summary says so).
    """
    have_span = {'start1', 'end1'}.issubset(sp.columns)
    sp['cbr_pos'] = np.nan
    for chrom, block in sp.groupby('chrom1', sort=False):
        c = cbrs.get(chrom)
        if c is None or c.size == 0:
            continue
        idx, _ = nearest_cbr(block['midpoint1'].to_numpy(float), c)
        sp.loc[block.index, 'cbr_pos'] = c[idx]
    for col in ('span_crosses_cbr', 'span_far_side_bp', 'gap_to_cbr_kb'):
        sp[col] = np.nan
    if not have_span:
        logging.warning('  the table has no start1/end1: spans cannot be '
                        'placed against the boundary, and the "clear" gene '
                        'set is not formed')
        return sp, False
    s = pd.to_numeric(sp['start1'], errors='coerce').to_numpy(float)
    e = pd.to_numeric(sp['end1'], errors='coerce').to_numpy(float)
    c = sp['cbr_pos'].to_numpy(float)
    crosses = (s <= c) & (c <= e)
    sp['span_crosses_cbr'] = crosses
    # bp of the span on the far side of the nearest CBR (from the midpoint):
    # the CBR sits at a half-base junction, so e - c + 0.5 counts the bases
    # beyond a CBR to the right of the midpoint, c - s + 0.5 those to the left.
    x = sp['x_kb'].to_numpy(float)
    mid = sp['midpoint1'].to_numpy(float)
    far = np.where(c >= mid, np.clip(e - c + 0.5, 0, None),
                   np.clip(c - s + 0.5, 0, None))
    sp['span_far_side_bp'] = np.where(crosses, far, 0.0)
    gap = np.minimum(np.abs(s - c), np.abs(e - c)) / 1000.0
    sp['gap_to_cbr_kb'] = np.where(crosses, 0.0, gap)
    # retained genes: bp of the span inside any eliminated tract (not only the
    # nearest one); eliminated genes: bp outside every tract.
    elim_bp = np.zeros(len(sp))
    for i, (chrom, a, b) in enumerate(zip(sp['chrom1'], s, e)):
        st, en = merged.get(chrom, (np.array([]), np.array([])))
        if np.isfinite(a) and np.isfinite(b):
            elim_bp[i] = overlap_bp(a, b, st, en)
    span = np.where(np.isfinite(s) & np.isfinite(e), e - s + 1, np.nan)
    sp['span_bp'] = span
    sp['span_other_region_bp'] = np.where(x < 0, elim_bp, span - elim_bp)
    for name, fm in flag_sets.items():
        fl = np.zeros(len(sp))
        for i, (chrom, a, b) in enumerate(zip(sp['chrom1'], s, e)):
            st, en = fm.get(chrom, (np.array([]), np.array([])))
            if np.isfinite(a) and np.isfinite(b):
                fl[i] = overlap_bp(a, b, st, en)
        sp[f'{name}_bp'] = fl
    return sp, True


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def fmt(v, nd=3):
    return 'NA' if v is None or not np.isfinite(v) else f'{v:.{nd}f}'


def ptxt(p):
    """A P or q value as plain text for the summary: 0.677, 2.1e-09.

    region_lib.fmt_p writes mathtext for figures, which reads badly in a text
    file."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p == 0:
        return '< 1e-300'
    return f'{p:.3g}' if p >= 0.001 else f'{p:.1e}'


def write_summary(path, args, filt, arms, stats_lines, bins, tests, genes,
                  bin_width, have_span, flag_names, offax=None, cap_src=''):
    L = []
    L.append('Fig. 3D edge check: the bin next to the CBR against the rest of '
             'its arm')
    L.append('=' * 78)
    L.append(f'table        {args.annotated}')
    L.append(f'window       +/- {args.slope_window:g} kb, {args.n_bins} bins of '
             f'{bin_width:g} kb, drawn where n >= {args.min_bin}')
    L.append(f'dN/dS        {args.omega_min:g} to {args.omega_max:g} '
             f'({args.omega_col}); log10 for means and fits')
    L.append('')
    L.append('Genes (as cbr_changepoint_relax.py counts them)')
    for r in REGIONS:
        f = filt.loc[r]
        L.append(f'  {r:<10} {int(f.in_table):,} in table, {int(f.usable):,} '
                 f'usable')
    L.append('')
    L.append('Arm fits (fit_arms, imported from cbr_changepoint_relax.py)')
    for r in REGIONS:
        a = arms.get(r)
        if a:
            L.append(f'  {r:<10} slope {a["slope_per_100kb"]:+.4f} log10 per '
                     f'100 kb, P = {ptxt(a["p"])}, n = {a["n"]:,}')
    if stats_lines:
        L.append('  against the figure run (--stats):')
        L.extend(stats_lines)
    else:
        L.append('  (no --stats: not compared with the figure run)')
    L.append('')
    L.append('Bins nearest the CBR')
    near = bins[(bins['hi_kb'] > -4 * bin_width) & (bins['lo_kb'] < 4 * bin_width)]
    L.append(f'  {"bin (kb)":<18}{"arm":<12}{"n":>5}{"mean log10":>12}'
             f'{"SE":>8}{"median":>9}  drawn')
    for _, b in near.iterrows():
        L.append(f'  [{b.lo_kb:+7.2f}, {b.hi_kb:+7.2f}){"":2}{b.arm:<12}'
                 f'{int(b.n):>5}{fmt(b.mean_log10_dnds):>12}'
                 f'{fmt(b.se_log10_dnds):>8}{fmt(b.median_dnds):>9}  '
                 f'{"yes" if b.drawn else "no"}')
    L.append('')
    L.append('Edge against the rest of the arm (q: BH over every row)')
    for _, t in tests.iterrows():
        what = ('crossing genes' if t.gene_set == 'crossing' else
                 f'edge {t.width_kb:>5g} kb  {t.gene_set:<5}')
        head = (f'  {t.arm:<10} {what:<20} n {int(t.n_edge):>3} vs '
                f'{int(t.n_rest):>4}')
        if isinstance(t.get('note'), str) and t.note:
            L.append(head + f'  {t.note}')
            continue
        L.append(head + f'  mean log10 {t.mean_log10_edge:+.3f} vs '
                 f'{t.mean_log10_rest:+.3f} (diff {t.diff_log10:+.3f}, '
                 f'Welch P = {ptxt(t.p_welch)}, q = {ptxt(t.q)})')
        L.append(f'  {"":<41}median dN/dS {t.median_dnds_edge:.3f} vs '
                 f'{t.median_dnds_rest:.3f}, ratio {t.median_ratio:.2f} '
                 f'[{fmt(t.ratio_lo95, 2)}, {fmt(t.ratio_hi95, 2)}], '
                 f'MWU P = {ptxt(t.p_mwu)}')
        L.append(f'  {"":<41}off the arm\'s line by {t.resid_mean_log10:+.3f} '
                 f'(P = {ptxt(t.p_resid)})')
    L.append('')
    L.append('Genes in the edge bins')
    for arm in REGIONS:
        g = genes[genes['arm'] == arm]
        L.append(f'  {arm}: {len(g)} genes from {g["cbr_key"].nunique()} CBRs '
                 f'(most from one CBR: '
                 f'{int(g["cbr_key"].value_counts().max()) if len(g) else 0})')
        if have_span and len(g):
            n_cross = int(g['span_crosses_cbr'].astype(bool).sum())
            L.append(f'    {n_cross} reach across the CBR; '
                     f'{int((g["span_other_region_bp"] > 0).sum())} have part '
                     f'of their span in the other region class')
            for fn in flag_names:
                L.append(f'    {int((g[f"{fn}_bp"] > 0).sum())} overlap {fn}')
    L.append('  (cbr_edge_genes.tsv lists them)')
    L.append('')
    if offax is not None:
        lo, hi = offax['ylim_log10']
        L.append(f'Fig. 3D y-axis (plot_slope_step rule: percentiles '
                 f'{100 - offax["cap_pct"]:g} and {offax["cap_pct"]:g} of the '
                 f'{offax["n_drawn_window"]} drawn genes, each widened by '
                 f'{100 * YLIM_PAD:g}% of their range; cap from {cap_src})')
        L.append(f'  limits     log10 dN/dS {lo:+.3f} to {hi:+.3f} '
                 f'(dN/dS {10 ** lo:.3g} to {10 ** hi:.3g})')
        for name in ('below', 'above'):
            parts = ', '.join(f'{offax[f"n_{name}_{r}"]} {r}' for r in REGIONS)
            L.append(f'  {name:<10} {offax[f"n_{name}"]} genes ({parts})')
        L.append('  (not drawn; cbr_edge_offaxis_genes.tsv lists them)')
        L.append('')
    # What to say: the drawn retained bin, all genes and clear genes.
    L.append('What to say')
    sel = tests[(tests['arm'] == 'retained') & np.isclose(tests['width_kb'],
                                                          bin_width)]
    for _, t in sel.iterrows():
        if isinstance(t.get('note'), str) and t.note:
            L.append(f'  {t.gene_set}: not tested ({t.note})')
            continue
        # One pre-specified question (the drawn retained bin), so its own P
        # decides; q, over every row of the table, is printed beside it.
        verdict = ('differs from' if t.p_welch < 0.05 else
                   'does not differ from')
        L.append(f'  {t.gene_set}: the retained bin within {bin_width:g} kb of '
                 f'the CBR ({int(t.n_edge)} genes) {verdict} the rest of the '
                 f'retained arm (Welch P = {ptxt(t.p_welch)}, q = '
                 f'{ptxt(t.q)}; median dN/dS {t.median_ratio:.2f}-fold '
                 f'[{fmt(t.ratio_lo95, 2)}, {fmt(t.ratio_hi95, 2)}]).')
    cr = tests[(tests['arm'] == 'retained') & (tests['gene_set'] == 'crossing')]
    for _, t in cr.iterrows():
        if isinstance(t.get('note'), str) and t.note:
            L.append(f'  crossing: not tested ({t.note})')
            continue
        verdict = 'evolve faster than' if (t.p_welch < 0.05 and
                                            t.diff_log10 > 0) else \
            'evolve slower than' if t.p_welch < 0.05 else \
            'do not differ from'
        L.append(f'  crossing: the {int(t.n_edge)} retained genes whose span '
                 f'reaches across the CBR into eliminated DNA {verdict} the '
                 f'other retained genes of the window (Welch P = '
                 f'{ptxt(t.p_welch)}; median dN/dS {t.median_ratio:.2f}-fold '
                 f'[{fmt(t.ratio_lo95, 2)}, {fmt(t.ratio_hi95, 2)}]).')
    if not have_span:
        L.append('  No spans in the table, so partly eliminated genes could '
                 'not be checked.')
    if offax is not None:
        n, tot = offax['n_off'], offax['n_drawn_window']
        if n == 0:
            L.append(f'  3D legend: all {tot} genes lie within the y-axis, so '
                     'every gene is drawn.')
        else:
            by = [f'{offax[f"n_{s}_{r}"]} {r} {s}' for s in ('below', 'above')
                  for r in REGIONS if offax[f'n_{s}_{r}']]
            L.append(f'  3D legend: {n} of the {tot} genes lie beyond the y-axis '
                     f'and are not drawn ({", ".join(by)}); every fit and mean '
                     'uses them.')
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    return L


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Test the fig. 3D bins next to the CBR against the rest of '
                    'their arm, and list their genes.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', required=True,
                   help='kaks_annotated.tsv from annotate_regions.py, the table '
                        'fig. 3D was drawn from.')
    p.add_argument('--eliminated-bed', required=True,
                   help='Eliminated regions BED (as for cbr_changepoint_relax.py). '
                        'Required: it places each gene against its boundary.')
    p.add_argument('--chrom-sizes', default=None)
    p.add_argument('--terminal-tol', type=int, default=1000,
                   help='Must match annotate_regions.py (default 1000).')
    p.add_argument('--check-tol-bp', type=float, default=10.0)
    p.add_argument('--stats', default=None,
                   help='cbr_changepoint_stats.tsv of the figure run; the arm '
                        'fits must match it.')
    p.add_argument('--dist-col', default='signed_cbr_dist1')
    p.add_argument('--omega-col', default='Ka_Ks')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--slope-window', type=float, default=150.0)
    p.add_argument('--omega-min', type=float, default=1e-4)
    p.add_argument('--omega-max', type=float, default=10.0)
    p.add_argument('--n-bins', type=int, default=24)
    p.add_argument('--min-bin', type=int, default=3,
                   help='Genes a bin needs to be drawn (the figure uses 3).')
    p.add_argument('--edge-widths', default='12.5,25,50',
                   help='Comma-separated edge widths in kb; the drawn bin width '
                        'is added if missing.')
    p.add_argument('--flag-bed', action='append', default=[],
                   help='BED of intervals to flag (e.g. CBR zones). Genes '
                        'overlapping one leave the "clear" set. Repeatable.')
    p.add_argument('--cap-pct', type=float, default=None,
                   help='The y-axis percentile cap of the figure run (default: '
                        'read from the provenance.json beside --stats, else 99, '
                        "cbr_changepoint_relax.py's default).")
    p.add_argument('--n-boot', type=int, default=2000)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None)
    a = p.parse_args(argv)
    if a.cap_pct is not None and not 50.0 < a.cap_pct <= 100.0:
        p.error('--cap-pct must lie in (50, 100]')
    if a.outdir is None:
        a.outdir = os.path.join('results', 'cbr_edge_check_' +
                                datetime.now().strftime('%Y%m%d_%H%M%S'))
    if a.n_bins % 2:
        p.error('--n-bins must be even, so that no bin straddles the CBR')
    try:
        a.widths = sorted({float(v) for v in a.edge_widths.split(',') if v.strip()})
    except ValueError:
        p.error('--edge-widths takes numbers, e.g. 12.5,25,50')
    if any(v <= 0 or v >= a.slope_window for v in a.widths):
        p.error('every edge width must lie between 0 and --slope-window')
    return a


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'cbr_edge_check')
    logging.info('cbr_edge_check.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    rng = np.random.default_rng(args.seed)

    df, sp, filt = load_genes(args)
    logging.info('Cross-checking %s against a fresh computation', args.dist_col)
    cbrs, merged, check_diag = cross_check(df, args)

    x = sp['x_kb'].to_numpy(float)
    y = sp['logw'].to_numpy(float)
    arms, step = fit_arms(x, y, args.slope_window)
    for r in REGIONS:
        a = arms.get(r)
        if a:
            logging.info('  %s arm: slope %+.4f per 100 kb, n = %d', r,
                         a['slope_per_100kb'], a['n'])
    stats_lines = (check_against_stats(arms, args.stats) if args.stats
                   else [])

    bin_width = 2.0 * args.slope_window / args.n_bins
    widths = sorted(set(args.widths) | {bin_width})
    bins, _ = figure_bins(x, y, 10.0 ** y, args.slope_window, args.n_bins,
                          args.min_bin)

    flag_sets = {}
    for path in args.flag_bed:
        name = os.path.splitext(os.path.basename(path))[0].replace('.', '_')
        flag_sets[name] = merge_intervals(parse_bed(path))
    sp, have_span = place_genes(sp, cbrs, merged, flag_sets, args)
    clear = np.ones(len(sp), dtype=bool)
    if have_span:
        clear &= ~sp['span_crosses_cbr'].astype(bool).to_numpy()
        clear &= ~(sp['span_other_region_bp'].to_numpy(float) > 0)
        for name in flag_sets:
            clear &= ~(sp[f'{name}_bp'].to_numpy(float) > 0)
    sets = GENE_SETS if have_span else ('all',)

    rows = []
    for arm in REGIONS:
        for wdt in widths:
            for gs in sets:
                rows.append(edge_test(sp, arm, wdt, args.slope_window, arms,
                                      clear, gs, rng, args.n_boot))
        if have_span:
            rows.append(crossing_test(sp, arm, args.slope_window, arms, clear,
                                      rng, args.n_boot))
    tests = pd.DataFrame(rows)
    tests['q'] = benjamini_hochberg(tests['p_welch'].to_numpy(float)
                                    if 'p_welch' in tests else
                                    np.full(len(tests), np.nan))
    tests['drawn_bin'] = np.isclose(tests['width_kb'].to_numpy(float),
                                    bin_width)

    # genes of the two edge bins
    keep_cols = [c for c in ('gene1', 'chrom1', 'start1', 'end1', 'midpoint1')
                 if c in sp.columns]
    rate_cols = [c for c in ('Ka', 'Ks', args.omega_col) if c in sp.columns]
    extra = [c for c in ('cbr_pos', 'span_bp', 'span_crosses_cbr',
                         'span_far_side_bp', 'span_other_region_bp',
                         'gap_to_cbr_kb') if c in sp.columns]
    extra += [f'{n}_bp' for n in flag_sets]
    parts = []
    for arm in REGIONS:
        edge, _ = arm_masks(x, arm, bin_width, args.slope_window)
        g = sp.loc[edge, keep_cols + rate_cols + ['x_kb', 'logw'] + extra].copy()
        g.insert(0, 'arm', arm)
        parts.append(g)
    genes = pd.concat(parts, ignore_index=True)
    genes['cbr_key'] = genes['chrom1'].astype(str) + ':' + \
        genes['cbr_pos'].map(lambda v: 'NA' if not np.isfinite(v) else f'{v:.1f}')
    genes = genes.sort_values(['arm', 'x_kb'], ascending=[False, True])

    # genes beyond the fig. 3D y-axis
    cap_pct, cap_src = resolve_cap_pct(args)
    offax, offax_genes = off_axis(sp, x, y, args.slope_window, cap_pct)
    logging.info('Fig. 3D y-axis %+.3f to %+.3f (cap %g, from %s): %d genes '
                 'below, %d above', *offax['ylim_log10'], cap_pct, cap_src,
                 offax['n_below'], offax['n_above'])

    bins.to_csv(os.path.join(args.outdir, 'cbr_edge_bins.tsv'), sep='\t',
                index=False, float_format='%.6g')
    tests.to_csv(os.path.join(args.outdir, 'cbr_edge_tests.tsv'), sep='\t',
                 index=False, float_format='%.6g')
    # Coordinates at full precision ('%.6g' would print 1.36e+07); only the
    # derived columns are rounded.
    for col, nd in (('x_kb', 3), ('logw', 4), ('gap_to_cbr_kb', 3)):
        if col in genes.columns:
            genes[col] = genes[col].round(nd)
    genes.to_csv(os.path.join(args.outdir, 'cbr_edge_genes.tsv'), sep='\t',
                 index=False)
    # Coordinates at full precision, as cbr_edge_genes.tsv (derived columns are
    # rounded in off_axis()).
    offax_genes.to_csv(os.path.join(args.outdir, 'cbr_edge_offaxis_genes.tsv'),
                       sep='\t', index=False)
    lines = write_summary(os.path.join(args.outdir, 'cbr_edge_summary.txt'),
                          args, filt, arms, stats_lines, bins, tests, genes,
                          bin_width, have_span, list(flag_sets), offax, cap_src)
    for line in lines:
        logging.info(line)

    prov = {
        'script': os.path.abspath(__file__), 'git': git_commit(_HERE),
        'run_utc': datetime.now(timezone.utc).isoformat(),
        'argv': sys.argv, 'args': {k: v for k, v in vars(args).items()},
        'inputs': {'annotated': fingerprint(args.annotated),
                   'eliminated_bed': fingerprint(args.eliminated_bed),
                   'chrom_sizes': fingerprint(args.chrom_sizes),
                   'stats': fingerprint(args.stats),
                   'flag_beds': [fingerprint(f) for f in args.flag_bed]},
        'cross_check': check_diag,
        'bin_width_kb': bin_width, 'widths_kb': widths,
        'gene_sets': list(sets),
        'off_axis_3d': {**offax, 'cap_source': cap_src},
        'versions': {'python': platform.python_version(),
                     'numpy': np.__version__, 'pandas': pd.__version__,
                     'scipy': scipy.__version__}}
    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump(prov, fh, indent=2, default=str)
    logging.info('Done.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
