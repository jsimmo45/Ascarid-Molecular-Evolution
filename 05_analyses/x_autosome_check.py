#!/usr/bin/env python3
"""
x_autosome_check.py -- Ascaris X chromosomes versus autosomes, and the
eliminated-versus-retained contrast within each.

Two questions from the 2026-10-10 review of the optional sex-chromosome
sentence (results draft, change note 6). The pairwise Ascaris x P. univalens
rates are classified by the Ascaris germline chromosome each gene sits on.

1. DO X-LINKED GENES DIFFER FROM AUTOSOMAL GENES IN dS, dN AND dN/dS?
   The draft says both dS and dN/dS are lower on the five X chromosomes, but
   no test was run and the per-chromosome dN/dS ranges overlap. Tested here on
   RETAINED genes only, so the comparison cannot be produced by eliminated genes
   (dN/dS 1.8x higher) sitting in different proportions on the X and the
   autosomes. Two tests per metric:

     gene level        Mann-Whitney, X-linked genes vs autosomal genes.
     chromosome level  exact permutation over every way of labelling k of the
                       N chromosomes "X" (C(24, 5) = 42,504 for Ascaris), on
                       the difference in mean per-chromosome median.

   Genes on one chromosome share its mutational environment -- base
   composition, recombination, replication timing -- so they are not
   independent draws, and the gene-level P overstates the evidence. THE
   CHROMOSOME-LEVEL TEST IS THE ONE TO QUOTE. Its floor is 1 / C(N, k), about
   2e-5 two-sided for Ascaris, and it is unchanged by how many genes each
   chromosome carries.

   The X/A ratio of median dS, r, is converted to the male-to-female mutation
   ratio it implies under male-driven evolution: with every X carried once in
   males and twice in females, r = 2(2 + a) / (3(1 + a)), so
   a = (4 - 3r) / (3r - 2). That assumes the whole dS gap is mutational. The
   interval comes from resampling chromosomes within each class, the same unit
   as the permutation test, so it is wide with only five X chromosomes, and
   honestly so.

2. DOES THE ELIMINATED-VERSUS-RETAINED CONTRAST HOLD WITHIN THE AUTOSOMES?
   If X-linked genes have lower dN/dS and eliminated genes sat mostly on
   autosomes, part of the genome-wide 1.77-fold contrast could be chromosome
   composition. The contrast is recomputed within the autosomes and within the
   X, beside the genome-wide value (all pairs, as Results quotes it), and the
   X-linked share of each region is compared with Fisher's exact test.

No figure: the outputs are numbers for one Results sentence and, if it is
kept, a supplementary table (per_chromosome.tsv).

Outputs (in --outdir)
---------------------
  x_autosome_summary.txt   the tests, and what the text can and cannot claim
  x_autosome_tests.tsv     every test, with BH q within each family
  per_chromosome.tsv       per-chromosome n and medians (retained genes)
  x_autosome_check.log
  provenance.json

Usage
-----
  python x_autosome_check.py \\
    --annotated results/annotate_regions_20260817/kaks_annotated.tsv \\
    --outdir results/x_autosome_check_<stamp>

Dependencies: numpy, pandas, scipy; region_lib and rate_tables from lib/.
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import math
import os
import platform
import re
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

from rate_tables import read_rate_table
from region_lib import (benjamini_hochberg, fingerprint, git_commit,
                        setup_logging)

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Metrics in manuscript notation (rate_tables maps Ka/Ks/Ka_Ks onto these).
METRICS = ('dS', 'dN', 'dN_dS')
METRIC_LABEL = {'dS': 'dS', 'dN': 'dN', 'dN_dS': 'dN/dS'}

#: Above this many label assignments, the chromosome-level test switches from
#: exact enumeration to Monte Carlo. Ascaris (C(24, 5) = 42,504) is exact.
MAX_EXACT = 2_000_000


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def p_text(p):
    """Plain-text P for the summary (region_lib.fmt_p writes mathtext)."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p == 0:
        return '< 1e-300'
    return f'{p:.2g}' if p < 0.001 else f'{p:.3f}'


def classify_chromosomes(chroms, x_pattern, autosome_pattern):
    """'X', 'A' or 'other' for each chromosome name.

    'other' covers unplaced scaffolds and anything the two patterns miss; those
    genes are excluded from the X/A comparisons but kept in the genome-wide
    contrast, which is computed on all pairs exactly as Results quotes it.
    """
    xr, ar = re.compile(x_pattern), re.compile(autosome_pattern)
    out = []
    for c in chroms:
        s = '' if pd.isna(c) else str(c)
        out.append('X' if xr.search(s) else 'A' if ar.search(s) else 'other')
    return np.array(out, dtype=object)


def male_bias_alpha(r):
    """Male-to-female mutation ratio implied by an X/A neutral-rate ratio r.

    With every X carried once in males and twice in females, a site on the X
    spends 1/3 of generations in males and an autosomal site 1/2, so
    r = 2(2 + a) / (3(1 + a)) and a = (4 - 3r) / (3r - 2).
    r = 1 gives a = 1 (no male excess); r falls towards 2/3 as a grows without
    bound. r above 4/3 has no non-negative solution.
    """
    if r is None or not np.isfinite(r):
        return np.nan
    if r <= 2.0 / 3.0:
        return np.inf
    a = (4.0 - 3.0 * r) / (3.0 * r - 2.0)
    return a if a >= 0 else np.nan


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def gene_level(x_vals, a_vals):
    """Mann-Whitney X vs autosomes on genes; medians and the X/A ratio."""
    out = {'n_x': len(x_vals), 'n_a': len(a_vals),
           'median_x': np.median(x_vals) if len(x_vals) else np.nan,
           'median_a': np.median(a_vals) if len(a_vals) else np.nan,
           'p': np.nan}
    out['ratio_x_a'] = (out['median_x'] / out['median_a']
                        if out['median_a'] else np.nan)
    if len(x_vals) >= 5 and len(a_vals) >= 5:
        out['p'] = mannwhitneyu(x_vals, a_vals, alternative='two-sided')[1]
    return out


def chromosome_permutation(medians, is_x, seed=0, n_mc=200_000):
    """Exact (or Monte Carlo) permutation test on per-chromosome medians.

    Statistic: mean median of the chromosomes labelled X minus mean median of
    the rest. Every way of choosing which k of the N chromosomes are X is
    enumerated, so the P is exact and needs no distributional assumption.

    Returns dict: stat, p_two, p_x_lower, n_assignments, exact, n_x, n_a.
    """
    m = np.asarray(medians, dtype=float)
    is_x = np.asarray(is_x, dtype=bool)
    n, k = len(m), int(is_x.sum())
    res = {'n_x': k, 'n_a': n - k, 'stat': np.nan, 'p_two': np.nan,
           'p_x_lower': np.nan, 'n_assignments': 0, 'exact': True}
    if k < 2 or n - k < 2:
        return res
    total = m.sum()
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
    stats = sel / k - (total - sel) / (n - k)
    tol = 1e-12 * max(1.0, abs(obs))
    res.update(stat=obs, n_assignments=n_assign,
               p_two=float(np.mean(np.abs(stats) >= abs(obs) - tol)),
               p_x_lower=float(np.mean(stats <= obs + tol)))
    return res


def ratio_bootstrap(groups_x, groups_a, n_boot=2000, seed=0):
    """Percentile CI for median(X genes) / median(autosomal genes), resampling
    whole chromosomes with replacement within each class.

    groups_x, groups_a : lists of per-chromosome value arrays.
    """
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    for b in range(n_boot):
        xs = np.concatenate([groups_x[i] for i in
                             rng.integers(0, len(groups_x), len(groups_x))])
        as_ = np.concatenate([groups_a[i] for i in
                              rng.integers(0, len(groups_a), len(groups_a))])
        ma = np.median(as_)
        out[b] = np.median(xs) / ma if ma else np.nan
    out = out[np.isfinite(out)]
    if not len(out):
        return np.nan, np.nan
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def region_contrast(df, mask, metric):
    """Eliminated vs retained within `mask`: n, medians, fold, Mann-Whitney P."""
    sub = df[mask]
    r = sub.loc[sub['region'] == 'retained', metric].dropna().to_numpy(float)
    e = sub.loc[sub['region'] == 'eliminated', metric].dropna().to_numpy(float)
    out = {'n_ret': len(r), 'n_elim': len(e),
           'median_ret': np.median(r) if len(r) else np.nan,
           'median_elim': np.median(e) if len(e) else np.nan, 'p': np.nan}
    out['fold'] = (out['median_elim'] / out['median_ret']
                   if out['median_ret'] else np.nan)
    if len(r) >= 5 and len(e) >= 5:
        out['p'] = mannwhitneyu(r, e, alternative='two-sided')[1]
    return out


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def per_chromosome_table(df, min_genes):
    """Per-chromosome counts and medians of retained genes, X and autosomes."""
    rows = []
    for chrom, sub in df[df['chrom_class'].isin(['X', 'A'])].groupby('chrom'):
        ret = sub[sub['region'] == 'retained']
        row = {'chrom': chrom, 'class': sub['chrom_class'].iloc[0],
               'n_retained': len(ret),
               'n_eliminated': int((sub['region'] == 'eliminated').sum())}
        for m in METRICS:
            v = ret[m].dropna()
            row[f'n_{m}'] = len(v)
            row[f'median_{m}'] = v.median() if len(v) else np.nan
        row['used'] = len(ret) >= min_genes
        rows.append(row)
    tab = pd.DataFrame(rows)
    if len(tab):
        key = tab['chrom'].map(lambda c: (int(re.sub(r'\D', '', c) or 10 ** 6)))
        tab = (tab.assign(_x=tab['class'].eq('X'), _k=key)
                  .sort_values(['_x', '_k', 'chrom'])
                  .drop(columns=['_x', '_k']).reset_index(drop=True))
    return tab


def run(df, args):
    """All tests. Returns (tests DataFrame, per-chromosome table, notes dict)."""
    tests, notes = [], {}
    ret = df[(df['region'] == 'retained') & df['chrom_class'].isin(['X', 'A'])]
    chrom_tab = per_chromosome_table(df, args.min_genes)
    used = chrom_tab[chrom_tab['used']]
    dropped = chrom_tab[~chrom_tab['used']]
    if len(dropped):
        logging.warning('  chromosomes with < %d retained genes, left out of '
                        'the chromosome-level tests: %s', args.min_genes,
                        ', '.join(dropped['chrom']))
    logging.info('  chromosome-level units: %d X, %d autosomes',
                 int(used['class'].eq('X').sum()),
                 int(used['class'].eq('A').sum()))

    # --- 1. X vs autosomes, retained genes ----------------------------------
    for m in METRICS:
        xv = ret.loc[ret['chrom_class'] == 'X', m].dropna().to_numpy(float)
        av = ret.loc[ret['chrom_class'] == 'A', m].dropna().to_numpy(float)
        g = gene_level(xv, av)
        tests.append({'family': 'x_vs_a_gene', 'metric': m, 'subset': 'retained',
                      'group1': 'autosome', 'group2': 'X',
                      'n1': g['n_a'], 'n2': g['n_x'],
                      'median1': g['median_a'], 'median2': g['median_x'],
                      'ratio': g['ratio_x_a'], 'p': g['p'],
                      'p_one_sided': np.nan, 'units': 'genes'})

        cm = used[f'median_{m}'].to_numpy(float)
        ok = np.isfinite(cm)
        perm = chromosome_permutation(cm[ok], used['class'].eq('X').to_numpy()[ok],
                                      seed=args.seed)
        tests.append({'family': 'x_vs_a_chromosome', 'metric': m,
                      'subset': 'retained', 'group1': 'autosome', 'group2': 'X',
                      'n1': perm['n_a'], 'n2': perm['n_x'],
                      'median1': float(np.mean(cm[ok][~used['class'].eq('X')
                                                      .to_numpy()[ok]])),
                      'median2': float(np.mean(cm[ok][used['class'].eq('X')
                                                      .to_numpy()[ok]])),
                      'ratio': np.nan, 'p': perm['p_two'],
                      'p_one_sided': perm['p_x_lower'],
                      'units': (f"chromosomes ({perm['n_assignments']:,} "
                                f"{'exact' if perm['exact'] else 'Monte Carlo'} "
                                'assignments; medians are means of '
                                'per-chromosome medians)')})
        logging.info('  %-6s X/A %.3f  gene P %s  chromosome P %s (X lower: %s)',
                     METRIC_LABEL[m], g['ratio_x_a'], p_text(g['p']),
                     p_text(perm['p_two']), p_text(perm['p_x_lower']))

    # --- male mutation bias from dS -----------------------------------------
    groups = {cls: [ret.loc[(ret['chrom'] == c), 'dS'].dropna().to_numpy(float)
                    for c in used.loc[used['class'] == cls, 'chrom']]
              for cls in ('X', 'A')}
    groups = {k: [v for v in vs if len(v)] for k, vs in groups.items()}
    r = next(t['ratio'] for t in tests
             if t['family'] == 'x_vs_a_gene' and t['metric'] == 'dS')
    if len(groups['X']) >= 2 and len(groups['A']) >= 2:
        lo, hi = ratio_bootstrap(groups['X'], groups['A'], args.n_boot, args.seed)
    else:
        lo = hi = np.nan
    notes['dS_ratio'] = (r, lo, hi)
    notes['alpha'] = (male_bias_alpha(r), male_bias_alpha(hi),
                      male_bias_alpha(lo))
    logging.info('  dS X/A %.3f (chromosome bootstrap 95%% CI %.3f-%.3f) -> '
                 'alpha %.2f (%.2f-%.2f)', r, lo, hi, *notes['alpha'])

    # --- 2. eliminated vs retained within autosomes and within the X --------
    strata = (('all pairs', np.ones(len(df), dtype=bool)),
              ('autosomes', df['chrom_class'].eq('A').to_numpy()),
              ('X', df['chrom_class'].eq('X').to_numpy()))
    for m in METRICS:
        for name, mask in strata:
            c = region_contrast(df, mask, m)
            tests.append({'family': 'region_within_class', 'metric': m,
                          'subset': name, 'group1': 'retained',
                          'group2': 'eliminated', 'n1': c['n_ret'],
                          'n2': c['n_elim'], 'median1': c['median_ret'],
                          'median2': c['median_elim'], 'ratio': c['fold'],
                          'p': c['p'], 'p_one_sided': np.nan, 'units': 'genes'})
            if m == 'dN_dS':
                logging.info('  dN/dS eliminated/retained, %-9s %.2f-fold '
                             '(n = %s / %s), P %s', name, c['fold'],
                             f"{c['n_elim']:,}", f"{c['n_ret']:,}",
                             p_text(c['p']))

    # --- composition ----------------------------------------------------------
    known = df[df['region'].isin(['retained', 'eliminated'])
               & df['chrom_class'].isin(['X', 'A'])]
    tab = pd.crosstab(known['region'], known['chrom_class']).reindex(
        index=['retained', 'eliminated'], columns=['A', 'X'], fill_value=0)
    odds, p_comp = fisher_exact(tab.to_numpy())
    notes['composition'] = (tab, odds, p_comp)
    logging.info('  X-linked share: retained %.1f%%, eliminated %.1f%% '
                 '(Fisher P %s)', 100 * tab.loc['retained', 'X'] /
                 max(tab.loc['retained'].sum(), 1),
                 100 * tab.loc['eliminated', 'X'] /
                 max(tab.loc['eliminated'].sum(), 1), p_text(p_comp))

    out = pd.DataFrame(tests)
    out['q'] = np.nan
    for (fam, met), block in out.groupby(['family', 'metric']):
        ok = block['p'].notna()
        if ok.any():
            out.loc[block.index[ok], 'q'] = benjamini_hochberg(
                block.loc[ok, 'p'].to_numpy(float))
    return out, chrom_tab, notes


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def write_summary(path, tests, chrom_tab, notes, args, n_rows, class_counts):
    """Human-readable results, ending with what the text can claim."""
    def row(fam, metric, subset=None):
        t = tests[(tests['family'] == fam) & (tests['metric'] == metric)]
        if subset is not None:
            t = t[t['subset'] == subset]
        return t.iloc[0] if len(t) else None

    L = ['X chromosomes vs autosomes (Ascaris), and the main contrast within '
         'each', '=' * 72, '',
         f'Input              {args.annotated}',
         f'Pairs              {n_rows:,}',
         f'Chromosome classes X {class_counts.get("X", 0):,}   autosome '
         f'{class_counts.get("A", 0):,}   other {class_counts.get("other", 0):,}'
         '  (other = unplaced; excluded from X/A tests)',
         f'X pattern          {args.x_pattern}    autosome pattern '
         f'{args.autosome_pattern}',
         f'Chromosome units   >= {args.min_genes} retained genes', '']

    L += ['-' * 72, '1. X vs autosomes, RETAINED genes only', '-' * 72,
          f"{'metric':<7} {'median A':>9} {'median X':>9} {'X/A':>6} "
          f"{'gene P':>9} {'chrom P':>9} {'X lower':>9}"]
    for m in METRICS:
        g, c = row('x_vs_a_gene', m), row('x_vs_a_chromosome', m)
        L.append(f"{METRIC_LABEL[m]:<7} {g['median1']:>9.4f} {g['median2']:>9.4f} "
                 f"{g['ratio']:>6.3f} {p_text(g['p']):>9} {p_text(c['p']):>9} "
                 f"{p_text(c['p_one_sided']):>9}")
    c0 = row('x_vs_a_chromosome', 'dS')
    L += ['', 'gene P: Mann-Whitney on genes. chrom P: permutation of the X label '
          f"over {c0['units']},", "two-sided; 'X lower' is one-sided.",
          'Genes on one chromosome are not independent: quote the chromosome '
          'P.', '']

    r, lo, hi = notes['dS_ratio']
    a, a_lo, a_hi = notes['alpha']
    L += ['Male-driven evolution (assumes every X is carried once in males and '
          'twice', 'in females, and that the whole dS gap is mutational):',
          f'  X/A median dS r = {r:.3f}  (95% CI {lo:.3f}-{hi:.3f}, chromosomes '
          'resampled within class)',
          f'  alpha = (4 - 3r) / (3r - 2) = {a:.2f}  (95% CI {a_lo:.2f}-'
          f'{a_hi:.2f})', '']

    L += ['-' * 72, '2. Eliminated vs retained, within each chromosome class',
          '-' * 72,
          f"{'metric':<7} {'subset':<10} {'n elim':>7} {'n ret':>7} "
          f"{'med elim':>9} {'med ret':>9} {'fold':>6} {'P':>9}"]
    for m in METRICS:
        for s in ('all pairs', 'autosomes', 'X'):
            t = row('region_within_class', m, s)
            L.append(f"{METRIC_LABEL[m]:<7} {s:<10} {int(t['n2']):>7,} "
                     f"{int(t['n1']):>7,} {t['median2']:>9.4f} "
                     f"{t['median1']:>9.4f} {t['ratio']:>6.2f} "
                     f"{p_text(t['p']):>9}")
    tab, odds, p_comp = notes['composition']
    sh_r = tab.loc['retained', 'X'] / max(tab.loc['retained'].sum(), 1)
    sh_e = tab.loc['eliminated', 'X'] / max(tab.loc['eliminated'].sum(), 1)
    L += ['', "'all pairs' is the genome-wide contrast as Results quotes it "
          '(1.77-fold,', 'P = 2.06e-40 for dN/dS); it includes unplaced '
          'scaffolds.', '',
          f'X-linked share: retained {100 * sh_r:.1f}% '
          f"({tab.loc['retained', 'X']:,}/{tab.loc['retained'].sum():,}), "
          f'eliminated {100 * sh_e:.1f}% '
          f"({tab.loc['eliminated', 'X']:,}/{tab.loc['eliminated'].sum():,}); "
          f'Fisher P {p_text(p_comp)}.', '']

    # --- what the text can claim ------------------------------------------
    L += ['-' * 72, 'What the text can claim (chromosome-level P < 0.05 for the '
          'X/A rows)', '-' * 72]
    for m in ('dS', 'dN_dS'):
        g, c = row('x_vs_a_gene', m), row('x_vs_a_chromosome', m)
        direction = 'lower' if g['ratio'] < 1 else 'higher'
        if np.isfinite(c['p']) and c['p'] < 0.05:
            L.append(f"  {METRIC_LABEL[m]}: X-linked genes {direction} "
                     f"(X/A {g['ratio']:.2f}; chromosome-level P = "
                     f"{p_text(c['p'])}). Can be stated.")
        else:
            L.append(f"  {METRIC_LABEL[m]}: X/A {g['ratio']:.2f}, chromosome-"
                     f"level P = {p_text(c['p'])}. Not supported: do not "
                     'claim a difference.')
    a_txt = (f'alpha ~ {a:.1f} (95% CI {a_lo:.1f}-{a_hi:.1f})'
             if np.isfinite(a) else 'alpha not defined for this r')
    L.append(f'  If the dS gap is stated: male-biased mutation as a hypothesis, '
             f'{a_txt}.')
    auto = row('region_within_class', 'dN_dS', 'autosomes')
    allp = row('region_within_class', 'dN_dS', 'all pairs')
    L.append(f"  Main contrast on autosomes only: {auto['ratio']:.2f}-fold "
             f"(P = {p_text(auto['p'])}), against {allp['ratio']:.2f}-fold on "
             'all pairs.')
    L.append('')

    L += ['-' * 72, 'Per-chromosome medians, retained genes (per_chromosome.tsv)',
          '-' * 72,
          chrom_tab.to_string(index=False, float_format=lambda v: f'{v:.4f}'),
          '']
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Ascaris X vs autosomes, and eliminated vs retained within '
                    'each chromosome class.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', required=True,
                   help='kaks_annotated.tsv from annotate_regions.py '
                        '(ANNOTATED_TABLE in config.sh)')
    p.add_argument('--chrom-col', default='chrom1',
                   help='Ascaris chromosome column (default: chrom1)')
    p.add_argument('--region-col', default='region1',
                   help='Ascaris region column (default: region1)')
    p.add_argument('--x-pattern', default=r'^chrX',
                   help=r"Regex for X chromosomes (default: '^chrX')")
    p.add_argument('--autosome-pattern', default=r'^chr\d+$',
                   help=r"Regex for autosomes (default: '^chr\d+$')")
    p.add_argument('--min-genes', type=int, default=30,
                   help='Retained genes a chromosome needs to be a unit of the '
                        'chromosome-level test (default: 30)')
    p.add_argument('--n-boot', type=int, default=2000,
                   help='Chromosome bootstrap replicates (default: 2000)')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.min_genes < 1 or args.n_boot < 100:
        p.error('--min-genes must be >= 1 and --n-boot >= 100')
    if args.outdir is None:
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.outdir = os.path.join('results', f'x_autosome_check_{stamp}')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'x_autosome_check')
    logging.info('X vs autosomes, and the main contrast within each class')
    logging.info('  input:  %s', os.path.abspath(args.annotated))
    logging.info('  outdir: %s', os.path.abspath(args.outdir))

    if not os.path.isfile(args.annotated):
        sys.exit(f'ERROR: annotated table not found: {args.annotated}')
    df = read_rate_table(args.annotated, require=METRICS)
    for col in (args.chrom_col, args.region_col):
        if col not in df.columns:
            sys.exit(f'ERROR: column {col!r} not in {args.annotated}. '
                     f'Columns: {list(df.columns)}')
    df = df.rename(columns={args.chrom_col: 'chrom', args.region_col: 'region'})
    for m in METRICS:
        df[m] = pd.to_numeric(df[m], errors='coerce')
    df['chrom_class'] = classify_chromosomes(df['chrom'], args.x_pattern,
                                             args.autosome_pattern)

    class_counts = df['chrom_class'].value_counts().to_dict()
    logging.info('Loaded %s pairs: X %s, autosome %s, other %s', f'{len(df):,}',
                 f"{class_counts.get('X', 0):,}", f"{class_counts.get('A', 0):,}",
                 f"{class_counts.get('other', 0):,}")
    x_chroms = sorted(df.loc[df['chrom_class'] == 'X', 'chrom'].unique())
    logging.info('  X chromosomes: %s', ', '.join(map(str, x_chroms)) or 'none')
    if not x_chroms:
        sys.exit(f'ERROR: no chromosome matched --x-pattern {args.x_pattern!r}. '
                 f"Names seen: {sorted(map(str, df['chrom'].dropna().unique()))[:12]}")
    other = sorted(map(str, df.loc[df['chrom_class'] == 'other', 'chrom']
                       .dropna().unique()))
    if other:
        logging.info('  other (excluded from X/A tests): %d sequences, e.g. %s',
                     len(other), ', '.join(other[:5]))
    logging.info('  regions: %s', df['region'].value_counts().to_dict())

    tests, chrom_tab, notes = run(df, args)
    tests.to_csv(os.path.join(args.outdir, 'x_autosome_tests.tsv'), sep='\t',
                 index=False, float_format='%.6g')
    chrom_tab.to_csv(os.path.join(args.outdir, 'per_chromosome.tsv'), sep='\t',
                     index=False, float_format='%.6g')
    write_summary(os.path.join(args.outdir, 'x_autosome_summary.txt'), tests,
                  chrom_tab, notes, args, len(df), class_counts)

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(),
                   'python': platform.python_version(),
                   'argv': sys.argv, 'n_pairs': int(len(df)),
                   'inputs': {'annotated': fingerprint(args.annotated)}},
                  fh, indent=2)
    logging.info('Summary: %s', os.path.join(args.outdir,
                                             'x_autosome_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
