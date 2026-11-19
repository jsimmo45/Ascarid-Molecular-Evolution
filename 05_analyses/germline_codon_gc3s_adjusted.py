#!/usr/bin/env python3
"""
germline_codon_gc3s_adjusted.py - the codon-bias contrasts of germline-restricted
genes with base composition (GC3s) controlled using every gene.

Why
---
germline_restricted_codon.py (run 2026-10-29) found that eliminated and
retained germline-restricted genes carry the same codon bias (dENC 2.20 vs
2.14, P = 0.86) and that germline restriction alone lowers it among retained
genes (2.14 vs 2.69, P = 9.8e-6). Both comparisons are confounded by base
composition. Eliminated germline-restricted genes are far richer in GC3s than
retained ones (median 0.475 vs 0.397), and retained germline-restricted genes
are poorer in it than somatic-peak genes (0.397 vs 0.427). dENC removes
composition only through Wright's curve, which is not a complete correction
in this genome: genome-wide, matching on GC3s enlarged the eliminated deficit
from -0.43 to -0.72 (codon_usage_bias.py). The GC3s-matched elimination
contrast of that run kept 138 of 213 genes and was too imprecise to settle
the question (+0.09, 95% CI -0.71 to +0.50), and the restriction contrast was
not matched on GC3s at all.

What the script does
--------------------
It reads the joined per-gene table of that run (germline_codon_genes.tsv:
the group of every gene at every definition of germline restriction, its
stage of peak expression, log2 expression, dENC, GC3s and the ENC-ceiling
flag), so the genes, groups and metric are exactly those already reported.

1. Regression on every expressed gene (groups R_S, R_Gbroad, R_GR, E_GR,
   E_other, as S12B-C), with the group as a factor (reference R_GR) and GC3s
   as a restricted cubic spline (Harrell's knots), and in the fuller model
   log2 expression as a second spline. The shape of the GC3s dependence is
   estimated from every gene rather than from matched pairs, so no gene is
   dropped. The two contrasts:
       elimination  E_GR - R_GR    (the E_GR coefficient)
       restriction  R_GR - R_S     (minus the R_S coefficient)
   Models: group only (the unadjusted means, to show what the adjustment does
   within one estimator), + GC3s, + GC3s + expression. OLS with HC3 standard
   errors (Wald P) and a stratified bootstrap interval (genes resampled within
   each group), at every definition. At the default definition also: median
   regression (the earlier contrasts were differences of medians; IRLS,
   checked against an exact linear-programming fit), without genes at the
   ENC ceiling, with 7 GC3s knots instead of --knots, and a test of whether
   the GC3s slope differs for R_GR and E_GR (if it does, one adjusted
   difference averages over composition).
2. The same models genome-wide (every gene with an expression value,
   eliminated vs retained): the size against which the elimination contrast
   is judged.
3. Matched contrasts, as differences of medians, using the matching of
   germline_restricted_codon.matched_difference (imported): restriction (each
   R_GR gene paired with an unused R_S gene) and elimination (each E_GR gene
   with an unused R_GR gene), on GC3s alone within --gc3s-caliper (default
   0.01, the caliper at which codon_usage_bias.py calibrated its matched
   test) and on GC3s plus expression level within --caliper. No stage
   matching: R_GR and R_S differ in peak stage by definition.
4. Descriptive: n and median dENC of R_S, R_GR and E_GR in fixed GC3s bins at
   the default definition (no test: codon_usage_bias.py found stratified
   tests miscalibrated with strata wider than 0.01).

Reading
-------
  elimination interval below 0: the conclusion changes. At equal composition,
      eliminated genes carry weaker codon bias than equally restricted
      retained genes, so codon usage, like dN/dS, separates them.
  interval spanning 0 but excluding the genome-wide adjusted difference: the
      conclusion holds and is firmer (no extra deficit, and one as large as
      the genome-wide difference is excluded).
  interval spanning both: unresolved.
  restriction interval below 0: germline restriction lowers codon bias at
      equal composition; otherwise the restriction contrast may be
      compositional.

Writes, in --outdir
-------------------
  gc3s_models.tsv     every model x contrast x definition: estimate, SE, P,
                      bootstrap interval, n
  gc3s_matched.tsv    the matched contrasts
  gc3s_strata.tsv     median dENC per group in GC3s bins (default definition)
  gc3s_adjusted_summary.txt   for reading, ending with what to say
  germline_codon_gc3s_adjusted.log
  provenance.json

Usage
-----
    python germline_codon_gc3s_adjusted.py \\
      --genes results/germline_restricted_codon_20261029_143308/germline_codon_genes.tsv \\
      --outdir results/germline_codon_gc3s_adjusted_20261029_170000

Dependencies: numpy, pandas, scipy, germline_restricted_codon (and what it
imports: germline_restricted_dnds, functional_enrichment, paper_colors,
region_lib)
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
from scipy import sparse
from scipy.optimize import linprog
from scipy.stats import chi2, norm
from scipy import __version__ as scipy_version

# One definition of the groups, the joined table, the matching and the text
# formatting with germline_restricted_codon.py (and through it S12B-C).
from germline_restricted_codon import (definitions, matched_difference,  # noqa: E402
                                       ptxt, ci_txt)
from germline_restricted_dnds import GROUP_MEMBERS  # noqa: E402
from functional_enrichment import PEAK_COL  # noqa: E402
from region_lib import fingerprint, git_commit, setup_logging  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Groups of the model, as S12B-C; R_GR is the reference level.
MODEL_GROUPS = ['R_S', 'R_Gbroad', 'R_GR', 'E_GR', 'E_other']
REFERENCE = 'R_GR'
#: Harrell's knot quantiles for a restricted cubic spline with k knots.
HARRELL_Q = {3: [0.10, 0.50, 0.90],
             4: [0.05, 0.35, 0.65, 0.95],
             5: [0.05, 0.275, 0.50, 0.725, 0.95],
             6: [0.05, 0.23, 0.41, 0.59, 0.77, 0.95],
             7: [0.025, 0.1833, 0.3417, 0.50, 0.6583, 0.8167, 0.975]}
#: Fixed GC3s bins of the descriptive table.
STRATA_EDGES = np.round(np.arange(0.20, 0.7001, 0.05), 2)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def load_genes(path, metric):
    """The joined per-gene table of germline_restricted_codon.py, checked."""
    df = pd.read_csv(path, sep='\t', dtype={'gene1': str})
    need = ['gene1', 'region1', PEAK_COL, 'log_expr', 'group_default', metric,
            'GC3s', 'ENC_at_ceiling']
    miss = [c for c in need if c not in df.columns]
    if miss:
        sys.exit(f'ERROR: {path} lacks {miss}; is it germline_codon_genes.tsv '
                 f'from germline_restricted_codon.py?\n  columns: {list(df.columns)}')
    if df['gene1'].duplicated().any():
        sys.exit(f'ERROR: duplicated gene IDs in {path}.')
    bad = ~df['region1'].isin(['retained', 'eliminated'])
    if bad.any():
        sys.exit(f'ERROR: {int(bad.sum())} genes with a region other than '
                 f'retained/eliminated in {path}.')
    df['ENC_at_ceiling'] = df['ENC_at_ceiling'].astype(str).str.lower().isin(
        ['true', '1', '1.0'])
    for c in (metric, 'GC3s', 'log_expr'):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    logging.info('Genes: %s (%s retained, %s eliminated) from %s',
                 f'{len(df):,}', f"{int((df['region1'] == 'retained').sum()):,}",
                 f"{int((df['region1'] == 'eliminated').sum()):,}",
                 os.path.basename(path))
    return df


def model_frame(df, labels, metric, drop_ceiling=False):
    """Genes of the five model groups with the metric, GC3s and expression.

    Returns (frame with y, gc3s, x, group; accounting dict). The accounting
    names how many genes each filter removes, so a reader can see where the
    expressed genes went.
    """
    fine_to_group = {f: g for g in MODEL_GROUPS for f in GROUP_MEMBERS[g]}
    grp = labels.map(fine_to_group)
    in_group = grp.notna()
    have = df[metric].notna() & df['GC3s'].notna() & df['log_expr'].notna()
    keep = in_group & have
    acct = {'in_model_groups': int(in_group.sum()),
            'missing_value': int((in_group & ~have).sum())}
    if drop_ceiling:
        acct['at_ceiling_dropped'] = int((keep & df['ENC_at_ceiling']).sum())
        keep &= ~df['ENC_at_ceiling']
    f = pd.DataFrame({'y': df.loc[keep, metric].to_numpy(float),
                      'gc3s': df.loc[keep, 'GC3s'].to_numpy(float),
                      'x': df.loc[keep, 'log_expr'].to_numpy(float),
                      'stage': df.loc[keep, PEAK_COL].astype(str).to_numpy(object),
                      'group': grp[keep].to_numpy(object)})
    acct['used'] = len(f)
    acct['per_group'] = {g: int((f['group'] == g).sum()) for g in MODEL_GROUPS}
    return f, acct


# ---------------------------------------------------------------------------
# Splines and fits
# ---------------------------------------------------------------------------

def rcs_knots(x, k):
    """Harrell's knot positions (quantiles) for a k-knot restricted spline."""
    if k not in HARRELL_Q:
        raise ValueError(f'knots must be one of {sorted(HARRELL_Q)}')
    t = np.unique(np.quantile(np.asarray(x, float), HARRELL_Q[k]))
    if len(t) < k:
        raise ValueError(f'only {len(t)} distinct knot positions for {k} knots')
    return t


def rcs_basis(x, knots):
    """Restricted (natural) cubic spline basis in Harrell's parameterization.

    Columns: x, then k - 2 nonlinear terms; the fitted curve is cubic between
    the knots and linear beyond the outer two. The nonlinear terms are divided
    by (t_k - t_1)^2 so they stay of the order of x and the design is well
    conditioned.
    """
    x = np.asarray(x, float)
    t = np.asarray(knots, float)
    k = len(t)
    norm2 = (t[-1] - t[0]) ** 2
    cols = [x]
    for j in range(k - 2):
        a = np.clip(x - t[j], 0, None) ** 3
        b = np.clip(x - t[k - 2], 0, None) ** 3 * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
        c = np.clip(x - t[k - 1], 0, None) ** 3 * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2])
        cols.append((a - b + c) / norm2)
    return np.column_stack(cols)


def design(frame, spec, knots):
    """Design matrix and column names for one model.

    spec: {'gc3s', 'expr', 'stage', 'interact': bool}. knots: arrays for
    'gc3s' and 'expr' (fixed from the full sample, so bootstrap fits use the
    same basis) and 'centre' (the GC3s at which interacting groups'
    contrasts are evaluated). Group dummies only for groups present.
    """
    g = frame['group'].to_numpy(object)
    cols, names = [np.ones(len(frame))], ['intercept']
    for grp in MODEL_GROUPS:
        if grp != REFERENCE and (g == grp).any():
            cols.append((g == grp).astype(float))
            names.append(f'group[{grp}]')
    if spec.get('gc3s'):
        B = rcs_basis(frame['gc3s'], knots['gc3s'])
        cols += list(B.T)
        names += [f'gc3s_{j}' for j in range(B.shape[1])]
    if spec.get('expr'):
        B = rcs_basis(frame['x'], knots['expr'])
        cols += list(B.T)
        names += [f'expr_{j}' for j in range(B.shape[1])]
    if spec.get('stage'):
        # Peak-stage dummies, the most common stage as reference. Only for
        # germline-peak genes: among them stage is not tied to group.
        st = frame['stage'].to_numpy(object)
        levels = pd.Series(st).value_counts().index.tolist()
        for lev in levels[1:]:
            cols.append((st == lev).astype(float))
            names.append(f'stage[{lev}]')
    if spec.get('interact'):
        gc = frame['gc3s'].to_numpy(float) - knots['centre']
        for grp in ('R_GR', 'E_GR'):
            cols.append((g == grp).astype(float) * gc)
            names.append(f'gc3s_x_{grp}')
    return np.column_stack(cols), names


def contrast_vectors(names):
    """Coefficient vectors of the two contrasts (None if a group is absent)."""
    out = {}
    if 'group[E_GR]' in names:
        c = np.zeros(len(names))
        c[names.index('group[E_GR]')] = 1.0
        out['elimination'] = c
    if 'group[R_S]' in names:
        c = np.zeros(len(names))
        c[names.index('group[R_S]')] = -1.0      # R_GR - R_S
        out['restriction'] = c
    return out


def ols_hc3(X, y):
    """OLS with HC3 covariance (MacKinnon & White 1985), as
    toxocara_annotation_diagnostic.ols_hc3. Raises on a rank-deficient design."""
    n, k = X.shape
    if np.linalg.matrix_rank(X) < k:
        raise np.linalg.LinAlgError('rank-deficient design (an empty group or a '
                                    'constant covariate)')
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    resid = y - X @ beta
    h = np.clip(np.einsum('ij,jk,ik->i', X, XtX_inv, X), 0.0, 1 - 1e-8)
    omega = (resid / (1.0 - h)) ** 2
    cov = XtX_inv @ (X.T * omega) @ X @ XtX_inv
    return beta, cov


def ols_beta(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]


def median_irls(X, y, beta0=None, max_iter=1000, tol=1e-9):
    """Median (least absolute deviation) regression by iteratively reweighted
    least squares. Weights 1/max(|r|, eps) with eps a millionth of the
    response's MAD. Returns (beta, iterations, converged)."""
    beta = ols_beta(X, y) if beta0 is None else np.array(beta0, float)
    mad = float(np.median(np.abs(y - np.median(y)))) or 1.0
    eps = 1e-6 * mad
    for it in range(1, max_iter + 1):
        r = y - X @ beta
        sw = 1.0 / np.sqrt(np.maximum(np.abs(r), eps))
        new = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        if np.max(np.abs(new - beta)) < tol * (1.0 + np.max(np.abs(beta))):
            return new, it, True
        beta = new
    return beta, max_iter, False


def median_lp(X, y):
    """Exact median regression as a linear programme (HiGHS); a check on
    median_irls, run once per model at the full sample."""
    n, k = X.shape
    c = np.concatenate([np.zeros(k), 0.5 * np.ones(n), 0.5 * np.ones(n)])
    A = sparse.hstack([sparse.csr_matrix(X), sparse.identity(n, format='csr'),
                       -sparse.identity(n, format='csr')], format='csr')
    res = linprog(c, A_eq=A, b_eq=y, bounds=[(None, None)] * k + [(0, None)] * (2 * n),
                  method='highs')
    if res.status != 0:
        raise RuntimeError(f'linear programme failed: {res.message}')
    return res.x[:k]


def strat_index(groups, rng):
    """Bootstrap row indices, resampled within each group (sizes preserved)."""
    parts = []
    for g in pd.unique(groups):
        ix = np.flatnonzero(groups == g)
        parts.append(ix[rng.integers(0, len(ix), len(ix))])
    return np.concatenate(parts)


# ---------------------------------------------------------------------------
# One model, fitted and bootstrapped
# ---------------------------------------------------------------------------

def fit_model(frame, spec, knots, n_boot, rng, median=False, lp_check=False,
              level=0.95):
    """Fit one model and bootstrap it. Returns {contrast: row dict}, info.

    OLS: estimate, HC3 SE, Wald z and P, bootstrap percentile interval.
    Median regression: IRLS estimate (the LP estimate if lp_check and the two
    disagree), bootstrap SE, z, P and percentile interval.
    """
    X, names = design(frame, spec, knots)
    y = frame['y'].to_numpy(float)
    cvec = contrast_vectors(names)
    info = {'n': len(y), 'k': X.shape[1], 'names': names}
    if not cvec:
        return {}, info
    if median:
        beta, iters, conv = median_irls(X, y)
        info.update(irls_iterations=iters, irls_converged=conv)
        if not conv:
            logging.warning('    IRLS did not converge in %d iterations', iters)
        if lp_check:
            b_lp = median_lp(X, y)
            gap = max(abs(float(c @ (beta - b_lp))) for c in cvec.values())
            info['irls_vs_lp_max_gap'] = gap
            logging.info('    median regression: IRLS vs exact LP, largest '
                         'contrast gap %.2e', gap)
            if gap > 0.01:
                logging.warning('    IRLS and LP disagree by %.3f; using the LP '
                                'estimate (intervals still from IRLS fits)', gap)
                beta = b_lp
        cov = None
    else:
        beta, cov = ols_hc3(X, y)
    groups = frame['group'].to_numpy(object)
    boots = {k: np.empty(n_boot) for k in cvec}
    for b in range(n_boot):
        ix = strat_index(groups, rng)
        if median:
            bb, _, _ = median_irls(X[ix], y[ix], beta0=beta, max_iter=200, tol=1e-6)
        else:
            bb = ols_beta(X[ix], y[ix])
        for k, c in cvec.items():
            boots[k][b] = c @ bb
    alpha = (1 - level) / 2
    out = {}
    for k, c in cvec.items():
        est = float(c @ beta)
        lo, hi = (float(v) for v in np.quantile(boots[k], [alpha, 1 - alpha]))
        if median:
            se = float(np.std(boots[k], ddof=1))
            se_kind = 'bootstrap'
        else:
            se = float(np.sqrt(max(c @ cov @ c, 0.0)))
            se_kind = 'HC3'
        z = est / se if se > 0 else np.nan
        out[k] = {'estimate': est, 'se': se, 'se_kind': se_kind, 'z': z,
                  'p': float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan,
                  'ci_lo': lo, 'ci_hi': hi, 'n_boot': n_boot}
    info['beta'] = beta
    info['cov'] = cov
    return out, info


def interaction_test(frame, knots):
    """Do R_GR and E_GR have their own GC3s slopes? HC3 Wald test of the two
    group x (GC3s - centre) terms added to the GC3s + expression model."""
    spec = {'gc3s': True, 'expr': True, 'interact': True}
    X, names = design(frame, spec, knots)
    beta, cov = ols_hc3(X, frame['y'].to_numpy(float))
    idx = [names.index(n) for n in ('gc3s_x_R_GR', 'gc3s_x_E_GR') if n in names]
    if not idx:
        return None
    b, V = beta[idx], cov[np.ix_(idx, idx)]
    stat = float(b @ np.linalg.pinv(V) @ b)
    return {'chi2': stat, 'df': len(idx), 'p': float(chi2.sf(stat, len(idx))),
            'slopes': {names[i]: float(beta[i]) for i in idx},
            'slope_se': {names[i]: float(np.sqrt(cov[i, i])) for i in idx}}


# ---------------------------------------------------------------------------
# Genome-wide reference and matched contrasts
# ---------------------------------------------------------------------------

def genome_wide(df, metric, spec, knots_k, expr_k, n_boot, rng, median=False):
    """Eliminated vs retained over every gene with an expression value,
    adjusted as `spec`. Returns (row dict, n)."""
    keep = df[metric].notna() & df['GC3s'].notna() & df['log_expr'].notna()
    f = pd.DataFrame({'y': df.loc[keep, metric].to_numpy(float),
                      'gc3s': df.loc[keep, 'GC3s'].to_numpy(float),
                      'x': df.loc[keep, 'log_expr'].to_numpy(float),
                      'group': df.loc[keep, 'region1'].to_numpy(object)})
    cols, names = [np.ones(len(f)), (f['group'] == 'eliminated').to_numpy(float)], \
        ['intercept', 'eliminated']
    if spec.get('gc3s'):
        B = rcs_basis(f['gc3s'], rcs_knots(f['gc3s'], knots_k))
        cols += list(B.T)
        names += [f'gc3s_{j}' for j in range(B.shape[1])]
    if spec.get('expr'):
        B = rcs_basis(f['x'], rcs_knots(f['x'], expr_k))
        cols += list(B.T)
        names += [f'expr_{j}' for j in range(B.shape[1])]
    X, y = np.column_stack(cols), f['y'].to_numpy(float)
    c = np.zeros(len(names))
    c[1] = 1.0
    if median:
        beta, _, _ = median_irls(X, y)
    else:
        beta, cov = ols_hc3(X, y)
    groups = f['group'].to_numpy(object)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        ix = strat_index(groups, rng)
        bb = (median_irls(X[ix], y[ix], beta0=beta, max_iter=200, tol=1e-6)[0]
              if median else ols_beta(X[ix], y[ix]))
        boots[b] = c @ bb
    est = float(c @ beta)
    se = float(np.std(boots, ddof=1)) if median else float(np.sqrt(c @ cov @ c))
    z = est / se if se > 0 else np.nan
    lo, hi = (float(v) for v in np.quantile(boots, [0.025, 0.975]))
    return {'estimate': est, 'se': se, 'se_kind': 'bootstrap' if median else 'HC3',
            'z': z, 'p': float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan,
            'ci_lo': lo, 'ci_hi': hi, 'n_boot': n_boot,
            'n_num': int((groups == 'eliminated').sum()),
            'n_den': int((groups == 'retained').sum())}, len(f)


def matched_rows(frame, metric_col, gc3s_caliper, caliper, n_match, n_match_boot,
                 rng):
    """Matched restriction and elimination contrasts (differences of medians).

    matched_difference() pairs each numerator gene, in random order, with an
    unused denominator gene of the same 'peak' within `caliper` on x and
    within `gc3s_caliper` on GC3s. Peak is held constant here (no stage
    matching), and caliper = inf turns expression matching off.
    """
    def pick(group):
        m = frame['group'] == group
        return pd.DataFrame({'peak': 'all', 'x': frame.loc[m, 'x'].to_numpy(float),
                             'gc3s': frame.loc[m, 'gc3s'].to_numpy(float),
                             'v': frame.loc[m, metric_col].to_numpy(float)})
    rows, reps_all = [], []
    specs = [('restriction', 'R_GR', 'R_S', np.inf, 'GC3s'),
             ('restriction', 'R_GR', 'R_S', caliper, 'GC3s + expression'),
             ('elimination', 'E_GR', 'R_GR', np.inf, 'GC3s'),
             ('elimination', 'E_GR', 'R_GR', caliper, 'GC3s + expression')]
    for contrast, num, den, cal, label in specs:
        e, r = pick(num), pick(den)
        if len(e) < 10 or len(r) < 10:
            logging.warning('  matched %s (%s): %d / %d genes; not run', contrast,
                            label, len(e), len(r))
            continue
        res, reps = matched_difference(e, r, cal, n_match, n_match_boot, rng,
                                       gc3s_caliper=gc3s_caliper)
        row = {'contrast': contrast, 'numerator': num, 'denominator': den,
               'matched_on': label, 'gc3s_caliper': gc3s_caliper,
               'expr_caliper': cal if np.isfinite(cal) else np.nan,
               'n_numerator': res['n_eliminated'], 'n_pool': res['n_retained_pool'],
               'n_no_candidate': res['n_no_candidate'],
               'median_n_matched': res['median_n_matched'], 'diff': res['diff'],
               'match_lo': res['match_lo'], 'match_hi': res['match_hi'],
               'ci_lo': res['ci_lo'], 'ci_hi': res['ci_hi'],
               'median_p': res['median_p'],
               'pct_rep_p_lt_0.05': res['pct_rep_p_lt_0.05'],
               'median_gc3s_balance_p': res['median_gc3s_balance_p'],
               'median_expr_balance_p': res['median_level_balance_p']}
        rows.append(row)
        reps.insert(0, 'matched_on', label)
        reps.insert(0, 'contrast', contrast)
        reps_all.append(reps)
        logging.info('  matched %-11s on %-17s %+.3f %s; %s of %s matched, %s '
                     'without a candidate; GC3s balance P %s', contrast, label,
                     row['diff'], ci_txt(row['ci_lo'], row['ci_hi']),
                     f"{row['median_n_matched']:.0f}", row['n_numerator'],
                     row['n_no_candidate'], ptxt(row['median_gc3s_balance_p']))
    return pd.DataFrame(rows), reps_all


def strata_table(frame):
    """n and median of y for R_S, R_GR and E_GR in fixed GC3s bins."""
    rows = []
    for lo, hi in zip(STRATA_EDGES[:-1], STRATA_EDGES[1:]):
        m = (frame['gc3s'] >= lo) & (frame['gc3s'] < hi)
        row = {'gc3s_lo': lo, 'gc3s_hi': hi}
        for g in ('R_S', 'R_GR', 'E_GR'):
            v = frame.loc[m & (frame['group'] == g), 'y']
            row[f'n_{g}'] = len(v)
            row[f'median_{g}'] = float(v.median()) if len(v) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def elimination_verdict(row, gw):
    """The reading of one adjusted elimination contrast (see the docstring).

    `gw`: the genome-wide row of the same model (estimate and interval), the
    size of the eliminated deficit against which the contrast is judged."""
    lo, hi = row['ci_lo'], row['ci_hi']
    ref = gw['estimate']
    if hi < 0:
        return ('CHANGES: at equal composition, eliminated genes carry weaker '
                'codon bias than equally restricted retained genes')
    if lo > 0:
        return ('eliminated genes carry STRONGER codon bias than equally '
                'restricted retained genes at equal composition')
    if gw['ci_hi'] >= 0:
        return ('NO ADJUSTED DEFICIT TO EXPLAIN: genome-wide, eliminated genes do '
                f'not differ at equal composition either ({ref:+.2f} '
                f"{ci_txt(gw['ci_lo'], gw['ci_hi'])})")
    if lo > ref:
        return ('HOLDS, firmer: no extra deficit, and the interval excludes the '
                f'genome-wide adjusted difference ({ref:+.2f})')
    return ('UNRESOLVED: the interval spans 0 and the genome-wide adjusted '
            f'difference ({ref:+.2f})')


def restriction_verdict(row):
    if row['ci_hi'] < 0:
        return 'germline restriction lowers codon bias at equal composition'
    if row['ci_lo'] > 0:
        return 'germline restriction RAISES codon bias at equal composition'
    return ('NOT SIGNIFICANT at equal composition: the restriction contrast may '
            'be compositional')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Codon-bias contrasts of germline-restricted genes with GC3s '
                    'controlled by regression on every gene and by matching.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--genes', required=True,
                   help='germline_codon_genes.tsv from germline_restricted_codon.py')
    p.add_argument('--metric', default='dENC', choices=['dENC', 'ENC'],
                   help='Codon metric (default dENC, as fig. 2C)')
    p.add_argument('--knots', type=int, default=5, choices=sorted(HARRELL_Q),
                   help='Knots of the GC3s spline (default 5; 7 is also fitted '
                        'at the default definition)')
    p.add_argument('--expr-knots', type=int, default=4, choices=sorted(HARRELL_Q),
                   help='Knots of the expression spline (default 4)')
    p.add_argument('--gc3s-caliper', type=float, default=0.01,
                   help='GC3s caliper of the matched contrasts (default 0.01, as '
                        'codon_usage_bias.py)')
    p.add_argument('--caliper', type=float, default=0.25,
                   help='log2 expression caliper of the GC3s + expression match '
                        '(default 0.25, as S12C)')
    p.add_argument('--n-boot', type=int, default=2000,
                   help='Bootstrap replicates of the OLS models (default 2000)')
    p.add_argument('--n-boot-median', type=int, default=1000,
                   help='Bootstrap replicates of the median regressions '
                        '(default 1000)')
    p.add_argument('--n-match', type=int, default=500)
    p.add_argument('--n-match-boot', type=int, default=1000)
    p.add_argument('--min-group', type=int, default=10,
                   help='A contrast is reported only if both of its groups hold '
                        'at least this many genes (default 10, as S12)')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', required=True)
    a = p.parse_args(argv)
    if a.gc3s_caliper <= 0 or a.caliper <= 0:
        p.error('calipers must be positive')
    if min(a.n_boot, a.n_boot_median) < 100:
        p.error('use at least 100 bootstrap replicates')
    return a


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'germline_codon_gc3s_adjusted')
    logging.info('germline_codon_gc3s_adjusted.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    M = args.metric
    # Separate random streams per analysis, so changing one part (e.g. the
    # number of median-regression replicates) leaves the others unchanged.
    rng_ols, rng_med, rng_gw, rng_match = (np.random.default_rng(s) for s in
                                           np.random.SeedSequence(args.seed).spawn(4))

    df = load_genes(args.genes, M)
    defs = definitions(df)
    d0 = defs[0]['key']
    logging.info('Definitions: %s (default %s)', ', '.join(d['key'] for d in defs), d0)

    model_rows, accounting = [], {}

    def add_rows(res, definition, label, model, frame, extra=None):
        for contrast, r in res.items():
            num, den = (('E_GR', 'R_GR') if contrast == 'elimination'
                        else ('R_GR', 'R_S'))
            n_num = int((frame['group'] == num).sum())
            n_den = int((frame['group'] == den).sum())
            row = {'definition': definition, 'label': label, 'model': model,
                   'contrast': contrast, 'numerator': num, 'denominator': den,
                   'n_num': n_num, 'n_den': n_den, 'n_model': len(frame),
                   **r, **(extra or {}), 'note': ''}
            if min(n_num, n_den) < args.min_group:
                row['note'] = f'fewer than {args.min_group} genes in a group'
            model_rows.append(row)

    # --- 1. models at every definition ----------------------------------------
    base_specs = [('group only', {}), ('GC3s', {'gc3s': True}),
                  ('GC3s + expression', {'gc3s': True, 'expr': True})]
    frames = {}
    for d in defs:
        frame, acct = model_frame(df, df[f"group_{d['key']}"], M)
        accounting[d['key']] = acct
        frames[d['key']] = frame
        logging.info('[%s] %s genes in the model (%s)', d['label'], f'{len(frame):,}',
                     ', '.join(f'{g} {n:,}' for g, n in acct['per_group'].items()))
        if len(frame) < 50:
            logging.warning('  too few genes; skipped')
            continue
        knots = {'gc3s': rcs_knots(frame['gc3s'], args.knots),
                 'expr': rcs_knots(frame['x'], args.expr_knots)}
        for name, spec in base_specs:
            res, info = fit_model(frame, spec, knots, args.n_boot, rng_ols)
            add_rows(res, d['key'], d['label'], name, frame)
            if name == 'group only':
                # Self-check: with only group dummies the contrasts are
                # differences of group means.
                gm = frame.groupby('group')['y'].mean()
                for c, (a, b) in (('elimination', ('E_GR', 'R_GR')),
                                  ('restriction', ('R_GR', 'R_S'))):
                    if c in res and a in gm and b in gm:
                        if not np.isclose(res[c]['estimate'], gm[a] - gm[b],
                                          atol=1e-9):
                            sys.exit(f'ERROR: group-only {c} contrast '
                                     f'{res[c]["estimate"]} is not the difference '
                                     f'of group means {gm[a] - gm[b]}')
        el = [r for r in model_rows if r['definition'] == d['key']
              and r['model'] == 'GC3s + expression' and r['contrast'] == 'elimination']
        if el:
            logging.info('  elimination, GC3s + expression: %+.3f %s, P = %s',
                         el[0]['estimate'], ci_txt(el[0]['ci_lo'], el[0]['ci_hi']),
                         ptxt(el[0]['p']))

    # --- default definition: robustness models --------------------------------
    f0 = frames[d0]
    if len(f0) < 50 or min((f0['group'] == g).sum() for g in ('R_S', 'R_GR', 'E_GR')) \
            < args.min_group:
        sys.exit(f'ERROR: the default definition ({d0}) leaves too few genes in '
                 f'R_S, R_GR or E_GR for the models.')
    k0 = {'gc3s': rcs_knots(f0['gc3s'], args.knots),
          'expr': rcs_knots(f0['x'], args.expr_knots)}
    logging.info('[%s] robustness models', defs[0]['label'])
    for name, spec in (('GC3s, median regression', {'gc3s': True}),
                       ('GC3s + expression, median regression',
                        {'gc3s': True, 'expr': True})):
        res, info = fit_model(f0, spec, k0, args.n_boot_median, rng_med,
                              median=True, lp_check=True)
        add_rows(res, d0, defs[0]['label'], name, f0,
                 {'irls_vs_lp_max_gap': info.get('irls_vs_lp_max_gap', np.nan)})
    fnc, acct_nc = model_frame(df, df[f'group_{d0}'], M, drop_ceiling=True)
    accounting[f'{d0}_no_ceiling'] = acct_nc
    knc = {'gc3s': rcs_knots(fnc['gc3s'], args.knots),
           'expr': rcs_knots(fnc['x'], args.expr_knots)}
    res, _ = fit_model(fnc, {'gc3s': True, 'expr': True}, knc, args.n_boot, rng_ols)
    add_rows(res, d0, defs[0]['label'], 'GC3s + expression, no ENC-ceiling genes', fnc)
    k7 = {'gc3s': rcs_knots(f0['gc3s'], 7), 'expr': k0['expr']}
    res, _ = fit_model(f0, {'gc3s': True, 'expr': True}, k7, args.n_boot, rng_ols)
    add_rows(res, d0, defs[0]['label'], 'GC3s (7 knots) + expression', f0)
    # As S12C matched on peak stage: germline-peak genes only (R_GR, E_GR,
    # R_Gbroad), with stage dummies; gives the elimination contrast only.
    fg = f0[f0['group'].isin(['R_GR', 'E_GR', 'R_Gbroad'])].reset_index(drop=True)
    kg = {'gc3s': rcs_knots(fg['gc3s'], args.knots),
          'expr': rcs_knots(fg['x'], args.expr_knots)}
    res, _ = fit_model(fg, {'gc3s': True, 'expr': True, 'stage': True}, kg,
                       args.n_boot, rng_ols)
    add_rows(res, d0, defs[0]['label'],
             'GC3s + expression + stage, germline-peak genes', fg)
    k0['centre'] = float(np.median(f0.loc[f0['group'].isin(['R_GR', 'E_GR']), 'gc3s']))
    inter = interaction_test(f0, k0)
    if inter:
        logging.info('  own GC3s slopes for R_GR and E_GR: chi2 = %.2f, df = %d, '
                     'P = %s', inter['chi2'], inter['df'], ptxt(inter['p']))

    models = pd.DataFrame(model_rows)

    # --- 2. genome-wide reference ------------------------------------------------
    logging.info('Genome-wide reference (eliminated vs retained)')
    gw_rows = []
    for name, spec, med in (('group only', {}, False),
                            ('GC3s', {'gc3s': True}, False),
                            ('GC3s + expression', {'gc3s': True, 'expr': True}, False),
                            ('GC3s + expression, median regression',
                             {'gc3s': True, 'expr': True}, True)):
        r, n = genome_wide(df, M, spec, args.knots, args.expr_knots,
                           args.n_boot_median if med else args.n_boot, rng_gw,
                           median=med)
        gw_rows.append({'definition': 'none', 'label': 'All genes', 'model': name,
                        'contrast': 'genome_wide', 'numerator': 'eliminated',
                        'denominator': 'retained', 'n_model': n, **r, 'note': ''})
        logging.info('  %-38s %+.3f %s', name, r['estimate'],
                     ci_txt(r['ci_lo'], r['ci_hi']))
    models = pd.concat([models, pd.DataFrame(gw_rows)], ignore_index=True)

    # --- 3. matched contrasts --------------------------------------------------
    logging.info('Matched contrasts (GC3s within %g; expression within %g where '
                 'named)', args.gc3s_caliper, args.caliper)
    matched, reps_all = matched_rows(f0.rename(columns={'y': 'v_metric'}),
                                     'v_metric', args.gc3s_caliper, args.caliper,
                                     args.n_match, args.n_match_boot, rng_match)

    # --- 4. strata ----------------------------------------------------------------
    strata = strata_table(f0)

    # --- tables -------------------------------------------------------------------
    out = args.outdir
    models.to_csv(os.path.join(out, 'gc3s_models.tsv'), sep='\t', index=False,
                  float_format='%.6g')
    matched.to_csv(os.path.join(out, 'gc3s_matched.tsv'), sep='\t', index=False,
                   float_format='%.6g')
    if reps_all:
        pd.concat(reps_all).to_csv(os.path.join(out, 'gc3s_match_replicates.tsv'),
                                   sep='\t', index=False, float_format='%.6g')
    strata.to_csv(os.path.join(out, 'gc3s_strata.tsv'), sep='\t', index=False,
                  float_format='%.6g')

    # --- summary --------------------------------------------------------------------
    def get(definition, model, contrast):
        r = models[(models['definition'] == definition) & (models['model'] == model)
                   & (models['contrast'] == contrast)]
        return r.iloc[0] if len(r) else None

    def line(r):
        return (f"{r['estimate']:+.3f} {ci_txt(r['ci_lo'], r['ci_hi'])}, "
                f"P = {ptxt(r['p'])}")

    L = [f'Codon bias ({M}) of germline-restricted genes, GC3s controlled',
         '=' * 78, f'genes   {args.genes}',
         f"metric  {M} ({'higher' if M == 'dENC' else 'lower'} = more codon bias); "
         f'GC3s spline {args.knots} knots, '
         f'expression spline {args.expr_knots} knots',
         'OLS: estimate [95% bootstrap CI], HC3 Wald P. Median regression: '
         'bootstrap SE and CI.', '',
         f'Genes in the model at the default definition ({d0}):']
    a0 = accounting[d0]
    L.append(f"  in the five groups {a0['in_model_groups']:,}; missing a value "
             f"{a0['missing_value']:,}; used {a0['used']:,} ("
             + ', '.join(f'{g} {n:,}' for g, n in a0['per_group'].items()) + ')')
    L.append(f"  without ENC-ceiling genes: {acct_nc['used']:,} "
             f"({acct_nc['at_ceiling_dropped']:,} dropped)")
    gmed = f0.groupby('group')['gc3s'].median()
    L.append('  median GC3s: ' + ', '.join(f'{g} {gmed[g]:.3f}' for g in MODEL_GROUPS
                                          if g in gmed))
    hi_e = f0.loc[f0['group'] == 'E_GR', 'gc3s']
    for g in ('R_GR', 'R_S'):
        top = f0.loc[f0['group'] == g, 'gc3s'].quantile(0.99)
        L.append(f'  E_GR genes above the 99th percentile of {g} GC3s ({top:.3f}): '
                 f'{int((hi_e > top).sum())} of {len(hi_e)}')

    for contrast in ('elimination', 'restriction'):
        title = ('Elimination, E_GR - R_GR (eliminated vs retained, both '
                 'germline-restricted)' if contrast == 'elimination' else
                 'Restriction, R_GR - R_S (retained: germline-restricted vs '
                 'somatic peak)')
        L += ['', title, f'  default definition ({defs[0]["label"]}):']
        for model in ('group only', 'GC3s', 'GC3s + expression',
                      'GC3s, median regression',
                      'GC3s + expression, median regression',
                      'GC3s + expression, no ENC-ceiling genes',
                      'GC3s (7 knots) + expression',
                      'GC3s + expression + stage, germline-peak genes'):
            r = get(d0, model, contrast)
            if r is not None:
                L.append(f'    {model:<42} {line(r)}' + (f"  [{r['note']}]"
                                                         if r['note'] else ''))
        L.append('  every definition, GC3s + expression (OLS):')
        for d in defs:
            r = get(d['key'], 'GC3s + expression', contrast)
            if r is not None:
                L.append(f"    {d['label']:<18} {line(r)}  (n {int(r['n_num']):,} / "
                         f"{int(r['n_den']):,})" + (f"  [{r['note']}]" if r['note'] else ''))
    if inter:
        L += ['', 'Own GC3s slopes for R_GR and E_GR (GC3s + expression model, at '
              f"GC3s {k0['centre']:.3f}): chi2 = {inter['chi2']:.2f}, df = "
              f"{inter['df']}, P = {ptxt(inter['p'])}"]
        for n, s in inter['slopes'].items():
            L.append(f"  {n.replace('gc3s_x_', '')}: GC3s slope {s:+.2f} per unit "
                     f"GC3s away from the common curve (SE {inter['slope_se'][n]:.2f})")
    L += ['', 'Genome-wide, eliminated - retained (every gene with an expression '
          'value):']
    for r in gw_rows:
        L.append(f"  {r['model']:<42} {line(r)}  (n {r['n_num']:,} / {r['n_den']:,})")
    L += ['', f'Matched (differences of medians; GC3s within {args.gc3s_caliper:g}, '
          f'expression within {args.caliper:g} where named):']
    for r in matched.itertuples():
        L.append(f'  {r.contrast:<11} on {r.matched_on:<17} {r.diff:+.3f} '
                 f'{ci_txt(r.ci_lo, r.ci_hi)} (matching alone '
                 f'{ci_txt(r.match_lo, r.match_hi)}); {r.median_n_matched:.0f} of '
                 f'{r.n_numerator} matched; GC3s balance P '
                 f'{ptxt(r.median_gc3s_balance_p)}, expression balance P '
                 f'{ptxt(r.median_expr_balance_p)}')
    L += ['', f'Median {M} in GC3s bins (default definition; descriptive):',
          '  GC3s         R_S (n)          R_GR (n)         E_GR (n)']
    for r in strata.itertuples():
        if r.n_R_S + r.n_R_GR + r.n_E_GR == 0:
            continue
        cells = []
        for g in ('R_S', 'R_GR', 'E_GR'):
            n, m = getattr(r, f'n_{g}'), getattr(r, f'median_{g}')
            cells.append(f'{m:6.2f} ({n:>5,})' if n else f'{"-":>6} ({0:>5})')
        L.append(f'  {r.gc3s_lo:.2f}-{r.gc3s_hi:.2f}  ' + '   '.join(cells))

    # What to say
    el = get(d0, 'GC3s + expression', 'elimination')
    el_med = get(d0, 'GC3s + expression, median regression', 'elimination')
    rs = get(d0, 'GC3s', 'restriction')
    gw = [r for r in gw_rows if r['model'] == 'GC3s + expression'][0]
    gw_med = [r for r in gw_rows if r['model'] ==
              'GC3s + expression, median regression'][0]
    L += ['', 'What to say']
    if el is not None:
        L.append(f"  Elimination at equal GC3s and expression (OLS, {int(el['n_model']):,} "
                 f"genes): {line(el)}; median regression {line(el_med) if el_med is not None else 'NA'}. "
                 f"Genome-wide adjusted difference: OLS {gw['estimate']:+.2f}, median "
                 f"{gw_med['estimate']:+.2f}.")
        L.append(f"    reading (OLS): {elimination_verdict(el, gw)}")
        if el_med is not None:
            L.append(f"    reading (median): {elimination_verdict(el_med, gw_med)}")
        allel = models[(models['model'] == 'GC3s + expression') &
                       (models['contrast'] == 'elimination') & (models['note'] == '')]
        if len(allel):
            L.append(f"    every definition: {allel['estimate'].min():+.2f} to "
                     f"{allel['estimate'].max():+.2f}; smallest P "
                     f"{ptxt(allel['p'].min())}; intervals below 0 at "
                     f"{int((allel['ci_hi'] < 0).sum())} of {len(allel)}")
    if rs is not None:
        L.append(f"  Restriction at equal GC3s (OLS): {line(rs)}; reading: "
                 f"{restriction_verdict(rs)}.")
        rs2 = get(d0, 'GC3s + expression', 'restriction')
        if rs2 is not None:
            L.append(f"    with expression level also held equal: {line(rs2)}; "
                     f"{restriction_verdict(rs2)}.")
        allrs = models[(models['model'] == 'GC3s') & (models['contrast'] == 'restriction')
                       & (models['note'] == '')]
        if len(allrs):
            L.append(f"    every definition (GC3s): {allrs['estimate'].min():+.2f} to "
                     f"{allrs['estimate'].max():+.2f}; largest P {ptxt(allrs['p'].max())}")
    for r in matched.itertuples():
        if r.matched_on == 'GC3s':
            L.append(f'  Matched on GC3s, {r.contrast}: {r.diff:+.2f} '
                     f'{ci_txt(r.ci_lo, r.ci_hi)} ({r.median_n_matched:.0f} of '
                     f'{r.n_numerator} matched).')
    with open(os.path.join(out, 'gc3s_adjusted_summary.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    for ln in L:
        logging.info(ln)

    prov = {'script': os.path.abspath(__file__), 'git': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).isoformat(), 'argv': sys.argv,
            'args': vars(args), 'inputs': {'genes': fingerprint(args.genes)},
            'definitions': [d['key'] for d in defs], 'accounting': accounting,
            'interaction_test': inter,
            'knots_default': {k: [float(v) for v in np.atleast_1d(x)]
                              for k, x in k0.items()},
            'versions': {'python': platform.python_version(),
                         'numpy': np.__version__, 'pandas': pd.__version__,
                         'scipy': scipy_version}}
    with open(os.path.join(out, 'provenance.json'), 'w') as fh:
        json.dump(prov, fh, indent=2, default=str)
    logging.info('Done.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
