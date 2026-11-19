#!/usr/bin/env python3
"""
germline_restricted_codon.py - is the weaker codon bias of eliminated genes a
germline-restriction effect?

The question (Jianbin Wang, meeting 2026-10-28)
-----------------------------------------------
Would codon usage help distinguish eliminated genes from genes that are only
repressed in the soma? A retained gene silenced in somatic cells escapes
somatic selection just as an eliminated gene does, and supplementary fig. S12C
compares their dN/dS. Codon bias asks the same question of a different kind of
selection. It is set by translational selection where a gene is translated,
which is the germline for both groups. If not being used in the soma were the
whole story, eliminated and retained genes with the same germline-restricted
expression should carry the same codon bias. Weaker bias in eliminated genes
would mean that selection is less effective in eliminated DNA itself (fig. 2C
shows weaker bias for eliminated genes overall).

What the script does
--------------------
It reads two tables that already exist and joins them by gene:

  genes_classified.tsv    from germline_restricted_dnds.py (S12B-C): each
                          gene's expression group at every definition of
                          germline restriction (group_<definition> columns),
                          its stage of peak expression and expression level
  codon_usage_per_gene.tsv  from codon_usage_bias.py (fig. 2B-C): ENC, GC3s,
                          dENC (ENC expected from GC3s minus observed; higher
                          = more codon bias) and the ENC ceiling flag

so the groups are exactly those of S12B-C and the codon metric exactly that of
fig. 2C. The groups (R_S, R_Gbroad, R_GR, E_GR, E_other, R_expr, E_expr), the
contrasts, the matching (candidate_lists, match_once) and the bootstrap of
medians (boot_medians) are imported from germline_restricted_dnds.py, and the
figure helpers from it and from functional_enrichment.py.

For every definition of germline restriction, every S12 contrast is computed
on --metric (default dENC) as a DIFFERENCE of medians (numerator minus
denominator; dENC can be zero or negative, so a ratio means nothing), with a
95% percentile bootstrap interval (each group resampled on its own),
Mann-Whitney P and BH q over every test. The two that answer the question:

  elimination  E_GR - R_GR   same expression domain, eliminated or not
  restriction  R_GR - R_S    does germline restriction alone weaken codon bias?

The elimination contrast is repeated with each eliminated germline-restricted
gene matched to a retained germline-restricted gene of the same peak stage and
expression level (log2 peak TPM + 1 within --caliper, as S12C), and again with
GC3s also matched (within --gc3s-caliper): dENC corrects for base composition
through Wright's curve, which is not a complete correction in this genome
(codon_usage_bias.py: the matched dENC contrast was stronger than the raw one).
At the default definition it is also repeated without genes at the ENC
ceiling of 61. For reference, the contrast over every joined gene (fig. 2C's
comparison) is reported first.

Writes, in --outdir
-------------------
  germline_codon_contrasts.{svg,pdf,png}   89 mm, the S12C layout: elimination
                          (filled circles) and restriction (open squares)
                          differences at every definition, matched rows below
  germline_codon_groups.tsv      n, median dENC [CI], GC3s, ENC, % at the
                                 ceiling, expression level, per definition x group
  germline_codon_contrasts.tsv   every contrast at every definition
  germline_codon_matched.tsv     the matched elimination contrasts
  germline_codon_match_replicates.tsv   one row per matching replicate
  germline_codon_genes.tsv       the joined per-gene table
  germline_codon_summary.txt     for reading, ending with what to say
  germline_restricted_codon.log
  provenance.json

Usage
-----
    python germline_restricted_codon.py \\
      --classified results/germline_restricted_dnds_20261019_153637/genes_classified.tsv \\
      --codon results/codon_usage_ascaris_20261027/codon_usage_per_gene.tsv \\
      --outdir results/germline_restricted_codon_20261029_120000

Dependencies: numpy, pandas, scipy, matplotlib, paper_colors, region_lib,
germline_restricted_dnds (and what it imports), functional_enrichment
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
from matplotlib.lines import Line2D  # noqa: E402  (after the backend is set)
from scipy.stats import mannwhitneyu  # noqa: E402
from scipy import __version__ as scipy_version  # noqa: E402

# The shared style must import, or the script stops: a figure drawn without it
# would look plausible and be wrong.
from paper_colors import (NEUTRAL, FONT, LINEWIDTH, apply_style,  # noqa: E402
                          make_outdir)
from region_lib import (benjamini_hochberg, fingerprint, git_commit,  # noqa: E402
                        setup_logging)
# One set of groups, contrasts, matching and figure layout with S12B-C.
from germline_restricted_dnds import (GROUP_MEMBERS, GROUP_DESC,  # noqa: E402
                                      SUMMARY_GROUPS, CONTRASTS, DRAWN,
                                      boot_medians, candidate_lists,
                                      match_once, contrast_height_mm,
                                      finish_figure)
from functional_enrichment import (blank_figure, add_axes_mm, side_text,  # noqa: E402
                                   place_with_table, PEAK_COL)

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Axis label per metric, for the figure.
METRIC_LABEL = {'dENC': 'ΔENC difference (difference of medians)',
                'ENC': 'ENC difference (difference of medians)',
                'GC3s': 'GC3s difference (difference of medians)'}


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def definitions(classified):
    """The definitions of germline restriction the S12 run used, from its
    group_<key> columns, in S12C's order and with its labels: floors, then tau
    thresholds, then 'calls agree'."""
    defs = []
    for col in classified.columns:
        if not col.startswith('group_') or col == 'group_default':
            continue
        key = col[len('group_'):]
        m_floor = re.fullmatch(r'floor_q([0-9.]+)', key)
        m_tau = re.fullmatch(r'tau_([0-9.]+)', key)
        m_conc = re.fullmatch(r'floor_q([0-9.]+)_concordant', key)
        if m_floor:
            q = float(m_floor.group(1))
            defs.append({'key': key, 'block': 'floor', 'order': (0, -q),
                         'label': f'Floor {round(100 * q)}th pct'})
        elif m_tau:
            t = float(m_tau.group(1))
            defs.append({'key': key, 'block': 'tau', 'order': (1, t),
                         'label': f'τ ≥ {t:g}'})
        elif m_conc:
            defs.append({'key': key, 'block': 'robust', 'order': (2, 0),
                         'label': 'Calls agree'})
        else:
            logging.warning('  column %s is not a definition this script knows; '
                            'skipped', col)
    defs.sort(key=lambda d: d['order'])
    if not defs:
        sys.exit('ERROR: the classified table has no group_<definition> columns; '
                 'is it genes_classified.tsv from germline_restricted_dnds.py?')
    return defs


def load_tables(args):
    """Join the classified genes to the codon table by gene, checking that the
    two agree on every gene's region."""
    cl = pd.read_csv(args.classified, sep='\t', dtype={'gene1': str})
    cu = pd.read_csv(args.codon, sep='\t', dtype={'gene_id': str})
    need_cl = ['gene1', 'region1', PEAK_COL, 'log_expr', 'group_default']
    need_cu = ['gene_id', 'region', args.metric, 'GC3s', 'ENC', 'ENC_at_ceiling']
    for frame, need, path in ((cl, need_cl, args.classified),
                              (cu, need_cu, args.codon)):
        miss = [c for c in need if c not in frame.columns]
        if miss:
            sys.exit(f'ERROR: {path} lacks {miss}.\n  columns: {list(frame.columns)}')
    for frame, col, path in ((cl, 'gene1', args.classified),
                             (cu, 'gene_id', args.codon)):
        dup = frame[col].duplicated()
        if dup.any():
            sys.exit(f'ERROR: {int(dup.sum())} duplicated gene IDs in {path}.')
    logging.info('Classified genes: %s (%s)', f'{len(cl):,}',
                 os.path.basename(args.classified))
    logging.info('Codon table: %s genes (%s)', f'{len(cu):,}',
                 os.path.basename(args.codon))
    cu = cu.rename(columns={'gene_id': 'gene1', 'region': 'codon_region'})
    keep_cu = ['gene1', 'codon_region', args.metric, 'GC3s', 'ENC',
               'ENC_at_ceiling'] + (['total_codons'] if 'total_codons' in cu
                                    else [])
    keep_cu = list(dict.fromkeys(keep_cu))
    df = cl.merge(cu[keep_cu], on='gene1', how='inner')
    both = df['codon_region'].notna()
    disagree = both & (df['codon_region'] != df['region1'])
    if disagree.any() and not args.allow_region_mismatch:
        ex = df.loc[disagree, ['gene1', 'region1', 'codon_region']].head(5)
        sys.exit(f'ERROR: {int(disagree.sum())} genes carry a different region in '
                 f'the two tables, so they were built from different annotated '
                 f'tables:\n{ex.to_string(index=False)}\n  Rerun one of them, or '
                 f'pass --allow-region-mismatch to drop those genes.')
    df = df[~disagree].copy()
    df['ENC_at_ceiling'] = df['ENC_at_ceiling'].astype(str).str.lower().isin(
        ['true', '1', '1.0'])
    acct = {
        'classified': {r: int((cl['region1'] == r).sum())
                       for r in ('retained', 'eliminated')},
        'codon': {r: int((cu['codon_region'] == r).sum())
                  for r in ('retained', 'eliminated')},
        'joined': {r: int((df['region1'] == r).sum())
                   for r in ('retained', 'eliminated')},
        'region_disagreements_dropped': int(disagree.sum()),
    }
    logging.info('  joined: %s retained, %s eliminated (%s dropped for a region '
                 'disagreement)', f"{acct['joined']['retained']:,}",
                 f"{acct['joined']['eliminated']:,}",
                 acct['region_disagreements_dropped'])
    return df, acct


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def group_mask(df, labels, group, metric):
    """Genes of `group` (S12's members) with a value of the metric."""
    return labels.isin(GROUP_MEMBERS[group]) & df[metric].notna()


def diff_contrast(a, b, rng, n_boot, min_n, level=0.95):
    """Difference of medians a - b, percentile bootstrap interval (groups
    resampled independently), Mann-Whitney P."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    out = {'n_num': len(a), 'n_den': len(b),
           'median_num': float(np.median(a)) if len(a) else np.nan,
           'median_den': float(np.median(b)) if len(b) else np.nan,
           'diff': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan,
           'U': np.nan, 'p': np.nan, 'note': ''}
    if len(a) < min_n or len(b) < min_n:
        out['note'] = f'fewer than {min_n} genes in a group; not tested'
        return out
    out['diff'] = out['median_num'] - out['median_den']
    d = boot_medians(a, rng, n_boot) - boot_medians(b, rng, n_boot)
    alpha = (1 - level) / 2
    out['ci_lo'], out['ci_hi'] = (float(x) for x in
                                  np.quantile(d, [alpha, 1 - alpha]))
    u, p = mannwhitneyu(a, b, alternative='two-sided')
    out['U'], out['p'] = float(u), float(p)
    return out


def matched_difference(e, r, caliper, n_rep, n_boot, rng, gc3s_caliper=None,
                       level=0.95):
    """Elimination contrast on matched genes, as a difference of medians.

    The matching of germline_restricted_dnds.matched_contrast(): each
    eliminated gene, in random order, takes an unused retained gene of the
    same peak stage with expression level within `caliper`
    (candidate_lists, match_once, imported); with `gc3s_caliper`, the
    candidate must also lie within that GC3s. `e`, `r`: frames with columns
    peak, x (log2 peak TPM + 1), gc3s and v (the metric). Observed replicates
    give the estimate and the spread due to matching alone; resample-then-match
    replicates give the interval.
    """
    codes, _ = pd.factorize(pd.concat([e['peak'], r['peak']], ignore_index=True))
    e_stage, r_stage = codes[:len(e)], codes[len(e):]
    e_x, r_x = e['x'].to_numpy(float), r['x'].to_numpy(float)
    e_g, r_g = e['gc3s'].to_numpy(float), r['gc3s'].to_numpy(float)
    e_v, r_v = e['v'].to_numpy(float), r['v'].to_numpy(float)

    def cands(ie, ir):
        c = candidate_lists(e_stage[ie], e_x[ie], r_stage[ir], r_x[ir], caliper)
        if gc3s_caliper is not None:
            rg = r_g[ir]
            c = [ci[np.abs(rg[ci] - g) <= gc3s_caliper]
                 for ci, g in zip(c, e_g[ie])]
        return c

    all_e, all_r = np.arange(len(e_v)), np.arange(len(r_v))
    c0 = cands(all_e, all_r)
    n_none = int(sum(len(c) == 0 for c in c0))
    rows = []
    for rep in range(n_rep):
        pick = match_once(c0, len(r_v), rng)
        ok = pick >= 0
        if ok.sum() < 5:
            continue
        ev, rv = e_v[ok], r_v[pick[ok]]
        rows.append({'kind': 'observed', 'replicate': rep,
                     'n_matched': int(ok.sum()),
                     'median_eliminated': float(np.median(ev)),
                     'median_matched_retained': float(np.median(rv)),
                     'diff': float(np.median(ev) - np.median(rv)),
                     'p': float(mannwhitneyu(ev, rv, alternative='two-sided')[1]),
                     'level_balance_p': float(mannwhitneyu(
                         e_x[ok], r_x[pick[ok]], alternative='two-sided')[1]),
                     'gc3s_balance_p': float(mannwhitneyu(
                         e_g[ok], r_g[pick[ok]], alternative='two-sided')[1])})
    for rep in range(n_boot):
        ie = rng.integers(0, len(e_v), len(e_v))
        ir = rng.integers(0, len(r_v), len(r_v))
        pick = match_once(cands(ie, ir), len(ir), rng)
        ok = pick >= 0
        if ok.sum() < 5:
            continue
        ev, rv = e_v[ie][ok], r_v[ir][pick[ok]]
        rows.append({'kind': 'bootstrap', 'replicate': rep,
                     'n_matched': int(ok.sum()),
                     'median_eliminated': float(np.median(ev)),
                     'median_matched_retained': float(np.median(rv)),
                     'diff': float(np.median(ev) - np.median(rv)),
                     'p': np.nan, 'level_balance_p': np.nan,
                     'gc3s_balance_p': np.nan})
    reps = pd.DataFrame(rows)
    alpha = (1 - level) / 2
    out = {'n_eliminated': len(e_v), 'n_retained_pool': len(r_v),
           'n_no_candidate': n_none, 'caliper_log2': caliper,
           'gc3s_caliper': gc3s_caliper if gc3s_caliper is not None else np.nan,
           'n_rep_observed': 0, 'n_rep_bootstrap': 0, 'diff': np.nan,
           'match_lo': np.nan, 'match_hi': np.nan, 'ci_lo': np.nan,
           'ci_hi': np.nan, 'median_n_matched': np.nan,
           'pct_rep_p_lt_0.05': np.nan, 'median_p': np.nan,
           'median_level_balance_p': np.nan, 'median_gc3s_balance_p': np.nan}
    if len(reps):
        obs = reps[reps['kind'] == 'observed']
        bts = reps[reps['kind'] == 'bootstrap']
        out['n_rep_observed'], out['n_rep_bootstrap'] = len(obs), len(bts)
        if len(obs):
            out.update(diff=float(obs['diff'].median()),
                       match_lo=float(obs['diff'].quantile(alpha)),
                       match_hi=float(obs['diff'].quantile(1 - alpha)),
                       median_n_matched=float(obs['n_matched'].median()),
                       **{'pct_rep_p_lt_0.05': float(100 * (obs['p'] < 0.05).mean())},
                       median_p=float(obs['p'].median()),
                       median_level_balance_p=float(obs['level_balance_p'].median()),
                       median_gc3s_balance_p=float(obs['gc3s_balance_p'].median()))
        if len(bts):
            out.update(ci_lo=float(bts['diff'].quantile(alpha)),
                       ci_hi=float(bts['diff'].quantile(1 - alpha)))
    return out, reps


# ---------------------------------------------------------------------------
# Figure: the S12C layout with a linear axis
# ---------------------------------------------------------------------------

def draw_difference_forest(ax, rows, metric, key_ncol=1, ax_t=None):
    """Elimination and restriction differences at every definition.

    As germline_restricted_dnds.draw_forest (S12C), but on a linear axis with
    the rule at 0, since the contrasts are differences of medians. Filled
    circles: eliminated - retained, both germline-restricted; open squares:
    germline-restricted - somatic peak, retained. n in the table axes.
    """
    ys, y, prev = [], 0.0, None
    for r in rows:
        if prev is not None and r['block'] != prev:
            y -= 0.45
        ys.append(y)
        prev = r['block']
        y -= 1.0
    ys = np.array(ys)
    ends = [v for r in rows for c in DRAWN if r.get(c)
            for v in (r[c][0], r[c][1], r[c][2])
            if v is not None and np.isfinite(v)]
    lo, hi = min(ends + [-0.1]), max(ends + [0.1])
    pad = 0.06 * (hi - lo)
    x_lo, x_hi = lo - pad, hi + pad
    style = {'elimination': (+0.17, 'o', NEUTRAL['ink']),
             'restriction': (-0.17, 's', 'white')}
    for r, yy in zip(rows, ys):
        for c in DRAWN:
            v = r.get(c)
            if not v or v[0] is None or not np.isfinite(v[0]):
                continue
            est, a, b = v
            dy, marker, face = style[c]
            if a is not None and np.isfinite(a) and np.isfinite(b):
                ax.plot([a, b], [yy + dy] * 2, color=NEUTRAL['ink'],
                        lw=LINEWIDTH['thin'], solid_capstyle='butt', zorder=2)
            ax.plot([est], [yy + dy], ls='none', marker=marker, ms=3.4,
                    mfc=face, mec=NEUTRAL['ink'], mew=0.7, zorder=4)
        if ax_t is not None:
            side_text(ax_t, yy, r['n_text'], dx_pt=0.0, x=0.0)
        else:
            side_text(ax, yy, r['n_text'])
    if ax_t is not None:
        side_text(ax_t, ys[0] + 0.62, 'n', dx_pt=0.0, x=0.0, va='bottom')
    else:
        side_text(ax, ys[0] + 0.62, 'n', va='bottom')
    # The rule at 0 in data coordinates (axvline's y is in axes units, which
    # audit() reads as data units).
    y_lo, y_hi = ys[-1] - 0.7, ys[0] + 0.7
    ax.plot([0.0, 0.0], [y_lo, y_hi], color=NEUTRAL['rule'],
            lw=LINEWIDTH['rule'], ls=(0, (4, 2)), zorder=1)
    ax.set_xlim(x_lo, x_hi)
    ax.set_yticks(ys)
    ax.set_yticklabels([r['label'] for r in rows])
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(y_lo, y_hi)
    if ax_t is not None:
        ax_t.set_ylim(y_lo, y_hi)
    ax.set_xlabel(METRIC_LABEL.get(metric, f'{metric} difference'))
    h = [Line2D([], [], ls='none', marker='o', ms=3.4, mfc=NEUTRAL['ink'],
                mec=NEUTRAL['ink'], mew=0.7),
         Line2D([], [], ls='none', marker='s', ms=3.4, mfc='white',
                mec=NEUTRAL['ink'], mew=0.7)]
    ax.legend(h, ['Eliminated − retained, germline only',
                  'Germline only − somatic peak, retained'],
              loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=key_ncol,
              frameon=False, fontsize=FONT['legend'], borderaxespad=0.4,
              handletextpad=0.3, labelspacing=0.3, columnspacing=1.6)
    return {'x_range': (float(x_lo), float(x_hi))}


def plot_difference_panel(rows, outdir, stem, metric, height_mm=None):
    """One 89 mm column, the S12C geometry: the forest and its n column, laid
    out from the drawn text by place_with_table."""
    height = height_mm or contrast_height_mm(rows, 2)
    fig = blank_figure('single', height)
    top = 4.0 + 3.6 * 2
    ax = add_axes_mm(fig, left=26.0, bottom=11.0, right=17.0, top=top)
    ax_t = add_axes_mm(fig, left=89.0 - 16.0, bottom=11.0, right=1.0, top=top)
    ax_t.axis('off')
    info = draw_difference_forest(ax, rows, metric, key_ncol=1, ax_t=ax_t)
    place_with_table(fig, ax, ax_t)
    paths, problems = finish_figure(fig, outdir, stem)
    info['height_mm'] = height
    return paths, problems, info


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def ptxt(p):
    """P or q as plain text for the summary."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p == 0:
        return '< 1e-300'
    return f'{p:.3g}' if p >= 0.001 else f'{p:.1e}'


def ci_txt(lo, hi, nd=2):
    if lo is None or not (np.isfinite(lo) and np.isfinite(hi)):
        return '[NA]'
    return f'[{lo:+.{nd}f}, {hi:+.{nd}f}]'


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Codon bias (dENC) of eliminated vs retained '
                    'germline-restricted genes, on the groups of S12B-C.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--classified', required=True,
                   help='genes_classified.tsv from germline_restricted_dnds.py')
    p.add_argument('--codon', required=True,
                   help='codon_usage_per_gene.tsv from codon_usage_bias.py '
                        '(the Ascaris run of fig. 2B-C)')
    p.add_argument('--metric', default='dENC', choices=['dENC', 'ENC', 'GC3s'],
                   help='Codon metric to compare (default dENC, as fig. 2C)')
    p.add_argument('--caliper', type=float, default=0.25,
                   help='Matching caliper on log2(peak TPM + 1), as S12C '
                        '(default 0.25)')
    p.add_argument('--gc3s-caliper', type=float, default=0.02,
                   help='GC3s caliper for the GC3s-matched row (default 0.02)')
    p.add_argument('--n-match', type=int, default=500,
                   help='Matching replicates on the observed genes (default 500)')
    p.add_argument('--n-match-boot', type=int, default=1000,
                   help='Resample-then-match replicates (default 1000)')
    p.add_argument('--n-boot', type=int, default=2000)
    p.add_argument('--min-group', type=int, default=10,
                   help='A contrast is tested only if both groups hold at least '
                        'this many genes (default 10, as S12)')
    p.add_argument('--allow-region-mismatch', action='store_true',
                   help='Drop genes whose region differs between the two '
                        'tables instead of stopping')
    p.add_argument('--fig-height-mm', type=float, default=None,
                   help='Figure height (default: from the number of rows, as '
                        'S12C)')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None,
                   help='Default: results/germline_restricted_codon_<YYYYmmdd_HHMMSS>')
    a = p.parse_args(argv)
    if a.caliper <= 0 or a.gc3s_caliper <= 0:
        p.error('--caliper and --gc3s-caliper must be positive')
    if a.outdir is None:
        a.outdir = make_outdir('germline_restricted_codon', base='results',
                               date=datetime.now().strftime('%Y%m%d_%H%M%S'))
    return a


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'germline_restricted_codon')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('germline_restricted_codon.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()
    # Separate streams for the bootstrap and the matching, as S12.
    rng_boot, rng_match = (np.random.default_rng(s) for s in
                           np.random.SeedSequence(args.seed).spawn(2))
    M = args.metric

    df, acct = load_tables(args)
    defs = definitions(df)
    d0 = defs[0]['key']
    if d0.startswith('floor_q') and not df['group_default'].fillna('').equals(
            df[f'group_{d0}'].fillna('')):
        logging.warning('  group_default differs from group_%s; the default '
                        'definition here is %s', d0, d0)

    # --- reference: every joined gene, eliminated vs retained (fig. 2C) ------
    contrast_rows = []
    ref = diff_contrast(df.loc[df['region1'] == 'eliminated', M].dropna(),
                        df.loc[df['region1'] == 'retained', M].dropna(),
                        rng_boot, args.n_boot, args.min_group)
    contrast_rows.append({'definition': 'none', 'label': 'All joined genes',
                          'contrast': 'all_genes', 'numerator': 'eliminated',
                          'denominator': 'retained', 'gene_set': 'all',
                          'description': 'every joined gene, as fig. 2C',
                          **ref})

    # --- groups and contrasts at every definition ----------------------------
    group_rows = []
    for d in defs:
        labels = df[f"group_{d['key']}"]
        frame = df[labels.notna()]
        labels = labels[labels.notna()]
        for g in SUMMARY_GROUPS:
            m = group_mask(frame, labels, g, M)
            v = frame.loc[m, M].to_numpy(float)
            lo = hi = np.nan
            if len(v) >= args.min_group:
                bm = boot_medians(v, rng_boot, args.n_boot)
                lo, hi = (float(x) for x in np.quantile(bm, [0.025, 0.975]))
            group_rows.append({
                'definition': d['key'], 'group': g, 'description': GROUP_DESC[g],
                'n': len(v), f'median_{M}': float(np.median(v)) if len(v) else np.nan,
                'ci_lo': lo, 'ci_hi': hi,
                'median_GC3s': float(frame.loc[m, 'GC3s'].median()) if len(v) else np.nan,
                'median_ENC': float(frame.loc[m, 'ENC'].median()) if len(v) else np.nan,
                'pct_at_ceiling': (float(100 * frame.loc[m, 'ENC_at_ceiling'].mean())
                                   if len(v) else np.nan),
                'median_log2_expr': (float(frame.loc[m, 'log_expr'].median())
                                     if len(v) else np.nan)})
        for ckey, num, den, desc in CONTRASTS:
            a = frame.loc[group_mask(frame, labels, num, M), M]
            b = frame.loc[group_mask(frame, labels, den, M), M]
            res = diff_contrast(a, b, rng_boot, args.n_boot, args.min_group)
            contrast_rows.append({'definition': d['key'], 'label': d['label'],
                                  'contrast': ckey, 'numerator': num,
                                  'denominator': den, 'gene_set': 'all',
                                  'description': desc, **res})
    # The ENC ceiling: the default elimination contrast without ceiling genes.
    lab0 = df[f'group_{d0}']
    fr0 = df[lab0.notna() & ~df['ENC_at_ceiling']]
    l0 = lab0[fr0.index]
    res = diff_contrast(fr0.loc[group_mask(fr0, l0, 'E_GR', M), M],
                        fr0.loc[group_mask(fr0, l0, 'R_GR', M), M],
                        rng_boot, args.n_boot, args.min_group)
    contrast_rows.append({'definition': d0, 'label': defs[0]['label'],
                          'contrast': 'elimination', 'numerator': 'E_GR',
                          'denominator': 'R_GR', 'gene_set': 'no_ceiling',
                          'description': 'eliminated vs retained, both '
                                         'germline-restricted; ENC below the '
                                         'ceiling', **res})
    contrasts = pd.DataFrame(contrast_rows)
    contrasts['q_bh'] = benjamini_hochberg(contrasts['p'].to_numpy(float))
    groups_tab = pd.DataFrame(group_rows)

    # --- matched elimination contrasts ---------------------------------------
    match_specs = [(d0, None, 'Matched')]
    conc = [d['key'] for d in defs if d['block'] == 'robust']
    if conc:
        match_specs.append((conc[0], None, 'Agree, matched'))
    match_specs.append((d0, args.gc3s_caliper, 'Matched, GC3s'))
    matched_rows, rep_frames = [], []
    for key, gcal, label in match_specs:
        lab = df[f'group_{key}']
        frame = df[lab.notna()]
        lab = lab[lab.notna()]

        def pick(group):
            m = group_mask(frame, lab, group, M)
            return pd.DataFrame({'peak': frame.loc[m, PEAK_COL].to_numpy(),
                                 'x': frame.loc[m, 'log_expr'].to_numpy(float),
                                 'gc3s': frame.loc[m, 'GC3s'].to_numpy(float),
                                 'v': frame.loc[m, M].to_numpy(float)})
        e, r = pick('E_GR'), pick('R_GR')
        if len(e) < args.min_group or len(r) < args.min_group:
            logging.warning('  %s (%s): %d eliminated, %d retained germline-'
                            'restricted genes; not run', label, key, len(e), len(r))
            continue
        res, reps = matched_difference(e, r, args.caliper, args.n_match,
                                       args.n_match_boot, rng_match,
                                       gc3s_caliper=gcal)
        res.update(definition=key, label=label)
        matched_rows.append(res)
        reps.insert(0, 'label', label)
        rep_frames.append(reps)
        logging.info('  %-15s diff %+.3f, CI %s, %s of %s matched (median); '
                     '%s without a candidate', label, res['diff'],
                     ci_txt(res['ci_lo'], res['ci_hi']),
                     res['median_n_matched'], res['n_eliminated'],
                     res['n_no_candidate'])
    matched = pd.DataFrame(matched_rows)

    # --- tables ----------------------------------------------------------------
    out = args.outdir
    groups_tab.to_csv(os.path.join(out, 'germline_codon_groups.tsv'), sep='\t',
                      index=False, float_format='%.6g')
    contrasts.to_csv(os.path.join(out, 'germline_codon_contrasts.tsv'), sep='\t',
                     index=False, float_format='%.6g')
    matched.to_csv(os.path.join(out, 'germline_codon_matched.tsv'), sep='\t',
                   index=False, float_format='%.6g')
    if rep_frames:
        pd.concat(rep_frames).to_csv(
            os.path.join(out, 'germline_codon_match_replicates.tsv'), sep='\t',
            index=False, float_format='%.6g')
    gcols = (['gene1', 'region1', PEAK_COL, 'log_expr', 'group_default'] +
             [f"group_{d['key']}" for d in defs] +
             list(dict.fromkeys([M, 'GC3s', 'ENC', 'ENC_at_ceiling'])) +
             [c for c in ('total_codons',) if c in df.columns])
    df[gcols].to_csv(os.path.join(out, 'germline_codon_genes.tsv'), sep='\t',
                     index=False, float_format='%.6g')

    # --- figure ------------------------------------------------------------------
    def look(key, ckey):
        r = contrasts[(contrasts['definition'] == key) &
                      (contrasts['contrast'] == ckey) &
                      (contrasts['gene_set'] == 'all')]
        return r.iloc[0] if len(r) and np.isfinite(r.iloc[0]['diff']) else None
    forest_rows = []
    for d in defs:
        el, rs = look(d['key'], 'elimination'), look(d['key'], 'restriction')
        forest_rows.append({
            'label': d['label'], 'block': d['block'],
            'elimination': (el['diff'], el['ci_lo'], el['ci_hi']) if el is not None else None,
            'restriction': (rs['diff'], rs['ci_lo'], rs['ci_hi']) if rs is not None else None,
            'n_text': (f"{int(el['n_num']):,} / {int(el['n_den']):,}"
                       if el is not None else '')})
    for m in matched_rows:
        forest_rows.append({
            'label': m['label'], 'block': 'matched',
            'elimination': (m['diff'], m['ci_lo'], m['ci_hi']),
            'restriction': None,
            'n_text': (f"{int(round(m['median_n_matched'])):,} pairs"
                       if np.isfinite(m['median_n_matched']) else '')})
    _, fig_problems, fig_info = plot_difference_panel(
        forest_rows, out, 'germline_codon_contrasts', M, args.fig_height_mm)

    # --- summary -------------------------------------------------------------------
    L = [f'Codon bias ({M}) of germline-restricted genes, eliminated vs retained',
         '=' * 78,
         f'classified  {args.classified}',
         f'codon       {args.codon}',
         f'metric      {M} ({"higher = more codon bias" if M == "dENC" else ""})',
         '', 'Genes                retained  eliminated']
    for k in ('classified', 'codon', 'joined'):
        L.append(f'  {k:<18} {acct[k]["retained"]:>8,}  {acct[k]["eliminated"]:>10,}')
    L.append(f'  region disagreements dropped: {acct["region_disagreements_dropped"]}')
    L += ['', f'Groups at the default definition ({d0}): median {M} [95% CI], n']
    for r in groups_tab[groups_tab['definition'] == d0].itertuples():
        med = getattr(r, f'median_{M}')
        ci = (f'[{r.ci_lo:.2f}, {r.ci_hi:.2f}]' if np.isfinite(r.ci_lo) else '[NA]')
        L.append(f'  {r.group:<9} {med:7.3f} {ci:<16} n = {r.n:<6,} GC3s '
                 f'{r.median_GC3s:.3f}  ceiling {r.pct_at_ceiling:4.1f}%  '
                 f'log2 expr {r.median_log2_expr:.2f}')
    L += ['', f'Contrasts (difference of medians of {M}, numerator minus '
          'denominator; q: BH over every row)']
    for r in contrasts.itertuples():
        head = f'  {r.label:<17} {r.contrast:<17} {r.gene_set:<10}'
        if r.note:
            L.append(head + f' {r.note}')
            continue
        L.append(head + f' {r.diff:+.3f} {ci_txt(r.ci_lo, r.ci_hi)}  '
                 f'({r.median_num:.3f} vs {r.median_den:.3f}; n {r.n_num:,} / '
                 f'{r.n_den:,}; P = {ptxt(r.p)}, q = {ptxt(r.q_bh)})')
    L += ['', 'Matched elimination contrasts (stage exact, log2 expression within '
          f'{args.caliper:g}; GC3s within {args.gc3s_caliper:g} where named)']
    for m in matched_rows:
        L.append(f"  {m['label']:<15} {m['diff']:+.3f} {ci_txt(m['ci_lo'], m['ci_hi'])} "
                 f"(matching alone {ci_txt(m['match_lo'], m['match_hi'])}); "
                 f"{m['median_n_matched']:.0f} of {m['n_eliminated']} matched, "
                 f"{m['n_no_candidate']} without a candidate; "
                 f"{m['pct_rep_p_lt_0.05']:.0f}% of replicates P < 0.05; "
                 f"balance P: expression {ptxt(m['median_level_balance_p'])}, "
                 f"GC3s {ptxt(m['median_gc3s_balance_p'])}")
    if fig_problems:
        L += ['', 'FIGURE AUDIT: ' + '; '.join(fig_problems)]
    # What to say
    el = look(d0, 'elimination')
    rs = look(d0, 'restriction')
    L += ['', 'What to say']
    if ref['note'] == '':
        L.append(f"  All joined genes: eliminated {ref['median_num']:.2f} vs "
                 f"retained {ref['median_den']:.2f} ({ref['diff']:+.2f} "
                 f"{ci_txt(ref['ci_lo'], ref['ci_hi'])}; P = {ptxt(ref['p'])}).")
    def verb(row):
        if row['p'] >= 0.05:
            return 'do not differ in codon bias from'
        return ('carry weaker codon bias than' if row['diff'] < 0 else
                'carry stronger codon bias than')
    if el is not None:
        L.append(f"  Germline-restricted genes ({defs[0]['label']}): eliminated "
                 f"genes {verb(el)} equally restricted retained "
                 f"genes ({M} {el['median_num']:.2f} vs {el['median_den']:.2f}, "
                 f"{el['diff']:+.2f} {ci_txt(el['ci_lo'], el['ci_hi'])}; "
                 f"P = {ptxt(el['p'])}; n {int(el['n_num'])} / {int(el['n_den'])}).")
    for m in matched_rows:
        L.append(f"  {m['label']}: {m['diff']:+.2f} {ci_txt(m['ci_lo'], m['ci_hi'])}.")
    if rs is not None:
        L.append(f"  Retained genes: germline-restricted genes {verb(rs)} "
                 f"somatic-peak genes ({rs['diff']:+.2f} "
                 f"{ci_txt(rs['ci_lo'], rs['ci_hi'])}; P = {ptxt(rs['p'])}).")
    L.append('  Reading: an elimination difference that holds after matching, '
             'with no restriction difference, points to weaker selection in '
             'eliminated DNA itself; a restriction difference with no '
             'elimination difference means the lower codon bias of eliminated '
             'genes follows their expression domain.')
    with open(os.path.join(out, 'germline_codon_summary.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    for line in L:
        logging.info(line)

    prov = {'script': os.path.abspath(__file__), 'git': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).isoformat(), 'argv': sys.argv,
            'args': vars(args),
            'inputs': {'classified': fingerprint(args.classified),
                       'codon': fingerprint(args.codon)},
            'accounting': acct, 'definitions': [d['key'] for d in defs],
            'figure': {'problems': fig_problems, **fig_info},
            'versions': {'python': platform.python_version(),
                         'numpy': np.__version__, 'pandas': pd.__version__,
                         'scipy': scipy_version,
                         'matplotlib': matplotlib.__version__}}
    with open(os.path.join(out, 'provenance.json'), 'w') as fh:
        json.dump(prov, fh, indent=2, default=str)
    logging.info('Done.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
