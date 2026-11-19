#!/usr/bin/env python3
"""
germline_restricted_dnds.py - is the dN/dS excess of eliminated genes a
germline-restriction effect?

The question (M. Labrador, 2026-10-16)
--------------------------------------
Germline genes that are silenced in somatic tissue -- in animals without
programmed DNA elimination, or retained Ascaris genes expressed only in the
germline -- also escape somatic constraint. If narrow, germline-only expression
is what raises dN/dS, eliminated genes are simply germline-restricted genes and
elimination adds nothing. Figure 6B does not settle this: it groups genes by the
stage where expression peaks, and a testis-peak gene can still be expressed in
somatic tissue.

What the script does
--------------------
1. Joins the annotated rate table to the expression workbook through the same
   position bridge as Figure 6 (`--bridge workbook_model_bridge.tsv`, read by
   `workbook_model_bridge.read_bridge`); `load_workbook` and the stage
   definitions are imported from functional_enrichment.py, so the two analyses
   share one join and one set of stage definitions. `--name-join` reproduces
   the superseded Target-name join of the 2026-10-16 run (results draft note
   60), which matched each model to a different gene.

2. Defines germline-restricted expression. A gene is germline-restricted when
   its expression peaks in a germline stage (testis, ovary, zygote or early
   embryo) and every somatic stage (late embryo, larvae, muscle, intestine,
   carcass) is at or below a noise floor. The floor for each somatic stage is a
   quantile (default the 95th) of that stage's TPM among genes the workbook
   itself calls eliminated. Their DNA is absent from somatic cells, so the TPM
   they show there is mapping noise, plus whatever the few germline precursor
   cells in whole-embryo and larval samples contribute. A retained gene at or
   below the floor is indistinguishable from a gene that is not there. Tissue
   specificity (tau on log2(TPM + 1) across the nine stages, with a germline
   peak) gives a second, floor-free definition.

3. Compares dN/dS between groups at each definition:
     R_S       retained, expression peaks in a somatic stage
     R_Gbroad  retained, germline peak, detected in soma above the floor
     R_GR      retained, germline-restricted
     E_GR      eliminated, germline-restricted
     E_other   eliminated, any other expressed gene (somatic peak, or germline
               peak with somatic signal above the floor)
   The two contrasts the question turns on, and the only ones drawn:
     elimination   E_GR / R_GR   same expression domain, different elimination
     restriction   R_GR / R_S    does germline restriction alone raise dN/dS?
   Folds are ratios of medians with percentile bootstrap intervals (each group
   resampled independently), P from two-sided Mann-Whitney U, and BH q across
   every test in contrasts.tsv.

4. Repeats the elimination contrast with every eliminated germline-restricted
   gene matched to a retained germline-restricted gene of the same peak stage
   and similar expression level (caliper on log2 of the peak TPM + 1). Matching
   is without replacement in random order, as in the GC3s match of S2. The point
   estimate is the median fold over matching replicates on the observed genes;
   the interval comes from replicates that resample both groups before matching,
   so it carries sampling as well as matching variability and can be read beside
   the bootstrap intervals of the unmatched rows.

5. Fits the Figure 6E model family (`fit_model`, imported from
   multivariate_insulation.py: standardized OLS, 95% CI = beta +/- t x SE) with
   tissue specificity, germline peak and germline restriction added, on raw and
   on ranked dN/dS, to see whether they absorb the elimination coefficient.

6. Checks the elimination labels against the workbook. The rate table's region1
   comes from the v50 eliminated-region BED; the workbook carries its own
   elimination call (on v37) for the gene whose expression is joined. With the
   name join the two calls agreed at chance, which is how note 60 was found;
   with the position bridge they agree for 95% of eliminated and 99.9% of
   retained genes, the rest being boundary differences between the two
   assemblies' calls. The table keeps each model's miniprot Rank, and the
   concordant subset, where both calls agree, is analysed separately.

Supplementary figure S12 (2026-10-19)
-------------------------------------
This control is panels B and C of S12; panel A is the dN/dS-against-expression
scatter, `S6_dnds_vs_expression` from functional_enrichment.py, and panel D the
within-class control of dnds_within_class.py. The panels are saved separately
and placed at assembly, as Figure 6 is: A | B on the first row, at one height
(`--b-height-mm`, which run_germline_restricted_dnds.sh reads from the newest
S6_dnds_vs_expression.svg), and C | D on the second, C at `--c-width` (default
89 mm since 2026-10-19, so D fits beside it) and `--c-height-mm`, which D's
wrapper reads. C's n column is a table axes of its own, as 6C's is: hung off
the plot's edge at 89 mm it left the plot 52% of the width, under audit()'s
55%.
No panel letters are drawn; they are added at assembly at 10 pt, and each key
starts at its axes' left edge, clear of the letter's corner. The 185 mm
two-panel figure of the 2026-10-17 run is no longer written.

Writes (into --outdir)
----------------------
    S12B_germline_groups.{svg,pdf,png}     89 mm wide, --b-height-mm tall:
                                  dN/dS by group at the default definition
    S12C_germline_contrasts.{svg,pdf,png}  --c-width wide (default 89 mm):
                                  the elimination and restriction contrasts at
                                  every definition
    genes_classified.tsv          one row per analysed gene, every derived column
    somatic_floor.tsv             floor per somatic stage, by source and quantile
    group_summary.tsv             n, median and CI, quartiles, tau, expression
                                  level and peak-stage mix, per definition x group
    contrasts.tsv                 every contrast at every definition
    matched_contrast.tsv          one row per matching replicate
    models.tsv                    every coefficient of every model
    concordance.tsv               region1 x workbook elimination call
    germline_restricted_summary.txt
    provenance.json
    germline_restricted_dnds.log

Usage
-----
    python germline_restricted_dnds.py \\
      --kaks results/annotate_regions_20260817/kaks_annotated.tsv \\
      --annotation Ascaris_genes_15714_annotation_expression.xlsx \\
      --bridge results/workbook_model_bridge_20261017_103009/workbook_model_bridge.tsv \\
      --gff1 As_miniprot.gff \\
      --outdir results/germline_restricted_dnds_20261017_160000

Dependencies: numpy, pandas (+ openpyxl), scipy, matplotlib, paper_colors,
region_lib, rate_tables, functional_enrichment, multivariate_insulation,
workbook_model_bridge
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, NullLocator
from scipy.stats import mannwhitneyu, rankdata
from scipy import __version__ as scipy_version

# A script that cannot import the shared style must stop, not fall back to
# defaults: a figure drawn without it would look plausible and be wrong.
from paper_colors import (REGION, NEUTRAL, FONT, LINEWIDTH, apply_style, audit,
                          save, make_outdir)
from region_lib import (parse_gff_models, select_primary_alignment,
                        benjamini_hochberg, fingerprint, git_commit,
                        setup_logging)
from rate_tables import read_rate_table
# One join, one set of stage definitions and one figure toolkit with Figure 6.
from functional_enrichment import (target_bridge, load_workbook, STAGES,
                                   GERMLINE, SOMATIC, LOW, PEAK_COL,
                                   CATEGORY_OF, blank_figure, add_axes_mm,
                                   style_box, side_text, fit_to_canvas,
                                   offcanvas_text, place_with_table)
# The Figure 6E model, unchanged: standardized OLS, t-based 95% intervals.
from multivariate_insulation import fit_model
# The position join of Figure 6, read in one place.
from workbook_model_bridge import read_bridge, shared_gene_counts

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Somatic stages that must sit at or below the floor. Adult muscle, intestine
#: and carcass are dissected, wholly somatic tissue; late embryo and larvae are
#: whole animals and carry a few germline precursor cells, which is why the
#: floor is set per stage from genes that are absent from the soma rather than
#: at zero.
FLOOR_STAGES = list(SOMATIC)
ADULT_SOMATIC = ['muscle', 'intestine', 'carcass']

#: Groups drawn in S12B, left to right: retained first, eliminated second.
GROUP_ORDER = ['R_S', 'R_Gbroad', 'R_GR', 'E_GR', 'E_other']
GROUP_REGION = {'R_S': 'retained', 'R_Gbroad': 'retained', 'R_GR': 'retained',
                'E_GR': 'eliminated', 'E_other': 'eliminated',
                'R_expr': 'retained', 'E_expr': 'eliminated'}
GROUP_TICK = {'R_S': 'Somatic\npeak', 'R_Gbroad': 'Germline,\nbroad',
              'R_GR': 'Germline\nonly', 'E_GR': 'Germline\nonly',
              'E_other': 'Other'}
GROUP_DESC = {
    'R_S': 'retained, peak in a somatic stage',
    'R_Gbroad': 'retained, germline peak, soma above the floor',
    'R_GR': 'retained, germline-restricted',
    'E_GR': 'eliminated, germline-restricted',
    'E_other': 'eliminated, other expressed',
    'R_expr': 'retained, every expressed gene',
    'E_expr': 'eliminated, every expressed gene',
}
#: Fine labels written by classify() -> analysis groups.
GROUP_MEMBERS = {
    'R_S': ['R_S'], 'R_Gbroad': ['R_Gbroad'], 'R_GR': ['R_GR'],
    'E_GR': ['E_GR'], 'E_other': ['E_S', 'E_Gbroad'],
    'R_expr': ['R_S', 'R_Gbroad', 'R_GR'],
    'E_expr': ['E_S', 'E_Gbroad', 'E_GR'],
}
SUMMARY_GROUPS = GROUP_ORDER + ['R_expr', 'E_expr']

#: (key, numerator, denominator, description). The first two answer the
#: question and are the only ones drawn.
CONTRASTS = [
    ('elimination', 'E_GR', 'R_GR',
     'eliminated vs retained, both germline-restricted'),
    ('restriction', 'R_GR', 'R_S',
     'germline-restricted vs somatic peak, retained'),
    ('germline_broad', 'R_Gbroad', 'R_S',
     'germline peak detected in soma vs somatic peak, retained'),
    ('eliminated_split', 'E_GR', 'E_other',
     'germline-restricted vs other, eliminated'),
    ('all_expressed', 'E_expr', 'R_expr',
     'every expressed gene, eliminated vs retained'),
]
DRAWN = ['elimination', 'restriction']

#: Model family: (name, continuous predictors, binary predictors).
MODEL_SPECS = [
    ('base', ['log_expr'], ['is_eliminated']),
    ('breadth', ['log_expr', 'tau'], ['is_eliminated', 'germline_peak']),
    ('breadth_restricted', ['log_expr', 'tau'],
     ['is_eliminated', 'germline_peak', 'germline_only']),
]


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def gff_rank(gff_path):
    """Rank and Identity of each model's best alignment, keyed on model ID.

    Read separately from any join, so the Rank column is the same whichever
    join is used. A Rank above 1 marks a secondary alignment of its query
    protein.
    """
    best, _ = select_primary_alignment(parse_gff_models(gff_path))
    ids = best['gene_id'].astype(str)
    return dict(zip(ids, best['rank'])), dict(zip(ids, best['identity']))


def tau_index(tpm):
    """Tissue specificity tau (Yanai et al. 2005) on log2(TPM + 1).

    tau = sum_i (1 - x_i / max x) / (n - 1): 0 for equal expression in every
    sample, 1 for expression in one. The log transform is the form
    Kryuchkova-Mostacci and Robinson-Rechavi (2017) found most robust. NaN for
    a gene with a missing stage or no expression anywhere.
    """
    x = np.log2(np.clip(np.asarray(tpm, float), 0.0, None) + 1.0)
    mx = x.max(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        tau = (1.0 - x / mx[:, None]).sum(axis=1) / (x.shape[1] - 1)
    ok = np.isfinite(x).all(axis=1) & (mx > 0)
    return np.where(ok, tau, np.nan)


def stage_floors(frame, quantiles):
    """{quantile: Series(stage -> TPM)} over the rows of `frame`, and the
    number of reference genes with a value at every floor stage."""
    sub = frame[FLOOR_STAGES].apply(pd.to_numeric, errors='coerce').dropna()
    out = {q: sub.quantile(q) for q in quantiles} if len(sub) else {}
    return out, len(sub)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def restricted_by_floor(frame, floor_vec):
    """(restricted, classifiable): every floor stage <= its floor, and every
    floor stage known. A gene missing a stage is not classified, not passed."""
    S = frame[FLOOR_STAGES]
    classifiable = S.notna().all(axis=1)
    below = S.le(floor_vec[FLOOR_STAGES]).all(axis=1)
    return below & classifiable, classifiable


def classify(frame, restricted, classifiable, germline_stages):
    """Fine group label per gene: {R,E}_{S,GR,Gbroad,low} or 'unclassified'.

    Unclassified: no stage of peak expression (not in the workbook), a germline
    peak at a stage left out by --germline-stages, or a germline peak with a
    somatic stage missing.
    """
    pre = pd.Series(np.where(frame['region1'].eq('eliminated'), 'E', 'R'),
                    index=frame.index)
    peak = frame[PEAK_COL]
    gpk = peak.isin(germline_stages)
    lab = pd.Series('unclassified', index=frame.index, dtype=object)
    lab[peak.isin(SOMATIC)] = pre + '_S'
    lab[peak.eq(LOW)] = pre + '_low'
    lab[gpk & classifiable & restricted] = pre + '_GR'
    lab[gpk & classifiable & ~restricted] = pre + '_Gbroad'
    return lab


def group_mask(frame, labels, group):
    """Genes of `group` with a dN/dS estimate."""
    return labels.isin(GROUP_MEMBERS[group]) & frame['dN_dS'].notna()


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def boot_medians(v, rng, n_boot, chunk=250):
    """Medians of `n_boot` resamples of `v`, drawn in chunks to bound memory
    (2,000 x 6,000 indices at once would be ~100 MB)."""
    v = np.asarray(v, float)
    out = np.empty(n_boot)
    for s in range(0, n_boot, chunk):
        k = min(chunk, n_boot - s)
        out[s:s + k] = np.median(v[rng.integers(0, len(v), size=(k, len(v)))],
                                 axis=1)
    return out


def ratio_contrast(a, b, rng, n_boot, min_n, level=0.95):
    """Ratio of medians a/b, percentile bootstrap interval, Mann-Whitney P.

    The groups are resampled independently. A replicate whose denominator
    median is 0 has no ratio and is dropped; the count is returned, since
    dN/dS = 0 is common enough in a small group for it to happen.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    out = {'n_num': len(a), 'n_den': len(b),
           'median_num': float(np.median(a)) if len(a) else np.nan,
           'median_den': float(np.median(b)) if len(b) else np.nan,
           'fold': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan,
           'U': np.nan, 'p': np.nan, 'n_boot_dropped': 0, 'note': ''}
    if len(a) < min_n or len(b) < min_n:
        out['note'] = f'fewer than {min_n} genes in a group; not tested'
        return out
    if out['median_den'] > 0:
        out['fold'] = out['median_num'] / out['median_den']
    ba, bb = boot_medians(a, rng, n_boot), boot_medians(b, rng, n_boot)
    ok = bb > 0
    out['n_boot_dropped'] = int((~ok).sum())
    if ok.sum() >= 0.9 * n_boot:
        alpha = (1 - level) / 2
        out['ci_lo'], out['ci_hi'] = (float(x) for x in
                                      np.quantile(ba[ok] / bb[ok],
                                                  [alpha, 1 - alpha]))
    else:
        out['note'] = 'denominator median 0 in >10% of replicates; no interval'
    u, p = mannwhitneyu(a, b, alternative='two-sided')
    out['U'], out['p'] = float(u), float(p)
    return out


def match_once(e_cands, n_ret, rng):
    """One matching pass: eliminated genes in random order, each given a
    retained partner drawn at random from its candidates not yet taken.
    Returns the partner index per eliminated gene (-1 for none)."""
    used = np.zeros(n_ret, dtype=bool)
    pick = np.full(len(e_cands), -1)
    for i in rng.permutation(len(e_cands)):
        c = e_cands[i]
        if len(c):
            c = c[~used[c]]
        if len(c):
            j = c[rng.integers(len(c))]
            used[j] = True
            pick[i] = j
    return pick


def candidate_lists(e_stage, e_x, r_stage, r_x, caliper):
    """Retained candidates for each eliminated gene: same peak stage, expression
    level within `caliper` (log2 units)."""
    return [np.flatnonzero((r_stage == s) & (np.abs(r_x - x) <= caliper))
            for s, x in zip(e_stage, e_x)]


def matched_contrast(e, r, caliper, n_rep, n_boot, rng, level=0.95):
    """Elimination contrast with stage-exact, level-calipered matching.

    `e`, `r`: frames of eliminated and retained germline-restricted genes with
    columns peak, x (log2 peak TPM + 1) and dN_dS. Returns (summary dict,
    per-replicate frame). Observed replicates (`n_rep`) give the point estimate
    and the spread due to matching alone; resample-then-match replicates
    (`n_boot`) give the interval.
    """
    # Integer stage codes: comparing object arrays of strings in the inner loop
    # is an order of magnitude slower than comparing integers.
    codes, _ = pd.factorize(pd.concat([e['peak'], r['peak']], ignore_index=True))
    e_stage, r_stage = codes[:len(e)], codes[len(e):]
    e_x, e_v = e['x'].to_numpy(float), e['dN_dS'].to_numpy(float)
    r_x, r_v = r['x'].to_numpy(float), r['dN_dS'].to_numpy(float)
    cands = candidate_lists(e_stage, e_x, r_stage, r_x, caliper)
    n_none = int(sum(len(c) == 0 for c in cands))
    rows = []
    for rep in range(n_rep):
        pick = match_once(cands, len(r_v), rng)
        ok = pick >= 0
        if ok.sum() < 5:
            continue
        me, mr = np.median(e_v[ok]), np.median(r_v[pick[ok]])
        rows.append({'kind': 'observed', 'replicate': rep,
                     'n_matched': int(ok.sum()),
                     'median_eliminated': me, 'median_matched_retained': mr,
                     'fold': me / mr if mr > 0 else np.nan,
                     'p': mannwhitneyu(e_v[ok], r_v[pick[ok]],
                                       alternative='two-sided')[1],
                     'level_balance_p': mannwhitneyu(
                         e_x[ok], r_x[pick[ok]], alternative='two-sided')[1]})
    for rep in range(n_boot):
        ie = rng.integers(0, len(e_v), len(e_v))
        ir = rng.integers(0, len(r_v), len(r_v))
        bc = candidate_lists(e_stage[ie], e_x[ie], r_stage[ir], r_x[ir], caliper)
        pick = match_once(bc, len(ir), rng)
        ok = pick >= 0
        if ok.sum() < 5:
            continue
        me, mr = np.median(e_v[ie][ok]), np.median(r_v[ir][pick[ok]])
        rows.append({'kind': 'bootstrap', 'replicate': rep,
                     'n_matched': int(ok.sum()), 'median_eliminated': me,
                     'median_matched_retained': mr,
                     'fold': me / mr if mr > 0 else np.nan,
                     'p': np.nan, 'level_balance_p': np.nan})
    reps = pd.DataFrame(rows)
    obs = reps[reps['kind'] == 'observed'] if len(reps) else reps
    bts = reps[reps['kind'] == 'bootstrap'] if len(reps) else reps
    alpha = (1 - level) / 2
    out = {'n_eliminated': len(e_v), 'n_retained_pool': len(r_v),
           'n_no_candidate': n_none, 'caliper_log2': caliper,
           'n_rep_observed': len(obs), 'n_rep_bootstrap': len(bts),
           'fold': np.nan, 'match_lo': np.nan, 'match_hi': np.nan,
           'ci_lo': np.nan, 'ci_hi': np.nan, 'median_n_matched': np.nan,
           'pct_rep_p_lt_0.05': np.nan, 'median_p': np.nan,
           'median_level_balance_p': np.nan}
    if len(obs):
        f = obs['fold'].dropna()
        out.update(fold=float(f.median()),
                   match_lo=float(f.quantile(alpha)),
                   match_hi=float(f.quantile(1 - alpha)),
                   median_n_matched=float(obs['n_matched'].median()),
                   **{'pct_rep_p_lt_0.05': float(100 * (obs['p'] < 0.05).mean())},
                   median_p=float(obs['p'].median()),
                   median_level_balance_p=float(obs['level_balance_p'].median()))
    if len(bts):
        f = bts['fold'].dropna()
        if len(f):
            out.update(ci_lo=float(f.quantile(alpha)),
                       ci_hi=float(f.quantile(1 - alpha)))
    return out, reps


def fit_models(frame, label):
    """The MODEL_SPECS family on raw and ranked dN/dS. Tidy coefficient rows."""
    cols = ['dN_dS', 'log_expr', 'tau', 'is_eliminated', 'germline_peak',
            'germline_only']
    mf = frame[cols].dropna().copy()
    mf['dN_dS_rank'] = rankdata(mf['dN_dS'])
    rows = []
    for response in ('dN_dS', 'dN_dS_rank'):
        for name, cont, binary in MODEL_SPECS:
            m = fit_model(mf, response, cont, binary,
                          f'{label}:{response}:{name}')
            if m is None:
                continue
            for nm, b, se, lo, hi, p in zip(m['names'], m['beta'], m['se'],
                                            m['ci_lo'], m['ci_hi'], m['p']):
                rows.append({'subset': label, 'response': response,
                             'model': name, 'term': nm, 'beta': b, 'se': se,
                             'ci_lo': lo, 'ci_hi': hi, 'p': p, 'n': m['n'],
                             'n_retained': m['n_retained'],
                             'n_eliminated': m['n_eliminated'], 'r2': m['r2']})
    return pd.DataFrame(rows)


def concordance_table(df, have_rank):
    """region1 x workbook elimination call, for genes with a workbook match.

    Per cell: peak-stage mix, the share peaking in adult somatic tissue (where
    an eliminated gene cannot be expressed), the share of primary (Rank 1)
    alignments and median dN/dS.
    """
    rows = []
    m_wb = df['wb_eliminated'].notna()
    for reg in ('retained', 'eliminated'):
        for call in (False, True):
            sub = df[m_wb & df['region1'].eq(reg)
                     & df['wb_eliminated'].eq(call)]
            pk = sub[PEAK_COL].dropna()
            n_pk = len(pk)

            def pct(mask_n):
                return 100.0 * mask_n / n_pk if n_pk else np.nan
            row = {'region1': reg,
                   'workbook_call': 'eliminated' if call else 'retained',
                   'concordant': (reg == 'eliminated') == call,
                   'n': len(sub), 'n_with_peak': n_pk,
                   'pct_germline_peak': pct(pk.isin(GERMLINE).sum()),
                   'pct_somatic_peak': pct(pk.isin(SOMATIC).sum()),
                   'pct_adult_somatic_peak': pct(pk.isin(ADULT_SOMATIC).sum()),
                   'pct_low_expression': pct(pk.eq(LOW).sum()),
                   'median_dN_dS': float(sub['dN_dS'].median())
                   if sub['dN_dS'].notna().any() else np.nan,
                   'median_somatic_max_tpm': float(sub['somatic_max'].median())
                   if sub['somatic_max'].notna().any() else np.nan}
            if have_rank:
                rk = sub['gff_rank'].dropna()
                row['n_with_rank'] = len(rk)
                row['pct_rank1'] = 100.0 * (rk == 1).mean() if len(rk) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def finish_figure(fig, outdir, stem):
    """Audit, then save in every format through paper_colors.save()."""
    problems = list(audit(fig, verbose=False)) + offcanvas_text(fig)
    for prob in problems:
        logging.warning('  audit %s: %s', stem, prob)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths, problems


def draw_groups(ax, df, labels, rng, cap_pct, jitter_max, min_box):
    """S12B: dN/dS by group, retained left, eliminated right.

    Boxes are quartiles of every gene in the group, whiskers the 5th-95th
    percentile, no fliers; ax.boxplot is handed every value, so the box does not
    depend on what is drawn. At most `jitter_max` points per group are jittered
    on top. The y-axis stops at the `cap_pct` percentile of all values in the
    drawn groups, widened to every whisker end; jittered points above it are
    open triangles on the top edge, counted per group.
    """
    xs = {'R_S': 0.0, 'R_Gbroad': 1.0, 'R_GR': 2.0, 'E_GR': 3.35, 'E_other': 4.35}
    vals = {g: df.loc[group_mask(df, labels, g), 'dN_dS'].to_numpy(float)
            for g in GROUP_ORDER}
    if not any(len(v) for v in vals.values()):
        sys.exit('ERROR: S12B has no gene in any group; check the join counts '
                 'in the log.')
    pooled = np.concatenate([v for v in vals.values() if len(v)])
    ycap = float(np.percentile(pooled, cap_pct))
    whisk = [np.percentile(v, 95) for v in vals.values() if len(v) >= min_box]
    ycap = max([ycap] + whisk)
    y_lo, y_hi = -0.03 * ycap, 1.05 * ycap
    info = {}
    for g in GROUP_ORDER:
        v, x, reg = vals[g], xs[g], GROUP_REGION[g]
        rec = {'n': len(v), 'n_jitter_drawn': 0, 'n_above_cap': int((v > ycap).sum()),
               'n_triangles_drawn': 0, 'drawn_box': len(v) >= min_box}
        if len(v) >= min_box:
            bp = ax.boxplot(v, positions=[x], widths=0.56, patch_artist=True,
                            showfliers=False, whis=(5, 95), zorder=2)
            style_box(bp, REGION[reg], 0.40)
        if len(v):
            show = v if len(v) <= jitter_max else rng.choice(v, jitter_max,
                                                             replace=False)
            inside, above = show[show <= ycap], show[show > ycap]
            ax.scatter(x + rng.uniform(-0.19, 0.19, len(inside)), inside,
                       s=1.3, linewidths=0, color=to_rgba(REGION[reg], 0.35),
                       zorder=3, rasterized=True)
            if len(above):
                ax.scatter(x + rng.uniform(-0.19, 0.19, len(above)),
                           np.full(len(above), y_hi), marker='^', s=6,
                           facecolors='none', edgecolors=REGION[reg],
                           linewidths=0.4, zorder=4, clip_on=False)
            rec.update(n_jitter_drawn=int(len(inside)),
                       n_triangles_drawn=int(len(above)))
        info[g] = rec
    ax.set_xticks([xs[g] for g in GROUP_ORDER])
    ax.set_xticklabels([f'{GROUP_TICK[g]}\nn = {len(vals[g]):,}'
                        for g in GROUP_ORDER], linespacing=1.05)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.55, 4.9)
    ax.set_ylim(y_lo, y_hi)
    ax.set_ylabel('dN/dS')
    handles = [Patch(facecolor=to_rgba(REGION[r], 0.40), edgecolor='black',
                     linewidth=LINEWIDTH['thin']) for r in ('retained', 'eliminated')]
    ax.legend(handles, ['Retained', 'Eliminated'], loc='lower left',
              bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False,
              fontsize=FONT['legend'], borderaxespad=0.4, handlelength=1.2,
              handletextpad=0.4, columnspacing=1.2)
    return {'y_cap': ycap, 'cap_pct': cap_pct, 'groups': info,
            'jitter_max': jitter_max}


def draw_forest(ax, rows, key_ncol=1, ax_t=None):
    """S12C: elimination and restriction folds at every definition.

    `rows`: dicts with label, block, and per contrast (fold, lo, hi), plus the n
    of the elimination contrast. Filled circles: eliminated / retained, both
    germline-restricted. Open squares: germline-restricted / somatic peak,
    retained. Log2 axis, dashed rule at 1; an interval running off the axis
    ends in an arrowhead. With `ax_t`, the n column is written in that table
    axes (as 6C's columns) instead of hung off the plot's right edge.
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
            for v in (r[c][1], r[c][2]) if v is not None and np.isfinite(v) and v > 0]
    # Limits snap outward to the nearest tick, so the axis spans the intervals
    # and no more; ends beyond 1/8 or 8 get arrowheads.
    nice = [0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
    lo, hi = min(ends + [0.9]), max(ends + [1.1])
    x_lo = max([t for t in nice if t <= lo] or [nice[0]])
    x_hi = min([t for t in nice if t >= hi] or [nice[-1]])
    # (vertical offset within the row, marker, face colour) per contrast
    style = {'elimination': (+0.17, 'o', NEUTRAL['ink']),
             'restriction': (-0.17, 's', 'white')}
    n_clipped = 0
    for r, yy in zip(rows, ys):
        for c in DRAWN:
            v = r.get(c)
            if not v or v[0] is None or not np.isfinite(v[0]):
                continue
            fold, a, b = v
            dy, marker, face = style[c]
            if a is not None and np.isfinite(a) and np.isfinite(b):
                ax.plot([max(a, x_lo), min(b, x_hi)], [yy + dy] * 2,
                        color=NEUTRAL['ink'], lw=LINEWIDTH['thin'],
                        solid_capstyle='butt', zorder=2)
                for end, lim, mk in ((a, x_lo, '<'), (b, x_hi, '>')):
                    if (mk == '<' and end < lim) or (mk == '>' and end > lim):
                        ax.plot([lim], [yy + dy], marker=mk, ms=2.8,
                                color=NEUTRAL['ink'], zorder=3, clip_on=False)
                        n_clipped += 1
            ax.plot([min(max(fold, x_lo), x_hi)], [yy + dy], ls='none',
                    marker=marker, ms=3.4, mfc=face, mec=NEUTRAL['ink'],
                    mew=0.7, zorder=4)
        if ax_t is not None:
            side_text(ax_t, yy, r['n_text'], dx_pt=0.0, x=0.0)
        else:
            side_text(ax, yy, r['n_text'])
    if ax_t is not None:
        side_text(ax_t, ys[0] + 0.62, 'n', dx_pt=0.0, x=0.0, va='bottom')
    else:
        side_text(ax, ys[0] + 0.62, 'n', va='bottom')
    # The rule at 1 in data coordinates, not axvline: axvline's y runs 0-1 in
    # axes units, which audit() reads as data units and so reports the rule as
    # running up into the legend above the panel.
    y_lo, y_hi = ys[-1] - 0.7, ys[0] + 0.7
    ax.plot([1.0, 1.0], [y_lo, y_hi], color=NEUTRAL['rule'],
            lw=LINEWIDTH['rule'], ls=(0, (4, 2)), zorder=1)
    ax.set_xscale('log', base=2)
    ax.set_xlim(x_lo, x_hi)
    ticks = [t for t in nice if x_lo <= t <= x_hi]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f'{t:g}' for t in ticks])
    ax.set_yticks(ys)
    ax.set_yticklabels([r['label'] for r in rows])
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(y_lo, y_hi)
    if ax_t is not None:
        ax_t.set_ylim(y_lo, y_hi)
    ax.set_xlabel('dN/dS fold change (ratio of medians)')
    h = [Line2D([], [], ls='none', marker='o', ms=3.4, mfc=NEUTRAL['ink'],
                mec=NEUTRAL['ink'], mew=0.7),
         Line2D([], [], ls='none', marker='s', ms=3.4, mfc='white',
                mec=NEUTRAL['ink'], mew=0.7)]
    ax.legend(h, ['Eliminated / retained, germline only',
                  'Germline only / somatic peak, retained'],
              loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=key_ncol,
              frameon=False, fontsize=FONT['legend'], borderaxespad=0.4,
              handletextpad=0.3, labelspacing=0.3, columnspacing=1.6)
    return {'x_range': (x_lo, x_hi), 'n_interval_ends_clipped': n_clipped}


def forest_span(rows, block_gap=0.45, pad=0.7):
    """Height of the forest in row units, as draw_forest lays it out: one per
    row, `block_gap` between blocks, `pad` above the first and below the last."""
    gaps = sum(1 for a, b in zip(rows, rows[1:]) if a['block'] != b['block'])
    return (len(rows) - 1) + block_gap * gaps + 2 * pad


def contrast_height_mm(rows, key_rows):
    """Canvas height for S12C: 4.2 mm per row unit, as 6C, plus the key, the
    x-axis and its label."""
    return float(min(150.0, 17.0 + 3.6 * key_rows + 4.2 * forest_span(rows)))


def plot_groups_panel(df, labels, outdir, stem, rng, args):
    """S12B: dN/dS by expression group, one 89 mm column as tall as S12A.

    The margins are starting values; fit_to_canvas measures the tick labels,
    the key above the axes and the y label, and moves only the edges that
    overflow, so the panel fills the column in whatever font the node has.
    """
    fig = blank_figure('single', args.b_height_mm)
    ax = add_axes_mm(fig, left=13.0, bottom=14.0, right=1.5, top=8.0)
    info = draw_groups(ax, df, labels, rng, args.cap_pct, args.jitter_max,
                       args.min_box)
    fit_to_canvas(fig, ax)
    paths, problems = finish_figure(fig, outdir, stem)
    return paths, problems, info


def plot_contrast_panel(rows, outdir, stem, args):
    """S12C: the elimination and restriction contrasts at every definition.

    At 136 or 185 mm the key fits on one row; at 89 mm it takes two. The height
    follows the number of rows (contrast_height_mm) unless --c-height-mm is
    given. The n column is a table axes beside the plot, and place_with_table
    sets both boxes from the drawn text, as in 6C.
    """
    key_ncol = 1 if args.c_width == 'single' else 2
    height = args.c_height_mm or contrast_height_mm(rows, 2 // key_ncol)
    width_mm = {'single': 89.0, 'onehalf': 136.0, 'double': 185.0}[args.c_width]
    fig = blank_figure(args.c_width, height)
    # Starting boxes only: place_with_table sets both from the drawn text.
    top = 4.0 + 3.6 * (2 // key_ncol)
    ax = add_axes_mm(fig, left=26.0, bottom=11.0, right=17.0, top=top)
    ax_t = add_axes_mm(fig, left=width_mm - 16.0, bottom=11.0, right=1.0, top=top)
    ax_t.axis('off')
    info = draw_forest(ax, rows, key_ncol=key_ncol, ax_t=ax_t)
    place_with_table(fig, ax, ax_t)
    paths, problems = finish_figure(fig, outdir, stem)
    info['height_mm'] = height
    return paths, problems, info


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_float_list(s, name, lo=0.0, hi=1.0):
    try:
        vals = [float(x) for x in str(s).split(',') if x.strip()]
    except ValueError:
        sys.exit(f'ERROR: --{name} must be comma-separated numbers, got {s!r}')
    bad = [v for v in vals if not lo < v < hi]
    if not vals or bad:
        sys.exit(f'ERROR: --{name} values must lie in ({lo:g}, {hi:g}); got {s!r}')
    return vals


def parse_args():
    p = argparse.ArgumentParser(
        description='dN/dS of eliminated vs retained germline-restricted genes '
                    '(reviewer control for Figure 6).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--kaks', required=True,
                   help='Annotated rate table (ANNOTATED_TABLE, from '
                        'annotate_regions.py)')
    p.add_argument('--annotation', required=True,
                   help='Ascaris_genes_15714_annotation_expression.xlsx')
    p.add_argument('--bridge', default=None,
                   help='workbook_model_bridge.tsv (workbook_model_bridge.py): '
                        'model ID -> workbook gene, by position, as Figure 6. '
                        'The join to use')
    p.add_argument('--name-join', action='store_true',
                   help='Join by the GFF Target= name instead, as the '
                        '2026-10-16 run did. Wrong (results draft note 60); '
                        'only to reproduce that run')
    p.add_argument('--gff1', required=True,
                   help='Ascaris miniprot GFF: each model\'s miniprot Rank for '
                        'the concordance table, and the names for --name-join')
    p.add_argument('--gene-col', default='gene1')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--wb-elim-col', default='eliminated',
                   help="Workbook column with its own elimination call "
                        "(default 'eliminated'); '' to run without it, which "
                        "drops the concordance analyses and needs "
                        "--floor-source analysed")
    p.add_argument('--wb-elim-pattern', default=r'^eliminated',
                   help='Regex marking an eliminated gene in --wb-elim-col '
                        "(default '^eliminated'; retained genes read '.')")
    p.add_argument('--germline-stages', default=','.join(GERMLINE),
                   help='Peak stages that count as germline (default: '
                        f"{','.join(GERMLINE)}, as Figure 6). Dropping "
                        'early_embryo tests whether pre-elimination embryos '
                        'drive the result')
    p.add_argument('--floor-source', choices=['workbook', 'analysed', 'concordant'],
                   default='workbook',
                   help='Genes that set the somatic noise floor: every gene the '
                        'workbook calls eliminated (default; independent of the '
                        'rate table), the analysed eliminated genes, or those '
                        'both annotations call eliminated')
    p.add_argument('--floor-quantiles', default='0.95,0.90,0.75',
                   help='Floor quantiles, one definition each (default '
                        '0.95,0.90,0.75; the first is the default definition)')
    p.add_argument('--tau-thresholds', default='0.7,0.8,0.9',
                   help='tau thresholds, one definition each (default 0.7,0.8,0.9)')
    p.add_argument('--caliper', type=float, default=0.25,
                   help='Matching caliper on log2(peak TPM + 1) (default 0.25)')
    p.add_argument('--n-match', type=int, default=500,
                   help='Matching replicates on the observed genes (default 500)')
    p.add_argument('--n-match-boot', type=int, default=1000,
                   help='Resample-then-match replicates for the matched '
                        'interval (default 1000)')
    p.add_argument('--n-boot', type=int, default=2000,
                   help='Bootstrap replicates per contrast (default 2000)')
    p.add_argument('--min-group', type=int, default=10,
                   help='A contrast is tested only if both groups hold at least '
                        'this many genes (default 10)')
    p.add_argument('--min-box', type=int, default=5,
                   help='S12B draws a box for a group of at least this many '
                        'genes (default 5)')
    p.add_argument('--cap-pct', type=float, default=99.0,
                   help='S12B y-axis stops at this percentile, widened to '
                        'every whisker end; drawing only (default 99)')
    p.add_argument('--jitter-max', type=int, default=400,
                   help='S12B jitters at most this many points per group; '
                        'boxes and tests use every gene (default 400)')
    p.add_argument('--b-height-mm', type=float, default=72.0,
                   help='Height of S12B, to match S12A beside it (default 72, '
                        "functional_enrichment.py's S6 height; the wrapper reads "
                        'it from the newest S6_dnds_vs_expression.svg)')
    p.add_argument('--c-width', choices=['single', 'onehalf', 'double'],
                   default='single',
                   help='Width of S12C: 89, 136 or 185 mm (default single, '
                        'beside S12D on the bottom row)')
    p.add_argument('--c-height-mm', type=float, default=None,
                   help='Height of S12C (default: from the number of rows)')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None,
                   help='Default: results/germline_restricted_dnds_<YYYYmmdd_HHMMSS>')
    args = p.parse_args()
    if bool(args.bridge) == bool(args.name_join):
        sys.exit('ERROR: give --bridge workbook_model_bridge.tsv (the position '
                 'join), or --name-join to reproduce the superseded name join; '
                 'exactly one of the two.')
    args.floor_q = parse_float_list(args.floor_quantiles, 'floor-quantiles')
    args.tau_t = parse_float_list(args.tau_thresholds, 'tau-thresholds')
    args.germ = [s.strip() for s in args.germline_stages.split(',') if s.strip()]
    bad = [s for s in args.germ if s not in GERMLINE]
    if not args.germ or bad:
        sys.exit(f'ERROR: --germline-stages must be a subset of {GERMLINE}; '
                 f'got {args.germline_stages!r}')
    if args.caliper <= 0:
        sys.exit('ERROR: --caliper must be positive')
    if not 40.0 <= args.b_height_mm <= 120.0:
        sys.exit(f'ERROR: --b-height-mm {args.b_height_mm:g} is outside 40-120 mm')
    if args.c_height_mm is not None and not 30.0 <= args.c_height_mm <= 150.0:
        sys.exit(f'ERROR: --c-height-mm {args.c_height_mm:g} is outside 30-150 mm')
    if not args.wb_elim_col and args.floor_source != 'analysed':
        sys.exit("ERROR: without --wb-elim-col the floor can only come from the "
                 "analysed genes; add --floor-source analysed")
    if args.outdir is None:
        # make_outdir with a seconds stamp: two runs on one day must not share
        # a folder (script_to_figure_map, global rules).
        args.outdir = make_outdir('germline_restricted_dnds', base='results',
                                  date=datetime.now().strftime('%Y%m%d_%H%M%S'))
    return args


def counts_by_region(df, mask):
    return (int((mask & df['region1'].eq('retained')).sum()),
            int((mask & df['region1'].eq('eliminated')).sum()))


def fmt_fold(r):
    """'1.80 [1.52, 2.11]' or 'NA (reason)'."""
    if r is None or not np.isfinite(r.get('fold', np.nan)):
        return 'NA' + (f" ({r['note']})" if r is not None and r.get('note') else '')
    ci = (f" [{r['ci_lo']:.2f}, {r['ci_hi']:.2f}]"
          if np.isfinite(r.get('ci_lo', np.nan)) else ' [no interval]')
    return f"{r['fold']:.2f}{ci}"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'germline_restricted_dnds')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('germline_restricted_dnds.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()
    # Separate streams, so a change to the figure's jitter can never move a
    # statistic, and the matching does not depend on how many bootstraps ran.
    rng_boot, rng_match, rng_fig = (np.random.default_rng(s) for s in
                                    np.random.SeedSequence(args.seed).spawn(3))

    # --- rate table -----------------------------------------------------------
    rates = read_rate_table(args.kaks, require=('dN_dS',))
    for col in (args.gene_col, args.region_col):
        if col not in rates.columns:
            sys.exit(f'ERROR: column {col!r} not in {args.kaks}\n'
                     f'  columns: {list(rates.columns)}')
    rates = rates.rename(columns={args.gene_col: 'gene1',
                                  args.region_col: 'region1'})
    rates['gene1'] = rates['gene1'].astype(str)
    for c in ('dN', 'dS', 'dN_dS'):
        if c in rates.columns:
            rates[c] = pd.to_numeric(rates[c], errors='coerce')
    n_rows = len(rates)
    df = rates[rates['region1'].isin(['retained', 'eliminated'])].copy()
    logging.info('Rate table: %s rows, %s retained or eliminated',
                 f'{n_rows:,}', f'{len(df):,}')
    if df.empty or df['region1'].nunique() < 2:
        sys.exit('ERROR: the rate table needs both retained and eliminated rows.')

    # --- join, as Figure 6 ---------------------------------------------------
    align_diag, bridge_info = None, None
    if args.bridge:
        name_of, bridge_info = read_bridge(args.bridge)
        df['annot_name'] = df['gene1'].map(name_of)
        join_step, join_col = 'given a workbook gene (bridge)', 'annot_name'
    else:
        logging.warning('--name-join: the superseded Target-name join, wrong '
                        'for most genes (results draft note 60)')
        name_of, target_of, align_diag = target_bridge(args.gff1)
        df['target'] = df['gene1'].map(target_of)
        df['annot_name'] = df['gene1'].map(name_of)
        join_step, join_col = 'with a GFF Target', 'target'
    rank_of, ident_of = gff_rank(args.gff1)
    df['gff_rank'] = pd.to_numeric(df['gene1'].map(rank_of), errors='coerce')
    df['gff_identity'] = pd.to_numeric(df['gene1'].map(ident_of), errors='coerce')
    have_rank = bool(df['gff_rank'].notna().any())
    if not have_rank:
        logging.warning('  the GFF carries no Rank= attribute; the concordance '
                        'table has no rank column')

    wb = load_workbook(args.annotation)
    have_wb_elim = bool(args.wb_elim_col)
    if have_wb_elim:
        if args.wb_elim_col not in wb.columns:
            sys.exit(f'ERROR: workbook has no column {args.wb_elim_col!r}; '
                     f"pass --wb-elim-col '' --floor-source analysed to run "
                     f'without it.\n  columns: {list(wb.columns)}')
        lab = wb[args.wb_elim_col].astype(str).str.strip()
        wb['wb_eliminated_label'] = lab
        wb['wb_eliminated'] = lab.str.contains(args.wb_elim_pattern, regex=True)
        logging.info('  workbook elimination call: %s eliminated, %s not '
                     "(pattern %r on %r)", f"{int(wb['wb_eliminated'].sum()):,}",
                     f"{int((~wb['wb_eliminated']).sum()):,}",
                     args.wb_elim_pattern, args.wb_elim_col)
        if not wb['wb_eliminated'].any():
            sys.exit(f'ERROR: no workbook gene matches {args.wb_elim_pattern!r} '
                     f'in {args.wb_elim_col!r}; values: '
                     f'{sorted(lab.unique())[:8]}')
    keep = ['name', PEAK_COL] + STAGES + \
        (['wb_eliminated', 'wb_eliminated_label'] if have_wb_elim else [])
    df = df.merge(wb[keep], left_on='annot_name', right_on='name', how='left')
    if not have_wb_elim:
        df['wb_eliminated'] = np.nan
        df['wb_eliminated_label'] = np.nan

    # --- per-gene expression features ------------------------------------------
    for s in STAGES:
        df[s] = pd.to_numeric(df[s], errors='coerce')
    df['category'] = df[PEAK_COL].map(CATEGORY_OF)
    df['germline_max'] = df[args.germ].max(axis=1, skipna=False)
    df['somatic_max'] = df[FLOOR_STAGES].max(axis=1, skipna=False)
    df['adult_somatic_max'] = df[ADULT_SOMATIC].max(axis=1, skipna=False)
    df['log_expr'] = np.log2(df[STAGES].max(axis=1, skipna=False) + 1.0)
    df['tau'] = tau_index(df[STAGES].to_numpy(float))
    df['is_eliminated'] = df['region1'].eq('eliminated').astype(int)
    df['germline_peak'] = df[PEAK_COL].isin(args.germ).astype(int)
    df['expressed'] = df[PEAK_COL].isin(args.germ + SOMATIC)
    df['concordant'] = (df['wb_eliminated'].notna()
                        & (df['wb_eliminated'].astype(object) ==
                           df['region1'].eq('eliminated')))

    n_steps = {
        'analysed (retained or eliminated)': counts_by_region(df, df['gene1'].notna()),
        join_step: counts_by_region(df, df[join_col].notna()),
        'matched to the workbook': counts_by_region(df, df['name'].notna()),
        'with a stage of peak expression': counts_by_region(df, df['category'].notna()),
        'expressed (germline or somatic peak)': counts_by_region(df, df['expressed']),
        'expressed, every somatic stage known': counts_by_region(
            df, df['expressed'] & df[FLOOR_STAGES].notna().all(axis=1)),
        'expressed, with dN/dS': counts_by_region(
            df, df['expressed'] & df['dN_dS'].notna()),
    }
    if have_wb_elim:
        n_steps['matched, both elimination calls agree'] = counts_by_region(
            df, df['concordant'])
    for k, (r, e) in n_steps.items():
        logging.info('  %-40s retained %s, eliminated %s', k, f'{r:,}', f'{e:,}')
    if n_steps['matched to the workbook'] == (0, 0):
        sys.exit('ERROR: no rate-table gene matched a workbook name.\n'
                 f"  bridged names: {df['annot_name'].dropna().head(3).tolist()}\n"
                 f"  workbook:      {wb['name'].head(3).tolist()}")
    # Two analysed models can take one workbook gene; each keeps its values, as
    # in Figure 6, and the count goes in the summary.
    shared = shared_gene_counts(df['annot_name'], df['region1'])
    if sum(shared):
        logging.info('  %s analysed genes share their workbook gene with another '
                     'analysed gene (retained %s, eliminated %s)',
                     f'{sum(shared):,}', f'{shared[0]:,}', f'{shared[1]:,}')

    # --- somatic noise floors ---------------------------------------------------
    floor_frames = {'analysed': df[df['region1'].eq('eliminated')
                                   & df['name'].notna()]}
    if have_wb_elim:
        floor_frames['workbook'] = wb[wb['wb_eliminated']]
        floor_frames['concordant'] = df[df['region1'].eq('eliminated')
                                        & df['wb_eliminated'].eq(True)]
    floors, floor_n, floor_rows = {}, {}, []
    for src, frame in floor_frames.items():
        floors[src], floor_n[src] = stage_floors(frame, args.floor_q)
        for q, vec in floors[src].items():
            for s in FLOOR_STAGES:
                floor_rows.append({'source': src, 'quantile': q, 'stage': s,
                                   'floor_tpm': float(vec[s]),
                                   'n_reference': floor_n[src],
                                   'used': src == args.floor_source})
    pd.DataFrame(floor_rows).to_csv(os.path.join(args.outdir, 'somatic_floor.tsv'),
                                    sep='\t', index=False, float_format='%.6g')
    if floor_n.get(args.floor_source, 0) < 20:
        sys.exit(f'ERROR: only {floor_n.get(args.floor_source, 0)} reference genes '
                 f'for --floor-source {args.floor_source}; a floor needs at least 20')
    use_floor = floors[args.floor_source]
    for q in args.floor_q:
        logging.info('  floor (%s, q = %g): %s', args.floor_source, q,
                     ', '.join(f'{s} {use_floor[q][s]:.2f}' for s in FLOOR_STAGES))

    # --- definitions --------------------------------------------------------------
    q0 = args.floor_q[0]
    # Labels are short because they sit beside the S12C forest; the caption defines them.
    defs = [{'key': f'floor_q{q:g}', 'kind': 'floor', 'q': q, 'subset': 'all',
             'block': 'floor', 'label': f'Floor {round(100 * q)}th pct'}
            for q in args.floor_q]
    defs += [{'key': f'tau_{t:g}', 'kind': 'tau', 't': t, 'subset': 'all',
              'block': 'tau', 'label': f'τ ≥ {t:g}'} for t in args.tau_t]
    if have_wb_elim:
        defs.append({'key': f'floor_q{q0:g}_concordant', 'kind': 'floor', 'q': q0,
                     'subset': 'concordant', 'block': 'robust',
                     'label': 'Calls agree'})

    labels_by_def, group_rows, contrast_rows = {}, [], []
    for d in defs:
        frame = df if d['subset'] == 'all' else df[df['concordant']]
        if d['kind'] == 'floor':
            restricted, classifiable = restricted_by_floor(frame, use_floor[d['q']])
        else:
            classifiable = frame['tau'].notna()
            restricted = frame['tau'].ge(d['t'])
        labels = classify(frame, restricted, classifiable, args.germ)
        labels_by_def[d['key']] = labels
        for g in SUMMARY_GROUPS:
            m = group_mask(frame, labels, g)
            v = frame.loc[m, 'dN_dS'].to_numpy(float)
            lo = hi = np.nan
            if len(v) >= args.min_group:
                bm = boot_medians(v, rng_boot, args.n_boot)
                lo, hi = (float(x) for x in np.quantile(bm, [0.025, 0.975]))
            mix = frame.loc[m, PEAK_COL].value_counts()
            group_rows.append({
                'definition': d['key'], 'subset': d['subset'], 'group': g,
                'description': GROUP_DESC[g], 'n': len(v),
                'median': float(np.median(v)) if len(v) else np.nan,
                'median_ci_lo': lo, 'median_ci_hi': hi,
                'q1': float(np.percentile(v, 25)) if len(v) else np.nan,
                'q3': float(np.percentile(v, 75)) if len(v) else np.nan,
                'median_tau': float(frame.loc[m, 'tau'].median()) if len(v) else np.nan,
                'median_log2_expr': float(frame.loc[m, 'log_expr'].median())
                if len(v) else np.nan,
                'peak_stage_mix': '; '.join(f'{k} {n}' for k, n in mix.items())})
        for ckey, num, den, desc in CONTRASTS:
            a = frame.loc[group_mask(frame, labels, num), 'dN_dS'].to_numpy(float)
            b = frame.loc[group_mask(frame, labels, den), 'dN_dS'].to_numpy(float)
            res = ratio_contrast(a, b, rng_boot, args.n_boot, args.min_group)
            contrast_rows.append({'definition': d['key'], 'subset': d['subset'],
                                  'label': d['label'], 'contrast': ckey,
                                  'numerator': num, 'denominator': den,
                                  'description': desc, **res})
    groups_tab = pd.DataFrame(group_rows)
    contrasts = pd.DataFrame(contrast_rows)
    contrasts['q_bh'] = benjamini_hochberg(contrasts['p'].to_numpy(float))
    groups_tab.to_csv(os.path.join(args.outdir, 'group_summary.tsv'), sep='\t',
                      index=False, float_format='%.6g')
    contrasts.to_csv(os.path.join(args.outdir, 'contrasts.tsv'), sep='\t',
                     index=False, float_format='%.6g')

    # --- matched elimination contrast -------------------------------------------
    match_defs = [(f'floor_q{q0:g}', 'all', 'Matched')]
    if have_wb_elim:
        match_defs.append((f'floor_q{q0:g}_concordant', 'concordant',
                           'Agree, matched'))
    matched, match_reps = {}, []
    for key, subset, label in match_defs:
        frame = df if subset == 'all' else df[df['concordant']]
        lab = labels_by_def[key]

        def pick(group):
            m = group_mask(frame, lab, group)
            return pd.DataFrame({'peak': frame.loc[m, PEAK_COL].to_numpy(),
                                 'x': frame.loc[m, 'log_expr'].to_numpy(float),
                                 'dN_dS': frame.loc[m, 'dN_dS'].to_numpy(float)})
        e, r = pick('E_GR'), pick('R_GR')
        if len(e) < args.min_group or len(r) < args.min_group:
            logging.warning('  matched %s: %d eliminated, %d retained '
                            'germline-restricted genes; not run', key, len(e), len(r))
            continue
        res, reps = matched_contrast(e, r, args.caliper, args.n_match,
                                     args.n_match_boot, rng_match)
        res.update(definition=key, subset=subset, label=label)
        matched[key] = res
        reps.insert(0, 'definition', key)
        match_reps.append(reps)
        logging.info('  matched %s: fold %.3f [%.3f, %.3f], %s of %d matched '
                     '(median), %d without a candidate', key, res['fold'],
                     res['ci_lo'], res['ci_hi'], res['median_n_matched'],
                     res['n_eliminated'], res['n_no_candidate'])
    if match_reps:
        pd.concat(match_reps).to_csv(os.path.join(args.outdir,
                                                  'matched_contrast.tsv'),
                                     sep='\t', index=False, float_format='%.6g')

    # --- models ---------------------------------------------------------------------
    lab0 = labels_by_def[f'floor_q{q0:g}']
    df['group_default'] = lab0
    df['germline_only'] = lab0.isin(['R_GR', 'E_GR']).astype(int)
    model_frames = [('all', df[df['expressed']])]
    if have_wb_elim:
        model_frames.append(('concordant', df[df['expressed'] & df['concordant']]))
    models = pd.concat([fit_models(f, name) for name, f in model_frames],
                       ignore_index=True)
    models.to_csv(os.path.join(args.outdir, 'models.tsv'), sep='\t', index=False,
                  float_format='%.6g')

    # --- label check ------------------------------------------------------------------
    conc = concordance_table(df, have_rank) if have_wb_elim else pd.DataFrame()
    if len(conc):
        conc.to_csv(os.path.join(args.outdir, 'concordance.tsv'), sep='\t',
                    index=False, float_format='%.6g')

    # --- per-gene table -------------------------------------------------------------
    for d in defs:
        df[f"group_{d['key']}"] = labels_by_def[d['key']].reindex(df.index)
    out_cols = (['gene1'] + [c for c in ('gene2',) if c in df.columns] +
                ['region1'] + [c for c in ('dN', 'dS', 'dN_dS') if c in df.columns] +
                [c for c in ('target',) if c in df.columns] +
                ['annot_name', 'gff_rank', 'gff_identity',
                 'wb_eliminated_label', 'wb_eliminated', 'concordant', PEAK_COL,
                 'category'] + STAGES +
                ['germline_max', 'somatic_max', 'adult_somatic_max', 'log_expr',
                 'tau', 'group_default'] + [f"group_{d['key']}" for d in defs])
    df[out_cols].to_csv(os.path.join(args.outdir, 'genes_classified.tsv'),
                        sep='\t', index=False, float_format='%.6g')

    # --- figure --------------------------------------------------------------------
    def look(defkey, ckey):
        r = contrasts[(contrasts['definition'] == defkey)
                      & (contrasts['contrast'] == ckey)]
        return r.iloc[0].to_dict() if len(r) else None

    forest_rows = []
    for d in defs:
        el, rs = look(d['key'], 'elimination'), look(d['key'], 'restriction')
        forest_rows.append({
            'label': d['label'], 'block': d['block'], 'definition': d['key'],
            'elimination': (el['fold'], el['ci_lo'], el['ci_hi']) if el else None,
            'restriction': (rs['fold'], rs['ci_lo'], rs['ci_hi']) if rs else None,
            'n_text': f"{el['n_num']:,} / {el['n_den']:,}" if el else ''})
    for key, subset, label in match_defs:
        if key in matched:
            m = matched[key]
            forest_rows.append({
                'label': label, 'block': 'robust', 'definition': key + '_matched',
                'elimination': (m['fold'], m['ci_lo'], m['ci_hi']),
                'restriction': None,
                'n_text': f"{int(round(m['median_n_matched'])):,} pairs"})
    logging.info('S12B: dN/dS by expression group')
    _, b_problems, a_info = plot_groups_panel(
        df[df['expressed']], lab0[df['expressed']], args.outdir,
        'S12B_germline_groups', rng_fig, args)
    logging.info('S12C: contrasts by definition')
    _, c_problems, b_info = plot_contrast_panel(
        forest_rows, args.outdir, 'S12C_germline_contrasts', args)
    fig_problems = ([f'S12B: {x}' for x in b_problems]
                    + [f'S12C: {x}' for x in c_problems])

    # --- summary ---------------------------------------------------------------------
    if args.bridge:
        join_line = (f"Join         position bridge {os.path.basename(args.bridge)}"
                     f": {bridge_info['assigned']:,} of {bridge_info['models']:,} "
                     f"models given a workbook gene (as Figure 6)")
    else:
        join_line = (f"Join         SUPERSEDED name join (results draft note 60): "
                     f"GFF Target= of {os.path.basename(args.gff1)}, "
                     f"{align_diag['n_ids']:,} models")
    L = ['Germline-restricted expression and dN/dS, eliminated vs retained', '',
         f'Rate table   {os.path.basename(args.kaks)}, {n_rows:,} rows',
         f'Workbook     {os.path.basename(args.annotation)}, {len(wb):,} genes',
         join_line,
         f"Germline stages {', '.join(args.germ)}; somatic stages "
         f"{', '.join(FLOOR_STAGES)}", '',
         'Genes at each step                         retained  eliminated']
    for k, (r, e) in n_steps.items():
        L.append(f'  {k:<40} {r:>8,}  {e:>10,}')
    L.append(f'  Sharing their workbook gene with another analysed gene: retained '
             f'{shared[0]:,}, eliminated {shared[1]:,} (each keeps its values, '
             f'as in Figure 6).')

    has_peak = df[df['category'].notna()]
    L += ['', 'CHECK AGAINST FIGURE 6A (functional_enrichment.py run with the same '
          'join must print the same counts)',
          '  Superseded name join, as in the Results draft: 8,059 retained and 310 '
          'eliminated genes with a stage of peak expression;',
          '           germline 35.2% / 41.3%, somatic 35.2% / 33.5%, low '
          '29.7% / 25.2% (retained / eliminated).']
    pk_r, pk_e = n_steps['with a stage of peak expression']
    shares = []
    for cat in ('germline', 'somatic', 'low'):
        pr = 100 * (has_peak['region1'].eq('retained') & has_peak['category'].eq(cat)).sum() / max(pk_r, 1)
        pe = 100 * (has_peak['region1'].eq('eliminated') & has_peak['category'].eq(cat)).sum() / max(pk_e, 1)
        shares.append(f'{cat} {pr:.1f}% / {pe:.1f}%')
    L += [f'  this run: {pk_r:,} retained and {pk_e:,} eliminated; '
          + ', '.join(shares) + '.',
          '  (With --germline-stages other than the default, the germline share '
          'here still uses Figure 6\'s categories.)']

    if len(conc):
        L += ['', 'Elimination calls: rate table (region1, v50 BED) x workbook '
              '(its own call for the joined gene)',
              '  region1     workbook    n      w/ peak  germline  somatic  '
              'adult soma  low     rank 1   median dN/dS']
        for r in conc.itertuples():
            rk = (f'{r.pct_rank1:5.1f}%' if have_rank and np.isfinite(r.pct_rank1)
                  else '   NA ')
            L.append(f'  {r.region1:<11} {r.workbook_call:<10} {r.n:6,}  '
                     f'{r.n_with_peak:7,}  {r.pct_germline_peak:6.1f}%  '
                     f'{r.pct_somatic_peak:6.1f}%  {r.pct_adult_somatic_peak:8.1f}%  '
                     f'{r.pct_low_expression:5.1f}%  {rk}   {r.median_dN_dS:.4f}')
        disc = conc[(conc['region1'] == 'eliminated') & ~conc['concordant']]
        if len(disc) and disc['n'].iloc[0]:
            L.append(f"  Eliminated in region1 but retained in the workbook: "
                     f"{int(disc['n'].iloc[0]):,} genes; their expression belongs "
                     f"to a gene the workbook places in retained DNA.")
        L.append('  Percentages of genes with a peak stage; adult soma = muscle, '
                 'intestine or carcass.')

    L += ['', f'Somatic noise floor: {args.floor_source} genes '
          f'(n = {floor_n[args.floor_source]:,} with every somatic stage known); '
          'TPM at each quantile']
    L.append('  quantile  ' + '  '.join(f'{s:>11}' for s in FLOOR_STAGES))
    for q in args.floor_q:
        L.append(f'  {q:<8g}  ' + '  '.join(f'{use_floor[q][s]:11.2f}'
                                            for s in FLOOR_STAGES))
    other_src = [s for s in floors if s != args.floor_source]
    if other_src:
        L.append(f"  Other sources ({', '.join(other_src)}) are in somatic_floor.tsv.")

    L += ['', f'Groups (dN/dS; median [95% bootstrap CI], n), default definition '
          f'floor_q{q0:g}']
    g0 = groups_tab[groups_tab['definition'] == f'floor_q{q0:g}']
    for r in g0.itertuples():
        ci = (f'[{r.median_ci_lo:.3f}, {r.median_ci_hi:.3f}]'
              if np.isfinite(r.median_ci_lo) else '[NA]')
        L.append(f'  {r.group:<9} {r.median:.4f} {ci:<16} n = {r.n:<6,} '
                 f'tau {r.median_tau:.2f}  log2 expr {r.median_log2_expr:.2f}  '
                 f'({r.description})')
    L.append('  Peak-stage mix of each group is in group_summary.tsv.')

    L += ['', 'Contrasts (fold = ratio of medians [95% bootstrap CI]; Mann-Whitney '
          'P; BH q over every row)']
    for d in defs:
        L.append(f"  {d['key']} ({d['label']})")
        for ckey, *_ in CONTRASTS:
            r = look(d['key'], ckey)
            if r is None:
                continue
            L.append(f"    {ckey:<17} {r['numerator']:>8} / {r['denominator']:<8} "
                     f"n {r['n_num']:>5,} / {r['n_den']:<6,} "
                     f"{fmt_fold(r):<24} P = {r['p']:.3g}  q = {r['q_bh']:.3g}")

    if matched:
        L += ['', f'Matched elimination contrast (exact on peak stage, caliper '
              f'{args.caliper:g} on log2 peak TPM + 1, without replacement)']
        for key, m in matched.items():
            L += [f"  {key}: fold {m['fold']:.3f}; interval {m['ci_lo']:.3f}-"
                  f"{m['ci_hi']:.3f} ({m['n_rep_bootstrap']} resample-then-match "
                  f"replicates); matching alone {m['match_lo']:.3f}-{m['match_hi']:.3f} "
                  f"({m['n_rep_observed']} replicates)",
                  f"    {m['median_n_matched']:.0f} of {m['n_eliminated']} "
                  f"eliminated genes matched (median); {m['n_no_candidate']} have "
                  f"no candidate; pool {m['n_retained_pool']:,} retained",
                  f"    Mann-Whitney P < 0.05 in {m['pct_rep_p_lt_0.05']:.0f}% of "
                  f"replicates (median P {m['median_p']:.3g}); expression-level "
                  f"balance, median P {m['median_level_balance_p']:.3g}"]

    if len(models):
        L += ['', 'Models (standardized OLS, as Figure 6E; expressed genes; the '
              'elimination coefficient is in SD units of the response)',
              '  subset      response    model                beta_elim [95% CI]'
              '        P          other terms']
        for (sub, resp, mod), g in models.groupby(['subset', 'response', 'model'],
                                                  sort=False):
            el = g[g['term'] == 'is_eliminated']
            if not len(el):
                continue
            el = el.iloc[0]
            others = '; '.join(f"{t.term} {t.beta:+.3f} (P {t.p:.2g})"
                               for t in g.itertuples()
                               if t.term not in ('intercept', 'is_eliminated'))
            L.append(f'  {sub:<11} {resp:<11} {mod:<20} {el.beta:+.3f} '
                     f'[{el.ci_lo:+.3f}, {el.ci_hi:+.3f}]  {el.p:<10.3g} {others}')
        n0 = models.iloc[0]
        L.append(f"  n = {int(n0['n']):,} in the first model "
                 f"({int(n0['n_retained']):,} retained, "
                 f"{int(n0['n_eliminated']):,} eliminated); every row in models.tsv.")

    # Outcome flags: computed from the default definition, not written in advance.
    el0, rs0 = look(f'floor_q{q0:g}', 'elimination'), look(f'floor_q{q0:g}', 'restriction')

    def where(r):
        """Position of a 95% interval relative to 1."""
        if r is None or not np.isfinite(r.get('ci_lo', np.nan)):
            return 'none'
        if r['ci_lo'] > 1:
            return 'above'
        if r['ci_hi'] < 1:
            return 'below'
        return 'spans'
    w_el, w_rs = where(el0), where(rs0)
    reading = {
        ('above', 'spans'): 'elimination raises dN/dS among germline-restricted '
                            'genes; germline restriction alone does not',
        ('above', 'above'): 'elimination and germline restriction both raise dN/dS',
        ('spans', 'above'): 'germline restriction raises dN/dS among retained '
                            'genes, and eliminated genes are no higher than '
                            'restricted retained genes',
        ('below', 'above'): 'germline restriction raises dN/dS among retained '
                            'genes, and eliminated genes are lower than '
                            'restricted retained genes',
        ('spans', 'spans'): 'neither contrast differs from 1',
    }.get((w_el, w_rs), 'no single reading; see both intervals')
    L += ['', 'Outcome at the default definition (from the intervals above)',
          f'  elimination {fmt_fold(el0)}, interval {w_el} 1; restriction '
          f'{fmt_fold(rs0)}, interval {w_rs} 1',
          f'  reading: {reading}.']

    L += ['', 'Figure audit (paper_colors.audit, plus text past the canvas edge)',
          '  ' + ('clean' if not fig_problems else '; '.join(fig_problems)), '',
          'Legend notes (S12B and S12C; S12A is functional_enrichment.py\'s '
          'S6_dnds_vs_expression)', '',
          f"  Panels: S12B 89 x {args.b_height_mm:.1f} mm, beside S12A; S12C "
          f"{ {'single': 89, 'onehalf': 136, 'double': 185}[args.c_width] } x "
          f"{b_info['height_mm']:.1f} mm, "
          + ("bottom row, beside S12D (dnds_within_class.py, drawn at this "
             "height). " if args.c_width == 'single' else "below. ")
          + "Panel letters are added at assembly.",
          f"  S12B: default definition (floor at the {round(100 * q0)}th percentile of "
          f"{args.floor_source} genes, every somatic stage). Boxes are quartiles of "
          f"every gene, whiskers the 5th-95th percentile. y-axis stops at "
          f"{a_info['y_cap']:.3f} (the {a_info['cap_pct']:g}th percentile of the "
          f"drawn groups, widened to every whisker end).",
          '  Per group: n, jittered points drawn, genes above the axis (all), '
          'triangles drawn (from the jitter sample):']
    for g in GROUP_ORDER:
        rec = a_info['groups'][g]
        L.append(f"    {g:<9} n {rec['n']:>6,}  drawn {rec['n_jitter_drawn']:>4}  "
                 f"above axis {rec['n_above_cap']:>4}  triangles "
                 f"{rec['n_triangles_drawn']:>3}" +
                 ('' if rec['drawn_box'] else '  (too few for a box)'))
    L += [f"  S12C: folds with 95% bootstrap intervals ({args.n_boot:,} resamples per "
          f"group); matched rows: median over {args.n_match} matchings, interval "
          f"from {args.n_match_boot:,} resample-then-match replicates. The n column "
          f"is eliminated / retained germline-only genes (pairs for matched rows). "
          f"{b_info['n_interval_ends_clipped']} interval ends run off the axis "
          f"(arrowheads).", '']
    summary = '\n'.join(L)
    with open(os.path.join(args.outdir, 'germline_restricted_summary.txt'),
              'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(), 'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'scipy': scipy_version,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv, 'args': vars(args),
            'inputs': {k: fingerprint(v) for k, v in
                       (('kaks', args.kaks), ('annotation', args.annotation),
                        ('bridge', args.bridge), ('gff1', args.gff1))},
            'join': {'mode': 'position bridge' if args.bridge
                     else 'name (superseded)', 'bridge': bridge_info,
                     'alignment': align_diag,
                     'shared_workbook_genes': {'retained': shared[0],
                                               'eliminated': shared[1]}},
            'gene_counts': {k: {'retained': r, 'eliminated': e}
                            for k, (r, e) in n_steps.items()},
            'germline_stages': args.germ, 'floor_stages': FLOOR_STAGES,
            'floor_source': args.floor_source,
            'floors': {src: {str(q): {s: float(v[s]) for s in FLOOR_STAGES}
                             for q, v in fl.items()} for src, fl in floors.items()},
            'floor_reference_n': floor_n,
            'definitions': defs,
            'matched': matched,
            'figure_audit': fig_problems,
        }, fh, indent=2, default=str)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
