#!/usr/bin/env python3
"""
make_supplementary_tables.py -- assemble Supplementary Tables S1-S11 of the
dN/dS paper from the analysis outputs: one Excel workbook for the journal, one
TSV per table part for the repository, and the captions.

Every value in these tables was computed by an analysis script. This script
selects, orders, renames and documents columns, and checks that tables which
report the same quantity agree with each other and with the per-gene data
behind them. It computes values itself in four places only, each deterministic
and each verified against an analysis output before it is used:

  * S2B  the per-branch `masked` flag (dN/dS > cap, no substitution, or no
         estimate). Every n, median and P of S2A is recomputed from S2B with
         the same rule.
  * S4   GC3 over all codons, and the ENC-ceiling counts and tests, from the
         per-gene tables, with the call codon_usage_bias.py uses. The same call
         is first shown to reproduce that script's summary for ENC, GC3s,
         dENC and CAI, and its printed ceiling-excluded results.
  * S6B  Tajima ratios from the substitution counts at full precision. The
         source rounds them to two decimals, and its `fold_faster` column is
         1 / the smallest of three rounded ratios, which is not the all-sites
         value the text quotes; it is dropped. Each chi-square is checked
         against (n1 - n2)^2 / (n1 + n2).
  * S9C  the median insulation score of each region and their two-sided
         Mann-Whitney test at each Hi-C stage, from that run's per-gene table
         with the run's own core-set rule (labelled region, finite dN/dS, and
         insulation, expression and GC present). Each is checked against the
         value the run printed in its summary, and the gene counts against
         its fitted models.
  * percentages from proportions (x 100).


Tables (numbered by first citation in the text, 2026-10-12)
------------------------------------------------------------
If an edit to the text changes the order of first citation, edit NUMBER below
(and FIGURE, for the supplementary figures the captions cite) and rerun;
nothing else in this file depends on the numbers. Sheets and TSVs follow NUMBER.

  S1   rates_pairwise        Ascaris x P. univalens rates, regions, gene structure
  S2   rate_heterogeneity    A tree additivity, B Tajima tests, C clock LRT
  S3   omega_above_one       pairs with dN/dS > 1: A genes, B region test
  S4   codon_usage           A all genes, B ENC ceiling, C GC3s-matched
  S5   microsynteny          A counts, B tests, C threshold sweep
  S6   x_autosome            per chromosome: retained-gene rates, Hi-C dosage
  S7   branch_models         codeml branch dN, dS, dN/dS: A tests, B per quartet
  S8   toxocara_models       A gene-model quality, B strata, C regression
  S9   insulation_models     A elimination coefficient by compartment control,
                             B every coefficient, C by Hi-C stage
  S10  popgen_cohort         the population-sample germline screen
  S11  popgen_alternatives   alternatives to relaxed constraint (curated file,
                             plus the McDonald-Kreitman row from mk_summary)

Until 2026-10-12.3 the numbers were the frozen figure map's: rates_pairwise S1,
branch_models S2, omega_above_one S3, codon_usage S4, insulation_models S5,
rate_heterogeneity S6, toxocara_models S7, microsynteny S8, popgen_cohort S9,
popgen_alternatives S10, x_autosome S11.


Inputs
------
A manifest, `supplementary_table_sources.tsv` beside this script: key<TAB>path
per line. Relative paths resolve against --base-dir; ${VAR} expands from the
environment. `--source key=path` overrides one line. Every source a selected
table needs is checked before anything is read, and all problems are reported
together. `cohort_audit` is optional (a cross-check for S9 only).


Checks
------
S3B's odds ratios are recomputed from its counts (conditional MLE with an
exact interval, as Methods states and microsynteny.py computes it) to confirm
which estimator the file holds; the table always shows the file's values.

Each check is logged and listed in the summary; a failed check is a WARNING,
the tables are still written, and the exit status is 2 so that the batch log
shows it. Missing inputs or columns stop the run (exit 1) before anything is
written.


Outputs (in --outdir)
---------------------
  Supplementary_Tables.xlsx          Contents sheet, one sheet per table part
                                     (title and caption above the header row),
                                     and a Columns sheet defining every column
  tsv/Table_S<n>[<part>]_<key>.tsv   the same tables, machine-readable
  table_captions.txt                 titles and captions for the Supplementary
                                     Material document
  supplementary_tables_summary.txt   rows per table, every check, every column
                                     left out and why, numbers for the text
  make_supplementary_tables.log
  provenance.json                    inputs (with SHA-256), versions, checks


Usage
-----
  python make_supplementary_tables.py \\
    --manifest 05_analyses/supplementary_table_sources.tsv \\
    --base-dir $PROJECT_DIR \\
    --outdir results/supplementary_tables_<stamp>

  # one or two tables only, e.g. after rerunning one analysis
  python make_supplementary_tables.py ... --only S2,codon_usage

Dependencies: numpy, pandas, scipy, openpyxl; region_lib and rate_tables from
lib/.


Changes 2026-11-05.3
--------------------
  * S7 caption: the primary rows are the values in fig. 4, the branch figure
    since the renumbering of 2026-10-27 (it still said fig. 5). The S7 check
    and comment that name the branch panels now say fig. 4A and 4B-C.

Changes 2026-11-05.2
--------------------
  * The first Isaac run (job 6550984) gave eliminated alpha 0.228, 0.183 and
    0.225 at the three floors, so the Results no longer say "a fifth to a
    quarter" but give the range, 0.18 to 0.23; the check now compares the
    rounded range with that text.

Changes 2026-11-05
------------------
  * S11: a row for positive selection, built from pi_by_region.py's
    mk_summary.tsv (new manifest key mk_summary, primary callset): DoS with
    its 95% interval and alpha in each class at minor-allele-frequency floors
    of 0.10, 0.15 and 0.20, placed after the site-counting row, with a verdict
    computed from the values. The Results cite these floors ("at every
    threshold") and until now no table showed the 0.15 and 0.20 runs.
  * Checks: every floor present for both classes, on the same genes; the
    0.10 floor reproduces the values the Results and fig. 6D quote (DoS
    +0.102 / +0.064, alpha 0.35 / 0.23, 6,882 / 264 genes); eliminated DoS at
    or below retained DoS at every floor; eliminated alpha positive at every
    floor, with the range the Results give (2026-11-05.2: 0.18 to 0.23).
  * S11 title and caption name dN/dS as well as piN/piS, since positive
    selection is an alternative explanation for the dN/dS contrast.
  * New check: the curated shown_in column cites main figures 1 to 6 only
    (it still cited the pre-2026-10-27 Fig. 7, now fixed in the file).

Changes 2026-11-01
------------------
  * S9C: the insulation result at each of five Hi-C stages (testis, ovary,
    zygote, 48 h and 4-cell (60 h) embryos), each run without an eigenvector
    and with the expression of its own stage (early embryo for 48 h and 60 h):
    region medians with the Mann-Whitney test, the elimination coefficient of
    both models, and the insulation coefficient of the dN/dS model. Testis is
    the no-eigenvector run of A and B. New manifest keys: stage_pergene_* and
    stage_summary_* for every stage, stage_models_* for the four new runs.
  * S9B now has a Hi-C stage column and every coefficient of the four new runs.
  * Checks: each stage run used the insulation track and expression column of
    its stage and no eigenvector; its recomputed medians, P and counts match
    what it printed and its models; every stage fits the same genes with the
    same dN/dS.
  * S10 caption: the population sample (BioProject PRJNA1195031) is Landeryou
    et al. 2025, worms from Ethiopia, as the Methods state; the caption cited
    another study and country (fixed until now only in the clean repository).

Changes 2026-10-12.3
--------------------
  * Tables renumbered by first citation in the text (NUMBER), sheets and TSVs in
    that order; supplementary figure numbers in captions come from FIGURE;
    check names take their table numbers from NUMBER.
  * S2 caption (rate heterogeneity): each Tajima contrast is shown with both
    outgroups; the caption now says that the outgroup on the longer path pulls
    the ratio towards 1, gives both values, says which one the text quotes, and
    limits the stability-across-sites statement to the contrast it holds for.
    New check: the nearer outgroup gives the larger ratio in every site subset.
  * Summary note lines print confidence bounds below 0.1 to two significant
    digits (S3: [0.015, 3.69], not [0.02, 3.69]).
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import platform
import re
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from scipy import __version__ as scipy_version

from rate_tables import read_rate_table
from region_lib import fingerprint, git_commit, setup_logging

try:
    import openpyxl
    from openpyxl.styles import Alignment, Border, Font, Side
    from openpyxl.utils import get_column_letter
except ImportError:                                  # fail loudly, never fall back
    sys.exit('ERROR: openpyxl is required (it is in the kaks environment).')

_HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT_VERSION = '2026-11-05.3'

# ---------------------------------------------------------------------------
# Numbering and labels
# ---------------------------------------------------------------------------

#: Table key -> supplementary table number, by first citation in the text
#: (2026-10-12). Renumber here, and only here.
NUMBER = {
    'rates_pairwise': 'S1',
    'rate_heterogeneity': 'S2',
    'omega_above_one': 'S3',
    'codon_usage': 'S4',
    'microsynteny': 'S5',
    'x_autosome': 'S6',
    'branch_models': 'S7',
    'toxocara_models': 'S8',
    'insulation_models': 'S9',
    'popgen_cohort': 'S10',
    'popgen_alternatives': 'S11',
}
#: Sheets, TSVs and the summary follow the table numbers.
TABLE_ORDER = sorted(NUMBER, key=lambda k: int(NUMBER[k][1:]))
#: Supplementary figures that captions cite, by first citation in the text.
FIGURE = {
    'branch_tox_region': 'S8',      # branch dN/dS with T. canis labels (was S14)
}

#: Manifest keys each table needs, and keys it uses when present.
REQUIRES = {
    'rates_pairwise': ['annotated', 'gene_structure'],
    'branch_models': ['branch_stats', 'branch_quartets', 'branch_decomposition'],
    'omega_above_one': ['omega_genes', 'omega_stats', 'annotated'],
    'codon_usage': [f'codon_{k}_{sp}' for sp in ('ascaris', 'parascaris', 'toxocara')
                    for k in ('summary', 'report', 'genes')],
    'insulation_models': ['models_noEV', 'models_EVtestis', 'models_EV4cell',
                          'stage_pergene_testis', 'stage_summary_testis'] +
                         [f'stage_{k}_{s}' for s in ('ovary', 'zygote', '48h', '60h')
                          for k in ('models', 'pergene', 'summary')],
    'rate_heterogeneity': ['rate_table', 'rate_caption', 'clock_sites12',
                           'clock_sites3', 'clock_sitesall'],
    'toxocara_models': ['tox_premise', 'tox_strata', 'tox_regression'],
    'microsynteny': ['synteny_summary', 'synteny_tests', 'synteny_sweep'],
    'popgen_cohort': ['cohort_metadata'],
    'popgen_alternatives': ['popgen_alternatives', 'mk_summary'],
    'x_autosome': ['x_autosome', 'x_dosage'],
}
OPTIONAL = {
    'toxocara_models': ['branch_decomposition'],   # cross-check of S7B against S2A
    'popgen_cohort': ['cohort_audit'],             # cross-check of S9 values
}

#: Species as the tables name them: genus-only Ascaris, abbreviated binomials.
SPECIES_ORDER = ['Ascaris', 'Parascaris', 'Baylisascaris', 'Toxocara']
SPECIES = {'Ascaris': 'Ascaris', 'Parascaris': 'P. univalens',
           'Baylisascaris': 'B. schroederi', 'Toxocara': 'T. canis'}
#: ...and as a column-name prefix.
SPECIES_COL = {'Ascaris': 'Ascaris', 'Parascaris': 'P_univalens',
               'Baylisascaris': 'B_schroederi', 'Toxocara': 'T_canis'}
LABELLING = {'as_region': 'Ascaris', 'tox_region': 'T. canis'}
LABELLING_ORDER = ['as_region', 'tox_region']
METRIC = {'omega': 'dN/dS', 'dN': 'dN', 'dS': 'dS'}
METRIC_ORDER = ['omega', 'dN', 'dS']

#: Written in place of a P value that underflowed to exactly 0.
P_UNDERFLOW = '<1e-300'

#: Relative tolerance for comparing with values the analysis scripts printed to six
#: significant digits (%.6g): half a unit in the sixth digit is up to 5e-6.
PRINTED_RTOL = 1e-5


# ---------------------------------------------------------------------------
# Table model
# ---------------------------------------------------------------------------

@dataclass
class Col:
    """One published column.

    name  published column name
    src   column of the part's frame it is taken from
    desc  definition, for the Columns sheet
    fmt   display format: int, dist (bp, one decimal), f1-f4, v4 (up to four
          decimals), pct, p, g, text, bool. The TSVs keep ten significant digits
          for every float except dist; the formats set what the workbook shows.
    optional  drop the column (logged) when every value is missing
    values    append the values that occur to the definition
    """
    name: str
    src: str
    desc: str
    fmt: str = 'g'
    optional: bool = False
    values: bool = False


@dataclass
class Part:
    suffix: str                 # '' for a one-part table, else 'A', 'B', ...
    title: str                  # heading of the part
    frame: pd.DataFrame         # internal column names
    cols: list
    drop: dict = field(default_factory=dict)   # source column -> why it is left out
    freeze_cols: int = 0
    wrap: bool = False          # wrap long text in data cells
    sources: tuple = ()


@dataclass
class Table:
    key: str
    title: str
    caption: str
    parts: list
    notes: list = field(default_factory=list)  # numbers for the text, for the summary


class Checks:
    """Every consistency check, in order; failures are warnings, not errors."""

    def __init__(self):
        self.items = []

    def add(self, table, name, ok, detail=''):
        ok = bool(ok)
        self.items.append({'table': NUMBER.get(table, table), 'check': name,
                           'ok': ok, 'detail': detail})
        msg = f"[{NUMBER.get(table, table)}] {'OK' if ok else 'MISMATCH'}: {name}"
        if detail:
            msg += f' ({detail})'
        (logging.info if ok else logging.warning)(msg)
        return ok

    @property
    def n_failed(self):
        return sum(not c['ok'] for c in self.items)


class Context:
    def __init__(self, sources, args):
        self.src = sources
        self.args = args
        self.checks = Checks()
        self._cache = {}

    def tsv(self, key, rate=False, **kw):
        """Read a manifest source once; `rate` reads it through rate_tables."""
        path = self.src[key]
        ck = (path, rate)
        if ck not in self._cache:
            if rate:
                self._cache[ck] = read_rate_table(path, require=('dN', 'dS', 'dN_dS'))
            else:
                self._cache[ck] = pd.read_csv(path, sep='\t', **kw)
            logging.info('  read %s: %s rows (%s)', key, f'{len(self._cache[ck]):,}',
                         os.path.basename(path))
        return self._cache[ck].copy()


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def need(df, cols, label):
    """Stop with the columns that are present when expected ones are missing."""
    missing = [c for c in cols if c not in df.columns]
    if missing:
        sys.exit(f'ERROR: {label}: missing column(s) {missing}. '
                 f'Columns present: {list(df.columns)}')


def to_bool(s):
    """True/False from bools or their string forms; anything else is NA."""
    m = {True: True, False: False, 'True': True, 'False': False, 'TRUE': True,
         'FALSE': False, 'true': True, 'false': False}
    return s.map(lambda v: m.get(v, pd.NA) if not (isinstance(v, float) and np.isnan(v))
                 else pd.NA).astype('boolean')


def close(a, b, rtol=1e-9, atol=0.0):
    """Elementwise closeness with NaN == NaN. Returns a numpy bool array."""
    a = np.asarray(pd.to_numeric(pd.Series(a), errors='coerce'), dtype=float)
    b = np.asarray(pd.to_numeric(pd.Series(b), errors='coerce'), dtype=float)
    both_nan = np.isnan(a) & np.isnan(b)
    with np.errstate(invalid='ignore'):
        ok = np.abs(a - b) <= atol + rtol * np.maximum(np.abs(a), np.abs(b))
    return ok | both_nan


def max_rel(a, b):
    """Largest relative difference, for check details."""
    a = np.asarray(pd.to_numeric(pd.Series(a), errors='coerce'), dtype=float)
    b = np.asarray(pd.to_numeric(pd.Series(b), errors='coerce'), dtype=float)
    ok = ~(np.isnan(a) | np.isnan(b))
    if not ok.any():
        return 0.0
    d = np.abs(a[ok] - b[ok]) / np.maximum(np.maximum(np.abs(a[ok]), np.abs(b[ok])), 1e-300)
    return float(d.max())


def parse_p(s):
    """'5.06 x 10^-10', '0.005' or '< 1e-300' as printed in a summary -> float."""
    s = str(s).strip()
    m = re.match(r'^([0-9.]+)\s*[x×]\s*10\^\s*([+-]?\d+)$', s)
    if m:
        return float(m.group(1)) * 10.0 ** int(m.group(2))
    if s.startswith('<'):
        return 0.0
    return float(s)


def p_words(p):
    """P for a caption or a summary line, in plain text."""
    if p is None or not np.isfinite(p):
        return 'NA'
    if p == 0:
        return '< 1e-300'
    return f'{p:.2g}' if p < 0.001 else f'{p:.3f}'


def fmt_num(v, dp):
    """A number for a summary note: dp decimals, or two significant digits below
    0.1, so that a confidence bound of 0.0155 prints as 0.015 rather than 0.02."""
    v = float(v)
    if not np.isfinite(v):
        return 'NA'
    if v != 0 and abs(v) < 0.1:
        return f'{v:.2g}'
    return f'{v:.{dp}f}'


def mwu(ret, elim):
    """Two-sided Mann-Whitney U, exactly as the analysis scripts call it:
    retained first, so U is the retained sample's statistic."""
    res = mannwhitneyu(ret, elim, alternative='two-sided')
    return float(res.statistic), float(res.pvalue)


def describe(vals):
    """n, median, Q1, Q3 as pandas computes them (the analysis scripts' way)."""
    v = pd.Series(vals, dtype=float).dropna()
    return {'n': int(len(v)), 'median': v.median(), 'Q1': v.quantile(0.25),
            'Q3': v.quantile(0.75)}


def _conditional_or_fallback(a, b, c, d, level):
    """Conditional MLE odds ratio and exact CI for scipy < 1.10 (no
    scipy.stats.contingency.odds_ratio): `a` follows Fisher's noncentral
    hypergeometric distribution given the margins; the limits are where the
    observed count sits in either (1 - level)/2 tail, as in R's fisher.test.
    The same code as microsynteny.py, so every odds ratio in the paper is one
    estimator."""
    from scipy.optimize import brentq
    from scipy.stats import nchypergeom_fisher
    M, n_col, N = a + b + c + d, a + c, a + b
    lo_sup, hi_sup = max(0, N - (M - n_col)), min(N, n_col)
    alpha = (1 - level) / 2

    def dist(log_or):
        return nchypergeom_fisher(M, n_col, N, np.exp(log_or))

    bound = 20.0
    mle = np.nan
    if a == lo_sup:
        mle, lo = 0.0, 0.0
    else:
        lo = np.exp(brentq(lambda x: dist(x).sf(a - 1) - alpha, -bound, bound))
    if a == hi_sup:
        mle, hi = np.inf, np.inf
    else:
        hi = np.exp(brentq(lambda x: dist(x).cdf(a) - alpha, -bound, bound))
    if lo_sup < a < hi_sup:
        mle = np.exp(brentq(lambda x: dist(x).mean() - a, -bound, bound))
    return mle, lo, hi


def conditional_or(a, b, c, d, level=0.95):
    """Conditional MLE odds ratio and exact CI for [[a, b], [c, d]]."""
    try:
        from scipy.stats.contingency import odds_ratio
    except ImportError:
        return tuple(float(x) for x in _conditional_or_fallback(a, b, c, d, level))
    res = odds_ratio(np.array([[a, b], [c, d]]), kind='conditional')
    ci = res.confidence_interval(confidence_level=level)
    return float(res.statistic), float(ci.low), float(ci.high)


def species_pair(s):
    """'Ascaris x Parascaris' -> 'Ascaris × P. univalens'."""
    parts = [p.strip() for p in re.split(r'\s+x\s+', str(s))]
    return ' × '.join(SPECIES.get(p, p) for p in parts)


# ---------------------------------------------------------------------------
# S1  Pairwise rates
# ---------------------------------------------------------------------------

def build_rates_pairwise(ctx):
    T = 'rates_pairwise'
    ck = ctx.checks
    tab = ctx.tsv('annotated', rate=True)
    need(tab, ['gene1', 'gene2', 'align_len', 'status', 'chrom1', 'start1', 'end1',
               'midpoint1', 'region1', 'elim_type1', 'cbr_proximal1', 'min_cbr_dist1',
               'signed_cbr_dist1', 'chrom2', 'start2', 'end2', 'midpoint2', 'region2',
               'elim_type2', 'cbr_proximal2', 'min_cbr_dist2', 'signed_cbr_dist2',
               'consensus_region', 'cbr_proximal_either', 'proximity1'],
         ctx.src['annotated'])
    n_dup = int(tab['gene1'].duplicated().sum())
    ck.add(T, 'one row per Ascaris gene', n_dup == 0, f'{n_dup} repeated')

    gs = ctx.tsv('gene_structure')
    need(gs, ['gene_id', 'region', 'chrom', 'start', 'end', 'strand', 'genomic_span',
              'cds_length', 'intron_fraction', 'intergenic_nearest'],
         ctx.src['gene_structure'])
    gs_dup = int(gs['gene_id'].duplicated().sum())
    ck.add(T, 'one gene-structure row per gene', gs_dup == 0, f'{gs_dup} repeated')
    gs = gs.drop_duplicates('gene_id').add_prefix('gs_')
    m = tab.merge(gs, left_on='gene1', right_on='gs_gene_id', how='left')
    have = m['gs_gene_id'].notna()
    ck.add(T, 'gene structure for every pair', have.all(),
           f'{int(have.sum()):,} of {len(m):,}')
    reg_bad = int((m.loc[have, 'region1'].astype(str) != m.loc[have, 'gs_region'].astype(str)).sum())
    ck.add(T, 'region1 equals the gene-structure region', reg_bad == 0, f'{reg_bad} differ')
    chr_bad = int((m.loc[have, 'chrom1'].astype(str) != m.loc[have, 'gs_chrom'].astype(str)).sum())
    pos_bad = int((~close(m.loc[have, 'start1'], m.loc[have, 'gs_start'])
                   | ~close(m.loc[have, 'end1'], m.loc[have, 'gs_end'])).sum())
    ck.add(T, 'coordinates equal the gene-structure coordinates', chr_bad + pos_bad == 0,
           f'{chr_bad} chromosomes, {pos_bad} positions differ')
    span_bad = int((~close(m.loc[have, 'gs_genomic_span'],
                           m.loc[have, 'end1'] - m.loc[have, 'start1'] + 1)).sum())
    ck.add(T, 'genomic span = end - start + 1', span_bad == 0, f'{span_bad} differ')

    drop = {'midpoint1': 'derivable from start and end',
            'midpoint2': 'derivable from start and end',
            'proximity1': 'same information as CBR_proximal_ascaris',
            'gs_gene_id': 'join key', 'gs_region': 'equals region_ascaris (checked)',
            'gs_chrom': 'equals chromosome_ascaris (checked)',
            'gs_start': 'equals start_ascaris (checked)',
            'gs_end': 'equals end_ascaris (checked)',
            'gs_intergenic_up': 'nearest distance kept; the two sides are not used',
            'gs_intergenic_down': 'nearest distance kept; the two sides are not used',
            'gs_locus_models': 'pipeline bookkeeping'}
    status_vals = sorted(m['status'].dropna().astype(str).unique())
    cols = [
        Col('gene_ascaris', 'gene1', 'Ascaris gene model: miniprot protein-to-genome '
            'alignment on the germline assembly AG.v50, longest isoform per locus', 'text'),
        Col('gene_parascaris', 'gene2', 'P. univalens gene model: miniprot alignment on '
            'the germline assembly pu_v1.13', 'text'),
        Col('dN', 'dN', 'Nonsynonymous substitutions per nonsynonymous site (yn00)', 'f4'),
        Col('dS', 'dS', 'Synonymous substitutions per synonymous site (yn00)', 'f4'),
        Col('dN_dS', 'dN_dS', 'dN/dS', 'f4'),
        Col('alignment_length_codons', 'align_len', 'Length of the pairwise codon '
            'alignment, codons', 'int'),
    ]
    if len(status_vals) <= 1:
        drop['status'] = f"constant ({status_vals[0] if status_vals else 'empty'})"
    else:
        cols.append(Col('filter_status', 'status', 'yn00 filter outcome', 'text', values=True))
    for i, sp, lab in ((1, 'ascaris', 'Ascaris'), (2, 'parascaris', 'P. univalens')):
        cols += [
            Col(f'chromosome_{sp}', f'chrom{i}', f'{lab} germline sequence', 'text'),
            Col(f'start_{sp}', f'start{i}', f'Start of the {lab} gene model (1-based)', 'int'),
            Col(f'end_{sp}', f'end{i}', f'End of the {lab} gene model (1-based, inclusive)', 'int'),
            Col(f'region_{sp}', f'region{i}', f'{lab} region containing the gene midpoint, '
                f'from the {lab} elimination annotation', 'text', values=True),
            Col(f'eliminated_tract_type_{sp}', f'elim_type{i}', f'For eliminated genes, '
                f'the {lab} eliminated tract containing the gene: end (terminal) or '
                'internal', 'text', optional=True, values=True),
            Col(f'distance_to_nearest_CBR_{sp}_bp', f'min_cbr_dist{i}', 'Distance from '
                'the gene midpoint to the nearest chromosomal breakage region (CBR) on the '
                'same sequence, bp. CBRs are the boundaries of merged eliminated tracts, '
                'boundaries at sequence ends excluded, placed at the half-base junction '
                '(so distances end in .5). Empty on sequences without a CBR', 'dist'),
            Col(f'signed_distance_to_CBR_{sp}_bp', f'signed_cbr_dist{i}', 'The same '
                'distance, positive on the eliminated side of the CBR', 'dist'),
            Col(f'CBR_proximal_{sp}', f'cbr_proximal{i}', 'Gene midpoint within 100 kb '
                'of a CBR', 'bool'),
        ]
    cols += [
        Col('region_consensus', 'consensus_region', 'The region where the two '
            'annotations agree; discordant where they disagree; unannotated where '
            'either lacks a retained or eliminated call', 'text', values=True),
        Col('CBR_proximal_either', 'cbr_proximal_either', 'CBR-proximal in either '
            'species', 'bool'),
        Col('strand_ascaris', 'gs_strand', 'Strand of the Ascaris gene model', 'text'),
        Col('genomic_span_bp', 'gs_genomic_span', 'Span of the Ascaris gene model, '
            'first to last base, bp', 'int'),
        Col('CDS_length_bp', 'gs_cds_length', 'Length of the Ascaris coding sequence, bp',
            'int'),
        Col('intron_fraction', 'gs_intron_fraction', '1 − CDS length / genomic span', 'f4'),
        Col('intergenic_distance_nearest_bp', 'gs_intergenic_nearest', 'Distance from '
            'the Ascaris gene model to the nearest neighbouring gene model, measured '
            'against every model in the annotation, bp (negative where models overlap)',
            'int'),
    ]
    n = len(m)
    n_ret = int((m['region1'] == 'retained').sum())
    n_elim = int((m['region1'] == 'eliminated').sum())
    title = (f'Pairwise substitution rates for {n:,} Ascaris × P. univalens one-to-one '
             'orthologue pairs, with region annotations and gene structure.')
    caption = (
        'dN, dS and dN/dS were estimated with yn00 (Yang and Nielsen 2000) from the codon '
        'alignment of each pair, after excluding estimates with dS > 5, dS = 0 or dN/dS > 50 '
        'and alignments shorter than 100 codons (Materials and Methods). Region fields come '
        "from each species' own elimination annotation, in germline coordinates: a gene is "
        'assigned to the region that contains its midpoint, and distances are measured from '
        'the midpoint to the nearest chromosomal breakage region (CBR). All principal '
        f'analyses use the Ascaris assignment (region_ascaris: {n_ret:,} retained, '
        f'{n_elim:,} eliminated). Gene-structure fields describe the Ascaris gene model.')
    part = Part('', '', m, cols, drop=drop, freeze_cols=1,
                sources=('annotated', 'gene_structure'))
    notes = [f'{n:,} pairs: {n_ret:,} retained, {n_elim:,} eliminated (region_ascaris)',
             'consensus: ' + ', '.join(f'{k} {v:,}' for k, v in
                                       m['consensus_region'].value_counts().items())]
    return Table(T, title, caption, [part], notes)


# ---------------------------------------------------------------------------
# S2  Four-species branch models
# ---------------------------------------------------------------------------

def usable(q, sp, metric, cap, apply_cap):
    """Genes usable for one metric on one branch -- branch_lib.usable_mask.

    Capped (dN/dS, dN): dN/dS, dN and dS finite, dN/dS <= cap, dS > 0.
    Uncapped (dS): the metric finite. Negative values (codeml's -1 for an
    estimate it could not obtain) are always excluded.
    """
    om = pd.to_numeric(q[f'{sp}_omega'], errors='coerce')
    dn = pd.to_numeric(q[f'{sp}_dN'], errors='coerce')
    ds = pd.to_numeric(q[f'{sp}_dS'], errors='coerce')
    vals = {'omega': om, 'dN': dn, 'dS': ds}[metric]
    keep = ((om.notna() & dn.notna() & ds.notna() & (om <= cap) & (ds > 0))
            if apply_cap else vals.notna())
    return (keep & (vals >= 0)).fillna(False).astype(bool)


def build_branch_models(ctx):
    T = 'branch_models'
    ck = ctx.checks
    cap = ctx.args.omega_cap
    bs = ctx.tsv('branch_stats')
    dec = ctx.tsv('branch_decomposition')
    q4 = ctx.tsv('branch_quartets')
    need(bs, ['branch', 'grouping', 'metric', 'ret_n', 'elim_n', 'ret_median',
              'elim_median', 'fold_change', 'p_value', 'q_value', 'ret_ci_lo', 'ret_ci_hi',
              'elim_ci_lo', 'elim_ci_hi'], ctx.src['branch_stats'])
    need(dec, ['branch', 'grouping', 'metric', 'masked', 'ret_n', 'elim_n', 'ret_median',
               'elim_median', 'fold_change', 'p_value', 'q_value', 'q_star'],
         ctx.src['branch_decomposition'])
    rate_cols = [f'{sp}_{m}' for sp in SPECIES_ORDER + ['internal'] for m in METRIC_ORDER]
    need(q4, ['as_gene', 'pu_gene', 'bsc_gene', 'tox_gene', 'align_len', 'as_region',
              'tox_region'] + rate_cols, ctx.src['branch_quartets'])
    dec['masked'] = to_bool(dec['masked'])
    if dec['masked'].isna().any():
        sys.exit(f"ERROR: {ctx.src['branch_decomposition']}: 'masked' is not True/False")
    dec['masked'] = dec['masked'].astype(bool)

    # --- S2B: the quartets -----------------------------------------------------
    if 'status' in q4.columns:
        st = sorted(q4['status'].dropna().astype(str).unique())
        ck.add(T, 'every quartet has status OK', st == ['OK'], ', '.join(st))
        q_ok = q4[q4['status'].astype(str) == 'OK'].copy()
    else:
        q_ok = q4.copy()
    if 'internal_subtends' in q4.columns:
        splits = sorted(set(q4['internal_subtends'].dropna().astype(str)) - {''})
        ck.add(T, 'unrooted tree: the internal branch splits (Ascaris, P. univalens) '
               'from (B. schroederi, T. canis) in every quartet',
               splits == ['Ascaris|Parascaris'], ', '.join(splits) or 'none recorded')
    for sp in SPECIES_ORDER:
        q4[f'{sp}_masked'] = ~usable(q4, sp, 'omega', cap, True)

    # --- recompute every decomposition row from the quartets ---------------------
    worst = {'n': 0, 'median': 0.0, 'p': 0.0}
    n_rows = 0
    for _, r in dec.iterrows():
        g, sp, metric, masked = r['grouping'], r['branch'], r['metric'], bool(r['masked'])
        if g not in q_ok.columns or sp not in SPECIES_ORDER or metric not in METRIC_ORDER:
            ck.add(T, f'decomposition row understood: {g} {sp} {metric}', False)
            continue
        keep = usable(q_ok, sp, metric, cap, masked)
        sub = q_ok[keep]
        v = pd.to_numeric(sub[f'{sp}_{metric}'], errors='coerce')
        ret = v[sub[g] == 'retained'].dropna()
        elim = v[sub[g] == 'eliminated'].dropna()
        _, p = mwu(ret, elim)
        worst['n'] = max(worst['n'], abs(len(ret) - int(r['ret_n'])),
                         abs(len(elim) - int(r['elim_n'])))
        worst['median'] = max(worst['median'], max_rel(ret.median(), r['ret_median']),
                              max_rel(elim.median(), r['elim_median']))
        worst['p'] = max(worst['p'], max_rel(p, r['p_value']))
        n_rows += 1
    S = NUMBER[T]
    ck.add(T, f'{S}A n and medians reproduced from {S}B with the masking rule (cap {cap:g})',
           worst['n'] == 0 and worst['median'] < 1e-9,
           f"{n_rows} rows; largest n difference {worst['n']}, "
           f"median rel. diff {worst['median']:.1e}")
    ck.add(T, f'{S}A P values reproduced from {S}B', worst['p'] < PRINTED_RTOL,
           f"largest rel. diff {worst['p']:.1e}")

    # --- branch_stats (fig. 4A) against the decomposition's dN/dS rows -----------
    om = dec[(dec['metric'] == 'omega') & dec['masked']]
    j = bs.merge(om, on=['branch', 'grouping'], suffixes=('_bs', '_dec'), how='outer',
                 indicator=True)
    both = j['_merge'] == 'both'
    ck.add(T, 'branch_stats and decomposition cover the same branches and labellings',
           both.all(), f'{int((~both).sum())} unmatched rows')
    agree = np.ones(int(both.sum()), dtype=bool)
    for c in ('ret_n', 'elim_n', 'ret_median', 'elim_median', 'fold_change', 'p_value'):
        agree &= close(j.loc[both, f'{c}_bs'], j.loc[both, f'{c}_dec'], rtol=PRINTED_RTOL)
    ck.add(T, 'dN/dS rows identical in branch_stats (fig. 4A) and the decomposition '
           '(fig. 4B, 4C)', agree.all(), f'{int((~agree).sum())} differ')

    # --- S2A -------------------------------------------------------------------------
    a = dec.merge(bs[['branch', 'grouping', 'ret_ci_lo', 'ret_ci_hi', 'elim_ci_lo',
                      'elim_ci_hi']].assign(metric='omega', masked=True),
                  on=['branch', 'grouping', 'metric', 'masked'], how='left')
    a['primary'] = (((a['metric'] == 'omega') & a['masked'])
                    | ((a['metric'] == 'dN') & a['masked'])
                    | ((a['metric'] == 'dS') & ~a['masked'])).astype(bool)
    a['labelling'] = a['grouping'].map(LABELLING)
    a['branch_label'] = a['branch'].map(SPECIES)
    a['metric_label'] = a['metric'].map(METRIC)
    a['mask_label'] = np.where(a['masked'], 'applied', 'none')
    for c in ('labelling', 'branch_label', 'metric_label'):
        if a[c].isna().any():
            sys.exit(f'ERROR: unexpected grouping, branch or metric in '
                     f"{ctx.src['branch_decomposition']}")
    a = a.assign(_p=~a['primary'],
                 _g=a['grouping'].map(LABELLING_ORDER.index),
                 _m=a['metric'].map(METRIC_ORDER.index),
                 _b=a['branch'].map(SPECIES_ORDER.index)
                 ).sort_values(['_p', '_g', '_m', '_b']).reset_index(drop=True)
    cols_a = [
        Col('labelling', 'labelling', 'Elimination annotation used to label genes',
            'text', values=True),
        Col('branch', 'branch_label', 'Terminal branch', 'text'),
        Col('metric', 'metric_label', 'Branch estimate compared', 'text', values=True),
        Col('masking', 'mask_label', 'applied: masked branches excluded; none: every '
            'branch with an estimate', 'text', values=True),
        Col('primary', 'primary', 'The comparison reported in the text and figures (dN/dS '
            'and dN with masking, dS without); the other rows apply the other rule', 'bool'),
        Col('n_retained', 'ret_n', 'Retained genes compared', 'int'),
        Col('n_eliminated', 'elim_n', 'Eliminated genes compared', 'int'),
        Col('median_retained', 'ret_median', 'Median over retained genes', 'f4'),
        Col('median_retained_CI_low', 'ret_ci_lo', 'dN/dS rows: lower 95% bootstrap '
            'bound of the retained median', 'f4'),
        Col('median_retained_CI_high', 'ret_ci_hi', 'dN/dS rows: upper bound', 'f4'),
        Col('median_eliminated', 'elim_median', 'Median over eliminated genes', 'f4'),
        Col('median_eliminated_CI_low', 'elim_ci_lo', 'dN/dS rows: lower 95% bootstrap '
            'bound of the eliminated median', 'f4'),
        Col('median_eliminated_CI_high', 'elim_ci_hi', 'dN/dS rows: upper bound', 'f4'),
        Col('fold_change', 'fold_change', 'Eliminated median / retained median', 'f3'),
        Col('P', 'p_value', 'Two-sided Mann–Whitney U test', 'p'),
        Col('q', 'q_star', 'Benjamini–Hochberg q within each metric and masking rule, '
            'over the four branches and two labellings (eight tests)', 'p'),
    ]
    drop_a = {'grouping': 'shown as labelling', 'branch': 'shown as branch',
              'metric': 'shown as metric', 'masked': 'shown as masking',
              'q_value': 'BH within one labelling (four tests); the text and figures use '
                         'q over both labellings (eight tests)',
              '_p': 'sort key', '_g': 'sort key', '_m': 'sort key', '_b': 'sort key'}

    # --- S2B -------------------------------------------------------------------------
    cols_b = [
        Col('gene_ascaris', 'as_gene', 'Ascaris gene model', 'text'),
        Col('gene_parascaris', 'pu_gene', 'P. univalens gene model', 'text'),
        Col('gene_baylisascaris', 'bsc_gene', 'B. schroederi gene (Augustus prediction)',
            'text'),
        Col('gene_toxocara', 'tox_gene', 'T. canis open reading frame (TransDecoder, on '
            'StringTie transcripts)', 'text'),
        Col('region_ascaris_annotation', 'as_region', 'Region of the Ascaris gene in the '
            'Ascaris annotation', 'text', values=True),
        Col('region_toxocara_annotation', 'tox_region', 'Region of the T. canis gene in '
            'the T. canis annotation (Wang et al. 2017)', 'text', values=True),
        Col('alignment_length_codons', 'align_len', 'Length of the four-species codon '
            'alignment, codons', 'int'),
    ]
    for sp in SPECIES_ORDER:
        p_, lab = SPECIES_COL[sp], SPECIES[sp]
        cols_b += [
            Col(f'{p_}_dN_dS', f'{sp}_omega', f'dN/dS (ω) on the {lab} branch', 'f4'),
            Col(f'{p_}_dN', f'{sp}_dN', f'dN on the {lab} branch', 'f4'),
            Col(f'{p_}_dS', f'{sp}_dS', f'dS on the {lab} branch', 'f4'),
            Col(f'{p_}_masked', f'{sp}_masked', f'True when the {lab} branch is masked '
                f'(dN/dS > {cap:g}, no substitution, or no estimate) and excluded from the '
                'dN/dS and dN comparisons', 'bool'),
        ]
    cols_b += [
        Col('internal_dN_dS', 'internal_omega', 'dN/dS on the internal branch (not '
            'interpreted)', 'f4'),
        Col('internal_dN', 'internal_dN', 'dN on the internal branch', 'f4'),
        Col('internal_dS', 'internal_dS', 'dS on the internal branch', 'f4'),
    ]
    drop_b = {'status': 'constant (OK; checked)',
              'internal_subtends': 'constant (Ascaris|Parascaris; checked)'}

    n_q = len(q4)
    title = ('Branch-specific dN, dS and dN/dS from codeml free-ratio models of '
             f'{n_q:,} four-species orthologue quartets.')
    caption = (
        'Quartets of one-to-one orthologues of Ascaris, P. univalens, B. schroederi and '
        'T. canis were analysed with the codeml free-ratio model on the unrooted tree '
        '((Ascaris, P. univalens), B. schroederi, T. canis). '
        '(A) Each row compares eliminated with retained genes on one terminal branch, with '
        'genes labelled by the Ascaris or by the T. canis elimination annotation. A branch '
        f'is masked when its dN/dS exceeds {cap:g} or when it carries no substitution '
        '(dN = dS = 0); masking removes that branch only. dN/dS and dN are compared with '
        'masked branches excluded and dS on every branch, because a cap on dN/dS is also a '
        'filter on dS (Materials and Methods). These primary rows are the values in fig. 4 '
        f"and supplementary fig. {FIGURE['branch_tox_region']}; the remaining two rows per "
        'branch apply the other '
        'rule to dN and to dS. Fold change, eliminated median divided by retained median; '
        'P, two-sided Mann–Whitney U test; q, Benjamini–Hochberg correction within each '
        'metric and masking rule over the four branches and two labellings (eight tests). '
        'Intervals on dN/dS medians are percentile bootstrap 95% intervals (2,000 '
        'resamples). (B) Per-quartet estimates. The internal branch separates '
        '(Ascaris, P. univalens) from (B. schroederi, T. canis); it is not an ancestral '
        'branch of either pair and is not interpreted.')
    prim = a[a['primary']]
    notes = []
    for _, r in prim.iterrows():
        notes.append(f"{r['labelling']:<8} {r['metric_label']:<5} {r['branch_label']:<13} "
                     f"{int(r['ret_n']):>5}/{int(r['elim_n']):<4} fold {r['fold_change']:.2f}  "
                     f"P {p_words(r['p_value'])}  q {p_words(r['q_star'])}")
    parts = [Part('A', '(A) Eliminated versus retained genes on each terminal branch.',
                  a, cols_a, drop=drop_a,
                  sources=('branch_decomposition', 'branch_stats', 'branch_quartets')),
             Part('B', '(B) Per-quartet branch estimates.', q4, cols_b, drop=drop_b,
                  freeze_cols=1, sources=('branch_quartets',))]
    return Table(T, title, caption, parts, notes)


# ---------------------------------------------------------------------------
# S3  dN/dS above 1
# ---------------------------------------------------------------------------

def build_omega_above_one(ctx):
    T = 'omega_above_one'
    ck = ctx.checks
    og = ctx.tsv('omega_genes', rate=True)
    st = ctx.tsv('omega_stats')
    tab = ctx.tsv('annotated', rate=True)
    need(og, ['gene1', 'gene2', 'region1', 'consensus_region', 'elim_type1', 'chrom1',
              'start1', 'end1', 'signed_cbr_dist1', 'cbr_proximal1', 'align_len'],
         ctx.src['omega_genes'])
    need(st, ['analysis', 'min_ds', 'threshold', 'n_eliminated', 'k_eliminated',
              'n_retained', 'k_retained', 'prop_eliminated', 'prop_retained',
              'odds_ratio', 'or_lo95', 'or_hi95', 'p'], ctx.src['omega_stats'])
    thr = sorted(pd.to_numeric(st['threshold'], errors='coerce').dropna().unique())
    if len(thr) != 1:
        sys.exit(f"ERROR: {ctx.src['omega_stats']}: expected one dN/dS threshold, got {thr}")
    thr = float(thr[0])

    # --- the gene list is exactly the S1 pairs above the threshold ---------------
    above = tab[pd.to_numeric(tab['dN_dS'], errors='coerce') > thr]
    same = set(og['gene1']) == set(above['gene1']) and len(og) == len(above)
    S, S1 = NUMBER[T], NUMBER['rates_pairwise']
    ck.add(T, f'{S}A lists exactly the {S1} pairs with dN/dS > {thr:g}', same,
           f'{len(og)} listed, {len(above)} in {S1}')
    jj = og.merge(above[['gene1', 'gene2', 'dN', 'dS', 'dN_dS', 'region1']], on='gene1',
                  suffixes=('', '_t1'))
    vals_ok = (len(jj) == len(og)
               and close(jj['dN'], jj['dN_t1'], rtol=PRINTED_RTOL).all()
               and close(jj['dS'], jj['dS_t1'], rtol=PRINTED_RTOL).all()
               and close(jj['dN_dS'], jj['dN_dS_t1'], rtol=PRINTED_RTOL).all()
               and (jj['gene2'] == jj['gene2_t1']).all()
               and (jj['region1'] == jj['region1_t1']).all())
    ck.add(T, f'{S}A values identical to {S1}', vals_ok)

    # --- S3B -----------------------------------------------------------------------
    st = st.copy()
    st['min_ds'] = pd.to_numeric(st['min_ds'], errors='coerce')
    prim = st[st['analysis'] == 'primary']
    sens = st[st['analysis'] != 'primary'].copy()
    floor0 = sens[sens['min_ds'] == 0]
    dup = (len(prim) == 1 and len(floor0) == 1 and all(
        close(prim[c].values, floor0[c].values).all()
        for c in ('n_eliminated', 'k_eliminated', 'n_retained', 'k_retained',
                  'odds_ratio', 'or_lo95', 'or_hi95', 'p')))
    ck.add(T, "the 'primary' row equals the no-floor sensitivity row", dup)
    b = (sens if dup else st).sort_values('min_ds', kind='mergesort').reset_index(drop=True)
    b['pct_eliminated'] = 100.0 * pd.to_numeric(b['prop_eliminated'])
    b['pct_retained'] = 100.0 * pd.to_numeric(b['prop_retained'])
    # Which odds ratio is in the file? Methods reports the conditional MLE with an
    # exact interval for every Fisher test (decided 2026-10-12 for this one too;
    # positive_selection.py wrote the Haldane-Anscombe estimate before that).
    ke, ne = b['k_eliminated'].astype(float), b['n_eliminated'].astype(float)
    kr, nr = b['k_retained'].astype(float), b['n_retained'].astype(float)
    hal = ((ke + .5) * (nr - kr + .5)) / ((ne - ke + .5) * (kr + .5))
    se = np.sqrt(1 / (ke + .5) + 1 / (ne - ke + .5) + 1 / (kr + .5) + 1 / (nr - kr + .5))
    is_hal = (close(hal, b['odds_ratio'], rtol=1e-4).all()
              and close(np.exp(np.log(hal) - 1.959964 * se), b['or_lo95'], rtol=1e-3).all()
              and close(np.exp(np.log(hal) + 1.959964 * se), b['or_hi95'], rtol=1e-3).all())
    cond = pd.DataFrame([conditional_or(int(k1), int(n1 - k1), int(k2), int(n2 - k2))
                         for k1, n1, k2, n2 in zip(ke, ne, kr, nr)],
                        columns=['or', 'lo', 'hi'])
    is_cond = (close(cond['or'], b['odds_ratio'], rtol=1e-4).all()
               and close(cond['lo'], b['or_lo95'], rtol=1e-3).all()
               and close(cond['hi'], b['or_hi95'], rtol=1e-3).all())
    or_method = 'conditional' if is_cond else 'haldane' if is_hal else 'unknown'
    ck.add(T, 'odds ratios are conditional MLE with exact intervals, as Methods states',
           is_cond, {'conditional': 'reproduced from the counts',
                     'haldane': 'the file holds Haldane–Anscombe estimates with Woolf '
                                'intervals: rerun positive_selection.py with the '
                                'conditional estimate',
                     'unknown': 'reproduced by neither estimator'}[or_method])
    n_e_all = int((tab['region1'] == 'eliminated').sum())
    k_e_all = int(((tab['region1'] == 'eliminated') & (tab['dN_dS'] > thr)).sum())
    n_r_all = int((tab['region1'] == 'retained').sum())
    k_r_all = int(((tab['region1'] == 'retained') & (tab['dN_dS'] > thr)).sum())
    r0 = b[pd.to_numeric(b['min_ds']) == 0]
    cnt_ok = (len(r0) == 1 and int(r0['n_eliminated'].iloc[0]) == n_e_all
              and int(r0['k_eliminated'].iloc[0]) == k_e_all
              and int(r0['n_retained'].iloc[0]) == n_r_all
              and int(r0['k_retained'].iloc[0]) == k_r_all)
    ck.add(T, f'{S}B counts without a dS floor equal the {S1} counts', cnt_ok,
           f'{S1}: {k_e_all}/{n_e_all} eliminated, {k_r_all}/{n_r_all} retained')

    cols_a = [
        Col('gene_ascaris', 'gene1', 'Ascaris gene model', 'text'),
        Col('gene_parascaris', 'gene2', 'P. univalens gene model', 'text'),
        Col('dN', 'dN', 'Pairwise dN (yn00)', 'f4'),
        Col('dS', 'dS', 'Pairwise dS (yn00)', 'f4'),
        Col('dN_dS', 'dN_dS', 'Pairwise dN/dS', 'f3'),
        Col('alignment_length_codons', 'align_len', 'Pairwise codon alignment length, '
            'codons', 'int'),
        Col('chromosome_ascaris', 'chrom1', 'Ascaris germline sequence', 'text'),
        Col('start_ascaris', 'start1', 'Start (1-based)', 'int'),
        Col('end_ascaris', 'end1', 'End (1-based, inclusive)', 'int'),
        Col('region_ascaris', 'region1', 'Region in the Ascaris annotation', 'text',
            values=True),
        Col('eliminated_tract_type_ascaris', 'elim_type1', 'Eliminated tract type: end '
            '(terminal) or internal', 'text', optional=True, values=True),
        Col('region_consensus', 'consensus_region', 'Region where the Ascaris and '
            'P. univalens annotations agree, else discordant', 'text', values=True),
        Col('signed_distance_to_CBR_ascaris_bp', 'signed_cbr_dist1', 'Distance from the '
            'gene midpoint to the nearest CBR, bp, positive on the eliminated side', 'dist'),
        Col('CBR_proximal_ascaris', 'cbr_proximal1', 'Gene midpoint within 100 kb of a CBR',
            'bool'),
    ]
    cols_b = [] if dup else [Col('analysis', 'analysis', 'Analysis as labelled by '
                                 'positive_selection.py', 'text', values=True)]
    cols_b += [
        Col('dS_floor', 'min_ds', 'Pairs with dS below this value excluded (0: none)', 'g'),
        Col('n_eliminated', 'n_eliminated', 'Eliminated pairs', 'int'),
        Col('n_eliminated_dN_dS_above_1', 'k_eliminated', 'Eliminated pairs with '
            f'dN/dS > {thr:g}', 'int'),
        Col('pct_eliminated', 'pct_eliminated', 'Percentage of eliminated pairs above 1',
            'f2'),
        Col('n_retained', 'n_retained', 'Retained pairs', 'int'),
        Col('n_retained_dN_dS_above_1', 'k_retained', f'Retained pairs with dN/dS > {thr:g}',
            'int'),
        Col('pct_retained', 'pct_retained', 'Percentage of retained pairs above 1', 'f2'),
        Col('odds_ratio', 'odds_ratio', 'Odds of dN/dS > 1, eliminated relative to retained'
            + {'conditional': ' (conditional maximum-likelihood estimate)',
               'haldane': ', with 0.5 added to every cell', 'unknown': ''}[or_method], 'f3'),
        Col('odds_ratio_CI_low', 'or_lo95', 'Lower 95% confidence bound'
            + {'conditional': ' (exact)', 'haldane': ' (Woolf)', 'unknown': ''}[or_method],
            'f3'),
        Col('odds_ratio_CI_high', 'or_hi95', 'Upper 95% confidence bound', 'f3'),
        Col('P', 'p', "Two-sided Fisher's exact test", 'p'),
    ]
    drop_b = {'threshold': f'constant ({thr:g}; in the caption)',
              'prop_eliminated': 'shown as a percentage',
              'prop_retained': 'shown as a percentage',
              'n_total': 'n_eliminated + n_retained',
              'or_method': 'which implementation computed the estimator (named in the '
                           'caption and checked against the counts)'}
    if dup:
        drop_b['analysis'] = "the 'primary' row repeats the no-floor row (checked)"
    qv = pd.to_numeric(b.get('q', pd.Series(dtype=float)), errors='coerce')
    if qv.notna().any() and not (qv.dropna() >= 1).all():
        cols_b.append(Col('q', 'q', 'Benjamini–Hochberg q over the dS floors', 'p'))
    else:
        drop_b['q'] = 'equal to 1 in every row'
    og = og.sort_values('dN_dS', ascending=False, kind='mergesort').reset_index(drop=True)
    n = len(og)
    title = f'Ascaris × P. univalens orthologue pairs with pairwise dN/dS above {thr:g}.'
    or_txt = {'conditional': 'Odds ratios compare eliminated with retained pairs: '
                             'conditional maximum-likelihood estimates with exact 95% '
                             'confidence intervals.',
              'haldane': 'Odds ratios compare eliminated with retained pairs, with 0.5 '
                         'added to every cell (Haldane–Anscombe correction), and carry '
                         'Woolf 95% confidence intervals.',
              'unknown': 'Odds ratios compare eliminated with retained pairs.'}[or_method]
    caption = (
        f'(A) The {n} pairs of supplementary table {NUMBER["rates_pairwise"]} with '
        f'dN/dS > {thr:g}, in decreasing order of dN/dS. (B) Their distribution between '
        "regions (two-sided Fisher's exact test), repeated after excluding pairs below "
        'increasing dS floors, where pairwise dN/dS is least reliable. ' + or_txt)
    notes = [f"{int(r['k_eliminated'])}/{int(r['n_eliminated'])} eliminated vs "
             f"{int(r['k_retained'])}/{int(r['n_retained'])} retained above {thr:g} at dS "
             f"floor {r['min_ds']:g}: OR {fmt_num(r['odds_ratio'], 2)} "
             f"[{fmt_num(r['or_lo95'], 2)}, {fmt_num(r['or_hi95'], 2)}], P {p_words(r['p'])}"
             + ('' if is_cond else f"; conditional MLE {fmt_num(c['or'], 2)} "
                                   f"[{fmt_num(c['lo'], 2)}, {fmt_num(c['hi'], 2)}]")
             for (_, r), (_, c) in zip(b.iterrows(), cond.iterrows())]
    parts = [Part('A', f'(A) Pairs with dN/dS > {thr:g}.', og, cols_a,
                  drop={'Ka_Ks': 'renamed'}, sources=('omega_genes',)),
             Part('B', '(B) Distribution between regions, with dS floors.', b, cols_b,
                  drop=drop_b, sources=('omega_stats',))]
    return Table(T, title, caption, parts, notes)


# ---------------------------------------------------------------------------
# S4  Codon usage
# ---------------------------------------------------------------------------

CODON_SPECIES = [('ascaris', 'Ascaris'), ('parascaris', 'P. univalens'),
                 ('toxocara', 'T. canis')]
#: (per-gene column, label), in table order. GC3_all is not in the summary file.
CODON_METRICS = [('GC3s', 'GC3s'), ('GC3_all', 'GC3 (all codons)'), ('ENC', 'ENC'),
                 ('dENC', 'ΔENC'), ('CAI', 'CAI')]
ENC_CEILING = 61.0


def read_codon_summary(path):
    """codon_usage_summary.tsv -> {metric: {'retained': {...}, 'eliminated': {...},
    'U': float, 'p': float}}.

    The test row reuses the columns of the region rows. codon_usage_bias.py
    writes it as n empty, P in `median` and U in `mean`; a table with P in `n`
    and U in `median` is read too. Either way P must lie in [0, 1], and the S4
    check then recomputes both from the per-gene table.
    """
    df = pd.read_csv(path, sep='\t')
    need(df, ['metric', 'region', 'n', 'median', 'mean', 'Q1', 'Q3'], path)
    out = {}
    for _, r in df.iterrows():
        d = out.setdefault(str(r['metric']), {})
        reg = str(r['region'])
        if reg.startswith('MWU'):
            v = {c: pd.to_numeric(pd.Series([r[c]]), errors='coerce').iloc[0]
                 for c in ('n', 'median', 'mean')}
            if np.isfinite(v['median']) and np.isfinite(v['mean']):
                p, U = float(v['median']), float(v['mean'])
            elif np.isfinite(v['n']) and np.isfinite(v['median']):
                p, U = float(v['n']), float(v['median'])
            else:
                sys.exit(f"ERROR: {path}: cannot find P and U in the {r['metric']} test row")
            if not 0.0 <= p <= 1.0:
                sys.exit(f"ERROR: {path}: {r['metric']} test row gives P = {p}, outside "
                         '[0, 1]; the column layout is not the expected one')
            d['p'], d['U'] = p, U
        else:
            d[reg] = {'n': int(float(r['n'])), 'median': float(r['median']),
                      'Q1': float(r['Q1']), 'Q3': float(r['Q3'])}
    return out


def read_codon_report(path):
    """The ceiling-sensitivity and GC3s-matched blocks of codon_usage_summary.txt."""
    with open(path, encoding='utf-8') as fh:
        text = fh.read()
    rep = {'sens': {}, 'matched': {}}
    m = re.search(r'Sensitivity: genes at the ENC ceiling excluded\n(.*?)\n\s*\n', text, re.S)
    if m:
        for line in m.group(1).splitlines():
            mm = re.match(r'^\s*(ENC|dENC)\s+([-+0-9.]+)\s+([-+0-9.]+)\s+(.+?)\s+\(n = '
                          r'([\d,]+) / ([\d,]+)\)\s*$', line)
            if mm:
                rep['sens'][mm.group(1)] = {
                    'ret_median': float(mm.group(2)), 'elim_median': float(mm.group(3)),
                    'p': parse_p(mm.group(4)), 'n_ret': int(mm.group(5).replace(',', '')),
                    'n_elim': int(mm.group(6).replace(',', ''))}
    mc = re.search(r'within ([0-9.]+) GC3s, repeated (\d+) times', text)
    if mc:
        rep['caliper'], rep['replicates'] = float(mc.group(1)), int(mc.group(2))
    for line in text.splitlines():
        mm = re.match(r'^\s*(CAI|ENC|dENC)\s+([-+]?[0-9.]+)\s+\[\s*([-+]?[0-9.]+),\s*'
                      r'([-+]?[0-9.]+)\]\s+(.+?)\s+(\d+)%\s+([0-9.]+)\s*$', line)
        if mm:
            rep['matched'][mm.group(1)] = {
                'diff': float(mm.group(2)), 'lo': float(mm.group(3)),
                'hi': float(mm.group(4)), 'p': parse_p(mm.group(5)),
                'pct_sig': float(mm.group(6)), 'gc3s_p': float(mm.group(7))}
    mn = re.search(r'Matched n = ([\d,]+) per group', text)
    if mn:
        rep['matched_n'] = int(mn.group(1).replace(',', ''))
    ms = re.search(r'surviving matching:\s*(.+?)\.\s*$', text, re.M)
    if ms:
        rep['share'] = {k: float(v) for k, v in
                        re.findall(r'(CAI|ENC|dENC)\s+(-?\d+(?:\.\d+)?)%', ms.group(1))}
    return rep


def build_codon_usage(ctx):
    T = 'codon_usage'
    ck = ctx.checks
    rows_a, rows_b, rows_c = [], [], []
    calipers, reps = set(), set()
    notes = []
    for sp, lab in CODON_SPECIES:
        summ = read_codon_summary(ctx.src[f'codon_summary_{sp}'])
        rep = read_codon_report(ctx.src[f'codon_report_{sp}'])
        g = ctx.tsv(f'codon_genes_{sp}')
        need(g, ['gene_id', 'region', 'ENC', 'GC3s', 'GC3_all', 'ENC_at_ceiling', 'dENC',
                 'CAI'], ctx.src[f'codon_genes_{sp}'])
        regs = sorted(g['region'].dropna().unique())
        ck.add(T, f'{lab}: per-gene regions are retained and eliminated only',
               regs == ['eliminated', 'retained'], ', '.join(map(str, regs)))
        ret_g, elim_g = g[g['region'] == 'retained'], g[g['region'] == 'eliminated']
        # --- A: every metric; summary-file values, verified against the per-gene table
        for col, mlab in CODON_METRICS:
            dr, de = describe(ret_g[col]), describe(elim_g[col])
            U, p = mwu(ret_g[col].dropna(), elim_g[col].dropna())
            if col in summ:
                s = summ[col]
                bad = [f'{r} {k}' for r, d in (('retained', dr), ('eliminated', de))
                       for k in ('n', 'median', 'Q1', 'Q3')
                       if not close(s[r][k], d[k], rtol=PRINTED_RTOL).all()]
                if not close(s['U'], U).all():
                    bad.append('U')
                if not close(s['p'], p, rtol=PRINTED_RTOL).all():
                    bad.append(f"P ({p_words(s['p'])} vs {p_words(p)})")
                ck.add(T, f'{lab} {mlab}: summary reproduced from the per-gene table',
                       not bad, ('differs: ' + ', '.join(bad)) if bad else
                       f"P {p_words(p)}")
                dr = dict(dr, **{k: s['retained'][k] for k in ('median', 'Q1', 'Q3')})
                de = dict(de, **{k: s['eliminated'][k] for k in ('median', 'Q1', 'Q3')})
                U, p = s['U'], s['p']
            elif col != 'GC3_all':
                sys.exit(f"ERROR: {ctx.src[f'codon_summary_{sp}']}: no {col} rows")
            rows_a.append({'species': lab, 'metric': mlab, 'n_ret': dr['n'],
                           'n_elim': de['n'], 'med_ret': dr['median'], 'q1_ret': dr['Q1'],
                           'q3_ret': dr['Q3'], 'med_elim': de['median'],
                           'q1_elim': de['Q1'], 'q3_elim': de['Q3'], 'U': U, 'p': p})
            if col == 'GC3s':
                notes.append(f"{lab}: GC3s {de['median']:.3f} vs {dr['median']:.3f}, "
                             f"P {p_words(p)}")
        # --- B: the ENC ceiling
        at = to_bool(g['ENC_at_ceiling'])
        enc = pd.to_numeric(g['ENC'], errors='coerce')
        ck.add(T, f'{lab}: ENC_at_ceiling flags exactly the genes at ENC = 61',
               bool(((enc >= ENC_CEILING - 1e-9) == at.fillna(False)).all()))
        n_at_r = int(at[g['region'] == 'retained'].sum())
        n_at_e = int(at[g['region'] == 'eliminated'].sum())
        below = g[~at.fillna(False).astype(bool)]
        for col, mlab in (('ENC', 'ENC'), ('dENC', 'ΔENC')):
            r_v = below.loc[below['region'] == 'retained', col].dropna()
            e_v = below.loc[below['region'] == 'eliminated', col].dropna()
            U, p = mwu(r_v, e_v)
            s = rep['sens'].get(col)
            if s is None:
                ck.add(T, f'{lab} {mlab}: ceiling-excluded result found in the report', False)
            else:
                ok = (s['n_ret'] == len(r_v) and s['n_elim'] == len(e_v)
                      and abs(s['ret_median'] - r_v.median()) <= 5.1e-5
                      and abs(s['elim_median'] - e_v.median()) <= 5.1e-5
                      and close(s['p'], p, rtol=1e-2, atol=1e-3 if s['p'] >= 1e-3 else 0).all())
                ck.add(T, f'{lab} {mlab}: ceiling-excluded test reproduces the report', ok,
                       f"P {p_words(s['p'])} vs {p_words(p)}")
            rows_b.append({'species': lab, 'metric': mlab, 'n_all_ret': len(ret_g),
                           'n_at_ret': n_at_r, 'pct_at_ret': 100 * n_at_r / len(ret_g),
                           'n_all_elim': len(elim_g), 'n_at_elim': n_at_e,
                           'pct_at_elim': 100 * n_at_e / len(elim_g),
                           'n_ret': len(r_v), 'n_elim': len(e_v), 'med_ret': r_v.median(),
                           'med_elim': e_v.median(), 'U': U, 'p': p})
        # --- C: GC3s-matched (resampling results, from the report)
        if 'caliper' in rep:
            calipers.add(rep['caliper'])
            reps.add(rep['replicates'])
        for col, mlab in (('CAI', 'CAI'), ('ENC', 'ENC'), ('dENC', 'ΔENC')):
            mt = rep['matched'].get(col)
            if mt is None:
                sys.exit(f"ERROR: {ctx.src[f'codon_report_{sp}']}: no GC3s-matched row "
                         f'for {col}')
            raw = summ[col]['eliminated']['median'] - summ[col]['retained']['median']
            share = rep.get('share', {}).get(col, np.nan)
            ck.add(T, f'{lab} {mlab}: reported share surviving matching = matched / raw '
                   'difference', abs(100 * mt['diff'] / raw - share) <= 1.0 if raw else False,
                   f"{100 * mt['diff'] / raw:.1f}% vs {share:g}%")
            rows_c.append({'species': lab, 'metric': mlab,
                           'n_matched': rep.get('matched_n', np.nan),
                           'raw_diff': raw, 'diff': mt['diff'], 'lo': mt['lo'],
                           'hi': mt['hi'], 'p': mt['p'], 'pct_sig': mt['pct_sig'],
                           'gc3s_p': mt['gc3s_p'], 'share': share})
        ck.add(T, f'{lab}: matched n equals the number of eliminated genes',
               rep.get('matched_n') == len(elim_g), f"{rep.get('matched_n')} vs {len(elim_g)}")
    ck.add(T, 'one caliper and one replicate count across species',
           len(calipers) == 1 and len(reps) == 1,
           f'calipers {sorted(calipers)}, replicates {sorted(reps)}')
    caliper = next(iter(calipers)) if len(calipers) == 1 else np.nan
    nrep = next(iter(reps)) if len(reps) == 1 else 0

    stat_cols = [
        Col('n_retained', 'n_ret', 'Retained genes', 'int'),
        Col('n_eliminated', 'n_elim', 'Eliminated genes', 'int'),
    ]
    cols_a = [Col('species', 'species', 'Species; genes labelled by its own annotation',
                  'text'),
              Col('metric', 'metric', 'Codon-usage statistic', 'text', values=True)] \
        + stat_cols + [
        Col('median_retained', 'med_ret', 'Median, retained genes', 'f4'),
        Col('Q1_retained', 'q1_ret', 'First quartile, retained', 'f4'),
        Col('Q3_retained', 'q3_ret', 'Third quartile, retained', 'f4'),
        Col('median_eliminated', 'med_elim', 'Median, eliminated genes', 'f4'),
        Col('Q1_eliminated', 'q1_elim', 'First quartile, eliminated', 'f4'),
        Col('Q3_eliminated', 'q3_elim', 'Third quartile, eliminated', 'f4'),
        Col('U', 'U', 'Mann–Whitney U of the retained sample', 'g'),
        Col('P', 'p', 'Two-sided Mann–Whitney U test', 'p')]
    cols_b = [
        Col('species', 'species', 'Species', 'text'),
        Col('metric', 'metric', 'Statistic tested after excluding genes at the ceiling',
            'text', values=True),
        Col('n_retained_all', 'n_all_ret', 'Retained genes analysed', 'int'),
        Col('n_retained_at_ceiling', 'n_at_ret', 'Retained genes with ENC at the '
            'ceiling of 61', 'int'),
        Col('pct_retained_at_ceiling', 'pct_at_ret', 'Percentage at the ceiling', 'f1'),
        Col('n_eliminated_all', 'n_all_elim', 'Eliminated genes analysed', 'int'),
        Col('n_eliminated_at_ceiling', 'n_at_elim', 'Eliminated genes at the ceiling',
            'int'),
        Col('pct_eliminated_at_ceiling', 'pct_at_elim', 'Percentage at the ceiling', 'f1'),
        Col('n_retained', 'n_ret', 'Retained genes below the ceiling', 'int'),
        Col('n_eliminated', 'n_elim', 'Eliminated genes below the ceiling', 'int'),
        Col('median_retained', 'med_ret', 'Median, retained genes below the ceiling', 'f4'),
        Col('median_eliminated', 'med_elim', 'Median, eliminated genes below the ceiling',
            'f4'),
        Col('U', 'U', 'Mann–Whitney U of the retained sample', 'g'),
        Col('P', 'p', 'Two-sided Mann–Whitney U test', 'p')]
    cols_c = [
        Col('species', 'species', 'Species', 'text'),
        Col('metric', 'metric', 'Statistic compared on GC3s-matched gene sets', 'text',
            values=True),
        Col('n_matched_per_group', 'n_matched', 'Genes per matched set (every eliminated '
            'gene, each with one retained partner)', 'int'),
        Col('difference_unmatched', 'raw_diff', 'Eliminated minus retained median, all '
            'genes', 'f4'),
        Col('difference_matched_median', 'diff', 'Eliminated minus retained median on the '
            'matched sets, median over replicates', 'f4'),
        Col('difference_matched_2.5pct', 'lo', '2.5th percentile over replicates', 'f4'),
        Col('difference_matched_97.5pct', 'hi', '97.5th percentile over replicates', 'f4'),
        Col('P_median', 'p', 'Two-sided Mann–Whitney P on the matched sets, median over '
            'replicates', 'p'),
        Col('pct_replicates_P_below_0.05', 'pct_sig', 'Percentage of replicates with '
            'P < 0.05', 'f1'),
        Col('GC3s_balance_P', 'gc3s_p', 'Mann–Whitney P for GC3s between the matched sets '
            '(median over replicates); far from significance when the matching worked', 'f2'),
        Col('pct_difference_surviving', 'share', 'Matched difference as a percentage of '
            'the unmatched difference', 'f1'),
    ]
    title = 'Codon usage of retained and eliminated genes in Ascaris, P. univalens and T. canis.'
    caption = (
        "Genes of each species are labelled by that species' own elimination annotation; "
        'genes with fewer than 100 codons, or lacking any of the two-, four- and six-fold '
        'degenerate amino-acid classes from which ENC is estimated, are excluded. GC3s, GC '
        'content at synonymous third '
        'positions (every sense codon except ATG and TGG); GC3, the same with ATG and TGG '
        'included; ENC, effective number of codons (Wright 1990), whose ceiling of 61 is '
        'reached by genes with no detectable bias; ΔENC, the ENC expected from GC3s minus '
        'the observed ENC, so that positive values denote more bias than composition '
        'predicts; CAI, codon adaptation index (Sharp and Li 1987), computed for Ascaris '
        'against a reference set defined by somatic expression and for P. univalens and '
        'T. canis against the most codon-biased retained genes, which makes the contrast '
        'partly circular in those two species, so it is not interpreted there; reference '
        'sets differ, so CAI is not compared across species. (A) Medians and quartiles, '
        'two-sided Mann–Whitney U tests. (B) Genes at the ENC ceiling, and ENC and ΔENC '
        'tested after excluding them. (C) GC3s-matched comparison: each eliminated gene is '
        f'paired with an unused retained gene within {caliper:g} GC3s, and the matching is '
        f'repeated {nrep} times.')
    parts = [
        Part('A', '(A) All analysed genes.', pd.DataFrame(rows_a), cols_a,
             sources=tuple(f'codon_{k}_{sp}' for sp, _ in CODON_SPECIES
                           for k in ('summary', 'genes'))),
        Part('B', '(B) The ENC ceiling.', pd.DataFrame(rows_b), cols_b,
             sources=tuple(f'codon_{k}_{sp}' for sp, _ in CODON_SPECIES
                           for k in ('genes', 'report'))),
        Part('C', '(C) GC3s-matched comparison.', pd.DataFrame(rows_c), cols_c,
             sources=tuple(f'codon_report_{sp}' for sp, _ in CODON_SPECIES)),
    ]
    return Table(T, title, caption, parts, notes)


# ---------------------------------------------------------------------------
# S9  Insulation and dN/dS models (S5 before the 2026-10-12 renumbering)
# ---------------------------------------------------------------------------

MODEL_RUNS = [('models_noEV', 'none'), ('models_EVtestis', 'testis'),
              ('models_EV4cell', '4-cell embryo (60 h)')]
#: Part C, one run per Hi-C stage, each without an eigenvector and with that
#: stage's insulation track and expression:
#: (stage label, models key, per-gene key, summary key, token in the insulation
#:  track's file name, expression column of the workbook).
#: Testis is the no-eigenvector run of parts A and B (models_noEV); the tokens
#: and columns are those of run_multivariate_insulation.sh.
STAGE_RUNS = [
    ('Testis', 'models_noEV', 'stage_pergene_testis', 'stage_summary_testis',
     'teste', 'testis'),
    ('Ovary', 'stage_models_ovary', 'stage_pergene_ovary', 'stage_summary_ovary',
     'ovary', 'ovary'),
    ('Zygote (0 h)', 'stage_models_zygote', 'stage_pergene_zygote',
     'stage_summary_zygote', '0hr', 'zygote'),
    ('48 h embryo', 'stage_models_48h', 'stage_pergene_48h', 'stage_summary_48h',
     '48hr', 'early_embryo'),
    ('4-cell embryo (60 h)', 'stage_models_60h', 'stage_pergene_60h',
     'stage_summary_60h', '60hr', 'early_embryo'),
]
EXPR_LABEL = {'testis': 'testis', 'ovary': 'ovary', 'zygote': 'zygote',
              'early_embryo': 'early embryo'}
#: Per-gene columns the core-set rule and the region comparison need.
PERGENE_COLS = ['gene', 'region', 'dN_dS', 'insulation', 'log_expr', 'gc']
TERM_LABEL = {'intercept': 'Intercept', 'insulation': 'Insulation score',
              'log_expr': 'Expression, log2(x + 1)', 'gc': 'Local GC content',
              'log_intergenic': 'Nearest intergenic distance, log10(d + 1)',
              'ev1': 'Compartment eigenvector (EV1)', 'is_eliminated': 'Eliminated (0/1)'}
RESPONSE_LABEL = {'dN_dS': 'dN/dS', 'insulation': 'Insulation score'}
MODEL_COLS = ['model', 'response', 'term', 'beta', 'se', 't', 'p', 'ci_lo', 'ci_hi', 'n',
              'n_retained', 'n_eliminated', 'r2', 'df']
_NUM = r'([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)'


def same_fit(a, b):
    """Two model blocks (rows of one model) are the same fit."""
    if len(a) != len(b) or list(a['term']) != list(b['term']):
        return False
    return all(close(a[c].values, b[c].values, rtol=1e-9).all()
               for c in ('beta', 'se', 'p', 'n', 'r2'))


def kept_models(d):
    """Model names of one run to publish, in file order, and the refits on the
    intergenic-distance subset left out because they equal the base fit."""
    names = list(dict.fromkeys(d['model']))
    keep, dropped = [], []
    for mname in names:
        mm = re.match(r'^(Model \d), IG set$', mname)
        if mm and same_fit(d[d['model'] == mname], d[d['model'] == mm.group(1)]):
            dropped.append(mname)
            continue
        keep.append(mname)
    n_ig = sum(bool(re.match(r'^Model \d, IG set$', m)) for m in names)
    return keep, dropped, n_ig


def read_mv_summary(path):
    """The header and the region comparison of a multivariate_insulation.py
    summary: the insulation track file, the expression column, the eigenvector
    line, and the printed medians, counts and P (used to check the values
    recomputed here, never published). Keys are None when a line is absent."""
    with open(path, encoding='utf-8') as fh:
        text = fh.read()
    out = {'track': None, 'expr_col': None, 'eigenvector': None, 'med_ret': None,
           'n_ret': None, 'med_elim': None, 'n_elim': None, 'p': None}
    m = re.search(r'^Insulation\s+(\S+)', text, re.M)
    if m:
        out['track'] = m.group(1)
    m = re.search(r"^Expression\s+\S+ column '([^']+)'", text, re.M)
    if m:
        out['expr_col'] = m.group(1)
    m = re.search(r'^Eigenvector\s+(.+?)\s*$', text, re.M)
    if m:
        out['eigenvector'] = m.group(1)
    m = re.search(r'retained median ' + _NUM + r' \(n = ([\d,]+)\); eliminated median '
                  + _NUM + r' \(n = ([\d,]+)\); P = ' + _NUM, text)
    if m:
        out.update(med_ret=float(m.group(1)), n_ret=int(m.group(2).replace(',', '')),
                   med_elim=float(m.group(3)), n_elim=int(m.group(4).replace(',', '')),
                   p=float(m.group(5)))
    return out


def model_term(d, model, term, label):
    """The single row of one term of one model, or stop."""
    r = d[(d['model'] == model) & (d['term'] == term)]
    if len(r) != 1:
        sys.exit(f'ERROR: {label}: {model} has {len(r)} {term!r} rows, expected 1')
    return r.iloc[0]


def stage_comparison(ctx, T, stage, mkey, gkey, skey, tok, expr, mods):
    """Part C row for one Hi-C stage, with its checks. Returns (row, core genes
    as a gene -> dN/dS frame)."""
    ck = ctx.checks
    summ = read_mv_summary(ctx.src[skey])
    ck.add(T, f'{stage}: the run used the {tok} insulation track',
           bool(summ['track']) and f'.{tok}.' in summ['track'],
           summ['track'] or 'no Insulation line in the summary')
    ck.add(T, f"{stage}: the run used expression column '{expr}'",
           summ['expr_col'] == expr, f"summary: {summ['expr_col']!r}")
    ck.add(T, f'{stage}: the run has no eigenvector',
           (summ['eigenvector'] or '').startswith('none'),
           summ['eigenvector'] or 'no Eigenvector line in the summary')
    ev_models = sorted(m for m in set(mods['model']) if m.endswith('+EV1'))
    ck.add(T, f'{stage}: no +EV1 model in the run', not ev_models, ', '.join(ev_models))

    g = ctx.tsv(gkey)
    need(g, PERGENE_COLS, ctx.src[gkey])
    dnds = pd.to_numeric(g['dN_dS'], errors='coerce')
    core_mask = (g['region'].isin(['retained', 'eliminated']) & np.isfinite(dnds)
                 & g['insulation'].notna() & g['log_expr'].notna() & g['gc'].notna())
    core = g.loc[core_mask].assign(dN_dS=dnds[core_mask])
    ret = core.loc[core['region'] == 'retained', 'insulation'].to_numpy(float)
    elim = core.loc[core['region'] == 'eliminated', 'insulation'].to_numpy(float)
    if len(ret) == 0 or len(elim) == 0:
        sys.exit(f'ERROR: {ctx.src[gkey]}: a region is empty in the core set')
    _u, p = mwu(ret, elim)
    med_r, med_e = float(np.median(ret)), float(np.median(elim))

    e2 = model_term(mods, 'Model 2', 'is_eliminated', ctx.src[mkey])
    e1 = model_term(mods, 'Model 1', 'is_eliminated', ctx.src[mkey])
    i1 = model_term(mods, 'Model 1', 'insulation', ctx.src[mkey])
    for nm, e in (('Model 1', e1), ('Model 2', e2)):
        ck.add(T, f'{stage}: core-set counts equal those of {nm}',
               int(e['n_retained']) == len(ret) and int(e['n_eliminated']) == len(elim),
               f"per-gene {len(ret):,} / {len(elim):,}; model "
               f"{int(e['n_retained']):,} / {int(e['n_eliminated']):,}")
    if summ['p'] is None:
        ck.add(T, f'{stage}: region comparison found in the summary', False,
               ctx.src[skey])
    else:
        # The summary prints medians to 4 decimals and P to 3 significant digits.
        ok = (summ['n_ret'] == len(ret) and summ['n_elim'] == len(elim)
              and abs(med_r - summ['med_ret']) <= 5.1e-5
              and abs(med_e - summ['med_elim']) <= 5.1e-5
              and bool(close(p, summ['p'], rtol=6e-3)[0]))
        ck.add(T, f'{stage}: recomputed medians, counts and P match the summary', ok,
               f"medians {med_r:+.4f} / {med_e:+.4f}, P {p:.3g}; printed "
               f"{summ['med_ret']:+.4f} / {summ['med_elim']:+.4f}, P {summ['p']:.3g}")
    row = {'stage': stage, 'expr': EXPR_LABEL.get(expr, expr),
           'n_ret': len(ret), 'n_elim': len(elim),
           'med_ret': med_r, 'med_elim': med_e, 'p_mwu': p,
           'b2': e2['beta'], 'lo2': e2['ci_lo'], 'hi2': e2['ci_hi'], 'p2': e2['p'],
           'r2_2': e2['r2'],
           'b1': e1['beta'], 'lo1': e1['ci_lo'], 'hi1': e1['ci_hi'], 'p1': e1['p'],
           'bi1': i1['beta'], 'pi1': i1['p'], 'r2_1': e1['r2'], 'n_fit': int(e1['n'])}
    return row, core[['gene', 'dN_dS']]


def build_insulation_models(ctx):
    T = 'insulation_models'
    ck = ctx.checks
    runs = {}
    for key, spec in MODEL_RUNS:
        d = ctx.tsv(key)
        need(d, MODEL_COLS, ctx.src[key])
        d['spec'] = spec
        runs[key] = d
    base = runs['models_noEV']
    base_models = list(dict.fromkeys(base['model']))
    for need_m in ('Model 1', 'Model 2'):
        if need_m not in base_models:
            sys.exit(f"ERROR: {ctx.src['models_noEV']}: no '{need_m}' rows")
    # Models fitted without an eigenvector must be identical in every run.
    for key, _spec in MODEL_RUNS[1:]:
        d = runs[key]
        for mname in base_models:
            if mname in set(d['model']):
                ck.add(T, f"'{mname}' identical in the {key} and models_noEV runs",
                       same_fit(base[base['model'] == mname], d[d['model'] == mname]))
        for mname in ('Model 1 +EV1', 'Model 2 +EV1'):
            if mname not in set(d['model']):
                sys.exit(f'ERROR: {ctx.src[key]}: no {mname!r} rows; was it run with '
                         '--eigenvector?')
    # A refit on the intergenic subset that equals the base fit is not repeated.
    keep_models, dropped_ig, n_ig = kept_models(base)
    ck.add(T, 'refits on the intergenic-distance subset equal the base fits',
           len(dropped_ig) == n_ig, f"dropped: {', '.join(dropped_ig) or 'none'}")

    # The four further stage runs (testis is models_noEV).
    stage_runs = {}
    for stage, mkey, _g, _s, _tok, _expr in STAGE_RUNS[1:]:
        d = ctx.tsv(mkey)
        need(d, MODEL_COLS, ctx.src[mkey])
        d['spec'] = 'none'
        for need_m in ('Model 1', 'Model 2'):
            if need_m not in set(d['model']):
                sys.exit(f"ERROR: {ctx.src[mkey]}: no '{need_m}' rows")
        keep_s, dropped_s, n_ig_s = kept_models(d)
        ck.add(T, f'{stage}: refits on the intergenic-distance subset equal the base fits',
               len(dropped_s) == n_ig_s, f"dropped: {', '.join(dropped_s) or 'none'}")
        stage_runs[mkey] = (d, keep_s)

    # --- A: the elimination coefficient across compartment specifications -------
    rows = []
    for key, spec in MODEL_RUNS:
        d = runs[key]
        for resp, mname in (('insulation', 'Model 2'), ('dN_dS', 'Model 1')):
            mod = mname if spec == 'none' else f'{mname} +EV1'
            blk = d[d['model'] == mod]
            e = blk[blk['term'] == 'is_eliminated']
            if len(e) != 1:
                sys.exit(f'ERROR: {ctx.src[key]}: {mod} has no single is_eliminated row')
            e = e.iloc[0]
            ev = blk[blk['term'] == 'ev1']
            rows.append({'spec': spec, 'response': RESPONSE_LABEL.get(resp, resp),
                         'model': mod, 'n': e['n'], 'n_ret': e['n_retained'],
                         'n_elim': e['n_eliminated'], 'r2': e['r2'], 'beta': e['beta'],
                         'se': e['se'], 'lo': e['ci_lo'], 'hi': e['ci_hi'], 'p': e['p'],
                         'beta_ev': ev['beta'].iloc[0] if len(ev) else np.nan,
                         'p_ev': ev['p'].iloc[0] if len(ev) else np.nan})
            if spec != 'none' and len(ev) != 1:
                ck.add(T, f'{mod} in {key} carries an ev1 term', False)
    a = pd.DataFrame(rows)
    a = a.assign(_r=a['response'].map({'Insulation score': 0, 'dN/dS': 1})).sort_values(
        ['_r'], kind='mergesort').reset_index(drop=True)

    # --- B: every coefficient, testis runs first, then the other stages ------------
    blocks = [base[base['model'].isin(keep_models)].assign(stage=STAGE_RUNS[0][0])]
    for key, _spec in MODEL_RUNS[1:]:
        d = runs[key]
        blocks.append(d[d['model'].str.endswith('+EV1')].assign(stage=STAGE_RUNS[0][0]))
    for stage, mkey, *_rest in STAGE_RUNS[1:]:
        d, keep_s = stage_runs[mkey]
        blocks.append(d[d['model'].isin(keep_s)].assign(stage=stage))
    b = pd.concat(blocks, ignore_index=True)
    b['term_label'] = b['term'].map(TERM_LABEL).fillna(b['term'])
    b['response_label'] = b['response'].map(RESPONSE_LABEL).fillna(b['response'])
    unknown = sorted(set(b['term']) - set(TERM_LABEL))
    ck.add(T, 'every model term has a label', not unknown, ', '.join(unknown))

    # --- C: the five Hi-C stages ------------------------------------------------------
    rows_c, cores = [], {}
    for stage, mkey, gkey, skey, tok, expr in STAGE_RUNS:
        mods = base if mkey == 'models_noEV' else stage_runs[mkey][0]
        row, core_genes = stage_comparison(ctx, T, stage, mkey, gkey, skey, tok, expr,
                                           mods)
        rows_c.append(row)
        cores[stage] = core_genes
    c = pd.DataFrame(rows_c)
    # Every stage reads the same rate table: the same core genes, the same dN/dS.
    ref_stage = STAGE_RUNS[0][0]
    ref = cores[ref_stage]
    for stage, *_rest in STAGE_RUNS[1:]:
        cur = cores[stage]
        ck.add(T, f'{stage}: the same core genes as {ref_stage}',
               set(cur['gene']) == set(ref['gene']),
               f"{len(cur):,} and {len(ref):,} genes, "
               f"{len(set(cur['gene']) ^ set(ref['gene'])):,} in one only")
        if ref['gene'].is_unique and cur['gene'].is_unique:
            mg = ref.merge(cur, on='gene', suffixes=('_ref', '_cur'))
            ck.add(T, f'{stage}: dN/dS of every shared gene equals {ref_stage}',
                   bool(close(mg['dN_dS_ref'], mg['dN_dS_cur'], rtol=PRINTED_RTOL).all()),
                   f'{len(mg):,} genes; largest relative difference '
                   f"{max_rel(mg['dN_dS_ref'], mg['dN_dS_cur']):.2g}")
        else:
            ck.add(T, f'{stage}: dN/dS of every shared gene equals {ref_stage}', False,
                   'gene IDs are not unique in a per-gene table; not compared')

    core = base[(base['model'] == 'Model 1') & (base['term'] == 'is_eliminated')].iloc[0]
    n_core, n_ret, n_elim = int(core['n']), int(core['n_retained']), int(core['n_eliminated'])
    ck.add(T, f'every stage fits the {n_core:,} genes of the testis models',
           bool((c['n_fit'] == n_core).all()),
           ', '.join(f"{r['stage']} {int(r['n_fit']):,}" for _, r in c.iterrows()))
    cols_a = [
        Col('compartment_eigenvector', 'spec', 'Pre-elimination compartment eigenvector '
            'added to the model', 'text', values=True),
        Col('response', 'response', 'Response variable', 'text', values=True),
        Col('model', 'model', 'Model (Materials and Methods)', 'text'),
        Col('n', 'n', 'Genes in the fit', 'int'),
        Col('n_retained', 'n_ret', 'Retained genes in the fit', 'int'),
        Col('n_eliminated', 'n_elim', 'Eliminated genes in the fit', 'int'),
        Col('R2', 'r2', 'Coefficient of determination', 'f4'),
        Col('beta_eliminated', 'beta', 'Elimination coefficient: eliminated minus retained, '
            'in standard deviations of the response', 'f4'),
        Col('SE', 'se', 'Standard error of the coefficient', 'f4'),
        Col('CI_low', 'lo', 'Lower 95% bound, β − t(0.975, df) × SE', 'f4'),
        Col('CI_high', 'hi', 'Upper 95% bound', 'f4'),
        Col('P', 'p', 'Two-sided t test', 'p'),
        Col('beta_EV1', 'beta_ev', 'Coefficient of the compartment eigenvector (per SD)',
            'f4'),
        Col('P_EV1', 'p_ev', 'Two-sided t test for the eigenvector coefficient', 'p'),
    ]
    cols_b = [
        Col('hic_stage', 'stage', 'Hi-C stage of the insulation track; expression is '
            'from the same stage (early embryo for the 48 h and 60 h embryos)', 'text',
            values=True),
        Col('compartment_eigenvector', 'spec', 'Eigenvector added (none for models '
            'without one)', 'text', values=True),
        Col('model', 'model', 'Model: 1 (dN/dS), 2 (insulation); +IG adds intergenic '
            'distance, +EV1 the eigenvector', 'text', values=True),
        Col('response', 'response_label', 'Response variable', 'text', values=True),
        Col('term', 'term_label', 'Model term; continuous terms are z-scores within the '
            'fitted set', 'text', values=True),
        Col('beta', 'beta', 'Coefficient', 'f4'),
        Col('SE', 'se', 'Standard error', 'f4'),
        Col('t', 't', 't statistic', 'f3'),
        Col('P', 'p', 'Two-sided t test', 'p'),
        Col('CI_low', 'ci_lo', 'Lower 95% bound', 'f4'),
        Col('CI_high', 'ci_hi', 'Upper 95% bound', 'f4'),
        Col('n', 'n', 'Genes in the fit', 'int'),
        Col('n_retained', 'n_retained', 'Retained genes', 'int'),
        Col('n_eliminated', 'n_eliminated', 'Eliminated genes', 'int'),
        Col('R2', 'r2', 'Coefficient of determination', 'f4'),
        Col('df', 'df', 'Residual degrees of freedom', 'int'),
    ]
    cols_c = [
        Col('hic_stage', 'stage', 'Hi-C stage of the insulation track', 'text', values=True),
        Col('expression_stage', 'expr', 'Expression sample used in both models', 'text',
            values=True),
        Col('n_retained', 'n_ret', 'Retained genes, in the comparison and in both fits',
            'int'),
        Col('n_eliminated', 'n_elim', 'Eliminated genes, in the comparison and in both '
            'fits', 'int'),
        Col('median_insulation_retained', 'med_ret', 'Median insulation score of retained '
            'genes (lower = stronger boundary)', 'f4'),
        Col('median_insulation_eliminated', 'med_elim', 'Median insulation score of '
            'eliminated genes', 'f4'),
        Col('P_Mann_Whitney', 'p_mwu', 'Two-sided Mann-Whitney U, retained against '
            'eliminated insulation scores (unadjusted)', 'p'),
        Col('insulation_model_beta_eliminated', 'b2', 'Model 2 (insulation): elimination '
            'coefficient, in standard deviations of insulation, at equal expression and '
            'GC', 'f4'),
        Col('insulation_model_CI_low', 'lo2', 'Lower 95% bound', 'f4'),
        Col('insulation_model_CI_high', 'hi2', 'Upper 95% bound', 'f4'),
        Col('insulation_model_P', 'p2', 'Two-sided t test', 'p'),
        Col('insulation_model_R2', 'r2_2', 'Coefficient of determination of model 2',
            'f4'),
        Col('dNdS_model_beta_eliminated', 'b1', 'Model 1 (dN/dS): elimination coefficient, '
            'in standard deviations of dN/dS, at equal insulation, expression and GC',
            'f4'),
        Col('dNdS_model_CI_low', 'lo1', 'Lower 95% bound', 'f4'),
        Col('dNdS_model_CI_high', 'hi1', 'Upper 95% bound', 'f4'),
        Col('dNdS_model_P', 'p1', 'Two-sided t test', 'p'),
        Col('dNdS_model_beta_insulation', 'bi1', 'Model 1: insulation coefficient (per '
            'standard deviation of insulation)', 'f4'),
        Col('dNdS_model_P_insulation', 'pi1', 'Two-sided t test for the insulation '
            'coefficient', 'p'),
        Col('dNdS_model_R2', 'r2_1', 'Coefficient of determination of model 1', 'f4'),
    ]
    title = ('Elimination status in regression models of chromatin insulation and dN/dS, '
             'with and without control for compartment, at five developmental stages.')
    ig_note = (' Refits restricted to genes with an intergenic distance are identical to '
               'the base fits, because every gene has one, and are not repeated.'
               if dropped_ig else '')
    caption = (
        f'Ordinary least squares on {n_core:,} Ascaris genes ({n_ret:,} retained, '
        f'{n_elim:,} eliminated) with dN/dS, insulation score, expression and local GC '
        'content. Model 1: dN/dS ~ insulation + expression + GC + eliminated. '
        'Model 2: insulation ~ expression + GC + eliminated. The response and the continuous '
        'predictors are z-scored within each fitted set, so the elimination coefficient is '
        'the eliminated-minus-retained difference in standard deviations of the response; '
        'intervals are β ± t(0.975, df) × SE. Insulation scores are from 100 kb windows on '
        '20 kb bins of pre-elimination Hi-C, and expression is log2(x + 1) in the sample of '
        'the same stage (the early-embryo sample for the 48 h and 60 h embryos); GC is '
        'measured in a 20 kb window centred on each gene. (A) Testis: the elimination '
        'coefficient '
        'without a compartment term, and with the pre-elimination compartment eigenvector '
        '(EV1) from testis or from 4-cell (60 h) embryos added, each fitted on the genes '
        'with an EV1 value; EV1 is oriented on each chromosome to correlate positively with '
        'GC. (B) Every coefficient of every model, including models with nearest intergenic '
        'distance (+IG), which is a mediator rather than a confounder and is therefore '
        'reported both with and without.' + ig_note + ' (C) Each of the five stages, '
        'without a compartment term: the median insulation score of retained and '
        'eliminated genes, compared by two-sided Mann-Whitney U without adjustment (fig. 5D '
        'for testis), and the elimination coefficients of both models and the insulation '
        'coefficient of model 1, each adjusted for the other terms of its model.')
    notes = [f"{r['response']:<16} EV1 {r['spec']:<22} β {r['beta']:+.4f} "
             f"[{r['lo']:+.4f}, {r['hi']:+.4f}]  P {p_words(r['p'])}  n {int(r['n']):,}"
             for _, r in a.iterrows()]
    notes += [f"C {r['stage']:<21} medians {r['med_ret']:+.4f} / {r['med_elim']:+.4f} "
              f"(ret / elim), P {p_words(r['p_mwu'])};  insulation model β "
              f"{r['b2']:+.4f} P {p_words(r['p2'])};  dN/dS model β {r['b1']:+.4f} "
              f"P {p_words(r['p1'])}, insulation {r['bi1']:+.4f} P {p_words(r['pi1'])}"
              for _, r in c.iterrows()]
    stage_sources = tuple(k for _st, mk, gk, sk, _t, _e in STAGE_RUNS
                          for k in (mk, gk, sk))
    parts = [Part('A', '(A) Elimination coefficient by compartment control.', a, cols_a,
                  drop={'_r': 'sort key'}, sources=tuple(k for k, _ in MODEL_RUNS)),
             Part('B', '(B) Every coefficient.', b, cols_b,
                  drop={'term': 'shown as a label', 'response': 'shown as a label'},
                  sources=tuple(k for k, _ in MODEL_RUNS) +
                  tuple(mk for _st, mk, *_r in STAGE_RUNS[1:])),
             Part('C', '(C) By developmental stage.', c, cols_c,
                  drop={'n_fit': 'checked: equals the testis fit'},
                  sources=stage_sources)]
    return Table(T, title, caption, parts, notes)


# ---------------------------------------------------------------------------
# S6  Rate heterogeneity
# ---------------------------------------------------------------------------

CLOCK_SUBSETS = [('clock_sites12', 'codon positions 1 and 2', '12'),
                 ('clock_sites3', 'codon position 3', '3'),
                 ('clock_sitesall', 'all sites', 'all')]


def read_blocks(path):
    """supp_table_rates.tsv: '# Panel X ...' comment lines start TSV blocks."""
    blocks, cur, name = {}, [], None
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                if name and cur:
                    blocks[name] = cur
                m = re.match(r'#\s*Panel\s+(\w+)', line)
                name, cur = (m.group(1) if m else line.strip('# \n')), []
            elif line.strip():
                cur.append(line.rstrip('\n').split('\t'))
    if name and cur:
        blocks[name] = cur
    out = {}
    for k, rows in blocks.items():
        out[k] = pd.DataFrame(rows[1:], columns=rows[0])
    return out


def outgroup_caption(T, ck, pa, b):
    """Caption text on the outgroup dependence of the Tajima ratios, and a check.

    Each contrast appears with two outgroups. The distance of an outgroup from
    the path joining the two contrasted lineages is taken from the fitted tree
    (Panel A path lengths): (d(l1, og) + d(l2, og) - d(l1, l2)) / 2. Homoplasy on
    the longer outgroup path should pull the ratio towards 1; the check asks
    whether the nearer outgroup gives the larger ratio in every site subset.
    """
    tree_d = {}
    for _, r in pa.iterrows():
        sp = [SPECIES.get(x.strip(), x.strip()) for x in re.split(r'\s+x\s+', str(r['pair']))]
        if len(sp) == 2:
            tree_d[frozenset(sp)] = float(r['tree_predicted'])
    contrasts, ok_all, pos3_all = [], True, True
    for (l1, l2), g in b.groupby(['l1', 'l2'], sort=False):
        ogs = list(dict.fromkeys(g['og']))
        need_d = [frozenset(x) for og in ogs for x in ((l1, og), (l2, og))] + [frozenset((l1, l2))]
        if len(ogs) != 2 or any(k not in tree_d for k in need_d):
            continue
        dist = {og: (tree_d[frozenset((l1, og))] + tree_d[frozenset((l2, og))]
                     - tree_d[frozenset((l1, l2))]) / 2 for og in ogs}
        near, far = sorted(ogs, key=dist.get)
        fold = g.set_index(['og', 'subset'])['fold']
        ratio = g.set_index(['og', 'subset'])['ratio']
        subsets = list(dict.fromkeys(g['subset']))
        ok = all(fold[(near, sb)] > fold[(far, sb)] for sb in subsets)
        lab12, lab3 = CLOCK_SUBSETS[0][1], CLOCK_SUBSETS[1][1]
        att = {sb: fold[(far, sb)] / fold[(near, sb)] for sb in (lab12, lab3)}
        ok_all &= ok
        pos3_all &= att[lab3] < att[lab12]
        r12, r3 = ratio[(near, lab12)], ratio[(near, lab3)]
        contrasts.append({'l1': l1, 'l2': l2, 'near': near, 'far': far,
                          'f_near': fold[(near, 'all sites')], 'f_far': fold[(far, 'all sites')],
                          'r12': r12, 'r3': r3,
                          'stable': abs(r12 - r3) / ((r12 + r3) / 2) < 0.05})
    if not contrasts:
        ck.add(T, 'nearer outgroup gives the larger Tajima ratio in every site subset', False,
               'no contrast with two outgroups and tree distances')
        return 'Ratios depend on the outgroup.'
    ck.add(T, 'nearer outgroup gives the larger Tajima ratio in every site subset', ok_all,
           '; '.join(f"{x['l2']} / {x['l1']}: {x['f_near']:.2f} with {x['near']}, "
                     f"{x['f_far']:.2f} with {x['far']} (all sites)" for x in contrasts))
    vals = ', and '.join(f"{x['l2']} / {x['l1']} is {x['f_near']:.2f}-fold"
                         + (' at all sites' if i == 0 else '')
                         + f" with {x['near']}" + (' as outgroup' if i == 0 else '')
                         + f" and {x['f_far']:.2f}-fold with {x['far']}"
                         for i, x in enumerate(contrasts))
    if ok_all:
        txt = ('Each contrast is shown with both possible outgroups. An outgroup on a longer '
               'path from the two lineages adds homoplasies that pull the ratio towards 1'
               + (', most strongly at codon position 3' if pos3_all else '') + ': '
               + vals + '. The main text quotes the nearer outgroup.')
    else:
        txt = ('Each contrast is shown with both possible outgroups, and the ratio depends on '
               'the outgroup: ' + vals + '. The main text quotes the outgroup nearer to the two '
               'lineages.')
    stable = [x for x in contrasts if x['stable']]
    if stable:
        txt += ' ' + ' '.join(
            f"With it, the {x['l1']} / {x['l2']} ratio is {x['r12']:.2f} at codon positions 1 "
            f"and 2 and {x['r3']:.2f} at position 3, so that asymmetry is not an artefact of "
            'saturation at synonymous sites.' for x in stable)
    return txt


def build_rate_heterogeneity(ctx):
    T = 'rate_heterogeneity'
    ck = ctx.checks
    blocks = read_blocks(ctx.src['rate_table'])
    if 'A' not in blocks or 'B' not in blocks:
        sys.exit(f"ERROR: {ctx.src['rate_table']}: expected '# Panel A' and '# Panel B' "
                 f'blocks, found {sorted(blocks)}')
    pa, pb = blocks['A'], blocks['B']
    need(pa, ['pair', 'median_Ks', 'n_orthologues', 'tree_predicted', 'deviation_pct'],
         ctx.src['rate_table'] + ' (Panel A)')
    need(pb, ['contrast'] + [f'{c}_{s}' for s in ('12', '3', 'all')
                             for c in ('ratio', 'n1', 'n2')], ctx.src['rate_table'] + ' (Panel B)')
    for c in ('median_Ks', 'n_orthologues', 'tree_predicted', 'deviation_pct'):
        pa[c] = pd.to_numeric(pa[c])
    pa['pair_label'] = pa['pair'].map(species_pair)
    dev = 100 * (pa['median_Ks'] - pa['tree_predicted']) / pa['tree_predicted']
    ck.add(T, 'deviation = observed / predicted − 1', close(dev, pa['deviation_pct'],
                                                            atol=1e-3).all())
    off = pa[pa['deviation_pct'].abs() > 1e-3]
    ck.add(T, 'exactly one pair departs from the tree (the pair left out of the exact fit)',
           len(off) == 1, f'{len(off)} pairs')
    with open(ctx.src['rate_caption'], encoding='utf-8') as fh:
        cap_txt = fh.read()
    m_sites = re.search(r'With ([\d,]+) aligned sites', cap_txt)
    m_aln = re.search(r'Based on ([\d,]+) four-way', cap_txt)
    m_model = re.search(r'baseml ([^,;)]+)', cap_txt)
    ck.add(T, 'site and alignment counts found in caption.txt', bool(m_sites and m_aln))

    # --- B: Tajima tests, counts from the rate table, chi-square from the clock runs
    clock = {}
    for key, lab, suf in CLOCK_SUBSETS:
        c = ctx.tsv(key)
        need(c, ['test', 'comparison', 'statistic', 'value', 'df', 'p_value'], ctx.src[key])
        clock[suf] = (lab, c)
    rows_b, worst_chi, worst_ratio = [], 0.0, 0.0
    for _, r in pb.iterrows():
        mm = re.match(r'^(\w+) vs (\w+) \(out: (\w+)\)$', str(r['contrast']).strip())
        if not mm:
            sys.exit(f"ERROR: cannot parse contrast {r['contrast']!r}")
        l1, l2, og = mm.groups()
        for suf in ('12', '3', 'all'):
            lab, c = clock[suf]
            n1, n2 = float(r[f'n1_{suf}']), float(r[f'n2_{suf}'])
            hit = c[(c['test'] == 'tajima_RRT')
                    & (c['comparison'].astype(str).str.strip() == str(r['contrast']).strip())]
            if len(hit) != 1:
                sys.exit(f"ERROR: {r['contrast']!r} not found once in the {lab} clock results")
            h = hit.iloc[0]
            chi = (n1 - n2) ** 2 / (n1 + n2)
            worst_chi = max(worst_chi, max_rel(chi, float(h['value'])))
            worst_ratio = max(worst_ratio, abs(float(r[f'ratio_{suf}']) - n1 / n2))
            rows_b.append({'l1': SPECIES.get(l1, l1), 'l2': SPECIES.get(l2, l2),
                           'og': SPECIES.get(og, og), 'subset': lab, 'n1': n1, 'n2': n2,
                           'ratio': n1 / n2, 'fold': n2 / n1, 'chi2': float(h['value']),
                           'df': int(float(h['df'])), 'p': float(h['p_value'])})
    ck.add(T, 'each Tajima χ² equals (n1 − n2)² / (n1 + n2)', worst_chi < 1e-6,
           f'largest rel. diff {worst_chi:.1e}')
    ck.add(T, 'rounded ratios in the rate table agree with the counts', worst_ratio <= 0.0051,
           f'largest abs. diff {worst_ratio:.4f}')
    rows_c = []
    for _, lab, suf in CLOCK_SUBSETS:
        lab, c = clock[suf]
        g = c[c['test'] == 'global_clock_LRT']
        if len(g) != 1:
            sys.exit(f'ERROR: no single global_clock_LRT row in the {lab} clock results')
        g = g.iloc[0]
        rows_c.append({'subset': lab, 'lrt': float(g['value']), 'df': int(float(g['df'])),
                       'p': float(g['p_value'])})
    b, c = pd.DataFrame(rows_b), pd.DataFrame(rows_c)
    o = off.iloc[0] if len(off) == 1 else None
    og_txt = outgroup_caption(T, ck, pa, b)

    cols_a = [
        Col('species_pair', 'pair_label', 'Species pair', 'text'),
        Col('median_dS', 'median_Ks', 'Median pairwise dS (yn00) over the orthologue pairs',
            'f4'),
        Col('n_orthologue_pairs', 'n_orthologues', 'Orthologue pairs after filtering', 'int'),
        Col('tree_predicted_dS', 'tree_predicted', 'Path length between the two species on '
            'the four-taxon tree', 'f4'),
        Col('deviation_pct', 'deviation_pct', '100 × (observed / predicted − 1)', 'f2'),
    ]
    cols_b = [
        Col('lineage_1', 'l1', 'First lineage', 'text'),
        Col('lineage_2', 'l2', 'Second lineage', 'text'),
        Col('outgroup', 'og', 'Outgroup', 'text'),
        Col('sites', 'subset', 'Alignment sites used', 'text', values=True),
        Col('substitutions_lineage_1', 'n1', 'Sites at which lineage 1 alone differs from '
            'the other two sequences', 'int'),
        Col('substitutions_lineage_2', 'n2', 'Sites at which lineage 2 alone differs', 'int'),
        Col('ratio_1_to_2', 'ratio', 'Lineage 1 / lineage 2 (1 for equal rates)', 'f3'),
        Col('fold_2_over_1', 'fold', 'Lineage 2 / lineage 1', 'f3'),
        Col('chi2', 'chi2', "Tajima's χ², (n1 − n2)² / (n1 + n2)", 'f1'),
        Col('df', 'df', 'Degrees of freedom', 'int'),
        Col('P', 'p', 'P from χ²', 'p'),
    ]
    cols_c = [
        Col('sites', 'subset', 'Alignment sites used', 'text', values=True),
        Col('two_delta_lnL', 'lrt', '2ΔlnL, unconstrained versus global clock', 'f1'),
        Col('df', 'df', 'Degrees of freedom', 'int'),
        Col('P', 'p', 'P from χ²', 'p'),
    ]
    n_sites = m_sites.group(1) if m_sites else '[VERIFY]'
    n_aln = m_aln.group(1) if m_aln else '[VERIFY]'
    model = m_model.group(1).strip() if m_model else 'GTR+G'
    title = 'Substitution-rate heterogeneity among the four lineages.'
    caption = (
        '(A) Median pairwise dS (yn00) for each species pair, against the path length '
        'between the two species on the four-taxon tree. The five branch lengths are solved '
        'exactly from five of the six pairs, which the tree therefore reproduces by '
        'construction; the sixth, '
        + (f"{o['pair_label']}, is the test of additivity and lies "
           f"{abs(o['deviation_pct']):.1f}% {'below' if o['deviation_pct'] < 0 else 'above'} "
           'prediction, as expected for the longest and most saturated path. '
           if o is not None else 'left out of the fit, is the test of additivity. ')
        + f"(B) Tajima's relative-rate test on {n_aln} four-species codon alignments, at "
        'codon positions 1 and 2, position 3 and all sites: substitutions inferred on each '
        'of two lineages against an outgroup, their ratio (1 for equal rates), and χ² '
        '(1 df). No correction for multiple hits is applied, which compresses the ratio '
        'towards 1 for the faster lineage. ' + og_txt + ' (C) '
        f'Likelihood-ratio test of a global molecular clock (baseml, {model}; 2 df) for '
        f'each site subset. With {n_sites} aligned sites the test has essentially unlimited '
        'power, so the size of the rate differences in B, not the P value, is the '
        f'informative result. P values that underflow are shown as {P_UNDERFLOW}.')
    notes = [f"{r['l2']} / {r['l1']} (out {r['og']}), {r['subset']}: "
             f"{r['fold']:.2f}-fold" for _, r in b.iterrows() if r['subset'] == 'all sites']
    parts = [Part('A', '(A) Pairwise dS against the four-taxon tree.', pa, cols_a,
                  drop={'pair': 'shown with species names', 'source_column': 'constant'},
                  sources=('rate_table',)),
             Part('B', "(B) Tajima's relative-rate tests.", b, cols_b,
                  sources=('rate_table',) + tuple(k for k, _, _ in CLOCK_SUBSETS)),
             Part('C', '(C) Global molecular clock.', c, cols_c,
                  sources=tuple(k for k, _, _ in CLOCK_SUBSETS))]
    return Table(T, title, caption, parts, notes)


# ---------------------------------------------------------------------------
# S7  T. canis gene-model diagnostic
# ---------------------------------------------------------------------------

TOX_QUANTITY = {'orf_len': 'T. canis ORF length (amino acids)',
                'align_len': 'Quartet codon alignment length (codons)',
                'orf_complete_rate': 'Complete T. canis ORFs (fraction)'}
TOX_TEST = {'mannwhitney': 'Mann–Whitney U', 'fisher': "Fisher's exact"}
TOX_STRATUM = {'all': 'all quartets', 'complete_orf': 'complete T. canis ORF'}
TOX_COVARIATE = {'orf_len': 'ORF length', 'align_len': 'alignment length',
                 'orf_5prime_partial': "5′-partial ORF", 'orf_3prime_partial': "3′-partial ORF",
                 'orf_internal': 'internal ORF', '-': 'none'}


def build_toxocara_models(ctx):
    T = 'toxocara_models'
    ck = ctx.checks
    t1, t2, t3 = ctx.tsv('tox_premise'), ctx.tsv('tox_strata'), ctx.tsv('tox_regression')
    need(t1, ['grouping', 'metric', 'test', 'ret_n', 'elim_n', 'ret_value', 'elim_value',
              'ratio', 'p_value'], ctx.src['tox_premise'])
    need(t2, ['grouping', 'metric', 'capped', 'branch', 'stratum', 'ret_n', 'elim_n',
              'ret_median', 'elim_median', 'fold_change', 'p_value', 'q_value'],
         ctx.src['tox_strata'])
    need(t3, ['grouping', 'branch', 'model', 'n', 'r2', 'pseudocount', 'covariates',
              'beta_eliminated', 'se_eliminated', 't_eliminated', 'p_eliminated',
              'fold_from_beta', 'pct_attenuation'], ctx.src['tox_regression'])

    # --- A ---------------------------------------------------------------------------
    t1['labelling'] = t1['grouping'].map(LABELLING)
    t1['quantity'] = t1['metric'].map(TOX_QUANTITY).fillna(t1['metric'])
    t1['test_label'] = t1['test'].map(TOX_TEST).fillna(t1['test'])
    drop_a = {'grouping': 'shown as labelling', 'metric': 'shown as quantity',
              'test': 'shown with its full name'}
    cols_a = [
        Col('labelling', 'labelling', 'Elimination annotation used to label quartets',
            'text', values=True),
        Col('quantity', 'quantity', 'Gene-model property compared', 'text', values=True),
        Col('test', 'test_label', 'Two-sided test', 'text', values=True),
        Col('n_retained', 'ret_n', 'Retained quartets', 'int'),
        Col('n_eliminated', 'elim_n', 'Eliminated quartets', 'int'),
        Col('value_retained', 'ret_value', 'Median (lengths) or fraction (complete ORFs), '
            'retained', 'v4'),
        Col('value_eliminated', 'elim_value', 'The same, eliminated', 'v4'),
        Col('ratio', 'ratio', 'Eliminated / retained', 'f3'),
        Col('P', 'p_value', 'Two-sided P', 'p'),
    ]
    if 'divergence_confounded' in t1.columns:
        dc = to_bool(t1['divergence_confounded'])
        if not dc.fillna(False).any():
            drop_a['divergence_confounded'] = 'False in every row (no divergence-derived ' \
                                              'covariate among these properties)'
        else:
            cols_a.append(Col('divergence_confounded', 'divergence_confounded', 'Property '
                              'derived from sequence divergence, so not a clean measure of '
                              'model quality', 'bool'))

    # --- B ---------------------------------------------------------------------------
    t2['capped'] = to_bool(t2['capped'])
    t2['labelling'] = t2['grouping'].map(LABELLING)
    t2['branch_label'] = t2['branch'].map(SPECIES)
    t2['metric_label'] = t2['metric'].map(METRIC)
    t2['mask_label'] = np.where(t2['capped'].fillna(False), 'applied', 'none')
    t2['stratum_label'] = t2['stratum'].map(TOX_STRATUM).fillna(t2['stratum'])
    bad = t2[['labelling', 'branch_label', 'metric_label']].isna().any(axis=1)
    if bad.any():
        sys.exit(f"ERROR: {ctx.src['tox_strata']}: unexpected grouping, branch or metric")
    exp_mask = ((t2['metric'] != 'dS') == t2['capped'].fillna(False)).all()
    ck.add(T, f"{NUMBER[T]}B masks dN/dS and dN and leaves dS unmasked, as in "
           f"{NUMBER['branch_models']}", bool(exp_mask))
    drop_b = {'grouping': 'shown as labelling', 'branch': 'shown as branch',
              'metric': 'shown as metric', 'capped': 'shown as masking',
              'stratum': 'shown as stratum'}
    cols_b = [
        Col('labelling', 'labelling', 'Elimination annotation used to label quartets',
            'text', values=True),
        Col('metric', 'metric_label', 'Branch estimate compared', 'text', values=True),
        Col('masking', 'mask_label', f"Masking as in supplementary table "
            f"{NUMBER['branch_models']}", 'text', values=True),
        Col('branch', 'branch_label', 'Terminal branch', 'text'),
        Col('stratum', 'stratum_label', 'Quartets compared', 'text', values=True),
        Col('n_retained', 'ret_n', 'Retained quartets', 'int'),
        Col('n_eliminated', 'elim_n', 'Eliminated quartets', 'int'),
        Col('median_retained', 'ret_median', 'Median, retained', 'f4'),
        Col('median_eliminated', 'elim_median', 'Median, eliminated', 'f4'),
        Col('fold_change', 'fold_change', 'Eliminated median / retained median', 'f3'),
        Col('P', 'p_value', 'Two-sided Mann–Whitney U test', 'p'),
        Col('q', 'q_value', 'Benjamini–Hochberg q within each labelling and metric, over '
            'branches and strata', 'p'),
    ]
    if 'underpowered' in t2.columns:
        up = to_bool(t2['underpowered'])
        if not up.fillna(False).any():
            drop_b['underpowered'] = 'False in every row (every stratum has at least the ' \
                                     'minimum number of eliminated quartets)'
        else:
            cols_b.append(Col('underpowered', 'underpowered', 'Fewer eliminated quartets '
                              'than the minimum for a stratum', 'bool'))
    # cross-check: the 'all quartets' rows are the S2A rows
    if 'branch_decomposition' in ctx.src and os.path.isfile(ctx.src['branch_decomposition']):
        dec = ctx.tsv('branch_decomposition')
        dec['masked'] = to_bool(dec['masked']).fillna(False).astype(bool)
        allq = t2[t2['stratum'] == 'all'].copy()
        allq['capped'] = allq['capped'].fillna(False).astype(bool)
        al = allq.merge(
            dec, left_on=['grouping', 'branch', 'metric', 'capped'],
            right_on=['grouping', 'branch', 'metric', 'masked'], suffixes=('', '_s2'),
            how='left')
        ok = (al['ret_n_s2'].notna().all()
              and close(al['ret_n'], al['ret_n_s2']).all()
              and close(al['elim_n'], al['elim_n_s2']).all()
              and close(al['ret_median'], al['ret_median_s2']).all()
              and close(al['elim_median'], al['elim_median_s2']).all()
              and close(al['p_value'], al['p_value_s2'], rtol=PRINTED_RTOL).all())
        ck.add(T, f"{NUMBER[T]}B 'all quartets' rows equal the matching "
               f"{NUMBER['branch_models']}A rows", ok,
               f'{len(al)} rows')

    # --- C ---------------------------------------------------------------------------
    t3['labelling'] = t3['grouping'].map(LABELLING)
    t3['branch_label'] = t3['branch'].map(SPECIES)
    t3['model_label'] = t3['model'].map({'unadjusted': 'unadjusted',
                                         'covariate_adjusted': 'adjusted'}).fillna(t3['model'])
    t3['covariate_label'] = t3['covariates'].astype(str).map(
        lambda s: ', '.join(TOX_COVARIATE.get(x.strip(), x.strip()) for x in s.split(',')))
    ck.add(T, 'fold = exp(β)', close(t3['fold_from_beta'], np.exp(t3['beta_eliminated']),
                                     rtol=3 * PRINTED_RTOL).all())
    piv = t3.pivot_table(index=['grouping', 'branch'], columns='model',
                         values='beta_eliminated')
    if {'unadjusted', 'covariate_adjusted'} <= set(piv.columns):
        att = (1 - piv['covariate_adjusted'] / piv['unadjusted']) * 100
        chk = t3.merge(att.rename('att').reset_index(), on=['grouping', 'branch'])
        # β is printed to six digits, so a small difference of two β carries an
        # absolute error of a few 1e-4 percentage points; 0.005 tests the definition.
        ck.add(T, 'attenuation = 100 × (1 − β adjusted / β unadjusted)',
               close(chk['pct_attenuation'], chk['att'], rtol=0, atol=0.005).all())
    t3['attenuation'] = np.where(t3['model'] == 'covariate_adjusted',
                                 t3['pct_attenuation'], np.nan)
    drop_c = {'grouping': 'shown as labelling', 'branch': 'shown as branch',
              'model': 'shown as model', 'covariates': 'shown with readable names',
              'pct_attenuation': 'shown on the adjusted row only'}
    cols_c = [
        Col('labelling', 'labelling', 'Elimination annotation used to label quartets',
            'text', values=True),
        Col('branch', 'branch_label', 'Branch whose dN is the response', 'text'),
        Col('model', 'model_label', 'Without (unadjusted) or with (adjusted) the T. canis '
            'gene-model covariates', 'text', values=True),
        Col('covariates', 'covariate_label', 'Covariates; ORF type enters as indicators '
            'against complete ORFs', 'text'),
        Col('n', 'n', 'Quartets in the fit', 'int'),
        Col('R2', 'r2', 'Coefficient of determination', 'f4'),
        Col('pseudocount', 'pseudocount', 'c in ln(dN + c): half the smallest positive dN on '
            'that branch', 'g'),
        Col('beta_eliminated', 'beta_eliminated', 'Elimination coefficient on ln(dN + c)',
            'f4'),
        Col('SE_HC3', 'se_eliminated', 'Heteroscedasticity-consistent (HC3) standard error',
            'f4'),
        Col('t', 't_eliminated', 't statistic', 'f3'),
        Col('P', 'p_eliminated', 'Two-sided P', 'p'),
        Col('fold', 'fold_from_beta', 'exp(β): eliminated / retained dN', 'f3'),
        Col('attenuation_pct', 'attenuation', 'Adjusted rows: percentage reduction of β on '
            'adding the covariates', 'f2'),
    ]
    title = 'Assessment of T. canis gene-model quality as a source of the elevated dN of ' \
            'eliminated genes.'
    caption = (
        'Quartets of supplementary table ' + NUMBER['branch_models'] + ', labelled by the '
        'Ascaris or by the T. canis elimination annotation. (A) Gene-model properties of '
        'eliminated and retained quartets: TransDecoder ORF length and codon alignment '
        "length (Mann–Whitney U tests) and the fraction of complete ORFs (Fisher's exact "
        'test). (B) Branch comparisons repeated on quartets whose T. canis ORF is complete, '
        'beside all quartets, with dN/dS and dN compared with masked branches excluded and '
        'dS on every branch. (C) Ordinary least squares of ln(dN + c) on elimination status, '
        'without and with the T. canis gene-model covariates, with HC3 standard errors. The '
        'same T. canis covariates are applied to all four branches as a control: '
        'attenuation that appears on every branch reflects covariates tracking gene class, '
        'not annotation noise.')
    notes = []
    for _, r in t1.iterrows():
        notes.append(f"{r['labelling']:<8} {r['quantity']}: {r['elim_value']:.4g} vs "
                     f"{r['ret_value']:.4g}, P {p_words(r['p_value'])}")
    parts = [Part('A', '(A) Gene-model properties.', t1, cols_a, drop=drop_a,
                  sources=('tox_premise',)),
             Part('B', '(B) Branch comparisons on complete T. canis ORFs.', t2, cols_b,
                  drop=drop_b, sources=('tox_strata',)),
             Part('C', '(C) Regression of dN on elimination status with gene-model '
                  'covariates.', t3, cols_c, drop=drop_c, sources=('tox_regression',))]
    return Table(T, title, caption, parts, notes)


# ---------------------------------------------------------------------------
# S8  Micro-synteny
# ---------------------------------------------------------------------------

SYN_CLASS = {'retained': 'retained', 'boundary': 'CBR-spanning', 'eliminated': 'eliminated'}
SYN_ORDER = ['retained', 'boundary', 'eliminated']


def build_microsynteny(ctx):
    T = 'microsynteny'
    ck = ctx.checks
    s, t, w = ctx.tsv('synteny_summary'), ctx.tsv('synteny_tests'), ctx.tsv('synteny_sweep')
    need(s, ['pair_class', 'max_distance', 'n_pairs', 'n_conserved', 'pct_conserved',
             'pct_lo', 'pct_hi', 'n_same_seq_beyond', 'n_different_seq',
             'pct_different_seq'], ctx.src['synteny_summary'])
    need(t, ['max_distance', 'group1', 'group2', 'conserved1', 'n1', 'conserved2', 'n2',
             'odds_ratio', 'or_lo', 'or_hi', 'p', 'q'], ctx.src['synteny_tests'])
    need(w, ['max_distance', 'n_pairs', 'n_conserved', 'pct_conserved', 'pair_class',
             'odds_ratio', 'or_lo', 'or_hi', 'p'], ctx.src['synteny_sweep'])
    d = sorted(s['max_distance'].unique())
    if len(d) != 1:
        sys.exit(f"ERROR: {ctx.src['synteny_summary']}: one threshold expected, got {d}")
    d = int(d[0])
    ck.add(T, 'tests use the summary threshold', set(t['max_distance']) == {d})
    parts_sum = s['n_conserved'] + s['n_same_seq_beyond'] + s['n_different_seq']
    ck.add(T, 'conserved + same sequence beyond threshold + different sequence = all pairs',
           bool((parts_sum == s['n_pairs']).all()))
    s = s.assign(_o=s['pair_class'].map(SYN_ORDER.index)).sort_values('_o').reset_index(drop=True)
    s['class_label'] = s['pair_class'].map(SYN_CLASS)
    w_d = w[w['max_distance'] == d].set_index('pair_class')
    same = all(int(w_d.loc[c, 'n_conserved']) == int(r['n_conserved'])
               and int(w_d.loc[c, 'n_pairs']) == int(r['n_pairs'])
               for c, r in s.set_index('pair_class').iterrows())
    ck.add(T, f'sweep at {d} genes equals the summary', same)
    t['comparison'] = t['group1'].map(SYN_CLASS) + ' vs ' + t['group2'].map(SYN_CLASS)
    re_ = t[(t['group1'] == 'retained') & (t['group2'] == 'eliminated')]
    ck.add(T, f'sweep test at {d} genes equals the retained-versus-eliminated test',
           len(re_) == 1
           and close(re_['p'].values, [w_d.loc['retained', 'p']], rtol=PRINTED_RTOL).all()
           and close(re_['odds_ratio'].values, [w_d.loc['retained', 'odds_ratio']],
                     rtol=PRINTED_RTOL).all())
    if 'or_method' in t.columns:
        ck.add(T, 'odds ratios are conditional MLE (scipy)',
               set(t['or_method'].astype(str)) == {'scipy'})
    # C: one row per threshold
    rows = []
    for dist, g in w.groupby('max_distance', sort=True):
        g = g.set_index('pair_class')
        row = {'d': int(dist)}
        for c in SYN_ORDER:
            if c in g.index:
                row[f'{c}_k'] = int(g.loc[c, 'n_conserved'])
                row[f'{c}_n'] = int(g.loc[c, 'n_pairs'])
                row[f'{c}_pct'] = float(g.loc[c, 'pct_conserved'])
        for k in ('odds_ratio', 'or_lo', 'or_hi', 'p'):
            row[k] = float(g.loc['retained', k]) if 'retained' in g.index else np.nan
        rows.append(row)
    c = pd.DataFrame(rows)

    cols_a = [
        Col('pair_class', 'class_label', 'Regions of the two adjacent Ascaris genes: both '
            'retained, both eliminated, or one of each (CBR-spanning)', 'text', values=True),
        Col('n_pairs', 'n_pairs', 'Adjacent pairs whose two genes both have a '
            'P. univalens orthologue', 'int'),
        Col('n_conserved', 'n_conserved', f'Orthologues on one P. univalens sequence within '
            f'{d} genes', 'int'),
        Col('pct_conserved', 'pct_conserved', 'Percentage conserved', 'f1'),
        Col('pct_CI_low', 'pct_lo', 'Lower Wilson 95% bound', 'f1'),
        Col('pct_CI_high', 'pct_hi', 'Upper Wilson 95% bound', 'f1'),
        Col('n_same_sequence_beyond', 'n_same_seq_beyond', f'Orthologues on one sequence '
            f'but more than {d} genes apart', 'int'),
        Col('n_different_sequences', 'n_different_seq', 'Orthologues on different '
            'P. univalens sequences', 'int'),
        Col('pct_different_sequences', 'pct_different_seq', 'Percentage on different '
            'sequences', 'f2'),
    ]
    cols_b = [
        Col('comparison', 'comparison', 'Pair classes compared', 'text'),
        Col('conserved_1', 'conserved1', 'Conserved pairs, first class', 'int'),
        Col('n_1', 'n1', 'Pairs, first class', 'int'),
        Col('conserved_2', 'conserved2', 'Conserved pairs, second class', 'int'),
        Col('n_2', 'n2', 'Pairs, second class', 'int'),
        Col('odds_ratio', 'odds_ratio', 'Conditional maximum-likelihood odds ratio of '
            'conservation, first relative to second class', 'f2'),
        Col('odds_ratio_CI_low', 'or_lo', 'Lower exact 95% bound', 'f2'),
        Col('odds_ratio_CI_high', 'or_hi', 'Upper exact 95% bound', 'f2'),
        Col('P', 'p', "Two-sided Fisher's exact test", 'p'),
        Col('q', 'q', 'Benjamini–Hochberg q over the three comparisons', 'p'),
    ]
    cols_c = [Col('max_rank_distance', 'd', 'Largest P. univalens rank distance counted '
                  'as conserved', 'int')]
    for cl in SYN_ORDER:
        lab = SYN_CLASS[cl].replace('-', '_')
        cols_c += [Col(f'{lab}_conserved', f'{cl}_k', f'Conserved {SYN_CLASS[cl]} pairs',
                       'int'),
                   Col(f'{lab}_n', f'{cl}_n', f'{SYN_CLASS[cl].capitalize()} pairs', 'int'),
                   Col(f'{lab}_pct', f'{cl}_pct', 'Percentage conserved', 'f1')]
    cols_c += [Col('odds_ratio_retained_vs_eliminated', 'odds_ratio', 'Conditional '
                   'maximum-likelihood odds ratio', 'f2'),
               Col('odds_ratio_CI_low', 'or_lo', 'Lower exact 95% bound', 'f2'),
               Col('odds_ratio_CI_high', 'or_hi', 'Upper exact 95% bound', 'f2'),
               Col('P', 'p', "Two-sided Fisher's exact test", 'p')]
    title = 'Conservation of local gene order between Ascaris and P. univalens.'
    caption = (
        'Adjacent Ascaris gene pairs in which both genes have a P. univalens orthologue, '
        'classed by the regions of the two genes. A pair is conserved when its two '
        f'orthologues lie on one P. univalens sequence within {d} genes of each other, with '
        'every gene model ranked by its midpoint. (A) Counts at the primary threshold of '
        f'{d} genes, with Wilson 95% intervals; pairs that are not conserved are split into '
        'those on one sequence beyond the threshold and those on different sequences. (B) '
        "Two-sided Fisher's exact tests, with conditional maximum-likelihood odds ratios, "
        'exact 95% intervals and Benjamini–Hochberg q over the three comparisons. (C) The '
        'comparison at rank-distance thresholds from '
        f"{int(c['d'].min())} to {int(c['d'].max())} genes.")
    r0 = s.set_index('pair_class')
    notes = [f"{SYN_CLASS[k]}: {int(r0.loc[k, 'n_conserved']):,}/{int(r0.loc[k, 'n_pairs']):,}"
             f" ({r0.loc[k, 'pct_conserved']:.1f}%)" for k in SYN_ORDER if k in r0.index]
    notes += [f"{r['comparison']}: OR {fmt_num(r['odds_ratio'], 1)} "
              f"[{fmt_num(r['or_lo'], 1)}, {fmt_num(r['or_hi'], 1)}], P {p_words(r['p'])}"
              for _, r in t.iterrows()]
    parts = [Part('A', f'(A) Conserved adjacent pairs at {d} genes.', s, cols_a,
                  drop={'pair_class': 'shown with readable names', 'max_distance':
                        f'constant ({d})', '_o': 'sort key'}, sources=('synteny_summary',)),
             Part('B', '(B) Tests.', t, cols_b,
                  drop={'max_distance': f'constant ({d})', 'group1': 'in comparison',
                        'group2': 'in comparison', 'odds_ratio_sample': 'the sample odds '
                        'ratio; the conditional MLE is reported', 'or_method': 'constant'},
                  sources=('synteny_tests',)),
             Part('C', '(C) Threshold sweep.', c, cols_c, sources=('synteny_sweep',))]
    return Table(T, title, caption, parts, notes)


# ---------------------------------------------------------------------------
# S9  Population-sample screen
# ---------------------------------------------------------------------------

SEX_LABEL = {'M': 'male', 'F': 'female', 'unknown': 'unassigned'}
REASON_LABEL = {'somatic': 'germline fraction below threshold',
                'low_mapped_depth': 'too few reads mapped to assess',
                'control': 'tissue control'}


def build_popgen_cohort(ctx):
    T = 'popgen_cohort'
    ck = ctx.checks
    thr = ctx.args.cohort_threshold
    md = ctx.tsv('cohort_metadata')
    need(md, ['run', 'sex', 'x_to_autosome', 'germline_fraction', 'germline_rescaled',
              'depth_retained_median', 'depth_eliminated_median', 'pass', 'fail_reason',
              'is_control'], ctx.src['cohort_metadata'])
    md['pass'] = to_bool(md['pass']).fillna(False).astype(bool)
    md['is_control'] = to_bool(md['is_control']).fillna(False).astype(bool)
    ck.add(T, 'one row per run', not md['run'].duplicated().any())
    worms, ctrl = md[~md['is_control']], md[md['is_control']]
    sel = md[md['pass']]
    ck.add(T, f'every selected library has germline fraction ≥ {thr:g}',
           bool((sel['germline_fraction'] >= thr).all()))
    som = worms[worms['fail_reason'].astype(str) == 'somatic']
    ck.add(T, f"every 'somatic' library is below {thr:g}",
           bool((som['germline_fraction'] < thr).all()))
    ck.add(T, 'no control is selected', not ctrl['pass'].any())
    resc = np.clip((md['germline_fraction'] - 0.03) / 1.09, 0, 1)
    ck.add(T, 'rescaled fraction = (fraction − 0.03) / 1.09, clipped to [0, 1]',
           close(resc, md['germline_rescaled'], rtol=1e-6, atol=1e-9).all())
    if 'cohort_audit' in ctx.src and os.path.isfile(ctx.src['cohort_audit']):
        au = ctx.tsv('cohort_audit')
        need(au, ['sample', 'depth_retained_median', 'depth_eliminated_median',
                  'germline_fraction', 'x_to_autosome'], ctx.src['cohort_audit'])
        au['run'] = au['sample'].astype(str).str.replace(r'\.sub\d+$', '', regex=True)
        j = md.merge(au, on='run', how='left', suffixes=('', '_audit'))
        found = j['sample'].notna()
        ok = (found.all()
              and all(close(j.loc[found, c], j.loc[found, f'{c}_audit']).all()
                      for c in ('depth_retained_median', 'depth_eliminated_median',
                                'germline_fraction', 'x_to_autosome')))
        ck.add(T, 'screen values identical to the audit table', ok,
               f'{int(found.sum())} of {len(md)} runs found')
        subs = sorted(set(au['sample'].astype(str).str.extract(r'\.sub(\d+)$')[0].dropna()))
        sub_txt = (f'subsampled to {int(subs[0]) / 1e6:g} million read pairs'
                   if len(subs) == 1 else 'subsampled to a common number of read pairs')
    else:
        sub_txt = 'subsampled to a common number of read pairs'
        logging.info('  cohort_audit not given: values not cross-checked')
    md['sex_label'] = md['sex'].map(SEX_LABEL).fillna(md['sex'])
    md['reason_label'] = np.where(md['pass'], '',
                                  md['fail_reason'].map(REASON_LABEL).fillna(md['fail_reason']))
    ranges = []
    for code in ('F', 'M', 'unknown'):
        v = worms.loc[worms['sex'] == code, 'x_to_autosome']
        if len(v):
            ranges.append(f'{SEX_LABEL[code]} {v.min():.2f}–{v.max():.2f} (n = {len(v)})')
    sel_sex = ', '.join(f'{n} {SEX_LABEL.get(k, k)}' for k, n in
                        sel['sex'].value_counts().reindex(['M', 'F', 'unknown']).dropna()
                        .astype(int).items())
    reasons = worms.loc[~worms['pass'], 'fail_reason'].value_counts().to_dict()

    cols = [
        Col('run', 'run', 'Sequence Read Archive run', 'text'),
        Col('control', 'is_control', 'Tissue control (BioProject PRJNA62057) rather than '
            'a population-sample worm', 'bool'),
        Col('sex', 'sex_label', 'Sex called from the X-to-autosome depth ratio', 'text',
            values=True),
        Col('X_to_autosome_depth_ratio', 'x_to_autosome', 'Median depth on the X '
            'chromosomes / median depth on the autosomes (retained windows)', 'f3'),
        Col('germline_fraction', 'germline_fraction', 'Median depth in eliminated windows / '
            'median depth in retained autosomal windows', 'f3'),
        Col('germline_fraction_rescaled', 'germline_rescaled', '(germline fraction − 0.03) '
            '/ 1.09, clipped to 0–1: the somatic floor at 0 and the germline control at 1',
            'f3'),
        Col('median_depth_retained', 'depth_retained_median', 'Median depth, retained '
            'autosomal windows, in the subsampled screen', 'f2'),
        Col('median_depth_eliminated', 'depth_eliminated_median', 'Median depth, eliminated '
            'windows, in the subsampled screen', 'f2'),
        Col('selected', 'pass', 'Carried forward to genotype calling', 'bool'),
        Col('reason_not_selected', 'reason_label', 'Why a library was not carried forward',
            'text', values=True),
    ]
    title = 'Population-sample libraries screened for germline DNA.'
    caption = (
        f'{len(worms)} whole-genome libraries of adult worms from an Ethiopian population '
        "(Landeryou et al. 2025; BioProject PRJNA1195031) and "
        f'{len(ctrl)} tissue controls (BioProject PRJNA62057: SRR502751, male spermatids, '
        'germline; SRR502756, female carcass; SRR5573034, male intestine), each '
        f'{sub_txt} and mapped to the Ascaris germline assembly. Depths are medians over '
        '10 kb windows at mapping quality ≥ 20. Germline fraction, median depth in '
        'eliminated windows divided by median depth in retained autosomal windows; '
        'rescaled, (fraction − 0.03)/1.09 clipped to 0–1, which places the somatic floor at 0 '
        'and the germline control at 1. Sex is called from the X-to-autosome '
        'depth ratio; observed ranges among worms: ' + '; '.join(ranges) + '. Libraries '
        f'with a germline fraction of at least {thr:.2f} were carried forward ({len(sel)}: '
        f'{sel_sex}).')
    notes = [f'{len(worms)} worms, {len(ctrl)} controls; {len(sel)} selected ({sel_sex}); '
             'not selected: ' + ', '.join(f'{v} {k}' for k, v in reasons.items()),
             'sex calls among worms: ' + '; '.join(ranges)]
    part = Part('', '', md, cols,
                drop={'sex': 'shown as a word', 'fail_reason': 'shown as a phrase'},
                freeze_cols=1,
                sources=('cohort_metadata', 'cohort_audit'))
    return Table(T, title, caption, [part], notes)


# ---------------------------------------------------------------------------
# S11 Alternatives tested (curated, plus the McDonald-Kreitman row)
# ---------------------------------------------------------------------------

#: McDonald-Kreitman variant sets shown in S11, as pi_by_region.py labels them
#: (f'MAF>={t:g}' for its --mk-maf floors, default 0.10,0.15,0.20), and as the
#: table prints them.
MK_FLOORS = {'MAF>=0.1': 'MAF ≥ 0.10', 'MAF>=0.15': 'MAF ≥ 0.15',
             'MAF>=0.2': 'MAF ≥ 0.20'}
#: What the Results and fig. 6D quote for the first floor (DoS to three
#: decimals, alpha to two), and the genes the Results give.
MK_TEXT = {'retained': {'DoS': 0.102, 'alpha': 0.35, 'genes': 6882},
           'eliminated': {'DoS': 0.064, 'alpha': 0.23, 'genes': 264}}
#: The range of eliminated alpha across the floors, as the Results give it
#: (two decimals).
MK_ALPHA_TEXT = (0.18, 0.23)


def signed(v, nd=3):
    """Fixed decimals with an explicit sign and a true minus."""
    return f'{v:+.{nd}f}'.replace('-', '\u2212')


def unsigned(v, nd=3):
    """Fixed decimals, true minus for negatives."""
    return f'{v:.{nd}f}'.replace('-', '\u2212')


def mk_alternative_row(mk, ck, T, src_label):
    """The positive-selection row of S11 from pi_by_region.py's mk_summary.tsv.

    Positive selection on eliminated genes, rather than relaxed constraint,
    would raise their Dn relative to Pn, so their DoS and alpha would exceed
    those of retained genes. Pooled DoS (with its gene-bootstrap interval) and
    alpha are shown for each class at every MAF floor. Returns the row (a dict
    of the curated columns) and lines for the summary.
    """
    need(mk, ['variant_set', 'region', 'genes', 'Pn', 'Ps', 'Dn', 'Ds', 'alpha',
              'DoS', 'DoS_lo', 'DoS_hi'], src_label)
    mk = mk.copy()
    mk['variant_set'] = mk['variant_set'].astype(str).str.strip()
    mk['region'] = mk['region'].astype(str).str.strip()
    have = sorted(set(mk['variant_set']))
    missing = [f for f in MK_FLOORS if f not in have]
    if missing:
        sys.exit(f'ERROR: {src_label}: no McDonald-Kreitman rows for {missing} '
                 f'(present: {have}). Rerun pi_by_region.py with its default '
                 f'--mk-maf 0.10,0.15,0.20.')
    rows = {}
    for f in MK_FLOORS:
        for reg in ('retained', 'eliminated'):
            r = mk[(mk['variant_set'] == f) & (mk['region'] == reg)]
            if len(r) != 1:
                sys.exit(f'ERROR: {src_label}: {len(r)} rows for {f} / {reg}, '
                         f'expected 1')
            rows[(f, reg)] = {c: float(r.iloc[0][c]) for c in
                              ('genes', 'Pn', 'Ps', 'Dn', 'Ds', 'alpha', 'DoS',
                               'DoS_lo', 'DoS_hi')}
    regions = ('retained', 'eliminated')

    genes = {reg: sorted({int(rows[(f, reg)]['genes']) for f in MK_FLOORS})
             for reg in regions}
    ck.add(T, 'McDonald-Kreitman: the same genes at every MAF floor',
           all(len(v) == 1 for v in genes.values()),
           '; '.join(f'{reg} {v}' for reg, v in genes.items()))
    f0 = next(iter(MK_FLOORS))
    got = {reg: (round(rows[(f0, reg)]['DoS'], 3), round(rows[(f0, reg)]['alpha'], 2),
                 int(rows[(f0, reg)]['genes'])) for reg in regions}
    want = {reg: (MK_TEXT[reg]['DoS'], MK_TEXT[reg]['alpha'], MK_TEXT[reg]['genes'])
            for reg in regions}
    ck.add(T, 'McDonald-Kreitman at MAF >= 0.10 reproduces the Results and fig. 6D '
           '(DoS, alpha, genes per class)', got == want,
           f'file {got}, text {want}')
    excess = [MK_FLOORS[f] for f in MK_FLOORS
              if rows[(f, 'eliminated')]['DoS'] > rows[(f, 'retained')]['DoS']]
    ck.add(T, 'eliminated DoS at or below retained DoS at every MAF floor (the Results: '
           '"at every threshold", no excess of adaptive substitution)', not excess,
           'excess at ' + ', '.join(excess) if excess else 'every floor')
    a_el = [rows[(f, 'eliminated')]['alpha'] for f in MK_FLOORS]
    a_txt = ', '.join(f'{MK_FLOORS[f]} {unsigned(a, 3)}' for f, a in zip(MK_FLOORS, a_el))
    ck.add(T, 'eliminated alpha positive at every MAF floor (the Results)',
           all(a > 0 for a in a_el), a_txt)
    lo, hi = MK_ALPHA_TEXT
    ck.add(T, f'eliminated alpha ranges from {lo:.2f} to {hi:.2f} across the MAF floors '
           f'(the Results)', (round(min(a_el), 2), round(max(a_el), 2)) == (lo, hi),
           a_txt)

    parts = []
    for f, label in MK_FLOORS.items():
        rt, el = rows[(f, 'retained')], rows[(f, 'eliminated')]
        parts.append(
            f"{label}: DoS {signed(rt['DoS'])} (95% CI {unsigned(rt['DoS_lo'])} to "
            f"{unsigned(rt['DoS_hi'])}) retained and {signed(el['DoS'])} "
            f"({unsigned(el['DoS_lo'])} to {unsigned(el['DoS_hi'])}) eliminated, "
            f"α {unsigned(rt['alpha'], 2)} and {unsigned(el['alpha'], 2)}")
    n_rt, n_el = genes['retained'][0], genes['eliminated'][0]
    # one sentence, floors separated by semicolons, no final full stop, as in the
    # curated rows
    result = '; '.join(parts + [f'pooled over {n_rt:,} retained and {n_el:,} eliminated '
                                f'genes'])
    if excess:
        verdict = ('Not rejected at every floor: eliminated DoS exceeds retained DoS at '
                   + ', '.join(excess))
    else:
        verdict = ('Rejected: no excess of adaptive substitution in eliminated genes at '
                   f'any floor; eliminated α {unsigned(min(a_el), 2)} to '
                   f'{unsigned(max(a_el), 2)}')
    row = {
        'alternative': 'Positive selection on eliminated genes',
        'prediction': ('Adaptive substitution, not relaxed constraint, raises dN/dS: '
                       'eliminated genes would show a higher direction of selection '
                       '(DoS) and α than retained genes'),
        'test': ('McDonald–Kreitman test pooled within each class, with polymorphism '
                 'counted from SNPs at minor-allele frequencies of at least 0.10, 0.15 '
                 'and 0.20'),
        'result': result,
        'verdict': verdict,
        'shown_in': 'Fig. 6D (MAF ≥ 0.10)',
    }
    notes = [f"MK {MK_FLOORS[f]}: DoS {rows[(f, 'retained')]['DoS']:+.4f} / "
             f"{rows[(f, 'eliminated')]['DoS']:+.4f}, alpha "
             f"{rows[(f, 'retained')]['alpha']:.3f} / {rows[(f, 'eliminated')]['alpha']:.3f} "
             f"(retained / eliminated)" for f in MK_FLOORS]
    return row, notes



def build_popgen_alternatives(ctx):
    T = 'popgen_alternatives'
    ck = ctx.checks
    alt = ctx.tsv('popgen_alternatives', keep_default_na=False)
    need(alt, ['alternative', 'prediction', 'test', 'result', 'verdict', 'shown_in'],
         ctx.src['popgen_alternatives'])
    empty = int((alt[['alternative', 'test', 'result', 'verdict']] == '').any(axis=1).sum())
    ck.add(T, 'every row of the curated file has an alternative, test, result and verdict',
           empty == 0, f'{empty} incomplete rows')
    figs = sorted({int(n) for s in alt['shown_in'] for n in re.findall(r'\bFig\. (\d+)', str(s))})
    ck.add(T, 'the curated shown_in column cites main figures 1 to 6 only (numbering of '
           '2026-10-27)', all(1 <= n <= 6 for n in figs),
           'main figures cited: ' + (', '.join(map(str, figs)) or 'none'))
    if (alt['alternative'].str.startswith('Positive selection')).any():
        sys.exit(f"ERROR: {ctx.src['popgen_alternatives']}: the positive-selection row is "
                 f"built from mk_summary; remove it from the curated file.")

    # The McDonald-Kreitman row, from the analysis output rather than typed in,
    # placed after the site-counting row (the order of the Results).
    row, mk_notes = mk_alternative_row(ctx.tsv('mk_summary'), ck, T, ctx.src['mk_summary'])
    after = [i for i, a in enumerate(alt['alternative']) if a.startswith('Nei–Gojobori')]
    ck.add(T, 'one site-counting row in the curated file, for the McDonald-Kreitman row to '
           'follow', len(after) == 1, f'{len(after)} found')
    pos = after[0] + 1 if len(after) == 1 else len(alt)
    alt = pd.concat([alt.iloc[:pos], pd.DataFrame([row], columns=alt.columns),
                     alt.iloc[pos:]], ignore_index=True)
    cols = [
        Col('alternative', 'alternative', 'Explanation other than relaxed purifying '
            'selection, or technical bias', 'text'),
        Col('prediction', 'prediction', 'What it predicts if it produced the contrast', 'text'),
        Col('test', 'test', 'How it was tested or controlled', 'text'),
        Col('result', 'result', 'Result on the primary callset; folds are eliminated '
            'relative to retained', 'text'),
        Col('verdict', 'verdict', 'Conclusion', 'text', values=True),
        Col('shown_in', 'shown_in', 'Figure showing the result (— if text only)', 'text'),
    ]
    title = ('Alternative explanations for the elevated πN/πS and dN/dS of eliminated genes, '
             'and the tests applied to them.')
    caption = (
        'Each row names an explanation other than relaxed purifying selection, or a '
        'technical bias, that could produce or inflate the eliminated-versus-retained '
        'πN/πS contrast in the population sample, or, for positive selection, the dN/dS '
        'contrast; how it was tested or controlled; the result on the primary callset; and '
        'the conclusion. For positive selection, the direction of selection (DoS) is given '
        'with its 95% interval from 2,000 resamples of genes, and α = 1 − NI is the estimated '
        'fraction of nonsynonymous substitutions driven by positive selection '
        '(Supplementary Methods).')
    part = Part('', '', alt, cols, wrap=True, sources=('popgen_alternatives', 'mk_summary'))
    return Table(T, title, caption, [part], [f'{len(alt)} rows'] + mk_notes)


# ---------------------------------------------------------------------------
# S11 X chromosomes versus autosomes
# ---------------------------------------------------------------------------

def build_x_autosome(ctx):
    T = 'x_autosome'
    ck = ctx.checks
    xa, hx = ctx.tsv('x_autosome'), ctx.tsv('x_dosage')
    need(xa, ['chrom', 'class', 'n_retained', 'n_eliminated', 'median_dS', 'median_dN',
              'median_dN_dS'], ctx.src['x_autosome'])
    need(hx, ['chrom', 'class', 'length_mb', 'n_bins', 'n_usable', 'n_eliminated', 'n_low',
              'n_high', 'median_ratio', 'q25', 'q75', 'pooled_ratio', 'frac_deviant', 'call',
              'male_copies'], ctx.src['x_dosage'])
    j = xa.merge(hx, on='chrom', how='outer', suffixes=('', '_hic'), indicator=True)
    ck.add(T, 'the two tables cover the same chromosomes', (j['_merge'] == 'both').all(),
           ', '.join(j.loc[j['_merge'] != 'both', 'chrom'].astype(str)) or 'all matched')
    j = j[j['_merge'] == 'both'].copy()
    ck.add(T, 'X/autosome classes agree', bool((j['class'] == j['class_hic']).all()))
    drop = {'class_hic': 'equals class (checked)', '_merge': 'join indicator'}
    for c in ('n_dS', 'n_dN', 'n_dN_dS'):
        if c in j.columns:
            same = bool((j[c] == j['n_retained']).all())
            ck.add(T, f'{c} equals n_retained', same)
            drop[c] = 'equals n_retained (checked)'
    if 'used' in j.columns:
        u = to_bool(j['used'])
        ck.add(T, "every chromosome 'used' in the permutation test",
               bool(u.fillna(False).all()))
        drop['used'] = 'True for every chromosome (checked)'
    for c in [c for c in j.columns if c.startswith('ratio_')]:
        ck.add(T, f'{c} equals the pooled ratio',
               close(j[c], j['pooled_ratio'], rtol=PRINTED_RTOL).all())
        drop[c] = 'equals the pooled ratio (one testis library; checked)'
    ck.add(T, 'male copies = 2 × median ratio',
           close(j['male_copies'], 2 * j['median_ratio'], rtol=PRINTED_RTOL).all())
    j['class_label'] = j['class'].map({'A': 'autosome', 'X': 'X'}).fillna(j['class'])
    j['pct_deviant'] = 100 * j['frac_deviant']
    order = xa['chrom'].tolist()
    present = set(j['chrom'])
    j = j.set_index('chrom').loc[[c for c in order if c in present]].reset_index()
    drop.update({'class': 'shown as a word', 'frac_deviant': 'shown as a percentage'})
    cols = [
        Col('chromosome', 'chrom', 'Ascaris germline chromosome', 'text'),
        Col('class', 'class_label', 'Autosome or X chromosome', 'text', values=True),
        Col('length_Mb', 'length_mb', 'Chromosome length, Mb', 'f3'),
        Col('n_retained_genes', 'n_retained', f"Retained genes (supplementary table "
            f"{NUMBER['rates_pairwise']})", 'int'),
        Col('n_eliminated_genes', 'n_eliminated', 'Eliminated genes', 'int'),
        Col('median_dS_retained', 'median_dS', 'Median dS over retained genes', 'f4'),
        Col('median_dN_retained', 'median_dN', 'Median dN over retained genes', 'f4'),
        Col('median_dN_dS_retained', 'median_dN_dS', 'Median dN/dS over retained genes',
            'f4'),
        Col('HiC_bins', 'n_bins', '100 kb bins on the chromosome', 'int'),
        Col('HiC_bins_usable', 'n_usable', 'Bins used for the ratio', 'int'),
        Col('HiC_bins_eliminated', 'n_eliminated_hic', 'Bins overlapping eliminated DNA '
            '(excluded)', 'int'),
        Col('HiC_bins_low_density', 'n_low', 'Bins below 0.2 × the median autosomal '
            'read-end density in either library (excluded)', 'int'),
        Col('HiC_bins_high_density', 'n_high', 'Bins above 5 × that density in either '
            'library (excluded)', 'int'),
        Col('testis_ovary_ratio_median', 'median_ratio', 'Median over usable bins of '
            '(testis read ends / testis autosomal total) / (ovary read ends / ovary '
            'autosomal total)', 'f3'),
        Col('testis_ovary_ratio_Q1', 'q25', 'First quartile of the per-bin ratio', 'f3'),
        Col('testis_ovary_ratio_Q3', 'q75', 'Third quartile', 'f3'),
        Col('testis_ovary_ratio_pooled', 'pooled_ratio', 'Ratio of summed read ends over '
            'usable bins, each library scaled by its autosomal total', 'f3'),
        Col('pct_bins_deviant', 'pct_deviant', 'Percentage of usable bins at least 0.5 log2 '
            'units from the class expectation (0.5 for X, 1 for autosomes)', 'f1'),
        Col('male_copy_number', 'male_copies', '2 × median ratio', 'f3'),
        Col('call', 'call', 'Copy number in males, from the median ratio on a log2 scale',
            'text', values=True),
    ]
    a, x = j[j['class'] == 'A'], j[j['class'] == 'X']
    title = ('Ascaris germline chromosomes: rates of retained genes, and male-to-female '
             'Hi-C read depth.')
    caption = (
        'Per chromosome, counts of retained and eliminated genes among the Ascaris × '
        f"P. univalens pairs (supplementary table {NUMBER['rates_pairwise']}) and medians "
        'over retained genes. The Hi-C columns compare read-end coverage of testis and '
        'ovary Hi-C libraries from dissected gonads (Simmons et al. 2026; read ends at '
        'mapping quality ≥ 30) in 100 kb bins, each library scaled by its own usable '
        'autosomal total; bins overlapping eliminated DNA, or below 0.2 or above 5 times the '
        'median autosomal density in either library, are excluded. A ratio of 1 is expected '
        'where males and females carry the same number of copies and 0.5 where males carry '
        'one and females two. Each chromosome is called by the expectation nearer its median '
        'ratio on a log2 scale; a median within 0.25 log2 units of the midpoint would be '
        'called intermediate.')
    notes = [f"autosome median ratios {a['median_ratio'].min():.3f}–{a['median_ratio'].max():.3f}"
             f"; X {x['median_ratio'].min():.3f}–{x['median_ratio'].max():.3f}",
             'calls: ' + ', '.join(f'{k} {v}' for k, v in j['call'].value_counts().items())]
    part = Part('', '', j, cols, drop=drop, sources=('x_autosome', 'x_dosage'))
    return Table(T, title, caption, [part], notes)


BUILDERS = {
    'rates_pairwise': build_rates_pairwise,
    'branch_models': build_branch_models,
    'omega_above_one': build_omega_above_one,
    'codon_usage': build_codon_usage,
    'insulation_models': build_insulation_models,
    'rate_heterogeneity': build_rate_heterogeneity,
    'toxocara_models': build_toxocara_models,
    'microsynteny': build_microsynteny,
    'popgen_cohort': build_popgen_cohort,
    'popgen_alternatives': build_popgen_alternatives,
    'x_autosome': build_x_autosome,
}


# ---------------------------------------------------------------------------
# From parts to published frames
# ---------------------------------------------------------------------------

def is_blank(s):
    return s.isna() | (s.astype(str).str.strip() == '')


def finalize(table, part, dropped_log):
    """Select, order, rename and type the columns of one part.

    Returns (frame with published names, list of the Col specs kept). Stops if a
    spec names a column the builder did not make; logs every source column that
    is left out, with the reason when the builder gave one.
    """
    f = part.frame.reset_index(drop=True)
    label = NUMBER[table.key] + part.suffix
    missing = [c.src for c in part.cols if c.src not in f.columns]
    if missing:
        sys.exit(f'ERROR: [{label}] columns {missing} were not built (internal error)')
    used = {c.src for c in part.cols}
    for c in f.columns:
        if c in used:
            continue
        why = part.drop.get(c)
        if why is None:
            logging.warning('[%s] source column left out without a stated reason: %s',
                            label, c)
            why = 'not used'
        dropped_log.append((label, c, why))
    kept, data = [], {}
    for c in part.cols:
        s = f[c.src]
        if c.optional and is_blank(s).all():
            dropped_log.append((label, c.src, 'empty in every row'))
            continue
        if c.fmt == 'int':
            v = pd.to_numeric(s, errors='coerce')
            frac = v.dropna() % 1
            if (frac != 0).any():
                sys.exit(f'ERROR: [{label}] {c.name} has non-integer values')
            s = v.round().astype('Int64')
        elif c.fmt == 'bool':
            s = to_bool(s)
        elif c.fmt in ('text',):
            s = s.where(~is_blank(s), None)
        kept.append(c)
        data[c.name] = s.values
    out = pd.DataFrame(data, columns=[c.name for c in kept])
    if out.columns.duplicated().any():
        sys.exit(f'ERROR: [{label}] duplicate published column names')
    return out, kept


def definition(col, series):
    d = col.desc
    if col.values:
        vals = sorted({str(v) for v in series.dropna()})
        if vals:
            d += '. Values: ' + ', '.join(vals)
    return d


def tsv_cell(v, fmt):
    """One TSV field. Floats keep ten significant digits (more than any source
    carries), distances one decimal (they end in .5), P values six digits, with
    an underflowed P written as P_UNDERFLOW. Negative zero prints as 0."""
    if v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)):
        return ''
    if fmt == 'p':
        v = float(v)
        return P_UNDERFLOW if v == 0 else f'{v:.6g}'
    if fmt == 'int':
        return str(int(v))
    if fmt == 'bool':
        return 'True' if bool(v) else 'False'
    if fmt == 'dist':
        return f'{float(v) + 0.0:.1f}'
    if isinstance(v, (float, np.floating)):
        v = float(v)
        return f'{(0.0 if v == 0 else v):.10g}'
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return str(v)


def write_tsv(out, cols, path):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(out.columns) + '\n')
        fmts = [c.fmt for c in cols]
        for row in out.itertuples(index=False, name=None):
            fh.write('\t'.join(tsv_cell(v, f) for v, f in zip(row, fmts)) + '\n')


NUMFMT = {'int': '0', 'dist': '0.0', 'f1': '0.0', 'f2': '0.00', 'f3': '0.000', 'f4': '0.0000',
          'v4': '0.0###', 'pct': '0.0', 'p': '[<0.001]0.00E+00;0.000', 'g': 'General',
          'text': 'General', 'bool': 'General'}


def xl_value(v, fmt):
    if v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)):
        return None
    if fmt == 'p' and float(v) == 0:
        return P_UNDERFLOW
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        v = float(v)
        return 0.0 if v == 0 else v
    return v


def display_text(v, fmt):
    """Roughly what the workbook shows for a value, to size its column."""
    if v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)):
        return ''
    try:
        if fmt == 'p':
            v = float(v)
            return P_UNDERFLOW if v == 0 else (f'{v:.2E}' if v < 0.001 else f'{v:.3f}')
        if fmt == 'int':
            return str(int(v))
        if fmt in ('dist', 'f1', 'pct'):
            return f'{float(v):.1f}'
        if fmt in ('f2', 'f3', 'f4'):
            return f'{float(v):.{fmt[1]}f}'
        if fmt == 'v4':
            return f'{float(v):.4f}'.rstrip('0')
        if fmt == 'bool':
            return 'FALSE'
        if isinstance(v, (float, np.floating)):
            return f'{float(v):.10g}'
    except (TypeError, ValueError):
        pass
    return str(v)


def text_width(series):
    lens = series.dropna().astype(str).str.len()
    return int(lens.quantile(0.95)) if len(lens) else 0


def write_workbook(entries, path, contents_rows, column_rows):
    """entries: list of (sheet name, heading, caption, frame, cols, part)."""
    wb = openpyxl.Workbook()
    bold = Font(bold=True)
    head_border = Border(bottom=Side(style='thin'))
    wrap_top = Alignment(wrap_text=True, vertical='top')

    ws = wb.active
    ws.title = 'Contents'
    ws.append(['Supplementary tables'])
    ws['A1'].font = Font(bold=True, size=12)
    ws.append([])
    hdr = ['Table', 'Sheet', 'Title', 'Rows', 'Columns', 'Caption']
    ws.append(hdr)
    for cell in ws[3]:
        cell.font, cell.border = bold, head_border
    for r in contents_rows:
        ws.append(r)
    for row in ws.iter_rows(min_row=4):
        for cell in row:
            cell.alignment = wrap_top
    for col, w in zip('ABCDEF', (7, 7, 60, 8, 8, 110)):
        ws.column_dimensions[col].width = w
    ws.freeze_panes = 'A4'

    for sheet, heading, caption, out, cols, part in entries:
        ws = wb.create_sheet(sheet)
        ncol = max(len(cols), 1)
        ws.cell(row=1, column=1, value=heading).font = Font(bold=True, size=11)
        ws.cell(row=2, column=1, value=caption).alignment = wrap_top
        # widths first: the caption's row height depends on them
        widths = []
        for i, c in enumerate(cols, 1):
            shown = out[c.name].map(lambda v, f=c.fmt: display_text(v, f))
            w = max(len(c.name), text_width(shown) if len(out) else 0) + 2
            w = min(max(w, 8), 60 if part.wrap else 45)
            widths.append(w)
            ws.column_dimensions[get_column_letter(i)].width = w
        span = max(min(ncol, 12), 6)
        ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=span)
        total = sum(widths[:span]) + 8 * max(0, span - len(widths))
        lines = max(1, math.ceil(len(caption) / (1.2 * max(total, 40))))
        ws.row_dimensions[2].height = 15 * lines + 4
        for i, c in enumerate(cols, 1):
            cell = ws.cell(row=4, column=i, value=c.name)
            cell.font, cell.border = bold, head_border
            cell.alignment = Alignment(wrap_text=True, vertical='bottom')
        fmts = [c.fmt for c in cols]
        for r_i, row in enumerate(out.itertuples(index=False, name=None), start=5):
            for c_i, (v, f) in enumerate(zip(row, fmts), start=1):
                cell = ws.cell(row=r_i, column=c_i, value=xl_value(v, f))
                if f not in ('g', 'text', 'bool') and cell.value is not None \
                        and not isinstance(cell.value, str):
                    cell.number_format = NUMFMT[f]
                if part.wrap:
                    cell.alignment = wrap_top
        last = get_column_letter(ncol)
        ws.freeze_panes = ws.cell(row=5, column=1 + part.freeze_cols)
        if len(out):
            ws.auto_filter.ref = f'A4:{last}{4 + len(out)}'

    ws = wb.create_sheet('Columns')
    ws.append(['Table', 'Sheet', 'Column', 'Definition'])
    for cell in ws[1]:
        cell.font, cell.border = bold, head_border
    for r in column_rows:
        ws.append(r)
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = wrap_top
    for col, w in zip('ABCD', (7, 7, 36, 110)):
        ws.column_dimensions[col].width = w
    ws.freeze_panes = 'A2'
    tmp = path + '.tmp'
    wb.save(tmp)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def read_manifest(path, base_dir, overrides):
    """key -> absolute path, after ${VAR} expansion and --source overrides."""
    raw = {}
    with open(path, encoding='utf-8') as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = [p.strip() for p in s.split('\t') if p.strip()]
            if len(parts) != 2:
                sys.exit(f'ERROR: {path}:{ln}: expected key<TAB>path, got {s!r}')
            if parts[0] in raw:
                sys.exit(f'ERROR: {path}:{ln}: key {parts[0]!r} given twice')
            raw[parts[0]] = parts[1]
    for ov in overrides or []:
        key, sep, p = ov.partition('=')
        if not sep or not key or not p:
            sys.exit(f'ERROR: --source expects key=path, got {ov!r}')
        raw[key.strip()] = p.strip()
    out, bad = {}, []
    for key, p in raw.items():
        e = os.path.expanduser(os.path.expandvars(p))
        if '$' in e:
            bad.append(f'  {key}: unexpanded variable in {p!r}')
            continue
        if not os.path.isabs(e):
            e = os.path.join(base_dir, e)
        out[key] = os.path.normpath(e)
    if bad:
        sys.exit('ERROR: export the variable(s) or pass --source key=path:\n' + '\n'.join(bad))
    return out


def select_tables(only):
    if not only:
        return list(TABLE_ORDER)
    by_num = {v.upper(): k for k, v in NUMBER.items()}
    keys = []
    for tok in [t.strip() for t in only.split(',') if t.strip()]:
        k = tok if tok in NUMBER else by_num.get(tok.upper())
        if k is None:
            sys.exit(f'ERROR: --only: unknown table {tok!r}; use a key ({", ".join(NUMBER)}) '
                     f'or a number ({", ".join(NUMBER.values())})')
        if k not in keys:
            keys.append(k)
    return [k for k in TABLE_ORDER if k in keys]


def check_sources(keys, src):
    problems = []
    for k in keys:
        for s in REQUIRES[k]:
            if s not in src:
                problems.append(f'  {NUMBER[k]} {k}: {s!r} not in the manifest')
            elif not os.path.isfile(src[s]) or os.path.getsize(src[s]) == 0:
                problems.append(f'  {NUMBER[k]} {k}: {s}: missing or empty: {src[s]}')
        for s in OPTIONAL.get(k, []):
            if s in src and not os.path.isfile(src[s]):
                logging.warning('optional source %s for %s not found (%s); its cross-check '
                                'is skipped', s, NUMBER[k], src[s])
    if problems:
        sys.exit('ERROR: inputs not ready:\n' + '\n'.join(problems))


def rel(path, base):
    try:
        r = os.path.relpath(path, base)
        return path if r.startswith('..') else r
    except ValueError:
        return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Assemble the supplementary tables (S1-S11) from analysis outputs.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument('--manifest', default=os.path.join(_HERE, 'supplementary_table_sources.tsv'),
                   help='key<TAB>path source list (default: beside this script)')
    p.add_argument('--base-dir', default=os.getcwd(),
                   help='Directory relative manifest paths resolve against '
                        '(default: the current directory)')
    p.add_argument('--source', action='append', default=[], metavar='KEY=PATH',
                   help='Override one manifest entry; repeatable')
    p.add_argument('--only', default=None,
                   help='Comma-separated table keys or numbers to build (default: all)')
    p.add_argument('--omega-cap', type=float, default=5.0,
                   help='Per-branch dN/dS cap used by the branch models (default 5)')
    p.add_argument('--cohort-threshold', type=float, default=0.40,
                   help='Germline-fraction threshold used to select the population '
                        'sample (default 0.40); checked against the selection')
    p.add_argument('--workbook-name', default='Supplementary_Tables.xlsx')
    p.add_argument('--outdir', required=True)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(os.path.join(args.outdir, 'tsv'), exist_ok=True)
    setup_logging(args.outdir, 'make_supplementary_tables')
    logging.info('make_supplementary_tables.py %s', SCRIPT_VERSION)
    base_dir = os.path.abspath(args.base_dir)
    logging.info('manifest: %s', os.path.abspath(args.manifest))
    logging.info('base dir: %s', base_dir)
    src = read_manifest(args.manifest, base_dir, args.source)
    keys = select_tables(args.only)
    check_sources(keys, src)
    ctx = Context(src, args)

    tables = []
    for k in keys:
        logging.info('--- %s %s ---', NUMBER[k], k)
        tables.append(BUILDERS[k](ctx))

    entries, contents, column_rows, dropped, written = [], [], [], [], []
    caption_blocks, summary_rows = [], []
    for t in tables:
        num = NUMBER[t.key]
        caption_blocks.append(f'Supplementary Table {num}. {t.title}\n{t.caption}\n')
        for i, part in enumerate(t.parts):
            out, cols = finalize(t, part, dropped)
            sheet = num + part.suffix
            heading = f'Supplementary Table {num}. {t.title}' + \
                (f' {part.title}' if part.title else '')
            entries.append((sheet, heading, t.caption, out, cols, part))
            contents.append([num if i == 0 else '', sheet,
                             t.title + (f' {part.title}' if part.title else ''),
                             len(out), len(cols), t.caption if i == 0 else ''])
            for c in cols:
                column_rows.append([num, sheet, c.name, definition(c, out[c.name])])
            fn = os.path.join(args.outdir, 'tsv', f'Table_{sheet}_{t.key}.tsv')
            write_tsv(out, cols, fn)
            written.append(fn)
            summary_rows.append((sheet, len(out), len(cols),
                                 t.title + (f' {part.title}' if part.title else '')))
            logging.info('  %s: %s rows x %d columns', sheet, f'{len(out):,}', len(cols))

    wb_path = os.path.join(args.outdir, args.workbook_name)
    write_workbook(entries, wb_path, contents, column_rows)
    written.append(wb_path)
    cap_path = os.path.join(args.outdir, 'table_captions.txt')
    with open(cap_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(caption_blocks))
    written.append(cap_path)

    # --- summary ---------------------------------------------------------------------
    ck = ctx.checks
    used_keys = sorted({s for k in keys for s in REQUIRES[k] + OPTIONAL.get(k, [])
                        if s in src})
    L = [f'make_supplementary_tables.py {SCRIPT_VERSION}   '
         f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
         f'Output: {os.path.abspath(args.outdir)}', '',
         f'CHECKS: {len(ck.items) - ck.n_failed} passed, {ck.n_failed} failed']
    for c in ck.items:
        if not c['ok']:
            L.append(f"  FAILED [{c['table']}] {c['check']}" +
                     (f" ({c['detail']})" if c['detail'] else ''))
    L += ['', 'Numbering follows the order of first citation in the text (2026-10-12). If the '
          'text changes that order, edit NUMBER (and FIGURE) in this script and rerun.', '',
          f"{'Sheet':<6} {'Rows':>7} {'Cols':>5}  Title"]
    L += [f'{s:<6} {r:>7,} {c:>5}  {ti}' for s, r, c, ti in summary_rows]
    L += ['', 'Numbers for the text and legends']
    for t in tables:
        if t.notes:
            L.append(f'  {NUMBER[t.key]}')
            L += [f'    {n}' for n in t.notes]
    L += ['', 'Source columns left out of the tables']
    for label, c, why in dropped:
        L.append(f'  {label:<5} {c:<26} {why}')
    L += ['', 'Checks']
    for c in ck.items:
        L.append(f"  [{'OK' if c['ok'] else 'MISMATCH'}] {c['table']:<4} {c['check']}"
                 + (f"  ({c['detail']})" if c['detail'] else ''))
    L += ['', 'Inputs']
    fps = {}
    for k in used_keys:
        fp = fingerprint(src[k])
        fps[k] = fp
        L.append(f"  {k:<26} {rel(src[k], base_dir)}"
                 + (f"  sha256:{fp['sha256_16']}" if fp else '  (not found)'))
    sum_path = os.path.join(args.outdir, 'supplementary_tables_summary.txt')
    with open(sum_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    written.append(sum_path)

    prov = {
        'script': os.path.abspath(__file__), 'script_version': SCRIPT_VERSION,
        'git_commit': git_commit(_HERE),
        'run_utc': datetime.now(timezone.utc).isoformat(),
        'host': platform.node(), 'python': platform.python_version(),
        'versions': {'numpy': np.__version__, 'pandas': pd.__version__,
                     'scipy': scipy_version, 'openpyxl': openpyxl.__version__},
        'args': vars(args), 'base_dir': base_dir,
        'numbering': {k: NUMBER[k] for k in keys},
        'inputs': fps,
        'tables': {s: {'rows': r, 'columns': c} for s, r, c, _ in summary_rows},
        'checks': ck.items,
        'outputs': [os.path.abspath(p) for p in written],
    }
    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump(prov, fh, indent=2, default=str)
    logging.info('Wrote %s and %d TSVs', wb_path, len(written) - 3)
    if ck.n_failed:
        logging.warning('%d check(s) failed; see %s', ck.n_failed, sum_path)
        return 2
    logging.info('All %d checks passed', len(ck.items))
    return 0


if __name__ == '__main__':
    sys.exit(main())
