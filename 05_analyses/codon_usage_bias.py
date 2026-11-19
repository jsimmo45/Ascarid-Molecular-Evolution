#!/usr/bin/env python3
"""
codon_usage_bias.py - codon usage bias in retained versus eliminated genes.

Tests whether genes removed by programmed DNA elimination show relaxed
translational selection, using three metrics computed per gene from CDS:

    ENC    effective number of codons (Wright 1990). 20 = one codon per amino
           acid, 61 = every synonymous codon used equally. Lower is more biased.
    GC3s   GC content at synonymous third positions.
    dENC   ENC expected from GC3s alone, minus observed ENC. Positive means more
           bias than base composition explains, which is the signature of
           translational selection rather than mutational pressure.
    CAI    codon adaptation index (Sharp and Li 1987), relative to a reference
           set of highly expressed genes.

**dENC is the metric to lead with.** ENC and GC3s are not independent: eliminated
genes have higher GC3s, and GC3s alone shifts ENC, so raw ENC confounds
composition with selection. dENC removes the composition term.

Reads a CDS FASTA and a table carrying per-gene region calls, and writes:

    codon_usage_panel.{svg,pdf,png}     Figure 2B-C strip, 185 x 60 mm
    enc_vs_gc3.{svg,pdf,png}            the ENC-GC3 panels alone
    cai_boxplot.{svg,pdf,png}           }
    denc_boxplot.{svg,pdf,png}          } single metrics, for supplements
    enc_boxplot.{svg,pdf,png}           }
    gc3_boxplot.{svg,pdf,png}           }
    codon_usage_per_gene.tsv            per-gene ENC, GC3s, dENC, CAI
    codon_usage_summary.tsv             medians, IQRs and tests
    codon_usage_summary.txt             the same, formatted for reading
    provenance.json                     inputs, arguments and versions

Notes on method
---------------
**GC3s, not GC3 over all codons.** Wright's expected-ENC curve is defined on GC
content at *synonymous* third positions. Methionine (ATG) and tryptophan (TGG)
have no synonymous third position and always carry G there, so including them
inflates GC3 by roughly two to three points and, because Met and Trp content
varies between genes, adds noise rather than a constant offset. On the steep part
of the curve, where these genomes sit, that is worth around one ENC unit in dENC
per gene, and it does not cancel between groups because the two groups differ in
GC3s. Both quantities are written to the per-gene table; GC3s is what is plotted
and what dENC uses.

**The CAI reference set should come from expression, not from codon bias.** If
the reference is chosen as the most codon-biased retained genes, those genes have
high CAI by construction and the retained-versus-eliminated contrast is partly
circular. Pass `--ref-genes` with a list built from expression data
(`make_cai_reference.py` does this). Without it the script falls back to the most
biased retained genes, prints a warning, and marks the summary accordingly --
usable as a relative check within one species, not as a headline result.

**ENC has a ceiling.** Genes with no detectable bias can return an estimate above
61, which is clamped. The clamped count is reported per group; where it is large,
compare medians rather than means, since a ceiling truncates the upper tail.

**The GC3s-controlled block answers the obvious objection to CAI and ENC.** Both
move with GC3s, and the two regions differ in GC3s, so a raw contrast in either
cannot separate translational selection from base composition. Two controls are
run: matched sampling, which pairs each eliminated gene with a retained gene at
the same GC3s and repeats the draw many times, and a stratified rank test, which
keeps every gene. They fail in opposite directions, which is why both are here.

Calibrated on synthetic data in both directions. On genes built so that the CAI
difference is *entirely* a GC3s artifact, the raw test gives P = 3 x 10^-39 while
the matched test gives median P = 0.96 with 0% of replicates significant. On
genes built with a real difference in preference for optimal codons, chosen
GC3-balanced so it is orthogonal to composition, the matched test gives
P = 2 x 10^-46 with 100% of replicates significant. So it is neither blind nor
trigger-happy.

The stratified test needed the same treatment and initially failed it: with
quantile strata it returned P = 0.001 on the null data, because quantile bins are
widest where the data are densest and the eliminated genes then sit at the upper
edge of every bin. Strata are therefore fixed-width and tied to the caliper.
Even so its calibration degrades as the caliper widens -- at 0.02 it returns
P = 2 x 10^-4 on null data -- so the matched result is the primary one and the
script warns above 0.01.

Changes 2026-10-12 (figure review of Fig 2B-C and supplementary fig. S4)
-------------------------------------------------------------------------
* The corner value of each ENC-GC3s sub-panel is now the **median dENC over
  genes**, the quantity the boxplot beside it shows, the summary tests and the
  text quotes. It used to print the offset of the median centroid from Wright's
  curve (expected ENC at the median GC3s minus the median ENC), which is a
  different number because the curve is not linear (P. univalens retained 3.13
  against a median dENC of 2.64), so one figure carried two different values
  labelled dENC for the same genes. The dotted segment still shows the
  centroid's offset; the summary gives both numbers.
* Axis limits are shared by every species unless --auto-limits is given
  (--gc3s-lim, --enc-lim, --denc-lim), so the Figure 2 strip and the two rows of
  the supplementary figure are drawn on one scale. Points outside the limits are
  counted in the summary, per panel and group, rather than dropped silently.
* The boxplot bracket carries significance markers by default (--bracket-label
  stars): the rest of Figure 2 uses markers, and the figure specification allows
  one convention per figure. --bracket-label p restores the printed P.
No statistic changes: every test, the matched replicates and the centroids are
computed before any figure is drawn, from the same random stream as before.

Changes 2026-10-14 (boxes from every gene, as in Fig 2D)
--------------------------------------------------------
* Boxes and whiskers are computed from every gene (box_stats) and drawn with
  ax.bxp. ax.boxplot computed them from the drawn genes, after those outside
  the central 99% of the pooled values (and outside --denc-lim) had been
  removed, so the drawn box was that of a trimmed sample and its median line
  was not the median in the summary.
* The drawing range is widened to take in every whisker end, so no whisker is
  cut. With the shared axis (--denc-lim), a whisker end outside the limits
  widens the axis for that species and logs a warning; widen --denc-lim for
  every species to keep the Figure 2 strip and the supplementary rows on one
  scale.
* The summary lists each metric's quartiles, whisker ends and genes beyond the
  whiskers, and the axis range of each boxplot with the genes not drawn.
* --denc-lim defaults to -6.5 16.5 (was 16). At 16 the bracket raised the top
  to 16.09 and 16.10 for P. univalens and T. canis, so the three species were
  not quite on one scale.
No statistic changes: tests, medians and the matched replicates are as before.

Changes 2026-10-26 (Jianbin Wang's review of Fig 2)
---------------------------------------------------
* The Figure 2 strip reads left to right as ENC against GC3s (panel B, the two
  sub-panels) and then the dENC boxplot (panel C), which summarises them; the
  boxplot used to come first.
* Both groups' median centroids are filled circles in their region colour on
  every sub-panel. Colour identifies the group, so the open ring that marked
  the other group's centroid is gone, and the legend has one centroid entry.
* The dotted segment from each centroid to Wright's curve is gone. It showed
  the centroid's offset, a different number from the median dENC that the
  corner prints and the boxplot shows; the summary still reports the offset.
No statistic changes: every test, the centroids and the matched replicates are
computed as before.

Usage
-----
    python codon_usage_bias.py \\
      --cds-fasta As_cds.fa \\
      --rates-table kaks_regions/kaks_annotated.tsv \\
      --gene-col gene1 --region-col region1 \\
      --ref-genes cai_reference_ascaris.txt \\
      --outdir results/codon_usage_ascaris

Requires numpy, pandas, scipy and matplotlib, plus paper_colors.py, which is
located automatically in lib/.
"""

import argparse
import hashlib
import json
import logging
import os
import platform
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from scipy import stats

_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, 'lib'), os.path.join(os.path.dirname(_HERE), 'lib')):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, _cand)

from paper_colors import (COL, FONT, LINEWIDTH, NEUTRAL, REGION, apply_style,
                          audit, mm, panel_letter, save)

apply_style()

EXTS = ('.svg', '.pdf', '.png')

#: Box order, left to right: retained blue/left, eliminated red/right.
REGIONS = ('retained', 'eliminated')

#: The box_stats keys ax.bxp() reads.
BXP_KEYS = ('med', 'q1', 'q3', 'whislo', 'whishi', 'fliers')

# ---------------------------------------------------------------------------
# Genetic code (NCBI translation table 1)
# ---------------------------------------------------------------------------

STOP_CODONS = {'TAA', 'TAG', 'TGA'}

CODON_TABLE = {
    'F': ['TTT', 'TTC'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'],
    'M': ['ATG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'],
    'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'],
    'K': ['AAA', 'AAG'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'C': ['TGT', 'TGC'],
    'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
}

CODON_TO_AA = {c: aa for aa, cods in CODON_TABLE.items() for c in cods}

#: Amino acids with no synonymous third position. Excluded from GC3s and CAI.
NON_DEGENERATE = {'M', 'W'}

#: How many amino acid families fall in each degeneracy class, which is what
#: Wright's ENC expression sums over: 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6.
CLASS_SIZE = {2: 9, 3: 1, 4: 5, 6: 3}

ENC_MIN, ENC_MAX = 20.0, 61.0


# ---------------------------------------------------------------------------
# Sequence handling
# ---------------------------------------------------------------------------

def parse_fasta(path):
    """Read a FASTA file into {id: sequence}, keyed on the first token of the
    header. Sequences are upper-cased; line wrapping is irrelevant."""
    seqs, gid, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if gid is not None:
                    seqs[gid] = ''.join(chunks).upper()
                gid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if gid is not None:
        seqs[gid] = ''.join(chunks).upper()
    logging.info(f'Parsed {len(seqs):,} sequences from {path}')
    return seqs


def count_codons(seq):
    """Codon counts for an in-frame CDS, excluding stops and any codon
    containing an ambiguity character."""
    counts = Counter()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if codon in CODON_TO_AA:
            counts[codon] += 1
    return counts


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def compute_enc(codon_counts):
    """Effective number of codons (Wright 1990), homozygosity method.

    For each amino acid family of k synonymous codons with n observed codons,
    F = (n * sum(p_i^2) - 1) / (n - 1), and

        ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6

    where Fk is the mean F over families of degeneracy k.

    Returns (enc, clamped). ``clamped`` is True when the raw estimate exceeded
    61, which happens for genes with no detectable bias and is reported rather
    than hidden, since a ceiling truncates the upper tail of the distribution.

    A missing degeneracy class is interpolated, not assumed. Wright's own
    treatment of an unestimable F3 -- the three-fold class holds isoleucine
    alone, so a gene short on isoleucine loses it -- is to take the mean of F2
    and F4. Substituting the class size instead, which is the tempting
    shortcut, is equivalent to asserting F = 1 for that class: maximum bias, and
    a downward push on ENC of up to 1.6 units for the single amino acid it is
    meant to stand in for.
    """
    family_F = defaultdict(list)
    for aa, codons in CODON_TABLE.items():
        deg = len(codons)
        if deg == 1:
            continue
        n_total = sum(codon_counts.get(c, 0) for c in codons)
        if n_total <= 1:
            continue
        homo = sum((codon_counts.get(c, 0) / n_total) ** 2 for c in codons)
        family_F[deg].append((n_total * homo - 1) / (n_total - 1))

    mean_F = {k: float(np.mean(v)) for k, v in family_F.items() if v}

    # F3 rests on one amino acid, so interpolate it when absent.
    if 3 not in mean_F and 2 in mean_F and 4 in mean_F:
        mean_F[3] = (mean_F[2] + mean_F[4]) / 2.0
    # Without the 2-, 4- or 6-fold classes there is not enough information for
    # an estimate at all; say so rather than invent one.
    if not all(k in mean_F for k in (2, 4, 6)):
        return None, False

    enc = 2.0
    for deg, size in CLASS_SIZE.items():
        f = mean_F.get(deg)
        if f is None or f <= 0:
            # F <= 0 means usage even beyond the unbiased expectation, i.e. the
            # class contributes its full codon count.
            enc += size * deg
        else:
            enc += size / f

    clamped = enc > ENC_MAX
    return float(min(max(enc, ENC_MIN), ENC_MAX)), bool(clamped)


def compute_gc3(codon_counts, synonymous_only=True):
    """GC at third codon positions.

    With ``synonymous_only`` (the default, and what Wright's expected-ENC curve
    is defined on) methionine and tryptophan are excluded: they have no
    synonymous third position and always carry G there, so including them
    inflates the value by two to three points and by an amount that varies with
    each gene's Met and Trp content.
    """
    gc = total = 0
    for codon, n in codon_counts.items():
        if synonymous_only and CODON_TO_AA.get(codon) in NON_DEGENERATE:
            continue
        total += n
        if codon[2] in ('G', 'C'):
            gc += n
    return gc / total if total else None


def enc_expected(gc3s):
    """ENC expected from GC3s alone under mutational bias (Wright 1990)."""
    s = np.asarray(gc3s, dtype=float)
    return 2.0 + s + 29.0 / (s ** 2 + (1.0 - s) ** 2)


def compute_rscu(pooled):
    """Relative synonymous codon usage from pooled counts."""
    rscu = {}
    for aa, codons in CODON_TABLE.items():
        k = len(codons)
        total = sum(pooled.get(c, 0) for c in codons)
        for c in codons:
            rscu[c] = 1.0 if total == 0 else (pooled.get(c, 0) * k) / total
    return rscu


def compute_w_table(pooled):
    """Relative adaptiveness w = RSCU / max(RSCU) within each family.

    Codons unobserved in the reference get a half-count rather than zero. Sharp
    and Li's own recommendation: a literal zero makes log(w) undefined, and
    substituting an arbitrarily tiny number turns one rare codon into an
    unbounded penalty on the gene's whole CAI.
    """
    pooled = Counter({c: pooled.get(c, 0) for c in CODON_TO_AA})
    for c in CODON_TO_AA:
        if pooled[c] == 0:
            pooled[c] = 0.5
    rscu = compute_rscu(pooled)
    w = {}
    for aa, codons in CODON_TABLE.items():
        peak = max(rscu[c] for c in codons) or 1.0
        for c in codons:
            w[c] = rscu[c] / peak
    return w


def compute_cai(codon_counts, w_table):
    """CAI: geometric mean of w over a gene's codons, excluding Met and Trp."""
    log_sum, n = 0.0, 0
    for codon, count in codon_counts.items():
        if CODON_TO_AA.get(codon) in NON_DEGENERATE:
            continue
        log_sum += count * np.log(w_table[codon])
        n += count
    return float(np.exp(log_sum / n)) if n else None


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def fmt_p(p):
    """Journal-style p-value: 1.9 x 10^-68, never 1.90e-68."""
    if p is None or not np.isfinite(p):
        return 'n/a'
    if p <= 0:
        return 'P < 10$^{-300}$'
    if p >= 1e-3:
        return f'P = {p:.3f}'
    exp = int(np.floor(np.log10(p)))
    return f'P = {p / 10 ** exp:.1f} × 10$^{{{exp}}}$'


def stars(p):
    """Significance marker, the convention of Figure 2 (dnds_decomposition.py)."""
    if p is None or not np.isfinite(p):
        return 'n/a'
    for cut, mark in ((1e-4, '****'), (1e-3, '***'), (1e-2, '**'), (5e-2, '*')):
        if p < cut:
            return mark
    return 'ns'


def fmt_p_plain(p):
    if p is None or not np.isfinite(p):
        return 'n/a'
    if p <= 0:
        # Underflowed to exactly zero. log10(0) is -inf and int() of that
        # raises, so a large-n test with a strong effect would crash the run
        # at the formatting step. Report the floor honestly instead.
        return '< 10^-300'
    if p >= 1e-3:
        return f'{p:.3f}'
    exp = int(np.floor(np.log10(p)))
    return f'{p / 10 ** exp:.2f} x 10^{exp}'


def van_elteren(df, metric, strata_col, min_per_stratum=5):
    """Stratified Wilcoxon rank-sum test (van Elteren 1960).

    Ranks within each GC3s stratum and combines the per-stratum statistics with
    the locally-best weights 1/(N+1), so a difference in CAI can only register
    if it holds at matched GC3s. Uses every gene, unlike the matched test, at
    the cost of assuming the effect has one direction across strata.

    Returns (Z, p, n_strata_used) or None.
    """
    sub = df.dropna(subset=[metric, strata_col])
    num = den = 0.0
    used = 0
    for _, block in sub.groupby(strata_col, observed=True):
        ret = block.loc[block.region == 'retained', metric].to_numpy()
        eli = block.loc[block.region == 'eliminated', metric].to_numpy()
        n1, n2 = len(ret), len(eli)
        if n1 < min_per_stratum or n2 < min_per_stratum:
            continue
        N = n1 + n2
        ranks = stats.rankdata(np.concatenate([ret, eli]))
        W = ranks[n1:].sum()
        expected = n2 * (N + 1) / 2.0
        _, counts = np.unique(np.concatenate([ret, eli]), return_counts=True)
        tie = (counts ** 3 - counts).sum()
        var = (n1 * n2 / 12.0) * ((N + 1) - tie / (N * (N - 1.0)))
        if var <= 0:
            continue
        w = 1.0 / (N + 1.0)
        num += w * (W - expected)
        den += w * w * var
        used += 1
    if used == 0 or den <= 0:
        return None
    z = num / np.sqrt(den)
    return float(z), float(2 * stats.norm.sf(abs(z))), used


def gc3s_matched_test(df, metric, rng, caliper=0.01, reps=500,
                      min_matched=50):
    """Compare a metric between regions at matched GC3s.

    Each eliminated gene is paired with an unused retained gene whose GC3s is
    within ``caliper``, and the comparison is run on the matched sets. Which
    retained genes get picked is arbitrary -- there are 9,085 of them and 322
    slots -- so a single draw is one sample from a large space of equally valid
    matched sets, and quoting its p-value would be quoting a lottery. The whole
    procedure is therefore repeated and the distribution of outcomes reported.

    The post-match GC3s test is reported alongside: if the match worked it is
    non-significant, and if it is not, nothing below it means anything.
    """
    ret = df[df.region == 'retained'].dropna(subset=['GC3s', metric])
    eli = df[df.region == 'eliminated'].dropna(subset=['GC3s', metric])
    if len(eli) < min_matched or len(ret) < min_matched:
        return None

    order = np.argsort(ret['GC3s'].to_numpy())
    r_gc = ret['GC3s'].to_numpy()[order]
    r_val = ret[metric].to_numpy()[order]
    e_gc = eli['GC3s'].to_numpy()
    e_val = eli[metric].to_numpy()

    p_metric, p_gc3s, diffs, n_matched = [], [], [], []
    for _ in range(reps):
        used = np.zeros(len(r_gc), dtype=bool)
        pick_r, pick_e = [], []
        for i in rng.permutation(len(e_gc)):
            lo = np.searchsorted(r_gc, e_gc[i] - caliper, 'left')
            hi = np.searchsorted(r_gc, e_gc[i] + caliper, 'right')
            if hi <= lo:
                continue
            free = np.flatnonzero(~used[lo:hi])
            if not len(free):
                continue
            j = lo + free[rng.integers(len(free))]
            used[j] = True
            pick_r.append(j)
            pick_e.append(i)
        if len(pick_r) < min_matched:
            continue
        rv, ev = r_val[pick_r], e_val[pick_e]
        _, p = stats.mannwhitneyu(rv, ev, alternative='two-sided')
        _, pg = stats.mannwhitneyu(r_gc[pick_r], e_gc[pick_e],
                                   alternative='two-sided')
        p_metric.append(p)
        p_gc3s.append(pg)
        diffs.append(float(np.median(ev) - np.median(rv)))
        n_matched.append(len(pick_r))
    if not p_metric:
        return None
    p_metric = np.asarray(p_metric)
    return {'median_p': float(np.median(p_metric)),
            'p_lo': float(np.percentile(p_metric, 2.5)),
            'p_hi': float(np.percentile(p_metric, 97.5)),
            'frac_sig': float((p_metric < 0.05).mean()),
            'median_diff': float(np.median(diffs)),
            'diff_lo': float(np.percentile(diffs, 2.5)),
            'diff_hi': float(np.percentile(diffs, 97.5)),
            'median_gc3s_p': float(np.median(p_gc3s)),
            'n_matched': int(np.median(n_matched)),
            'reps': len(p_metric)}


def centroid_ci(x, y, rng, n_boot=2000, level=95.0):
    """Median centroid of a point cloud with a percentile bootstrap interval on
    each axis independently.

    The intervals describe how well each coordinate is located, which is what a
    reader needs to judge whether a centroid has been placed meaningfully. They
    are not a test: ENC depends on GC3s, so the two axes are not independent and
    a rectangle drawn from them is not a joint confidence region. The comparison
    that is licensed here is the vertical offset from the expectation curve --
    that is dENC, and it is tested separately on the per-gene values.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return None
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    bx = np.median(x[idx], axis=1)
    by = np.median(y[idx], axis=1)
    lo, hi = (100 - level) / 2, 100 - (100 - level) / 2
    return {
        'gc3s': float(np.median(x)), 'enc': float(np.median(y)), 'n': int(len(x)),
        'gc3s_lo': float(np.percentile(bx, lo)), 'gc3s_hi': float(np.percentile(bx, hi)),
        'enc_lo': float(np.percentile(by, lo)), 'enc_hi': float(np.percentile(by, hi)),
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def save_figure(fig, outdir, stem):
    for problem in audit(fig, verbose=False):
        logging.warning(f'FIGURE AUDIT [{stem}]: {problem}')
    save(fig, outdir, stem, formats=[e.lstrip('.') for e in EXTS])
    plt.close(fig)
    logging.info(f"Saved: {stem}{'/'.join(EXTS)}")


def box_stats(values, whis=1.5):
    """Box and whisker statistics for one group, in the form ``ax.bxp`` draws.

    Computed from every value passed, never from the values drawn. Median and
    quartiles are ``np.percentile`` with linear interpolation, which is what
    pandas ``median()`` and ``quantile()`` give in the summary table, so the
    box on the page is the box in the table. Whiskers follow Tukey,
    matplotlib's default rule: each ends at the most extreme value within
    ``whis`` x IQR of the box, or at the box edge if no value lies beyond it.
    The same function as gene_structure.py (Figure 2D) and
    dnds_decomposition.py (Figure 2A), without the log-axis option, which these
    linear axes do not need.
    """
    v = np.asarray(values, dtype=float)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    lo_lim = q1 - whis * (q3 - q1)
    hi_lim = q3 + whis * (q3 - q1)
    up = v[v <= hi_lim]
    whishi = float(up.max()) if len(up) and up.max() >= q3 else float(q3)
    down = v[v >= lo_lim]
    whislo = float(down.min()) if len(down) and down.min() <= q1 else float(q1)
    return {'med': float(med), 'q1': float(q1), 'q3': float(q3),
            'whislo': whislo, 'whishi': whishi,
            # Read unconditionally by bxp() in older matplotlib releases.
            'fliers': np.array([]),
            'beyond_whiskers': int(((v < whislo) | (v > whishi)).sum())}


def centroid_legend(show_ci):
    """Handles for the shared centroid legend.

    The centroid entry is a *pair* of markers, blue then red, so the key uses
    the same two colours as the panels. A single grey stand-in marker is smaller
    but it invents a third colour that appears nowhere in the figure, which
    makes the reader look for grey data. Both groups' centroids are filled on
    every sub-panel, so one entry covers them: colour identifies the group.
    """
    centroids = tuple(
        Line2D([], [], marker='o', linestyle='none', markersize=3.4,
               markerfacecolor=REGION[k], markeredgecolor=NEUTRAL['ink'],
               markeredgewidth=0.7)
        for k in ('retained', 'eliminated'))

    curve = Line2D([], [], linestyle='--', color=NEUTRAL['ink'],
                   linewidth=LINEWIDTH['rule'])
    handles = [curve, centroids]
    labels = ['ENC expected from GC3s',
              'Median centroid, 95% CI' if show_ci else 'Median centroid']
    hmap = {tuple: HandlerTuple(ndivide=None, pad=0.4)}
    return handles, labels, hmap


def _boxplot_axes(ax, df, metric, ylabel, rng, max_points=1200, cap_pct=99.5,
                  box_alpha=0.40, ylim=None, bracket_label='stars'):
    """Retained versus eliminated for one metric. Returns the test result and
    the accounting of what is drawn.

    The test, the box and the whiskers use every gene (box_stats). The drawing
    range is the central ``cap_pct`` percent of the pooled values, widened to
    take in every whisker end. ``ylim`` fixes the y-range (shared across
    species); a whisker end outside it widens it for this panel, with a
    warning, so no whisker is cut. Genes outside the drawing range are not
    drawn, and are counted in ``res['n_not_drawn']``.
    """
    ret = df.loc[df['region'] == 'retained', metric].dropna()
    eli = df.loc[df['region'] == 'eliminated', metric].dropna()
    if len(ret) < 5 or len(eli) < 5:
        ax.text(0.5, 0.5, 'insufficient data', transform=ax.transAxes,
                ha='center', va='center', fontsize=FONT['annotation'])
        ax.set_axis_off()
        return None

    _, p = stats.mannwhitneyu(ret, eli, alternative='two-sided')
    res = dict(n_ret=len(ret), n_elim=len(eli), med_ret=float(ret.median()),
               med_elim=float(eli.median()), p=float(p))

    # Box and whiskers from every gene (box_stats), never from the genes drawn:
    # ax.boxplot() computed them from the trimmed values below.
    box = {key: box_stats(vals.to_numpy(dtype=float))
           for key, vals in zip(REGIONS, (ret, eli))}
    wlo = min(b['whislo'] for b in box.values())
    whi = max(b['whishi'] for b in box.values())
    res['box'] = {key: {k: b[k] for k in ('q1', 'med', 'q3', 'whislo', 'whishi',
                                          'beyond_whiskers')}
                  for key, b in box.items()}

    # Drawing range: the central mass, widened to every whisker end, then the
    # fixed range if given. The fixed range is itself widened if a whisker end
    # lies outside it: a cut whisker misreports the spread, while a wider axis
    # only costs a shared scale, and the warning says which bound to change.
    lo, hi = np.percentile(pd.concat([ret, eli]), [100 - cap_pct, cap_pct])
    lo, hi = min(lo, wlo), max(hi, whi)
    if ylim is not None:
        ylim = [float(ylim[0]), float(ylim[1])]
        if wlo < ylim[0]:
            ylim[0] = wlo - (ylim[1] - ylim[0]) * 0.03
            logging.warning('  %s: --denc-lim lower bound lowered to %.2f so the '
                            'lower whisker end (%.2f) is not cut; give every '
                            'species the same --denc-lim to keep one scale',
                            metric, ylim[0], wlo)
        if whi > ylim[1]:
            ylim[1] = whi               # the bracket headroom is added below
            logging.warning('  %s: --denc-lim upper bound raised to %.2f so the '
                            'upper whisker end is not cut; give every species '
                            'the same --denc-lim to keep one scale', metric, whi)
        lo, hi = max(lo, ylim[0]), min(hi, ylim[1])
    shown = [ret[(ret >= lo) & (ret <= hi)], eli[(eli >= lo) & (eli <= hi)]]
    res['n_not_drawn'] = {'retained': int(len(ret) - len(shown[0])),
                          'eliminated': int(len(eli) - len(shown[1]))}
    res['draw_range'] = (float(lo), float(hi))

    # The box sits *under* the jitter. Drawn over it at high opacity the fill
    # washes out exactly the points the strip plot exists to show -- the ones in
    # the interquartile range, which is most of them. Region colour still reads
    # at 0.40 because the box is a large solid area next to its partner.
    # bxp() draws the precomputed statistics; boxplot() would recompute them
    # from whatever it is handed.
    bp = ax.bxp([{k: box[key][k] for k in BXP_KEYS} for key in REGIONS],
                positions=[0, 1], widths=0.5, patch_artist=True,
                showfliers=False, zorder=2)
    for box_patch, key in zip(bp['boxes'], REGIONS):
        # Alpha goes in the face colour, NOT via set_alpha(). Patch-level alpha
        # dims every component of the artist, so set_alpha(0.4) also fades the
        # edge from black to mid-grey -- the outline stops reading as an outline
        # and the box loses its boundary at print size.
        box_patch.set_facecolor(to_rgba(REGION[key], box_alpha))
        box_patch.set_alpha(None)
        box_patch.set_edgecolor('black')
        box_patch.set_linewidth(LINEWIDTH['thin'])
    # Whiskers, caps and the median go back on top, so the summary statistics
    # stay legible through the cloud.
    for part, z in (('whiskers', 4), ('caps', 4), ('medians', 5)):
        for artist in bp[part]:
            artist.set_linewidth(LINEWIDTH['data'] if part == 'medians'
                                 else LINEWIDTH['thin'])
            artist.set_color('black')
            artist.set_zorder(z)

    for i, (vals, key) in enumerate(zip(shown, REGIONS)):
        v = vals.to_numpy()
        if len(v) > max_points:
            v = v[rng.choice(len(v), max_points, replace=False)]
        ax.scatter(np.full(len(v), i) + rng.normal(0, 0.06, len(v)), v,
                   s=1.2, alpha=0.28, color=REGION[key], edgecolors='none',
                   zorder=3, rasterized=True)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Retained\n(n = {len(ret):,})',
                        f'Eliminated\n(n = {len(eli):,})'])
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.6, 1.6)

    # Bracket above the highest drawn point, never through the cloud. Both
    # whisker ends are genes' values inside the drawing range, so the bracket
    # also clears the whiskers.
    top = max(s.max() for s in shown)
    bottom = min(s.min() for s in shown)
    if ylim is None:
        span = top - bottom
        y = top + span * 0.06
        y0, y1 = bottom - span * 0.05, y + span * 0.16
    else:
        span = ylim[1] - ylim[0]
        y = top + span * 0.04
        y0, y1 = ylim
        if y + span * 0.13 > y1:            # bracket and label must fit
            y1 = y + span * 0.13
            logging.warning('  %s: --denc-lim upper bound raised to %.2f so the '
                            'bracket clears the data', metric, y1)
    ax.set_ylim(y0, y1)
    res['axis'] = (float(y0), float(y1))
    ax.plot([0, 0, 1, 1], [y, y + span * 0.03, y + span * 0.03, y],
            color=NEUTRAL['ink'], linewidth=LINEWIDTH['thin'],
            solid_joinstyle='miter')
    ax.text(0.5, y + span * 0.045, stars(p) if bracket_label == 'stars' else fmt_p(p),
            ha='center', va='bottom', fontsize=FONT['annotation'],
            color=NEUTRAL['ink'])
    return res


def _enc_gc3_axes(ax, df, focus, centroids, xlim, ylim, show_ylabel,
                  show_ci=False, denc_median=None):
    """One ENC-versus-GC3 sub-panel, drawn for ``focus`` with both centroids.

    The focus group's points are shown; both groups' centroids are drawn on
    every sub-panel, filled in their region colour, so the two can be compared
    without switching between panels. The corner value is the focus group's
    median dENC over genes (``denc_median``), the number the boxplot and the
    text use. Returns how many of the focus group's genes fall outside the axis
    limits.
    """
    sub = df[df['region'] == focus].dropna(subset=['GC3s', 'ENC'])
    n_out = int((~(sub['GC3s'].between(*xlim) & sub['ENC'].between(*ylim))).sum())
    ax.scatter(sub['GC3s'], sub['ENC'], s=1.2, alpha=0.22, color=REGION[focus],
               edgecolors='none', zorder=2, rasterized=True)

    grid = np.linspace(max(xlim[0], 0.01), min(xlim[1], 0.99), 300)
    ax.plot(grid, enc_expected(grid), linestyle='--', color=NEUTRAL['ink'],
            linewidth=LINEWIDTH['rule'], zorder=4)

    for key in ('retained', 'eliminated'):
        cen = centroids.get(key)
        if cen is None:
            continue
        own = key == focus
        # Both groups filled in their region colour: colour identifies the
        # group. (Until 2026-10-26 the other group was an open ring and a dotted
        # drop marked each centroid's offset from the curve; see Changes.)
        ax.plot([cen['gc3s']], [cen['enc']], marker='o', markersize=3.4,
                markerfacecolor=REGION[key], markeredgecolor=NEUTRAL['ink'],
                markeredgewidth=0.7, zorder=6)
        # Interval on top of the marker, not behind it. With n = 9,085 the
        # retained interval is around 0.3 mm wide at print size -- narrower than
        # the marker -- so drawn underneath it is invisible, and an invisible
        # interval that the legend claims is worse than no interval at all.
        if show_ci:
            ax.errorbar(cen['gc3s'], cen['enc'],
                        xerr=[[cen['gc3s'] - cen['gc3s_lo']], [cen['gc3s_hi'] - cen['gc3s']]],
                        yerr=[[cen['enc'] - cen['enc_lo']], [cen['enc_hi'] - cen['enc']]],
                        fmt='none', ecolor=NEUTRAL['ink'],
                        elinewidth=LINEWIDTH['hairline'], capsize=1.4,
                        capthick=LINEWIDTH['hairline'], zorder=8)
        if own and denc_median is not None and focus in denc_median:
            # Corner, not adjacent to the marker. Beside the centroid the label
            # lands in the densest part of the cloud and needs a white box to
            # stay readable, which then hides data. In the corner it is clean,
            # and putting the retained label top-left and the eliminated label
            # top-right follows the house left/right region convention, so which
            # group it belongs to is legible without reading it. The value is
            # the median over genes, as in the boxplot and the text.
            corner = 'left' if focus == 'retained' else 'right'
            ax.text(0.03 if corner == 'left' else 0.97, 0.97,
                    f'median ΔENC = {denc_median[focus]:.2f}'.replace('-', '\u2212'),
                    transform=ax.transAxes, ha=corner, va='top',
                    fontsize=FONT['annotation'], color=NEUTRAL['ink'],
                    zorder=9)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel('GC3s')
    if show_ylabel:
        ax.set_ylabel('ENC')
    else:
        ax.tick_params(labelleft=False)
    ax.set_title(f'{focus.capitalize()} (n = {len(sub):,})',
                 fontsize=FONT['panel_title'])
    return n_out


def _limits(df, limits):
    """(xlim, ylim) for the ENC-GC3s panels: the shared limits, or from the data.

    The data-driven version (--auto-limits) takes the lower ENC bound from the
    0.1% quantile of the drawn points, not their minimum: the minimum can belong
    to a gene outside the GC3s range, which then sets the axis without being
    visible.
    """
    if limits.get('gc3s') and limits.get('enc'):
        return tuple(limits['gc3s']), tuple(limits['enc'])
    pts = df.dropna(subset=['GC3s', 'ENC'])
    xlim = (max(0.0, pts['GC3s'].quantile(0.001) - 0.03),
            min(1.0, pts['GC3s'].quantile(0.999) + 0.03))
    inx = pts[pts['GC3s'].between(*xlim)]
    # Headroom above the ENC ceiling so the corner label lands on blank canvas.
    # Genes pinned at 61 form a dense row along the top of the retained panel,
    # and 1.5 units is not enough clearance for a 7 pt label.
    ylim = (max(ENC_MIN - 1, float(inx['ENC'].quantile(0.001)) - 2), ENC_MAX + 3.5)
    return xlim, ylim


def plot_panel_strip(df, centroids, outdir, stem, rng, panel_b='cai',
                     width='double', height_mm=60.0, letters=False,
                     box_alpha=0.40, show_ci=False, limits=None,
                     denc_median=None, bracket_label='stars'):
    """Figure 2B-C: the two ENC-GC3 sub-panels (B), then one metric boxplot (C).

    Laid out with explicit margins rather than by a layout solver, so the shared
    legend below cannot renegotiate the size of the panels. The scatter pair
    comes first because the boxplot summarises it (median dENC per group).
    """
    metric, ylabel = {
        'cai': ('CAI', 'Codon adaptation index'),
        'denc': ('dENC', 'ΔENC (expected − observed)'),
        'enc': ('ENC', 'Effective number of codons'),
    }[panel_b]

    w_mm, h_mm = COL[width] * 25.4, height_mm
    fig = plt.figure(figsize=(COL[width], mm(height_mm)))
    fig.set_layout_engine('none')

    left, right, bottom, top, gutter = 13.0, 3.0, 16.5, 5.0, 13.0
    usable = w_mm - left - right - 2 * gutter
    w_box, w_sc = usable * 0.26, usable * 0.37
    y0, hgt = bottom / h_mm, (h_mm - bottom - top) / h_mm

    # Left to right: retained scatter, eliminated scatter, boxplot. Same widths
    # and gutters as the earlier box-first strip; the gutter before the boxplot
    # holds its tick labels and y-label, as the left margin did.
    ax_ret = fig.add_axes([left / w_mm, y0, w_sc / w_mm, hgt])
    ax_eli = fig.add_axes([(left + w_sc + gutter) / w_mm, y0, w_sc / w_mm, hgt])
    ax_box = fig.add_axes([(left + w_sc + gutter + w_sc + gutter) / w_mm, y0,
                           w_box / w_mm, hgt])

    limits = limits or {}
    res = _boxplot_axes(ax_box, df, metric, ylabel, rng, box_alpha=box_alpha,
                        ylim=limits.get('denc') if metric == 'dENC' else None,
                        bracket_label=bracket_label)

    xlim, ylim = _limits(df, limits)
    out = {'retained': _enc_gc3_axes(ax_ret, df, 'retained', centroids, xlim, ylim,
                                     True, show_ci, denc_median),
           'eliminated': _enc_gc3_axes(ax_eli, df, 'eliminated', centroids, xlim,
                                       ylim, False, show_ci, denc_median)}
    if res is not None:
        res['enc_gc3s_outside'] = out
        res['xlim'], res['ylim'] = xlim, ylim

    # Legend centred under the scatter pair it explains.
    pair_centre = (left + w_sc + gutter / 2) / w_mm
    handles, labels, hmap = centroid_legend(show_ci)
    fig.legend(handles=handles, labels=labels, handler_map=hmap,
               loc='lower center', bbox_to_anchor=(pair_centre, 0.015), ncol=2,
               frameon=False, fontsize=FONT['legend'], handlelength=1.7,
               columnspacing=1.4, handletextpad=0.5)

    if letters:
        panel_letter(ax_ret, 'B')
        panel_letter(ax_box, 'C')
    save_figure(fig, outdir, stem)
    return res


def plot_single_boxplot(df, metric, ylabel, outdir, stem, rng,
                        width='single', height_mm=62.0, box_alpha=0.40,
                        limits=None, bracket_label='stars'):
    w_mm, h_mm = COL[width] * 25.4, height_mm
    fig = plt.figure(figsize=(COL[width], mm(height_mm)))
    fig.set_layout_engine('none')
    left, right, bottom, top = 16.0, 3.0, 16.0, 4.0
    ax = fig.add_axes([left / w_mm, bottom / h_mm,
                       (w_mm - left - right) / w_mm, (h_mm - bottom - top) / h_mm])
    limits = limits or {}
    res = _boxplot_axes(ax, df, metric, ylabel, rng, box_alpha=box_alpha,
                        ylim=limits.get('denc') if metric == 'dENC' else None,
                        bracket_label=bracket_label)
    save_figure(fig, outdir, stem)
    return res


def plot_enc_gc3(df, centroids, outdir, stem, width='double', height_mm=62.0,
                 show_ci=False, limits=None, denc_median=None):
    w_mm, h_mm = COL[width] * 25.4, height_mm
    fig = plt.figure(figsize=(COL[width], mm(height_mm)))
    fig.set_layout_engine('none')
    left, right, bottom, top, gutter = 14.0, 3.0, 16.5, 5.0, 10.0
    each = (w_mm - left - right - gutter) / 2
    y0, hgt = bottom / h_mm, (h_mm - bottom - top) / h_mm
    ax1 = fig.add_axes([left / w_mm, y0, each / w_mm, hgt])
    ax2 = fig.add_axes([(left + each + gutter) / w_mm, y0, each / w_mm, hgt])

    xlim, ylim = _limits(df, limits or {})
    _enc_gc3_axes(ax1, df, 'retained', centroids, xlim, ylim, True, show_ci,
                  denc_median)
    _enc_gc3_axes(ax2, df, 'eliminated', centroids, xlim, ylim, False, show_ci,
                  denc_median)

    handles, labels, hmap = centroid_legend(show_ci)
    fig.legend(handles=handles, labels=labels, handler_map=hmap,
               loc='lower center', bbox_to_anchor=(0.5, 0.015), ncol=2,
               frameon=False, fontsize=FONT['legend'], handlelength=1.7,
               columnspacing=1.4, handletextpad=0.5)
    save_figure(fig, outdir, stem)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def fingerprint(path):
    st = os.stat(path)
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return {'path': os.path.abspath(path), 'bytes': st.st_size,
            'mtime_utc': datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
            'sha256_16': h.hexdigest()[:16]}


def git_commit():
    try:
        out = subprocess.run(['git', '-C', _HERE, 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            dirty = subprocess.run(['git', '-C', _HERE, 'status', '--porcelain'],
                                   capture_output=True, text=True, timeout=5)
            return out.stdout.strip() + ('+dirty' if dirty.stdout.strip() else '')
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Codon usage bias (ENC, GC3s, dENC, CAI) for retained '
                    'versus eliminated genes.')
    p.add_argument('--cds-fasta', required=True,
                   help='CDS FASTA, nucleotide and in-frame')
    p.add_argument('--rates-table', required=True,
                   help='Table carrying gene IDs and per-gene region calls')
    p.add_argument('--gene-col', default='gene1',
                   help='Gene ID column in the table (default: gene1)')
    p.add_argument('--region-col', default='region1',
                   help='Region column in the table (default: region1)')
    p.add_argument('--ref-genes', default=None,
                   help='One gene ID per line, the CAI reference set. Build it '
                        'from expression with make_cai_reference.py. Without '
                        'this the fallback is the most codon-biased retained '
                        'genes, which makes the CAI contrast partly circular.')
    p.add_argument('--min-codons', type=int, default=100,
                   help='Minimum informative codons per gene (default: 100)')
    p.add_argument('--panel-b', choices=['cai', 'denc', 'enc'], default='cai',
                   help="Metric for the boxplot in the Figure 2 strip "
                        "(default: cai, matching the frozen figure map; 'denc' "
                        'is the stronger choice when no expression-based '
                        'reference set is available)')
    p.add_argument('--n-boot', type=int, default=2000,
                   help='Bootstrap replicates for the centroid intervals')
    p.add_argument('--match-reps', type=int, default=500,
                   help='Replicates of the GC3s-matched comparison. 0 skips '
                        'the whole GC3s-controlled block.')
    p.add_argument('--match-caliper', type=float, default=0.01,
                   help='Maximum GC3s difference within a matched pair '
                        '(default: 0.01)')
    p.add_argument('--box-alpha', type=float, default=0.40,
                   help='Opacity of the box fill. The jitter is drawn on top of '
                        'the box, so a high value hides the points inside the '
                        'interquartile range (default: 0.40)')
    p.add_argument('--centroid-ci', action='store_true',
                   help='Draw bootstrap interval bars on the centroids. Off by '
                        'default: at n = 9,085 the retained interval is about '
                        '0.2 mm at print size, narrower than the marker itself, '
                        'so the bars read as clutter rather than as uncertainty. '
                        'The intervals are in the summary either way, and the '
                        'legend only claims a CI when this is set.')
    p.add_argument('--gc3s-lim', type=float, nargs=2, default=[0.20, 0.70],
                   metavar=('LO', 'HI'),
                   help='GC3s axis of the ENC-GC3s panels, shared by every species '
                        '(default: 0.20 0.70)')
    p.add_argument('--enc-lim', type=float, nargs=2, default=[36.0, 64.5],
                   metavar=('LO', 'HI'),
                   help='ENC axis of the ENC-GC3s panels; the top leaves room for '
                        'the corner label above the ceiling of 61 (default: 36 64.5)')
    p.add_argument('--denc-lim', type=float, nargs=2, default=[-6.5, 16.5],
                   metavar=('LO', 'HI'),
                   help='dENC axis of the boxplot, shared by every species '
                        '(default: -6.5 16.5; 16 until 2026-10-14, which the '
                        'bracket raised to 16.09 and 16.10 for P. univalens and '
                        'T. canis). Widened, with a warning, if a whisker end or '
                        'the bracket would not fit')
    p.add_argument('--auto-limits', action='store_true',
                   help='Set every axis from this species\' data instead of the '
                        'shared limits')
    p.add_argument('--bracket-label', choices=['stars', 'p'], default='stars',
                   help='Boxplot bracket label: significance markers, as in the '
                        'rest of Figure 2 (default), or the printed P')
    p.add_argument('--min-per-stratum', type=int, default=5,
                   help='Minimum genes of each region in a stratum for it to '
                        'contribute to the stratified test (default: 5)')
    p.add_argument('--seed', type=int, default=0,
                   help='RNG seed for jitter, subsampling and the bootstrap. '
                        'Fixed by default so a rerun reproduces the figure.')
    p.add_argument('--panel-letters', action='store_true',
                   help='Draw panel letters. Leave off if they are added during '
                        'figure assembly.')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'codon_usage_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    # mode='w': FileHandler appends by default, so rerunning into an
    # existing directory concatenates logs while the figures and
    # provenance.json are overwritten -- the log then describes more
    # runs than the directory contains.
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(os.path.join(args.outdir, 'codon_usage.log'),
                                     mode='w'),
                  logging.StreamHandler(sys.stdout)])
    for noisy in ('matplotlib.font_manager', 'fontTools', 'PIL'):
        logging.getLogger(noisy).setLevel(logging.ERROR)

    logging.info('=' * 60)
    logging.info('Codon usage bias')
    logging.info('=' * 60)

    for path in (args.cds_fasta, args.rates_table):
        if not os.path.isfile(path):
            sys.exit(f'ERROR: input not found: {path}')

    rng = np.random.default_rng(args.seed)
    sequences = parse_fasta(args.cds_fasta)

    table = pd.read_csv(args.rates_table, sep='\t')
    for col in (args.gene_col, args.region_col):
        if col not in table.columns:
            sys.exit(f"ERROR: no '{col}' column in {args.rates_table}. "
                     f'Columns present: {table.columns.tolist()}')

    region_map = dict(zip(table[args.gene_col], table[args.region_col]))
    counts = pd.Series(list(region_map.values())).value_counts()
    logging.info('Region calls in %s: %s', args.region_col,
                 ', '.join(f'{k} {v:,}' for k, v in counts.items()))

    # --- per-gene metrics ---------------------------------------------------
    rows, per_gene_counts, n_clamped = [], {}, Counter()
    for gene_id, seq in sequences.items():
        region = region_map.get(gene_id)
        if region not in ('retained', 'eliminated'):
            continue
        if len(seq) % 3:
            seq = seq[:len(seq) - len(seq) % 3]
        cc = count_codons(seq)
        total = sum(cc.values())
        if total < args.min_codons:
            continue
        enc, clamped = compute_enc(cc)
        if enc is None:
            continue
        if clamped:
            n_clamped[region] += 1
        per_gene_counts[gene_id] = cc
        gc3s = compute_gc3(cc, synonymous_only=True)
        rows.append({'gene_id': gene_id, 'region': region, 'ENC': enc,
                     'GC3s': gc3s, 'GC3_all': compute_gc3(cc, synonymous_only=False),
                     'ENC_at_ceiling': clamped, 'total_codons': total})

    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit('ERROR: no genes passed filters. Check that the FASTA headers '
                 f"and the '{args.gene_col}' column use the same ID convention.")
    logging.info('Metrics for %s genes (%s retained, %s eliminated)',
                 f'{len(df):,}', f'{(df.region == "retained").sum():,}',
                 f'{(df.region == "eliminated").sum():,}')
    for region in ('retained', 'eliminated'):
        n = int((df.region == region).sum())
        if n and n_clamped[region]:
            logging.info('  %s at the ENC ceiling of 61: %d of %d (%.1f%%)',
                         region, n_clamped[region], n, 100 * n_clamped[region] / n)

    df['dENC'] = enc_expected(df['GC3s']) - df['ENC']

    # --- CAI reference set --------------------------------------------------
    if args.ref_genes:
        with open(args.ref_genes) as fh:
            ref_ids = {ln.strip() for ln in fh if ln.strip()}
        ref_source = f'expression-based list ({os.path.basename(args.ref_genes)})'
        circular = False
        logging.info('CAI reference: %s IDs from %s', f'{len(ref_ids):,}', args.ref_genes)
    else:
        retained = df[df.region == 'retained']
        cutoff = retained['ENC'].quantile(0.05)
        ref_ids = set(retained.loc[retained['ENC'] <= cutoff, 'gene_id'])
        ref_source = 'fallback: most codon-biased 5% of retained genes'
        circular = True
        logging.warning('No --ref-genes given. Falling back to the most '
                        'codon-biased 5%% of retained genes (ENC <= %.1f, n = %d). '
                        'This selects the reference by the very property CAI '
                        'measures, so the retained-vs-eliminated CAI contrast is '
                        'partly circular -- lead with dENC, which is '
                        'reference-free.', cutoff, len(ref_ids))

    pooled = Counter()
    n_ref = 0
    for gid in ref_ids:
        if gid in per_gene_counts:
            pooled += per_gene_counts[gid]
            n_ref += 1
    if n_ref == 0:
        # Show both ID spaces. An annotation workbook usually carries several
        # naming conventions, and the only thing that identifies the mismatch is
        # seeing the two side by side.
        theirs = sorted(ref_ids)[:3]
        ours = sorted(per_gene_counts)[:3]
        sys.exit(
            'ERROR: none of the reference gene IDs appear in the analysed set.\n'
            f'  --ref-genes IDs look like : {", ".join(theirs)}\n'
            f'  analysed gene IDs look like: {", ".join(ours)}\n'
            '  These are different ID conventions. Rebuild the reference with\n'
            '  make_cai_reference.py --probe --target-table <the rates table>\n'
            '  to see which expression column uses the analysed convention.')
    if n_ref < 30:
        logging.warning('Only %d reference genes matched. A CAI w table built '
                        'on that few genes is unstable; widen the reference or '
                        'check the ID convention.', n_ref)
    logging.info('CAI reference: %s genes pooled, %s codons',
                 f'{n_ref:,}', f'{sum(pooled.values()):,}')

    w_table = compute_w_table(pooled)
    df['CAI'] = [compute_cai(per_gene_counts[g], w_table) for g in df['gene_id']]

    # --- centroids ----------------------------------------------------------
    centroids = {}
    for region in ('retained', 'eliminated'):
        sub = df[df.region == region]
        cen = centroid_ci(sub['GC3s'], sub['ENC'], rng, n_boot=args.n_boot)
        if cen:
            cen['enc_expected'] = float(enc_expected(cen['gc3s']))
            cen['dENC_at_centroid'] = cen['enc_expected'] - cen['enc']
            centroids[region] = cen
            logging.info('%s centroid: GC3s %.4f [%.4f-%.4f], ENC %.2f [%.2f-%.2f], '
                         'dENC at centroid %.2f', region, cen['gc3s'],
                         cen['gc3s_lo'], cen['gc3s_hi'], cen['enc'],
                         cen['enc_lo'], cen['enc_hi'], cen['dENC_at_centroid'])

    # --- summary ------------------------------------------------------------
    METRICS = [('ENC', 'Effective number of codons'),
               ('GC3s', 'GC3s'),
               ('dENC', 'ΔENC (expected − observed)'),
               ('CAI', 'Codon adaptation index')]
    summary_rows, tests = [], {}
    for metric, _ in METRICS:
        ret = df.loc[df.region == 'retained', metric].dropna()
        eli = df.loc[df.region == 'eliminated', metric].dropna()
        for region, vals in (('retained', ret), ('eliminated', eli)):
            summary_rows.append({'metric': metric, 'region': region, 'n': len(vals),
                                 'median': vals.median(), 'mean': vals.mean(),
                                 'Q1': vals.quantile(0.25), 'Q3': vals.quantile(0.75)})
        if len(ret) and len(eli):
            u, p = stats.mannwhitneyu(ret, eli, alternative='two-sided')
            tests[metric] = {'U': float(u), 'p': float(p),
                             'med_ret': float(ret.median()), 'med_elim': float(eli.median())}
            summary_rows.append({'metric': metric, 'region': 'MWU', 'n': None,
                                 'median': p, 'mean': u, 'Q1': None, 'Q3': None})

    pd.DataFrame(summary_rows).to_csv(
        os.path.join(args.outdir, 'codon_usage_summary.tsv'), sep='\t', index=False)
    df.to_csv(os.path.join(args.outdir, 'codon_usage_per_gene.tsv'),
              sep='\t', index=False)

    lines = ['Codon usage bias', '=' * 60, '',
             f'CDS         {args.cds_fasta}',
             f'Regions     {args.region_col} in {args.rates_table}',
             f'CAI ref     {ref_source}', '',
             f"{'metric':<8}{'retained':>11}{'eliminated':>12}{'P':>15}", '']
    for metric, _ in METRICS:
        t = tests.get(metric)
        if t:
            lines.append(f"  {metric:<6}{t['med_ret']:>11.4f}{t['med_elim']:>12.4f}"
                         f'{fmt_p_plain(t["p"]):>15}')
    # A gene at the ENC ceiling has its dENC computed against a clamped ENC, so
    # its dENC is pushed down. Eliminated genes hit the ceiling more often, so
    # the ceiling could manufacture part of the dENC difference on its own.
    # Repeating the test without those genes is the first thing a referee will
    # ask for, and it costs nothing to answer here.
    clear = df[~df['ENC_at_ceiling'].astype(bool)]
    lines += ['', 'Sensitivity: genes at the ENC ceiling excluded']
    for metric in ('ENC', 'dENC'):
        r = clear.loc[clear.region == 'retained', metric].dropna()
        e = clear.loc[clear.region == 'eliminated', metric].dropna()
        if len(r) > 5 and len(e) > 5:
            _, p = stats.mannwhitneyu(r, e, alternative='two-sided')
            lines.append(f'  {metric:<6}{r.median():>11.4f}{e.median():>12.4f}'
                         f'{fmt_p_plain(p):>15}   '
                         f'(n = {len(r):,} / {len(e):,})')
            logging.info('Ceiling-excluded %s: retained %.4f, eliminated %.4f, '
                         'P = %s', metric, r.median(), e.median(), fmt_p_plain(p))

    # --- GC3s-controlled comparison -----------------------------------------
    # CAI and ENC both move with GC3s, and the two groups differ in GC3s, so a
    # raw contrast in either cannot separate translational selection from base
    # composition. Two independent controls: matched sampling (discards data,
    # assumes nothing) and a stratified rank test (uses everything, assumes one
    # direction across strata). Agreement between them is the result; either
    # alone is weaker. dENC is included as a positive control -- it already
    # removes composition analytically, so matching should barely move it.
    matched, strat = {}, {}
    if args.match_reps > 0:
        # Fixed-width strata, one caliper wide -- NOT quantiles. Quantile strata
        # are wide where the data are dense, which is exactly where both groups
        # live, and the eliminated genes then sit systematically at the upper
        # edge of every stratum. On a synthetic set built with no region effect
        # beyond GC3s, decile strata returned P = 0.001 for CAI; caliper-width
        # strata return the correct null. Width is tied to the caliper so the
        # stratified test cannot be more permissive than the matched one.
        gc = df['GC3s'].dropna()
        width = max(args.match_caliper, 1e-4)
        edges = np.arange(np.floor(gc.min() / width) * width,
                          gc.max() + width, width)
        df['_gc3s_bin'] = pd.cut(df['GC3s'], bins=edges, include_lowest=True)
        logging.info('GC3s-controlled comparison: %d matched replicates, '
                     '%d fixed-width strata of %.3f GC3s, caliper %.3f',
                     args.match_reps, len(edges) - 1, width, args.match_caliper)
        if args.match_caliper > 0.01:
            logging.warning('Caliper %.3f is wider than the 0.01 this was '
                            'calibrated at. On synthetic null data the '
                            'stratified test returned P = 2 x 10^-4 at a 0.02 '
                            'caliper -- a false positive. Treat the stratified '
                            'result as unreliable and use the matched one.',
                            args.match_caliper)
        for metric in ('CAI', 'ENC', 'dENC'):
            m = gc3s_matched_test(df, metric, rng, caliper=args.match_caliper,
                                  reps=args.match_reps)
            if m:
                matched[metric] = m
                logging.info('  matched %-5s diff %+.4f [%+.4f, %+.4f], '
                             'median P = %s, %.0f%% of replicates P < 0.05 '
                             '(n = %d per group, GC3s check P = %.2f)',
                             metric, m['median_diff'], m['diff_lo'], m['diff_hi'],
                             fmt_p_plain(m['median_p']), 100 * m['frac_sig'],
                             m['n_matched'], m['median_gc3s_p'])
            v = van_elteren(df, metric, '_gc3s_bin',
                             min_per_stratum=args.min_per_stratum)
            if v:
                strat[metric] = v
                logging.info('  stratified %-5s Z = %+.2f, P = %s (%d strata)',
                             metric, v[0], fmt_p_plain(v[1]), v[2])

    if matched or strat:
        lines += ['', 'GC3s-controlled comparison',
                  '  Matched sampling: each eliminated gene paired with an unused',
                  f'  retained gene within {args.match_caliper:g} GC3s, repeated '
                  f'{args.match_reps} times.',
                  '',
                  f"  {'metric':<8}{'median diff':>13}{'95% range':>22}"
                  f"{'median P':>16}{'% reps P<.05':>14}{'GC3s check':>12}"]
        for metric in ('CAI', 'ENC', 'dENC'):
            m = matched.get(metric)
            if not m:
                continue
            rng_txt = f"[{m['diff_lo']:+.4f}, {m['diff_hi']:+.4f}]"
            pct_txt = f"{100 * m['frac_sig']:.0f}%"
            lines.append(
                f"  {metric:<8}{m['median_diff']:>+13.4f}{rng_txt:>22}"
                f"{fmt_p_plain(m['median_p']):>16}{pct_txt:>14}"
                f"{m['median_gc3s_p']:>12.2f}")
        if matched:
            # How much of each raw contrast is left once GC3s is equalised.
            # Under 100% means composition was carrying part of it; over 100%
            # means the raw comparison was understating the difference because
            # the groups' GC3s difference pushed the metric the other way.
            share = []
            for metric in ('CAI', 'ENC', 'dENC'):
                m, t = matched.get(metric), tests.get(metric)
                if m and t:
                    raw = t['med_elim'] - t['med_ret']
                    if abs(raw) > 1e-12:
                        share.append(f'{metric} {100 * m["median_diff"] / raw:.0f}%')
            lines += ['',
                      f"  Matched n = {list(matched.values())[0]['n_matched']} per group. "
                      'The GC3s check is the',
                      '  post-matching GC3s test: it has to be non-significant for',
                      '  anything else in this block to mean what it says.']
            if share:
                lines.append('  Share of the raw difference surviving matching: '
                             + ', '.join(share) + '.')
        if strat:
            lines += ['', '  Stratified rank test (van Elteren, GC3s strata), all genes:',
                      '']
            for metric in ('CAI', 'ENC', 'dENC'):
                v = strat.get(metric)
                if v:
                    lines.append(f'  {metric:<7}Z = {v[0]:+7.2f}   '
                                 f'P = {fmt_p_plain(v[1])}   ({v[2]} strata)')

    denc_median = {r: float(df.loc[df.region == r, 'dENC'].median())
                   for r in ('retained', 'eliminated')
                   if (df.region == r).any()}
    lines += ['', 'Centroids (median, 95% bootstrap CI)']
    for region, cen in centroids.items():
        lines.append(f"  {region:<11} GC3s {cen['gc3s']:.4f} "
                     f"[{cen['gc3s_lo']:.4f}-{cen['gc3s_hi']:.4f}]   "
                     f"ENC {cen['enc']:.2f} [{cen['enc_lo']:.2f}-{cen['enc_hi']:.2f}]")
        lines.append(f"  {'':<11} median dENC over genes {denc_median[region]:.2f} "
                     f"(printed in the panel corner); offset of the centroid from "
                     f"the curve {cen['dENC_at_centroid']:.2f} (not drawn)")
    lines += ['  The two differ because the expected-ENC curve is not linear: the',
              '  offset at the median GC3s is not the median of the per-gene offsets.',
              '  The text and the boxplot use the median over genes.']
    lines += ['', 'Reporting notes',
              '  dENC is the metric to lead with: ENC alone confounds',
              '  translational selection with GC3s composition, and the two',
              '  groups differ in GC3s.',
              '  GC3s excludes Met and Trp, which have no synonymous third',
              '  position. Wright\'s expected-ENC curve is defined on GC3s;',
              '  GC3 over all codons is also in the per-gene table, for',
              '  comparison with analyses that use it.',
              '  Centroid intervals are per-axis percentile bootstraps. They',
              '  locate each coordinate; they are not a joint region and not a',
              '  test. The licensed comparison is the vertical offset from the',
              '  expectation curve, which is dENC and is tested per gene.']
    if circular:
        lines += ['  CAI used the fallback reference set, so its',
                  '  retained-vs-eliminated contrast is partly circular.',
                  '  Supporting evidence only.']
    for region in ('retained', 'eliminated'):
        n = int((df.region == region).sum())
        if n and n_clamped[region]:
            lines.append(f'  {region}: {n_clamped[region]} of {n} genes '
                         f'({100 * n_clamped[region] / n:.1f}%) sit at the ENC '
                         f'ceiling of 61.')
    # --- figures ------------------------------------------------------------
    # Drawn after every statistic, so the jitter draws cannot touch a result.
    limits = ({} if args.auto_limits else
              {'gc3s': args.gc3s_lim, 'enc': args.enc_lim, 'denc': args.denc_lim})
    strip = plot_panel_strip(df, centroids, args.outdir, 'codon_usage_panel', rng,
                             panel_b=args.panel_b, letters=args.panel_letters,
                             box_alpha=args.box_alpha, show_ci=args.centroid_ci,
                             limits=limits, denc_median=denc_median,
                             bracket_label=args.bracket_label)
    plot_enc_gc3(df, centroids, args.outdir, 'enc_vs_gc3',
                 show_ci=args.centroid_ci, limits=limits, denc_median=denc_median)
    single = {}
    for metric, ylabel in METRICS:
        single[metric] = plot_single_boxplot(df, metric, ylabel, args.outdir,
                                             f'{metric.lower()}_boxplot', rng,
                                             box_alpha=args.box_alpha, limits=limits,
                                             bracket_label=args.bracket_label)

    # What the figures leave out, per panel and group: nothing is cut silently.
    lines += ['', 'Figures',
              f"  Axis limits: {'from this species' + chr(39) + ' data (--auto-limits)' if args.auto_limits else 'shared across species'}"]
    if strip:
        xl, yl = strip['xlim'], strip['ylim']
        o = strip['enc_gc3s_outside']
        lines.append(f'  ENC-GC3s panels: GC3s {xl[0]:.2f}-{xl[1]:.2f}, ENC {yl[0]:.1f}-'
                     f"{yl[1]:.1f}; genes outside: retained {o['retained']}, "
                     f"eliminated {o['eliminated']}")
    drawn = {m: r for m, r in single.items() if r}
    strip_metric = {'cai': 'CAI', 'denc': 'dENC', 'enc': 'ENC'}[args.panel_b]
    lines += ['  Boxes: quartiles of every gene; centre line, the median. Whiskers:',
              '  Tukey, to the most extreme gene within 1.5 x IQR of the box, as in',
              f'  Figure 2A and D. The strip\'s boxplot is the {strip_metric} row.',
              f"    {'metric':<8}{'group':<12}{'Q1':>9}{'median':>9}{'Q3':>9}"
              f"{'whisker lo':>12}{'whisker hi':>12}{'beyond':>8}"]
    for metric, r in drawn.items():
        for key in REGIONS:
            b = r['box'][key]
            lines.append(f"    {metric:<8}{key:<12}{b['q1']:>9.4f}{b['med']:>9.4f}"
                         f"{b['q3']:>9.4f}{b['whislo']:>12.4f}{b['whishi']:>12.4f}"
                         f"{b['beyond_whiskers']:>8,}")
    lines += ['  Boxplot axes, and genes not drawn (outside the central 99% of the',
              '  pooled values widened to the whisker ends, or outside the axis);',
              '  the boxes, whiskers and tests use every gene:']
    for metric, r in drawn.items():
        nd, (a0, a1) = r['n_not_drawn'], r['axis']
        lines.append(f'    {metric:<8}axis {a0:.4g} to {a1:.4g}; not drawn: retained '
                     f'{nd["retained"]}, eliminated {nd["eliminated"]}')
    lines.append(f'  Bracket label: {"significance markers" if args.bracket_label == "stars" else "printed P"}'
                 + (' (**** P < 0.0001, *** P < 0.001, ** P < 0.01, * P < 0.05)'
                    if args.bracket_label == 'stars' else ''))
    summary_txt = '\n'.join(lines) + '\n'
    with open(os.path.join(args.outdir, 'codon_usage_summary.txt'), 'w') as fh:
        fh.write(summary_txt)

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(),
            'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv,
            'args': vars(args),
            'inputs': {k: fingerprint(v) for k, v in
                       (('cds_fasta', args.cds_fasta),
                        ('rates_table', args.rates_table),
                        ('ref_genes', args.ref_genes)) if v},
            'cai_reference': {'source': ref_source, 'circular': circular,
                              'n_genes_pooled': n_ref,
                              'n_codons': int(sum(pooled.values()))},
            'genes_analysed': int(len(df)),
            'enc_at_ceiling': dict(n_clamped),
        }, fh, indent=2)

    print()
    print(summary_txt)
    logging.info('Summary: %s', os.path.join(args.outdir, 'codon_usage_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
