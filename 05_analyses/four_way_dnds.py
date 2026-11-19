#!/usr/bin/env python3
"""
four_way_dnds.py - four-taxon dS distance tree and pairwise dN/dS comparison.

Reads the six pairwise substitution-rate tables for Ascaris, Parascaris,
Baylisascaris and Toxocara and writes:

    ds_tree_four_way.{svg,pdf,png}   dS distance tree with gene-resampling support
    ds_overlay.{svg,pdf,png}         dS distributions, all six comparisons
    dnds_six_way_DIAGNOSTIC.*        pipeline check, not a manuscript figure
    four_way_summary.txt             every value, plus drafted legend text
    provenance.json                  inputs, arguments and versions for the run


Method
------
**Topology.** The four-point condition selects the unrooted topology from the
4x4 matrix of median dS. It is checked against the accepted ascaridid topology
(((Ascaris, Parascaris), Baylisascaris), Toxocara); a disagreement is reported
as a warning and branch lengths are still fitted to the accepted topology.

**Branch lengths.** An unrooted four-taxon tree has 2n-3 = 5 branches but six
pairwise distances, so the system is over-determined by one degree of freedom.
Both fits are always reported. `--fit exact` (default) solves five distances
exactly and leaves d(Parascaris, Toxocara) as the residual; `--fit ols` fits all
six by least squares. Least squares is the conventional choice and is the worse
one here: the Toxocara distances lie near saturation, and least squares
propagates that error inward, moving the Ascaris branch by 36%. The exact solve
insulates the well-behaved ingroup triplet, which is exactly determined by three
unsaturated distances, and sacrifices the single most saturated distance of the
six.

**The root.** Pairwise distances constrain five branches, not six, so they
cannot locate the root along the Toxocara branch. The tree is drawn unrooted
with a basal trifurcation, and the reported Toxocara value is the full unrooted
branch. Rooting with Toxocara as outgroup is a statement for the figure legend,
not something these data show.

**Support.** Each pairwise dS is a median over thousands of orthologues, so the
sampling unit is the gene. `--bootstrap` resamples every table with replacement,
recomputes all six medians, re-tests the topology and refits the branches,
reporting topology recovery and percentile intervals. A site-wise phylogenetic
bootstrap does not apply to a tree built from summary distances. A four-taxon
unrooted tree has one non-trivial split, so there is exactly one support value;
terminal branches are 100% by construction and are not annotated.

**Saturation.** Toxocara dS is 1.39-1.43, close to the ~1.5 threshold beyond
which yn00 estimates become unreliable. That branch is therefore drawn shorter
than the true divergence and the six distances are not perfectly additive; the
residual table quantifies by how much, and the discrepancy falls on the longest
path. This is a distance summary mapped onto an accepted topology, not a
phylogenetic inference. Quantitative statements about lineage-specific selection
come from the codeml four-way branch model, not from this tree.

**The diagnostic figure.** `dnds_six_way_DIAGNOSTIC` shows dN/dS across all six
comparisons and is not a manuscript figure. Absolute dN/dS is not comparable
between comparisons of different divergence depth, because dN/dS is
time-dependent and dS saturates: the Toxocara comparisons land near 0.10 and the
ascarid-ascarid ones near 0.18-0.21. Side by side that reads as Toxocara genes
being more constrained, which is the opposite of what the branch models show.
It is kept as a sanity check on the pipeline; `--no-dnds-figure` skips it.

**Figures.** Built at journal column width with absolute point sizes and audited
before saving. They must be placed at 100% and never rescaled: point sizes are
absolute, so scaling a figure scales its type below the journal's minimum.


Usage
-----
    python four_way_dnds.py \
      --as-pu   kaks_results/kaks_filtered.tsv \
      --as-bsc  kaks_As_Bsc_v2/kaks_filtered.tsv \
      --pu-bsc  kaks_Pu_Bsc_v2/kaks_filtered.tsv \
      --tox-as  kaks_Tox_As/kaks_filtered.tsv \
      --tox-pu  kaks_Tox_Pu/kaks_filtered.tsv \
      --tox-bsc kaks_Tox_Bsc_v2/kaks_filtered.tsv \
      --bootstrap 1000

Input tables need columns for the synonymous rate, the nonsynonymous rate and
their ratio, under either the dN/dS or the legacy Ka/Ks column names.

Requires numpy, pandas and matplotlib, plus paper_colors.py and rate_tables.py,
which are located automatically in lib/.
"""

import argparse
import hashlib
import json
import logging
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

# Shared modules (paper_colors, rate_tables) may sit beside this script in a flat
# layout, or in a sibling lib/ directory when scripts are filed per figure.
# Python only puts the script's own directory on sys.path, so find lib/ if it
# exists. A no-op in the flat layout, and independent of the current working
# directory, so the script runs the same under sbatch, from the shell, or from a
# clone of the repository.
_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, 'lib'),
              os.path.join(os.path.dirname(_HERE), 'lib')):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, _cand)

from paper_colors import (COL, FONT, LINEWIDTH, NEUTRAL, apply_style, audit,
                          headroom, mm, pair_ramp, panel_letter, save)
from rate_tables import read_rate_table

apply_style()


SPECIES = ['Ascaris', 'Parascaris', 'Baylisascaris', 'Toxocara']
AS, PU, BSC, TOX = 0, 1, 2, 3

# Tip labels. Genus-only for Ascaris per the manuscript convention; full
# binomials for the rest. --label-style genus overrides to genus for all four.
BINOMIAL = {
    'Ascaris': 'Ascaris',
    'Parascaris': 'P. univalens',
    'Baylisascaris': 'B. schroederi',
    'Toxocara': 'T. canis',
}

# The accepted ascaridid topology, as an unrooted split. Rooting with Toxocara
# as outgroup turns this into (((Ascaris, Parascaris), Baylisascaris), Toxocara).
REFERENCE_SPLIT = frozenset({frozenset({AS, PU}), frozenset({BSC, TOX})})

# The three possible unrooted topologies for four taxa
ALL_SPLITS = [
    frozenset({frozenset({AS, PU}), frozenset({BSC, TOX})}),
    frozenset({frozenset({AS, BSC}), frozenset({PU, TOX})}),
    frozenset({frozenset({AS, TOX}), frozenset({PU, BSC})}),
]

BRANCH_ORDER = ['Ascaris', 'Parascaris', 'Baylisascaris', 'Toxocara', 'internal']

EXTS = ('.pdf', '.svg', '.png')


def split_label(split):
    """Human-readable name for an unrooted four-taxon split."""
    a, b = sorted(split, key=lambda s: min(s))
    name = lambda s: '+'.join(SPECIES[i][:3] for i in sorted(s))
    return f'{name(a)} | {name(b)}'


def save_figure(fig, outpath):
    """Save at exact built size in every export format, after auditing.

    No bbox_inches: the figure was built at final printed size and trimming the
    bounding box would change it. Anything the audit flags is logged loudly
    rather than raised, so a run still completes and the problem is visible.
    """
    for problem in audit(fig, verbose=False):
        logging.warning(f"FIGURE AUDIT [{os.path.basename(outpath)}]: {problem}")
    stem = outpath[:-4] if outpath.endswith('.pdf') else outpath
    # Through save(), not fig.savefig: it suppresses the per-run timestamp each
    # format embeds, so a rerun on unchanged input is byte-identical and `diff`
    # can answer whether a figure actually moved.
    save(fig, os.path.dirname(stem) or '.', os.path.basename(stem),
         formats=[e.lstrip('.') for e in EXTS])
    plt.close(fig)
    logging.info(f"Saved: {stem}{'/'.join(EXTS)}")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_rate_table(path, label, max_ds=None):
    """
    Load a filtered rate TSV.

    Returns the DataFrame restricted to usable rows and the corresponding Ks
    values as an array. The array is what the bootstrap resamples, so the
    filtering applied here defines the sampling universe -- and because the same
    row mask is applied to the returned DataFrame, the plotted distributions and
    the fitted medians describe the same set of genes.
    """
    # read_rate_table maps the pipeline's legacy column names onto dN/dS/dN_dS,
    # so this script speaks the manuscript's notation regardless of file vintage.
    # speaks the manuscript's notation regardless of the file's vintage.
    try:
        df = read_rate_table(path, require=('dS', 'dN_dS'))
    except (KeyError, ValueError) as exc:
        sys.exit(f"ERROR: {exc}")

    ks = pd.to_numeric(df['dS'], errors='coerce')
    keep = ks.notna() & (ks > 0)
    if max_ds is not None:
        keep &= ks <= max_ds
    n_drop = int((~keep).sum())

    kept = df.loc[keep].copy()
    med_dnds = pd.to_numeric(kept['dN_dS'], errors='coerce').median()
    logging.info(f"Loaded {label}: {len(df)} pairs, {int(keep.sum())} usable, "
                 f"median dS={ks[keep].median():.4f}, "
                 f"median dN/dS={med_dnds:.4f}"
                 + (f"  ({n_drop} rows dropped)" if n_drop else ""))
    return kept, ks[keep].to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# Topology and branch fitting  (unchanged -- validated)
# ---------------------------------------------------------------------------

def four_point_split(m):
    """
    Pick the unrooted topology by the four-point condition.

    For four taxa the correct topology is the one whose pair of "cross"
    distances sums to the least. Returns (split, {split: sum}).
    """
    sums = {
        ALL_SPLITS[0]: m[AS][PU] + m[BSC][TOX],
        ALL_SPLITS[1]: m[AS][BSC] + m[PU][TOX],
        ALL_SPLITS[2]: m[AS][TOX] + m[PU][BSC],
    }
    return min(sums, key=sums.get), sums


# Design matrix for the reference split {As,Pu} | {Bsc,Tox}.
# Columns: bA, bP, bB, bT, bI   (bI = the internal branch)
# Rows are the six pairwise distances in the order given by _PAIRS.
_PAIRS = [(AS, PU), (AS, BSC), (AS, TOX), (PU, BSC), (PU, TOX), (BSC, TOX)]
_DESIGN = np.array([
    [1, 1, 0, 0, 0],   # d(As,Pu)  = bA + bP
    [1, 0, 1, 0, 1],   # d(As,Bsc) = bA + bI + bB
    [1, 0, 0, 1, 1],   # d(As,Tox) = bA + bI + bT
    [0, 1, 1, 0, 1],   # d(Pu,Bsc) = bP + bI + bB
    [0, 1, 0, 1, 1],   # d(Pu,Tox) = bP + bI + bT
    [0, 0, 1, 1, 0],   # d(Bsc,Tox)= bB + bT
], dtype=float)

# Rows used by the exact fit: everything except d(Pu,Tox), left as the residual.
_EXACT_ROWS = [0, 1, 2, 3, 5]


def fit_branches(m, method='ols'):
    """
    Branch lengths for the reference topology from a 4x4 distance matrix.

    method='ols'   least squares over all six distances
    method='exact' determinate solve on five distances, leaving d(Pu,Tox) as
                   the residual
    """
    d = np.array([m[i][j] for i, j in _PAIRS], dtype=float)

    if method == 'exact':
        A = _DESIGN[_EXACT_ROWS]
        y = d[_EXACT_ROWS]
    elif method == 'ols':
        A, y = _DESIGN, d
    else:
        raise ValueError(f"unknown fit method: {method}")

    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return dict(zip(BRANCH_ORDER, beta))


def fitted_distances(branches):
    """Pairwise distances implied by a set of branch lengths."""
    beta = np.array([branches[b] for b in BRANCH_ORDER], dtype=float)
    return _DESIGN @ beta


def residual_table(m, branches):
    """Observed minus fitted, per pairwise distance."""
    obs = np.array([m[i][j] for i, j in _PAIRS], dtype=float)
    fit = fitted_distances(branches)
    rows = []
    for (i, j), o, f in zip(_PAIRS, obs, fit):
        rows.append({
            'pair': f'{SPECIES[i][:3]}-{SPECIES[j][:3]}',
            'observed': o, 'fitted': f, 'residual': o - f,
            'pct': 100 * (o - f) / f if f else np.nan,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Gene-resampling bootstrap  (unchanged -- validated)
# ---------------------------------------------------------------------------

def bootstrap_tree(ks_arrays, n_rep, method, seed=0):
    """
    Resample orthologues with replacement and refit the tree each time.

    The unit of resampling is the gene, because each pairwise dS in the matrix
    is a median over genes. Each replicate draws a new sample from every
    pairwise table independently (the tables do not share a gene set), rebuilds
    all six medians, re-tests the topology by the four-point condition, and
    refits the branches.

    Branch confidence intervals are computed only over replicates that recover
    the reference topology, since branch lengths are not comparable across
    topologies.
    """
    rng = np.random.default_rng(seed)
    split_counts = {s: 0 for s in ALL_SPLITS}
    kept = []

    for _ in range(n_rep):
        m = np.zeros((4, 4))
        for (i, j), v in ks_arrays.items():
            idx = rng.integers(0, len(v), len(v))
            med = float(np.median(v[idx]))
            m[i][j] = m[j][i] = med

        split, _ = four_point_split(m)
        split_counts[split] += 1
        if split == REFERENCE_SPLIT:
            kept.append(fit_branches(m, method=method))

    reps = pd.DataFrame(kept) if kept else pd.DataFrame(columns=BRANCH_ORDER)
    ci = {}
    for b in BRANCH_ORDER:
        if len(reps):
            ci[b] = (float(np.percentile(reps[b], 2.5)),
                     float(np.percentile(reps[b], 97.5)))
        else:
            ci[b] = (np.nan, np.nan)

    return {
        'n_rep': n_rep,
        'split_counts': split_counts,
        'support': 100.0 * split_counts[REFERENCE_SPLIT] / n_rep,
        'branches': reps,
        'ci': ci,
    }


# ---------------------------------------------------------------------------
# Figure 1D -- dS distance tree
# ---------------------------------------------------------------------------

def nice_scale(span):
    """A round scale-bar length roughly a fifth of the axis span."""
    target = span / 5.0
    for step in (0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0, 2.0):
        if step >= target:
            return step
    return 2.0


def plot_four_way_tree(branches, outpath, support=None, fit_method='exact',
                       break_long=True, break_ratio=2.5, width='single',
                       height_mm=72.0, label_style='binomial', letter=None):
    """
    Unrooted dS distance tree drawn with a basal trifurcation.

    Toxocara, Baylisascaris and the (Ascaris, Parascaris) clade radiate from one
    node. That is an honest depiction of what pairwise distances support: the
    topology is resolved, the position of the root along the Toxocara branch is
    not. Rooting with Toxocara as outgroup belongs in the legend.

    Branch lengths are labelled along branches, never at nodes -- a number at a
    node reads as a support value to any phylogeneticist. There is exactly ONE
    support value because a four-taxon unrooted tree has n-3 = 1 internal
    branch, hence one non-trivial split. Terminal branches define trivial splits
    (one taxon versus the other three) present in every possible tree, so they
    are 100% by construction and are never annotated.

    The Toxocara branch is ~8x the whole ingroup, which compresses Ascaris,
    Parascaris, Baylisascaris and the internal branch into about a tenth of the
    figure width and makes their labels collide. break_long draws that branch
    truncated, with a break symbol, so the ingroup gets legible space. The
    printed value is always the true length; the summary records the break.

    Explanatory notes for this panel are written to the summary file rather
    than into the artwork: they are caption material, and at final printed size
    they would fall below the journal's 6 pt floor.
    """
    bA = branches['Ascaris']
    bP = branches['Parascaris']
    bB = branches['Baylisascaris']
    bT = branches['Toxocara']
    bI = branches['internal']

    fig, ax = plt.subplots(figsize=(COL[width], mm(height_mm)))

    # y positions, top to bottom
    y_tox, y_bsc, y_as, y_pu = 3.0, 2.0, 1.0, 0.0
    y_n1 = (y_as + y_pu) / 2.0

    x_n2 = 0.0
    x_n1 = x_n2 + bI

    # Ingroup extent sets the drawing scale; Toxocara is truncated to it when it
    # would otherwise dominate.
    ingroup_max = max(bB, bI + bA, bI + bP)
    broken = bool(break_long and bT > break_ratio * ingroup_max)
    x_tox_draw = (ingroup_max * 1.75) if broken else bT

    tips = {
        'Toxocara':      (x_tox_draw, y_tox),
        'Baylisascaris': (x_n2 + bB, y_bsc),
        'Ascaris':       (x_n1 + bA, y_as),
        'Parascaris':    (x_n1 + bP, y_pu),
    }
    max_x = max(x for x, _ in tips.values())

    lw = LINEWIDTH['data'] * 1.4
    ink = NEUTRAL['ink']
    # basal trifurcation: vertical spine at n2 spanning its three children
    ax.plot([x_n2, x_n2], [y_n1, y_tox], '-', color=ink, lw=lw,
            solid_capstyle='round')

    if broken:
        # two solid segments with a gap and a pair of slashes marking the break
        xb = x_tox_draw * 0.60
        gap = x_tox_draw * 0.05
        ax.plot([x_n2, xb], [y_tox, y_tox], '-', color=ink, lw=lw)
        ax.plot([xb + gap, x_tox_draw], [y_tox, y_tox], '-', color=ink, lw=lw)
        for dx in (0.0, gap):
            ax.plot([xb + dx - 0.018 * x_tox_draw, xb + dx + 0.018 * x_tox_draw],
                    [y_tox - 0.16, y_tox + 0.16], '-', color=ink,
                    lw=lw * 0.7, solid_capstyle='round')
    else:
        ax.plot([x_n2, x_tox_draw], [y_tox, y_tox], '-', color=ink, lw=lw)

    ax.plot([x_n2, tips['Baylisascaris'][0]], [y_bsc, y_bsc], '-', color=ink, lw=lw)
    ax.plot([x_n2, x_n1], [y_n1, y_n1], '-', color=ink, lw=lw)
    # (Ascaris, Parascaris) clade
    ax.plot([x_n1, x_n1], [y_pu, y_as], '-', color=ink, lw=lw,
            solid_capstyle='round')
    ax.plot([x_n1, tips['Ascaris'][0]], [y_as, y_as], '-', color=ink, lw=lw)
    ax.plot([x_n1, tips['Parascaris'][0]], [y_pu, y_pu], '-', color=ink, lw=lw)

    label_x = max_x + 0.04 * max_x

    # dotted leaders so the taxon labels align despite unequal branch lengths
    for name, (x, y) in tips.items():
        text = name if label_style == 'genus' else BINOMIAL[name]
        ax.plot([x, label_x], [y, y], ls=':', color=NEUTRAL['ink_light'],
                lw=LINEWIDTH['hairline'], zorder=1)
        ax.text(label_x + 0.01 * max_x, y, text, fontsize=FONT['axis_label'],
                fontstyle='italic', va='center', color=ink)

    def label_branch(x1, x2, y, value, dy=0.15, ha='center'):
        """Label a branch, positioned over the branch it belongs to.

        ha='right' anchors the text at the branch's far end so it extends back
        along the branch -- needed for the internal branch, which starts at x=0
        and would otherwise overhang the basal spine.
        """
        pad = 0.008 * max_x   # keep right-anchored text clear of the clade spine
        x = {'center': (x1 + x2) / 2.0, 'right': x2 - pad, 'left': x1 + pad}[ha]
        ax.text(x, y + dy, f'{value:.4f}', ha=ha, va='bottom',
                fontsize=FONT['annotation'], color=ink)

    label_branch(x_n2, tips['Toxocara'][0], y_tox, bT)
    label_branch(x_n2, tips['Baylisascaris'][0], y_bsc, bB)
    label_branch(x_n1, tips['Ascaris'][0], y_as, bA)
    label_branch(x_n1, tips['Parascaris'][0], y_pu, bP, dy=-0.34)
    # internal: value above, support below, both right-anchored at the branch
    # end so the pair reads unambiguously as belonging to this branch
    # The internal branch is short -- its label is more than twice its drawn
    # length -- and it sits in the crook of the trifurcation, between the basal
    # spine and the (Ascaris, Parascaris) clade spine. On the branch it collides
    # with one or the other, and with bootstrap support drawn it would collide
    # with itself as well. Both values go in the empty margin left of the basal
    # spine, vertically centred on the branch they describe.
    x_margin = x_n2 - 0.035 * max_x
    ax.text(x_margin, y_n1 + 0.07, f'{bI:.4f}', ha='right', va='bottom',
            fontsize=FONT['annotation'], color=ink)
    if support is not None:
        # Support in ink rather than colour. Blue and red are reserved for
        # region identity throughout the paper and are not reused here.
        ax.text(x_margin, y_n1 - 0.07, f'{support:.0f}', ha='right', va='top',
                fontsize=FONT['annotation'], fontweight='bold', color=ink)

    # scale bar
    bar = nice_scale(ingroup_max)
    y_bar = -0.58
    ax.plot([0, bar], [y_bar, y_bar], '-', color=ink, lw=lw,
            solid_capstyle='butt')
    ax.text(bar / 2.0, y_bar - 0.16, f'{bar:g} substitutions per synonymous site',
            ha='center', va='top', fontsize=FONT['annotation'], color=ink)

    # Right margin sized for the longest tip label ('B. schroederi').
    # There is no bbox_inches='tight' to rescue an overhang -- the figure is
    # saved at exactly the size it was built at, so text that runs past the
    # canvas is simply lost.
    ax.set_xlim(-0.40 * max_x, label_x + 0.62 * max_x)
    ax.set_ylim(-1.05, 3.5)
    ax.axis('off')
    if letter:
        panel_letter(ax, letter, dx=-0.02, dy=0.97)

    save_figure(fig, outpath)
    return broken


# ---------------------------------------------------------------------------
# Figure 1C -- dS distribution overlay
# ---------------------------------------------------------------------------

def plot_ks_overlay(dfs, labels, outpath, width='single', height_mm=62.0,
                    xmax=1.8, binwidth=0.02, letter=None,
                    legend_loc='auto'):
    """
    dS distributions for all six comparisons, ordered by median.

    Six comparisons ordered by divergence are a magnitude rather than six
    identities, so they take a single-hue sequential ramp assigned in ascending
    median dS. The grouping that matters -- three ascarid-ascarid pairs at the
    light end, three Toxocara pairs at the dark end -- then emerges from the
    figure rather than being asserted by it. Blue and red are not used: they are
    reserved for retained and eliminated regions elsewhere in the paper.

    Distributions are drawn as step outlines. Six translucent filled histograms
    overlapping are unreadable.

    Counts rather than densities are appropriate because the six sample sizes
    lie within about 20% of one another; densities would be needed if they did
    not.
    """
    order = np.argsort([pd.to_numeric(d['dS'], errors='coerce').median()
                        for d in dfs])
    ramp = pair_ramp(len(dfs))

    fig, ax = plt.subplots(figsize=(COL[width], mm(height_mm)))
    bins = np.arange(0, xmax + binwidth, binwidth)

    for colour, k in zip(ramp, order):
        ks = pd.to_numeric(dfs[k]['dS'], errors='coerce').dropna()
        ax.hist(ks, bins=bins, histtype='step', color=colour,
                linewidth=LINEWIDTH['data'],
                label=f'{labels[k]}  ({len(ks):,}; {ks.median():.3f})')

    ax.set_xlabel('dS (synonymous substitutions per site)')
    ax.set_ylabel('Ortholog pairs')
    ax.set_xlim(0, xmax)
    # A legend above the axes costs three rows of panel height. That is
    # affordable when the panel is narrow and tall, and wasteful when it is wide
    # and short. At single-column width the space above the axes is the only
    # place the legend fits; at wider settings the upper right of the axes is
    # empty -- the distributions sit at the far left and along the lower right --
    # so the legend goes inside and the panel can be a short strip. headroom()
    # lifts the y ceiling first so the legend cannot land on a bar.
    #
    # Anchored at the axes left edge when above: the panel letter sits in the
    # y-label margin at negative x, so the two do not collide, whereas pushing
    # the legend right to make room runs its second column off the canvas.
    if legend_loc == 'auto':
        legend_loc = 'above' if width == 'single' else 'inside'
    common = dict(ncol=2, frameon=False, fontsize=FONT['legend'],
                  handlelength=1.2, columnspacing=1.0, borderaxespad=0.0,
                  title='comparison  (n; median dS)',
                  title_fontsize=FONT['legend'])
    if legend_loc == 'inside':
        headroom(ax, 0.30)
        ax.legend(loc='upper right', **common)
    else:
        ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), **common)
    if letter:
        panel_letter(ax, letter, dx=-0.13, dy=1.02)

    save_figure(fig, outpath)


# ---------------------------------------------------------------------------
# Diagnostic only -- NOT a manuscript figure
# ---------------------------------------------------------------------------

def plot_dnds_diagnostic(dfs, labels, outpath, rng, box_alpha=0.55,
                         width='double',
                         height_mm=70.0):
    """
    dN/dS across all six comparisons. INTERNAL DIAGNOSTIC, NOT A PAPER FIGURE.

    Absolute dN/dS is not comparable between comparisons of different divergence
    depth. dN/dS is time-dependent and dS saturates, so the Toxocara comparisons
    land near 0.10 while the ascarid-ascarid ones land near 0.18-0.21. Placed
    side by side this reads as Toxocara genes being under stronger purifying
    selection, which is the opposite of what the branch models show. Only the
    within-comparison retained-versus-eliminated contrast is interpretable, and
    that is what Figure 5 reports.

    It is kept as a sanity check on the pipeline: a comparison whose whole
    distribution shifts unexpectedly means something upstream has changed.
    """
    ramp = pair_ramp(len(dfs))
    order = np.argsort([pd.to_numeric(d['dS'], errors='coerce').median()
                        for d in dfs])

    fig, ax = plt.subplots(figsize=(COL[width], mm(height_mm)))
    data = [pd.to_numeric(dfs[k]['dN_dS'], errors='coerce').dropna().values
            for k in order]
    positions = range(len(order))

    bp = ax.boxplot(data, positions=positions, widths=0.5,
                    patch_artist=True, showfliers=False, zorder=2)
    for patch, colour in zip(bp['boxes'], ramp):
        # Alpha in the face colour, never set_alpha() -- that would fade the black
        # edge along with the fill. 0.55 rather than the 0.40 used on the
        # two-colour region panels: adjacent steps of the six-step ramp are only
        # dE 4.4 apart at 0.40 composited on white, which is below the ~5 needed
        # to tell them apart side by side in print. At 0.55 the minimum is 6.2.
        patch.set_facecolor(to_rgba(colour, box_alpha))
        patch.set_alpha(None)
        patch.set_edgecolor('black')
        patch.set_linewidth(LINEWIDTH['thin'])
    for part, z in (('whiskers', 4), ('caps', 4), ('medians', 5)):
        for artist in bp[part]:
            artist.set_color('black')
            artist.set_zorder(z)
            artist.set_linewidth(LINEWIDTH['data'] if part == 'medians'
                                 else LINEWIDTH['thin'])

    for i, (d, colour) in enumerate(zip(data, ramp)):
        d_plot = d[rng.choice(len(d), 1500, replace=False)] if len(d) > 1500 else d
        jitter = rng.normal(0, 0.06, size=len(d_plot))
        ax.scatter(np.full(len(d_plot), i) + jitter, d_plot,
                   s=1.5, alpha=0.15, color=colour,
                   edgecolors='none', zorder=3, rasterized=True)

    ax.set_xticks(list(positions))
    ax.set_xticklabels([f'{labels[k]}\n(n = {len(d):,})'
                        for k, d in zip(order, data)])
    ax.set_ylabel('dN/dS')
    # Neutrality line in ink, not red -- red is the eliminated-region colour.
    ax.axhline(1.0, color=NEUTRAL['rule'], linestyle='--',
               linewidth=LINEWIDTH['rule'])
    ax.set_ylim(0, 1.2)

    save_figure(fig, outpath)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def file_fingerprint(path):
    """Size, mtime and a SHA-256 prefix for one input file."""
    st = os.stat(path)
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return {
        'path': os.path.abspath(path),
        'bytes': st.st_size,
        'mtime_utc': datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
        'sha256_16': h.hexdigest()[:16],
    }


def git_commit():
    """Current commit of the repo containing this script, if there is one."""
    try:
        here = os.path.dirname(os.path.abspath(__file__)) or '.'
        out = subprocess.run(['git', '-C', here, 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            dirty = subprocess.run(['git', '-C', here, 'status', '--porcelain'],
                                   capture_output=True, text=True, timeout=5)
            return out.stdout.strip() + ('+dirty' if dirty.stdout.strip() else '')
    except Exception:
        pass
    return None


def write_provenance(outdir, args, inputs, row_counts):
    """Record exactly what produced the contents of this directory.

    Row counts and hashes let any figure be traced back to the exact tables it
    was built from, which is otherwise impossible to establish once several
    versions of an input exist.
    """
    rec = {
        'script': os.path.basename(__file__),
        'git_commit': git_commit(),
        'run_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'host': platform.node(),
        'python': sys.version.split()[0],
        'packages': {
            'numpy': np.__version__,
            'pandas': pd.__version__,
            'matplotlib': matplotlib.__version__,
        },
        'argv': sys.argv,
        'args': vars(args),
        'inputs': inputs,
        'usable_rows': row_counts,
    }
    path = os.path.join(outdir, 'provenance.json')
    with open(path, 'w') as fh:
        json.dump(rec, fh, indent=2)
    logging.info(f"Provenance: {path}")
    return rec


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Four-way dN/dS comparison and dS distance tree.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--as-pu', required=True,
                   help='Ascaris vs Parascaris kaks_filtered.tsv')
    p.add_argument('--as-bsc', required=True,
                   help='Ascaris vs Baylisascaris kaks_filtered.tsv (use _v2)')
    p.add_argument('--pu-bsc', required=True,
                   help='Parascaris vs Baylisascaris kaks_filtered.tsv (use _v2)')
    p.add_argument('--tox-as', required=True,
                   help='Toxocara vs Ascaris kaks_filtered.tsv')
    p.add_argument('--tox-pu', required=True,
                   help='Toxocara vs Parascaris kaks_filtered.tsv')
    p.add_argument('--tox-bsc', required=True,
                   help='Toxocara vs Baylisascaris kaks_filtered.tsv (use _v2)')
    p.add_argument('--fit', choices=['exact', 'ols'], default='exact',
                   help="Branch-length fitting. 'exact' (default) solves five "
                        "distances exactly and leaves d(Pu,Tox) -- the most "
                        "saturated -- as the residual, which keeps the saturated "
                        "Toxocara distances from perturbing the well-estimated "
                        "ingroup branches. 'ols' fits all six by least squares, "
                        "which is standard practice but here lets saturation "
                        "propagate inward. Both are always reported; this "
                        "chooses which the figure and CIs use.")
    p.add_argument('--no-break', dest='break_long', action='store_false',
                   help='Draw the Toxocara branch fully to scale instead of '
                        'truncating it. It is ~8x the whole ingroup, so this '
                        'compresses the other four branches into about a tenth '
                        'of the figure width and their labels crowd.')
    p.add_argument('--bootstrap', type=int, default=1000,
                   help='Gene-resampling replicates for topology support and '
                        'branch CIs (default: 1000; 0 disables)')
    p.add_argument('--max-ds', type=float, default=None,
                   help='Drop dS above this before taking medians and before '
                        'plotting, so the fitted distances and the plotted '
                        'distributions always describe the same genes. Default: '
                        'no upper bound. Try 3 as a saturation sensitivity check.')
    p.add_argument('--seed', type=int, default=0,
                   help='RNG seed for the bootstrap and scatter subsampling')
    p.add_argument('--tree-width', choices=['single', 'onehalf', 'double'],
                   default='single',
                   help='Journal column width for the tree (default: single, '
                        'i.e. 89 mm, for a Figure 1 with C and D side by side)')
    p.add_argument('--overlay-width', choices=['single', 'onehalf', 'double'],
                   default='single',
                   help='Journal column width for the Ks overlay (default: single)')
    p.add_argument('--tree-height-mm', type=float, default=72.0)
    p.add_argument('--overlay-height-mm', type=float, default=62.0)
    p.add_argument('--overlay-legend', choices=['auto', 'above', 'inside'],
                   default='auto',
                   help="Placement of the overlay legend. 'auto' (default) puts "
                        "it above the axes at single-column width, where there is "
                        "no room inside, and in the empty upper right at wider "
                        "settings, where a legend above would waste panel height.")
    p.add_argument('--overlay-xmax', type=float, default=1.8,
                   help='Upper x limit of the Ks overlay (default: 1.8, which '
                        'covers the Toxocara comparisons at 1.39-1.43)')
    p.add_argument('--label-style', choices=['binomial', 'genus'],
                   default='binomial',
                   help="Tree tip labels. 'binomial' gives Ascaris, "
                        "P. univalens, B. schroederi, T. canis, matching the "
                        "manuscript naming convention.")
    p.add_argument('--panel-letters', action='store_true',
                   help='Draw panel letters D and C on the tree and overlay. '
                        'Leave off if the letters are added during assembly.')
    p.add_argument('--box-alpha', type=float, default=0.55,
                   help='Opacity of the box fills in the six-way diagnostic. '
                        'The jitter is drawn over the boxes, so a high value '
                        'hides it; 0.55 rather than the 0.40 used on the '
                        'two-colour region panels because adjacent steps of the '
                        'six-step ramp fall to dE 4.4 at 0.40, below the ~5 '
                        'needed to separate them in print (default: 0.55)')
    p.add_argument('--no-dnds-figure', dest='dnds_figure', action='store_false',
                   help='Skip the six-way dN/dS boxplot. It is a pipeline '
                        'diagnostic, not a manuscript figure: absolute dN/dS is '
                        'not comparable across comparisons of different '
                        'divergence depth.')
    p.add_argument('--outdir', default=None,
                   help='Output directory. Created if absent. Defaults to '
                        'results/four_way_dnds_<YYYYMMDD> under the current '
                        'working directory, so outputs land next to the data '
                        'rather than in a distant home tree, and so '
                        'a rerun on a later day cannot overwrite the run a '
                        'figure was built from. Pass an explicit path to '
                        'control it; the SLURM wrapper does exactly that.')
    args = p.parse_args()
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'four_way_dnds_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    # mode='w': FileHandler appends by default, so rerunning into an
    # existing directory concatenates logs while the figures and
    # provenance.json are overwritten -- the log then describes more
    # runs than the directory contains.
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(os.path.join(args.outdir, 'four_way.log'),
                                     mode='w'),
                  logging.StreamHandler(sys.stdout)])
    for noisy in ('matplotlib.font_manager', 'fontTools', 'fontTools.subset',
                  'PIL'):
        logging.getLogger(noisy).setLevel(logging.ERROR)

    logging.info("=" * 60)
    logging.info("Four-Way dN/dS Analysis")
    logging.info("=" * 60)

    spec = [('As x Pu', args.as_pu, (AS, PU)),
            ('As x Bsc', args.as_bsc, (AS, BSC)),
            ('Pu x Bsc', args.pu_bsc, (PU, BSC)),
            ('Tox x As', args.tox_as, (AS, TOX)),
            ('Tox x Pu', args.tox_pu, (PU, TOX)),
            ('Tox x Bsc', args.tox_bsc, (BSC, TOX))]

    dfs, labels, ks_arrays, inputs, row_counts = [], [], {}, {}, {}
    for label, path, key in spec:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: input not found: {path}")
        df, arr = load_rate_table(path, label, max_ds=args.max_ds)
        dfs.append(df)
        labels.append(label)
        ks_arrays[key] = arr
        inputs[label] = file_fingerprint(path)
        row_counts[label] = int(len(arr))

    # --- distance matrix -----------------------------------------------------
    m = np.zeros((4, 4))
    for (i, j), arr in ks_arrays.items():
        m[i][j] = m[j][i] = float(np.median(arr))

    logging.info("\ndS distance matrix:")
    logging.info(f"{'':>15}" + ''.join(f'{s:>14}' for s in SPECIES))
    for i, sp in enumerate(SPECIES):
        vals = ''.join(f'{m[i][j]:>14.4f}' if j != i else f'{"—":>14}'
                       for j in range(4))
        logging.info(f"{sp:>15}{vals}")

    # --- topology ------------------------------------------------------------
    best, sums = four_point_split(m)
    logging.info("\nFour-point condition (smallest sum wins):")
    for s in ALL_SPLITS:
        mark = '  <- selected' if s == best else ''
        logging.info(f"  {split_label(s):<34} {sums[s]:.4f}{mark}")
    if best != REFERENCE_SPLIT:
        logging.warning(f"Selected topology {split_label(best)} differs from the "
                        f"accepted ascaridid topology {split_label(REFERENCE_SPLIT)}. "
                        f"Branch lengths below are still fitted to the accepted "
                        f"topology; treat with care.")

    # --- branch lengths, both fits ------------------------------------------
    fits = {meth: fit_branches(m, method=meth) for meth in ('ols', 'exact')}
    branches = fits[args.fit]

    logging.info("\nBranch-specific dS (topology: (((As,Pu),Bsc),Tox)):")
    logging.info(f"  {'branch':<16}{'OLS':>12}{'exact':>12}")
    for b in BRANCH_ORDER:
        logging.info(f"  {b:<16}{fits['ols'][b]:>12.4f}{fits['exact'][b]:>12.4f}")
    logging.info(f"  (figure uses --fit {args.fit})")
    logging.info("  The Toxocara branch is the FULL unrooted branch. Pairwise "
                 "distances cannot place the root along it; root with Toxocara "
                 "as outgroup in the legend, not in the data.")

    neg = {b: round(v, 4) for b, v in branches.items() if v < 0}
    if neg:
        logging.warning(f"NEGATIVE branch length(s): {neg}")
        logging.warning("  Indicates non-additive distances, expected as Ks "
                        "approaches saturation. Use the codeml four-way branch "
                        "model for quantitative branch dN/dS.")

    resid = residual_table(m, branches)
    logging.info("\nAdditivity residuals (observed - fitted):")
    for _, r in resid.iterrows():
        logging.info(f"  {r['pair']:<10} obs={r['observed']:.4f} "
                     f"fit={r['fitted']:.4f} resid={r['residual']:+.4f} "
                     f"({r['pct']:+.1f}%)")
    # A four-taxon tree has 6 distances and 5 branches, so exactly one degree of
    # freedom. Where that single discrepancy shows up is a property of the fit,
    # not of the data: OLS spreads it evenly over the four cross-distances and
    # fits the two sister-pair distances exactly, while the exact solve puts all
    # of it on the one distance it drops.
    worst = resid.loc[resid['residual'].abs().idxmax()]
    if args.fit == 'ols':
        logging.info(f"  One degree of freedom, so the discrepancy is fixed in "
                     f"total and least squares spreads it evenly: "
                     f"±{abs(worst['residual']):.4f} on each cross-distance, "
                     f"with Asc-Par and Bay-Tox fitting exactly.")
    else:
        logging.info(f"  The exact solve concentrates the entire discrepancy on "
                     f"{worst['pair']} ({worst['pct']:+.1f}%), the distance it "
                     f"drops. That is by design: it is the most saturated of the "
                     f"six, so it is the one to sacrifice.")
    logging.info("  Either way the discrepancy is non-additivity from Ks "
                 "approaching saturation on the Toxocara branch, not a fitting "
                 "error. Use the codeml four-way branch model for quantitative "
                 "branch dN/dS.")

    # --- bootstrap -----------------------------------------------------------
    boot = None
    if args.bootstrap > 0:
        logging.info(f"\nGene-resampling bootstrap: {args.bootstrap} replicates "
                     f"(resampling orthologues with replacement)")
        boot = bootstrap_tree(ks_arrays, args.bootstrap, args.fit, seed=args.seed)
        for s in ALL_SPLITS:
            pct = 100.0 * boot['split_counts'][s] / boot['n_rep']
            logging.info(f"  {split_label(s):<34} {pct:6.1f}%")
        logging.info(f"  Support for the accepted topology: {boot['support']:.1f}%")
        logging.info("  Branch 95% percentile intervals:")
        for b in BRANCH_ORDER:
            lo, hi = boot['ci'][b]
            logging.info(f"    {b:<16}{branches[b]:>9.4f}  [{lo:.4f}, {hi:.4f}]")

    # --- figures -------------------------------------------------------------
    rng = np.random.default_rng(args.seed)
    broken = plot_four_way_tree(
        branches, os.path.join(args.outdir, 'ds_tree_four_way.pdf'),
        support=boot['support'] if boot else None,
        fit_method=args.fit, break_long=args.break_long,
        width=args.tree_width, height_mm=args.tree_height_mm,
        label_style=args.label_style,
        letter='D' if args.panel_letters else None)

    plot_ks_overlay(dfs, labels, os.path.join(args.outdir, 'ds_overlay.pdf'),
                    width=args.overlay_width, height_mm=args.overlay_height_mm,
                    xmax=args.overlay_xmax, legend_loc=args.overlay_legend,
                    letter='C' if args.panel_letters else None)

    if args.dnds_figure:
        plot_dnds_diagnostic(
            dfs, labels,
            os.path.join(args.outdir, 'dnds_six_way_DIAGNOSTIC.pdf'), rng,
            box_alpha=args.box_alpha)

    write_provenance(args.outdir, args, inputs, row_counts)

    # --- summary -------------------------------------------------------------
    summary_path = os.path.join(args.outdir, 'four_way_summary.txt')
    with open(summary_path, 'w') as fh:
        fh.write("Four-Way dN/dS Summary\n")
        fh.write("=" * 60 + "\n\n")
        fh.write("Topology: (((Ascaris, Parascaris), Baylisascaris), Toxocara)\n")
        fh.write("This is a Ks distance summary mapped onto the accepted\n")
        fh.write("ascaridid topology, not a phylogenetic inference.\n\n")

        for label, df in zip(labels, dfs):
            ks = pd.to_numeric(df['dS'], errors='coerce')
            kk = pd.to_numeric(df['dN_dS'], errors='coerce')
            fh.write(f"{label}:\n")
            fh.write(f"  Ortholog pairs: {len(df)}\n")
            fh.write(f"  Median dS: {ks.median():.4f}\n")
            fh.write(f"  Median dN/dS: {kk.median():.4f}\n\n")

        fh.write("dS Distance Matrix:\n")
        fh.write(f"{'':>15}" + ''.join(f'{s:>14}' for s in SPECIES) + '\n')
        for i, sp in enumerate(SPECIES):
            vals = ''.join(f'{m[i][j]:>14.4f}' if j != i else f'{"—":>14}'
                           for j in range(4))
            fh.write(f"{sp:>15}{vals}\n")

        fh.write("\nFour-point condition (smallest sum selects the topology):\n")
        for s in ALL_SPLITS:
            mark = '  <- selected' if s == best else ''
            fh.write(f"  {split_label(s):<34} {sums[s]:.4f}{mark}\n")

        fh.write("\nBranch-specific dS:\n")
        fh.write(f"  {'branch':<16}{'OLS':>12}{'exact':>12}\n")
        for b in BRANCH_ORDER:
            fh.write(f"  {b:<16}{fits['ols'][b]:>12.4f}{fits['exact'][b]:>12.4f}\n")
        fh.write(f"  Figure and CIs use --fit {args.fit}.\n\n")
        fh.write("  'exact' solves five distances exactly and leaves d(Pu,Tox)\n")
        fh.write("  as the residual. 'ols' fits all six by least squares. LS is\n")
        fh.write("  the usual choice, but here it is the worse one: the Toxocara\n")
        fh.write("  distances are near saturation, and LS lets that error\n")
        fh.write("  propagate into the ingroup, moving the Ascaris branch by 36%.\n")
        fh.write("  The exact solve insulates the ingroup triplet (dS 0.20-0.37,\n")
        fh.write("  no saturation, exactly determined by three distances) and\n")
        fh.write("  sacrifices d(Pu,Tox), which is the most saturated of the six.\n")
        fh.write("  That is why 'exact' is the default.\n\n")
        fh.write("  The Toxocara value is the full unrooted branch. Pairwise\n")
        fh.write("  distances constrain five branches, not six, so they cannot\n")
        fh.write("  locate the root along it.\n\n")

        fh.write("Additivity residuals (observed - fitted):\n")
        fh.write(f"  {'pair':<10}{'observed':>11}{'fitted':>11}"
                 f"{'residual':>11}{'pct':>9}\n")
        for _, r in resid.iterrows():
            fh.write(f"  {r['pair']:<10}{r['observed']:>11.4f}{r['fitted']:>11.4f}"
                     f"{r['residual']:>+11.4f}{r['pct']:>+9.1f}\n")
        fh.write("  Toxocara dS is 1.39-1.43, near the ~1.5 saturation threshold,\n")
        fh.write("  so the longest path is compressed most and the distances are\n")
        fh.write("  not perfectly additive.\n\n")

        if boot:
            fh.write(f"Gene-resampling bootstrap ({boot['n_rep']} replicates,\n")
            fh.write("orthologues resampled with replacement within each pair):\n")
            for s in ALL_SPLITS:
                pct = 100.0 * boot['split_counts'][s] / boot['n_rep']
                fh.write(f"  {split_label(s):<34} {pct:6.1f}%\n")
            fh.write(f"\n  Support for the accepted topology: "
                     f"{boot['support']:.1f}%\n\n")
            fh.write("  Branch 95% percentile intervals:\n")
            for b in BRANCH_ORDER:
                lo, hi = boot['ci'][b]
                fh.write(f"    {b:<16}{branches[b]:>9.4f}  [{lo:.4f}, {hi:.4f}]\n")
            fh.write("\n  A site-wise phylogenetic bootstrap does not apply to a\n")
            fh.write("  tree built from summary distances; the gene is the\n")
            fh.write("  sampling unit here.\n")

        # Legend text for the tree panel, kept out of the artwork itself.
        fh.write("\n\nFOR THE FIGURE 1D LEGEND\n")
        fh.write("-" * 40 + "\n")
        fh.write(f"  Branch lengths fitted by "
                 f"{'least squares over all six distances' if args.fit == 'ols' else 'exact solve on five distances'}.\n")
        if broken:
            fh.write("  The Toxocara branch is drawn truncated (break marks); the\n")
            fh.write("  printed value is its true length and all other branches\n")
            fh.write("  are to scale.\n")
        if boot:
            fh.write("  The bold value on the internal branch is the percentage of\n")
            fh.write(f"  {boot['n_rep']} gene-resampling replicates recovering that\n")
            fh.write("  split -- the only non-trivial split in a four-taxon\n")
            fh.write("  unrooted tree. Terminal branches are 100% by construction\n")
            fh.write("  and are not annotated.\n")

        fh.write("\n\nNOTE ON dnds_six_way_DIAGNOSTIC\n")
        fh.write("-" * 40 + "\n")
        fh.write("  Internal diagnostic, not a manuscript figure. Absolute dN/dS\n")
        fh.write("  is not comparable between comparisons of different divergence\n")
        fh.write("  depth: dN/dS is time-dependent and dS saturates, so the\n")
        fh.write("  Toxocara comparisons sit near 0.10 and the ascarid-ascarid\n")
        fh.write("  ones near 0.18-0.21. That difference is not a difference in\n")
        fh.write("  constraint. Only the within-comparison retained-versus-\n")
        fh.write("  eliminated contrast is interpretable.\n")

    logging.info(f"\nSummary saved: {summary_path}")
    logging.info("Done.")


if __name__ == '__main__':
    main()
