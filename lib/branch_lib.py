#!/usr/bin/env python3
"""
branch_lib.py -- shared logic for the four-taxon branch-model analyses.

`branch_dnds_four.py` runs codeml and writes the annotated table;
`branch_dnds_decomposition.py` re-analyses that table. Before this module they
carried two copies of the tree constants, two copies of the significance-star
function, two different BH implementations and two different ideas of which
genes are usable on a branch -- so the two figures of one manuscript figure
could disagree about n, and did.

Everything here is the single definition of:

  * the tree, the tip-to-node mapping and how a branch string is interpreted
  * which genes are usable for one metric on one branch (`usable_mask`)
  * the retained-vs-eliminated comparison (`compare`, `split_and_compare`)
  * loading and validating the annotated table

THE TREE MUST BE WRITTEN UNROOTED -- read this before touching TREE_TOPOLOGY.

codeml with `clock = 0` does NOT unroot a rooted tree. Given a basal
bifurcation it prints "This is a rooted tree, without clock. Check." and goes
on to fit both branches at the root, although only their sum is estimable
(codeml.c; the clock-test Methods describe the same behaviour in baseml). Until
2026-10-10 TREE_TOPOLOGY was written rooted,
`(((Ascaris,Parascaris),Baylisascaris),Toxocara);`, and this docstring wrongly
said codeml would unroot it. It did not: the Toxocara branch and the
root-to-ascarid branch split one length arbitrarily. In the 2026-10-09 run,
1,268 genes had Toxocara dN = 0 with the whole path's dN on the other root
branch, and the "internal branch" recorded was that other root branch, not the
internal branch; the two summed gave a dS median of 1.277 against 1.267 on the
Ks tree. Run on a simulated quartet with the codeml build used for testing, the
rooted tree gives six branches and Toxocara dS 0.295 + 0.128, the unrooted tree
five branches and Toxocara dS 0.423, with the other three branches identical.
TREE_TOPOLOGY is now written unrooted, `check_unrooted()` lets a caller refuse
a rooted string before running anything, and `label_branches()` refuses a
table with more than 2n - 3 branches.

THE INTERNAL BRANCH IS NOT AN ANCESTOR BRANCH -- read this before quoting it.

An unrooted four-taxon tree has exactly ONE internal branch, and it separates

    {Ascaris, Parascaris}   from   {Baylisascaris, Toxocara}

It is not "the ascarid ancestor" and it is not "the branch leading to the
ascarid clade": those are not identifiable without a root, and a four-taxon
unrooted tree has none. The previous code stored `internal_omega_min` -- a
minimum over a list that in this topology always has exactly one element --
which invited exactly that misreading. `label_branches()` now returns the tip
set the internal branch actually subtends, computed from codeml's own branch
list rather than assumed, so a change of topology cannot silently relabel it.

ZERO ON A LOG AXIS. codeml bounds the free-ratio dN/dS to [0.0001, 999], so a
branch with no nonsynonymous change is reported as dN/dS = 0.0001. Its dN is
printed as (dN/dS) x dS to four decimals: 0.0000 on a short branch, but 0.0001
once dS reaches 0.5. On the Toxocara branch (median dS 1.27) most such genes
therefore show dN = 0.0001, not 0, so the dN/dS zeros (Figure 5A) and the dN
zeros (Figure 5C) are counted separately and need not agree (2026-10-10 run:
21 against 1 on that branch). `at_zero()` is the one definition of "zero" for
every metric, and `zero_band_axis()` gives zeros a band of their own below a
log axis, so both halves of Figure 5 show them the same way.

MARGINS. `fit_bottom_margin()` sets a figure's bottom margin from the drawn
extent of its labels rather than from a guessed fraction (2026-10-10; the
guessed fractions left 5-6 mm of empty figure under Figure 5A, 5B and 5D).

TICK LABELS (2026-10-14). The log axes are labelled 0.001, 0.01, 0.1, 1, 10
(`decade_label`), not 10^-3 ... 10^1. The power form is mathtext, which the
SVG carries in DejaVu Sans with 4.9 pt exponents; it reached the assembled
Figure 5 as a second typeface under the 6 pt floor.

Dependencies: numpy, pandas, scipy; matplotlib for zero_band_axis(),
mark_axis_break() and fit_bottom_margin() only, imported inside them or used
through the figure passed in.
"""

import logging
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


# ---------------------------------------------------------------------------
# Tree
# ---------------------------------------------------------------------------

# Written UNROOTED -- a trifurcation at the base. A rooted string makes codeml
# fit two unidentifiable root branches; see the module docstring.
#
# PAML numbers tips by their order in the SEQUENCE file, not the tree string.
# branch_dnds_four.py writes the alignment As, Pu, Bsc, Tox (MAFFT keeps the
# input order), which is what TIP_NUMBER assumes and is also the order here.
TREE_TOPOLOGY = "((Ascaris,Parascaris),Baylisascaris,Toxocara);"

# Tree tip labels, in the order panels are drawn.
SPECIES_ORDER = ['Ascaris', 'Parascaris', 'Baylisascaris', 'Toxocara']

# Tree tip label -> paper_colors species key, so SPECIES_LABEL gives the
# italic display name and SPECIES gives the colour. The tree needs plain ASCII
# tip names (PAML matches them literally); the figure needs italics.
SPECIES_KEY = {'Ascaris': 'ascaris', 'Parascaris': 'parascaris',
               'Baylisascaris': 'baylisascaris', 'Toxocara': 'toxocara'}

TIP_NUMBER = {'1': 'Ascaris', '2': 'Parascaris',
              '3': 'Baylisascaris', '4': 'Toxocara'}

METRICS = ('omega', 'dN', 'dS')

# Display names for the metrics. House rule: dN, dS and dN/dS in anything a
# reader sees; omega only as a model parameter in an equation. Plain text and
# upright, in the label typeface; paper_colors.house_text() would rewrite the
# old mathtext spellings anyway, but the source should say what is printed.
METRIC_LABEL = {'omega': 'Branch dN/dS',
                'dN': 'Branch dN',
                'dS': 'Branch dS'}


def root_degree(newick):
    """Number of children of the root of a Newick string (3 = unrooted, 2 = rooted)."""
    s = newick.strip()
    if not s.startswith('('):
        raise ValueError(f'not a Newick tree: {newick!r}')
    depth, commas = 0, 0
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
            if depth == 0:
                break
        elif ch == ',' and depth == 1:
            commas += 1
    if depth != 0:
        raise ValueError(f'unbalanced parentheses in {newick!r}')
    return commas + 1


def check_unrooted(newick=TREE_TOPOLOGY):
    """Raise ValueError if `newick` has a basal bifurcation; return the root degree.

    codeml with clock = 0 fits a rooted tree as given, with two root branches of
    which only the sum is estimable. Call this before starting any codeml job.
    """
    k = root_degree(newick)
    if k < 3:
        raise ValueError(
            f'{newick} is rooted (its root has {k} children). codeml with '
            'clock = 0 would fit both root branches although only their sum is '
            'estimable. Write the tree with a trifurcation at the base, e.g. '
            '((Ascaris,Parascaris),Baylisascaris,Toxocara);')
    return k


# ---------------------------------------------------------------------------
# codeml free-ratio output
# ---------------------------------------------------------------------------

def parse_branch_table(outfile):
    """Extract the 'dN & dS for each branch' table from a codeml output file.

    Returns {branch_string: {'omega': f, 'dN': f, 'dS': f}}, or {} if the table
    is absent. Raises nothing -- the caller decides what a short table means.
    """
    try:
        with open(outfile) as fh:
            lines = fh.readlines()
    except OSError:
        return {}

    in_table = False
    branch_data = {}
    for line in lines:
        if 'dN & dS for each branch' in line:
            in_table = True
            continue
        if not in_table:
            continue
        s = line.strip()
        if s.startswith('branch'):
            continue                      # header row
        if not s:
            if branch_data:               # blank line after the rows ends it
                break
            continue
        parts = s.split()
        if len(parts) >= 7 and '..' in parts[0]:
            try:
                branch_data[parts[0]] = {'omega': float(parts[4]),
                                         'dN': float(parts[5]),
                                         'dS': float(parts[6])}
            except (ValueError, IndexError):
                pass
    return branch_data


def descendant_tips(branch_data, n_tips=4):
    """Map every node to the set of tips reachable from it, away from its parent.

    Built from codeml's own branch list, so it survives a change in how PAML
    numbers internal nodes -- which is the whole reason terminal branches are
    identified by child tip number rather than by position in the table.
    """
    children = {}
    for branch in branch_data:
        nodes = branch.split('..')
        if len(nodes) != 2:
            continue
        children.setdefault(nodes[0], []).append(nodes[1])

    tips = {str(i) for i in range(1, n_tips + 1)}
    memo = {}

    def walk(node, seen):
        if node in memo:
            return memo[node]
        if node in tips:
            memo[node] = frozenset([node])
            return memo[node]
        if node in seen:                  # cycle guard; malformed input only
            return frozenset()
        out = set()
        for kid in children.get(node, []):
            out |= walk(kid, seen | {node})
        memo[node] = frozenset(out)
        return memo[node]

    return {node: walk(node, frozenset())
            for node in set(children) | set(sum(children.values(), []))}


def label_branches(branch_data, n_tips=4):
    """Split a codeml branch table into terminal branches and internal branches.

    Returns (terminal, internal) where

      terminal  {species_name: {'omega','dN','dS'}}
      internal  [{'branch','subtends','omega','dN','dS'}, ...]

    `subtends` is a sorted tuple of species names -- the smaller of the two
    sides of the split, so it is independent of which end codeml calls the
    parent. For the standard four-taxon tree it is
    ('Ascaris', 'Parascaris'), and the branch is the As/Pu-versus-Bsc/Tox
    split, NOT an ancestral ascarid branch. See the module docstring.

    Raises ValueError unless the table has exactly 2n - 3 branches, the count
    for an unrooted tree of n tips. One more than that is the signature of a
    rooted tree: codeml splits the root branch in two, and the two halves --
    for this tree the Toxocara branch and a second branch subtending only
    Toxocara -- take arbitrary shares of one length. Callers turn the error
    into a per-gene status.
    """
    expected = 2 * n_tips - 3
    if len(branch_data) != expected:
        raise ValueError(
            f'{len(branch_data)} branches in the codeml table; an unrooted '
            f'{n_tips}-taxon tree has {expected}'
            + (' -- codeml was given a rooted tree' if len(branch_data)
               == expected + 1 else ''))
    desc = descendant_tips(branch_data, n_tips)
    all_tips = {str(i) for i in range(1, n_tips + 1)}

    terminal, internal = {}, []
    for branch, data in sorted(branch_data.items()):
        nodes = branch.split('..')
        if len(nodes) != 2:
            continue
        child = nodes[1]
        if child in TIP_NUMBER:
            terminal[TIP_NUMBER[child]] = dict(data)
            continue
        below = desc.get(child, frozenset())
        other = all_tips - below
        side = below if len(below) <= len(other) else other
        internal.append({
            'branch': branch,
            'subtends': tuple(sorted(TIP_NUMBER[t] for t in side
                                     if t in TIP_NUMBER)),
            **data})
    return terminal, internal


def internal_branch_note(internal):
    """One-line description of what the internal branch(es) actually separate."""
    if not internal:
        return 'no internal branch recovered'
    if len(internal) == 1:
        b = internal[0]
        side = ', '.join(b['subtends'])
        rest = [s for s in SPECIES_ORDER if s not in b['subtends']]
        return (f"one internal branch ({b['branch']}), separating "
                f"[{side}] from [{', '.join(rest)}]. This is an UNROOTED "
                f"split, not an ancestral branch of either side.")
    return (f"{len(internal)} internal branches: "
            + '; '.join(f"{b['branch']} subtends [{', '.join(b['subtends'])}]"
                        for b in internal))


# ---------------------------------------------------------------------------
# Gene-level masking and comparison
# ---------------------------------------------------------------------------

def usable_mask(df, species, metric, omega_cap, apply_cap=True, ds_max=np.inf):
    """Boolean mask of genes usable for ONE metric on ONE branch.

    apply_cap=True: omega, dN and dS must all be finite on that branch, the
    branch's own omega must be <= omega_cap, and dS must be > 0. A gene dropped
    on one branch is still used on the others -- the dS ~ 0 artifact is
    branch-specific and dropping the whole gene would throw away three good
    estimates.

    Why dS > 0 (added 2026-10-10): on a branch with no substitution at all,
    codeml reports dN = dS = 0 and an omega that is wherever its optimiser
    stopped, because omega no longer affects the likelihood. codeml seeds its
    starting values from the clock, so the same gene came back with omega
    0.0001 in one test run and 0.79 or 3.48 in the next; whether such a gene
    passed the cap was a coin toss. dN > 0 with dS = 0 already gives omega =
    999 and was always masked, so this only removes the 0/0 genes, and removes
    them deterministically.

    apply_cap=False: only the requested metric must be finite. This is the
    correct set for testing dS, because the omega cap is itself a filter on dS
    (it removes genes with dS near zero, preferentially), so applying it before
    a dS test is circular.

    Negative values are always dropped: codeml writes -1 for an estimate it
    could not obtain, and that sentinel must not enter a median.
    """
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}, got {metric!r}")

    om = pd.to_numeric(df[f'{species}_omega'], errors='coerce')
    dn = pd.to_numeric(df[f'{species}_dN'], errors='coerce')
    ds = pd.to_numeric(df[f'{species}_dS'], errors='coerce')
    vals = {'omega': om, 'dN': dn, 'dS': ds}[metric]

    if apply_cap:
        keep = (om.notna() & dn.notna() & ds.notna() & (om <= omega_cap)
                & (ds > 0))
    else:
        keep = vals.notna()

    keep = keep & (vals >= 0)
    if metric == 'dS' and np.isfinite(ds_max):
        keep = keep & (ds <= ds_max)
    return keep.fillna(False).astype(bool)


def compare(vals_ret, vals_elim, min_n=5):
    """Mann-Whitney comparison returning medians, fold change and P."""
    out = {'ret_n': len(vals_ret), 'elim_n': len(vals_elim),
           'ret_median': vals_ret.median() if len(vals_ret) else np.nan,
           'elim_median': vals_elim.median() if len(vals_elim) else np.nan,
           'fold_change': np.nan, 'p_value': np.nan}
    if len(vals_ret) >= min_n and len(vals_elim) >= min_n:
        _, p = mannwhitneyu(vals_ret, vals_elim, alternative='two-sided')
        out['p_value'] = p
        rm = out['ret_median']
        if np.isfinite(rm) and rm > 0:
            out['fold_change'] = out['elim_median'] / rm
    return out


def split_and_compare(df, keep, region_col, col, min_n=5):
    """Apply a keep mask, split retained/eliminated, run compare()."""
    sub = df[keep]
    ret = pd.to_numeric(sub[sub[region_col] == 'retained'][col],
                        errors='coerce').dropna()
    elim = pd.to_numeric(sub[sub[region_col] == 'eliminated'][col],
                         errors='coerce').dropna()
    return compare(ret, elim, min_n), ret, elim


def bootstrap_median_ci(vals, n_boot=2000, seed=0, alpha=0.05):
    """Percentile bootstrap CI on a median.

    The previous bar panel used `Series.sem()` -- the standard error of the
    MEAN -- as an error bar on a median. That is a category error, not a
    conservative approximation: for a right-skewed distribution like branch
    omega the two are not even centred on the same quantity.
    """
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 5:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, v.size, size=(n_boot, v.size))
    meds = np.median(v[idx], axis=1)
    return (float(np.percentile(meds, 100 * alpha / 2)),
            float(np.percentile(meds, 100 * (1 - alpha / 2))))


# ---------------------------------------------------------------------------
# Zero on a log axis -- one definition for Figure 5A and 5C
# ---------------------------------------------------------------------------

#: codeml bounds the free-ratio dN/dS to [0.0001, 999] (SetxBound in codeml.c).
#: A branch with no nonsynonymous change is reported as dN/dS = 0.0001 with dN
#: printed as 0.0000, so 0.0001 is zero, not a small ratio. Testing `<= 0`
#: instead, as Figure 5A did until 2026-10-10, never fires.
CODEML_OMEGA_FLOOR = 1e-4


def at_zero(vals, metric):
    """Boolean array: estimates that mean no substitution of that class.

    dN/dS at codeml's floor (that is, dN = 0), or dN or dS printed as 0.
    """
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}, got {metric!r}")
    v = np.asarray(vals, dtype=float)
    if metric == 'omega':
        return v <= CODEML_OMEGA_FLOOR * (1 + 1e-9)
    return v <= 0


def as_plotted(vals, metric):
    """Values as drawn on a zero-band axis: every zero, including dN/dS at
    codeml's floor, set to exactly 0. Box statistics are taken on these, so a
    quartile that falls among the zeros is drawn at 0 rather than at 0.0001."""
    v = np.asarray(vals, dtype=float).copy()
    v[at_zero(v, metric)] = 0.0
    return v


def zero_band_limits(vals_list, metric, hi_cap=np.inf, phi=99.5, lo_min=1e-4):
    """Limits for zero_band_axis(): (lo, hi, n_above).

    lo is the decade at or below the smallest non-zero value, never below
    `lo_min` (codeml prints four decimals, so nothing non-zero is smaller), so
    no non-zero value falls into the zero band. hi is the decade at or above
    the `phi`-th percentile, capped at `hi_cap`'s decade. n_above counts the
    values above hi, which the axis clips and the legend must state.
    """
    v = (np.concatenate([np.asarray(a, dtype=float) for a in vals_list
                         if len(a)]) if vals_list else np.array([]))
    v = v[np.isfinite(v)]
    nz = v[~at_zero(v, metric)] if v.size else v
    if not nz.size:
        return lo_min, lo_min * 1e3, 0
    lo = max(lo_min, 10.0 ** np.floor(np.log10(nz.min())))
    top = min(float(np.percentile(nz, phi)), hi_cap)
    hi = max(10.0 ** np.ceil(np.log10(top)), lo * 10)
    return lo, hi, int((nz > hi).sum())


def decade_label(k, decimal=True):
    """Tick label for 10**k.

    decimal=True (the default since 2026-10-14) writes the plain number --
    0.001, 0.01, 0.1, 1, 10 -- in the label font. The power-of-ten form is
    mathtext: matplotlib writes it into the SVG in DejaVu Sans with the
    exponent at 70% of the tick size (4.9 pt at 7 pt, under the 6 pt floor),
    and Arial has no superscript minus to replace it with. audit() cannot see
    the exponent's size, because it is part of one mathtext string.
    """
    k = int(k)
    if decimal:
        return f'{10.0 ** k:.{max(0, -k)}f}'
    return r'$\mathdefault{10^{%d}}$' % k


def zero_band_axis(ax, lo, hi, band=0.3, label_lo=False, decimal=True):
    """Log y-axis from `lo` to `hi` with a band for 0 below it.

    symlog with linthresh = lo: [0, lo] is linear and gets `band` decades of
    height, everything above lo is logarithmic. Zero then sits inside the band
    instead of below the bottom of a log axis, so a box or whisker that reaches
    0 is drawn to 0, and a pile of zeros shows as a pile. The tick at lo is
    left unlabelled by default: it sits one band-height above "0", too close for
    two labels at 7 pt. Mark the scale change with mark_axis_break().
    Ticks are labelled as plain numbers (decade_label); decimal=False restores
    the 10^k form.
    Returns the axes-fraction height of lo, where the break belongs.
    """
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

    ax.set_yscale('symlog', linthresh=lo, linscale=band, base=10)
    ax.set_ylim(-0.6 * lo, hi)
    k0, k1 = int(round(np.log10(lo))), int(round(np.log10(hi)))
    ax.yaxis.set_major_locator(
        FixedLocator([0.0] + [10.0 ** k for k in range(k0, k1 + 1)]))

    def fmt(y, _pos):
        if y == 0:
            return '0'
        k = int(round(np.log10(y)))
        if k == k0 and not label_lo:
            return ''
        return decade_label(k, decimal)

    ax.yaxis.set_major_formatter(FuncFormatter(fmt))
    ax.yaxis.set_minor_locator(FixedLocator(
        [m * 10.0 ** k for k in range(k0, k1) for m in range(2, 10)]))
    ax.yaxis.set_minor_formatter(NullFormatter())
    to_axes = ax.transAxes.inverted()
    return float(to_axes.transform(ax.transData.transform((0, lo)))[1])


def mark_axis_break(ax, y, size=6.0):
    """Two short slants across the left spine at axes-fraction `y`, with the
    spine cut between them: the scale is not continuous there."""
    from matplotlib.path import Path

    gap = Path([(0, -0.5), (0, 0.5)], [Path.MOVETO, Path.LINETO])
    slants = Path([(-1, -1.0), (1, -0.2), (-1, 0.2), (1, 1.0)],
                  [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO])
    kw = dict(transform=ax.transAxes, clip_on=False, linestyle='none',
              markerfacecolor='none')
    ax.plot([0], [y], marker=gap, markersize=size * 0.55, zorder=9,
            markeredgecolor='white', markeredgewidth=2.0, **kw)
    ax.plot([0], [y], marker=slants, markersize=size, zorder=10,
            markeredgecolor='black', markeredgewidth=0.6, **kw)


def fit_bottom_margin(fig, pad_mm=1.0, passes=3, tol_mm=0.05):
    """Move the subplots' bottom edge so the lowest drawn text sits `pad_mm`
    above the bottom of the figure. Returns the change in mm (positive: empty
    space reclaimed; negative: text was cut off and the margin grew).

    A bottom margin given as a figure fraction is a guess at how tall the tick
    labels will be. Figure 5A (two-line tick labels plus a species label) and
    5B and 5D (one-line species labels) each left 5-6 mm of empty figure under
    the labels, which shows as a gap between stacked panels at assembly. This
    measures instead: it draws the figure, takes the tight bounding box of
    everything on it, and moves `bottom` by the difference.

    The figure keeps its size -- the panel footprint planned for assembly --
    and the axes grow into the reclaimed height. Nothing is scaled, so type
    stays at the house sizes. Positions given in axes fractions (the zero-band
    break, the stars above the panels) are unaffected, because an axes
    fraction does not depend on the axes' physical height. Repeats up to
    `passes` times, because a taller axes can move y tick labels.

    Call after subplots_adjust() and before audit() and save(). Requires a
    figure whose layout engine, if any, honours subplots_adjust() -- in
    practice layout='none', which every Figure 5 script uses. Constrained and
    tight layout ignore subplots_adjust(), so the margin would silently stay.
    """
    engine = getattr(fig, 'get_layout_engine', lambda: None)()
    if engine is not None and not getattr(engine, 'adjust_compatible', True):
        raise ValueError('fit_bottom_margin needs layout="none"; this figure '
                         f'uses {type(engine).__name__}, which ignores '
                         'subplots_adjust()')
    h_mm = fig.get_figheight() * 25.4
    moved = 0.0
    for _ in range(passes):
        fig.canvas.draw()
        # Figure.get_tightbbox() is in inches.
        lowest_mm = fig.get_tightbbox(fig.canvas.get_renderer()).y0 * 25.4
        excess_mm = lowest_mm - pad_mm
        if abs(excess_mm) < tol_mm:
            break
        sp = fig.subplotpars
        bottom = sp.bottom - excess_mm / h_mm
        if not 0.0 < bottom < sp.top:
            raise ValueError(f'fit_bottom_margin: bottom would be {bottom:.3f} '
                             f'with top at {sp.top:.3f}; the figure is too '
                             'short for its labels')
        fig.subplots_adjust(bottom=bottom)
        moved += excess_mm
    return moved


# ---------------------------------------------------------------------------
# Annotated table
# ---------------------------------------------------------------------------

def load_annotated(path, groupings=None):
    """Read and validate `branch_dnds_four_annotated.tsv`.

    Returns (DataFrame, groupings). Exits loudly rather than half-working: a
    missing rate column produced a KeyError deep inside a plotting function
    before, which reads as a plotting bug rather than a wrong input file.
    """
    if not os.path.isfile(path):
        raise SystemExit(f"ERROR: annotated table not found: {path}")

    df = pd.read_csv(path, sep='\t')
    logging.info('Loaded %s genes from %s', f'{len(df):,}',
                 os.path.basename(path))

    missing = [f'{sp}_{m}' for sp in SPECIES_ORDER for m in METRICS
               if f'{sp}_{m}' not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing rate columns: {missing}")

    # A table from the rooted tree (runs before 2026-10-10) records the second
    # root branch as "internal", subtending Toxocara alone, and its Toxocara
    # values are an arbitrary share of the real branch. Refuse it outright.
    if 'internal_subtends' in df.columns:
        splits = set(df['internal_subtends'].dropna().astype(str)) - {''}
        if splits and splits != {'Ascaris|Parascaris'}:
            raise SystemExit(
                f"ERROR: {os.path.basename(path)} records internal split(s) "
                f"{sorted(splits)}; the unrooted tree gives Ascaris|Parascaris. "
                "This table was made with the rooted tree (before 2026-10-10), "
                "so its Toxocara branch values are not estimates. Rerun "
                "branch_dnds_four.py and use its new annotated table.")

    if groupings:
        bad = [g for g in groupings if g not in df.columns]
        if bad:
            raise SystemExit(f"ERROR: grouping column(s) not in table: {bad}")
    else:
        groupings = [g for g in ('as_region', 'tox_region') if g in df.columns]
    if not groupings:
        raise SystemExit("ERROR: no region column found (looked for as_region, "
                         "tox_region). Pass --groupings explicitly.")

    for g in groupings:
        counts = df[g].value_counts().to_dict()
        logging.info('  %s: %s', g, counts)
        if not {'retained', 'eliminated'} <= set(counts):
            logging.warning('  %s lacks one of retained/eliminated -- panels '
                            'for this grouping will be empty', g)
    return df, groupings
