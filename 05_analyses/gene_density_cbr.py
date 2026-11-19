#!/usr/bin/env python3
"""
gene_density_cbr.py - gene density as a function of signed distance to the
nearest chromosomal breakage region (CBR).

Figure 3B (Figure 4B before the renumbering of 2026-10-27). Replaces
`replot_gene_density.py`, which reported a density roughly 600x too high (see
"The units" below).

Writes:

    gene_density_cbr.{svg,pdf,png}       Figure 3B panel
    gene_density_resolutions.*           the same profile at three bin widths,
                                         as a binning-artefact control
    gene_density_binned.tsv              every plotted number, per bin
    gene_density_summary.txt             densities, rate ratio, diagnostics
    provenance.json                      inputs, arguments and versions


Changes 2026-10-30 (laid-out Fig. 3B; the hand edits of 2026-10-15 in code)
---------------------------------------------------------------------------
* The y label is plain text, "Gene density (genes per kb)". It was mathtext
  (genes kb$^{-1}$), which opened in Illustrator in a second typeface with
  uneven letter spacing.
* The "mean 0.057" label of the genome-wide rule sits in the lower right corner
  inside the axes, where it was placed by hand, clear of the trace and ribbon
  (_label_mean_corner(); the lower left corner if the right one is not clear,
  and the old on-the-rule placement if neither is).
* The top margin is 3.0 mm, as in the profile panel beside it
  (cbr_meta_profile.py), so the two axes tops are level when the panels are
  placed top-aligned. The panel had sat 2.3 mm lower in the 2026-10 layouts.
* The region labels are upright, as in Fig. 3D: italic is for species names.
* Densities, the binned table and the summary are unchanged.
* 2026-10-31: fontTools' INFO logging is silenced, as in cbr_meta_profile.py;
  it added several hundred lines per PDF to the SLURM log.


The units
---------
The earlier version divided the pooled gene count in a distance bin by the bin
width in kb and called the result "genes per kb", giving a plateau near 5 against
a genome-wide average of 0.062 genes/kb - about 600-fold too high, not the 100x
recorded in the figure map.

The arithmetic `counts / bin_width_kb` is dimensionally correct but the numerator
is wrong. `np.histogram` pools genes from *every* CBR onto one distance axis, so
a bin at +25 kb holds the genes that are 25 kb from CBR 1, plus those 25 kb from
CBR 2, and so on across all ~600 CBRs, while the denominator is the width of a
single window. The reported value is therefore inflated by the number of CBRs
that contribute to the bin. That factor is not a constant: it falls off with
distance as CBRs run out of chromosome, which distorts the *shape* of the profile
and not only its scale.

Dividing by a fixed CBR count would fix the scale and leave the shape wrong. The
correct denominator is the number of base pairs of genome that actually lie at
that signed distance from their own nearest CBR:

    density(bin) = 1000 * genes_with_signed_distance_in_bin
                        / bp_of_genome_with_signed_distance_in_bin

This is a genuine local gene density in genes/kb. Averaged over all bins it
returns the genome-wide density by construction, which is the check that the
scale is right, and it removes the size bias in the same operation - a bin far
out on the eliminated side is divided by the small number of base pairs that
really sit that deep inside an eliminated tract, rather than being compared
against a full-width window that mostly does not exist.

The base-pair denominator is measured by tiling each chromosome at `--tile-bp`
and passing every tile centre through *the same* signed-distance function used
for the genes. Numerator and denominator therefore share one definition of
"nearest CBR" and one sign convention, and any quirk in that definition cancels
between them.


Sign convention
---------------
Eliminated is positive and drawn right; retained is negative and drawn left.
This is the house convention, and `cbr_changepoint_relax.py` already signs
distance eliminated-positive - so Figure 3B and Figure 3D of the same figure
disagreed with each other until this flip.


What counts as a CBR
--------------------
The boundaries of the eliminated regions, with two corrections the earlier
version did not make:

  - Overlapping or bookended eliminated intervals are merged first, so an
    internal junction between two adjacent eliminated tracts is not counted as a
    CBR and does not put a spurious zero point in the middle of eliminated
    sequence.
  - A boundary that sits at a chromosome end is dropped (`--terminal-tol`). A
    terminal eliminated tract runs to the telomere; its distal edge is the end of
    the assembly, not a breakage region. Keeping it inflates the CBR count and
    puts a phantom retained side beyond the end of the chromosome.


What counts as a gene
---------------------
Overlapping gene models are merged into loci before counting, and the locus is
counted once at its midpoint. A miniprot query protein can align more than once
(Rank=1, Rank=2, ...) and the secondary alignments usually sit on or beside the
primary, so counting models rather than loci inflates density wherever the
projection is redundant. Pass `--no-merge-loci` to count raw models.


Usage
-----
    python gene_density_cbr.py \\
      --gff As_miniprot.gff \\
      --eliminated-bed AG.v50.eliminated_strict.bed \\
      --chrom-sizes AG.v50.fa.fai \\
      --species Ascaris \\
      --outdir results/gene_density_20260817

Dependencies: numpy, pandas, matplotlib, scipy, paper_colors
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
from scipy.stats import chi2

from paper_colors import (REGION, REGION_SHADE, NEUTRAL, apply_style, audit,
                          save, figure, panel_letter, mm)
from plain_text import width_pt
from region_lib import (parse_gff_models, parse_bed, merge_intervals,
                        load_chrom_sizes, derive_cbrs, in_intervals,
                        signed_distance, fingerprint, git_commit, setup_logging)

_HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def merge_models_to_loci(models):
    """Collapse overlapping gene models into loci.

    Returns (loci_df, diag). `loci_df` has chrom, start, end, midpoint and the
    number of models merged into each locus. Coordinates stay 1-based inclusive.
    """
    rows, n_overlapping, same_target, diff_target = [], 0, 0, 0
    for chrom, block in models.groupby('chrom', sort=False):
        b = block.sort_values(['start', 'end'])
        cur_start = cur_end = None
        targets, n_members = [], 0
        for r in b.itertuples():
            if cur_start is None:
                cur_start, cur_end = r.start, r.end
                targets, n_members = [r.target], 1
                continue
            if r.start <= cur_end:                      # overlaps the open locus
                n_overlapping += 1
                if r.target is not None and r.target in targets:
                    same_target += 1
                else:
                    diff_target += 1
                cur_end = max(cur_end, r.end)
                targets.append(r.target)
                n_members += 1
            else:
                rows.append((chrom, cur_start, cur_end, n_members))
                cur_start, cur_end = r.start, r.end
                targets, n_members = [r.target], 1
        if cur_start is not None:
            rows.append((chrom, cur_start, cur_end, n_members))

    loci = pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'n_models'])
    loci['midpoint'] = (loci['start'] + loci['end']) // 2
    diag = {'n_models': int(len(models)), 'n_loci': int(len(loci)),
            'n_models_overlapping': int(n_overlapping),
            'overlap_same_target': int(same_target),
            'overlap_diff_target': int(diff_target)}
    logging.info('  %s models merge into %s loci (%s models overlap a neighbour)',
                 f"{diag['n_models']:,}", f"{diag['n_loci']:,}",
                 f"{diag['n_models_overlapping']:,}")
    return loci, diag


# ---------------------------------------------------------------------------
# CBR positions and signed distance
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# The profile
# ---------------------------------------------------------------------------

def tile_bp_per_bin(chrom_sizes, cbrs_by_chrom, elim_merged, bin_edges,
                    tile_bp=1000):
    """Base pairs of genome at each signed distance, and how many CBRs supply them.

    Each sequence carrying at least one CBR is tiled at `tile_bp`; each tile
    centre is pushed through `signed_distance` and contributes its own length to
    the bin it lands in. Passing the tiles through the same function as the genes
    is the point of doing it this way: the numerator and the denominator cannot
    disagree about which CBR is nearest or which side of it a position is on.

    Resolution error is bounded by one tile, so keep `tile_bp` well under the bin
    width. Returns (bp per bin, distinct CBRs per bin, total bp considered).
    """
    n_bins = len(bin_edges) - 1
    bp = np.zeros(n_bins)
    cbr_sets = [set() for _ in range(n_bins)]
    total_bp = 0

    for chrom, cbrs in cbrs_by_chrom.items():
        L = chrom_sizes.get(chrom)
        if not L:
            logging.warning('  no length for %s; its base pairs are omitted from '
                            'the denominator', chrom)
            continue
        n_tiles = int(np.ceil(L / tile_bp))
        left = np.arange(n_tiles, dtype=float) * tile_bp          # 0-based offset
        centres = left + tile_bp / 2.0 + 0.5                      # -> 1-based
        widths = np.full(n_tiles, float(tile_bp))
        widths[-1] = L - left[-1]                                 # partial last tile
        centres[-1] = left[-1] + widths[-1] / 2.0 + 0.5
        total_bp += L

        sd, idx = signed_distance(chrom, centres, cbrs_by_chrom, elim_merged)
        b = np.digitize(sd, bin_edges) - 1
        keep = (b >= 0) & (b < n_bins)
        np.add.at(bp, b[keep], widths[keep])
        for bi, ci in zip(b[keep], idx[keep]):
            cbr_sets[bi].add((chrom, int(ci)))

    return bp, np.array([len(s) for s in cbr_sets]), total_bp


def poisson_ci(counts, alpha=0.05):
    """Exact (Garwood) Poisson interval for a count.

    A bin holds a count, not a mean of per-gene values, so its uncertainty is
    Poisson rather than a bootstrap over genes. The exact interval is used
    because the far bins have counts in the single digits, where the normal
    approximation gives a lower limit below zero.
    """
    counts = np.asarray(counts, dtype=float)
    lo = np.where(counts > 0, chi2.ppf(alpha / 2, 2 * counts) / 2.0, 0.0)
    hi = chi2.ppf(1 - alpha / 2, 2 * (counts + 1)) / 2.0
    return lo, hi


def build_profile(loci, cbrs_by_chrom, elim_merged, chrom_sizes,
                  window_kb=25, max_dist_kb=500, tile_bp=1000,
                  min_bp_per_bin=0.0, min_cbr_per_bin=1):
    """Binned gene density against signed CBR distance.

    Returns a frame with one row per bin: centre, gene count, base pairs,
    contributing CBRs, density and its Poisson interval. Bins that fail either
    coverage floor are kept in the table with `used = False` and are not drawn.
    """
    w = window_kb * 1000
    max_dist = max_dist_kb * 1000
    edges = np.arange(-max_dist, max_dist + w, w, dtype=float)
    centres = (edges[:-1] + edges[1:]) / 2.0

    # --- numerator: gene loci ---
    sd_all, kept_chroms, dropped = [], 0, 0
    for chrom, block in loci.groupby('chrom', sort=False):
        if chrom not in cbrs_by_chrom:
            dropped += len(block)
            continue
        kept_chroms += 1
        sd, _ = signed_distance(chrom, block['midpoint'].to_numpy(),
                                cbrs_by_chrom, elim_merged)
        sd_all.append(sd)
    signed = np.concatenate(sd_all) if sd_all else np.array([])
    logging.info('  %s loci on %d sequences carry a CBR; %s loci sit on '
                 'sequences with none and are excluded from both the count and '
                 'the base-pair denominator',
                 f'{len(signed):,}', kept_chroms, f'{dropped:,}')
    # A sequence-name mismatch between the GFF and the BED is not an error in
    # either file, so nothing upstream catches it: the script would find no CBR
    # on any sequence carrying genes and draw an empty panel rather than a wrong
    # one. Fail here instead, with both name spaces on screen.
    if len(signed) == 0:
        gff_names = sorted(loci['chrom'].unique())[:5]
        bed_names = sorted(cbrs_by_chrom)[:5]
        sys.exit('ERROR: not one gene locus sits on a sequence carrying a CBR. '
                 'The sequence names in the GFF and the BED almost certainly do '
                 'not match.\n'
                 f'  GFF sequences: {gff_names}\n'
                 f'  BED sequences: {bed_names}')
    if dropped > 4 * len(signed):
        logging.warning('  %.0f%% of loci sit on sequences with no CBR. Check '
                        'that this is the expected assembly and not a partial '
                        'name match.', 100.0 * dropped / (dropped + len(signed)))

    counts, _ = np.histogram(signed, bins=edges)

    # --- denominator: base pairs at each signed distance ---
    bp, n_cbr, total_bp = tile_bp_per_bin(chrom_sizes, cbrs_by_chrom,
                                          elim_merged, edges, tile_bp)

    with np.errstate(divide='ignore', invalid='ignore'):
        density = np.where(bp > 0, 1000.0 * counts / bp, np.nan)
        lo_c, hi_c = poisson_ci(counts)
        lo = np.where(bp > 0, 1000.0 * lo_c / bp, np.nan)
        hi = np.where(bp > 0, 1000.0 * hi_c / bp, np.nan)

    used = (bp >= min_bp_per_bin) & (n_cbr >= min_cbr_per_bin) & (counts >= 0) & (bp > 0)

    df = pd.DataFrame({
        'bin_start_kb': edges[:-1] / 1000.0,
        'bin_end_kb': edges[1:] / 1000.0,
        'centre_kb': centres / 1000.0,
        'side': np.where(centres > 0, 'eliminated', 'retained'),
        'n_genes': counts,
        'bp': bp,
        'n_cbr': n_cbr,
        'density_per_kb': density,
        'density_lo95': lo,
        'density_hi95': hi,
        'used': used,
    })
    return df, signed, total_bp


def whole_genome_density(loci, elim_merged, chrom_sizes, cbrs_by_chrom,
                         cbr_only=True):
    """Gene density in eliminated versus retained sequence, ignoring CBR distance.

    The binned profile only covers sequence within `--max-dist-kb` of a CBR, so
    its two sides are "eliminated near a CBR" and "retained near a CBR", not the
    two compartments as wholes. This is the unconditioned comparison: every base
    pair on the sequence, split by whether it is inside an eliminated tract.

    It is the number that answers "are eliminated regions gene-poor?" — the
    profile answers the narrower question of whether density *steps* at the
    boundary. Exact, not tiled: interval lengths are summed directly.
    """
    chroms = [c for c in chrom_sizes
              if (not cbr_only or c in cbrs_by_chrom)]
    total_bp = sum(chrom_sizes[c] for c in chroms)
    elim_bp = 0.0
    for c in chroms:
        starts, ends = elim_merged.get(c, (np.array([]), np.array([])))
        elim_bp += float((ends - starts + 1).sum()) if starts.size else 0.0
    ret_bp = total_bp - elim_bp

    n_elim = 0
    keep = loci[loci['chrom'].isin(chroms)]
    for c, block in keep.groupby('chrom', sort=False):
        starts, ends = elim_merged.get(c, (np.array([]), np.array([])))
        n_elim += int(in_intervals(block['midpoint'].to_numpy(float),
                                     starts, ends).sum())
    n_ret = int(len(keep)) - n_elim

    out = {'n_sequences': len(chroms), 'total_bp': total_bp,
           'retained': {'n_genes': n_ret, 'bp': ret_bp,
                        'density_per_kb': 1000.0 * n_ret / ret_bp if ret_bp else np.nan},
           'eliminated': {'n_genes': n_elim, 'bp': elim_bp,
                          'density_per_kb': 1000.0 * n_elim / elim_bp if elim_bp else np.nan}}
    if n_ret > 0 and n_elim > 0:
        rr = out['eliminated']['density_per_kb'] / out['retained']['density_per_kb']
        se = np.sqrt(1.0 / n_elim + 1.0 / n_ret)
        out['ratio'] = {'eliminated_over_retained': rr,
                        'lo95': rr * np.exp(-1.96 * se),
                        'hi95': rr * np.exp(1.96 * se)}
    else:
        out['ratio'] = None
    return out


def side_summary(profile, within_kb=None):
    """Pooled density each side of the CBR, and the eliminated:retained ratio.

    Pooling counts and base pairs before dividing - rather than averaging the
    per-bin densities - keeps the estimate unbiased when bins differ in coverage,
    which they do by construction on the eliminated side.
    """
    d = profile[profile['used']]
    if within_kb is not None:
        d = d[d['centre_kb'].abs() <= within_kb]
    out = {}
    for side in ('retained', 'eliminated'):
        s = d[d['side'] == side]
        n, b = float(s['n_genes'].sum()), float(s['bp'].sum())
        out[side] = {'n_genes': n, 'bp': b,
                     'density_per_kb': 1000.0 * n / b if b else np.nan}
    r, e = out['retained'], out['eliminated']
    if r['density_per_kb'] and e['n_genes'] > 0 and r['n_genes'] > 0:
        # Poisson rate-ratio interval on the log scale.
        rr = e['density_per_kb'] / r['density_per_kb']
        se = np.sqrt(1.0 / e['n_genes'] + 1.0 / r['n_genes'])
        out['ratio'] = {'eliminated_over_retained': rr,
                        'lo95': rr * np.exp(-1.96 * se),
                        'hi95': rr * np.exp(1.96 * se)}
    else:
        out['ratio'] = None
    return out


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _shade_regions(ax, xlo, xhi):
    """Pale region tracts behind the data: retained left, eliminated right."""
    ax.axvspan(xlo, 0, color=REGION_SHADE['retained'], zorder=0, lw=0)
    ax.axvspan(0, xhi, color=REGION_SHADE['eliminated'], zorder=0, lw=0)
    ax.axvline(0, color=NEUTRAL['cbr'], lw=1.0, zorder=4)


def _label_rule(ax, d, y, xlo, xhi, text, edge_frac=0.30):
    """Put the reference-rule label wherever the trace is furthest from the rule.

    The label sits on the rule, so a fixed corner puts it straight through the
    data whenever the trace happens to run at the genome mean there - which is
    exactly what a retained-side plateau does. The clearance in each outer
    `edge_frac` of the axis is measured and the roomier end wins, with the label
    pushed to whichever side of the rule the trace is not on.
    """
    span = xhi - xlo
    ends = {'left': d[d['centre_kb'] <= xlo + edge_frac * span],
            'right': d[d['centre_kb'] >= xhi - edge_frac * span]}
    best, best_gap, best_sign = 'right', -np.inf, 1
    for side, sub in ends.items():
        if sub.empty:
            continue
        dev = sub['density_per_kb'] - y
        gap = float(np.nanmin(np.abs(dev)))
        if gap > best_gap:
            best, best_gap = side, gap
            best_sign = -1 if float(np.nanmedian(dev)) > 0 else 1
    x = xlo if best == 'left' else xhi
    ax.annotate(text, xy=(x, y),
                xytext=(3 if best == 'left' else -3, 3 * best_sign),
                textcoords='offset points',
                ha='left' if best == 'left' else 'right',
                va='bottom' if best_sign > 0 else 'top',
                color=NEUTRAL['rule'], fontsize=plt.rcParams['xtick.labelsize'])


def _label_mean_corner(ax, d, text, show_ci=True, pad_pt=3.0, clear_pt=2.0):
    """Put the genome-wide-rule label in a lower corner inside the axes.

    Lower right first, where the 2026-10-15 hand edit put it, then lower left.
    A corner is taken only if the trace (and its ribbon, when drawn) stays at
    least `clear_pt` above the label's top over the label's width, measured in
    data units once the limits are final. Returns True if a corner was used;
    the caller falls back to the on-the-rule placement otherwise.
    """
    fig = ax.figure
    size = float(plt.rcParams['xtick.labelsize'])
    w_pt = width_pt(fig, text, size)
    fig.canvas.draw()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    bb = ax.get_window_extent()
    px_per_pt = fig.dpi / 72.0
    dx_per_pt = (x1 - x0) / (bb.width / px_per_pt)
    dy_per_pt = (y1 - y0) / (bb.height / px_per_pt)
    # Cap height plus a little for the label's own box (digits, no descenders).
    top = y0 + (pad_pt + 0.75 * size) * dy_per_pt
    low = d['density_lo95'] if show_ci else d['density_per_kb']
    for side in ('right', 'left'):
        if side == 'right':
            xa, xb = x1 - (pad_pt + w_pt) * dx_per_pt, x1 - pad_pt * dx_per_pt
        else:
            xa, xb = x0 + pad_pt * dx_per_pt, x0 + (pad_pt + w_pt) * dx_per_pt
        under = (d['centre_kb'] >= xa) & (d['centre_kb'] <= xb)
        # The trace is drawn as a line between bin centres, so take the bins
        # just outside the label's span as well.
        idx = np.flatnonzero(under.to_numpy())
        if idx.size:
            idx = np.arange(max(idx[0] - 1, 0), min(idx[-1] + 2, len(d)))
        else:
            near = int(np.argmin(np.abs(d['centre_kb'].to_numpy()
                                        - 0.5 * (xa + xb))))
            idx = np.array([near])
        lowest = float(np.nanmin(low.to_numpy()[idx]))
        if lowest >= top + clear_pt * dy_per_pt:
            ax.annotate(text, xy=(1.0 if side == 'right' else 0.0, 0.0),
                        xycoords='axes fraction',
                        xytext=(-pad_pt if side == 'right' else pad_pt, pad_pt),
                        textcoords='offset points',
                        ha='right' if side == 'right' else 'left', va='bottom',
                        color=NEUTRAL['rule'], fontsize=size)
            return True
    return False


def _region_labels(ax, xlo, xhi):
    """Name the two sides inside the axes, in the region colours. Upright, as
    in Fig. 3D: italic is for species names (2026-10-30)."""
    ax.text(xlo * 0.55, 0.94, 'Retained', transform=ax.get_xaxis_transform(),
            ha='center', va='top', color=REGION['retained'])
    ax.text(xhi * 0.55, 0.94, 'Eliminated', transform=ax.get_xaxis_transform(),
            ha='center', va='top', color=REGION['eliminated'])


def plot_profile(profile, genome_density, outdir, stem, species='',
                 width='single', height_mm=55.0, show_ci=True,
                 panel_letters=False, symmetric=False, max_dist_kb=None):
    """Figure 3B: one trace, Poisson ribbon, genome-wide reference rule.

    Layout is set with explicit margins rather than constrained_layout, which
    treats the axes as the flexible element and will silently shrink the panel
    toward zero to fit whatever else is on the canvas.
    """
    d = profile[profile['used']]
    if d.empty:
        logging.error('  every bin failed the coverage floors; nothing to plot')
        return None

    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax)
    w_in, h_in = fig.get_size_inches()
    # Top margin 3.0 mm, as cbr_meta_profile.plot_stack, so the two axes tops
    # are level when the panels are placed top-aligned (2026-10-30; was 3.5).
    left, right, bottom, top = mm(13.0), mm(3.0), mm(11.0), mm(3.0)
    ax = fig.add_axes([left / w_in, bottom / h_in,
                       (w_in - left - right) / w_in,
                       (h_in - bottom - top) / h_in])

    if symmetric and max_dist_kb:
        xlo, xhi = -max_dist_kb, max_dist_kb
    else:
        xlo, xhi = d['centre_kb'].min(), d['centre_kb'].max()
    _shade_regions(ax, xlo, xhi)

    if show_ci:
        ax.fill_between(d['centre_kb'], d['density_lo95'], d['density_hi95'],
                        color=NEUTRAL['ink'], alpha=0.18, lw=0, zorder=2)
    ax.plot(d['centre_kb'], d['density_per_kb'], color=NEUTRAL['ink'],
            lw=1.2, zorder=3)

    ax.axhline(genome_density, color=NEUTRAL['rule'], lw=0.8, ls=(0, (4, 2)),
               zorder=1)

    ax.set_xlim(xlo, xhi)
    ax.set_ylim(bottom=0)
    _region_labels(ax, xlo, xhi)
    ax.set_xlabel('Distance from nearest CBR (kb)')
    ax.set_ylabel('Gene density (genes per kb)')
    # The rule's label goes in once the limits are final: a lower corner if the
    # data leave room, else on the rule where the trace is furthest from it.
    mean_text = f'mean {genome_density:.3f}'
    if not _label_mean_corner(ax, d, mean_text, show_ci=show_ci):
        logging.info('  no clear lower corner for "%s"; placed on the rule',
                     mean_text)
        _label_rule(ax, d, genome_density, xlo, xhi, mean_text)
    if species:
        ax.set_title(species, fontstyle='italic', pad=2)
    if panel_letters:
        panel_letter(ax, 'B')

    for problem in audit(fig, verbose=False):
        logging.warning('  audit: %s', problem)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths


def plot_resolutions(profiles, genome_density, outdir, stem, species='',
                     width='single', height_mm=55.0, symmetric=False,
                     max_dist_kb=None):
    """The same profile at several bin widths - a binning-artefact control.

    Drawn in a light-to-dark neutral ramp rather than in colour: these are one
    quantity at three resolutions, not three categories, and blue and red are
    reserved for region.
    """
    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax)
    w_in, h_in = fig.get_size_inches()
    left, right, bottom, top = mm(13.0), mm(3.0), mm(11.0), mm(3.5)
    ax = fig.add_axes([left / w_in, bottom / h_in,
                       (w_in - left - right) / w_in,
                       (h_in - bottom - top) / h_in])

    greys = ['#BDBDBD', '#7A7A7A', NEUTRAL['ink']]
    xlo = xhi = 0.0
    finest = None
    for (wkb, prof), colour in zip(profiles, greys[-len(profiles):]):
        d = prof[prof['used']]
        if d.empty:
            continue
        finest = d if finest is None else finest
        xlo, xhi = min(xlo, d['centre_kb'].min()), max(xhi, d['centre_kb'].max())
        ax.plot(d['centre_kb'], d['density_per_kb'], color=colour, lw=1.0,
                label=f'{wkb:g} kb', zorder=3)

    if symmetric and max_dist_kb:
        xlo, xhi = -max_dist_kb, max_dist_kb
    _shade_regions(ax, xlo, xhi)
    ax.axhline(genome_density, color=NEUTRAL['rule'], lw=0.8, ls=(0, (4, 2)),
               zorder=1)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(bottom=0)
    _region_labels(ax, xlo, xhi)
    ax.set_xlabel('Distance from nearest CBR (kb)')
    ax.set_ylabel('Gene density (genes per kb)')
    if species:
        ax.set_title(species, fontstyle='italic', pad=2)
    ax.legend(frameon=False, loc='lower left', handlelength=1.4,
              handletextpad=0.5, borderaxespad=0.3, title='bin width',
              title_fontsize=plt.rcParams['legend.fontsize'])

    for problem in audit(fig, verbose=False):
        logging.warning('  audit: %s', problem)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Gene density against signed distance to the nearest CBR '
                    '(Figure 3B).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True,
                   help='Gene models. Every mRNA is used, not only the genes '
                        'under comparison elsewhere.')
    p.add_argument('--eliminated-bed', required=True,
                   help='Eliminated regions BED. Supplies both the eliminated '
                        'tracts and the CBR coordinates.')
    p.add_argument('--chrom-sizes', default=None,
                   help='.fai or two-column TSV of sequence lengths. Strongly '
                        'recommended: without it lengths are inferred from the '
                        'last annotated feature and the outermost retained bins '
                        'read slightly dense.')
    p.add_argument('--window-kb', type=float, default=25.0,
                   help='Bin width of the plotted profile (default: 25)')
    p.add_argument('--extra-windows', type=float, nargs='*',
                   default=[10.0, 50.0],
                   help='Additional bin widths for the resolution control '
                        '(default: 10 50). Pass with no values to skip it.')
    p.add_argument('--max-dist-kb', type=float, default=500.0,
                   help='Half-width of the distance axis (default: 500). The '
                        'base-pair denominator already removes the size bias, '
                        'so this is a legibility choice, not a correction.')
    p.add_argument('--tile-bp', type=int, default=1000,
                   help='Resolution of the base-pair denominator (default: '
                        '1000). Keep it well under --window-kb: the error in '
                        'the denominator is bounded by one tile.')
    p.add_argument('--min-bp-per-bin', type=float, default=0.0,
                   help='Drop bins backed by fewer base pairs than this')
    p.add_argument('--min-cbr-per-bin', type=int, default=1,
                   help='Drop bins to which fewer CBRs contribute (default: 1)')
    p.add_argument('--terminal-tol', type=int, default=1000,
                   help='Eliminated-tract boundaries within this many bp of a '
                        'sequence end are telomeric, not breakage regions, and '
                        'are dropped (default: 1000)')
    p.add_argument('--merge-loci', dest='merge_loci', action='store_true',
                   default=True,
                   help='Merge overlapping gene models into loci before '
                        'counting (default)')
    p.add_argument('--no-merge-loci', dest='merge_loci', action='store_false',
                   help='Count every gene model, including redundant miniprot '
                        'secondary alignments')
    p.add_argument('--no-ci', dest='show_ci', action='store_false',
                   help='Omit the Poisson ribbon')
    p.add_argument('--symmetric-axis', action='store_true',
                   help='Draw the full +/- --max-dist-kb axis. By default the '
                        'axis is trimmed to the bins that carry data, so the '
                        'eliminated side stops where the eliminated tracts run '
                        'out rather than trailing off into empty shading.')
    p.add_argument('--species', default='',
                   help='Species name for the panel title; omit for none')
    p.add_argument('--width', choices=['single', 'onehalf', 'double'],
                   default='single')
    p.add_argument('--height-mm', type=float, default=55.0)
    p.add_argument('--panel-letters', action='store_true',
                   help='Draw the panel letter. Off by default: panel letters '
                        'are added at assembly.')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if args.tile_bp >= args.window_kb * 1000:
        sys.exit(f'ERROR: --tile-bp {args.tile_bp} is not smaller than the '
                 f'{args.window_kb:g} kb bin width; the base-pair denominator '
                 f'would be no finer than the bins it fills.')
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'gene_density_cbr_' + datetime.now().strftime('%Y%m%d'))
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    # mode='w': FileHandler appends by default, so rerunning into an existing
    # directory concatenates logs while the figures are overwritten.
    setup_logging(args.outdir, 'gene_density_cbr')
    # matplotlib's PDF writer subsets fonts through fontTools, which logs every
    # table it prunes at INFO -- several hundred lines per figure, burying the
    # lines this script writes (as cbr_meta_profile.py). Warnings still show.
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.info('gene_density_cbr.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))

    apply_style()

    # --- inputs -----------------------------------------------------------
    models = parse_gff_models(args.gff)
    elim = parse_bed(args.eliminated_bed)
    chrom_sizes, inferred_sizes = load_chrom_sizes(args.chrom_sizes, models, elim)

    elim_merged = merge_intervals(elim)
    cbrs, cbr_diag = derive_cbrs(elim_merged, chrom_sizes, args.terminal_tol)
    if not cbrs:
        sys.exit('ERROR: no CBRs survive; check --eliminated-bed and '
                 '--terminal-tol.')

    if args.merge_loci:
        loci, locus_diag = merge_models_to_loci(models)
    else:
        loci = models.assign(midpoint=(models['start'] + models['end']) // 2,
                             n_models=1)
        locus_diag = {'n_models': int(len(models)), 'n_loci': int(len(models)),
                      'n_models_overlapping': None, 'overlap_same_target': None,
                      'overlap_diff_target': None}
        logging.info('  --no-merge-loci: counting all %s models',
                     f'{len(models):,}')

    # --- profiles ---------------------------------------------------------
    logging.info('Main profile at %g kb bins, +/- %g kb',
                 args.window_kb, args.max_dist_kb)
    profile, signed, total_bp = build_profile(
        loci, cbrs, elim_merged, chrom_sizes,
        window_kb=args.window_kb, max_dist_kb=args.max_dist_kb,
        tile_bp=args.tile_bp, min_bp_per_bin=args.min_bp_per_bin,
        min_cbr_per_bin=args.min_cbr_per_bin)

    # Reference density over the sequences that carry a CBR, not over the whole
    # assembly: sequences with no CBR contribute to neither side of the profile,
    # so including them would draw a reference rule the profile cannot be
    # compared against.
    ref_chroms = [c for c in loci['chrom'].unique() if c in cbrs]
    ref_bp = sum(chrom_sizes.get(c, 0) for c in ref_chroms)
    n_ref_loci = int(loci['chrom'].isin(ref_chroms).sum())
    genome_density = 1000.0 * n_ref_loci / ref_bp if ref_bp else np.nan
    logging.info('  reference density over CBR-carrying sequences: %.4f genes/kb '
                 '(%s loci over %.1f Mb, %d sequences)',
                 genome_density, f'{n_ref_loci:,}', ref_bp / 1e6, len(ref_chroms))

    drawn = profile[profile['used']]
    if drawn.empty:
        sys.exit('ERROR: no bin passes the coverage floors '
                 f'(--min-bp-per-bin {args.min_bp_per_bin:g}, '
                 f'--min-cbr-per-bin {args.min_cbr_per_bin}). Nothing to plot.')
    logging.info('  %d of %d bins pass the coverage floors; density range '
                 '%.4f - %.4f genes/kb', len(drawn), len(profile),
                 drawn['density_per_kb'].min(), drawn['density_per_kb'].max())

    # Scale check: the coverage-weighted mean of the drawn bins must land near
    # the genome-wide density. If it does not, the denominator is wrong.
    weighted = (1000.0 * drawn['n_genes'].sum() / drawn['bp'].sum())
    logging.info('  coverage-weighted mean of the profile: %.4f genes/kb '
                 '(%.2fx the reference)', weighted, weighted / genome_density)
    if not 0.3 < weighted / genome_density < 3.0:
        logging.warning('  the profile mean is far from the reference density; '
                        'the base-pair denominator is suspect')

    profile.to_csv(os.path.join(args.outdir, 'gene_density_binned.tsv'),
                   sep='\t', index=False, float_format='%.6g')

    extra = []
    for wkb in (args.extra_windows or []):
        logging.info('Control profile at %g kb bins', wkb)
        prof, _, _ = build_profile(
            loci, cbrs, elim_merged, chrom_sizes, window_kb=wkb,
            max_dist_kb=args.max_dist_kb, tile_bp=args.tile_bp,
            min_bp_per_bin=args.min_bp_per_bin,
            min_cbr_per_bin=args.min_cbr_per_bin)
        extra.append((wkb, prof))

    # --- figures ----------------------------------------------------------
    logging.info('Figure 3B')
    plot_profile(profile, genome_density, args.outdir, 'gene_density_cbr',
                 species=args.species, width=args.width,
                 height_mm=args.height_mm, show_ci=args.show_ci,
                 panel_letters=args.panel_letters,
                 symmetric=args.symmetric_axis, max_dist_kb=args.max_dist_kb)

    if extra:
        logging.info('Resolution control')
        series = sorted(extra + [(args.window_kb, profile)], key=lambda t: t[0])
        plot_resolutions(series, genome_density, args.outdir,
                         'gene_density_resolutions', species=args.species,
                         width=args.width, height_mm=args.height_mm,
                         symmetric=args.symmetric_axis,
                         max_dist_kb=args.max_dist_kb)

    # --- summary ----------------------------------------------------------
    overall = side_summary(profile)
    near = side_summary(profile, within_kb=100.0)
    whole = whole_genome_density(loci, elim_merged, chrom_sizes, cbrs,
                                 cbr_only=True)
    logging.info('Whole-sequence density (no distance conditioning): '
                 'retained %.4f, eliminated %.4f genes/kb',
                 whole['retained']['density_per_kb'],
                 whole['eliminated']['density_per_kb'])

    lines = ['Gene density at chromosomal breakage regions', '',
             f'Species              {args.species or "(unnamed)"}',
             f'Gene models          {locus_diag["n_models"]:,}',
             f'Gene loci counted    {locus_diag["n_loci"]:,}'
             + ('' if args.merge_loci else '  (--no-merge-loci)'),
             f'CBRs                 {cbr_diag["n_cbrs"]:,} on '
             f'{cbr_diag["n_sequences_with_cbrs"]} sequences',
             f'  terminal dropped   {cbr_diag["n_terminal_dropped"]:,} of '
             f'{cbr_diag["n_boundaries"]:,} tract boundaries',
             f'Bin width            {args.window_kb:g} kb',
             f'Denominator tile     {args.tile_bp:,} bp', '',
             f'Reference density    {genome_density:.4f} genes/kb  '
             f'({n_ref_loci:,} loci over {ref_bp/1e6:.1f} Mb on the '
             f'{len(ref_chroms)} sequences carrying a CBR)',
             f'Profile mean         {weighted:.4f} genes/kb  '
             f'({weighted / genome_density:.2f}x the reference)', '']

    lines.append('Whole sequence — no distance conditioning')
    lines.append('  the two compartments as wholes, not just the parts near a CBR')
    for side in ('retained', 'eliminated'):
        v = whole[side]
        lines.append(f'  {side:<12} {v["density_per_kb"]:.4f} genes/kb  '
                     f'({v["n_genes"]:,} genes over {v["bp"]/1e6:.1f} Mb)')
    if whole['ratio']:
        r = whole['ratio']
        lines.append(f'  ratio        {r["eliminated_over_retained"]:.3f} '
                     f'eliminated:retained [{r["lo95"]:.3f}, {r["hi95"]:.3f}]')
    lines.append('')

    for name, s in (('Whole axis', overall), ('Within 100 kb', near)):
        lines.append(f'{name}')
        for side in ('retained', 'eliminated'):
            v = s[side]
            lines.append(f'  {side:<12} {v["density_per_kb"]:.4f} genes/kb  '
                         f'({v["n_genes"]:,.0f} genes over {v["bp"]/1e6:.1f} Mb)')
        if s['ratio']:
            r = s['ratio']
            lines.append(f'  ratio        {r["eliminated_over_retained"]:.3f} '
                         f'eliminated:retained '
                         f'[{r["lo95"]:.3f}, {r["hi95"]:.3f}]')
        lines.append('')

    lines += [
        'Method', '',
        '  Density is genes per kb of sequence at that signed distance from its',
        '  own nearest CBR, not the pooled gene count divided by the bin width.',
        '  The denominator is measured by tiling each sequence and passing every',
        f'  tile centre through the same signed-distance function as the genes,',
        f'  at {args.tile_bp:,} bp resolution. Pooling over CBRs without this',
        '  inflates the density by the number of contributing CBRs, and by a',
        '  factor that falls with distance, so it distorts the shape and not',
        '  only the scale.', '',
        '  Positive is eliminated and drawn right; negative is retained and',
        '  drawn left, matching Figure 3D.', '',
        '  Ribbons are exact Poisson intervals on the bin count. A bin is a',
        '  count, so its uncertainty is Poisson rather than a bootstrap over',
        '  genes.', '']
    if inferred_sizes:
        lines += ['  Sequence lengths were inferred from the last annotated',
                  '  feature; every length is an underestimate and the outermost',
                  '  retained bins read slightly dense. Pass --chrom-sizes.', '']
    if args.merge_loci and locus_diag['n_models_overlapping']:
        lines += [f'  {locus_diag["n_models_overlapping"]:,} models overlap a '
                  f'neighbour and were merged',
                  f'  ({locus_diag["overlap_same_target"]:,} share a query '
                  f'protein, {locus_diag["overlap_diff_target"]:,} do not). A',
                  '  miniprot query can align more than once and the secondary',
                  '  alignments sit on or beside the primary, so counting models',
                  '  overstates density wherever the projection is redundant.', '']

    summary = '\n'.join(lines)
    with open(os.path.join(args.outdir, 'gene_density_summary.txt'), 'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(),
            'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv,
            'args': vars(args),
            'inputs': {'gff': fingerprint(args.gff),
                       'eliminated_bed': fingerprint(args.eliminated_bed),
                       'chrom_sizes': fingerprint(args.chrom_sizes)},
            'chrom_sizes_inferred': inferred_sizes,
            'cbrs': cbr_diag,
            'loci': locus_diag,
            'reference_density_per_kb': genome_density,
            'reference_sequences': len(ref_chroms),
            'reference_loci': n_ref_loci,
            'reference_bp': ref_bp,
            'profile_mean_per_kb': float(weighted),
            'side_summary': overall,
            'side_summary_within_100kb': near,
            'whole_sequence_density': whole,
        }, fh, indent=2, default=float)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
