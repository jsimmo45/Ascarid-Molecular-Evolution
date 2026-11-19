#!/usr/bin/env python3
"""
rbh_quality.py -- quality of the six reciprocal-best-hit ortholog sets
(supplementary fig. S18; S10 before the 2026-10-12 renumbering).

Reads what reciprocal_best_hits.py already wrote -- no BLAST is rerun -- and
draws the six pairs in one figure. The May figure was one two-panel file per
pair, three of them from before the Baylisascaris identifier fix.

Inputs, per pair (--pair LABEL RBH_TABLE [RATE_TABLE])
-----------------------------------------------------
  RBH_TABLE   rbh_full.tsv: gene1, gene2, pident_fwd/rev, bitscore_fwd/rev
              (after the script's identity >= 30% and length >= 50 aa filters).
              rbh_stats.txt beside it, if present, supplies the unfiltered count.
  RATE_TABLE  optional kaks_filtered.tsv of the same pair. When given, the
              figure shows the ortholog pairs that entered the rate analysis
              (the n Methods quotes) rather than every RBH, and the pair's
              median dS fixes its colour. With a rate table for every pair, the
              pairs are ordered by median dS and coloured with pair_ramp(), the
              rule Figure 1C uses, so each pair has the same colour in both
              figures. Otherwise the command-line order is used and the
              summary says so.

Panels (one 185 mm strip, figS10_rbh_quality)
---------------------------------------------
  A  Protein identity of each ortholog pair (mean of the forward and reverse
     BLAST identity), one stepped histogram per pair, as % of that pair's
     orthologs in 1% bins. The two clusters are named on the panel; each name
     sits at the lowest place over the body of its own curves where it clears
     every curve (label_groups()). The x-axis starts at the floor of the
     lowest drawn identity, so nothing is cut (reciprocal_best_hits.py keeps
     identity >= 30%).
  B  Reciprocity: log2(forward / reverse bitscore) per pair, horizontal boxes
     (quartiles; whiskers 5th-95th percentiles), with the number of pairs in a
     column to the right. Values near 0 mean the two searches agreed. The
     ascarid boxes are often too narrow to see, so the summary gives every
     pair's quartiles and whiskers for the legend.
  Panel letters A and B sit at the top of the canvas, each level with the
  left edge of its panel's labels (--no-panel-letters to omit them).

Outputs (in --outdir)
---------------------
  figS10_rbh_quality.{svg,pdf,png}   (stem kept from the frozen numbering;
                                      the assembly copies it as FigS18)
  rbh_quality_per_pair.tsv     counts at each step, medians, IQRs, and the
                               quartiles and whiskers of panel B
  rbh_quality_summary.txt
  rbh_quality.log
  provenance.json

Usage
-----
  python rbh_quality.py \\
    --pair 'Ascaris × P. univalens' rbh_results/rbh_full.tsv kaks_results/kaks_filtered.tsv \\
    --pair 'Ascaris × B. schroederi' rbh_As_Bsc_v2/rbh_full.tsv kaks_As_Bsc_v2/kaks_filtered.tsv \\
    ... --outdir results/rbh_quality_<stamp>

Dependencies: numpy, pandas, matplotlib; paper_colors, region_lib and
rate_tables from lib/.
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
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
from matplotlib.transforms import ScaledTranslation

# A script that cannot import the shared style must stop, not fall back to
# matplotlib defaults.
from paper_colors import (NEUTRAL, FONT, LINEWIDTH, apply_style, audit, figure,
                          save, mm, pair_ramp)

LOG2_RATIO_LABEL = 'log₂(forward / reverse bitscore)'
from region_lib import fingerprint, git_commit, setup_logging
from rate_tables import read_rate_table

_HERE = os.path.dirname(os.path.abspath(__file__))

RBH_COLS = ('gene1', 'gene2', 'pident_fwd', 'pident_rev', 'bitscore_fwd',
            'bitscore_rev')


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def raw_rbh_count(rbh_path):
    """Unfiltered RBH count from rbh_stats.txt beside the table, or None."""
    stats = os.path.join(os.path.dirname(os.path.abspath(rbh_path)),
                         'rbh_stats.txt')
    if not os.path.isfile(stats):
        return None
    with open(stats) as fh:
        for line in fh:
            m = re.match(r'\s*RBH pairs \(raw\):\s*([\d,]+)', line)
            if m:
                return int(m.group(1).replace(',', ''))
    return None


def load_pair(label, rbh_path, rate_path):
    """One pair: the RBH table (restricted to the analysed pairs when a rate
    table is given) and its counts."""
    rbh = pd.read_csv(rbh_path, sep='\t')
    miss = [c for c in RBH_COLS if c not in rbh.columns]
    if miss:
        sys.exit(f'ERROR: {rbh_path} lacks {miss}; is it rbh_full.tsv from '
                 f'reciprocal_best_hits.py? Columns: {list(rbh.columns)}')
    rbh['gene1'] = rbh['gene1'].astype(str)
    rbh['gene2'] = rbh['gene2'].astype(str)
    info = {'label': label, 'rbh_table': rbh_path, 'rate_table': rate_path,
            'n_rbh_raw': raw_rbh_count(rbh_path), 'n_rbh': len(rbh),
            'n_rate_rows': np.nan, 'n_analysed': np.nan,
            'n_rate_not_in_rbh': np.nan, 'median_dS': np.nan,
            'orientation': ''}
    use = rbh
    if rate_path:
        rate = read_rate_table(rate_path, require=('dS',))
        info['n_rate_rows'] = len(rate)
        ds = pd.to_numeric(rate['dS'], errors='coerce')
        info['median_dS'] = float(np.nanmedian(ds)) if ds.notna().any() else np.nan
        if {'gene1', 'gene2'} <= set(rate.columns):
            pairs = set(zip(rate['gene1'].astype(str), rate['gene2'].astype(str)))
            fwd = [(a, b) in pairs for a, b in zip(rbh['gene1'], rbh['gene2'])]
            rev = [(b, a) in pairs for a, b in zip(rbh['gene1'], rbh['gene2'])]
            n_f, n_r = int(np.sum(fwd)), int(np.sum(rev))
            keep = np.array(fwd if n_f >= n_r else rev)
            info['orientation'] = 'as listed' if n_f >= n_r else 'swapped'
            use = rbh[keep]
            info['n_analysed'] = int(keep.sum())
            info['n_rate_not_in_rbh'] = len(pairs) - int(keep.sum())
            if info['n_rate_not_in_rbh'] > 0.01 * len(pairs):
                logging.warning('  %s: %s rate-table pairs are not in the RBH '
                                'table -- are these the same run?', label,
                                f"{info['n_rate_not_in_rbh']:,}")
        else:
            logging.warning('  %s: rate table has no gene1/gene2 columns; all '
                            'RBH pairs are drawn', label)
    use = use.assign(
        identity=(use['pident_fwd'] + use['pident_rev']) / 2,
        log2_bitscore_ratio=np.log2(use['bitscore_fwd'] / use['bitscore_rev']))
    logging.info('  %-28s RBH %s%s, analysed %s, median dS %s', label,
                 f"{info['n_rbh']:,}",
                 f" (raw {info['n_rbh_raw']:,})" if info['n_rbh_raw'] else '',
                 f"{info['n_analysed']:,}" if np.isfinite(info['n_analysed'])
                 else 'n/a',
                 f"{info['median_dS']:.3f}" if np.isfinite(info['median_dS'])
                 else 'n/a')
    return use, info


def ratio_box(lr):
    """Panel B's box and whiskers for one pair, from the finite values the
    panel draws (np.percentile, as matplotlib's boxplot): median, quartiles
    (the box), 5th and 95th percentiles (the whiskers)."""
    v = lr.replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
    keys = ('log2_ratio_n_drawn', 'log2_ratio_median', 'log2_ratio_q1',
            'log2_ratio_q3', 'log2_ratio_iqr', 'log2_ratio_p5',
            'log2_ratio_p95')
    if not len(v):
        return dict.fromkeys(keys, np.nan)
    p5, q1, med, q3, p95 = np.percentile(v, [5, 25, 50, 75, 95])
    return dict(zip(keys, (len(v), float(med), float(q1), float(q3),
                           float(q3 - q1), float(p5), float(p95))))


# ---------------------------------------------------------------------------
# Figure helpers (shared pattern with the other converted scripts)
# ---------------------------------------------------------------------------

def blank_figure(width, height_mm):
    """A journal-width canvas with no axes and no layout engine; panels are
    placed explicitly, because constrained_layout shrinks a panel to nothing
    to make room for text it cannot fit and only warns."""
    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax)
    return fig


def drawn_texts(fig):
    """Every text object actually drawn (tick labels outside the view limits
    are returned by findobj but never drawn, so they are excluded)."""
    out = list(fig.texts)
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            tol = 1e-9 * max(abs(hi - lo), 1e-300)
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                loc = tick.get_loc()
                if loc is None or not (lo - tol <= loc <= hi + tol):
                    continue
                out += [lab for lab in (tick.label1, tick.label2)
                        if lab.get_visible()]
            out.append(axis.label)
            out.append(axis.get_offset_text())
        out += list(ax.texts)
        leg = ax.get_legend()
        if leg is not None:
            out += list(leg.get_texts())
    return [t for t in out if t.get_visible() and (t.get_text() or '').strip()]


def text_problems(fig, tol_px=1.0):
    """Drawn text running past the canvas edge, or overlapping other drawn
    text. Neither is checked by audit(), and both save without complaint."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    out = []
    texts = drawn_texts(fig)
    boxes = [t.get_window_extent(rend) for t in texts]
    for t, bb in zip(texts, boxes):
        if (bb.x0 < -tol_px or bb.y0 < -tol_px or bb.x1 > W + tol_px
                or bb.y1 > H + tol_px):
            out.append(f'text {t.get_text()[:28]!r} runs past the canvas edge')
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = boxes[i], boxes[j]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) > tol_px and
                    min(a.y1, b.y1) - max(a.y0, b.y0) > tol_px):
                out.append(f'text {texts[i].get_text()[:20]!r} overlaps '
                           f'{texts[j].get_text()[:20]!r}')
    return out


def place_row(fig, axes, shares, gaps_mm, pad_mm=0.8, iterations=4,
              top_extra_px=0.0):
    """Lay out panels side by side from their measured text; the row fills
    the canvas whatever font the node renders with. `top_extra_px` is kept
    free above everything the axes draw (the panel-letter row)."""
    W, H = fig.bbox.width, fig.bbox.height
    pad = mm(pad_mm) * fig.dpi
    gaps = [mm(g) * fig.dpi for g in gaps_mm]
    shares = np.asarray(shares, float) / np.sum(shares)
    for _ in range(iterations):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        tight = [ax.get_tightbbox(r) for ax in axes]
        boxes = [ax.get_window_extent(r) for ax in axes]
        lefts = [max(0.0, b.x0 - t.x0) for t, b in zip(tight, boxes)]
        rights = [max(0.0, t.x1 - b.x1) for t, b in zip(tight, boxes)]
        top = max(max(0.0, t.y1 - b.y1) for t, b in zip(tight, boxes))
        bottom = max(max(0.0, b.y0 - t.y0) for t, b in zip(tight, boxes))
        avail = W - 2 * pad - sum(lefts) - sum(rights) - sum(gaps)
        widths = avail * shares
        x = pad
        y0 = pad + bottom
        hh = H - 2 * pad - bottom - top - top_extra_px
        for i, ax in enumerate(axes):
            x += lefts[i]
            ax.set_position([x / W, y0 / H, widths[i] / W, hh / H])
            x += widths[i] + rights[i] + (gaps[i] if i < len(gaps) else 0.0)


def finish(fig, outdir, stem):
    """Audit, check text, save in every format, close."""
    problems = list(audit(fig, verbose=False)) + text_problems(fig)
    for prob in problems:
        logging.warning('  audit %s: %s', stem, prob)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths, problems


def add_panel_letters(fig, letters):
    """Bold panel letters, created before the layout so their height can be
    reserved; align_panel_letters() puts them in place afterwards."""
    return [fig.text(0.0, 1.0, s, ha='left', va='top',
                     fontsize=FONT['panel_letter'], fontweight='bold',
                     color=NEUTRAL['ink']) for s in letters]


def letter_room_px(fig, letters, gap_mm=0.8):
    """Height to keep free above the panels for the letter row."""
    if not letters:
        return 0.0
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    return (max(t.get_window_extent(r).height for t in letters)
            + mm(gap_mm) * fig.dpi)


def align_panel_letters(fig, axes, letters, pad_mm=0.8):
    """Each letter at the top of the canvas, level with the left edge of its
    panel's drawn extent (tick labels and axis label included). Call after
    the final layout. Returns the letters' x positions in mm, for the log."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    pad = mm(pad_mm) * fig.dpi
    out = []
    for ax, t in zip(axes, letters):
        x0 = max(ax.get_tightbbox(r).x0, pad)
        t.set_position((x0 / W, 1.0 - pad / H))
        out.append(round(x0 / fig.dpi * 25.4, 2))
    return out


def side_text(ax, y, s, dx_pt=3.0, **kw):
    """Text in a column just right of the axes, at data height `y`."""
    tr = ax.get_yaxis_transform() + ScaledTranslation(dx_pt / 72, 0,
                                                      ax.figure.dpi_scale_trans)
    opts = dict(ha='left', va='center_baseline', fontsize=FONT['annotation'],
                color=NEUTRAL['ink'], clip_on=False)
    opts.update(kw)
    return ax.text(1.0, y, s, transform=tr, **opts)


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

def draw_identity(ax, data, colors, lo=None):
    """A: identity histogram per pair, % of that pair's orthologs per 1% bin.

    No legend: the colour column in B keys each pair. The cluster names are
    added by label_groups() once the layout is final. Returns the panel
    facts and the drawn curves (bin edges, one array of % per pair).
    """
    allv = np.concatenate([d['identity'].to_numpy(float) for d in data])
    lo = np.floor(min(allv.min(), 100) / 5) * 5 if lo is None else lo
    if allv.min() < lo:                  # a fixed `lo` must not cut data
        logging.warning('  A: %d identities below the x-axis start %g%% are '
                        'not drawn', int((allv < lo).sum()), lo)
    edges = np.arange(lo, 101, 1.0)
    ymax = 0.0
    curves = []
    for d, c in zip(data, colors):
        v = d['identity'].to_numpy(float)
        h, _ = np.histogram(v, bins=edges)
        pct = 100 * h / max(len(v), 1)
        ymax = max(ymax, pct.max())
        curves.append(pct)
        try:
            ax.stairs(pct, edges, color=c, lw=LINEWIDTH['data'], zorder=3)
        except AttributeError:                       # matplotlib < 3.4
            ax.step(edges, np.append(pct, pct[-1]), where='post', color=c,
                    lw=LINEWIDTH['data'], zorder=3)
    ax.set_xlim(lo, 100)
    ax.set_ylim(0, ymax * 1.18)
    ax.set_xlabel('Protein identity of ortholog pair (%)')
    ax.set_ylabel('Ortholog pairs (% per 1% bin)')
    return ({'x_from': float(lo), 'min_identity': float(allv.min()),
             'n_below_axis': int((allv < lo).sum())}, edges, curves)


def label_groups(ax, edges, curves, groups, gap_pt=2.0, min_frac=0.5,
                 step=0.25):
    """Name each cluster of curves in panel A where the name clears every
    curve, its own group's and the other's.

    Every centre on a `step` grid is tried that keeps the name inside the
    x-range, over the body of its own group (the group's highest curve under
    the name reaches at least `min_frac` of the group's peak) and where its
    own group is the top curve at the centre (so the name is read with its
    own curves, not the other group's). The name's bottom goes `gap_pt`
    above the highest curve of ANY group under its width (plus a pad of
    gap_pt either side), and the lowest such place is kept, ties going to
    the one nearest the group's peak. Names already placed are obstacles for
    the next. The y-limit is raised if a name would not fit.

    Call after the final layout: a name's width and height in data units
    depend on the axes size. Returns the text artists (so a later layout
    pass can remove and redo them) and, per group, where the name went.
    """
    fig = ax.figure
    if not groups:
        return [], {}
    order = list(dict.fromkeys(groups))            # first-appearance order
    centres = (edges[:-1] + edges[1:]) / 2
    env_all = np.max(np.vstack(curves), axis=0)
    texts, placed = [], {}
    for _ in range(4):                             # ylim may grow; redo
        for t in texts:
            t.remove()
        texts, placed, boxes = [], {}, []
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        bb = ax.get_window_extent(r)
        (x_lo, x_hi), (y_lo, y_hi) = ax.get_xlim(), ax.get_ylim()
        px_x = (x_hi - x_lo) / bb.width            # data units per pixel
        px_y = (y_hi - y_lo) / bb.height
        gap_px = gap_pt / 72 * fig.dpi
        need_top = y_hi
        for g in order:
            mine = [c for c, gg in zip(curves, groups) if gg == g]
            rest = [c for c, gg in zip(curves, groups) if gg != g]
            own = np.max(np.vstack(mine), axis=0)
            other = (np.max(np.vstack(rest), axis=0) if rest
                     else np.zeros_like(own))
            t = ax.text(0, 0, g, ha='center', va='bottom',
                        fontsize=FONT['annotation'], color=NEUTRAL['ink'],
                        zorder=6)
            ext = t.get_window_extent(r)
            half_w = ext.width / 2 * px_x
            pad_x, h = gap_px * px_x, ext.height * px_y
            peak_x = float(centres[int(np.argmax(own))])
            best = None
            for xc in np.arange(x_lo + half_w, x_hi - half_w + 1e-9, step):
                under = (edges[1:] > xc - half_w - pad_x) & \
                        (edges[:-1] < xc + half_w + pad_x)
                if not under.any() or own[under].max() < min_frac * own.max():
                    continue
                k = int(np.clip(np.searchsorted(edges, xc, side='right') - 1,
                                0, len(own) - 1))
                if own[k] < other[k]:              # other group on top here
                    continue
                y = float(env_all[under].max()) + gap_px * px_y
                for (bx0, bx1, by0, by1) in boxes:   # names already placed
                    if bx0 < xc + half_w + pad_x and bx1 > xc - half_w - pad_x:
                        y = max(y, by1 + gap_px * px_y)
                cand = (round(y, 9), abs(xc - peak_x), float(xc), y)
                best = cand if best is None or cand < best else best
            if best is None:                       # no body wide enough
                under = (edges[1:] > peak_x - half_w) & \
                        (edges[:-1] < peak_x + half_w)
                best = (0, 0, peak_x,
                        float(env_all[under].max()) + gap_px * px_y)
            xc, y = best[2], best[3]
            t.set_position((xc, y))
            texts.append(t)
            boxes.append((xc - half_w, xc + half_w, y, y + h))
            placed[g] = {'x_centre': round(xc, 2), 'y_bottom': round(y, 3),
                         'own_peak_x': peak_x}
            need_top = max(need_top, y + h + 0.03 * (y_hi - y_lo))
        if need_top <= y_hi + 1e-12:
            break
        ax.set_ylim(y_lo, need_top)
    return texts, placed


def draw_reciprocity(ax, data, labels, colors, clip_pct=(0.5, 99.5),
                     box_alpha=0.40):
    """B: log2(forward / reverse bitscore) per pair, horizontal boxes."""
    vals = [d['log2_bitscore_ratio'].replace([np.inf, -np.inf], np.nan)
            .dropna().to_numpy(float) for d in data]
    ys = -np.arange(len(vals), dtype=float)
    kw = dict(positions=ys, widths=0.6, patch_artist=True, showfliers=False,
              whis=(5, 95), zorder=2)
    try:
        bp = ax.boxplot(vals, orientation='horizontal', **kw)
    except TypeError:                     # matplotlib < 3.10
        bp = ax.boxplot(vals, vert=False, **kw)
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(to_rgba(c, box_alpha))
        patch.set_edgecolor('black')
        patch.set_linewidth(LINEWIDTH['thin'])
    for part in ('whiskers', 'caps'):
        for art in bp[part]:
            art.set_color('black')
            art.set_linewidth(LINEWIDTH['thin'])
    for med in bp['medians']:
        med.set_color('black')
        med.set_linewidth(LINEWIDTH['data'])
    ax.axvline(0, color=NEUTRAL['rule'], lw=LINEWIDTH['rule'], ls=(0, (4, 2)),
               zorder=1)
    # x-range from the whiskers, symmetric about 0 so agreement reads as centre
    w = max(np.abs(np.percentile(v, [5, 95])).max() for v in vals if len(v))
    lim = max(w * 1.25, 0.05)
    ax.set_xlim(-lim, lim)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3, symmetric=True))
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _: f'{v:g}'.replace('-', '\u2212')))
    ax.set_ylim(ys[-1] - 0.6, ys[0] + 0.6)
    ax.yaxis.set_major_locator(FixedLocator(ys))
    ax.set_yticklabels(labels)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel(LOG2_RATIO_LABEL)
    # Right-hand column: the pair's colour, which keys panel A, then its n.
    for y, v, c in zip(ys, vals, colors):
        tr = ax.get_yaxis_transform() + ScaledTranslation(
            8.0 / 72, 0, ax.figure.dpi_scale_trans)
        ax.plot([1.0], [y], ls='none', marker='_', ms=9, mew=2.6, color=c,
                transform=tr, clip_on=False, zorder=5)
        side_text(ax, y, f'n = {len(v):,}', dx_pt=15.0)
    return {'x_limit': float(lim)}


def plot_s18(data, infos, colors, outdir, stem, args):
    """Supplementary fig. S18 (stem figS10_rbh_quality): the two panels in
    one strip, laid out from their measured text; cluster names placed on
    the final axes; panel letters level with each panel's left edge."""
    fig = blank_figure('double', args.height_mm)
    ax_a = fig.add_axes([0.07, 0.2, 0.40, 0.6])
    ax_b = fig.add_axes([0.62, 0.2, 0.28, 0.6])
    labels = [i['label'] for i in infos]
    groups = [args.outgroup_label if args.outgroup and args.outgroup in lab
              else args.ingroup_label for lab in labels] if args.outgroup \
        else None
    info_a, edges, curves = draw_identity(ax_a, data, colors)
    info_b = draw_reciprocity(ax_b, data, labels, colors)
    letters = ([] if args.no_panel_letters
               else add_panel_letters(fig, ['A', 'B']))
    extra = letter_room_px(fig, letters)
    axes = [ax_a, ax_b]
    place_row(fig, axes, shares=[0.52, 0.48], gaps_mm=[10.0],
              top_extra_px=extra)
    # The names can raise A's y-limit, which can change its tick labels and
    # so the layout; lay out again and redo the names on the final axes.
    names, _ = label_groups(ax_a, edges, curves, groups)
    place_row(fig, axes, shares=[0.52, 0.48], gaps_mm=[10.0],
              top_extra_px=extra)
    for t in names:
        t.remove()
    names, placed = label_groups(ax_a, edges, curves, groups)
    letter_x = align_panel_letters(fig, axes, letters)
    info_a.update({'names': placed, 'y_top': float(ax_a.get_ylim()[1])})
    paths, problems = finish(fig, outdir, stem)
    return paths, problems, {'A': info_a, 'B': info_b,
                             'letters_x_mm': letter_x}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Quality of the reciprocal-best-hit ortholog sets '
                    '(supplementary fig. S18; S10 before renumbering).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--pair', nargs='+', action='append', required=True,
                   metavar='LABEL RBH_TABLE [RATE_TABLE]',
                   help='One pair: a label, its rbh_full.tsv and optionally its '
                        'rate table. Repeat for every pair.')
    p.add_argument('--outgroup', default='T. canis',
                   help='Pairs whose label contains this are named as the '
                        'outgroup cluster in panel A (default "T. canis"; '
                        'empty string for no cluster names)')
    p.add_argument('--ingroup-label', default='Ascarid pairs')
    p.add_argument('--outgroup-label', default='Outgroup pairs')
    p.add_argument('--height-mm', type=float, default=58.0)
    p.add_argument('--no-panel-letters', action='store_true',
                   help='Omit the A and B panel letters (drawn by default, '
                        'since the legend refers to both panels)')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    for spec in args.pair:
        if len(spec) not in (2, 3):
            p.error(f'--pair takes LABEL RBH_TABLE [RATE_TABLE], got {spec}')
    if len(args.pair) > 6:
        p.error('at most six pairs (pair_ramp has six steps)')
    if args.outdir is None:
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.outdir = os.path.join('results', f'rbh_quality_{stamp}')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'rbh_quality')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('rbh_quality.py -- supplementary fig. S18 (S10 before '
                 'renumbering)')
    for spec in args.pair:
        for f in spec[1:]:
            if not os.path.isfile(f):
                sys.exit(f'ERROR: not found: {f}  (pair {spec[0]!r})')
    apply_style()

    loaded = [load_pair(s[0], s[1], s[2] if len(s) == 3 else None)
              for s in args.pair]
    data = [d for d, _ in loaded]
    infos = [i for _, i in loaded]
    have_ds = all(np.isfinite(i['median_dS']) for i in infos)
    if have_ds:
        order = list(np.argsort([i['median_dS'] for i in infos], kind='stable'))
        order_note = 'ordered by median dS (as Figure 1C)'
    else:
        order = list(range(len(infos)))
        order_note = ('command-line order: not every pair has a rate table, so '
                      'the colours may not match Figure 1C')
        logging.warning('  %s', order_note)
    data = [data[k] for k in order]
    infos = [infos[k] for k in order]
    colors = pair_ramp(len(infos)) if len(infos) >= 2 else [NEUTRAL['ink']]

    rows = []
    for d, i, c in zip(data, infos, colors):
        idv, lr = d['identity'], d['log2_bitscore_ratio']
        rows.append({**{k: i[k] for k in ('label', 'n_rbh_raw', 'n_rbh',
                                          'n_rate_rows', 'n_analysed',
                                          'n_rate_not_in_rbh', 'median_dS',
                                          'orientation')},
                     'n_drawn': len(d), 'colour': c,
                     'median_identity': float(idv.median()),
                     'identity_q1': float(idv.quantile(0.25)),
                     'identity_q3': float(idv.quantile(0.75)),
                     'pct_identity_below_50': float(100 * (idv < 50).mean()),
                     'median_abs_log2_bitscore_ratio': float(lr.abs().median()),
                     'pct_bitscores_within_5pct': float(
                         100 * (lr.abs() <= np.log2(1.05)).mean()),
                     **ratio_box(lr),
                     'rbh_table': i['rbh_table'], 'rate_table': i['rate_table']})
    per_pair = pd.DataFrame(rows)
    per_pair.to_csv(os.path.join(args.outdir, 'rbh_quality_per_pair.tsv'),
                    sep='\t', index=False, float_format='%.6g')

    logging.info('S18: six-pair strip')
    paths, problems, fig_info = plot_s18(data, infos, colors, args.outdir,
                                         'figS10_rbh_quality', args)

    # --- summary ----------------------------------------------------------------
    def n_txt(v):
        return f'{int(v):,}' if v is not None and np.isfinite(v) else 'n/a'
    L = ['Reciprocal-best-hit ortholog quality (supplementary fig. S18; '
         'S10 before renumbering)',
         '=' * 72, '', f'Pairs {order_note}.', '',
         f"{'pair':<30} {'RBH raw':>8} {'RBH':>7} {'analysed':>9} {'drawn':>7} "
         f"{'med dS':>7} {'med id %':>8} {'IQR':>12} {'|log2 bs| med':>13}"]
    for r in per_pair.itertuples():
        L.append(f"{r.label:<30} {n_txt(r.n_rbh_raw):>8} {r.n_rbh:>7,} "
                 f"{n_txt(r.n_analysed):>9} {r.n_drawn:>7,} "
                 f"{r.median_dS:>7.3f} {r.median_identity:>8.1f} "
                 f"{r.identity_q1:>5.1f}-{r.identity_q3:<6.1f} "
                 f"{r.median_abs_log2_bitscore_ratio:>13.4f}")
    L += ['',
          'RBH raw = before the identity >= 30% / length >= 50 aa filters '
          '(rbh_stats.txt);',
          'RBH = rbh_full.tsv; analysed = RBH pairs present in the rate table '
          '(the n Methods quotes);',
          'drawn = what the figure shows (analysed when a rate table is given, '
          'else RBH).', '']
    if per_pair['n_analysed'].notna().all():
        L.append('Analysed pairs, for Methods: ' + '; '.join(
            f"{r.label} {int(r.n_analysed):,}" for r in per_pair.itertuples()))
        stray = per_pair[per_pair['n_rate_not_in_rbh'] > 0]
        if len(stray):
            L.append('Rate-table pairs missing from the RBH table: ' + '; '.join(
                f"{r.label} {int(r.n_rate_not_in_rbh):,}"
                for r in stray.itertuples()))
        L.append('')
    # Panel B's boxes, for the legend: the ascarid boxes are too narrow to
    # see at the panel's scale, so the legend quotes them.
    L += ['Panel B, log2(forward / reverse bitscore) as drawn (for the legend)',
          f"{'pair':<30} {'n':>7} {'median':>9} {'IQR (box)':>21} "
          f"{'IQR width':>10} {'5th-95th (whiskers)':>21}"]
    for r in per_pair.itertuples():
        L.append(f"{r.label:<30} {n_txt(r.log2_ratio_n_drawn):>7} "
                 f"{r.log2_ratio_median:>9.4f} "
                 f"{r.log2_ratio_q1:>10.4f} to {r.log2_ratio_q3:<7.4f} "
                 f"{r.log2_ratio_iqr:>10.4f} "
                 f"{r.log2_ratio_p5:>10.4f} to {r.log2_ratio_p95:<7.4f}")
    widest = per_pair.loc[per_pair['log2_ratio_iqr'].idxmax()]
    L += [f"Widest IQR: {widest['label']}, {widest['log2_ratio_iqr']:.4f} "
          f"(a {100 * (2 ** widest['log2_ratio_iqr'] - 1):.2f}% bitscore "
          f"difference); x-axis +/-{fig_info['B']['x_limit']:.3f}.", '']
    fa = fig_info['A']
    names = '; '.join(f"{g} centred at {v['x_centre']:g}% identity"
                      for g, v in fa.get('names', {}).items())
    L += ['Drawn (figS10_rbh_quality; copied as FigS18 at assembly)',
          f"  A: x-axis from {fa['x_from']:g}% identity; lowest drawn identity "
          f"{fa['min_identity']:.1f}% (reciprocal_best_hits.py keeps identity "
          f">= 30%), so nothing is cut ({fa['n_below_axis']} below the axis). "
          "Colours keyed by the column in B, light to dark = increasing median "
          f"dS. y-axis to {fa['y_top']:.2f}.",
          f"  A cluster names: {names}." if names else
          '  A cluster names: none (--outgroup empty).',
          f"  B: x-axis +/-{fig_info['B']['x_limit']:.3f}; boxes = quartiles, "
          'whiskers = 5th-95th percentiles; dashed line = 0 (agreement).',
          ('  Panel letters A, B at x = ' + ', '.join(
              f'{v:g}' for v in fig_info['letters_x_mm']) + ' mm.')
          if fig_info['letters_x_mm'] else '  Panel letters: omitted.',
          '  AUDIT: ' + '; '.join(problems) if problems else '  audit: clean', '']
    with open(os.path.join(args.outdir, 'rbh_quality_summary.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(),
                   'python': platform.python_version(),
                   'versions': {'numpy': np.__version__,
                                'pandas': pd.__version__,
                                'matplotlib': matplotlib.__version__},
                   'argv': sys.argv,
                   'inputs': {i['label']: {'rbh': fingerprint(i['rbh_table']),
                                           'rate': fingerprint(i['rate_table'])}
                              for i in infos},
                   'order': order_note,
                   'figure': [os.path.basename(p) for p in paths],
                   'drawn': fig_info,
                   'audit': problems},
                  fh, indent=2, default=str)
    logging.info('Summary: %s', os.path.join(args.outdir,
                                             'rbh_quality_summary.txt'))
    logging.info('Done.')


if __name__ == '__main__':
    main()
