#!/usr/bin/env python3
"""
plain_text.py - a raised exponent as plain text, for figures that must open in
Illustrator in one typeface.

The house rule (script_to_figure_map.md, "text is plain", corrected 2026-10-15):
no mathtext and no Unicode sub- or superscript digits. Matplotlib writes
mathtext glyph by glyph in its math font (DejaVu Sans), and Arial has no
superscript digits, so Illustrator sets 10⁻²⁷ in Myriad Pro and will not switch
it back. The exponent is therefore its own text object, in the label typeface
at 6 pt (the journal floor), raised beside the run it belongs to: the fix made
by hand to Figure 6D on 2026-10-15. `region_lib.fmt_p` still writes mathtext;
`draw_p` here replaces it on a figure.

Joined at an anchor, not by measured width
------------------------------------------
Arial is not installed on the Isaac compute nodes, so matplotlib lays text out
there in DejaVu Sans, which is wider; the SVG still names Arial, and
Illustrator sets it in Arial. Two texts placed side by side by measured width
would open with a gap between them. Instead the run before the junction is
right-aligned and the run after it left-aligned at the same point, which
matplotlib's SVG backend writes as text-anchor end and start at one x, so the
two meet in any typeface. Only the position of the label as a whole (its centre
or edge) depends on widths, and those are taken from Helvetica's metrics, which
matplotlib ships and Arial shares, so the label lands where intended once
Illustrator sets it in Arial. A character Helvetica lacks (Greek) is measured
in the node's font.

Functions
---------
    p_parts(p)                       ('P = 6.1 × 10', '−27') or ('P = 0.034', None)
    draw_exponent_text(ax, x, y, base, exp, ...)
                                     `base` and a raised `exp`, joined at an anchor
    draw_p(ax, x, y, p, prefix='')   a P value in that form, after `prefix`
    draw_runs(ax, x, y, runs, ...)   a one-line label of runs that differ in size
                                     and baseline (a subscript), joined at anchors
    subscript_runs(before, sub, after, fontsize)
                                     the runs for "log2(...)" with a lowered 2
    replace_xlabel(ax, runs)         draw runs where the x label is, hide the label
    replace_ylabel(ax, runs)         the same for a left y label, read bottom to top
    width_pt(fig, s, size)           width of `s` as Arial sets it, in points

A subscript follows the exponent rule (2026-10-19, S12A's x label): its own
text at 6 pt, its baseline lowered SUB_DROP of the label's size below the
label's.

Rotated labels (2026-10-30, Fig. 3D's log10(dN/dS))
---------------------------------------------------
`draw_runs(..., rotation=90)` lays the runs out along the text direction, which
for a label read bottom to top is up the page: each run is drawn rotated about
its anchor (rotation_mode 'anchor'), the joins are made exactly as for a
horizontal label, and a run's rise moves it across the page instead of up it
(a lowered subscript moves right, towards the axis). `replace_ylabel` measures
the y label after layout and redraws it as runs on the label's own baseline.
"""

from __future__ import annotations

import os

import numpy as np
import matplotlib
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import ScaledTranslation

#: Exponent size in points: the journal floor, as in the 6D hand fix.
EXP_PT = 6.0
#: Raise of the exponent's baseline above the label's, in units of the label's
#: size (2.5 pt for a 7 pt label).
EXP_RAISE = 0.36
#: Subscript size in points, and the drop of its baseline below the label's, in
#: units of the label's size (1.2 pt for an 8 pt axis label).
SUB_PT = 6.0
SUB_DROP = 0.15
MINUS = '−'
TIMES = '×'

_HELVETICA = None          # (AFM, unicode -> glyph name), or False if unavailable


def _helvetica():
    """Helvetica's AFM metrics, shipped with matplotlib, loaded once."""
    global _HELVETICA
    if _HELVETICA is None:
        try:
            try:
                from matplotlib import _afm as afm_mod
            except ImportError:                  # matplotlib < 3.6
                from matplotlib import afm as afm_mod
            from fontTools.agl import UV2AGL     # a matplotlib dependency
            path = os.path.join(matplotlib.get_data_path(), 'fonts', 'afm',
                                'phvr8a.afm')
            with open(path, 'rb') as fh:
                _HELVETICA = (afm_mod.AFM(fh), UV2AGL)
        except Exception:                        # metrics are a refinement only
            _HELVETICA = False
    return _HELVETICA


def _measured(fig, s, size, weight='normal'):
    """(width, descent) of `s` in points, in the font the node actually has."""
    r = fig.canvas.get_renderer()
    w, _, d = r.get_text_width_height_descent(
        s, FontProperties(size=size, weight=weight), ismath=False)
    return w * 72.0 / fig.dpi, d * 72.0 / fig.dpi


def width_pt(fig, s, size, weight='normal'):
    """Advance width of `s` at `size` pt as Arial sets it.

    Helvetica's AFM widths, which Arial shares; a character Helvetica lacks is
    measured in the node's font. Bold, or no metrics available: measured.
    """
    h = _helvetica() if weight in ('normal', 'regular', 400) else False
    if not h:
        return _measured(fig, s, size, weight)[0]
    afm, agl = h
    total = 0.0
    for ch in s:
        try:
            total += afm.get_width_from_char_name(agl[ord(ch)]) * size / 1000.0
        except (KeyError, TypeError):
            total += _measured(fig, ch, size, weight)[0]
    return total


def p_parts(p, sig=3, digits=1):
    """A P value as (text, exponent).

    p >= 0.001: ('P = 0.034', None), `sig` significant digits (3, as
    region_lib.fmt_p). Smaller: ('P = 6.1 × 10', '−27'), mantissa to `digits`
    decimals, carried into the exponent if it rounds to 10. P = 0 (underflow):
    ('P < 10', '−300').
    """
    if p is None or not np.isfinite(p):
        return 'P = NA', None
    if p <= 0:
        return 'P < 10', MINUS + '300'
    if p >= 0.001:
        return f'P = {p:.{sig}g}', None
    e = int(np.floor(np.log10(p)))
    m = round(p / 10.0 ** e, digits)
    if m >= 10:
        m, e = m / 10.0, e + 1
    return f'P = {m:.{digits}f} {TIMES} 10', (MINUS if e < 0 else '') + str(abs(e))


def draw_exponent_text(ax, x, y, base, exp, transform=None, ha='center',
                       va='baseline', fontsize=7.0, exp_size=EXP_PT,
                       raise_frac=EXP_RAISE, dx_pt=0.0, dy_pt=0.0, **kw):
    """Draw `base` with `exp` raised after it, as two plain texts.

    (x, y) is in `transform` coordinates (default ax.transData), shifted by
    (dx_pt, dy_pt) points. `ha` aligns the whole label: 'left', 'center' or
    'right'. `va` is 'baseline', or 'bottom' to sit the label where a single
    text with va='bottom' would (matplotlib's line descent, from "lp").
    Further keywords (color, zorder, clip_on, fontweight) go to both texts.
    With no `exp`, one text is drawn. Returns the list of Text objects.
    """
    if ha not in ('left', 'center', 'right'):
        raise ValueError(f"ha must be 'left', 'center' or 'right', not {ha!r}")
    fig = ax.figure
    transform = ax.transData if transform is None else transform
    weight = kw.get('fontweight', 'normal')
    if va == 'bottom':
        dy_pt += max(_measured(fig, base, fontsize, weight)[1],
                     _measured(fig, 'lp', fontsize, weight)[1])
    elif va != 'baseline':
        raise ValueError(f"va must be 'baseline' or 'bottom', not {va!r}")

    def at(dx, dy):
        return transform + ScaledTranslation(dx / 72.0, dy / 72.0,
                                             fig.dpi_scale_trans)

    if not exp:
        return [ax.text(x, y, base, transform=at(dx_pt, dy_pt), ha=ha,
                        va='baseline', fontsize=fontsize, **kw)]
    w_base = width_pt(fig, base, fontsize, weight)
    w_exp = width_pt(fig, exp, exp_size, weight)
    # Offset of the junction from the anchor, for the alignment asked for.
    junction = {'left': w_base, 'center': (w_base - w_exp) / 2.0,
                'right': -w_exp}[ha]
    t_base = ax.text(x, y, base, transform=at(dx_pt + junction, dy_pt),
                     ha='right', va='baseline', fontsize=fontsize, **kw)
    t_exp = ax.text(x, y, exp,
                    transform=at(dx_pt + junction, dy_pt + raise_frac * fontsize),
                    ha='left', va='baseline', fontsize=exp_size, **kw)
    return [t_base, t_exp]


def draw_p(ax, x, y, p, prefix='', sig=3, **kw):
    """Draw `prefix` + a P value, with any exponent raised as its own text.

    Keywords as draw_exponent_text. Returns the list of Text objects.
    """
    base, exp = p_parts(p, sig=sig)
    return draw_exponent_text(ax, x, y, prefix + base, exp, **kw)


def draw_runs(ax, x, y, runs, transform=None, ha='center', va='baseline',
              fontsize=7.0, dx_pt=0.0, dy_pt=0.0, pivot=None, rotation=0,
              **kw):
    """Draw a one-line label made of runs that differ in size or baseline.

    `runs` is a list of (text, size_pt, rise_pt): size None means `fontsize`;
    rise_pt > 0 raises the run's baseline above the label's (an exponent),
    < 0 lowers it (a subscript). (x, y), `transform`, `ha`, `va`, `dx_pt` and
    `dy_pt` as draw_exponent_text; further keywords go to every text.

    Joined at anchors, as draw_exponent_text: the runs before `pivot` are
    right-aligned at their right edges, the pivot and the runs after it
    left-aligned at their left edges. The junction just before the pivot then
    meets in any typeface; every other junction rests on one run's Helvetica
    width, which Arial shares, so it is exact once Illustrator sets the text
    in Arial. The default pivot is the narrowest run, so the width a junction
    rests on is the smallest one: for a one-digit subscript, at most half a
    point off in the DejaVu Sans of a node without Arial. Returns the texts.

    `rotation` 90 draws a label read bottom to top (a y label): (x, y) is then
    a point on its baseline, `ha` is along the text (so up the page), each run
    is rotated about its anchor, and a rise moves a run across the page, left
    for an exponent and right for a subscript. Only va='baseline' is accepted
    with it; replace_ylabel finds the baseline.
    """
    if ha not in ('left', 'center', 'right'):
        raise ValueError(f"ha must be 'left', 'center' or 'right', not {ha!r}")
    if not runs:
        raise ValueError('draw_runs needs at least one run')
    if rotation not in (0, 90):
        raise ValueError(f'rotation must be 0 or 90, not {rotation!r}')
    if rotation and va != 'baseline':
        raise ValueError("a rotated label is placed by its baseline: "
                         "va must be 'baseline'")
    fig = ax.figure
    transform = ax.transData if transform is None else transform
    weight = kw.get('fontweight', 'normal')
    clean = []
    for run in runs:
        text, size, rise = run
        if not str(text):
            raise ValueError(f'empty run in {runs!r}')
        clean.append((str(text), float(fontsize if size is None else size),
                      float(rise)))
    if va == 'bottom':
        whole = ''.join(t for t, _, _ in clean)
        dy_pt += max(_measured(fig, whole, fontsize, weight)[1],
                     _measured(fig, 'lp', fontsize, weight)[1])
    elif va != 'baseline':
        raise ValueError(f"va must be 'baseline' or 'bottom', not {va!r}")
    widths = [width_pt(fig, t, size, weight) for t, size, _ in clean]
    start = {'left': 0.0, 'center': -sum(widths) / 2.0,
             'right': -sum(widths)}[ha]
    edges = start + np.concatenate([[0.0], np.cumsum(widths)[:-1]])
    if pivot is None:
        pivot = int(np.argmin(widths))
    if not 0 <= pivot < len(clean):
        raise ValueError(f'pivot {pivot} is not a run index (0-{len(clean) - 1})')

    def at(dx, dy):
        return transform + ScaledTranslation(dx / 72.0, dy / 72.0,
                                             fig.dpi_scale_trans)

    out = []
    for i, ((text, size, rise), edge, w) in enumerate(zip(clean, edges, widths)):
        off, align = (edge + w, 'right') if i < pivot else (edge, 'left')
        if rotation == 0:
            place, turn = at(dx_pt + off, dy_pt + rise), {}
        else:
            # Read bottom to top, the text runs up the page and its own "up"
            # points left: offsets along the text go to y, rises to -x.
            place = at(dx_pt - rise, dy_pt + off)
            turn = dict(rotation=90, rotation_mode='anchor')
        out.append(ax.text(x, y, text, transform=place, ha=align,
                           va='baseline', fontsize=size, **turn, **kw))
    return out


def subscript_runs(before, sub, after='', fontsize=8.0, sub_size=SUB_PT,
                   drop_frac=SUB_DROP):
    """Runs for draw_runs: `before`, then `sub` lowered at `sub_size`, then
    `after`. subscript_runs('Mean expression, log', '2', '(TPM + 0.1)') gives
    log2 with the 2 set as a subscript."""
    runs = [(before, None, 0.0), (sub, sub_size, -drop_frac * fontsize)]
    if after:
        runs.append((after, None, 0.0))
    return runs


def replace_xlabel(ax, runs, **kw):
    """Draw `runs` (draw_runs) where the x-axis label is, and hide the label.

    Set the plain label first, with the same words (log2 for log + subscript
    2): the layout measures it and reserves its room. Call once the layout is
    final, since the label's place depends on the tick labels. The runs take
    the label's size, colour, alignment and baseline, and are placed in axes
    coordinates plus a point offset, so they keep their place at any savefig
    dpi. Returns the texts.
    """
    fig = ax.figure
    lab = ax.xaxis.label
    text = lab.get_text()
    if not text.strip():
        raise ValueError('set the plain x label first, so the layout reserves '
                         'its room')
    if '\n' in text:
        raise ValueError('replace_xlabel handles a one-line label only')
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    size = lab.get_fontsize()
    weight = lab.get_fontweight()
    # The label is drawn va='top'; its box bottom is the baseline less the
    # line's descent, which matplotlib takes as the larger of the text's and
    # that of "lp".
    descent = max(_measured(fig, text, size, weight)[1],
                  _measured(fig, 'lp', size, weight)[1])
    box = lab.get_window_extent(rend)
    dy_pt = (box.y0 - ax.bbox.y0) * 72.0 / fig.dpi + descent
    opts = dict(color=lab.get_color(), clip_on=False)
    opts.update(kw)
    lab.set_visible(False)
    return draw_runs(ax, lab.get_position()[0], 0.0, runs,
                     transform=ax.transAxes,
                     ha=lab.get_horizontalalignment(), fontsize=size,
                     dy_pt=dy_pt, **opts)


def replace_ylabel(ax, runs, **kw):
    """Draw `runs` (draw_runs) where a left y-axis label is, and hide the label.

    The y-label counterpart of replace_xlabel (2026-10-30, for Fig. 3D's
    log10(dN/dS)), with the same contract: set the plain label first, with the
    same words, so the layout reserves its room, and call once the layout is
    final. Matplotlib draws a left y label rotated 90 degrees, read bottom to
    top, with its bottom towards the axis (va='bottom' in rotation mode
    'anchor'), so the right edge of its box is the line's bottom and the
    baseline lies one descent to the left of it. The runs are centred where the
    label was centred, take its size, colour and weight, and are placed in axes
    coordinates plus a point offset, so they keep their place at any savefig
    dpi. Returns the texts.
    """
    fig = ax.figure
    lab = ax.yaxis.label
    text = lab.get_text()
    if not text.strip():
        raise ValueError('set the plain y label first, so the layout reserves '
                         'its room')
    if '\n' in text:
        raise ValueError('replace_ylabel handles a one-line label only')
    if (ax.yaxis.get_label_position() != 'left'
            or round(float(lab.get_rotation())) % 360 != 90):
        raise ValueError('replace_ylabel handles a left y label read bottom to '
                         'top only')
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    size = lab.get_fontsize()
    weight = lab.get_fontweight()
    descent = max(_measured(fig, text, size, weight)[1],
                  _measured(fig, 'lp', size, weight)[1])
    box = lab.get_window_extent(rend)
    dx_pt = (box.x1 - ax.bbox.x0) * 72.0 / fig.dpi - descent
    y_mid = (0.5 * (box.y0 + box.y1) - ax.bbox.y0) / ax.bbox.height
    opts = dict(color=lab.get_color(), clip_on=False)
    if weight not in ('normal', 'regular', 400):
        opts['fontweight'] = weight
    opts.update(kw)
    lab.set_visible(False)
    return draw_runs(ax, 0.0, y_mid, runs, transform=ax.transAxes,
                     ha='center', fontsize=size, dx_pt=dx_pt, rotation=90,
                     **opts)
