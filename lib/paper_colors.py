"""
paper_colors.py - colour, type and page-geometry tokens for the manuscript figures.

Imported by every plotting script so that a design decision is made once and
applied everywhere:

    from paper_colors import REGION, SPECIES, apply_style, audit
    apply_style()

Colour is scoped
----------------
Blue and red encode genomic region -- retained and eliminated -- and nothing
else. Species, metrics, developmental stages and sub-categories draw from
separate scales that avoid those two hue families, so a reader who has learned
the region convention in one figure cannot misread it in another.

The scales are scoped rather than merged because they cannot be merged: six
mutually distinguishable categorical hues are not achievable under
red-green colour vision deficiency once blue and red are reserved. Region colour
and species colour therefore never appear together in one panel. Where both are
needed, region takes the colour and species identity is carried by the panel
title, tick label or facet.

Every scale here satisfies a lightness band, a chroma floor, all-pairs
separation under deuteranopia, protanopia and tritanopia, and 3:1 contrast
against white. `docs/colour_proof.pdf` shows each of them under normal vision,
the three dichromacies and greyscale.

Size is absolute
----------------
Figures are built at final printed width with font sizes in points, so a 7 pt
tick label prints at 7 pt. This holds only if the figure is placed at 100%:
scaling it afterwards scales its type with it. `figure()` sets a journal column
width, and `audit()` checks a finished figure for text below the 6 pt minimum, a
non-standard width, an over-height canvas, and legends overlapping data.

Text is plain
-------------
Species names print in italic and dN/dS prints upright, both in the same
typeface as every other label. Neither goes through mathtext: matplotlib
writes mathtext into the SVG glyph by glyph in its math font (DejaVu Sans on
the cluster), while plain text declares the Arial stack, so a mathtext label
opens in Illustrator in a second typeface. `house_text()` applies both
conventions to a finished figure, and `save()` and `audit()` call it, so a
script only has to pass `SPECIES_LABEL` and write "dN/dS".
"""

from __future__ import annotations

import datetime as _dt
import os as _os
import re as _re

import matplotlib as _mpl
import matplotlib.pyplot as _plt
import matplotlib.ticker as _mticker

__all__ = [
    "REGION", "REGION_ORDER", "REGION_SHADE", "SPECIES", "SPECIES_ORDER",
    "SPECIES_LABEL", "NEUTRAL", "PROXIMITY_HATCH", "pair_ramp", "apply_style",
    "make_outdir", "save", "COL", "COL_MM", "FONT", "LINEWIDTH", "MIN_FONT_PT",
    "mm", "pt", "figure", "panel_letter", "legend_outside", "headroom", "audit",
    "house_text", "ITALIC_NAMES", "NIGON", "NIGON_ORDER",
]

# ---------------------------------------------------------------------------
# 1. REGION -- the locked convention. Retained is drawn LEFT, eliminated RIGHT.
#    Validated as a pair: normal dE 22.8, protan 14.5, tritan 28.8, contrast OK.
# ---------------------------------------------------------------------------
REGION = {
    "retained":   "#4C72B0",
    "eliminated": "#C44E52",
}
REGION_ORDER = ["retained", "eliminated"]

# Background shading for region tracts in landscape/meta-profile panels.
# Pale enough to sit under data without competing with it.
REGION_SHADE = {
    "retained":   "#EAF0F8",
    "eliminated": "#FBECEC",
}

# ---------------------------------------------------------------------------
# 2. SPECIES -- categorical, four slots, fixed order, never cycled.
#    Use ONLY in panels that carry no region colour (Ks tree, Ks distributions,
#    ortholog-quality panels). Validated all-pairs: worst normal dE 15.4,
#    worst CVD dE 9.9 (protan), all >= 3:1 contrast on white.
#    Derived from ColorBrewer Dark2 + BrBG brown; deliberately avoids the blue
#    and red hue families so it cannot be confused with REGION across figures.
# ---------------------------------------------------------------------------
SPECIES = {
    "ascaris":       "#D95F02",  # orange
    "parascaris":    "#1B9E77",  # blue-green
    "baylisascaris": "#7570B3",  # slate purple
    "toxocara":      "#8C510A",  # brown
}
SPECIES_ORDER = ["ascaris", "parascaris", "baylisascaris", "toxocara"]

# Display names, as plain text. Genus-only for Ascaris per the manuscript
# convention. house_text(), which save() and audit() run, sets them in italic.
# These were mathtext until 2026-10-09: "$\it{P.\ univalens}$" printed in the
# math font rather than the label font, and mathtext's own space after the full
# stop plus "\ " gave a double-width gap (15.4 mm against 14.1 mm at 7 pt).
SPECIES_LABEL = {
    "ascaris":       "Ascaris",
    "parascaris":    "P. univalens",
    "baylisascaris": "B. schroederi",
    "toxocara":      "T. canis",
}

#: Every name house_text() sets in italic. Longest first, so that
#: "Parascaris univalens" is matched before "Parascaris".
ITALIC_NAMES = tuple(sorted(
    set(SPECIES_LABEL.values()) | {
        "Ascaris suum", "A. suum",
        "Parascaris univalens", "Parascaris",
        "Baylisascaris schroederi", "Baylisascaris",
        "Toxocara canis", "Toxocara",
    },
    key=lambda n: (-len(n), n)))

# ---------------------------------------------------------------------------
# 3. NEUTRAL -- for anything that is not region or species identity.
#    Metric traces (dN, dS, Ka/Ks) are all drawn in INK: the metric is named by
#    the axis label, and region is carried by the background shading. Drawing dN
#    red and dS blue collides with the region convention and must not be done.
# ---------------------------------------------------------------------------
NEUTRAL = {
    "ink":        "#2B2B2B",   # data lines, sliding-window medians
    "ink_light":  "#8A8A8A",   # secondary lines, per-gene points under a trace
    "rule":       "#6E6E6E",   # neutrality lines (omega = 1), reference rules
    "cbr":        "#E8B923",   # CBR marker lines -- yellow, as in the Fig 1A and 6F
                             # schematics (Fig 4 v6 hand fix, 2026-10-19; was
                             # green #1B7837)
    "sex_chrom":  "#B8860B",   # sex-chromosome accent in genome traces
    "grid":       "#DDDDDD",
}

# ---------------------------------------------------------------------------
# 4. SUB-CATEGORIES WITHIN A REGION (e.g. CBR-proximal vs interior, terminal vs
#    internal eliminated). Keep the region colour and vary the fill texture:
#    hue = region, texture = subclass, which also survives greyscale printing.
#    Tinting is not an option -- lightening blue drops its chroma below the point
#    at which it still reads as a hue.
# ---------------------------------------------------------------------------
PROXIMITY_HATCH = {
    "interior":     "",       # solid fill
    "cbr_proximal": "///",    # 45-degree hatch, same facecolor
}

# ---------------------------------------------------------------------------
# 5. PAIRWISE COMPARISONS -- six ortholog sets ordered by median Ks.
#    These are a MAGNITUDE, not an identity, so they take a single-hue
#    sequential ramp rather than six categorical hues (which cannot be made
#    CVD-safe). The three ascarid-ascarid pairs fall at the light end and the
#    three Toxocara pairs at the dark end automatically, which is the point.
# ---------------------------------------------------------------------------
_PAIR_RAMP = [
    "#CFC7E2", "#B3A6D4", "#9784C4", "#7B62B0", "#5F4396", "#432873",
]


def pair_ramp(n: int = 6) -> list[str]:
    """Return n sequential colours, light to dark, for magnitude-ordered series.

    Series MUST be sorted by the quantity being encoded (median Ks) before
    assignment, otherwise the ramp is meaningless.
    """
    if n < 2:
        raise ValueError("pair_ramp needs at least 2 steps")
    if n > len(_PAIR_RAMP):
        raise ValueError(
            f"pair_ramp is defined for up to {len(_PAIR_RAMP)} steps; "
            "more series than that should be faceted, not recoloured"
        )
    idx = [round(i * (len(_PAIR_RAMP) - 1) / (n - 1)) for i in range(n)]
    return [_PAIR_RAMP[i] for i in idx]


# ---------------------------------------------------------------------------
# 5b. NIGON ELEMENTS -- the seven ancestral linkage groups of Rhabditida.
#     A categorical identity, coloured as in the original Nigon-painting
#     figures. Element A is a dark red and C a light blue, so like SPECIES this
#     palette is SCOPED: use it only in panels that carry no region colour
#     (synteny / Nigon-painting figures). Marks drawn over it go in INK, since a
#     red mark reads as element A.
# ---------------------------------------------------------------------------
NIGON = {
    "A": "#932831",
    "B": "#F7941D",
    "C": "#50A9DD",
    "D": "#D7DF21",
    "E": "#52B747",
    "N": "#8780BC",
    "X": "#6D3285",
}
NIGON_ORDER = ["A", "B", "C", "D", "E", "N", "X"]


# ---------------------------------------------------------------------------
# 6. House style.
# ---------------------------------------------------------------------------
def apply_style(scale: float = 1.0) -> None:
    """Set the manuscript figure defaults, in absolute points.

    Sizes come from FONT and are absolute: a 7 pt tick label prints at 7 pt,
    provided the figure is built at final size and never scaled afterwards.
    `scale` exists only for making oversized reading copies -- leave it at 1.0
    for anything that goes to the journal.
    """
    f = {k: v * scale for k, v in FONT.items()}
    _params = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": f["tick"],
        "axes.titlesize": f["panel_title"],
        "axes.labelsize": f["axis_label"],
        "xtick.labelsize": f["tick"],
        "ytick.labelsize": f["tick"],
        "legend.fontsize": f["legend"],
        "figure.titlesize": f["suptitle"],

        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": LINEWIDTH["spine"],
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": LINEWIDTH["tick"],
        "ytick.major.width": LINEWIDTH["tick"],
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.minor.width": LINEWIDTH["hairline"],
        "ytick.minor.width": LINEWIDTH["hairline"],
        "axes.labelpad": 3.0,
        "axes.titlepad": 4.0,

        "axes.grid": False,
        "grid.color": NEUTRAL["grid"],
        "grid.linewidth": LINEWIDTH["hairline"],

        "lines.linewidth": LINEWIDTH["data"],
        "lines.markersize": 3.0,
        "patch.linewidth": LINEWIDTH["thin"],
        "boxplot.boxprops.linewidth": LINEWIDTH["thin"],
        "boxplot.whiskerprops.linewidth": LINEWIDTH["thin"],
        "boxplot.capprops.linewidth": LINEWIDTH["thin"],
        "boxplot.medianprops.linewidth": LINEWIDTH["data"],
        "boxplot.medianprops.color": "black",
        "boxplot.flierprops.markersize": 2.0,
        "boxplot.flierprops.markeredgewidth": LINEWIDTH["hairline"],

        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5,
        "legend.borderaxespad": 0.0,
        "legend.columnspacing": 1.0,

        "svg.fonttype": "none",
        # Element ids in an SVG are generated from a salt that is random per
        # process unless it is fixed. Fixing it makes two runs of the same
        # script byte-identical rather than differing in a hundred lines of
        # renamed <path id="..."> definitions that carry no information.
        "svg.hashsalt": "ascarid-molecular-evolution",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        # 600 ppi so any rasterized=True layer meets the line-art requirement
        # while text and axes stay vector.
        "savefig.dpi": 600,
        "figure.dpi": 600,
        "savefig.bbox": "standard",
        "savefig.pad_inches": 0.01,
        "figure.constrained_layout.use": True,
        "figure.constrained_layout.h_pad": 0.02,
        "figure.constrained_layout.w_pad": 0.02,
    }
    # Applied key by key rather than in one update(): rcParams.update() raises on
    # the first key the installed matplotlib does not recognise, which would take
    # the whole run down over a cosmetic setting. Cluster and laptop rarely have
    # the same version.
    unknown = []
    for key, value in _params.items():
        try:
            _mpl.rcParams[key] = value
        except (KeyError, ValueError):
            unknown.append(key)
    if unknown:
        import warnings
        warnings.warn("paper_colors.apply_style: matplotlib "
                      f"{_mpl.__version__} does not accept {unknown}; "
                      "figures will still render, check spacing.")


def region_colors(order=None) -> list[str]:
    """Region colours in canonical drawing order (retained first, i.e. left)."""
    return [REGION[k] for k in (order or REGION_ORDER)]


def species_colors(order=None) -> list[str]:
    """Species colours in canonical order."""
    return [SPECIES[k] for k in (order or SPECIES_ORDER)]


# ---------------------------------------------------------------------------
# 6b. TEXT CONVENTIONS -- species names italic, dN/dS upright, one typeface.
#     house_text() is the only public entry point; save() and audit() call it.
# ---------------------------------------------------------------------------

#: A species name as a whole word: not part of a longer word, and not directly
#: after a full stop (so "A. suum" is not matched inside some other token).
_NAME_RE = _re.compile(
    r"(?<![\w.])(?:" + "|".join(_re.escape(n) for n in ITALIC_NAMES) + r")(?!\w)")
_NAME_SET = frozenset(ITALIC_NAMES)

#: What may sit between names in a text that is still "only species names":
#: "Ascaris × P. univalens", "Ascaris vs T. canis", "Ascaris / Parascaris".
_JOINERS_RE = _re.compile(
    r"(?:\s+|[×/|,;&+()–—-]|\bx\b|\bvs\b\.?|\bversus\b|\band\b)*")

#: Mathtext spellings of dN/dS with italic N and S, rewritten to plain text.
#: Order matters: the paired form first, so "$dN/dS$" is not split in two.
_ITALIC_RATE = [
    # $dN/dS$, $d_N/d_S$, $\it{dN/dS}$                     -> dN/dS
    (_re.compile(r"\$\s*(?:\\(?:it|mathit)\s*)?\{?\s*d\s*_?\{?\s*([NS])\s*\}?"
                 r"\s*/\s*d\s*_?\{?\s*([NS])\s*\}?\s*\}?\s*\$"), r"d\1/d\2"),
    # d$N$, d$\it{N}$, d$\mathit{N}$, d$_N$, d$_{N}$         -> dN
    (_re.compile(r"(?<![A-Za-z\\])d\$\s*(?:\\(?:it|mathit)\s*)?\{?\s*_?\s*\{?"
                 r"\s*([NS])\s*\}?\s*\}?\s*\$"), r"d\1"),
    # $dN$, $d_N$, $d_{N}$                                  -> dN
    (_re.compile(r"\$\s*(?:\\(?:it|mathit)\s*)?\{?\s*d\s*_?\{?\s*([NS])\s*\}?"
                 r"\s*\}?\s*\$"), r"d\1"),
]
#: Anything still looking like a mathtext dN or dS after the rewrite.
_RATE_LEFTOVER = _re.compile(r"(?<![A-Za-z\\])d\$|\$[^$]*\bd_\{?[NS]")

#: A mathtext italic run, "$\it{...}$" or "$\mathit{...}$".
_MATH_ITALIC_RE = _re.compile(r"\$\s*\\(?:it|mathit|textit)\s*\{([^{}$]*)\}\s*\$")

_MIXED_NOTE = ("species name inside other text, {s!r}: set in mathtext italic so "
               "it stays italic, which prints it in the math font -- give the "
               "name its own text (axis label, title or tick label)")


def _split_math(s: str):
    """Split on unescaped '$' into [(is_math, text), ...]; None if the count is odd."""
    parts = _re.split(r"(?<!\\)\$", s)
    if len(parts) % 2 == 0:
        return None
    return [(i % 2 == 1, p) for i, p in enumerate(parts)]


def _unwrap_species(s: str) -> str:
    """Replace "$\\it{P.\\ univalens}$"-style runs around a known name with the name."""
    def repl(m):
        inner = m.group(1)
        for a, b in ((r"\ ", " "), (r"\,", " "), (r"\;", " "), (r"\:", " "),
                     ("~", " "), (r"\!", "")):
            inner = inner.replace(a, b)
        inner = " ".join(inner.split())
        return inner if inner in _NAME_SET else m.group(0)
    return _MATH_ITALIC_RE.sub(repl, s)


def _species_only(s: str) -> bool:
    """True if every line holds only species names and joiners, with >= 1 name."""
    if "$" in s:
        return False
    found = False
    for line in s.split("\n"):
        rest, n = _NAME_RE.subn(" ", line)
        found = found or n > 0
        if not _JOINERS_RE.fullmatch(rest):
            return False
    return found


def _tex_name(name: str) -> str:
    """A name for mathtext. Mathtext already puts a thin space after '.', so
    '\\,' (not '\\ ') makes the gap after an abbreviated genus one word space."""
    return name.replace(". ", r".\,").replace(" ", r"\ ")


def _house_string(s: str):
    """Apply the conventions to one string; returns (text, set_italic, notes)."""
    notes = []
    new = s
    for pat, rep in _ITALIC_RATE:
        new = pat.sub(rep, new)
    if _RATE_LEFTOVER.search(new):
        notes.append(f"italic dN/dS left in {s!r} -- write it as plain 'dN/dS'")

    unwrapped = _unwrap_species(new)
    if _species_only(unwrapped):
        return unwrapped, True, notes

    parts = _split_math(new)
    if parts is None:
        if _NAME_RE.search(new):
            notes.append(f"species name in {s!r} left upright: the text has an "
                         "odd number of '$'")
        return new, False, notes
    if any(_NAME_RE.search(seg) for is_math, seg in parts if not is_math):
        new = "".join(
            ("$" + seg + "$") if is_math else _NAME_RE.sub(
                lambda m: r"$\mathit{" + _tex_name(m.group(0)) + "}$", seg)
            for is_math, seg in parts)
        notes.append(_MIXED_NOTE.format(s=s))
    elif unwrapped != new:
        notes.append(_MIXED_NOTE.format(s=s))
    return new, False, notes


class _HouseFormatter(_mticker.Formatter):
    """Wrap a tick formatter so the labels it emits follow the house conventions.

    set_xticklabels() stores its strings inside the formatter (a FixedFormatter,
    or a FuncFormatter over a dict when the tick positions were fixed first),
    which rewrites every tick label at each draw. Wrapping survives that and
    does not depend on how a given matplotlib version stores the strings.
    """

    def __init__(self, inner):
        self.inner = inner

    def set_axis(self, axis):
        super().set_axis(axis)
        self.inner.set_axis(axis)

    def set_locs(self, locs):
        self.inner.set_locs(locs)

    def __call__(self, x, pos=None):
        return _house_string(str(self.inner(x, pos)))[0]

    def format_ticks(self, values):
        return [_house_string(str(v))[0] for v in self.inner.format_ticks(values)]

    def get_offset(self):
        return self.inner.get_offset()

    def format_data(self, value):
        return self.inner.format_data(value)

    def format_data_short(self, value):
        return self.inner.format_data_short(value)


def house_text(fig) -> list[str]:
    """Apply the house text conventions to every text in a figure.

    1. dN/dS is upright. Mathtext spellings with italic N and S -- d$N$/d$S$,
       d$\\it{N}$, $d_N$, $dN/dS$ -- are rewritten as plain dN/dS, dN and dS.
    2. Species names are italic in the label typeface. A text made only of
       names in ITALIC_NAMES, optionally joined by "×", "vs", "/", "," and the
       like, gets fontstyle='italic'; mathtext names in it are unwrapped to
       plain text first.
    3. A species name inside other text ("Ascaris\\n(n = 5,733)") cannot be
       partly italic in plain text. It is wrapped in mathtext italic, with the
       gap after an abbreviated genus corrected, and reported: it prints in
       the math font, and the fix is to give the name a text of its own.

    Tick labels are regenerated by their formatter at every draw, so editing
    the drawn text alone would be undone by the next draw. A formatter whose
    labels need changing is therefore wrapped so that it emits the converted
    labels itself.
    Idempotent. Returns notes for audit(); empty when nothing needs attention.
    """
    import matplotlib.text as _mtext

    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for get, put, locs in (
                    (axis.get_major_formatter, axis.set_major_formatter,
                     axis.get_majorticklocs),
                    (axis.get_minor_formatter, axis.set_minor_formatter,
                     axis.get_minorticklocs)):
                fmt = get()
                if isinstance(fmt, _HouseFormatter):
                    continue
                try:
                    labels = [str(v) for v in fmt.format_ticks(list(locs()))]
                except Exception:          # a formatter that cannot run detached
                    continue
                if any(_house_string(v)[0] != v for v in labels):
                    put(_HouseFormatter(fmt))

    notes = []
    for t in fig.findobj(_mtext.Text):
        s = t.get_text()
        if not s:
            continue
        new, italic, msgs = _house_string(s)
        if new != s:
            t.set_text(new)
        if italic and t.get_fontstyle() != "italic":
            t.set_fontstyle("italic")
        notes.extend(msgs)
    return list(dict.fromkeys(notes))


# ---------------------------------------------------------------------------
# 7. Run directories.
#    Every script writes into <base>/<analysis>_<YYYYMMDD>[_<tag>]/ so a rerun
#    cannot overwrite the run a figure was built from, and any stored figure can
#    be matched to the run that produced it.
# ---------------------------------------------------------------------------
def make_outdir(analysis: str, base: str = ".", tag: str | None = None,
                date: str | None = None, exist_ok: bool = True) -> str:
    """Create and return a dated output directory.

    make_outdir("dnds_decomposition", base="results", tag="as_region")
      -> results/dnds_decomposition_20260814_as_region/

    Pass date="20260812" to write into an existing run rather than today's.
    """
    stamp = date or _dt.date.today().strftime("%Y%m%d")
    name = f"{analysis}_{stamp}" + (f"_{tag}" if tag else "")
    path = _os.path.join(base, name)
    _os.makedirs(path, exist_ok=exist_ok)
    return path


#: Metadata keys carrying a wall-clock timestamp, per format. Set to None so two
#: runs of the same script on the same input produce identical files.
_NO_TIMESTAMP = {"svg": {"Date": None},
                 "pdf": {"CreationDate": None}}


def save(fig, outdir: str, stem: str, formats=("svg", "pdf", "png")) -> list[str]:
    """Save one figure in every export format and return the paths written.

    Timestamps are suppressed so a rerun is byte-identical to the run before it.
    Without this, two runs of the same script on the same input differ in the
    SVG's ``dc:date`` and the PDF's ``CreationDate``, which makes ``diff`` on a
    pair of runs useless for the one question worth asking of it -- whether the
    figure actually changed. ``svg.hashsalt``, set in ``apply_style()``, does the
    same job for the element identifiers matplotlib generates internally.
    """
    # FREEZE THE GEOMETRY BEFORE WRITING ANYTHING.
    #
    # Writing the PDF re-solves the layout and overwrites hand-set axes
    # positions with a constrained-layout solution, because
    # ``figure.constrained_layout.use`` is True in rcParams. The formats are
    # written svg -> pdf -> png, so the SVG came out with the geometry the
    # script asked for and the PDF and PNG came out with something else: a
    # panel set to bottom=0.20 landed at 0.054, and the PNG proof therefore did
    # not show the figure that the SVG placed in Illustrator. Verified
    # 2026-08-18 on a 2x4 grid -- ``get_position()`` reported 0.20 before the
    # PDF write and 0.054 after it.
    #
    # Draw once so any layout engine does its work, pin every axes to the
    # result, then write. Figures that rely on constrained_layout are
    # unaffected: they get their solved positions pinned, which is exactly what
    # audit() measured.
    #
    # house_text() first, so the text that is measured and pinned is the text
    # that is written: species in italic, dN/dS upright. A no-op when audit()
    # has already run it.
    house_text(fig)
    fig.canvas.draw()
    frozen = [(ax, ax.get_position().frozen()) for ax in fig.axes]
    fig.set_layout_engine("none")

    written = []
    for ext in formats:
        for ax, pos in frozen:
            ax.set_position(pos)
        p = _os.path.join(outdir, f"{stem}.{ext}")
        meta = _NO_TIMESTAMP.get(ext.lstrip("."))
        fig.savefig(p, **({"metadata": meta} if meta else {}))
        written.append(p)
    return written


# ---------------------------------------------------------------------------
# 8. PAGE GEOMETRY -- Molecular Biology and Evolution / Oxford University Press
#
#     Type area           234 mm high
#     Double column       185 mm wide
#     Single column        89 mm wide
#     Minimum label size     6 pt at final magnification
#     Fonts               sans serif, Arial or Helvetica, uniform within a figure
#     Resolution          600 ppi line art, 300 ppi colour/halftone
#                         (MBE asks 1200 ppi for pure line drawings in TIFF;
#                          submitting vector PDF/EPS avoids the question entirely)
#     Colour              CMYK for print; colour is free online
#
# Build every figure at its final printed size, in inches, with font sizes in
# points. Points are absolute, so a 7 pt tick label is 7 pt on paper -- provided
# the figure is placed at 100%. Scaled to 80%, that same label prints at 5.6 pt,
# below the journal minimum.
#
# Where one figure's panels come from several scripts, each script emits a
# full-width horizontal strip at COL["double"] and its own height. The strips
# stack without any scaling.
# ---------------------------------------------------------------------------
MM_PER_IN = 25.4
PT_PER_IN = 72.0


def mm(x: float) -> float:
    """Millimetres to inches, for figsize."""
    return x / MM_PER_IN


def pt(x: float) -> float:
    """Points to inches."""
    return x / PT_PER_IN


COL_MM = {"single": 89.0, "onehalf": 136.0, "double": 185.0}
COL = {k: mm(v) for k, v in COL_MM.items()}

MAX_HEIGHT_MM = 234.0
MAX_HEIGHT = mm(MAX_HEIGHT_MM)

#: Leave room for the caption on a figure that fills the page.
MAX_HEIGHT_WITH_CAPTION_MM = 200.0

#: Journal floor. Nothing in a figure may be smaller than this at final size.
MIN_FONT_PT = 6.0

#: Smallest plausible drawn panel. Below this an axes has not been designed
#: small, it has been squeezed by a layout solver making room for something that
#: does not fit. audit() reports it.
MIN_AXES_MM = 15.0

#: Fraction of the figure (width, height) the axes must collectively span.
#: Generous enough for a legend below the panels and an axis label beside them.
MIN_AXES_SPAN = (0.55, 0.45)

#: House sizes, in points at final printed size. Chosen so a 10% emergency
#: reduction still clears the 6 pt floor.
FONT = {
    "tick":        7.0,
    "axis_label":  8.0,
    "panel_title": 8.0,
    "annotation":  7.0,   # p-values, n =, fold changes drawn on the axes
    "legend":      7.0,
    "panel_letter": 10.0,  # A, B, C -- bold
    "suptitle":    9.0,
}

#: Stroke weights in points. 0.25 pt is the practical dropout floor in print.
LINEWIDTH = {
    "spine":   0.75,
    "tick":    0.75,
    "data":    1.2,
    "thin":    0.6,
    "rule":    0.9,   # dashed reference lines
    "hairline": 0.4,
}


def figure(width: str = "double", height_mm: float = 60.0, **kw):
    """Create a figure at exact journal width. width in {single, onehalf, double}."""
    import matplotlib.pyplot as plt
    if width not in COL:
        raise ValueError(f"width must be one of {sorted(COL)}")
    if height_mm > MAX_HEIGHT_MM:
        raise ValueError(f"height {height_mm} mm exceeds the {MAX_HEIGHT_MM} mm type area")
    return plt.subplots(figsize=(COL[width], mm(height_mm)), **kw)


def panel_letter(ax, letter: str, dx: float = -0.10, dy: float = 1.06):
    """Place a bold panel letter in axes coordinates, outside the axes box."""
    return ax.text(dx, dy, letter, transform=ax.transAxes,
                   fontsize=FONT["panel_letter"], fontweight="bold",
                   va="bottom", ha="left")


def legend_outside(ax, where: str = "right", **kw):
    """Put the legend outside the axes so it can never cover data.

    where='right'  -> to the right of the axes
    where='top'    -> above the axes, horizontal
    """
    opts = dict(frameon=False, fontsize=FONT["legend"], handlelength=1.4,
                handletextpad=0.5, borderaxespad=0.0)
    opts.update(kw)
    if where == "right":
        return ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), **opts)
    if where == "top":
        return ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02),
                         ncol=opts.pop("ncol", 3), **opts)
    raise ValueError("where must be 'right' or 'top'")


def headroom(ax, frac: float = 0.18):
    """Extend the y-axis upward so an in-axes legend sits over blank space.

    Use when the legend must be inside the axes: reserve the space first, then
    place the legend in it. Never drop a legend onto data and hope.
    """
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, lo + (hi - lo) * (1.0 + frac))


def audit(fig, verbose: bool = True):
    """Check a figure against the journal rules before saving.

    Reports: any text below MIN_FONT_PT, a width that is not a legal column
    width, a height over the type area, any legend whose box overlaps a data
    artist, and any species name or dN/dS that house_text() could not bring to
    the house convention. house_text() runs first, so every later check sees
    the text that save() will write. Returns a list of problem strings (empty
    means clean).
    """
    import matplotlib.text as _mtext
    import numpy as np
    problems = list(house_text(fig))

    w_in, h_in = fig.get_size_inches()
    w_mm, h_mm = w_in * MM_PER_IN, h_in * MM_PER_IN
    if not any(abs(w_mm - v) < 0.5 for v in COL_MM.values()):
        problems.append(
            f"width {w_mm:.1f} mm is not a column width "
            f"({', '.join(f'{k} {v:g}' for k, v in COL_MM.items())})")
    if h_mm > MAX_HEIGHT_MM:
        problems.append(f"height {h_mm:.1f} mm exceeds the {MAX_HEIGHT_MM:g} mm type area")

    small = {}
    for t in fig.findobj(_mtext.Text):
        s = (t.get_text() or "").strip()
        if not s:
            continue
        size = t.get_fontsize()
        if size < MIN_FONT_PT - 1e-6:
            small.setdefault(round(size, 2), []).append(s[:24])
    for size, items in sorted(small.items()):
        problems.append(f"{len(items)} text object(s) at {size} pt "
                        f"(floor is {MIN_FONT_PT:g} pt): {items[:3]}")

    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is None:
            continue
        lb = leg.get_window_extent(rend)
        hits = 0
        for art in list(ax.lines) + list(ax.collections) + list(ax.patches):
            # Test the drawn geometry, not the bounding box. A step histogram or
            # a long line is one artist whose box spans the whole panel, so a box
            # test reports a collision for any in-axes legend even when the ink
            # is nowhere near it. Sampling vertices asks the question that
            # matters: does anything actually drawn fall inside the legend?
            disp = None
            if hasattr(art, 'get_xydata'):                     # Line2D
                pts = np.asarray(art.get_xydata(), dtype=float)
                pts = pts[np.isfinite(pts).all(axis=1)]
                if len(pts):
                    disp = ax.transData.transform(pts)
            elif hasattr(art, 'get_offsets'):                  # scatter, collections
                pts = np.asarray(art.get_offsets(), dtype=float)
                pts = pts[np.isfinite(pts).all(axis=1)]
                if len(pts):
                    disp = ax.transData.transform(pts)
            elif hasattr(art, 'get_path'):                     # Polygon, Rectangle
                # hist(histtype='step') draws each series as one Polygon, whose
                # bounding box spans the panel; its vertices are what is drawn.
                #
                # Patch.get_transform() ALREADY composes get_patch_transform()
                # with the data transform. Applying get_patch_transform() first
                # and then get_transform() applies it twice, which is the
                # identity for a Polygon (patch transform = identity) but not
                # for a Rectangle from ax.bar(), whose patch transform maps the
                # unit square onto the bar. Every bar panel with a legend
                # therefore got garbage display coordinates -- reported as a
                # spurious overlap on 5D, and equally capable of hiding a real
                # one. Fixed 2026-08-18.
                try:
                    v = np.asarray(
                        art.get_path().transformed(art.get_transform()).vertices,
                        dtype=float)
                    disp = v[np.isfinite(v).all(axis=1)]
                except Exception:
                    disp = None
            if disp is not None and len(disp):
                # Densify: a sparse polyline can cross the legend box between two
                # vertices without either vertex being inside it.
                if len(disp) > 1 and len(disp) < 2000:
                    k = int(np.ceil(2000 / len(disp)))
                    t = np.linspace(0, 1, k, endpoint=False)[:, None, None]
                    segs = disp[:-1][None] + t * (disp[1:] - disp[:-1])[None]
                    disp = np.vstack([segs.reshape(-1, 2), disp[-1:]])
                inside = ((disp[:, 0] >= lb.x0) & (disp[:, 0] <= lb.x1) &
                          (disp[:, 1] >= lb.y0) & (disp[:, 1] <= lb.y1))
                if inside.any():
                    hits += 1
                continue
            try:
                ab = art.get_window_extent(rend)
            except Exception:
                continue
            if ab.width > 0 and ab.height > 0 and lb.overlaps(ab):
                hits += 1
        if hits:
            problems.append(f"legend on {ax.get_label() or 'an axes'} overlaps "
                            f"{hits} data artist(s) -- move it outside or call headroom()")

    # Collapsed axes. constrained_layout treats the panel as the flexible part
    # of the figure: asked for a legend or a label wider than the canvas, it
    # shrinks the axes to make room and emits at most a UserWarning, so a figure
    # can save successfully with a panel a few millimetres wide. Nothing else in
    # this audit catches that -- the width is legal, the type is legal, and the
    # legend overlaps nothing because there is nothing left to overlap.
    boxes = []
    for ax in fig.axes:
        if ax.get_label() == "<colorbar>" or getattr(ax, "_colorbar", None) is not None:
            continue
        pos = ax.get_position()
        boxes.append(pos)
        aw, ah = pos.width * w_mm, pos.height * h_mm
        if aw < MIN_AXES_MM or ah < MIN_AXES_MM:
            problems.append(
                f"axes drawn at {aw:.1f} x {ah:.1f} mm (floor is {MIN_AXES_MM:g} mm) -- "
                f"layout collapsed, usually a legend or label wider than the figure")

    # A milder squeeze passes the floor above and is still wrong: one panel
    # given 37 mm of an 89 mm figure because a legend claimed the rest is not a
    # small panel, it is a failed layout. Measuring the span of all the axes
    # together rather than each one separately keeps a legitimate three-panel
    # strip clean, since together those panels still cover the page.
    if boxes:
        x_span = max(b.x1 for b in boxes) - min(b.x0 for b in boxes)
        y_span = max(b.y1 for b in boxes) - min(b.y0 for b in boxes)
        if x_span < MIN_AXES_SPAN[0]:
            problems.append(
                f"axes span only {x_span * 100:.0f}% of the figure width "
                f"(expected >= {MIN_AXES_SPAN[0] * 100:.0f}%) -- something outside the "
                f"panels, usually a legend, is taking the page")
        if y_span < MIN_AXES_SPAN[1]:
            problems.append(
                f"axes span only {y_span * 100:.0f}% of the figure height "
                f"(expected >= {MIN_AXES_SPAN[1] * 100:.0f}%)")

    if verbose:
        if problems:
            print("FIGURE AUDIT -- " + str(len(problems)) + " problem(s):")
            for p_ in problems:
                print("  !", p_)
        else:
            print("FIGURE AUDIT -- clean")
    return problems
