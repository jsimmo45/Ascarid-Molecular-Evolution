#!/usr/bin/env python3
"""
figure_cohort_selection.py — supplementary fig. S14: how the cohort was chosen
and why its genotypes are comparable between region classes. (S16 before the
2026-10-12 renumbering; the output stem keeps S16_cohort_selection.)

What this assembles, and why it is one figure
---------------------------------------------
Four steps stand between "public Ascaris reads" and "a callset that can be
compared between retained and eliminated regions". Each is a separate script and
each, alone, looks like housekeeping:

    A  select_germline_samples.py   which libraries contain enough gonad?
    B  select_germline_samples.py   which sex is each worm?
    C  depth_match.py               is coverage equal across region classes?
    D  allele_balance_panel.py      what does the paralogy screen remove?

Together they are one argument, and the argument has a direction. Eliminated
regions are present only in germline DNA, so a somatic-heavy library genotypes
them at lower depth than retained regions *by construction* (A). Sex matters
because the X is retained and hemizygous in males, so an unbalanced sex
composition moves genome-wide heterozygosity (B). Lower depth in one class
lowers the heterozygote call rate in that class, which is a direct bias on
pi_N/pi_S — so depth is equalised per worm and the achieved depth verified, not
assumed (C). And what survives all of that can still be inflated by collapsed
paralogues, which the allele-balance screen removes several-fold more of in
eliminated regions than in retained ones (D).

A and B are selection. C and D are the two corrections that make the surviving
comparison legitimate. D is the load-bearing one: it is the control that would
most plausibly have produced the headline result as an artifact, and it does
not.

This is an ASSEMBLER. It reads summary tables and draws; it recomputes nothing.
The source scripts keep their own working figures. Nothing here can disagree
with them, because nothing here recalculates.

Finding the inputs
------------------
This project is split across two trees — code and cohort tables under
``$POPGEN_WORK_DIR``, results under ``$POPGEN_DIR``
— and make_outdir() puts each result in a DATED subdirectory, so no input is at
a path worth memorising. All three arguments therefore accept a directory to
search as well as a file:

    --metadata        a file, or a directory searched for sample_metadata.tsv
    --allele-balance  a file, or a directory searched for
                      allele_balance_screen.tsv
    --depth-match-glob  a quoted glob, or a directory searched recursively for
                      *.depth_match.tsv

Panel C needs one ``<sample>.depth_match.tsv`` per worm — depth_match.py writes
one file per BAM — which is why that argument is a glob rather than a path.

If you would rather look first:

    find $POPGEN_DIR \\
         $POPGEN_WORK_DIR \\
         \\( -name 'sample_metadata.tsv' -o -name 'allele_balance_screen.tsv' \\
            -o -name '*.depth_match.tsv' \\) -printf '%TF %10s  %p\\n' | sort

Example
-------
    python figure_cohort_selection.py \\
        --metadata         cohort_20260820 \\
        --depth-match-glob $POPGEN_DIR/bam_matched \\
        --allele-balance   $POPGEN_DIR/results/figures \\
        --outdir           $POPGEN_DIR/results/figures

Quote a glob if you pass one. An unquoted pattern is expanded by the shell into
many positional arguments and argparse takes only the first.

Changes 2026-10-12
------------------
Panel D's log axis ran to 10^4 % in the 2026-08-22 build. Its limits were taken
from every number stored for the panel, which includes ``no_ad_evaluable`` -- a
raw COUNT of sites, in the thousands -- so the top was six times that count.
The limits now come from the six drawn percentages only, the top stops at
100 %, and the log ticks print as plain numbers (0.1, 1, 10, 100). The same
limits and ticks are applied in allele_balance_panel.py, so the standalone
S16D and this panel still match.

Changes 2026-10-13
------------------
Panel D's bars were based at 0, which has no position on a log axis:
matplotlib maps it to 10^-1000, so each bar was written into the SVG and PDF
with its bottom edge ~54,000 pt (750 in) below the page, hidden only by the
axes clip path. Illustrator's canvas is 227 in, so it refused to move, scale
or paste any selection that held one of those bars ("the requested
transformation would make some objects too large"). The bars are now based on
the lower y-limit, so they end at the x-axis and stay on the page; the
rendered figure is pixel-identical. A percentage that is zero or missing has
no position on the axis either and would go off the page the same way, so it
is now left undrawn, logged, and listed under ``panel_d.not_drawn`` in
provenance.json. allele_balance_panel.py was not changed.

Changes 2026-10-31
------------------
Panel B's key calls the worm in the shaded band "Unassigned" (was "Unknown"),
the word the legend uses; the band itself keeps "Ambiguous". Nothing else
changes.
"""

from __future__ import annotations

import argparse
import glob as globmod
import json
import logging
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

from style_boot import RETAINED, ELIMINATED, paper_colors

if paper_colors is None:                # pragma: no cover
    raise SystemExit(
        "figure_cohort_selection.py builds a journal supplement and needs the "
        "full house style, not just the two region colours. Set "
        "PAPER_COLORS_DIR to the directory holding paper_colors.py (on Isaac: "
        ".../miniprot/Ascarid-Molecular-Evolution/lib).")

apply_style = paper_colors.apply_style
audit = paper_colors.audit
figure = paper_colors.figure
headroom = paper_colors.headroom
make_outdir = paper_colors.make_outdir
panel_letter = paper_colors.panel_letter
pc_save = paper_colors.save
FONT = paper_colors.FONT
LINEWIDTH = paper_colors.LINEWIDTH
INK = paper_colors.NEUTRAL["ink"]
INK_LIGHT = paper_colors.NEUTRAL["ink_light"]
RULE = paper_colors.NEUTRAL["rule"]

FIG_W = "double"
FIG_H_MM = 136.0

#: Panel D groups, kept identical to allele_balance_panel.py so the standalone
#: panel and the assembled one cannot drift apart.
AB_GROUPS = [("ab", "Allele balance\noutside window"),
             ("allhet", "All genotypes\nheterozygous"),
             ("dropped", "Either\n(total removed)")]
REGIONS = ["retained", "eliminated"]

#: Above this many bars, panel A drops per-library tick labels for a rank axis.
#: Fifty accessions at 7 pt across 89 mm overlap into a solid band; the audit
#: routinely covers that many libraries, so the switch has to be automatic
#: rather than a flag someone remembers to set.
MAX_TICK_LABELS = 20

#: Controls are not worms and must never be drawn as though they were; they are
#: the assay floor and ceiling that germline_rescaled is defined against.
CONTROL_HATCH = "///"

#: Panel D plots shares of variant sites, and a share cannot exceed 100 %, so
#: the log axis stops there. Same value in allele_balance_panel.py.
PCT_MAX = 100.0


def plain_log_ticks(axis) -> None:
    """Label a log axis 0.1, 1, 10, 100 rather than 10^-1 ... 10^2.

    The axis is a percentage and should read as one. The default exponent
    labels are also mathtext, which the SVG writes glyph by glyph in the math
    font rather than the label typeface. Minor ticks stay, unlabelled.
    Identical in allele_balance_panel.py.
    """
    axis.set_major_locator(LogLocator(base=10.0))
    axis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:g}"))
    axis.set_minor_formatter(NullFormatter())


def strip_suffix(name: str) -> str:
    """``SRR31692074.sub3000000`` -> ``SRR31692074``.

    The pilot audit ran on subsampled BAMs and the full pipeline on the whole
    ones, so the same worm appears under two stems depending on which stage
    wrote the table. Normalising here is what lets the cohort cross-check below
    be strict.
    """
    return re.sub(r"\.sub\d+$", "", str(name))


def short(name: str, n: int = 4) -> str:
    """Final n characters of a run accession.

    Full SRR accessions on the x-axis of an 88 mm panel leave no room for the
    panel. The last four digits are unique within this cohort; the caption
    states the convention.
    """
    return str(name)[-n:]


def load_depth_match(pattern: str) -> pd.DataFrame:
    """Concatenate the per-worm depth_match tables named by a glob.

    Accepts either a glob or a directory (in which case ``*.depth_match.tsv``
    inside it is assumed). Returns an empty frame with the right schema if
    nothing matches; the caller decides whether that is fatal.
    """
    pattern = os.path.expanduser(pattern)
    if Path(pattern).is_dir():
        # Recursive: depth_match.py is run per BAM and its --outdir is often a
        # per-sample subdirectory, so the files are one level down as often as
        # they are flat.
        paths = sorted(str(p) for p in Path(pattern).rglob("*.depth_match.tsv"))
    else:
        paths = sorted(globmod.glob(pattern))
    cols = ["sample", "depth_retained_in", "depth_eliminated_in",
            "frac_retained", "frac_eliminated", "depth_retained_out",
            "depth_eliminated_out", "target"]
    if not paths:
        return pd.DataFrame(columns=cols)
    frames = []
    for p in paths:
        try:
            frames.append(pd.read_csv(p, sep="\t"))
        except Exception as exc:
            logging.error("Could not read %s (%s) — skipped", p, exc)
    if not frames:
        return pd.DataFrame(columns=cols)
    df = pd.concat(frames, ignore_index=True)
    df["sample"] = df["sample"].map(strip_suffix)
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(
            f"depth_match tables are missing column(s) {missing}. These files "
            f"were not written by depth_match.py, or by a different version of "
            f"it — check {paths[0]}.")
    # A worm matched twice (a re-run left both files) would be double-plotted.
    dup = df["sample"].duplicated(keep="last")
    if dup.any():
        logging.warning("Duplicate depth_match rows for %s — keeping the last "
                        "file read for each. Delete stale files if that is not "
                        "the intended run.", sorted(df.loc[dup, "sample"]))
        df = df[~df["sample"].duplicated(keep="last")]
    return df.reset_index(drop=True)


def resolve_input(path: Path, filename: str) -> Path:
    """Return ``path`` if it is a file, else the newest ``filename`` beneath it.

    The inputs to this figure are written into dated output directories by
    make_outdir(), so their exact paths change every run and none of them is
    where a reader would guess. Accepting the PARENT and searching keeps the
    wrapper's defaults stable across re-runs. The newest match wins, and which
    file was taken is logged — a silent pick between two candidates is exactly
    the mismatched-run failure the cohort cross-check exists to catch.
    """
    path = path.expanduser()
    if path.is_file():
        return path
    if not path.is_dir():
        raise SystemExit(
            f"{path} is neither a file nor a directory. Pass {filename} "
            f"itself, or a directory to search for it.")
    hits = sorted(path.rglob(filename), key=lambda p: p.stat().st_mtime,
                  reverse=True)
    if not hits:
        raise SystemExit(
            f"No {filename} anywhere under {path}.\n"
            f"Locate it with:  find $POPGEN_DIR "
            f"$POPGEN_WORK_DIR "
            f"-name '{filename}'\n"
            f"If it does not exist yet, the script that writes it has not been "
            f"run against this callset.")
    if len(hits) > 1:
        logging.warning("%d copies of %s under %s — taking the newest:",
                        len(hits), filename, path)
        for h in hits[:5]:
            logging.warning("    %s", h)
    logging.info("Using %s", hits[0])
    return hits[0]


def panel_a(ax, meta: pd.DataFrame, min_germline: float) -> dict:
    """Germline fraction per library, with the selection threshold.

    Fill encodes selection outcome, not identity: pass is solid ink, fail is
    open, control is hatched. Region blue and red are reserved for retained
    versus eliminated everywhere in this paper and must not be spent on a
    distinction that is not that one.

    The x-axis is a RANK, not a list of names, whenever the audit covers more
    libraries than an 89 mm panel can label. Nearly fifty accessions printed at
    7 pt collide into a solid black band that carries no information at all, and
    the identity of a library that was screened out is not something a reader
    needs: what this panel has to show is the shape of the distribution and
    where the threshold falls in it. The ten libraries whose identity does
    matter are named in panels B and C.
    """
    d = meta
    x = np.arange(len(d), dtype=float)
    label_ticks = len(d) <= MAX_TICK_LABELS
    faces = np.where(d["pass"].to_numpy(), INK, "white")
    for i, r in enumerate(d.itertuples()):
        ax.bar(x[i], r.germline_fraction, 0.72,
               color=(INK if r.is_control else faces[i]),
               alpha=(0.35 if r.is_control else 1.0),
               hatch=(CONTROL_HATCH if r.is_control else ""),
               edgecolor=INK, linewidth=LINEWIDTH["thin"], zorder=2)
    ax.axhline(min_germline, color=RULE, ls="--", lw=LINEWIDTH["rule"],
               zorder=3, label=f"Threshold ({min_germline:.2f})")
    ax.set_ylabel("Germline fraction")
    ax.set_ylim(0, float(d.germline_fraction.max()))
    headroom(ax, 0.12)
    if label_ticks:
        ax.set_xticks(x)
        ax.set_xticklabels([short(s) for s in d["run"]], rotation=90)
        ax.set_xlabel("Library")
    else:
        # Rank ticks, spaced so the axis reads as an ordered axis rather than a
        # list. The controls sit at the right-hand end and are identified by
        # their hatch, not by a tick. Only the worms are ranked -- the controls
        # are appended after them, out of rank -- so the label says so.
        step = 10 if len(d) <= 60 else 20
        ticks = np.arange(0, len(d), step, dtype=float)
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(int(t) + 1) for t in ticks])
        n_ctrl = int(d["is_control"].sum())
        n_lib = int(len(d) - n_ctrl)
        ax.set_xlabel(f"{n_lib} libraries ranked, then {n_ctrl} tissue "
                      f"control{'s' if n_ctrl != 1 else ''}"
                      if n_ctrl else
                      f"Library, ranked by germline fraction (n = {len(d)})")
    ax.set_xlim(-0.7, len(d) - 0.3)
    ax.tick_params(axis="x", length=1.5, pad=1.5)

    worms = d[~d.is_control]
    ctrl = d[d.is_control]
    return {"n_libraries": int(len(d)),
            "n_worms": int(len(worms)),
            "n_pass": int(worms["pass"].sum()),
            "fail_reasons": worms.loc[~worms["pass"], "fail_reason"]
                                 .value_counts().to_dict(),
            "control_range": ([float(ctrl.germline_fraction.min()),
                               float(ctrl.germline_fraction.max())]
                              if len(ctrl) else None)}


def panel_b(ax, meta: pd.DataFrame, male_max: float, female_min: float) -> dict:
    """X:A per library, and the sex call it produces.

    The two thresholds are drawn as the boundaries of a shaded ambiguous band
    rather than as two lines, because what the panel has to show is that nothing
    in the analysis cohort lands inside it. Marker shape carries the call; shape
    survives greyscale and keeps the region hues free.

    ONLY THE SELECTED WORMS ARE DRAWN, and that is a claim about what the panel
    is for. Drawn over all ~50 audited libraries the panel is dominated by
    somatic-heavy preps whose X:A is estimated from almost no germline reads;
    those scatter into the ambiguous band and would make the sex call look
    unreliable when what is actually unreliable is the library. Every one of
    them has already been removed by panel A. The count of ambiguous calls among
    the excluded libraries is logged rather than drawn, so the number is on the
    record without being asserted as a property of the cohort.

    Panel A is the selection step; this panel and C describe what came out of
    it, in the same order, so a worm can be tracked between them.
    """
    sel = meta[(~meta.is_control) & meta["pass"]]
    amb_excluded = meta[(~meta.is_control) & (~meta["pass"])
                        & (meta.sex == "unknown")]
    if len(amb_excluded):
        logging.info("B  %d excluded librar(ies) had an ambiguous X:A; all were "
                     "already removed by the germline-fraction filter, and none "
                     "is drawn", len(amb_excluded))
    drop = sel[~np.isfinite(pd.to_numeric(sel.x_to_autosome, errors="coerce"))]
    if len(drop):
        logging.warning("B  %d SELECTED worm(s) have no X:A and cannot be "
                        "drawn: %s", len(drop),
                        ", ".join(short(s) for s in drop["run"]))
    d = sel[np.isfinite(pd.to_numeric(sel.x_to_autosome,
                                      errors="coerce"))].reset_index(drop=True)
    if d.empty:
        ax.text(0.5, 0.5, "no X:A values", ha="center", va="center",
                transform=ax.transAxes, fontsize=FONT["annotation"], color=INK)
        ax.set_xticks([])
        ax.set_yticks([])
        return {"calls_selected": {}, "calls_all_libraries": {},
                "n_unknown_selected": 0, "n_unknown_excluded": 0,
                "n_drawn": 0, "not_drawn": []}
    x = np.arange(len(d), dtype=float)
    ax.axhspan(male_max, female_min, color=RULE, alpha=0.12, lw=0, zorder=1,
               label="Ambiguous")
    markers = {"M": "o", "F": "s", "unknown": "X"}
    # "Unassigned", as the legend says: the X:A ratio puts the worm between
    # the thresholds, so no sex is assigned (the data value stays "unknown").
    labels = {"M": "Male", "F": "Female", "unknown": "Unassigned"}
    for call, mk in markers.items():
        sub = d[d.sex == call]
        if sub.empty:
            continue
        ax.scatter(x[d.sex.to_numpy() == call], sub.x_to_autosome, s=12.0,
                   marker=mk, color=INK, edgecolor="black",
                   linewidth=LINEWIDTH["hairline"], zorder=3,
                   label=labels[call])
    ax.set_ylabel("X : autosome depth")
    vals = d.x_to_autosome.dropna()
    lo = float(min(vals.min(), male_max)) if len(vals) else 0.0
    hi = float(max(vals.max(), female_min)) if len(vals) else 1.0
    pad = max(0.10 * (hi - lo), 0.05)
    ax.set_ylim(max(0.0, lo - pad), hi + pad)
    ax.set_xticks(x)
    ax.set_xticklabels([short(s) for s in d["run"]], rotation=90)
    ax.set_xlabel("Worm")
    ax.set_xlim(-0.7, len(d) - 0.3)
    ax.tick_params(axis="x", length=1.5, pad=1.5)
    return {"calls_selected": d.sex.value_counts().to_dict(),
            "calls_all_libraries": meta.sex.value_counts().to_dict(),
            "n_unknown_selected": int((d.sex == "unknown").sum()),
            "n_unknown_excluded": int(len(amb_excluded)),
            "n_drawn": int(len(d)),
            "not_drawn": [str(s) for s in drop["run"]]}


def panel_c(ax, dm: pd.DataFrame, order) -> dict:
    """Depth before and after matching, by region class.

    Bars are the ACHIEVED depth after subsampling — the number the caller
    actually saw — and the open markers are the input depth the worm arrived
    with. Drawing both is the point: the gap between marker and bar is how much
    was thrown away, and the flatness of the bars is the claim.

    This is a genuine region contrast, so it takes the region hues.
    """
    if dm.empty:
        ax.text(0.5, 0.5, "no depth_match tables found",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=FONT["annotation"], color=INK)
        ax.set_xticks([])
        ax.set_yticks([])
        return {}
    d = dm.set_index("sample").reindex([s for s in order
                                        if s in set(dm["sample"])]).reset_index()
    x = np.arange(len(d), dtype=float)
    w = 0.38
    ax.bar(x - w / 2, d.depth_retained_out, w, color=RETAINED,
           edgecolor="black", linewidth=LINEWIDTH["thin"], zorder=2,
           label="Retained")
    ax.bar(x + w / 2, d.depth_eliminated_out, w, color=ELIMINATED,
           edgecolor="black", linewidth=LINEWIDTH["thin"], zorder=2,
           label="Eliminated")
    ax.scatter(np.concatenate([x - w / 2, x + w / 2]),
               np.concatenate([d.depth_retained_in, d.depth_eliminated_in]),
               s=7.0, marker="_", color=INK, linewidth=LINEWIDTH["thin"],
               zorder=4, label="Before matching")
    target = float(pd.to_numeric(d.target, errors="coerce").median())
    if np.isfinite(target):
        ax.axhline(target, color=RULE, ls="--", lw=LINEWIDTH["rule"], zorder=3,
                   label=f"Target ({target:g}\u00d7)")
    ax.set_ylabel("Median depth (\u00d7)")
    # Limits from EVERY drawn value, including the pre-matching markers. Setting
    # them from the bars alone leaves the markers outside the axes, where
    # audit() cannot see them and a reader cannot either.
    drawn = np.concatenate([d.depth_retained_out, d.depth_eliminated_out,
                            d.depth_retained_in, d.depth_eliminated_in,
                            [target] if np.isfinite(target) else []])
    ax.set_ylim(0, float(np.nanmax(drawn)))
    headroom(ax, 0.12)
    ax.set_xticks(x)
    ax.set_xticklabels([short(s) for s in d["sample"]], rotation=90)
    ax.set_xlabel("Worm")
    ax.set_xlim(-0.7, len(d) - 0.3)
    ax.tick_params(axis="x", length=1.5, pad=1.5)

    ratio = (d.depth_eliminated_out / d.depth_retained_out).replace(
        [np.inf, -np.inf], np.nan)
    ratio_in = (d.depth_eliminated_in / d.depth_retained_in).replace(
        [np.inf, -np.inf], np.nan)
    # The RANGE of the pre-matching ratio, not just its median: the manuscript
    # quotes the span ("eliminated regions are covered at X% to Y% of retained
    # depth") and a median cannot source that sentence. Same for the achieved
    # depths — a worm whose input depth never reached the target is matched
    # BETWEEN classes but not TO the target, so "every worm at 13x" would be
    # wrong while "achieved ratio 1.00" stays right.
    short_of_target = d[d.depth_retained_out < 0.95 * target] \
        if np.isfinite(target) else d.iloc[0:0]
    return {"n_matched": int(len(d)),
            "target": target,
            "ratio_in_median": float(ratio_in.median()),
            "ratio_in_range": [float(ratio_in.min()), float(ratio_in.max())],
            "ratio_out_median": float(ratio.median()),
            "ratio_out_range": [float(ratio.min()), float(ratio.max())],
            "depth_out_range": [float(d.depth_retained_out.min()),
                                float(d.depth_retained_out.max())],
            "below_target": {str(smp): float(dep) for smp, dep in
                             zip(short_of_target["sample"],
                                 short_of_target["depth_retained_out"])},
            "worst_worm": str(d.loc[(ratio - 1).abs().idxmax(), "sample"])
                          if ratio.notna().any() else None}


def panel_d(ax, ab: pd.DataFrame) -> dict:
    """What the allele-balance screen removes, by region class.

    Grouped bars: x is the drop criterion, hue is region, y is log because the
    two classes differ by an order of magnitude and a linear axis would render
    the retained bars invisible. Identical layout to the standalone S16D so the
    two cannot disagree.

    EVERY BAR IS BASED ON THE AXIS FLOOR, NOT ON 0. Zero has no position on a
    log axis: matplotlib maps it to 10^-1000, so a bar based at 0 is written
    into the SVG and PDF with its bottom edge hundreds of inches off the page,
    hidden only by the axes clip path, and Illustrator then refuses to move or
    paste anything that includes it. Based on the lower y-limit, a bar ends at
    the x-axis and looks the same. The limits therefore have to be known before
    anything is drawn, so this works in two passes: read every region's row and
    set the limits from the six percentages, then draw. A percentage that is
    zero or missing has no position on the axis either; it is left undrawn,
    logged and listed in ``not_drawn``, rather than sent off the page.
    """
    x = np.arange(len(AB_GROUPS), dtype=float)
    bw = 0.34
    info = {}
    rows = {}

    # Pass 1 -- read. Nothing is drawn until the limits are known, because each
    # bar's base is the lower limit.
    for region in REGIONS:
        row = ab[ab.region == region]
        if row.empty:
            logging.error("allele_balance_screen.tsv has no '%s' row", region)
            continue
        row = row.iloc[0]
        rows[region] = row
        info[region] = {k: float(row[f"pct_{k}"]) for k, _lab in AB_GROUPS}
        info[region]["variant_sites"] = int(row.variant_sites)
        # Sites with no heterozygote carrying usable AD cannot be judged on
        # balance and sit in the denominator, so a high rate here means the
        # drawn percentages UNDERSTATE what the screen would remove. It is
        # higher in eliminated regions, which is the direction that matters.
        if "no_ad_evaluable" in row.index:
            n_no_ad = int(row.no_ad_evaluable)
            info[region]["no_ad_evaluable"] = n_no_ad
            info[region]["pct_no_ad"] = (100 * n_no_ad / int(row.variant_sites)
                                         if int(row.variant_sites) else float("nan"))

    # Limits from the six DRAWN percentages only. info[region] also carries
    # variant_sites, no_ad_evaluable (a raw count of sites, in the thousands)
    # and pct_no_ad, none of which is drawn; taking the limits from every value
    # in it let the no-AD count set the top, which is how the 2026-08-22 build
    # came to run to 10^4 %.
    drawn = [info[r][k] for r in REGIONS if r in info
             for k, _lab in AB_GROUPS
             if np.isfinite(info[r][k]) and info[r][k] > 0]
    limits = None
    if drawn:
        # Multiplicative headroom: headroom() does linear arithmetic and would
        # not reserve the space it claims to on a log axis. The top stops at
        # 100 %, the most a share can be.
        limits = (min(drawn) / 3.0, min(PCT_MAX, max(drawn) * 6.0))
    else:
        logging.error("D  allele_balance_screen.tsv has no positive "
                      "percentage to draw; panel D is empty")

    # Pass 2 -- draw. bar() measures height from `bottom`, so each bar is
    # (value - floor) tall and its top lands on the value, where the label sits.
    not_drawn = []
    for j, region in enumerate(REGIONS):
        if region not in rows:
            continue
        vals = np.array([info[region][k] for k, _lab in AB_GROUPS])
        pos = x + (j - 0.5) * bw
        ok = np.isfinite(vals) & (vals > 0)
        skipped = [f"{region} pct_{k} = {v:g}"
                   for (k, _lab), v, good in zip(AB_GROUPS, vals, ok)
                   if not good]
        if skipped:
            logging.warning("D  not drawn, no position on a log axis: %s",
                            ", ".join(skipped))
            not_drawn += skipped
        if limits is None or not ok.any():
            continue
        floor = limits[0]
        ax.bar(pos[ok], vals[ok] - floor, bw, bottom=floor,
               color=paper_colors.REGION[region],
               edgecolor="black", linewidth=LINEWIDTH["thin"], zorder=2,
               label=f"{region.capitalize()} "
                     f"(n = {int(rows[region].variant_sites):,})")
        for xi, v in zip(pos[ok], vals[ok]):
            ax.annotate(f"{v:.2f}", (xi, v), textcoords="offset points",
                        xytext=(0, 2), ha="center",
                        fontsize=FONT["annotation"], zorder=3)
    info["not_drawn"] = not_drawn

    ax.set_yscale("log")
    ax.set_ylabel("Variant sites removed (%)")
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _k, lab in AB_GROUPS])
    ax.set_xlim(-0.55, len(AB_GROUPS) - 0.45)
    if limits is not None:
        ax.set_ylim(*limits)
        info["y_limits_pct"] = [float(v) for v in ax.get_ylim()]
        info["bar_base_pct"] = float(limits[0])
    plain_log_ticks(ax.yaxis)
    if {"retained", "eliminated"} <= set(info):
        r, e = info["retained"]["dropped"], info["eliminated"]["dropped"]
        info["fold_eliminated_over_retained"] = (e / r if r else float("nan"))
    return info


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Supplementary fig. S14 (stem S16) — cohort selection and callset "
                    "corrections.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--metadata", required=True, type=Path,
                   help="sample_metadata.tsv from select_germline_samples.py, "
                        "or a directory to search for it")
    p.add_argument("--depth-match-glob", required=True,
                   help="Glob for the per-worm <sample>.depth_match.tsv files "
                        "written by depth_match.py, or a directory to search "
                        "recursively for them. QUOTE A GLOB.")
    p.add_argument("--allele-balance", required=True, type=Path,
                   help="allele_balance_screen.tsv from allele_balance_panel.py, "
                        "or a directory to search for it (it lands in a dated "
                        "make_outdir() subdirectory, so pass the parent)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-germline", type=float, default=0.40,
                   help="Match select_germline_samples.py; drawn in panel A")
    p.add_argument("--male-max", type=float, default=0.65,
                   help="Match select_germline_samples.py; lower edge of the "
                        "ambiguous band in panel B")
    p.add_argument("--female-min", type=float, default=0.90,
                   help="Match select_germline_samples.py; upper edge")
    p.add_argument("--tag", default=None)
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    for noisy in ("fontTools", "matplotlib.font_manager", "matplotlib"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    apply_style()

    meta_path = resolve_input(args.metadata, "sample_metadata.tsv")
    meta = pd.read_csv(meta_path, sep="\t")
    # select_germline_samples.py names the library column `run`; the audit table
    # it was built from names the same thing `sample`. Accept either rather than
    # failing on a table that has the column under the other name.
    if "run" not in meta.columns and "sample" in meta.columns:
        logging.info("No 'run' column; using 'sample' as the library name")
        meta = meta.rename(columns={"sample": "run"})
    need = {"run", "sex", "x_to_autosome", "germline_fraction", "pass",
            "is_control"}
    missing = need - set(meta.columns)
    if missing:
        logging.error("sample_metadata.tsv is missing column(s): %s",
                      sorted(missing))
        logging.error("Columns present: %s", list(meta.columns))
        return 2
    meta["run"] = meta["run"].map(strip_suffix)
    meta["is_control"] = meta["is_control"].astype(bool)
    meta["pass"] = meta["pass"].astype(bool)
    if "fail_reason" not in meta.columns:
        meta["fail_reason"] = ""
    meta["fail_reason"] = meta["fail_reason"].fillna("")
    # Worms first, ranked by the quantity the selection was made on, then the
    # controls. Every panel uses this order so a reader can follow one library
    # across A, B and C.
    meta = meta.sort_values(["is_control", "germline_fraction"],
                            ascending=[True, False]).reset_index(drop=True)

    dm = load_depth_match(args.depth_match_glob)
    if dm.empty:
        logging.error("No depth_match tables matched %r.",
                      args.depth_match_glob)
        logging.error("Find them with:  find $POPGEN_DIR -name "
                      "'*.depth_match.tsv' | head")
        logging.error("Panel C cannot be drawn without them, and a supplement "
                      "with an empty panel should not be assembled. Re-run "
                      "with the right --depth-match-glob.")
        return 1

    ab_path = resolve_input(args.allele_balance, "allele_balance_screen.tsv")
    ab = pd.read_csv(ab_path, sep="\t")
    if "region" not in ab.columns:
        logging.error("allele_balance_screen.tsv has no 'region' column; got %s",
                      list(ab.columns))
        return 2

    # The tables come from three separate runs in three output directories, so
    # assembling a stale one in is the obvious failure mode and it would be
    # completely invisible in the drawn figure.
    selected = set(meta.loc[(~meta.is_control) & meta["pass"], "run"])
    matched = set(dm["sample"])
    if matched != selected:
        logging.error("The selected cohort and the depth-matched set differ:\n"
                      "  selected      %d: %s\n  depth-matched %d: %s",
                      len(selected), sorted(selected),
                      len(matched), sorted(matched))
        logging.error("  selected but not matched: %s",
                      sorted(selected - matched) or "none")
        logging.error("  matched but not selected: %s",
                      sorted(matched - selected) or "none")
        logging.error("Re-run depth_match.py over the current selection rather "
                      "than assembling a figure from mismatched runs.")
        return 1
    logging.info("Cohort agrees between metadata and depth_match: %d worms",
                 len(selected))

    outdir = Path(make_outdir("s16_cohort_selection",
                              base=str(args.outdir.expanduser().resolve()),
                              tag=args.tag))
    logging.info("Output: %s", outdir)

    fig, axes = figure(FIG_W, height_mm=FIG_H_MM, nrows=2, ncols=2)
    (axA, axB), (axC, axD) = axes
    a_info = panel_a(axA, meta, args.min_germline)
    b_info = panel_b(axB, meta, args.male_max, args.female_min)
    c_info = panel_c(axC, dm, list(meta["run"]))
    d_info = panel_d(axD, ab)

    # Panel A's fill encoding has no drawn handle to build a legend from — the
    # bars ARE the categories — so it gets proxy handles. Everything else builds
    # its own.
    from matplotlib.patches import Patch
    a_handles, a_labels = axA.get_legend_handles_labels()
    a_handles = [Patch(facecolor=INK, edgecolor=INK,
                       linewidth=LINEWIDTH["thin"]),
                 Patch(facecolor="white", edgecolor=INK,
                       linewidth=LINEWIDTH["thin"]),
                 Patch(facecolor=INK, alpha=0.35, edgecolor=INK,
                       hatch=CONTROL_HATCH, linewidth=LINEWIDTH["thin"])
                 ] + a_handles
    a_labels = ["Selected", "Excluded", "Tissue control"] + a_labels
    axA.legend(a_handles, a_labels, loc="lower left", bbox_to_anchor=(0.0, 1.02),
               ncol=4, frameon=False, fontsize=FONT["legend"], handlelength=1.4,
               handletextpad=0.5, columnspacing=1.0, borderaxespad=0.0)

    # Every other legend goes ABOVE its axes too. Placing a legend in a corner
    # means guessing which corner the data leaves free, and that guess depends
    # on the data.
    for a, ncol in ((axB, 4), (axC, 4), (axD, 2)):
        if a.get_legend_handles_labels()[0]:
            a.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=ncol,
                     frameon=False, fontsize=FONT["legend"], handlelength=1.4,
                     handletextpad=0.5, columnspacing=1.0, borderaxespad=0.0)

    for a, letter in ((axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")):
        panel_letter(a, letter, dy=1.20)

    problems = audit(fig)
    if problems:
        logging.warning("S14 breaks %d journal rule(s) — listed above.",
                        len(problems))
    written = pc_save(fig, str(outdir), "S16_cohort_selection")
    plt.close(fig)

    (outdir / "provenance.json").write_text(json.dumps({
        "argv": sys.argv,
        "inputs": {"metadata": str(meta_path),
                   "depth_match_glob": args.depth_match_glob,
                   "allele_balance": str(ab_path)},
        "thresholds": {"min_germline": args.min_germline,
                       "male_max": args.male_max,
                       "female_min": args.female_min},
        "panel_a": a_info, "panel_b": b_info, "panel_c": c_info,
        "panel_d": d_info, "audit_problems": problems,
    }, indent=2, default=str))

    logging.info("=" * 74)
    logging.info("A  %d librar(ies): %d worms, %d selected; excluded %s",
                 a_info["n_libraries"], a_info["n_worms"], a_info["n_pass"],
                 a_info["fail_reasons"] or "none")
    if a_info["control_range"]:
        logging.info("   tissue controls span germline fraction %.3f-%.3f — "
                     "the assay floor and ceiling germline_rescaled is defined "
                     "against", *a_info["control_range"])
    over_one = meta[meta["germline_fraction"] > 1.0]
    if len(over_one):
        # For the legend: the fraction is a depth ratio, not a proportion.
        logging.info("   LEGEND: germline fraction above 1 in %s (%s) — it is a "
                     "depth ratio (eliminated : retained), not a proportion, so "
                     "a germline-only tissue can exceed 1",
                     ", ".join(short(r) for r in over_one["run"]),
                     ", ".join(f"{v:.2f}" for v in over_one["germline_fraction"]))
    logging.info("B  sex calls in the selected cohort: %s%s",
                 b_info["calls_selected"],
                 "" if not b_info["n_unknown_excluded"] else
                 f"  ({b_info['n_unknown_excluded']} further ambiguous call(s) "
                 f"among libraries panel A already excluded — not drawn)")
    if b_info.get("n_unknown_selected"):
        logging.warning("B  %d SELECTED worm(s) have an ambiguous X:A. A worm "
                        "carried into the analysis without a sex call is a "
                        "real loose end, not a cosmetic one: the X is retained "
                        "and hemizygous in males, so it sits in the retained "
                        "class of every pooled estimate.",
                        b_info["n_unknown_selected"])
    if c_info:
        logging.info("C  eliminated:retained depth ratio %.3f before matching "
                     "(range %.3f-%.3f), %.3f after (range %.3f-%.3f) across "
                     "%d worms",
                     c_info["ratio_in_median"], c_info["ratio_in_range"][0],
                     c_info["ratio_in_range"][1], c_info["ratio_out_median"],
                     c_info["ratio_out_range"][0], c_info["ratio_out_range"][1],
                     c_info["n_matched"])
        if c_info.get("below_target"):
            logging.info("   %d worm(s) never reached the %gx target and were "
                         "matched at their own depth instead: %s. Both classes "
                         "still agree in those worms — the ratio is what the "
                         "contrast depends on — but the cohort is NOT uniformly "
                         "at %gx, so do not write it that way.",
                         len(c_info["below_target"]), c_info["target"],
                         ", ".join(f"{k} at {v:.1f}x" for k, v in
                                   c_info["below_target"].items()),
                         c_info["target"])
    if "fold_eliminated_over_retained" in d_info:
        logging.info("D  the screen removes %.2f%% of retained and %.2f%% of "
                     "eliminated variant sites (%.1fx)",
                     d_info["retained"]["dropped"],
                     d_info["eliminated"]["dropped"],
                     d_info["fold_eliminated_over_retained"])
        no_ad = {r: d_info[r].get("pct_no_ad") for r in REGIONS
                 if r in d_info and d_info[r].get("pct_no_ad") is not None}
        if no_ad:
            logging.info("   %s of sites carried no heterozygote with usable "
                         "AD and could not be judged on balance. They are in "
                         "the denominator, so these percentages are a FLOOR on "
                         "what the screen would remove — and the rate is higher "
                         "in eliminated regions.",
                         ", ".join(f"{v:.1f}% {r}" for r, v in no_ad.items()))
    if d_info.get("not_drawn"):
        logging.warning("D  %d percentage(s) not drawn, being zero or missing: "
                        "%s", len(d_info["not_drawn"]),
                        ", ".join(d_info["not_drawn"]))
    logging.info("")
    logging.info("D is the load-bearing panel. A and B choose the cohort and C "
                 "removes the depth asymmetry, but D is the control that would "
                 "most plausibly have produced the headline as an artifact — "
                 "and the screen is applied, so what it removes is gone from "
                 "the callset the result is computed on.")
    logging.info("Figure: %s", ", ".join(Path(w).name for w in written))
    logging.info("Export from the SVG — Arial is absent on Isaac nodes.")
    logging.info("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
