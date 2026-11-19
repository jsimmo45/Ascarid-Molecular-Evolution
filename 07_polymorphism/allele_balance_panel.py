#!/usr/bin/env python3
"""
allele_balance_panel.py — supplementary fig. S16D: what the allele-balance
screen removes, by region class.

Why this deserves a panel
-------------------------
The allele-balance screen is the single most load-bearing control in the
polymorphism analysis, and it currently exists only as a number in a log file.
It removes variant sites whose heterozygous calls are not balanced, or where
every called genotype is heterozygous — both signatures of a duplication present
in these *A. lumbricoides* worms but collapsed in the *A. suum* reference. A
mappability mask cannot find these, because a mask knows only reference
paralogy.

The screen removes far more eliminated variants than retained ones: on the
2026-08-20 callset (at most one missing genotype) the filter's own log shows
18.5% of eliminated variant sites removed against 1.55% of retained, 12-fold.
Those variants would otherwise have entered piN/piS as polymorphism, inflating
eliminated diversity in exactly the direction that resembles the result. A
referee who does not see that number quantified will ask for it.

What this script does differently from reading the filter log
--------------------------------------------------------------
It recomputes the screen from the PRE-screen callset rather than parsing
`filter_genotypes.py`'s summary, so the panel is reproducible from data and the
two drop criteria can be reported separately. Sites failing both are counted and
reported rather than being silently assigned to one category. It then checks
every verdict against the screen's own record of what it removed,
`excluded_sites.txt` in the same directory as the callset, so the panel cannot
drift from the filter without saying so.

Figure conventions
------------------
Built to the house spec: `paper_colors.apply_style()`, exact column width in
millimetres, `audit()` before saving, dated output directory, provenance record.
Hue carries region, the locked blue/red convention, and the drop criterion runs
along the x-axis — so the two classes sit adjacent within each criterion and the
comparison the panel exists to make is a direct one.

The y-axis is log. On the 2026-08-20 callset retained lost 1.55% of variant
sites and eliminated 18.5%, so on a linear axis the retained bars are under a
millimetre tall at print size and the split by criterion is invisible. Bars are
grouped rather than stacked: the two criteria overlap, so their parts do not
sum, and stacking on a log axis would misrepresent them twice over.

Changes 2026-10-12
------------------
The log axis stops at 100 % and its ticks print as plain numbers (0.1, 1, 10,
100), matching panel D of figure_cohort_selection.py, whose own axis ran to
10^4 % because it took its limits from a site count as well as the drawn
percentages. This script was not affected by that; the change only keeps the
two panels identical. Nothing else moves, so the 2026-08-22 table and
standalone figure stay valid.

Changes 2026-10-14
------------------
The screen is now recomputed with filter_genotypes.py's rule term for term.
The 2026-08-22 table, which S14D drew, is superseded:

* A site is judged only if it has at least ``--min-het`` (3) heterozygous
  calls whose AD gives the two called alleles a non-zero total, and is kept
  otherwise. This script used to judge every site with one such
  heterozygote, and about 9% of true heterozygotes at 13 reads fall outside
  0.25-0.75 by binomial sampling alone, so a site with a single heterozygote,
  the commonest kind, was counted as removed about 9% of the time by chance.
  The 2026-08-22 table summed to about 11,600 removed sites; the applied
  screen removed at most 6,788 records.
* "All genotypes heterozygous" is the filter's test: usable heterozygotes are
  at least ``--max-het-frac`` (1.0) of called genotypes, where a call is any
  non-missing GT.
* Verdicts are checked against ``excluded_sites.txt``. Any disagreement is an
  error (exit 3) after the tables and figure are written, so a chained S14
  rebuild does not run on a panel that disagrees with the callset.
* The table drops ``no_ad_evaluable``. It counted sites with no usable
  heterozygote, which the screen keeps; figure_cohort_selection.py read it as
  the denominator of a "floor" that does not exist. Sites the screen does not
  judge are now ``not_judged``. The percentage denominators are unchanged:
  every biallelic SNP in classified CDS.
* The standalone bars start at the axis floor rather than at 0, which a log
  axis maps to 10^-1000; figure_cohort_selection.py's panel D was fixed the
  same way on 2026-10-13 (Illustrator could not move or paste it).

Example
-------
    python allele_balance_panel.py \\
        --vcf vcf/filtered_miss4/step3_gapfilter.bcf \\
        --gff As_miniprot.gff \\
        --gene-callable targets/gene_callable.tsv \\
        --outdir results/figures --tag filtered_miss4

``excluded_sites.txt`` is read from beside ``--vcf`` unless
``--excluded-sites`` names it. See run_allele_balance_panel.sh.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

# het_by_region sets its own rcParams at import time (it predates the house
# style), so paper_colors.apply_style() is called AFTER these imports and wins.
from het_by_region import RegionMap, vcf_samples
from pi_by_region import read_cds

def _find_paper_colors() -> Path:
    """
    Locate paper_colors.py.

    It lives in the Ka/Ks tree rather than beside these scripts, and its exact
    depth differs between checkouts, so guessing one relative path fails on a
    compute node with an import error that says nothing useful. Search the
    plausible locations, and if none has it, say exactly what was tried.
    """
    here = Path(__file__).resolve().parent
    candidates = []
    env = os.environ.get("PAPER_COLORS_DIR")
    if env:
        candidates.append(Path(env).expanduser())
    for base in (here, here.parent, here.parent.parent):
        candidates += [base, base / "lib", base / "scripts", base / "src",
                       base / "Ascarid-Molecular-Evolution" / "lib"]
    seen, ordered = set(), []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    for c in ordered:
        if (c / "paper_colors.py").is_file():
            return c
    raise SystemExit(
        "paper_colors.py not found. Set PAPER_COLORS_DIR to the directory "
        "holding it, e.g.\n"
        "    export PAPER_COLORS_DIR=/path/to/"
        "Ascarid-Molecular-Evolution/lib\n"
        "Searched:\n  " + "\n  ".join(str(c) for c in ordered))


sys.path.insert(0, str(_find_paper_colors()))
from paper_colors import (REGION, NEUTRAL, FONT, LINEWIDTH, apply_style, audit,
                          figure, legend_outside, make_outdir, panel_letter,
                          save)

REGIONS = ("retained", "eliminated")
RCODE = {"retained": 0, "eliminated": 1}
CODE_NAME = {v: k for k, v in RCODE.items()}

#: The panel plots shares of variant sites, and a share cannot exceed 100 %, so
#: the log axis stops there. Same value in figure_cohort_selection.py.
PCT_MAX = 100.0


def plain_log_ticks(axis) -> None:
    """Label a log axis 0.1, 1, 10, 100 rather than 10^-1 ... 10^2.

    The axis is a percentage and should read as one. The default exponent
    labels are also mathtext, which the SVG writes glyph by glyph in the math
    font rather than the label typeface. Minor ticks stay, unlabelled.
    Identical in figure_cohort_selection.py.
    """
    axis.set_major_locator(LogLocator(base=10.0))
    axis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:g}"))
    axis.set_minor_formatter(NullFormatter())


#: How many disagreements with excluded_sites.txt to print, per direction.
N_EXAMPLES = 10


def read_excluded(path: Path) -> set:
    """
    Load filter_genotypes.py's ``excluded_sites.txt`` as a set of (chrom, pos).

    The filter writes one ``CHROM<TAB>POS`` line (POS 1-based) per variant
    record its allele-balance screen dropped, and then removes every record at
    those positions with ``bcftools view -T ^excluded_sites.txt``. It is the
    only record of what the applied screen removed, and the filter's log was not
    kept for every callset, so it is the reference this script is checked
    against.
    """
    sites = set()
    n_lines = n_bad = 0
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            n_lines += 1
            f = line.rstrip("\n").split("\t")
            try:
                sites.add((f[0], int(f[1])))
            except (IndexError, ValueError):
                n_bad += 1
    if n_bad:
        logging.error("%s: %d of %d line(s) are not CHROM<TAB>POS; is this the "
                      "file filter_genotypes.py wrote?", path, n_bad, n_lines)
        sys.exit(2)
    logging.info("excluded_sites.txt: %d line(s), %d distinct position(s): %s",
                 n_lines, len(sites), path)
    return sites


def site_verdict(cells: list[tuple[str, str]], ab_lo: float, ab_hi: float,
                 min_het: int, max_het_frac: float) -> dict:
    """
    Apply filter_genotypes.py's allele-balance rule to one site.

    ``cells`` is one (GT, AD) pair per sample. Every step mirrors
    ``filter_genotypes.allele_balance_screen`` (2026-08-20), so the verdict is
    the one the callset was filtered with, not an approximation of it:

    * A genotype is CALLED unless its GT is ``./.``, ``.|.`` or ``.``.
    * A called genotype is a USABLE heterozygote if its two alleles differ, AD
      has a depth for each of them, and those two depths sum above zero. Its
      balance is the depth of the second called allele over the sum, which for
      a biallelic 0/1 call is alt / (ref + alt).
    * A site with fewer than ``min_het`` usable heterozygotes is NOT JUDGED and
      kept, whatever its genotypes look like.
    * A judged site fails on BALANCE if the mean over its usable heterozygotes
      lies outside [ab_lo, ab_hi], and on HETEROZYGOSITY if usable
      heterozygotes make up at least ``max_het_frac`` of its called genotypes.
      The filter tests balance first and logs a site failing both as "skewed";
      both verdicts are returned here so the overlap can be reported.
    """
    n_called = n_het = 0
    ab_sum = 0.0
    for gt, ad in cells:
        if gt in ("./.", ".|.", "."):
            continue
        n_called += 1
        a = gt.replace("|", "/").split("/")
        if len(a) != 2 or a[0] == a[1]:
            continue
        counts = [int(x) for x in ad.split(",") if x not in (".", "")]
        if len(counts) < 2 or sum(counts) == 0:
            continue
        try:
            d1, d2 = counts[int(a[0])], counts[int(a[1])]
        except (ValueError, IndexError):
            continue
        if d1 + d2 == 0:
            continue
        ab_sum += d2 / (d1 + d2)
        n_het += 1

    if n_het < min_het:
        return {"judged": False, "bad_ab": False, "all_het": False,
                "n_het": n_het, "n_called": n_called, "mean_ab": float("nan")}
    mean_ab = ab_sum / n_het
    return {"judged": True,
            "bad_ab": not (ab_lo <= mean_ab <= ab_hi),
            "all_het": n_called > 0 and (n_het / n_called) >= max_het_frac,
            "n_het": n_het, "n_called": n_called, "mean_ab": mean_ab}


def classify_sites(vcf: Path, rmap: RegionMap, n: int, ab_lo: float,
                   ab_hi: float, min_het: int, max_het_frac: float,
                   excluded: set | None) -> dict:
    """
    Stream the pre-screen callset and apply the screen to every biallelic SNP
    in classified CDS (the sites piN/piS is computed from).

    Returns per region code: variant sites, sites judged, sites failing the
    balance window, sites with every called genotype heterozygous, and sites
    failing both. The criteria overlap, so ``both`` is reported rather than
    absorbed into whichever was tested first.

    With ``excluded`` (the screen's excluded_sites.txt as a set), each site's
    verdict is compared with whether the screen listed it. The two disagree only
    if the rule here is not the rule the callset was filtered with, the
    parameters differ, or the two files come from different runs.
    """
    fmt = "%CHROM\\t%POS[\\t%GT\\t%AD]\\n"
    cmd = (f"bcftools query -f '{fmt}' "
           f"-i 'N_ALT=1 && TYPE=\"snp\"' {vcf}")
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 22)

    acc = {c: {"total": 0, "judged": 0, "ab": 0, "allhet": 0, "both": 0}
           for c in RCODE.values()}
    xc = {c: {"listed": 0, "agree_drop": 0, "agree_keep": 0,
              "panel_only": 0, "list_only": 0} for c in RCODE.values()}
    examples = {"panel_only": [], "list_only": []}
    n_lines = n_mapped = n_short = 0

    for line in proc.stdout:
        n_lines += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 + 2 * n:
            n_short += 1
            continue
        chrom, pos = f[0], int(f[1])
        code = rmap.lookup(chrom, pos - 1)
        if code < 0:
            continue
        n_mapped += 1
        a = acc[code]
        a["total"] += 1

        cells = [(f[2 + 2 * k], f[3 + 2 * k]) for k in range(n)]
        v = site_verdict(cells, ab_lo, ab_hi, min_het, max_het_frac)
        if v["judged"]:
            a["judged"] += 1
            if v["bad_ab"]:
                a["ab"] += 1
            if v["all_het"]:
                a["allhet"] += 1
            if v["bad_ab"] and v["all_het"]:
                a["both"] += 1
        dropped = v["bad_ab"] or v["all_het"]

        if excluded is not None:
            x = xc[code]
            listed = (chrom, pos) in excluded
            x["listed"] += listed
            if dropped and listed:
                x["agree_drop"] += 1
            elif not dropped and not listed:
                x["agree_keep"] += 1
            else:
                kind = "panel_only" if dropped else "list_only"
                x[kind] += 1
                if len(examples[kind]) < N_EXAMPLES:
                    examples[kind].append(
                        f"{chrom}:{pos} {CODE_NAME[code]} het={v['n_het']} "
                        f"called={v['n_called']} mean_ab={v['mean_ab']:.3f}")

    proc.wait()
    if proc.returncode != 0:
        logging.error("bcftools query failed:\n%s", proc.stderr.read()[-2000:])
        sys.exit(proc.returncode)

    logging.info("Streamed %d biallelic SNP record(s); %d in classified CDS",
                 n_lines, n_mapped)
    if n_short:
        logging.warning("%d record(s) had fewer fields than %d samples need "
                        "and were skipped", n_short, n)
    if n_mapped == 0:
        logging.error("No variants mapped into the region map — check that the "
                      "GFF contig names match the VCF")
        sys.exit(1)
    return {"counts": acc,
            "crosscheck": xc if excluded is not None else None,
            "examples": examples if excluded is not None else None}


def build_table(acc: dict, xc: dict | None, min_het: int) -> pd.DataFrame:
    """One row per region class; percentages are of every variant site."""
    rows = []
    for code, a in acc.items():
        total = a["total"]
        dropped = a["ab"] + a["allhet"] - a["both"]
        pct = (lambda k: 100 * k / total if total else np.nan)
        row = {
            "region": CODE_NAME[code],
            "variant_sites": total,
            "judged": a["judged"],
            "not_judged": total - a["judged"],
            "drop_ab": a["ab"],
            "drop_allhet": a["allhet"],
            "drop_both": a["both"],
            "dropped_total": dropped,
            "pct_ab": pct(a["ab"]),
            "pct_allhet": pct(a["allhet"]),
            "pct_dropped": pct(dropped),
            "pct_not_judged": pct(total - a["judged"]),
            "min_het": min_het,
        }
        if xc is not None:
            x = xc[code]
            row["listed_in_excluded_sites"] = x["listed"]
            row["disagree_panel_only"] = x["panel_only"]
            row["disagree_list_only"] = x["list_only"]
        rows.append(row)
    return pd.DataFrame(rows).set_index("region").reindex(REGIONS).reset_index()


#: x-axis groups: the drop criterion. Region is the hue, so the two classes sit
#: adjacent within each group and the comparison the panel exists to make is a
#: direct one. An earlier layout put region on the x-axis and carried the
#: criterion in hatch texture; it needed a two-entry legend of long labels that
#: `audit()` flagged as overlapping the bars in an 89 mm panel.
GROUPS = [("ab", "Allele balance\noutside window"),
          ("allhet", "All genotypes\nheterozygous"),
          ("dropped", "Either\n(total removed)")]


def draw(df: pd.DataFrame, width: str, height_mm: float, letter: str | None):
    """
    Grouped bars: x = drop criterion, hue = region, log y.

    The limits are set before anything is drawn, and each bar runs from the
    axis floor up to its value. A bar drawn from 0 on a log axis starts at
    10^-1000: matplotlib clips it, but the SVG path still runs tens of thousands
    of points off the page, and Illustrator then refuses to move or paste the
    panel. A share of 0 (or a missing one) cannot be placed on a log axis at
    all, so it is left undrawn and reported. Same rule as panel D of
    figure_cohort_selection.py (2026-10-13).

    Returns the figure, the axes and a dict for the provenance record:
    ``y_limits_pct``, ``bar_base_pct`` and ``not_drawn``.
    """
    fig, ax = figure(width, height_mm=height_mm)
    info = {"not_drawn": []}

    vals_by_region = {}
    for region in REGIONS:
        row = df[df.region == region].iloc[0]
        vals_by_region[region] = [float(row[f"pct_{k}"]) for k, _lab in GROUPS]
    pos_vals = [v for vals in vals_by_region.values() for v in vals
                if np.isfinite(v) and v > 0]

    ax.set_yscale("log")
    if pos_vals:
        # Multiplicative headroom: this axis is log, so paper_colors.headroom()'s
        # linear arithmetic would not reserve the space it claims to. The top
        # stops at 100 %, the most a share can be.
        ax.set_ylim(min(pos_vals) / 3.0, min(PCT_MAX, max(pos_vals) * 6.0))
    else:
        logging.warning("No positive share to draw; panel D will be empty.")
        ax.set_ylim(0.1, PCT_MAX)
    floor = ax.get_ylim()[0]
    info["y_limits_pct"] = [float(v) for v in ax.get_ylim()]
    info["bar_base_pct"] = float(floor)

    x = np.arange(len(GROUPS), dtype=float)
    bw = 0.34
    for j, region in enumerate(REGIONS):
        row = df[df.region == region].iloc[0]
        pos = x + (j - 0.5) * bw
        label = f"{region.capitalize()} (n = {int(row.variant_sites):,})"
        for i, (xi, v) in enumerate(zip(pos, vals_by_region[region])):
            key = GROUPS[i][0]
            if not (np.isfinite(v) and v > 0):
                info["not_drawn"].append({"region": region, "group": key,
                                          "pct": None if v != v else v})
                logging.warning("%s %s share is %s: not drawn (a log axis has "
                                "no place for it); state it in the legend.",
                                region, key, "missing" if v != v else f"{v:g}%")
                continue
            ax.bar(xi, v - floor, bw, bottom=floor, color=REGION[region],
                   edgecolor="black", linewidth=LINEWIDTH["thin"], zorder=2,
                   label=label)
            label = None          # one legend entry per region
            ax.annotate(f"{v:.2f}", (xi, v), textcoords="offset points",
                        xytext=(0, 2), ha="center",
                        fontsize=FONT["annotation"], zorder=3)

    ax.set_ylabel("Variant sites removed (%)")
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _k, lab in GROUPS])
    ax.set_xlim(-0.55, len(GROUPS) - 0.45)
    plain_log_ticks(ax.yaxis)

    # Legend above the axes, where no data can reach it. n goes in the legend
    # rather than the tick labels here because the ticks name the criterion,
    # not the gene class.
    legend_outside(ax, "top", ncol=2)
    if letter:
        panel_letter(ax, letter)
    return fig, ax, info


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Supplementary fig. S16D — allele-balance screen by region.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf", required=True, type=Path,
                   help="PRE-screen callset (step3_gapfilter.bcf). Passing the "
                        "post-screen file measures nothing — the sites this "
                        "panel counts have already been removed from it.")
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--gene-callable", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--excluded-sites", type=Path, default=None,
                   help="filter_genotypes.py's excluded_sites.txt for this "
                        "callset; every verdict is checked against it. Default: "
                        "excluded_sites.txt beside --vcf (skipped, with a "
                        "warning, if absent)")
    p.add_argument("--ab-min", type=float, default=0.25,
                   help="Lower edge of the accepted allele-balance window; "
                        "must match filter_genotypes.py")
    p.add_argument("--ab-max", type=float, default=0.75)
    p.add_argument("--min-het", type=int, default=3,
                   help="Judge a site only with at least this many usable "
                        "heterozygous calls; must match filter_genotypes.py")
    p.add_argument("--max-het-frac", type=float, default=1.0,
                   help="Drop a judged site when usable heterozygotes are at "
                        "least this fraction of called genotypes; must match "
                        "filter_genotypes.py")
    p.add_argument("--width", default="single",
                   choices=["single", "onehalf", "double"],
                   help="Column width. 'single' (89 mm) so two panels sit side "
                        "by side in the S16 2x2 assembly without scaling.")
    p.add_argument("--height-mm", type=float, default=58.0)
    p.add_argument("--panel-letter", default="D",
                   help="Empty string to omit")
    p.add_argument("--tag", default=None, help="Suffix for the output directory")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = p.parse_args(argv)
    if args.min_het < 1:
        p.error("--min-het must be at least 1")
    if not 0 < args.max_het_frac <= 1:
        p.error("--max-het-frac must be in (0, 1]")
    if not 0 <= args.ab_min < args.ab_max <= 1:
        p.error("need 0 <= --ab-min < --ab-max <= 1")
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    # matplotlib's PDF writer subsets fonts through fontTools, which logs a few
    # hundred INFO lines per save and buries the result table this script exists
    # to print.
    for noisy in ("fontTools", "matplotlib.font_manager", "matplotlib"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    apply_style()

    vcf = args.vcf.expanduser().resolve()
    for label, path in (("--vcf", vcf), ("--gff", args.gff.expanduser()),
                        ("--gene-callable", args.gene_callable.expanduser())):
        if not path.is_file():
            logging.error("%s not found: %s", label, path)
            return 2
    if "step3" not in vcf.name:
        logging.warning("%s is not named step3_*: this panel needs the callset "
                        "BEFORE the allele-balance screen. On the screened "
                        "callset every count below is ~0 by construction.",
                        vcf.name)

    # The screen's own record of what it removed. Look beside the callset unless
    # told otherwise; without it the panel still draws, but nothing confirms
    # that it reports the screen that was applied.
    excl_path = (args.excluded_sites.expanduser().resolve()
                 if args.excluded_sites else vcf.parent / "excluded_sites.txt")
    if excl_path.is_file():
        excluded = read_excluded(excl_path)
    elif args.excluded_sites:
        logging.error("--excluded-sites not found: %s", excl_path)
        return 2
    else:
        logging.warning("No excluded_sites.txt beside %s: the panel cannot be "
                        "checked against the applied screen.", vcf)
        excluded, excl_path = None, None

    outdir = make_outdir("s16d_allele_balance",
                         base=str(args.outdir.expanduser().resolve()),
                         tag=args.tag)
    logging.info("Output: %s", outdir)
    logging.info("Screen: judged at >= %d usable heterozygote(s); drop if mean "
                 "balance outside [%.2f, %.2f], or if usable heterozygotes are "
                 ">= %.2f of called genotypes", args.min_het, args.ab_min,
                 args.ab_max, args.max_het_frac)

    samples = vcf_samples(vcf)
    n = len(samples)
    logging.info("Samples: %d", n)
    if n == 0:
        logging.error("No samples in %s", vcf)
        return 2

    gc = pd.read_csv(args.gene_callable.expanduser(), sep="\t")
    gc = gc[gc["pass"].astype(bool) & gc.region.isin(REGIONS)]
    region_of = dict(zip(gc.gene, gc.region))
    logging.info("Callable genes: %s", gc.region.value_counts().to_dict())
    genes = read_cds(args.gff.expanduser(), set(gc.gene))
    rmap = RegionMap()
    if rmap.build(genes, region_of) == 0:
        logging.error("Empty region map")
        return 1

    res = classify_sites(vcf, rmap, n, args.ab_min, args.ab_max,
                         args.min_het, args.max_het_frac, excluded)
    df = build_table(res["counts"], res["crosscheck"], args.min_het)
    df.to_csv(Path(outdir) / "allele_balance_screen.tsv", sep="\t", index=False)

    logging.info("=" * 74)
    logging.info("%-11s %9s %9s %8s %8s %6s %9s", "region", "variants",
                 "judged", "AB fail", "all-het", "both", "% removed")
    for _, r in df.iterrows():
        logging.info("%-11s %9d %9d %8d %8d %6d %8.2f%%", r.region,
                     r.variant_sites, r.judged, r.drop_ab, r.drop_allhet,
                     r.drop_both, r.pct_dropped)
    for _, r in df.iterrows():
        logging.info("  %-10s removed %d of %d variant sites (%.2f%%): balance "
                     "%.2f%%, all heterozygous %.2f%%; %d (%.1f%%) had fewer "
                     "than %d usable heterozygotes, were not judged and stay "
                     "in the callset", r.region, r.dropped_total,
                     r.variant_sites, r.pct_dropped, r.pct_ab, r.pct_allhet,
                     r.not_judged, r.pct_not_judged, args.min_het)
    ret = df[df.region == "retained"].iloc[0]
    eli = df[df.region == "eliminated"].iloc[0]
    logging.info("")
    if ret.pct_dropped > 0:
        logging.info("Eliminated / retained removal rate = %.1fx",
                     eli.pct_dropped / ret.pct_dropped)
    else:
        logging.warning("No retained variant was removed, so the fold is "
                        "undefined. A zero here means the screen is not "
                        "firing, not that the callset is clean.")
    if eli.pct_dropped > 0:
        logging.info("One eliminated variant in %.0f was removed by the screen.",
                     100 / eli.pct_dropped)

    # --- agreement with the applied screen ---------------------------------
    n_disagree = 0
    xc = res["crosscheck"]
    if xc is not None:
        logging.info("")
        logging.info("Check against %s:", excl_path)
        for code in RCODE.values():
            x = xc[code]
            logging.info("  %-10s listed %d; agree %d removed + %d kept; "
                         "removed here only %d; listed only %d",
                         CODE_NAME[code], x["listed"], x["agree_drop"],
                         x["agree_keep"], x["panel_only"], x["list_only"])
            n_disagree += x["panel_only"] + x["list_only"]
        if n_disagree == 0:
            logging.info("  AGREE: every verdict here matches the applied "
                         "screen, so this table is what the callset lost.")
        else:
            logging.error("  DISAGREE at %d site(s): this panel does not "
                          "reproduce the screen that was applied. Check that "
                          "--ab-min/--ab-max/--min-het/--max-het-frac match the "
                          "filter_genotypes.py run, and that --vcf and "
                          "excluded_sites.txt come from the same directory.",
                          n_disagree)
            for kind, lab in (("panel_only", "removed here, not listed"),
                              ("list_only", "listed, kept here")):
                for e in res["examples"][kind]:
                    logging.error("    %s: %s", lab, e)

    fig, _ax, draw_info = draw(df, args.width, args.height_mm,
                               args.panel_letter or None)
    problems = audit(fig)
    paths = save(fig, outdir, "S16D_allele_balance")
    plt.close(fig)

    (Path(outdir) / "provenance.json").write_text(json.dumps({
        "argv": sys.argv,
        "vcf": str(vcf),
        "excluded_sites": str(excl_path) if excl_path else None,
        "samples": samples,
        "ab_window": [args.ab_min, args.ab_max],
        "min_het": args.min_het,
        "max_het_frac": args.max_het_frac,
        "counts": {CODE_NAME[k]: {kk: int(vv) for kk, vv in v.items()}
                   for k, v in res["counts"].items()},
        "crosscheck": ({CODE_NAME[k]: {kk: int(vv) for kk, vv in v.items()}
                        for k, v in xc.items()} if xc is not None else None),
        "crosscheck_examples": res["examples"],
        "panel": draw_info,
        "audit_problems": problems,
        "outputs": paths,
    }, indent=2))

    logging.info("Wrote %s", ", ".join(Path(p).name for p in paths))
    logging.info("Export from the SVG — Arial is absent on Isaac nodes.")
    logging.info("=" * 74)
    if n_disagree:
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
