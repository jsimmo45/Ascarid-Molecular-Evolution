#!/usr/bin/env python3
"""
reference_bias_check.py — quantify reference bias and test whether it can
account for the retained/eliminated piN/piS difference.

The concern
-----------
Reads carrying a non-reference allele map slightly worse than reads carrying
the reference allele, so heterozygotes lose alternate reads and are sometimes
miscalled homozygous-reference. Diversity is therefore underestimated, and the
underestimate grows with local divergence from the reference. Since these are
*A. lumbricoides* worms mapped to the *A. suum* v50 assembly, and eliminated
regions might be more diverged, this is a candidate explanation for the 3x
lower piS in eliminated genes.

(Easton et al. 2020 argue the two are one interbreeding complex rather than two
species, so this is intraspecific bias, not cross-species. Test 2 measures the
divergence directly rather than assuming either way.)

What it can and cannot do
-------------------------
Reference bias acts per READ, on total mismatch load. It has no way to know
whether the site a read covers is synonymous or nonsynonymous. So it can
depress piN and piS together — plausibly explaining the low piS — but has no
mechanism to change their RATIO. Test 3 checks that reasoning against the data
rather than trusting it.

Three tests
-----------
1. **Reference-allele fraction at heterozygous calls.** A true heterozygote
   gives 0.5. Systematic reference bias pushes it above 0.5. Compared between
   region classes. Run this on the PRE-allele-balance callset — the [0.25,0.75]
   filter trims exactly the tail this test is looking for.

2. **Divergence from the reference**, as the per-sample rate of
   homozygous-alternate genotypes by class. This is the quantity that sets how
   much bias is possible at all: if eliminated regions are no more diverged from
   the reference than retained ones, there is no room for a class-asymmetric
   effect.

3. **piN/piS stratified by local divergence.** For every site, count the variant
   sites within +/- ``--window`` bp as a proxy for the mismatch load a read
   spanning it carries. If reference bias drives the class difference, piN/piS
   should rise with local divergence. If piN/piS is flat across strata while
   piS falls, bias affects the denominator only — which is the same logic that
   defused the reduced-Ne alternative, applied to a different confound.

Example
-------
    python reference_bias_check.py \\
        --vcf vcf/filtered_miss4/step3_gapfilter.bcf \\
        --genome $DB_DIR/ascaris/v50/AG.v50.fa \\
        --gff ../As_miniprot.gff \\
        --gene-callable targets/gene_callable.tsv \\
        --outdir results/refbias

    # figure only, from an existing run (no VCF, seconds):
    python reference_bias_check.py --outdir results/refbias \\
        --replot-from results/refbias/refbias_20260822_filtered_miss4 \\
        --tag filtered_miss4_replot

Figure (supplementary fig. S16; file stem S17_reference_bias kept)
-----------------------------------------------------------------
Fixes of 2026-10-12, figure only (every table is unchanged): region colours
from paper_colors.REGION (the old ones came from pi_by_region); A as points
with intervals against 0.5 instead of bars on an axis starting at 0.48; plain
πS and πN/πS axis titles and "±" (no mathtext); strata ticks 2–3, 4–7, ≥8;
the eliminated / retained fold of each stratum in a row along the top of D,
with "×"; the class key above D, outside the axes; provenance.json, and
legend notes in the log (including sites per stratum).

Fix of 2026-10-13: pi_by_region is imported inside matplotlib.rc_context(),
so the rcParams it sets at import (Arial only, 15-20 pt) no longer override
the house style on the full-run path. The --replot-from path never imported
it and was not affected. A stratum under --min-s-sites in one class is now
named in the legend notes (it is drawn for the other class only, with no
fold), and the full run logs every stratum under the floor with its sites.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
# pi_by_region (codon index, FASTA and CDS readers) is imported inside main(),
# on the path that reads the VCF, so --replot-from needs only the plotting
# stack.

# style_boot locates lib/paper_colors.py and warns if it cannot.
from style_boot import paper_colors

if paper_colors is None:                # pragma: no cover
    raise SystemExit(
        "reference_bias_check.py builds supplementary fig. S16 and needs the "
        "full house style. Set PAPER_COLORS_DIR to the directory holding "
        "paper_colors.py.")

# Region colours from paper_colors.REGION itself. The 2026-08-22 figure took
# them from pi_by_region, which carried its own #2C6FBB / #C0392B, not the
# house #4C72B0 / #C44E52.
RETAINED = paper_colors.REGION["retained"]
ELIMINATED = paper_colors.REGION["eliminated"]
INK = paper_colors.NEUTRAL["ink"]

apply_style = paper_colors.apply_style
audit = paper_colors.audit
figure = paper_colors.figure
headroom = paper_colors.headroom
make_outdir = paper_colors.make_outdir
panel_letter = paper_colors.panel_letter
pc_save = paper_colors.save
FONT = paper_colors.FONT
LINEWIDTH = paper_colors.LINEWIDTH
RULE = paper_colors.NEUTRAL["rule"]

#: Supplementary fig. S16 (file stem S17_reference_bias, kept from an earlier
#: numbering): four panels, two rows, double-column width.
FIG_W = "double"
FIG_H_MM = 124.0


def query(vcf: Path, fmt: str):
    """Stream `bcftools query` line by line."""
    cmd = f"bcftools query -f '{fmt}' {vcf}"
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 22)
    for line in proc.stdout:
        yield line.rstrip("\n").split("\t")
    proc.wait()
    if proc.returncode != 0:
        logging.error("bcftools query failed:\n%s", proc.stderr.read()[-2000:])
        sys.exit(proc.returncode)


# ---------------------------------------------------------------------------
# Figure (supplementary fig. S16; file stem S17_reference_bias kept)
# ---------------------------------------------------------------------------

def stratum_display(label: str) -> str:
    """Tick text for a stratum label as stored in the table: "2-3" -> "2–3",
    ">=8" -> "≥8". The table keeps the ASCII form."""
    s = str(label)
    if s.startswith(">="):
        return "≥" + s[2:]
    return s.replace("-", "–")


def text_problems(fig, tol_px: float = 1.0) -> list[str]:
    """Drawn text overlapping other drawn text or running past the canvas
    edge. audit() checks neither, and both save without complaint."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    texts = list(fig.texts)
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks():
                loc = tick.get_loc()
                if loc is not None and lo - 1e-9 <= loc <= hi + 1e-9:
                    texts += [t for t in (tick.label1, tick.label2)
                              if t.get_visible()]
            texts.append(axis.label)
        texts += list(ax.texts)
        leg = ax.get_legend()
        if leg is not None:
            texts += list(leg.get_texts())
    texts = [t for t in texts if t.get_visible() and (t.get_text() or "").strip()]
    boxes = [t.get_window_extent(rend) for t in texts]
    out = []
    for t, b in zip(texts, boxes):
        if b.x0 < -tol_px or b.y0 < -tol_px or b.x1 > W + tol_px \
                or b.y1 > H + tol_px:
            out.append(f"text {t.get_text()[:24]!r} runs past the canvas edge")
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = boxes[i], boxes[j]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) > tol_px
                    and min(a.y1, b.y1) - max(a.y0, b.y0) > tol_px):
                out.append(f"text {texts[i].get_text()[:16]!r} overlaps "
                           f"{texts[j].get_text()[:16]!r}")
    return out


def plot_refbias(t1: list[dict], t2: list[dict], strat: pd.DataFrame,
                 window: int, outdir: Path, stem: str):
    """
    Supplementary fig. S16 (S18 in the frozen map; stem S17_reference_bias).

    A  mean reference-allele fraction at heterozygous calls per class, as a
       point with its bootstrap 95% interval, against the 0.5 no-bias line.
       Points rather than bars: the axis does not start at 0 (the panel is
       about distance from 0.5), and a bar cut off at 0.48 exaggerates the
       difference.
    B  homozygous-alternate genotype rate per class (bars from 0).
    C  piS and D piN/piS by local divergence (variants within +/- window bp),
       one line per class; the eliminated / retained piN/piS of each stratum
       in a row along the top of D, clear of the points; the class key above
       D. Colours are paper_colors REGION.

    Returns (files, problems, facts for the log and provenance.json).
    """
    t1d = {r["region"]: r for r in t1}
    t2d = {r["region"]: r for r in t2}
    fig, axes = figure(FIG_W, height_mm=FIG_H_MM, nrows=2, ncols=2)
    (axA, axB), (axC, axD) = axes
    order = [r for r in ("retained", "eliminated") if r in t1d or r in t2d]
    cols = {"retained": RETAINED, "eliminated": ELIMINATED}
    facts: dict = {}

    # A -- reference-allele fraction at heterozygous calls
    if t1d:
        oa = [r for r in order if r in t1d]
        x = np.arange(len(oa), dtype=float)
        vals = np.array([float(t1d[r]["mean_ref_frac"]) for r in oa])
        lo = np.array([float(t1d[r]["lo"]) for r in oa])
        hi = np.array([float(t1d[r]["hi"]) for r in oa])
        axA.axhline(0.5, ls="--", lw=LINEWIDTH["rule"], color=RULE, zorder=1)
        axA.errorbar(x, vals, yerr=[vals - lo, hi - vals], fmt="none",
                     ecolor=INK, elinewidth=LINEWIDTH["thin"], capsize=2.5,
                     zorder=3)
        for xi, v, r in zip(x, vals, oa):
            axA.plot([xi], [v], marker="o", ms=5.0, ls="none",
                     markerfacecolor=cols[r], markeredgecolor=INK,
                     markeredgewidth=0.6, zorder=4)
        span = max(float(hi.max()), 0.5) - min(float(lo.min()), 0.5)
        axA.set_ylim(min(float(lo.min()), 0.5) - 0.25 * span,
                     max(float(hi.max()), 0.5) + 0.25 * span)
        axA.set_xlim(-0.6, len(oa) - 0.4)
        axA.set_xticks(x)
        axA.set_xticklabels([f"{r.capitalize()}\n(n = {int(t1d[r]['n']):,})"
                             for r in oa])
        axA.set_ylabel("Reference-allele fraction\nat heterozygous calls")
        facts["A"] = {r: [round(float(t1d[r]["mean_ref_frac"]), 4),
                          round(float(t1d[r]["lo"]), 4),
                          round(float(t1d[r]["hi"]), 4)] for r in oa}

    # B -- divergence from the reference
    if t2d:
        ob = [r for r in order if r in t2d]
        x = np.arange(len(ob), dtype=float)
        vals = [float(t2d[r]["hom_alt_rate"]) for r in ob]
        axB.bar(x, vals, 0.6, color=[cols[r] for r in ob],
                edgecolor="black", linewidth=LINEWIDTH["thin"])
        for xi, v in zip(x, vals):
            axB.annotate(f"{v:.4f}", (xi, v), textcoords="offset points",
                         xytext=(0, 2), ha="center",
                         fontsize=FONT["annotation"])
        axB.set_xticks(x)
        axB.set_xticklabels([f"{r.capitalize()}\n"
                             f"(n = {int(t2d[r]['called']):,})" for r in ob])
        axB.set_ylabel("Homozygous-alternate rate")
        axB.set_ylim(0, max(vals) * 1.25)

    # C, D -- the stratified test
    folds = {}
    if not strat.empty:
        labs = list(dict.fromkeys(strat.stratum))
        xs = np.arange(len(labs), dtype=float)
        top_d = []
        for region in order:
            sdf = (strat[strat.region == region].set_index("stratum")
                   .reindex(labs))
            if not np.isfinite(sdf.piS.to_numpy(float)).any():
                # No stratum in this class met the site floor. Skip it rather
                # than drawing an invisible series with a legend key.
                logging.warning("S16 C/D: %s has no stratum above the site "
                                "floor and is omitted from the panel", region)
                continue
            axC.plot(xs, sdf.piS, "o-", color=cols[region],
                     lw=LINEWIDTH["data"], ms=3.0, label=region.capitalize())
            axD.plot(xs, sdf.piN_piS, "o-", color=cols[region],
                     lw=LINEWIDTH["data"], ms=3.0, label=region.capitalize())
            top_d += [v for v in sdf.piN_piS.to_numpy(float) if v == v]
        # The fold per stratum is the quantity the panel argues about, so it is
        # printed -- in one row along the top of D, in ink, so no label can sit
        # on a point (the 4-7 label did, beside the eliminated point).
        rr = strat[strat.region == "retained"].set_index("stratum").reindex(labs)
        ee = strat[strat.region == "eliminated"].set_index("stratum").reindex(labs)
        for i, lab in enumerate(labs):
            f = ee.piN_piS.get(lab, np.nan) / rr.piN_piS.get(lab, np.nan)
            if f == f:
                folds[str(lab)] = round(float(f), 3)
                axD.annotate(f"{f:.1f}×", (i, 0.955),
                             xycoords=axD.get_xaxis_transform(), ha="center",
                             va="top", fontsize=FONT["annotation"], color=INK)
        for a, ylab in ((axC, "πS"), (axD, "πN/πS")):
            a.set_xticks(xs)
            a.set_xticklabels([stratum_display(s) for s in labs])
            a.set_xlim(-0.5, len(labs) - 0.5)
            a.set_xlabel(f"Variants within ±{window} bp")
            a.set_ylabel(ylab)
            a.set_ylim(bottom=0)
        if top_d:
            # data in the lower 82% of D, the fold row above it
            axD.set_ylim(0, max(top_d) / 0.82)
        if axD.get_legend_handles_labels()[0]:
            axD.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2,
                       frameon=False, fontsize=FONT["legend"],
                       handlelength=1.4, handletextpad=0.5, columnspacing=1.2,
                       borderaxespad=0.0)
        facts["D_fold_by_stratum"] = folds
        # A stratum below the site floor in one class is drawn for the other
        # class only, with no fold above it; the legend has to say which.
        present = {lab: [rg for rg, df_ in (("retained", rr), ("eliminated", ee))
                         if np.isfinite(float(df_.piN_piS.get(lab, np.nan)))]
                   for lab in labs}
        facts["CD_one_class_only"] = {str(lab): rgs[0] for lab, rgs in
                                      present.items() if len(rgs) == 1}
        facts["CD_sites_by_stratum"] = {
            f"{r.region} {r.stratum}": [int(round(float(r.S_sites))),
                                        int(round(float(r.N_sites)))]
            for r in strat.itertuples()}

    for a, letter in ((axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")):
        panel_letter(a, letter)
    problems = list(audit(fig)) + text_problems(fig)
    if problems:
        logging.warning("S16 breaks %d rule(s): %s", len(problems),
                        "; ".join(problems))
    written = pc_save(fig, str(outdir), stem)
    plt.close(fig)
    logging.info("Figure: %s", ", ".join(Path(w).name for w in written))
    return written, problems, facts


def legend_notes(t1: list[dict], facts: dict, window: int,
                 min_s_sites: float | None = None) -> None:
    """What every mark of S16 is, and the counts the legend needs, to the log."""
    logging.info("")
    logging.info("LEGEND NOTES (supplementary fig. S16)")
    logging.info("  A: mean reference-allele fraction at heterozygous calls "
                 "(depth >= 6), bootstrap 95% interval; dashed line = 0.5, "
                 "no bias; n = heterozygous calls.")
    if len(t1) == 2:
        r_, e_ = float(t1[0]["mean_ref_frac"]), float(t1[1]["mean_ref_frac"])
        word = ("more" if e_ - r_ > 0.0005 else
                "less" if r_ - e_ > 0.0005 else "the same")
        logging.info("     retained %.3f, eliminated %.3f: %s reference allele "
                     "at eliminated heterozygous calls (Results cites only B "
                     "and D, so the legend must explain A).", r_, e_, word)
    logging.info("  B: homozygous-alternate genotypes / called genotypes; "
                 "n = called genotypes.")
    logging.info("  C, D: piS and piN/piS by the number of variant sites within "
                 "+/-%d bp of each site; numbers along the top of D = "
                 "eliminated / retained piN/piS in that stratum.", window)
    if facts.get("CD_sites_by_stratum"):
        logging.info("     synonymous / nonsynonymous sites per stratum: %s",
                     "; ".join(f"{k} {v[0]:,} / {v[1]:,}"
                               for k, v in facts["CD_sites_by_stratum"].items()))
    one = facts.get("CD_one_class_only") or {}
    if one:
        floor = (f"{min_s_sites:,.0f}" if min_s_sites is not None
                 else "--min-s-sites")
        logging.info("     drawn for one class only, no fold printed (the other "
                     "class has fewer than %s synonymous sites there): %s",
                     floor, "; ".join(
                         f"{stratum_display(k)} {v} only" for k, v in one.items()))


def replot(args) -> int:
    """Redraw S16 from an earlier run's two tables, without the VCF."""
    src = args.replot_from.expanduser()
    summ, strat_p = src / "refbias_summary.tsv", src / "pi_by_local_divergence.tsv"
    for f in (summ, strat_p):
        if not f.is_file():
            logging.error("--replot-from: %s not found", f)
            return 2
    summary = pd.read_csv(summ, sep="\t")
    t1 = summary[summary.test == "allele_balance"].to_dict("records")
    t2 = summary[summary.test == "divergence"].to_dict("records")
    strat = pd.read_csv(strat_p, sep="\t")
    old = src / "provenance.json"
    window = args.window
    min_s = None
    if old.is_file():
        try:
            prov = json.loads(old.read_text())
            window = int(prov.get("window", window))
            min_s = prov.get("min_s_sites")
        except (OSError, ValueError, TypeError):
            pass
    else:
        logging.warning("No provenance.json in %s, so the run's --window is not "
                        "recorded; the x label uses --window %d. Pass the value "
                        "the run used if it was not the default.", src, window)
    outdir = Path(make_outdir("refbias",
                              base=str(args.outdir.expanduser().resolve()),
                              tag=args.tag))
    logging.info("Output: %s (replot from %s; no VCF is read)", outdir, src)
    for f in (summ, strat_p):
        shutil.copyfile(f, outdir / f.name)
    written, problems, facts = plot_refbias(t1, t2, strat, window, outdir,
                                            "S17_reference_bias")
    legend_notes(t1, facts, window, min_s)
    (outdir / "provenance.json").write_text(json.dumps({
        "argv": sys.argv, "replot_from": str(src.resolve()),
        "min_s_sites": min_s,
        "sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                   for f in (summ, strat_p)},
        "window": window, "figure": [Path(w).name for w in written],
        "figure_facts": facts, "audit": problems}, indent=2, default=str))
    return 0


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Quantify reference bias and test it against the piN/piS result.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf", type=Path,
                   help="Prefer the PRE-allele-balance callset "
                        "(step3_gapfilter.bcf): the AB filter trims the tail "
                        "test 1 is looking for. Required unless --replot-from.")
    p.add_argument("--genome", type=Path)
    p.add_argument("--gff", type=Path)
    p.add_argument("--gene-callable", type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--replot-from", type=Path, default=None,
                   help="Redraw the figure from an earlier run's folder "
                        "(refbias_summary.tsv and pi_by_local_divergence.tsv) "
                        "without reading the VCF. Both tables are copied into "
                        "the new output folder.")
    p.add_argument("--window", type=int, default=150,
                   help="Half-window in bp for the local-divergence proxy; "
                        "match the read length")
    p.add_argument("--strata", default="0,1,2,4,8,1000",
                   help="Bin edges for the local variant count")
    p.add_argument("--min-s-sites", type=float, default=1000)
    p.add_argument("--tag", default=None,
                   help="Suffix for the dated output directory")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = p.parse_args(argv)
    if args.replot_from is None:
        missing = [f for f in ("vcf", "genome", "gff", "gene_callable")
                   if getattr(args, f) is None]
        if missing:
            p.error("required unless --replot-from: " + ", ".join(
                "--" + m.replace("_", "-") for m in missing))
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    for noisy in ("fontTools", "matplotlib.font_manager", "matplotlib"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    apply_style()
    if args.replot_from is not None:
        return replot(args)
    outdir = Path(make_outdir('refbias', base=str(args.outdir.expanduser().resolve()), tag=args.tag))
    logging.info('Output: %s', outdir)

    # pi_by_region sets its own rcParams when it is imported (Arial only,
    # 15-20 pt labels). Imported after apply_style(), it overrode the house
    # style on this path only: on Isaac, where Arial is absent, every text
    # fell back to DejaVu Sans at the wrong size and the S16 labels overran
    # the canvas (2026-10-13). rc_context() restores the rcParams in force
    # on entry, so whatever the module sets at import is discarded.
    with matplotlib.rc_context():
        from pi_by_region import (CodonIndex, read_fasta, read_cds,
                                  syn_fraction, is_synonymous, COMPLEMENT)
    gc = pd.read_csv(args.gene_callable.expanduser(), sep="\t")
    gc = gc[gc["pass"].astype(bool) & gc.region.isin(["retained", "eliminated"])]
    region_of = dict(zip(gc.gene, gc.region))

    fasta = read_fasta(args.genome.expanduser())
    genes = read_cds(args.gff.expanduser(), set(gc.gene))
    index = CodonIndex()
    index.build(genes, fasta)
    del fasta
    strands = [genes[g]["strand"] for g in index.gene_names]
    gene_region = np.array([region_of.get(g, "") for g in index.gene_names])

    # =====================================================================
    # PASS 1 — variant positions (for the local-divergence proxy) and
    #          tests 1 and 2
    # =====================================================================
    logging.info("Pass 1: variant sites, allele balance, divergence")
    var_pos: dict[str, list[int]] = {}
    ab_rows = []
    hom_alt = {"retained": 0, "eliminated": 0}
    called = {"retained": 0, "eliminated": 0}
    het_n = {"retained": 0, "eliminated": 0}

    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT,%AD]\\n"
    for f in query(args.vcf.expanduser(), fmt):
        if len(f) < 5:
            continue
        chrom, pos0, alt = f[0], int(f[1]) - 1, f[3]
        if alt in (".", "", "<*>"):
            continue
        hit = index.lookup(chrom, pos0)
        if hit is None:
            continue
        region = gene_region[hit[0]]
        if region not in ("retained", "eliminated"):
            continue
        var_pos.setdefault(chrom, []).append(pos0)

        for cell in f[4:]:
            if "," not in cell:
                continue
            gt, _, ad = cell.partition(",")
            g = gt.replace("|", "/")
            if g in ("./.", ".", ""):
                continue
            called[region] += 1
            a = g.split("/")
            if len(a) != 2:
                continue
            if a[0] == a[1] and a[0] != "0":
                hom_alt[region] += 1
            elif a[0] != a[1]:
                counts = [int(x) for x in ad.split(",") if x not in (".", "")]
                try:
                    i, j = int(a[0]), int(a[1])
                    d_ref, d_alt = counts[i], counts[j]
                except (ValueError, IndexError):
                    continue
                tot = d_ref + d_alt
                if tot >= 6:
                    het_n[region] += 1
                    ab_rows.append((region, d_ref / tot, tot))

    ab = pd.DataFrame(ab_rows, columns=["region", "ref_frac", "depth"])
    ab.to_csv(outdir / "het_allele_balance.tsv.gz", sep="\t", index=False)

    logging.info("=" * 70)
    logging.info("TEST 1 — reference-allele fraction at heterozygous calls")
    logging.info("  (0.5 = unbiased; >0.5 = reference bias)")
    t1 = []
    for region in ("retained", "eliminated"):
        s = ab[ab.region == region].ref_frac
        if s.empty:
            continue
        # bootstrap the mean
        rng = np.random.default_rng(0)
        v = s.to_numpy()
        boot = [rng.choice(v, min(len(v), 200000), replace=True).mean()
                for _ in range(200)]
        lo, hi = np.percentile(boot, [2.5, 97.5])
        logging.info("  %-11s n=%8d  mean=%.4f [%.4f, %.4f]  median=%.4f",
                     region, len(s), s.mean(), lo, hi, s.median())
        t1.append({"test": "allele_balance", "region": region, "n": len(s),
                   "mean_ref_frac": s.mean(), "lo": lo, "hi": hi,
                   "median_ref_frac": s.median()})
    if len(t1) == 2:
        d = t1[1]["mean_ref_frac"] - t1[0]["mean_ref_frac"]
        logging.info("  eliminated - retained = %+.4f", d)
        if abs(d) < 0.01:
            logging.info("  READ: bias is the same in both classes. It cannot "
                         "produce a class difference in anything.")
        else:
            logging.info("  READ: bias differs between classes; quantify its "
                         "effect on piS with test 3 before interpreting piS.")

    logging.info("")
    logging.info("TEST 2 — divergence from the reference (hom-alt genotype rate)")
    t2 = []
    for region in ("retained", "eliminated"):
        n, h = called[region], hom_alt[region]
        rate = h / n if n else np.nan
        logging.info("  %-11s hom-alt %d of %d called genotypes = %.5f",
                     region, h, n, rate)
        t2.append({"test": "divergence", "region": region, "called": n,
                   "hom_alt": h, "hom_alt_rate": rate})
    if len(t2) == 2 and t2[0]["hom_alt_rate"]:
        logging.info("  eliminated / retained = %.3f",
                     t2[1]["hom_alt_rate"] / t2[0]["hom_alt_rate"])
        logging.info("  This bounds how much class-asymmetric bias is even "
                     "possible: no excess divergence, no room for excess bias.")

    # =====================================================================
    # PASS 2 — piN/piS stratified by local divergence
    # =====================================================================
    logging.info("")
    logging.info("TEST 3 — piN/piS by local divergence (+/- %d bp)", args.window)
    var_arr = {c: np.array(sorted(v), dtype=np.int64) for c, v in var_pos.items()}
    edges = [float(x) for x in args.strata.split(",")]

    n_str = len(edges) - 1
    acc = {k: np.zeros((2, n_str)) for k in ("S", "N", "pi_s", "pi_n")}
    reg_idx = {"retained": 0, "eliminated": 1}

    fmt2 = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
    for f in query(args.vcf.expanduser(), fmt2):
        if len(f) < 5:
            continue
        chrom, pos0, ref, alt = f[0], int(f[1]) - 1, f[2], f[3]
        hit = index.lookup(chrom, pos0)
        if hit is None:
            continue
        gi, ci, cp = hit
        region = gene_region[gi]
        if region not in reg_idx:
            continue
        codon = index.gene_seq[gi][3 * ci: 3 * ci + 3]
        if len(codon) != 3:
            continue
        sf = syn_fraction(codon, cp)
        if sf != sf:
            continue

        # local divergence: variant sites within +/- window, excluding self
        arr = var_arr.get(chrom)
        if arr is None or arr.size == 0:
            local = 0
        else:
            lo = np.searchsorted(arr, pos0 - args.window, "left")
            hi = np.searchsorted(arr, pos0 + args.window, "right")
            local = int(hi - lo)
            # Exclude the site itself if it is one of the variants. Guard with
            # a bare index check: an earlier version compared against
            # arr[hi - 1], which wraps to the last element when hi == 0.
            i = int(np.searchsorted(arr, pos0, "left"))
            if i < arr.size and arr[i] == pos0:
                local -= 1

        k = int(np.digitize(local, edges[1:-1], right=False))
        k = min(max(k, 0), n_str - 1)
        r = reg_idx[region]
        acc["S"][r, k] += sf
        acc["N"][r, k] += 1.0 - sf

        if alt in (".", "", "<*>"):
            continue
        alts = [a for a in alt.split(",") if a not in (".", "<*>")]
        if len(alts) != 1 or len(alts[0]) != 1 or len(ref) != 1:
            continue
        n_called = n_alt = 0
        for gt in f[4:]:
            g = gt.replace("|", "/")
            if g in ("./.", ".", ""):
                continue
            for a in g.split("/"):
                if a == ".":
                    continue
                n_called += 1
                if a != "0":
                    n_alt += 1
        if n_called < 2:
            continue
        p = n_alt / n_called
        pi = (n_called / (n_called - 1)) * 2.0 * p * (1.0 - p)
        ab_base = alts[0].upper()
        if strands[gi] == "-":
            ab_base = ab_base.translate(COMPLEMENT)
        syn = is_synonymous(codon, cp, ab_base)
        if syn is None:
            continue
        if syn:
            acc["pi_s"][r, k] += pi
        else:
            acc["pi_n"][r, k] += pi

    rows = []
    labels = []
    for k in range(n_str):
        lo_e, hi_e = edges[k], edges[k + 1]
        labels.append(f"{lo_e:.0f}" if hi_e - lo_e == 1 else
                      f"{lo_e:.0f}-{hi_e - 1:.0f}" if hi_e < 1000 else f">={lo_e:.0f}")
    below: list[tuple[str, str, float]] = []    # strata under --min-s-sites
    for region, r in reg_idx.items():
        for k in range(n_str):
            S, N = acc["S"][r, k], acc["N"][r, k]
            if S < args.min_s_sites:
                below.append((region, labels[k], float(S)))
                continue
            piS = acc["pi_s"][r, k] / S
            piN = acc["pi_n"][r, k] / N if N else np.nan
            rows.append({"region": region, "stratum": labels[k],
                         "S_sites": S, "N_sites": N, "piS": piS, "piN": piN,
                         "piN_piS": piN / piS if piS else np.nan})
            logging.info("  %-11s local=%-8s S=%9.0f  piS=%.5f  piN/piS=%.4f",
                         region, labels[k], S, piS, piN / piS if piS else np.nan)
    for region, lab, S in below:
        logging.info("  %-11s local=%-8s S=%9.0f  below --min-s-sites %.0f; "
                     "not estimated", region, lab, S, args.min_s_sites)
    strat = pd.DataFrame(rows)
    if strat.empty:
        # An empty frame has no columns either, so every `strat.region` below
        # raises AttributeError rather than simply finding nothing. Give it the
        # schema so the summary and the figure degrade to "no rows" instead.
        strat = pd.DataFrame(columns=["region", "stratum", "S_sites",
                                      "N_sites", "piS", "piN", "piN_piS"])
        logging.warning("No stratum met --min-s-sites; TEST 3 has no rows and "
                        "panels C and D of S17 will be empty.")
    strat.to_csv(outdir / "pi_by_local_divergence.tsv", sep="\t", index=False)

    logging.info("")
    for region in ("retained", "eliminated"):
        s = strat[strat.region == region]
        if len(s) >= 2:
            logging.info("  %-11s piS varies %.2fx across strata, piN/piS "
                         "varies %.2fx", region,
                         s.piS.max() / s.piS.min(),
                         s.piN_piS.max() / s.piN_piS.min())
    logging.info("")
    logging.info("  READ: if piS falls with local divergence while piN/piS "
                 "stays flat, reference bias moves the denominator only and "
                 "cannot explain the class difference in the ratio.")

    pd.DataFrame(t1 + t2).to_csv(outdir / "refbias_summary.tsv", sep="\t",
                                 index=False)

    # --- supplementary fig. S16 --------------------------------------------
    # Four panels: the two bias measurements that BOUND the concern (A, B) and
    # the stratified test that would detect it if the bound were wrong (C, D).
    written, problems, facts = plot_refbias(t1, t2, strat, args.window, outdir,
                                            "S17_reference_bias")
    legend_notes(t1, facts, args.window, args.min_s_sites)
    (outdir / "provenance.json").write_text(json.dumps({
        "argv": sys.argv, "vcf": str(args.vcf), "window": args.window,
        "strata": args.strata, "min_s_sites": args.min_s_sites,
        "strata_below_floor": [{"region": g, "stratum": s, "S_sites": v}
                               for g, s, v in below],
        "figure": [Path(w).name for w in written], "figure_facts": facts,
        "audit": problems}, indent=2, default=str))
    logging.info("Figures written to %s", outdir)

    logging.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
