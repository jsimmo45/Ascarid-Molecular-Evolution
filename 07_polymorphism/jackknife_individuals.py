#!/usr/bin/env python3
"""
jackknife_individuals.py — is the piN/piS class contrast driven by one worm?

Why this exists
---------------
check_duplicate_individuals.py found no duplicates, but it did find that IBS0
is not a pair property in this cohort. Two libraries — SRR31692059 and
SRR31675198 in the primary set — show depressed IBS0 against EVERY other
library, which no pedigree can produce: kinship is pairwise, and a worm cannot
be a close relative of all nine others. A global IBS0 deficit is instead the
signature of a sample that rarely makes a confident homozygous-alternate call,
i.e. one carrying excess apparent heterozygosity.

That matters far more than the sample-independence question it came from.
Excess apparent heterozygosity from collapsed paralogy is already known to be
~12x worse in eliminated regions in this dataset, and paralog-derived "variants"
are fixed amino-acid differences between copies, so they land nonsynonymous far
more often than real polymorphism does. A single contaminated or paralogy-prone
worm could therefore inflate piN specifically in eliminated genes and
manufacture the entire result.

This script tests that directly: recompute the whole class contrast with each
worm held out in turn. If the ~2.2x fold survives dropping SRR31692059 and
SRR31675198, the result is not a one-worm artifact and the alternative is
closed. If it collapses, better to find out here.

Design
------
ONE pass over the callset, not one pass per held-out worm.

Every leave-one-out subset sees the same sites in the same order, so per site
the only thing that changes is the allele count. Tracking each sample's called
and alternate allele counts separately makes every subset's pi computable from
the full-set totals by subtraction, which turns eleven four-hour passes into
one. The bcftools stream and the codon index — the two expensive parts — are
paid for once.

Denominators (callable synonymous and nonsynonymous sites) are deliberately
NOT recomputed per subset. Callable-site status was fixed upstream by the
F_MISSING filter across the whole cohort, so it does not depend on which worm
is dropped, and holding it constant is what makes the subsets comparable.

Codon machinery, the reference-base self-check and the clustered-variant filter
are imported from pi_by_region.py rather than reimplemented, so the two scripts
cannot drift apart.

What this test CANNOT do
------------------------
Leave-one-out sees only an artifact carried by a single library. One shared by
two or three worms — a prep batch, a reference, an aligner version — survives
every single drop, because removing one leaves the rest still carrying it, and
the jackknife then looks flat and reassuring. This was verified against a
synthetic callset with a three-library confound planted in it: every
leave-one-out fold stayed within 3% of the full-cohort value while the
artifact was plainly there.

So a flat jackknife means "no ONE-LIBRARY artifact", never "no technical
artifact". --drop-set exists for that gap: name the suspected group and remove
it in one subset. het_by_region.py is how you find out which libraries belong
in that group.

Gene set (changed 2026-10-13)
-----------------------------
pi_by_region.py pools only genes with at least --min-sites (100 there)
callable coding sites. This script used to accept --min-sites and ignore it,
so its full-cohort fold pooled the smaller genes too: 1.887 on the primary
callset, against 1.884 from pi_by_region.py. The option is now applied: sites
are accumulated per gene during the pass, genes below the floor are dropped
afterwards, and the rest are pooled by class, as pi_by_region.py does. The
default is 0, which keeps every gene and reproduces the jackknife reported in
the paper; --min-sites 100 uses exactly the gene set of the main estimate.

Example
-------
    python jackknife_individuals.py \\
        --vcf vcf/filtered_miss4/snps_invariant.bcf \\
        --genome AG.v50.fa --gff As_miniprot.gff \\
        --gene-callable targets/gene_callable.tsv \\
        --max-local-variants 3 \\
        --drop-set SRR31692059,SRR31675198 \\
        --outdir results/jackknife
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pi_by_region import (COMPLEMENT, CodonIndex, collect_variant_positions,
                          index_strand_minus, is_synonymous,
                          local_variant_count, read_cds, read_fasta,
                          syn_fraction)

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
    from paper_colors import RETAINED, ELIMINATED       # type: ignore
except Exception:
    RETAINED, ELIMINATED = "#2C6FBB", "#C0392B"

plt.rcParams.update({
    "font.family": "Arial", "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelsize": 18, "xtick.labelsize": 12, "ytick.labelsize": 13,
    "legend.fontsize": 13, "figure.dpi": 150,
})

REGIONS = ("retained", "eliminated")


def vcf_samples(vcf: Path) -> list[str]:
    r = subprocess.run(["bash", "-c", f"bcftools query -l {vcf}"],
                       capture_output=True, text=True, check=True)
    return [s for s in r.stdout.split() if s]


def parse_gt(gt: str) -> tuple[int, int]:
    """(alleles called, alternate alleles) for one biallelic genotype string."""
    if not gt or gt in ("./.", ".", ".|."):
        return 0, 0
    called = alt = 0
    for a in gt.replace("|", "/").split("/"):
        if a == ".":
            continue
        called += 1
        if a != "0":
            alt += 1
    return called, alt


def pi_unbiased(n_called: np.ndarray, n_alt: np.ndarray) -> np.ndarray:
    """
    Nei's unbiased pi per site, vectorised over subsets.

    n/(n-1) * 2p(1-p). Subsets with fewer than two called alleles contribute
    zero rather than a divide-by-zero, matching pi_by_region.py, which skips
    such sites after having already counted them in the denominator.
    """
    out = np.zeros_like(n_called, dtype=np.float64)
    ok = n_called >= 2
    if not ok.any():
        return out
    n = n_called[ok].astype(np.float64)
    p = n_alt[ok] / n
    out[ok] = (n / (n - 1.0)) * 2.0 * p * (1.0 - p)
    return out


def stream_jackknife(vcf: Path, index: CodonIndex, strands: list[str],
                     region_of: dict, subsets: list[np.ndarray],
                     n_samples: int, var_pos=None, max_local: int = 0,
                     window: int = 150) -> dict:
    """
    One pass; per GENE accumulate site denominators (subset-invariant) and pi
    numerators for every subset. Pooling by class happens afterwards, in
    pool_class(), so genes can first be filtered on their callable sites.

    ``subsets`` is a list of boolean masks over samples, subset 0 being the
    full cohort. Everything else is a mask with one or more worms switched off.

    Returns per-gene arrays, in codon-index order: region, n_sites, S_sites,
    N_sites, S4_sites, N0_sites (length n_genes) and pi_syn, pi_non, pi_4,
    pi_0 (n_genes x n_subsets).
    """
    n_sub = len(subsets)
    # mask matrix (subsets x samples) as float, so a subset's allele totals are
    # one matrix-vector product rather than a Python loop over subsets
    M = np.array([s.astype(np.float64) for s in subsets])

    # Per-gene accumulators (pi_by_region.py's layout), so that --min-sites
    # can drop genes after the pass. Site counts are Python lists: one scalar
    # add per site is cheaper on a list than on a numpy array.
    n_genes = len(index.gene_names)
    n_sites = [0] * n_genes
    S_sites = [0.0] * n_genes
    N_sites = [0.0] * n_genes
    S4_sites = [0] * n_genes
    N0_sites = [0] * n_genes
    pi_syn = np.zeros((n_genes, n_sub))
    pi_non = np.zeros((n_genes, n_sub))
    pi_4 = np.zeros((n_genes, n_sub))
    pi_0 = np.zeros((n_genes, n_sub))

    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
    proc = subprocess.Popen(["bash", "-c", f"bcftools query -f '{fmt}' {vcf}"],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, bufsize=1 << 22)

    called = np.zeros(n_samples)
    alt_n = np.zeros(n_samples)
    n_lines = n_mapped = n_multi = n_skip = n_clustered = 0
    n_refcheck = n_refbad = 0

    for line in proc.stdout:
        n_lines += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 + n_samples:
            continue
        chrom, pos, ref, alt = f[0], int(f[1]) - 1, f[2], f[3]

        hit = index.lookup(chrom, pos)
        if hit is None:
            continue
        gi, ci, cp = hit
        region = region_of.get(index.gene_names[gi])
        if region not in REGIONS:
            continue
        n_mapped += 1

        codon = index.gene_seq[gi][3 * ci: 3 * ci + 3]
        if len(codon) != 3:
            n_skip += 1
            continue

        # Same REF/codon self-check as pi_by_region.py: it validates
        # coordinates, strand and reading frame at once against data this
        # script did not generate.
        exp = codon[cp]
        obs = ref.upper() if not index_strand_minus(strands, gi) \
            else ref.upper().translate(COMPLEMENT)
        if len(obs) == 1 and obs in "ACGT":
            n_refcheck += 1
            if obs != exp:
                n_refbad += 1

        if max_local > 0 and var_pos is not None:
            arr = var_pos.get(chrom)
            if arr is not None and arr.size and \
                    local_variant_count(arr, pos, window) > max_local:
                n_clustered += 1
                continue

        sf = syn_fraction(codon, cp)
        if sf != sf:
            n_skip += 1
            continue

        n_sites[gi] += 1
        S_sites[gi] += sf
        N_sites[gi] += 1.0 - sf
        fold4 = sf == 1.0
        fold0 = sf == 0.0
        if fold4:
            S4_sites[gi] += 1
        elif fold0:
            N0_sites[gi] += 1

        if alt in (".", "", "<*>"):
            continue
        alts = [x for x in alt.split(",") if x not in (".", "<*>")]
        if len(alts) != 1 or len(alts[0]) != 1 or len(ref) != 1:
            n_multi += 1
            continue

        for k, gt in enumerate(f[4:4 + n_samples]):
            called[k], alt_n[k] = parse_gt(gt)
        if called.sum() < 2:
            continue

        alt_base = alts[0].upper()
        if index_strand_minus(strands, gi):
            alt_base = alt_base.translate(COMPLEMENT)
        syn = is_synonymous(codon, cp, alt_base)
        if syn is None:
            continue

        # every subset's pi at this site, in one shot
        pis = pi_unbiased(M @ called, M @ alt_n)
        if syn:
            pi_syn[gi] += pis
            if fold4:
                pi_4[gi] += pis
        else:
            pi_non[gi] += pis
            if fold0:
                pi_0[gi] += pis

    proc.wait()
    if proc.returncode != 0:
        logging.error("bcftools query failed:\n%s", proc.stderr.read()[-2000:])
        sys.exit(proc.returncode)

    logging.info("Streamed %d records; %d mapped into classified coding sites",
                 n_lines, n_mapped)
    logging.info("  skipped: %d untranslatable codons, %d multiallelic",
                 n_skip, n_multi)
    if n_clustered:
        logging.info("  skipped: %d site(s) in variant clusters", n_clustered)
    if n_refcheck:
        bad = 100 * n_refbad / n_refcheck
        logging.info("  REF/codon agreement: %d of %d disagree (%.4f%%)",
                     n_refbad, n_refcheck, bad)
        if bad > 0.1:
            logging.error("REF bases disagree with the codon index at %.2f%% "
                          "of sites — a coordinate, strand or phase bug. Do "
                          "NOT use these results.", bad)
            sys.exit(1)
    return {"region": np.array([region_of.get(g, "") for g in index.gene_names],
                               dtype=object),
            "n_sites": np.array(n_sites, dtype=np.int64),
            "S_sites": np.array(S_sites), "N_sites": np.array(N_sites),
            "S4_sites": np.array(S4_sites, dtype=np.float64),
            "N0_sites": np.array(N0_sites, dtype=np.float64),
            "pi_syn": pi_syn, "pi_non": pi_non, "pi_4": pi_4, "pi_0": pi_0}


def pool_class(per_gene: dict, region: str, keep: np.ndarray) -> dict:
    """
    One region class pooled over the genes in ``keep``: the subset-invariant
    site denominators and, per subset, the pi numerators, in the layout
    ratios() reads. Pooling sums over genes, as pi_by_region.pooled() does.
    """
    m = (per_gene["region"] == region) & keep
    out = {k: float(per_gene[k][m].sum())
           for k in ("S_sites", "N_sites", "S4_sites", "N0_sites")}
    for k in ("pi_syn", "pi_non", "pi_4", "pi_0"):
        out[k] = per_gene[k][m].sum(axis=0)
    out["genes"] = int(m.sum())
    return out


def ratios(a: dict, k: int) -> dict:
    """piN/piS and pi0/pi4 for one region class under subset ``k``."""
    piS = a["pi_syn"][k] / a["S_sites"] if a["S_sites"] else np.nan
    piN = a["pi_non"][k] / a["N_sites"] if a["N_sites"] else np.nan
    pi4 = a["pi_4"][k] / a["S4_sites"] if a["S4_sites"] else np.nan
    pi0 = a["pi_0"][k] / a["N0_sites"] if a["N0_sites"] else np.nan
    return {"piS": piS, "piN": piN, "pi4": pi4, "pi0": pi0,
            "piN_piS": piN / piS if piS else np.nan,
            "pi0_pi4": pi0 / pi4 if pi4 else np.nan}


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Leave-one-worm-out jackknife of the piN/piS class "
                    "contrast, in a single pass.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf", required=True, type=Path,
                   help="snps_invariant.bcf — MUST contain invariant sites")
    p.add_argument("--genome", required=True, type=Path)
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--gene-callable", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--drop-set", action="append", default=[],
                   help="Comma-separated sample names to drop together, as an "
                        "extra subset beyond the single-worm ones. Repeatable. "
                        "Use it for the two libraries with a global IBS0 "
                        "deficit: if they share a cause, dropping them one at "
                        "a time can leave the other still carrying the "
                        "artifact.")
    p.add_argument("--max-local-variants", type=int, default=0,
                   help="Clustered-variant filter; match the pi_by_region run "
                        "being checked (3 for the pi_clustfilt results)")
    p.add_argument("--local-window", type=int, default=150)
    p.add_argument("--min-sites", type=int, default=0,
                   help="Leave out genes with fewer callable coding sites "
                        "(sites in the callset with a codon class), from every "
                        "subset, as pi_by_region.py does with its own "
                        "--min-sites (100 there). 0 keeps every gene: the "
                        "jackknife reported in the paper was run that way, so "
                        "its full-cohort fold (1.887) differs slightly from "
                        "pi_by_region.py's (1.884). 100 uses exactly the main "
                        "estimate's gene set. Before 2026-10-13 this option "
                        "was accepted and ignored.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    samples = vcf_samples(args.vcf.expanduser())
    n = len(samples)
    logging.info("Samples: %d", n)
    if n < 3:
        logging.error("Jackknife needs at least three samples")
        return 1

    gc = pd.read_csv(args.gene_callable.expanduser(), sep="\t")
    gc = gc[gc["pass"].astype(bool) & gc.region.isin(REGIONS)]
    keep = set(gc.gene)
    region_of = dict(zip(gc.gene, gc.region))
    logging.info("Candidate genes: %d (%s)", len(keep),
                 gc.region.value_counts().to_dict())

    fasta = read_fasta(args.genome.expanduser())
    genes = read_cds(args.gff.expanduser(), keep)
    index = CodonIndex()
    index.build(genes, fasta)
    del fasta
    strands = [genes[g]["strand"] for g in index.gene_names]

    # subset 0 = everyone; 1..n = drop one; then any explicit --drop-set
    labels = ["(none)"] + samples
    subsets = [np.ones(n, bool)]
    for k in range(n):
        m = np.ones(n, bool)
        m[k] = False
        subsets.append(m)
    for spec in args.drop_set:
        drop = [s.strip() for s in spec.split(",") if s.strip()]
        unknown = [s for s in drop if s not in samples]
        if unknown:
            logging.error("--drop-set names not in the VCF: %s", unknown)
            return 1
        m = np.ones(n, bool)
        for s in drop:
            m[samples.index(s)] = False
        if m.sum() < 2:
            logging.error("--drop-set '%s' leaves fewer than two samples", spec)
            return 1
        subsets.append(m)
        labels.append("+".join(drop))
    logging.info("Subsets: %d (full cohort, %d single drops, %d explicit)",
                 len(subsets), n, len(args.drop_set))

    var_pos = None
    if args.max_local_variants > 0:
        var_pos = collect_variant_positions(args.vcf.expanduser())
        logging.info("Clustered-variant filter: >%d variants within +/-%d bp",
                     args.max_local_variants, args.local_window)

    per_gene = stream_jackknife(args.vcf.expanduser(), index, strands,
                                region_of, subsets, n, var_pos=var_pos,
                                max_local=args.max_local_variants,
                                window=args.local_window)

    # Gene filter after the pass, on each gene's callable coding sites, then
    # pooling by class: the order pi_by_region.py uses for its main estimate.
    keep = per_gene["n_sites"] >= args.min_sites
    kept = {r: int(((per_gene["region"] == r) & keep).sum()) for r in REGIONS}
    dropped = {r: int(((per_gene["region"] == r) & ~keep).sum())
               for r in REGIONS}
    if args.min_sites > 0:
        logging.info("Genes with >= %d callable coding sites: %d (%s); left "
                     "out: %d (%s)", args.min_sites, sum(kept.values()), kept,
                     sum(dropped.values()), dropped)
    else:
        logging.info("--min-sites 0: every gene pooled, %d (%s). "
                     "pi_by_region.py's main estimate drops genes under 100 "
                     "callable coding sites; pass --min-sites 100 to match it.",
                     sum(kept.values()), kept)
    acc = {r: pool_class(per_gene, r, keep) for r in REGIONS}

    rows = []
    for k, lab in enumerate(labels):
        r = ratios(acc["retained"], k)
        e = ratios(acc["eliminated"], k)
        rows.append({
            "dropped": lab,
            "n_worms": int(subsets[k].sum()),
            "piS_ret": r["piS"], "piN_ret": r["piN"],
            "piS_eli": e["piS"], "piN_eli": e["piN"],
            "piN_piS_ret": r["piN_piS"], "piN_piS_eli": e["piN_piS"],
            "fold_piN_piS": e["piN_piS"] / r["piN_piS"]
            if r["piN_piS"] else np.nan,
            "pi0_pi4_ret": r["pi0_pi4"], "pi0_pi4_eli": e["pi0_pi4"],
            "fold_pi0_pi4": e["pi0_pi4"] / r["pi0_pi4"]
            if r["pi0_pi4"] else np.nan,
        })
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "jackknife_summary.tsv", sep="\t", index=False)

    full = df.iloc[0]
    loo = df.iloc[1:1 + n]
    logging.info("=" * 74)
    logging.info("Full cohort (n=%d): piN/piS fold = %.3f   pi0/pi4 fold = %.3f",
                 n, full.fold_piN_piS, full.fold_pi0_pi4)
    logging.info("")
    logging.info("%-16s %6s %11s %11s %11s", "dropped", "n", "piN/piS ret",
                 "piN/piS eli", "fold")
    for _, r in df.iloc[1:].iterrows():
        flag = ""
        if r.dropped != "(none)" and abs(r.fold_piN_piS - full.fold_piN_piS) \
                > 0.15 * full.fold_piN_piS:
            flag = "   <-- shifts the fold by >15%"
        logging.info("%-16s %6d %11.4f %11.4f %11.3f%s", r.dropped, r.n_worms,
                     r.piN_piS_ret, r.piN_piS_eli, r.fold_piN_piS, flag)

    lo, hi = loo.fold_piN_piS.min(), loo.fold_piN_piS.max()
    logging.info("")
    logging.info("Leave-one-out fold range: %.3f - %.3f (full %.3f)",
                 lo, hi, full.fold_piN_piS)

    # Tukey pseudo-values give a jackknife SE over individuals, which is a
    # different and more conservative uncertainty than the gene bootstrap in
    # pi_by_region.py: it asks how much the answer depends on WHICH WORMS were
    # sequenced rather than which genes were analysed. Reviewers asking about
    # n=10 are asking for this number.
    ps = n * full.fold_piN_piS - (n - 1) * loo.fold_piN_piS.to_numpy()
    jk_mean = float(ps.mean())
    jk_se = float(np.sqrt(ps.var(ddof=1) / n))
    logging.info("Jackknife over individuals: fold = %.3f +/- %.3f (SE), "
                 "95%% ~ [%.3f, %.3f]", jk_mean, jk_se,
                 jk_mean - 1.96 * jk_se, jk_mean + 1.96 * jk_se)

    logging.info("")
    if lo > 1.0:
        logging.info("VERDICT: every leave-one-out fold stays above 1.0. No "
                     "single worm produces the class contrast, so the result "
                     "is not an artifact of one contaminated or paralogy-prone "
                     "library. Quote the range above.")
        logging.info("")
        logging.info("LIMIT OF THIS TEST: leave-one-out can only see an "
                     "artifact carried by ONE library. An artifact shared by "
                     "two or three — same prep batch, same reference, same "
                     "aligner — survives every single drop, because dropping "
                     "one leaves the others still carrying it, and the fold "
                     "looks reassuringly flat. Do not read a flat jackknife as "
                     "'no technical artifact'; it means 'no ONE-LIBRARY "
                     "artifact'. Use --drop-set to remove a suspected group "
                     "together, and read het_by_region.py for whether such a "
                     "group exists.")
    else:
        logging.warning("VERDICT: dropping %s takes the fold to %.3f. The "
                        "contrast depends on that library — find out why "
                        "before going further.",
                        loo.loc[loo.fold_piN_piS.idxmin(), "dropped"], lo)
    if args.drop_set:
        for _, r in df.iloc[1 + n:].iterrows():
            logging.info("Combined drop %-24s fold = %.3f (n=%d)",
                         r.dropped, r.fold_piN_piS, r.n_worms)

    # --- figure -------------------------------------------------------------
    d = df.iloc[1:]
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(7.4, 0.42 * len(d) + 2.2))
    ax.axvline(full.fold_piN_piS, color="black", lw=1.6, ls="--",
               label=f"Full cohort ({full.fold_piN_piS:.2f})")
    ax.axvline(1.0, color="grey", lw=1.4, ls=":")
    ax.scatter(d.fold_piN_piS, y, s=90, color=ELIMINATED, edgecolor="black",
               linewidth=0.6, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(d.dropped)
    ax.set_xlabel("Eliminated / retained " + r"$\pi_N/\pi_S$")
    ax.set_ylabel("Worm held out")
    ax.set_xlim(left=min(0.95, float(d.fold_piN_piS.min()) - 0.1))
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(outdir / f"S_jackknife_individuals.{ext}",
                    bbox_inches="tight")
    plt.close(fig)
    logging.info("Figure: S_jackknife_individuals.{png,svg,pdf}")
    logging.info("Export from the SVG — Arial is absent on Isaac nodes.")
    logging.info("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
