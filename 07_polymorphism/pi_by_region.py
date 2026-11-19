#!/usr/bin/env python3
"""
pi_by_region.py — codon-aware piN/piS and the McDonald-Kreitman test, by region
class.

This is the analysis the whole pipeline exists to produce.

What it computes
----------------
For every callable coding site, the reference codon determines how much of that
site is "synonymous" and how much "nonsynonymous", by Nei-Gojobori counting: of
the three possible substitutions at that position, the fraction that leave the
amino acid unchanged contributes to S, the rest to N. Summing over the sites
that are **actually observed in the callset** gives per-gene S and N
denominators that reflect callability rather than the reference's gene models.
That is the entire reason for carrying invariant sites this far.

For every biallelic SNP, the same codon says whether the change is synonymous,
and per-site diversity is the unbiased estimator

    pi_i = n/(n-1) * 2 * p * (1-p)

with n the number of *called* alleles at that site and p the alternate
frequency. Unbiased for any n >= 2, which is why variable missingness costs
variance rather than accuracy.

    piN = sum(pi_i over nonsynonymous SNPs) / N_sites
    piS = sum(pi_i over synonymous SNPs)    / S_sites

The MK test then contrasts polymorphism against the divergence already measured
by the Ka/Ks pipeline. Ka and Ks are rates, so counts are recovered as
Dn = Ka * N and Ds = Ks * S using this script's own callable N and S. Reported
as the neutrality index NI = (Pn/Ps)/(Dn/Ds), alpha = 1 - NI, and DoS
(Stoletzki & Eyre-Walker) = Dn/(Dn+Ds) - Pn/(Pn+Ps). **DoS is the headline**:
NI is a ratio of ratios and is badly behaved when any cell is small, which it
will be for individual eliminated genes.

Why this matters for the paper's claim
--------------------------------------
Relaxed purifying selection predicts elevated piN/piS in eliminated genes AND
elevated dN/dS, with the two moving together, so alpha and DoS stay unchanged.
Adaptive evolution predicts elevated dN/dS *without* elevated piN/piS, raising
alpha. The dN/dS result alone cannot separate these. This can.

Caveats it enforces
-------------------
* Genes are pooled within class for the primary test. With 10 diploids, per-gene
  Pn and Ps are far too small to test individually; per-gene values are written
  out for inspection, not for significance.
* Codons containing an N, or not in the standard table, are skipped.
* Multiallelic sites are skipped (biallelic SNPs only).
* A position covered by more than one gene model is assigned to the first and
  counted once; the count of such positions is reported.
* Changes creating a stop codon are counted as nonsynonymous.
* Pn, Ps, P0 and P4 count only sites that segregate in the sample (both
  alleles among the called genotypes).

Changes 2026-10-14
------------------
The polymorphism counts (Pn, Ps, P0, P4, and so the "all variants" rows of
mk_summary.tsv and the Pn/Ps columns of pi_summary.tsv and per_gene.tsv)
included sites that do not segregate in the sample. Pn and Ps counted sites at
which every called allele is the alternate, which are fixed differences from
the reference; P0 and P4 counted those and also sites at which no called allele
is the alternate (the alternate survived only in masked genotypes). Both kinds
are now excluded and counted in the log.
Nothing else moves: pi is 0 at such sites, so piS, piN, piN/piS, pi0, pi4 and
pi0/pi4 are unchanged, and the MAF-floored MK rows already excluded them (their
minor-allele frequency is 0).

Example
-------
    python pi_by_region.py \\
        --vcf vcf/filtered_miss4/snps_invariant.bcf \\
        --genome $DB_DIR/ascaris/v50/AG.v50.fa \\
        --gff ../As_miniprot.gff \\
        --gene-callable targets/gene_callable.tsv \\
        --kaks ../kaks_regions/kaks_annotated.tsv \\
        --outdir results/pi_miss4
"""

from __future__ import annotations

import argparse
import gzip
import logging
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
    from paper_colors import RETAINED, ELIMINATED       # type: ignore
except Exception:
    RETAINED, ELIMINATED = "#2C6FBB", "#C0392B"

plt.rcParams.update({
    "font.family": "Arial", "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelsize": 20, "xtick.labelsize": 16, "ytick.labelsize": 16,
    "legend.fontsize": 15, "figure.dpi": 150,
})

# ---------------------------------------------------------------------------
# Genetic code
# ---------------------------------------------------------------------------

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODON_TABLE = dict(zip(CODONS, AAS))

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def syn_fraction(codon: str, pos: int) -> float:
    """
    Nei-Gojobori synonymous fraction at ``pos`` (0,1,2) of ``codon``.

    Of the three possible substitutions, the fraction leaving the amino acid
    unchanged. Changes to a stop codon count as nonsynonymous.
    """
    aa = CODON_TABLE.get(codon)
    if aa is None or aa == "*":
        return float("nan")
    same = 0
    for b in "ACGT":
        if b == codon[pos]:
            continue
        alt = codon[:pos] + b + codon[pos + 1:]
        alt_aa = CODON_TABLE.get(alt)
        if alt_aa is not None and alt_aa == aa:
            same += 1
    return same / 3.0


def is_synonymous(codon: str, pos: int, alt_base: str) -> bool | None:
    """True if substituting ``alt_base`` at ``pos`` is synonymous; None if the
    codon or the resulting codon is untranslatable or involves a stop."""
    aa = CODON_TABLE.get(codon)
    if aa is None or aa == "*":
        return None
    alt = codon[:pos] + alt_base + codon[pos + 1:]
    alt_aa = CODON_TABLE.get(alt)
    if alt_aa is None:
        return None
    if alt_aa == "*":
        return False
    return alt_aa == aa


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> dict[str, str]:
    """Whole-genome FASTA into {name: sequence}. ~276 Mb for AG.v50."""
    logging.info("Reading genome: %s", path)
    seqs: dict[str, list[str]] = {}
    name = None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    out = {k: "".join(v).upper() for k, v in seqs.items()}
    logging.info("  %d sequences, %.1f Mb", len(out),
                 sum(len(v) for v in out.values()) / 1e6)
    return out


def read_cds(gff: Path, keep: set[str]) -> dict[str, dict]:
    """CDS intervals per mRNA, 0-based half-open, restricted to ``keep``."""
    import re
    attr_id = re.compile(r"(?:^|;)ID=([^;]+)")
    attr_parent = re.compile(r"(?:^|;)Parent=([^;]+)")
    genes: dict[str, dict] = {}
    opener = gzip.open if str(gff).endswith(".gz") else open
    with opener(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                m = attr_id.search(f[8])
                if m and m.group(1) in keep:
                    genes.setdefault(m.group(1), {"chrom": f[0], "strand": f[6],
                                                  "cds": []})
            elif f[2] == "CDS":
                m = attr_parent.search(f[8])
                if m and m.group(1) in keep:
                    g = genes.setdefault(m.group(1), {"chrom": f[0],
                                                      "strand": f[6], "cds": []})
                    g["cds"].append((int(f[3]) - 1, int(f[4])))
    genes = {k: v for k, v in genes.items() if v["cds"]}
    for g in genes.values():
        ivs = sorted(g["cds"])
        merged = [ivs[0]]
        for s, e in ivs[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        g["cds"] = merged
    logging.info("GFF: %d of %d requested genes have CDS", len(genes), len(keep))
    return genes


class CodonIndex:
    """
    Position -> (gene, codon index, position in codon), built once per
    chromosome as sorted numpy arrays.

    A plain dict over ~14 M coding positions costs well over a gigabyte; four
    parallel arrays plus a binary search cost about 180 MB and are faster.
    """

    def __init__(self):
        self.pos: dict[str, np.ndarray] = {}
        self.gene: dict[str, np.ndarray] = {}
        self.codon: dict[str, np.ndarray] = {}
        self.cpos: dict[str, np.ndarray] = {}
        self.gene_names: list[str] = []
        self.gene_seq: list[str] = []
        self.n_overlap = 0

    def build(self, genes: dict[str, dict], fasta: dict[str, str]) -> None:
        by_chrom: dict[str, list] = defaultdict(list)
        seen: dict[str, set] = defaultdict(set)

        for gi, (gid, g) in enumerate(sorted(genes.items())):
            chrom, strand = g["chrom"], g["strand"]
            seq = fasta.get(chrom)
            if seq is None:
                continue
            positions: list[int] = []
            for s, e in g["cds"]:
                positions.extend(range(s, e))
            if strand == "-":
                positions.reverse()

            bases = [seq[p] for p in positions]
            if strand == "-":
                bases = [b.translate(COMPLEMENT) for b in bases]
            self.gene_names.append(gid)
            self.gene_seq.append("".join(bases))

            chrom_seen = seen[chrom]
            for i, p in enumerate(positions):
                if p in chrom_seen:      # position shared with another model
                    self.n_overlap += 1
                    continue
                chrom_seen.add(p)
                by_chrom[chrom].append((p, gi, i // 3, i % 3))

        for chrom, rows in by_chrom.items():
            rows.sort()
            arr = np.array(rows, dtype=np.int64)
            self.pos[chrom] = arr[:, 0].astype(np.int64)
            self.gene[chrom] = arr[:, 1].astype(np.int32)
            self.codon[chrom] = arr[:, 2].astype(np.int32)
            self.cpos[chrom] = arr[:, 3].astype(np.int8)

        n = sum(len(v) for v in self.pos.values())
        logging.info("Codon index: %d coding positions over %d sequences, "
                     "%d genes", n, len(self.pos), len(self.gene_names))
        if self.n_overlap:
            logging.info("  %d position(s) shared between gene models, counted "
                         "once", self.n_overlap)

    def lookup(self, chrom: str, pos0: int):
        arr = self.pos.get(chrom)
        if arr is None:
            return None
        i = np.searchsorted(arr, pos0)
        if i >= arr.size or arr[i] != pos0:
            return None
        return (int(self.gene[chrom][i]), int(self.codon[chrom][i]),
                int(self.cpos[chrom][i]))


# ---------------------------------------------------------------------------
# Streaming the callset
# ---------------------------------------------------------------------------

def collect_variant_positions(vcf: Path) -> dict[str, np.ndarray]:
    """
    Positions of every variant site, per chromosome, sorted.

    Used only for the clustered-variant filter. Costs one extra pass over the
    callset but only ~150k rows come back, so memory is trivial.
    """
    logging.info("Pre-pass: collecting variant positions for the "
                 "clustered-variant filter")
    cmd = f"bcftools query -f '%CHROM\\t%POS\\n' -i 'N_ALT>0' {vcf}"
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 22)
    pos: dict[str, list[int]] = {}
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) == 2:
            pos.setdefault(f[0], []).append(int(f[1]) - 1)
    proc.wait()
    if proc.returncode != 0:
        logging.error("bcftools query failed:\n%s", proc.stderr.read()[-2000:])
        sys.exit(proc.returncode)
    out = {c: np.array(sorted(v), dtype=np.int64) for c, v in pos.items()}
    logging.info("  %d variant site(s) over %d sequences",
                 sum(v.size for v in out.values()), len(out))
    return out


def local_variant_count(arr: np.ndarray, pos0: int, window: int) -> int:
    """Variant sites within +/- window of pos0, excluding pos0 itself."""
    lo = np.searchsorted(arr, pos0 - window, "left")
    hi = np.searchsorted(arr, pos0 + window, "right")
    n = int(hi - lo)
    i = int(np.searchsorted(arr, pos0, "left"))
    if i < arr.size and arr[i] == pos0:
        n -= 1
    return n


def stream_sites(vcf: Path, index: CodonIndex, strands: list[str],
                 mk_mafs: list[float] | None = None,
                 var_pos: dict | None = None,
                 max_local: int = 0, window: int = 150):
    """
    Walk the all-sites callset once, accumulating per gene:
      n_sites, S_sites, N_sites, pi_syn, pi_non, Ps, Pn, Ps_maf, Pn_maf

    ``mk_mafs`` gives one or more minor-allele-frequency floors for the MK
    counts only; each is tallied separately in a single pass.
    pi is always computed from every variant; the MK counts are additionally
    tallied with rare variants excluded, because slightly deleterious mutations
    segregate at low frequency, inflate Pn, and deflate alpha — and they do so
    MORE in a class under relaxed selection, which is exactly the class being
    tested. Without this the alpha comparison is partly circular.

    Genotypes are counted directly rather than trusting INFO/AC and AN, which
    are absent at invariant sites.
    """
    mk_mafs = mk_mafs or []
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
    cmd = f"bcftools query -f '{fmt}' {vcf}"
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 22)

    ng = len(index.gene_names)
    keys = ["n_sites", "S_sites", "N_sites", "pi_syn", "pi_non", "Ps", "Pn",
            "S4_sites", "N0_sites", "pi_4", "pi_0", "P4", "P0"]
    for i in range(len(mk_mafs)):
        keys += [f"Ps_maf{i}", f"Pn_maf{i}"]
    acc = {k: np.zeros(ng, dtype=np.float64) for k in keys}

    n_lines = n_mapped = n_skip_codon = n_multi = 0
    n_refcheck = n_refbad = 0
    n_clustered = 0
    n_foldbad = 0
    n_fixed_alt = n_alt_absent = 0
    for line in proc.stdout:
        n_lines += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        chrom, pos, ref, alt = f[0], int(f[1]) - 1, f[2], f[3]

        hit = index.lookup(chrom, pos)
        if hit is None:
            continue
        gi, ci, cp = hit
        n_mapped += 1

        seq = index.gene_seq[gi]
        codon = seq[3 * ci: 3 * ci + 3]
        if len(codon) != 3:
            n_skip_codon += 1
            continue
        # Self-check: the reference base reported by the caller must equal the
        # codon base the index predicts, after complementing on the minus
        # strand. Any coordinate, strand or phase error shows up here rather
        # than as a quietly wrong piN/piS.
        exp = codon[cp]
        obs = ref.upper() if not index_strand_minus(strands, gi) \
            else ref.upper().translate(COMPLEMENT)
        if len(obs) == 1 and obs in "ACGT":
            n_refcheck += 1
            if obs != exp:
                n_refbad += 1

        # Clustered-variant filter. A dense local pile-up of variants is the
        # classic misalignment / collapsed-paralog signature: such "variants"
        # are fixed differences between copies and are random with respect to
        # codon position, so they drag piN/piS toward the no-selection value.
        # Applied to EVERY site, variant and invariant, so the denominator
        # matches the numerator.
        if max_local > 0 and var_pos is not None:
            arr = var_pos.get(chrom)
            if arr is not None and arr.size:
                if local_variant_count(arr, pos, window) > max_local:
                    n_clustered += 1
                    continue

        sf = syn_fraction(codon, cp)
        if sf != sf:                      # NaN: stop codon or contains N
            n_skip_codon += 1
            continue

        acc["n_sites"][gi] += 1
        acc["S_sites"][gi] += sf
        acc["N_sites"][gi] += 1.0 - sf

        # Strictly 0-fold and 4-fold degenerate sites. At a 4-fold site every
        # substitution is synonymous; at a 0-fold site none is. pi0/pi4 is
        # therefore the same contrast as piN/piS but free of the Nei-Gojobori
        # weighting model — a useful check that the headline is not an artifact
        # of the counting scheme.
        fold4 = (sf == 1.0)
        fold0 = (sf == 0.0)
        if fold4:
            acc["S4_sites"][gi] += 1
        elif fold0:
            acc["N0_sites"][gi] += 1

        if alt in (".", "", "<*>"):
            continue                       # invariant: denominator only
        alts = [a for a in alt.split(",") if a not in (".", "<*>")]
        if len(alts) != 1 or len(alts[0]) != 1 or len(ref) != 1:
            n_multi += 1
            continue

        n_called = 0
        n_alt = 0
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
        # A polymorphism is a site that SEGREGATES in the sample: both alleles
        # among the called genotypes. Where every called allele is the
        # alternate, the site is a fixed difference from the reference (p = 1,
        # pi = 0), not a polymorphism, and must not enter Pn, Ps, P0 or P4.
        # Where no alternate allele survives genotype masking (p = 0), there
        # is nothing to count either.
        segregating = 0 < n_alt < n_called
        if n_alt == n_called:
            n_fixed_alt += 1
        elif n_alt == 0:
            n_alt_absent += 1

        alt_base = alts[0].upper()
        if index_strand_minus(strands, gi):
            alt_base = alt_base.translate(COMPLEMENT)
        syn = is_synonymous(codon, cp, alt_base)
        if syn is None:
            continue
        # Consistency check on the degeneracy classification: a variant at a
        # 4-fold site must be synonymous and one at a 0-fold site must not be.
        if fold4:
            if syn:
                acc["pi_4"][gi] += pi
                if segregating:
                    acc["P4"][gi] += 1
            else:
                n_foldbad += 1
        elif fold0:
            if not syn:
                acc["pi_0"][gi] += pi
                if segregating:
                    acc["P0"][gi] += 1
            else:
                n_foldbad += 1

        maf = min(p, 1.0 - p)
        pre = "Ps" if syn else "Pn"
        if syn:
            acc["pi_syn"][gi] += pi
        else:
            acc["pi_non"][gi] += pi
        if segregating:
            acc[pre][gi] += 1
            for i, t in enumerate(mk_mafs):
                if maf >= t:
                    acc[f"{pre}_maf{i}"][gi] += 1

    proc.wait()
    if proc.returncode != 0:
        logging.error("bcftools query failed:\n%s", proc.stderr.read()[-2000:])
        sys.exit(proc.returncode)

    logging.info("Streamed %d records; %d mapped to coding positions", n_lines,
                 n_mapped)
    logging.info("  skipped: %d untranslatable codons, %d multiallelic",
                 n_skip_codon, n_multi)
    logging.info("  not counted as polymorphisms: %d biallelic SNP(s) at which "
                 "every called allele is the alternate (fixed differences from "
                 "the reference) and %d at which none is; pi = 0 at both, so pi "
                 "is unaffected", n_fixed_alt, n_alt_absent)
    if n_foldbad:
        logging.error("%d variant(s) contradict their degeneracy class "
                      "(synonymous at a 0-fold site or vice versa). That is a "
                      "codon-table or frame bug.", n_foldbad)
        sys.exit(1)
    if n_clustered:
        logging.info("  skipped: %d site(s) in variant clusters "
                     "(>%d variants within +/-%d bp)", n_clustered, max_local,
                     window)

    # The reference-base check is the single most informative correctness
    # signal available: it validates coordinates, strand and reading frame
    # simultaneously, against data the script did not generate.
    if n_refcheck:
        bad = 100 * n_refbad / n_refcheck
        logging.info("  REF/codon agreement: %d of %d disagree (%.4f%%)",
                     n_refbad, n_refcheck, bad)
        if bad > 0.1:
            logging.error("REF bases disagree with the codon index at %.2f%% "
                          "of sites. That is a coordinate, strand or phase "
                          "bug, not noise — do NOT use these results.", bad)
            sys.exit(1)
    return acc


def index_strand_minus(strands: list[str], gi: int) -> bool:
    return strands[gi] == "-"


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def pooled(df: pd.DataFrame) -> dict:
    S, N = df.S_sites.sum(), df.N_sites.sum()
    piS = df.pi_syn.sum() / S if S else np.nan
    piN = df.pi_non.sum() / N if N else np.nan
    S4 = df.S4_sites.sum() if "S4_sites" in df else 0
    N0 = df.N0_sites.sum() if "N0_sites" in df else 0
    pi4 = df.pi_4.sum() / S4 if S4 else np.nan
    pi0 = df.pi_0.sum() / N0 if N0 else np.nan
    return {"genes": len(df), "S_sites": S, "N_sites": N,
            "Ps": df.Ps.sum(), "Pn": df.Pn.sum(),
            "piS": piS, "piN": piN,
            "piN_piS": piN / piS if piS else np.nan,
            "S4_sites": S4, "N0_sites": N0, "pi4": pi4, "pi0": pi0,
            "pi0_pi4": pi0 / pi4 if pi4 else np.nan,
            "Ds": df.Ds.sum() if "Ds" in df else np.nan,
            "Dn": df.Dn.sum() if "Dn" in df else np.nan}


def mk_stats(Pn, Ps, Dn, Ds) -> dict:
    out = {"Pn": Pn, "Ps": Ps, "Dn": Dn, "Ds": Ds}
    out["NI"] = ((Pn / Ps) / (Dn / Ds)) if Ps and Ds and Dn else np.nan
    out["alpha"] = 1 - out["NI"] if out["NI"] == out["NI"] else np.nan
    dos_d = (Dn / (Dn + Ds)) if (Dn + Ds) else np.nan
    dos_p = (Pn / (Pn + Ps)) if (Pn + Ps) else np.nan
    out["DoS"] = dos_d - dos_p
    return out


def bootstrap_ci(df: pd.DataFrame, fn, n_boot: int, seed: int = 0):
    rng = np.random.default_rng(seed)
    idx = np.arange(len(df))
    vals = []
    for _ in range(n_boot):
        take = rng.choice(idx, len(idx), replace=True)
        try:
            vals.append(fn(df.iloc[take]))
        except Exception:
            continue
    vals = np.array([v for v in vals if v == v])
    if vals.size == 0:
        return (np.nan, np.nan)
    return tuple(np.percentile(vals, [2.5, 97.5]))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Codon-aware piN/piS and MK test by region class.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf", required=True, type=Path,
                   help="snps_invariant.bcf — MUST contain invariant sites")
    p.add_argument("--genome", required=True, type=Path)
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--gene-callable", required=True, type=Path,
                   help="gene_callable.tsv from make_callable_targets.py")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--kaks", type=Path,
                   help="kaks_annotated.tsv; enables the MK test")
    p.add_argument("--min-sites", type=int, default=100,
                   help="Genes with fewer callable coding sites are excluded")
    p.add_argument("--max-local-variants", type=int, default=0,
                   help="Drop sites with more than this many variants within "
                        "+/- --local-window bp. 0 disables. Clustered variants "
                        "are the misalignment/paralogy signature; excluding "
                        "them STRENGTHENS the class contrast in this dataset "
                        "(2.07x -> 2.27x at a threshold of 3), so it is a "
                        "robustness check, not a rescue.")
    p.add_argument("--local-window", type=int, default=150,
                   help="Half-window for --max-local-variants; match read length")
    p.add_argument("--mk-maf", default="0.10,0.15,0.20",
                   help="Comma-separated minor-allele-frequency floors for the "
                        "MK counts only (pi always uses every variant). All are "
                        "tallied in one pass. Standard correction for slightly "
                        "deleterious polymorphism deflating alpha. With 20 "
                        "haplotypes, 0.10 only excludes singletons, so check "
                        "0.15 and 0.20 before concluding anything about a "
                        "difference in alpha between classes. Empty disables.")
    p.add_argument("--n-boot", type=int, default=2000)
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

    gc = pd.read_csv(args.gene_callable.expanduser(), sep="\t")
    gc = gc[gc["pass"].astype(bool) & gc.region.isin(["retained", "eliminated"])]
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

    mk_mafs = [float(x) for x in str(args.mk_maf).split(",") if x.strip()]
    logging.info("MK minor-allele-frequency floors: %s",
                 ", ".join(f"{t:g}" for t in mk_mafs) or "(none)")
    var_pos = None
    if args.max_local_variants > 0:
        var_pos = collect_variant_positions(args.vcf.expanduser())
        logging.info("Clustered-variant filter: drop sites with >%d variants "
                     "within +/-%d bp", args.max_local_variants,
                     args.local_window)
    acc = stream_sites(args.vcf.expanduser(), index, strands,
                       mk_mafs=mk_mafs, var_pos=var_pos,
                       max_local=args.max_local_variants,
                       window=args.local_window)

    df = pd.DataFrame({k: v for k, v in acc.items()})
    df.insert(0, "gene", index.gene_names)
    df["region"] = df.gene.map(region_of)
    df = df[df.n_sites >= args.min_sites].copy()
    logging.info("Genes with >= %d callable coding sites: %d (%s)",
                 args.min_sites, len(df), df.region.value_counts().to_dict())

    df["piS"] = np.where(df.S_sites > 0, df.pi_syn / df.S_sites, np.nan)
    df["piN"] = np.where(df.N_sites > 0, df.pi_non / df.N_sites, np.nan)

    # --- divergence, for MK -------------------------------------------------
    have_mk = False
    if args.kaks:
        k = pd.read_csv(args.kaks.expanduser(), sep="\t")
        k = k[["gene1", "Ka", "Ks"]].rename(columns={"gene1": "gene"})
        df = df.merge(k, on="gene", how="left")
        # Ka and Ks are rates; recover counts with this script's own callable
        # N and S so numerator and denominator refer to the same sites.
        df["Dn"] = df.Ka * df.N_sites
        df["Ds"] = df.Ks * df.S_sites
        have_mk = df.Ka.notna().sum() > 0
        logging.info("Genes with Ka/Ks: %d of %d", int(df.Ka.notna().sum()),
                     len(df))

    df.to_csv(outdir / "per_gene.tsv", sep="\t", index=False)

    # --- pooled by class ----------------------------------------------------
    logging.info("=" * 70)
    rows = []
    for region in ("retained", "eliminated"):
        sub = df[df.region == region]
        if sub.empty:
            continue
        r = pooled(sub)
        r["region"] = region
        ci = bootstrap_ci(sub, lambda d: (d.pi_non.sum() / d.N_sites.sum()) /
                          (d.pi_syn.sum() / d.S_sites.sum()), args.n_boot)
        r["piN_piS_lo"], r["piN_piS_hi"] = ci
        ci4 = bootstrap_ci(sub, lambda d: (d.pi_0.sum() / d.N0_sites.sum()) /
                           (d.pi_4.sum() / d.S4_sites.sum()), args.n_boot)
        r["pi0_pi4_lo"], r["pi0_pi4_hi"] = ci4
        rows.append(r)
        logging.info("%-11s genes=%5d  piN=%.5f  piS=%.5f  piN/piS=%.4f "
                     "[%.3f, %.3f]  Pn=%d Ps=%d",
                     region, r["genes"], r["piN"], r["piS"], r["piN_piS"],
                     ci[0], ci[1], int(r["Pn"]), int(r["Ps"]))
        logging.info("%-11s 0-fold/4-fold: pi0=%.5f  pi4=%.5f  pi0/pi4=%.4f "
                     "[%.3f, %.3f]  (%d 4-fold sites, %d 0-fold)",
                     "", r["pi0"], r["pi4"], r["pi0_pi4"], ci4[0], ci4[1],
                     int(r["S4_sites"]), int(r["N0_sites"]))
    summary = pd.DataFrame(rows)

    if len(summary) == 2:
        e = summary[summary.region == "eliminated"].iloc[0]
        t = summary[summary.region == "retained"].iloc[0]
        fold = e.piN_piS / t.piN_piS
        logging.info("")
        logging.info("piN/piS fold (eliminated / retained) = %.3f", fold)
        logging.info("  piS ratio  = %.3f   (mutation rate / Ne proxy)",
                     e.piS / t.piS)
        logging.info("  piN ratio  = %.3f", e.piN / t.piN)
        if t.pi0_pi4 == t.pi0_pi4 and e.pi0_pi4 == e.pi0_pi4:
            logging.info("")
            logging.info("  pi0/pi4 fold (eliminated / retained) = %.3f",
                         e.pi0_pi4 / t.pi0_pi4)
            logging.info("  pi4 ratio  = %.3f   (NG-weighting-free version of "
                         "the piS ratio)", e.pi4 / t.pi4)
            logging.info("  pi0/pi4 uses ONLY strictly 0-fold and 4-fold sites, "
                         "so it carries no Nei-Gojobori weighting assumption. "
                         "If it agrees with piN/piS, the counting model is not "
                         "driving the result; if pi4 shows the same depression "
                         "as piS, that depression is not a weighting artifact "
                         "either.")
        logging.info("Read these three together: relaxed purifying selection "
                     "predicts piN up with piS ~unchanged, exactly as dN rose "
                     "with dS unchanged in the divergence data.")

    # --- MK -----------------------------------------------------------------
    if have_mk:
        logging.info("")
        logging.info("McDonald-Kreitman, pooled within class")
        mk_rows = []
        variant_sets = [("all variants", "Pn", "Ps")]
        variant_sets += [(f"MAF>={t:g}", f"Pn_maf{i}", f"Ps_maf{i}")
                         for i, t in enumerate(mk_mafs)]
        for variant, pn_col, ps_col in variant_sets:
            logging.info("  -- %s", variant)
            for region in ("retained", "eliminated"):
                sub = df[(df.region == region) & df.Ka.notna()]
                if sub.empty:
                    continue
                m = mk_stats(sub[pn_col].sum(), sub[ps_col].sum(),
                             sub.Dn.sum(), sub.Ds.sum())
                m["region"] = region
                m["genes"] = len(sub)
                m["variant_set"] = variant
                lo, hi = bootstrap_ci(
                    sub, lambda d: mk_stats(d[pn_col].sum(), d[ps_col].sum(),
                                            d.Dn.sum(), d.Ds.sum())["DoS"],
                    args.n_boot)
                m["DoS_lo"], m["DoS_hi"] = lo, hi
                mk_rows.append(m)
                logging.info("     %-11s Pn=%5d Ps=%5d Dn=%8.1f Ds=%8.1f  "
                             "NI=%.3f  alpha=%.3f  DoS=%+.4f [%+.4f, %+.4f]",
                             region, int(m["Pn"]), int(m["Ps"]), m["Dn"],
                             m["Ds"], m["NI"], m["alpha"], m["DoS"], lo, hi)
        pd.DataFrame(mk_rows).to_csv(outdir / "mk_summary.tsv", sep="\t",
                                     index=False)
        logging.info("")
        logging.info("DoS is the headline, not NI: NI is a ratio of ratios and "
                     "misbehaves when any cell is small.")
        logging.info("If eliminated alpha keeps climbing toward retained as the "
                     "MAF floor rises, the two classes are indistinguishable in "
                     "alpha once the slightly-deleterious artifact is removed — "
                     "which is what relaxation predicts and adaptation does not. "
                     "The claim that survives every threshold is the weaker and "
                     "sufficient one: NO EXCESS of adaptive substitution in "
                     "eliminated genes.")
        logging.info("Quote the MAF-filtered row. Under relaxed selection the "
                     "eliminated class carries MORE low-frequency slightly "
                     "deleterious variants, which inflate Pn and deflate alpha "
                     "on their own — so an unfiltered alpha comparison is "
                     "partly guaranteed by the hypothesis it is testing.")
        logging.info("Overlapping DoS intervals between classes = relaxation, "
                     "not adaptation. That is the expected result, and it is "
                     "the claim the dN/dS data alone could not make.")

    summary.to_csv(outdir / "pi_summary.tsv", sep="\t", index=False)

    # --- figure -------------------------------------------------------------
    if len(summary) == 2:
        fig, ax = plt.subplots(figsize=(6.0, 5.0))
        vals = [summary[summary.region == "retained"].piN_piS.iloc[0],
                summary[summary.region == "eliminated"].piN_piS.iloc[0]]
        los = [summary[summary.region == "retained"].piN_piS_lo.iloc[0],
               summary[summary.region == "eliminated"].piN_piS_lo.iloc[0]]
        his = [summary[summary.region == "retained"].piN_piS_hi.iloc[0],
               summary[summary.region == "eliminated"].piN_piS_hi.iloc[0]]
        n = [int(summary[summary.region == "retained"].genes.iloc[0]),
             int(summary[summary.region == "eliminated"].genes.iloc[0])]
        ax.bar([0, 1], vals, color=[RETAINED, ELIMINATED], edgecolor="black",
               linewidth=0.6, width=0.6)
        ax.errorbar([0, 1], vals,
                    yerr=[np.array(vals) - np.array(los),
                          np.array(his) - np.array(vals)],
                    fmt="none", ecolor="black", capsize=6, linewidth=1.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Retained\n(n = {n[0]:,})",
                            f"Eliminated\n(n = {n[1]:,})"])
        ax.set_ylabel(r"$\pi_N/\pi_S$")
        for ext in ("png", "svg", "pdf"):
            fig.savefig(outdir / f"pin_pis_by_region.{ext}", bbox_inches="tight")
        plt.close(fig)
        logging.info("Figure: pin_pis_by_region.{png,svg,pdf}")

    logging.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
