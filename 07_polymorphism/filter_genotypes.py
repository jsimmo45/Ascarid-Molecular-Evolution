#!/usr/bin/env python3
"""
filter_genotypes.py — filter the all-sites callset and split it for downstream use.

Two passes, because one of the filters cannot be expressed in bcftools syntax.

Pass 1 (bcftools): per-genotype depth/quality masking, site-level mapping and
bias filters, an indel buffer, and a missingness threshold.

Pass 2 (python): per-site **allele balance** across heterozygous calls. A
collapsed paralog — two diverged copies in the worm that the reference holds as
one locus — produces a site where every worm looks heterozygous with a clean
50/50 read split, or a skewed split that no real heterozygote should show. This
is the failure mode flagged by the depth audit: two worms had eliminated-region
depth *exceeding* retained (ratios 1.09 and 1.12), which the genmap mask cannot
catch because the duplication is in the sample, not the reference.

Expect this filter to remove proportionally more eliminated sites than retained
ones. **That asymmetry is reported, not hidden** — if it is large, the callable
denominator differs by class and every rate must be per callable site.

Outputs
-------
    snps_invariant.bcf   SNPs + invariant sites  -> piN/piS
    indels.bcf           coding indels           -> segregating LoF load
    excluded_sites.txt   sites dropped by the allele-balance screen
    filter_summary.tsv   site counts by region class at each stage

Why invariant sites are carried the whole way: pi is a per-site rate, and the
denominator is the number of sites *observed*, not the number in the reference.
Dropping invariant sites here would reintroduce the bias the all-sites callset
exists to avoid.

Example
-------
    python filter_genotypes.py --vcf vcf/allsites.bcf \\
        --eliminated-bed $DB_DIR/ascaris/v50/AG.v50.eliminated_strict.bed \\
        --outdir vcf/filtered --min-dp 6 --max-dp 30 --min-gq 20 --threads 8
"""

from __future__ import annotations

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
from pathlib import Path


def require(tool: str) -> str:
    p = shutil.which(tool)
    if p is None:
        logging.error("%s not found on PATH", tool)
        sys.exit(2)
    return p


def sh(cmd: str, desc: str, capture: bool = False) -> str:
    logging.info("%s", desc)
    logging.debug("CMD: %s", cmd)
    proc = subprocess.run(["bash", "-o", "pipefail", "-c", cmd],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        logging.error("%s FAILED (rc=%d)\n%s", desc, proc.returncode,
                      proc.stderr.strip()[-3000:])
        sys.exit(proc.returncode)
    return proc.stdout


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    raw: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.split("\t")
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    out = {}
    for c, ivs in raw.items():
        ivs.sort()
        merged = [ivs[0]]
        for s, e in ivs[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        out[c] = merged
    return out


def in_intervals(ivs: list[tuple[int, int]], pos0: int) -> bool:
    lo, hi = 0, len(ivs) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid]
        if pos0 < s:
            hi = mid - 1
        elif pos0 >= e:
            lo = mid + 1
        else:
            return True
    return False


# ---------------------------------------------------------------------------
# Pass 2: allele balance
# ---------------------------------------------------------------------------

def allele_balance_screen(bcftools: str, vcf: Path, dest: Path,
                          min_ab: float, max_ab: float, min_het: int,
                          max_het_frac: float, elim: dict) -> dict:
    """
    Stream heterozygous calls and flag sites whose allele balance is
    inconsistent with real heterozygosity.

    Two independent signals, both hallmarks of a collapsed duplication:

    * **Skewed balance.** A true heterozygote gives ~0.5 alt fraction. A
      paralog collapse gives whatever the copy ratio is — commonly 0.33 or
      0.67 for a 3-copy locus. Sites whose mean het allele balance falls
      outside [min_ab, max_ab] are dropped.
    * **Everybody heterozygous.** With 10 diploids, all-het under HWE needs
      p ~ 0.5 and has probability ~0.001. Seeing it repeatedly means the
      "heterozygosity" is fixed difference between copies, not polymorphism.

    Only sites with at least ``min_het`` het calls are judged; below that the
    balance estimate is too noisy to act on and the site is kept.
    """
    logging.info("Allele-balance screen (streaming het calls)")
    fmt = "%CHROM\\t%POS[\\t%GT,%AD]\\n"
    cmd = f"{bcftools} query -f '{fmt}' -i 'N_ALT>0' {vcf}"
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 20)

    n_seen = n_drop_ab = n_drop_allhet = 0
    drop_by_region = {"retained": 0, "eliminated": 0}
    seen_by_region = {"retained": 0, "eliminated": 0}

    with dest.open("w") as out:
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom, pos = f[0], int(f[1])
            n_het = 0
            n_called = 0
            ab_sum = 0.0
            for cell in f[2:]:
                if "," not in cell:
                    continue
                gt, _, ad = cell.partition(",")
                if gt in ("./.", ".|.", "."):
                    continue
                n_called += 1
                a = gt.replace("|", "/").split("/")
                if len(a) != 2 or a[0] == a[1]:
                    continue
                counts = [int(x) for x in ad.split(",") if x not in (".", "")]
                if len(counts) < 2:
                    continue
                tot = sum(counts)
                if tot == 0:
                    continue
                # allele balance of the two called alleles
                try:
                    i, j = int(a[0]), int(a[1])
                    d1, d2 = counts[i], counts[j]
                except (ValueError, IndexError):
                    continue
                if d1 + d2 == 0:
                    continue
                ab_sum += d2 / (d1 + d2)
                n_het += 1

            n_seen += 1
            region = ("eliminated"
                      if in_intervals(elim.get(chrom, []), pos - 1)
                      else "retained")
            seen_by_region[region] += 1

            if n_het < min_het:
                continue
            mean_ab = ab_sum / n_het
            all_het = n_called > 0 and (n_het / n_called) >= max_het_frac

            if not (min_ab <= mean_ab <= max_ab):
                out.write(f"{chrom}\t{pos}\n")
                n_drop_ab += 1
                drop_by_region[region] += 1
            elif all_het:
                out.write(f"{chrom}\t{pos}\n")
                n_drop_allhet += 1
                drop_by_region[region] += 1

    proc.wait()
    if proc.returncode != 0:
        logging.error("bcftools query failed:\n%s", proc.stderr.read()[-2000:])
        sys.exit(proc.returncode)

    logging.info("  variant sites examined : %d", n_seen)
    logging.info("  dropped, skewed balance: %d", n_drop_ab)
    logging.info("  dropped, all-het       : %d", n_drop_allhet)
    for reg in ("retained", "eliminated"):
        s, d = seen_by_region[reg], drop_by_region[reg]
        logging.info("  %-11s %7d variant sites, %6d dropped (%.2f%%)",
                     reg, s, d, 100 * d / s if s else 0.0)
    if seen_by_region["retained"] and seen_by_region["eliminated"]:
        fr = drop_by_region["retained"] / seen_by_region["retained"]
        fe = drop_by_region["eliminated"] / seen_by_region["eliminated"]
        if fr > 0:
            logging.info("  eliminated dropped %.2fx as often as retained",
                         fe / fr)
        if fe > 3 * max(fr, 1e-9):
            logging.warning("Eliminated sites are dropped far more often. "
                            "Residual collapsed paralogy is real; report the "
                            "per-class callable denominators explicitly and do "
                            "not quote any per-bp rate.")
    return {"seen": seen_by_region, "dropped": drop_by_region}


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Filter the all-sites callset; split SNPs+invariant from indels.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf", required=True, type=Path, help="allsites.bcf")
    p.add_argument("--eliminated-bed", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--threads", type=int, default=8)

    p.add_argument("--min-dp", type=int, default=6,
                   help="Per-genotype FORMAT/DP floor. BAMs are matched to 13x, "
                        "so 6 is ~0.45x mean — low enough not to bite in either "
                        "class, high enough to drop dropout-prone calls.")
    p.add_argument("--max-dp", type=int, default=30,
                   help="Per-genotype FORMAT/DP ceiling (~2.3x the matched 13x). "
                        "Excess depth after matching means collapsed copies.")
    p.add_argument("--min-gq", type=int, default=20)
    p.add_argument("--max-missing", type=int, default=1,
                   help="Maximum samples with a missing genotype at a site")
    p.add_argument("--min-mq", type=int, default=40, help="INFO/MQ floor")
    p.add_argument("--max-mq0f", type=float, default=0.10,
                   help="INFO/MQ0F ceiling (fraction of MAPQ0 reads)")
    p.add_argument("--indel-buffer", type=int, default=5,
                   help="Drop SNPs within this many bp of an indel — "
                        "misalignment around indels is the commonest source of "
                        "false heterozygotes")
    p.add_argument("--min-ab", type=float, default=0.25,
                   help="Lower bound on mean het allele balance")
    p.add_argument("--max-ab", type=float, default=0.75)
    p.add_argument("--min-het", type=int, default=3,
                   help="Judge allele balance only at sites with at least this "
                        "many het calls")
    p.add_argument("--max-het-frac", type=float, default=1.0,
                   help="Drop sites where this fraction or more of called "
                        "genotypes are heterozygous (1.0 = only all-het sites)")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    bcftools = require("bcftools")

    vcf = args.vcf.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    elim = read_bed(args.eliminated_bed.expanduser())

    def count(path: Path) -> int:
        return int(sh(f"{bcftools} index -n {path}", f"count {path.name}").strip())

    # --- pass 1 ------------------------------------------------------------
    # setGT first: mask individual genotypes on depth/quality, so the
    # missingness filter that follows counts genotypes we actually trust.
    step1 = outdir / "step1_gtmask.bcf"
    if not step1.exists():
        sh(f"{bcftools} +setGT {vcf} -O u -- -t q "
           f"-i 'FMT/DP<{args.min_dp} | FMT/DP>{args.max_dp} | FMT/GQ<{args.min_gq}' "
           f"-n . "
           f"| {bcftools} view -O b --threads {args.threads} -o {step1}",
           "mask low-confidence genotypes")
        sh(f"{bcftools} index -f --threads {args.threads} {step1}", "index")

    # Site-level filters. MQ/MQ0F are only informative where reads carry a
    # variant signal; invariant sites often lack them, so the expression keeps
    # a site when the tag is absent rather than discarding it.
    n_samples = len([s for s in sh(f"{bcftools} query -l {vcf}",
                                   "list samples").split() if s])
    if n_samples == 0:
        logging.error("No samples in %s", vcf)
        return 2
    logging.info("%d sample(s) in the callset", n_samples)

    step2 = outdir / "step2_sitefilter.bcf"
    # F_MISSING is a fraction, so the allowed count must be divided by the
    # actual sample count, not a hardcoded cohort size.
    max_missing_frac = min(1.0, args.max_missing / n_samples)
    logging.info("Allowing up to %d missing genotype(s) of %d "
                 "(F_MISSING <= %.3f)", args.max_missing, n_samples,
                 max_missing_frac)
    expr = (f"(INFO/MQ>={args.min_mq} || INFO/MQ=\".\") && "
            f"(INFO/MQ0F<={args.max_mq0f} || INFO/MQ0F=\".\") && "
            f"F_MISSING<={max_missing_frac:.4f}")
    if not step2.exists():
        sh(f"{bcftools} view -i '{expr}' -O b --threads {args.threads} "
           f"-o {step2} {step1}", "site-level filters")
        sh(f"{bcftools} index -f --threads {args.threads} {step2}", "index")

    # --- indel buffer ------------------------------------------------------
    # bcftools filter --SnpGap removes SNPs near indels; --IndelGap thins
    # clustered indels, which are usually one misalignment reported twice.
    step3 = outdir / "step3_gapfilter.bcf"
    if not step3.exists():
        sh(f"{bcftools} filter --SnpGap {args.indel_buffer} "
           f"--IndelGap {args.indel_buffer} -O b --threads {args.threads} "
           f"-o {step3} {step2}", f"SnpGap/IndelGap {args.indel_buffer}")
        sh(f"{bcftools} index -f --threads {args.threads} {step3}", "index")

    # --- pass 2 ------------------------------------------------------------
    excl = outdir / "excluded_sites.txt"
    stats = allele_balance_screen(bcftools, step3, excl, args.min_ab,
                                  args.max_ab, args.min_het,
                                  args.max_het_frac, elim)

    step4 = outdir / "step4_ab.bcf"
    if excl.stat().st_size > 0:
        sh(f"{bcftools} view -T ^{excl} -O b --threads {args.threads} "
           f"-o {step4} {step3}", "drop allele-balance failures")
    else:
        logging.info("No sites failed the allele-balance screen")
        shutil.copy(step3, step4)
    sh(f"{bcftools} index -f --threads {args.threads} {step4}", "index")

    # --- split -------------------------------------------------------------
    snps = outdir / "snps_invariant.bcf"
    indels = outdir / "indels.bcf"
    # -V indels keeps SNPs AND invariant sites (invariant have no ALT, so they
    # are not typed as indels). This is the file piN/piS reads.
    sh(f"{bcftools} view -V indels -O b --threads {args.threads} -o {snps} {step4}",
       "split: SNPs + invariant")
    sh(f"{bcftools} index -f --threads {args.threads} {snps}", "index")
    sh(f"{bcftools} view -v indels -O b --threads {args.threads} -o {indels} {step4}",
       "split: indels")
    sh(f"{bcftools} index -f --threads {args.threads} {indels}", "index")

    # --- summary -----------------------------------------------------------
    rows = [("input", count(vcf)), ("gt_masked", count(step1)),
            ("site_filtered", count(step2)), ("gap_filtered", count(step3)),
            ("allele_balance", count(step4)),
            ("snps_invariant", count(snps)), ("indels", count(indels))]
    with (outdir / "filter_summary.tsv").open("w") as fh:
        fh.write("stage\tsites\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")

    logging.info("=" * 70)
    prev = rows[0][1]
    for k, v in rows:
        logging.info("  %-16s %10d  (%.1f%% of input)", k, v, 100 * v / prev)
    n_var = int(sh(f"{bcftools} view -H -v snps {snps} | wc -l",
                   "count SNPs").strip())
    logging.info("  variant SNPs in snps_invariant.bcf: %d (%.2f%%)",
                 n_var, 100 * n_var / max(count(snps), 1))
    logging.info("")
    logging.info("piN/piS reads snps_invariant.bcf; the invariant sites in it "
                 "ARE the denominator. Do not filter them out later.")
    logging.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
