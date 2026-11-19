#!/usr/bin/env python3
"""
select_germline_samples.py — turn the audit table into the analysis cohort.

Two jobs:

  1. **Select.** Keep worms whose gonad carryover is high enough that eliminated
     regions will genotype at usable depth, and whose retained-region depth
     shows the library mapped at all.
  2. **Annotate.** Emit the per-worm covariate table the downstream analysis
     needs: rescaled germline fraction and a sex call from X:A.

Both matter. Germline fraction is not just a filter — a worm at 0.72 and one at
0.36 have different power to call heterozygotes *in eliminated regions
specifically*, and that asymmetry falls directly on the retained-versus-
eliminated contrast if it is not carried forward and controlled.

Rescaling
---------
Raw germline fraction is bounded by the assay, not by biology. Wang-lab tissue
controls mapped through the identical pipeline give:

    floor   ~0.03   (female carcass 0.050, male intestine 0.012 — somatic)
    ceiling ~1.12   (male spermatids 1.119 — pure germline)

so ``rescaled = (raw - floor) / (ceiling - floor)``, clipped to [0, 1], is the
interpretable quantity: the fraction of the prep that is gonad.

Sex
---
X:A separates cleanly in this cohort (~0.5 male, ~1.0 female), but carries a
systematic upward bias of roughly +0.2 from mismapping and is read-length
sensitive. Thresholds are therefore deliberately wide, and anything landing
between them is called ``unknown`` rather than forced — an X:A of 0.82 has no
clean explanation and should not be silently binned.

Example
-------
    python select_germline_samples.py \
        --summary $POPGEN_DIR/audit_all_20260820/germline_fraction_summary.tsv \
        --min-germline 0.40 --outdir cohort_20260820
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import pandas as pd

#: Control runs — never part of the worm cohort, but reported for provenance.
CONTROL_RUNS = {"SRR502751", "SRR502756", "SRR5573034"}


def strip_suffix(name: str) -> str:
    """``SRR31692074.sub3000000`` -> ``SRR31692074``."""
    return re.sub(r"\.sub\d+$", "", str(name))


def call_sex(x_to_a: float, male_max: float, female_min: float) -> str:
    if pd.isna(x_to_a):
        return "unknown"
    if x_to_a <= male_max:
        return "M"
    if x_to_a >= female_min:
        return "F"
    return "unknown"


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Select the germline-positive cohort and annotate covariates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--summary", required=True, type=Path,
                   help="germline_fraction_summary.tsv from germline_fraction_audit.py")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-germline", type=float, default=0.40,
                   help="Minimum RAW germline fraction to include")
    p.add_argument("--min-retained-depth", type=float, default=1.0,
                   help="Minimum retained-window median depth in the pilot "
                        "subsample; catches libraries that barely mapped "
                        "(contamination) rather than libraries with no gonad")
    p.add_argument("--floor", type=float, default=0.03,
                   help="Somatic-control germline fraction (assay floor)")
    p.add_argument("--ceiling", type=float, default=1.12,
                   help="Spermatid-control germline fraction (assay ceiling)")
    p.add_argument("--male-max", type=float, default=0.65,
                   help="X:A at or below this is called male")
    p.add_argument("--female-min", type=float, default=0.90,
                   help="X:A at or above this is called female")
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

    df = pd.read_csv(args.summary.expanduser(), sep="\t")
    df["run"] = df["sample"].map(strip_suffix)
    df["is_control"] = df["run"].isin(CONTROL_RUNS)

    span = args.ceiling - args.floor
    if span <= 0:
        logging.error("--ceiling must exceed --floor")
        return 2
    df["germline_rescaled"] = ((df["germline_fraction"] - args.floor) / span).clip(0, 1)
    df["sex"] = df["x_to_autosome"].map(
        lambda v: call_sex(v, args.male_max, args.female_min))

    worms = df[~df["is_control"]].copy()
    logging.info("Worms: %d, controls: %d", len(worms), int(df["is_control"].sum()))

    # --- selection ---------------------------------------------------------
    enough_germline = worms["germline_fraction"] >= args.min_germline
    mapped = worms["depth_retained_median"] >= args.min_retained_depth
    worms["pass"] = enough_germline & mapped

    # Distinguish the two failure modes: they are different problems and a
    # reader of this table should not have to infer which applied.
    worms["fail_reason"] = ""
    worms.loc[~mapped, "fail_reason"] = "low_mapped_depth"
    worms.loc[mapped & ~enough_germline, "fail_reason"] = "somatic"

    sel = worms[worms["pass"]].sort_values("germline_fraction", ascending=False)

    logging.info("-" * 70)
    logging.info("Selected %d worm(s) at germline_fraction >= %.2f and "
                 "retained depth >= %.1fx", len(sel), args.min_germline,
                 args.min_retained_depth)
    for _, r in sel.iterrows():
        logging.info("  %-14s gf=%.3f  rescaled=%.2f  X:A=%.2f  sex=%s",
                     r["run"], r["germline_fraction"], r["germline_rescaled"],
                     r["x_to_autosome"], r["sex"])
    logging.info("Sex composition of the selection: %s",
                 sel["sex"].value_counts().to_dict())
    logging.info("Excluded: %s",
                 worms.loc[~worms["pass"], "fail_reason"].value_counts().to_dict())

    if len(sel) < 6:
        logging.warning("Fewer than 6 worms selected — pooled pi_N/pi_S and MK "
                        "will be badly underpowered. Consider lowering "
                        "--min-germline and reporting it as a sensitivity set.")

    # --- outputs -----------------------------------------------------------
    (outdir / "selected_runs.txt").write_text("\n".join(sel["run"]) + "\n")

    cols = ["run", "sex", "x_to_autosome", "germline_fraction",
            "germline_rescaled", "depth_retained_median",
            "depth_eliminated_median", "pass", "fail_reason", "is_control"]
    out = pd.concat([worms, df[df["is_control"]].assign(**{"pass": False,
                                                           "fail_reason": "control"})])
    out = out.sort_values(["is_control", "germline_fraction"],
                          ascending=[True, False])
    out[[c for c in cols if c in out.columns]].to_csv(
        outdir / "sample_metadata.tsv", sep="\t", index=False)

    logging.info("Wrote %s (%d runs) and %s",
                 outdir / "selected_runs.txt", len(sel),
                 outdir / "sample_metadata.tsv")
    logging.info("germline_rescaled is a COVARIATE for the genotyping stage, "
                 "not just a filter — carry it forward.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
