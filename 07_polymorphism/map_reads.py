#!/usr/bin/env python3
"""
map_reads.py — map one sequencing run to the Ascaris v50 germline assembly.

Scope
-----
One run per invocation, driven by the ENA manifest written by
``fetch_ena_reads.py``. Parallelism across runs is SLURM's job (see
``run_map_reads.sh``, an array job), not this script's.

Why the germline assembly
-------------------------
Everything downstream (retained vs eliminated contrasts) lives in germline
coordinates. Mapping to a post-PDE somatic assembly would silently discard
every read from an eliminated region — the exact reads the analysis is about.
The script refuses to run if the reference name looks post-PDE, because this
mistake is invisible in the output and expensive to discover later.

Pipeline
--------
    bwa-mem2 mem -> samtools fixmate -m -> samtools sort -> samtools markdup

Duplicates are *marked*, never removed, so a later step can choose. Coordinate
sorting happens once, in the stream, so no intermediate BAM is written.

Pilot mode
----------
``--subsample-pairs N`` takes the first N read pairs of each FASTQ. This is for
the coverage audit only — it is a head, not a random sample, so reads come from
a small number of flowcell tiles and the duplicate rate is not meaningful. It
is fine for asking "does this worm carry germline DNA at all", which is the
only question the pilot needs to answer.

Example
-------
    python map_reads.py --manifest reads/ena_manifest.tsv --run SRR31757523 \
        --reads-dir reads --genome $DB_DIR/ascaris/v50/AG.v50.fa \
        --outdir bam --threads 16 --subsample-pairs 3000000
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Substrings that indicate a post-PDE / somatic reference. Mapping to one of
# these would zero out every eliminated region without any error message.
POST_PDE_HINTS = ("postpde", "post_pde", "somatic", "soma")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def require(tool: str) -> str:
    """Return the path to ``tool`` or exit with a readable message."""
    path = shutil.which(tool)
    if path is None:
        logging.error("Required tool not found on PATH: %s", tool)
        logging.error("Activate the popgen env, and remember the mandatory "
                      "`export PATH=.../envs/popgen/bin:$PATH` line — Spack "
                      "Python on Isaac shadows conda without it.")
        sys.exit(2)
    return path


def run(cmd: str, desc: str) -> None:
    """Run a shell pipeline, logging it verbatim so the job log is reproducible."""
    logging.info("%s", desc)
    logging.debug("CMD: %s", cmd)
    proc = subprocess.run(["bash", "-o", "pipefail", "-c", cmd],
                          capture_output=True, text=True)
    if proc.stderr.strip():
        for line in proc.stderr.strip().splitlines()[-40:]:
            logging.debug("  | %s", line)
    if proc.returncode != 0:
        logging.error("%s FAILED (rc=%d)", desc, proc.returncode)
        logging.error("CMD: %s", cmd)
        logging.error("stderr tail:\n%s", proc.stderr.strip()[-4000:])
        sys.exit(proc.returncode)


def load_manifest(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_sample_map(path: Path) -> dict[str, str]:
    """
    Read a run -> sample table: TSV/CSV with a header containing a run column
    (``run_accession`` / ``run`` / ``srr``) and a sample column (``sample`` /
    ``sample_id`` / ``worm`` / ``worm_id``). Extra columns are ignored.
    """
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter=delim))
    if not rows:
        logging.error("Sample map is empty: %s", path)
        sys.exit(2)

    cols = {c.strip().lower(): c for c in rows[0].keys()}
    run_col = next((cols[c] for c in ("run_accession", "run", "srr", "run_id")
                    if c in cols), None)
    smp_col = next((cols[c] for c in ("sample", "sample_id", "worm", "worm_id",
                                      "sample_name", "specimen")
                    if c in cols), None)
    if run_col is None or smp_col is None:
        logging.error("Sample map needs a run column and a sample column; "
                      "found: %s", ", ".join(rows[0].keys()))
        sys.exit(2)

    mapping = {r[run_col].strip(): r[smp_col].strip()
               for r in rows if r.get(run_col) and r.get(smp_col)}
    logging.info("Sample map: %d runs -> %d distinct samples (%s -> %s)",
                 len(mapping), len(set(mapping.values())), run_col, smp_col)
    return mapping


def resolve_sample(row: dict, rows: list[dict], sample_map: dict[str, str] | None,
                   field: str) -> str:
    """
    Decide the ``SM:`` value for this run.

    ``SM`` is what a joint genotyper groups reads by, so getting it wrong
    merges distinct worms into one individual and silently destroys every
    population-genetic statistic downstream.

    PRJNA1195031 is submitted with a SINGLE BioSample (SAMN45204864, alias
    '0101a') for all 62 runs, so the ENA metadata cannot identify worms. The
    per-worm identities live in Supplementary Data 2 of Landeryou et al. 2025.
    Pass that table with ``--sample-map``. Absent it, this falls back to the
    run accession — safe for the coverage pilot, where SM is never read, and
    it keeps runs separate rather than collapsing them.
    """
    run_acc = row["run_accession"]

    if sample_map is not None:
        if run_acc in sample_map:
            return sample_map[run_acc]
        logging.error("%s is not in the sample map. Every run must be mapped "
                      "before genotyping.", run_acc)
        sys.exit(2)

    if field != "auto":
        value = row.get(field) or run_acc
        return value

    alias_values = {r.get("sample_alias", "") for r in rows}
    if len(rows) > 1 and len(alias_values) == 1:
        logging.warning("=" * 70)
        logging.warning("ENA reports ONE sample_alias (%r) for all %d runs in "
                        "this manifest.", next(iter(alias_values)), len(rows))
        logging.warning("Using the run accession as SM so that runs stay "
                        "distinct. This is fine for the coverage audit.")
        logging.warning("BEFORE JOINT GENOTYPING: supply the real worm IDs via "
                        "--sample-map (Supplementary Data 2 of Landeryou et al. "
                        "2025), or every worm will be genotyped as one "
                        "individual.")
        logging.warning("=" * 70)
        return run_acc

    return row.get("sample_alias") or row.get("sample_accession") or run_acc


def pick_run(rows: list[dict], run_acc: str | None, index: int | None) -> dict:
    """Select a manifest row by accession or by 1-based array index."""
    if run_acc:
        for r in rows:
            if r["run_accession"] == run_acc:
                return r
        logging.error("Run %s not in manifest", run_acc)
        sys.exit(2)
    if index is None:
        logging.error("Give --run or --index")
        sys.exit(2)
    if not 1 <= index <= len(rows):
        logging.error("--index %d out of range 1..%d", index, len(rows))
        sys.exit(2)
    return rows[index - 1]


def fastqs_for(row: dict, reads_dir: Path) -> list[Path]:
    """
    Resolve the FASTQ paths for a run, as laid down by fetch_ena_reads.py
    (``<reads_dir>/<run_accession>/<filename>``).

    Runs with three files (ENA emits an unpaired ``_ .fastq.gz`` alongside
    ``_1``/``_2`` for some submissions) are reduced to the _1/_2 pair; the
    orphan file is a small, non-representative subset and mixing it in would
    give those reads a different duplicate and insert-size profile.
    """
    run_acc = row["run_accession"]
    names = [Path(u).name for u in row["fastq_ftp"].split(";") if u]
    paths = [reads_dir / run_acc / n for n in names]

    paired = [p for p in paths if p.name.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
    if len(paired) == 2:
        paths = sorted(paired)
    elif len(paths) > 2:
        logging.warning("%s: %d FASTQs and no clean _1/_2 pair; using the first two",
                        run_acc, len(paths))
        paths = sorted(paths)[:2]

    missing = [p for p in paths if not p.exists()]
    if missing:
        logging.error("Missing FASTQ(s): %s", ", ".join(str(m) for m in missing))
        logging.error("Run fetch_ena_reads.py first.")
        sys.exit(2)
    return paths


def subsample(src: Path, dest: Path, n_pairs: int) -> None:
    """
    Write the first ``n_pairs`` records of a gzipped FASTQ to ``dest``.

    ``set +o pipefail`` is deliberate: ``head`` closes the pipe early and
    ``zcat`` exits 141 (SIGPIPE). That is expected, not an error.
    """
    n_lines = n_pairs * 4
    cmd = (f"set +o pipefail; zcat {src} | head -n {n_lines} | "
           f"gzip -1 > {dest}")
    logging.info("Subsampling %s -> %s (%d pairs)", src.name, dest.name, n_pairs)
    subprocess.run(["bash", "-c", cmd], check=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Map one run to the germline assembly (bwa-mem2 + markdup).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--manifest", required=True, type=Path,
                   help="ena_manifest.tsv from fetch_ena_reads.py")
    p.add_argument("--reads-dir", required=True, type=Path,
                   help="Directory holding <run_accession>/*.fastq.gz")
    p.add_argument("--genome", required=True, type=Path,
                   help="Germline reference FASTA, e.g. $DB_DIR/ascaris/v50/AG.v50.fa")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Destination for the sorted, duplicate-marked BAM")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--run", help="Run accession to map")
    group.add_argument("--index", type=int,
                       help="1-based row index into the manifest (for SLURM arrays)")
    p.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)),
                   help="Threads for bwa-mem2 and samtools sort")
    p.add_argument("--subsample-pairs", type=int, default=0,
                   help="Map only the first N read pairs (0 = all). Pilot use only.")
    p.add_argument("--tmpdir", type=Path, default=None,
                   help="Scratch for sort spill and subsampled FASTQs "
                        "(default: $TMPDIR, else <outdir>/tmp)")
    p.add_argument("--aligner", default="bwa-mem2", choices=["bwa-mem2", "bwa"],
                   help="Aligner binary")
    p.add_argument("--sample-map", type=Path,
                   help="TSV/CSV mapping run accession -> worm/sample ID, used "
                        "for the BAM SM tag. REQUIRED before joint genotyping: "
                        "PRJNA1195031 carries one BioSample for all 62 runs, so "
                        "ENA metadata cannot identify individuals.")
    p.add_argument("--sample-field", default="auto",
                   choices=["auto", "sample_alias", "sample_accession",
                            "run_accession"],
                   help="Manifest column to use as SM when no --sample-map is "
                        "given; 'auto' detects a degenerate alias column and "
                        "falls back to run_accession")
    p.add_argument("--allow-postpde-reference", action="store_true",
                   help="Override the post-PDE reference guard (you almost "
                        "certainly do not want this)")
    p.add_argument("--force", action="store_true",
                   help="Re-map even if the output BAM already exists")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    aligner = require(args.aligner)
    samtools = require("samtools")

    genome = args.genome.expanduser().resolve()
    if not genome.exists():
        logging.error("Reference not found: %s", genome)
        return 2

    # Guard: germline coordinates are load-bearing for every downstream step.
    lowered = str(genome).lower()
    if any(h in lowered for h in POST_PDE_HINTS) and not args.allow_postpde_reference:
        logging.error("Reference path looks post-PDE/somatic: %s", genome)
        logging.error("Eliminated regions are absent from a post-PDE assembly, so "
                      "every read from them would be lost silently. Use the v50 "
                      "germline FASTA, or pass --allow-postpde-reference if you "
                      "really mean it.")
        return 2

    # bwa-mem2 index files; build once, separately (needs ~10 GB RAM).
    idx_suffixes = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"] \
        if args.aligner == "bwa-mem2" else [".amb", ".ann", ".bwt", ".pac", ".sa"]
    missing_idx = [s for s in idx_suffixes if not Path(str(genome) + s).exists()]
    if missing_idx:
        logging.error("Reference index incomplete (missing %s)", ", ".join(missing_idx))
        logging.error("Build it once with:  %s index %s", args.aligner, genome)
        return 2

    rows = load_manifest(args.manifest.expanduser())
    row = pick_run(rows, args.run, args.index)
    run_acc = row["run_accession"]

    sample_map = (load_sample_map(args.sample_map.expanduser())
                  if args.sample_map else None)
    sample = resolve_sample(row, rows, sample_map, args.sample_field)
    logging.info("Read-group SM for %s: %s", run_acc, sample)

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    tag = f"{run_acc}.sub{args.subsample_pairs}" if args.subsample_pairs else run_acc
    bam = outdir / f"{tag}.bam"

    if bam.exists() and not args.force:
        logging.info("%s already exists; nothing to do (use --force to re-map)", bam)
        return 0

    tmp_root = args.tmpdir or Path(os.environ.get("TMPDIR", outdir / "tmp"))
    tmp_root.mkdir(parents=True, exist_ok=True)

    fq = fastqs_for(row, args.reads_dir.expanduser().resolve())
    logging.info("Run %s (sample %s): %s", run_acc, sample,
                 ", ".join(p.name for p in fq))

    with tempfile.TemporaryDirectory(dir=tmp_root, prefix=f"{tag}_") as tmp:
        tmp = Path(tmp)

        if args.subsample_pairs:
            sub = []
            for i, src in enumerate(fq, 1):
                dest = tmp / f"{run_acc}_{i}.sub.fastq.gz"
                subsample(src, dest, args.subsample_pairs)
                sub.append(dest)
            fq = sub

        # Read group: SM is what the joint genotyper groups by, so it must be
        # the sample, not the run. ID is per-run so multi-run samples stay
        # distinguishable for duplicate marking and per-lane diagnostics.
        rg = (f"@RG\\tID:{run_acc}\\tSM:{sample}\\tPL:"
              f"{row.get('instrument_platform', 'ILLUMINA')}\\tLB:{sample}")

        sort_threads = max(1, args.threads // 2)
        cmd = (
            f"{aligner} mem -t {args.threads} -R '{rg}' "
            f"{genome} {' '.join(str(p) for p in fq)} "
            f"| {samtools} fixmate -m -u -@ 2 - - "
            # -m 1G per sort thread: with 16 CPUs that is 8 sort threads = 8 GB,
            # on top of the ~10 GB bwa-mem2 index. Raising either without the
            # other is how this job gets OOM-killed at 90% completion.
            f"| {samtools} sort -@ {sort_threads} -m 1G -T {tmp}/sort -u - "
            f"| {samtools} markdup -@ 2 -T {tmp}/mkdup - {bam}"
        )
        run(cmd, f"Mapping {run_acc} -> {bam.name}")

    run(f"{samtools} index -@ {min(args.threads, 8)} {bam}", "Indexing BAM")
    run(f"{samtools} flagstat -@ 4 {bam} > {bam.with_suffix('.flagstat.txt')}",
        "flagstat")
    run(f"{samtools} idxstats {bam} > {bam.with_suffix('.idxstats.txt')}",
        "idxstats")

    logging.info("Done: %s", bam)
    return 0


if __name__ == "__main__":
    sys.exit(main())
