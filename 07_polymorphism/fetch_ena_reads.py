#!/usr/bin/env python3
"""
fetch_ena_reads.py — download raw FASTQs for an INSDC project from the ENA mirror.

Why ENA and not the SRA toolkit
-------------------------------
ENA mirrors every SRA submission and serves the original gzipped FASTQs over
plain HTTPS with published MD5 checksums. That means:

  * no ``prefetch`` / ``fasterq-dump`` cache to configure (the SRA toolkit's
    default cache lands in $HOME, which blows the Isaac home quota),
  * downloads are resumable with ``curl -C -`` after a walltime kill,
  * integrity is verifiable, so a truncated transfer is caught here rather than
    three steps later as a silent mapping-rate drop.

The script does two separable jobs, so a failed download never re-queries the
API and a re-run never re-downloads a file that is already complete:

  1. Query the ENA portal API and write a manifest TSV (one row per run).
  2. Download every FASTQ listed in the manifest, verifying MD5.

Typical use
-----------
    # 1. manifest only, so you can eyeball sample metadata and total size first
    python fetch_ena_reads.py --accession PRJNA1195031 \
        --outdir $POPGEN_DIR/reads --manifest-only

    # 2. subset to a pilot set, then pull everything
    python fetch_ena_reads.py --accession PRJNA1195031 \
        --outdir $POPGEN_DIR/reads --runs pilot_runs.txt --jobs 4

Exit status is non-zero if any run failed, so a SLURM wrapper can detect it.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import logging
import os
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

#: Metadata columns pulled from ENA. Keep ``fastq_ftp``/``fastq_md5``/
#: ``fastq_bytes`` — the downloader needs all three. The rest are for the
#: human reading the manifest.
ENA_FIELDS = [
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "sample_alias",
    "sample_title",
    "scientific_name",
    "library_strategy",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]

MANIFEST_NAME = "ena_manifest.tsv"


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def query_ena(accession: str, retries: int = 4, backoff: float = 5.0) -> str:
    """
    Fetch the ENA filereport for ``accession`` as TSV text.

    ENA occasionally returns 429/503 under load; retry with linear backoff
    rather than failing a queued job for a transient server hiccup.
    """
    params = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(ENA_FIELDS),
            "format": "tsv",
            "limit": 0,  # 0 = no cap
        }
    )
    url = f"{ENA_PORTAL}?{params}"
    logging.info("Querying ENA: %s", url)

    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                text = resp.read().decode("utf-8")
            if not text.strip():
                raise RuntimeError("ENA returned an empty report")
            return text
        except (urllib.error.URLError, RuntimeError, TimeoutError) as exc:
            last_err = exc
            logging.warning("ENA query attempt %d/%d failed: %s", attempt, retries, exc)
            if attempt < retries:
                time.sleep(backoff * attempt)
    raise RuntimeError(f"ENA query failed after {retries} attempts: {last_err}")


def parse_manifest(text: str) -> list[dict]:
    """Parse ENA TSV text into a list of row dicts, dropping runs with no FASTQ."""
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    kept = [r for r in rows if r.get("fastq_ftp")]
    dropped = len(rows) - len(kept)
    if dropped:
        logging.warning("%d run(s) have no fastq_ftp entry and were dropped", dropped)
    return kept


def write_manifest(rows: list[dict], path: Path) -> None:
    """Write the manifest TSV, preserving ENA column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ENA_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Manifest written: %s (%d runs)", path, len(rows))


def summarise(rows: list[dict]) -> None:
    """Log a short human summary: run count, total bases, total download size."""
    def _sum(field: str) -> int:
        total = 0
        for r in rows:
            for part in str(r.get(field, "")).split(";"):
                part = part.strip()
                if part.isdigit():
                    total += int(part)
        return total

    n_bases = _sum("base_count")
    n_bytes = _sum("fastq_bytes")
    species = sorted({r.get("scientific_name", "?") for r in rows})
    layouts = sorted({r.get("library_layout", "?") for r in rows})

    logging.info("Runs           : %d", len(rows))
    logging.info("Species        : %s", ", ".join(species))
    logging.info("Layout         : %s", ", ".join(layouts))
    logging.info("Total bases    : %.2f Gbp", n_bases / 1e9)
    logging.info("Download size  : %.1f GB", n_bytes / 1e9)


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def ftp_to_https(path: str) -> str:
    """
    Convert an ENA ``fastq_ftp`` path to an HTTPS URL.

    ENA reports paths like ``ftp.sra.ebi.ac.uk/vol1/fastq/...`` with no scheme.
    HTTPS is used rather than FTP because many HPC firewalls block passive FTP.
    """
    path = path.strip()
    if path.startswith(("http://", "https://")):
        return path
    if path.startswith("ftp://"):
        path = path[len("ftp://"):]
    return "https://" + path


def md5sum(path: Path, chunk: int = 8 << 20) -> str:
    """Streaming MD5 so a 5 GB FASTQ never lands in memory."""
    h = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def download_one(url: str, dest: Path, expected_md5: str | None,
                 retries: int, verify: bool) -> tuple[bool, str]:
    """
    Download ``url`` to ``dest`` with resume, then optionally verify MD5.

    Returns ``(ok, message)``. A file that already exists and checksums clean
    is skipped without re-downloading, which is what makes the whole script
    safe to re-run after a walltime kill.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)

    # Already present and correct? Do nothing.
    if dest.exists() and expected_md5 and verify:
        if md5sum(dest) == expected_md5:
            return True, f"{dest.name}: already complete (md5 ok)"
        logging.warning("%s exists but md5 mismatches — re-downloading", dest.name)
        dest.unlink()
    elif dest.exists() and not verify:
        return True, f"{dest.name}: already present (md5 not checked)"

    for attempt in range(1, retries + 1):
        # -C - resumes a partial file; -L follows redirects; --fail turns an
        # HTTP error page into a non-zero exit instead of a corrupt .gz.
        cmd = ["curl", "-sS", "--fail", "-L", "-C", "-",
               "--retry", "3", "--retry-delay", "5",
               "-o", str(dest), url]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode == 0:
            if not (expected_md5 and verify):
                return True, f"{dest.name}: downloaded"
            got = md5sum(dest)
            if got == expected_md5:
                return True, f"{dest.name}: downloaded (md5 ok)"
            logging.warning("%s md5 mismatch (%s != %s), attempt %d/%d",
                            dest.name, got, expected_md5, attempt, retries)
            dest.unlink(missing_ok=True)  # a resumed corrupt file never heals
        else:
            logging.warning("%s curl failed (rc=%d) attempt %d/%d: %s",
                            dest.name, proc.returncode, attempt, retries,
                            proc.stderr.strip()[:300])
        time.sleep(5 * attempt)

    return False, f"{dest.name}: FAILED after {retries} attempts"


def build_jobs(rows: list[dict], outdir: Path) -> list[tuple[str, Path, str | None]]:
    """
    Expand manifest rows into one (url, dest, md5) job per FASTQ file.

    ENA packs the mates of a paired run into semicolon-separated lists in
    ``fastq_ftp`` and ``fastq_md5``; the two lists are positionally matched.
    Files are laid out as ``<outdir>/<run_accession>/<filename>`` so a run is
    self-contained and easy to delete if it turns out to be junk.
    """
    jobs: list[tuple[str, Path, str | None]] = []
    for row in rows:
        run = row["run_accession"]
        urls = [u for u in row["fastq_ftp"].split(";") if u]
        md5s = [m for m in row.get("fastq_md5", "").split(";") if m]
        if md5s and len(md5s) != len(urls):
            logging.warning("%s: %d urls but %d md5s — checksums disabled for this run",
                            run, len(urls), len(md5s))
            md5s = []
        for i, url in enumerate(urls):
            https = ftp_to_https(url)
            dest = outdir / run / Path(urllib.parse.urlparse(https).path).name
            jobs.append((https, dest, md5s[i] if md5s else None))
    return jobs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download raw FASTQs for an INSDC project from ENA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--accession", required=True,
                   help="Project/study accession, e.g. PRJNA1195031")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Destination directory (reads land in <outdir>/<run>/)")
    p.add_argument("--manifest-only", action="store_true",
                   help="Write the manifest and stop; download nothing")
    p.add_argument("--runs", type=Path,
                   help="File with one run accession per line; restrict to these")
    p.add_argument("--jobs", type=int, default=4,
                   help="Concurrent downloads (ENA throttles above ~6)")
    p.add_argument("--retries", type=int, default=3,
                   help="Download attempts per file")
    p.add_argument("--no-verify", action="store_true",
                   help="Skip MD5 verification (faster, not recommended)")
    p.add_argument("--dry-run", action="store_true",
                   help="List what would be downloaded, then exit")
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

    if shutil.which("curl") is None:
        logging.error("curl not found on PATH — required for downloads")
        return 2

    outdir: Path = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / MANIFEST_NAME

    # --- 1. manifest -------------------------------------------------------
    if manifest_path.exists():
        logging.info("Reusing existing manifest: %s", manifest_path)
        rows = list(csv.DictReader(manifest_path.read_text().splitlines(),
                                   delimiter="\t"))
    else:
        rows = parse_manifest(query_ena(args.accession))
        if not rows:
            logging.error("No runs with FASTQ files found for %s", args.accession)
            return 1
        write_manifest(rows, manifest_path)

    if args.runs:
        wanted = {ln.strip() for ln in args.runs.read_text().splitlines() if ln.strip()}
        before = len(rows)
        rows = [r for r in rows if r["run_accession"] in wanted]
        logging.info("Restricted to %d of %d runs via %s", len(rows), before, args.runs)
        missing = wanted - {r["run_accession"] for r in rows}
        if missing:
            logging.warning("Not in manifest: %s", ", ".join(sorted(missing)))

    summarise(rows)

    if args.manifest_only:
        logging.info("--manifest-only set; stopping before download")
        return 0

    # --- 2. download -------------------------------------------------------
    jobs = build_jobs(rows, outdir)
    logging.info("%d FASTQ file(s) to fetch into %s", len(jobs), outdir)

    if args.dry_run:
        for url, dest, _ in jobs:
            print(f"{url}\t{dest}")
        return 0

    failures: list[str] = []
    verify = not args.no_verify
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(download_one, url, dest, md5, args.retries, verify): dest
            for url, dest, md5 in jobs
        }
        for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            ok, msg = fut.result()
            logging.info("[%d/%d] %s", i, len(jobs), msg)
            if not ok:
                failures.append(msg)

    if failures:
        logging.error("%d file(s) failed:", len(failures))
        for f in failures:
            logging.error("  %s", f)
        logging.error("Re-run the same command — completed files are skipped.")
        return 1

    logging.info("All %d file(s) downloaded and verified.", len(jobs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
