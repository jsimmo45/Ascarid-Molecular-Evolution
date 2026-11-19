#!/usr/bin/env python3
"""
Rebuild the four-way codon alignments that branch_dnds_four.py threw away.

WHY THIS EXISTS
---------------
branch_dnds_four.py aligns each quartet inside tempfile.mkdtemp(prefix=
'branch4_') and deletes the directory on exit unless --keep-tmp is passed. The
v2 run (SLURM 6089193, 2026-08-08) did not pass it, so the only quartet
alignments on disk are four_way_alignments/ -- 2,772 files, dated 07-31, built
BEFORE the Baylisascaris identifier fix. Those hold roughly half the
Baylisascaris genes and must not be used for anything.

This script regenerates them WITHOUT rerunning codeml, which is where the
1h32m went. Alignment alone is minutes.

WHERE THE QUARTETS COME FROM
----------------------------
branch_dnds_four_results_v2/branch_dnds_four_annotated.tsv, which already lists
the exact quartets the branch model used, one row per quartet with a gene ID
column per species. Taking the quartet list from there rather than re-deriving
it from the six rbh_*_v2/orthologs.tsv files guarantees that the clock test
runs on precisely the genes the dN/dS results came from. Re-deriving would risk
a set that differs by a few genes for no benefit.

Rows whose status column is not OK are skipped by default: those are quartets
where alignment or codeml failed in the original run, and they contribute
nothing here either.

OUTPUT
------
One FASTA per quartet, codon-aligned, four sequences each. Headers are written
as 'Species|gene_id' so that downstream species assignment is unambiguous --
molecular_clock_test.py then resolves them by exact token match rather than by
guessing at identifier formats.

USAGE
-----
Check the column detection first:

    python build_quartet_alignments.py \\
      --quartets branch_dnds_four_results_v2/branch_dnds_four_annotated.tsv \\
      --inspect

Then build:

    python build_quartet_alignments.py \\
      --quartets branch_dnds_four_results_v2/branch_dnds_four_annotated.tsv \\
      --outdir four_way_alignments_v2 \\
      --threads 24

GOTCHAS BAKED IN
----------------
  * Protein and CDS lengths must agree: a CDS shorter than 3 x (residues in the
    protein) means the two FASTAs disagree about that gene, and back-translation
    would silently frame-shift the tail. Such quartets are rejected, not
    patched.
  * Trailing stop codons are stripped before alignment. An internal stop means
    the CDS and protein are not the same gene model; those are rejected too.
  * MAFFT is run once per quartet in a private temp directory, so parallel
    workers cannot collide.
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, OrderedDict
from multiprocessing import Pool

import pandas as pd

SPECIES_ORDER = ['Ascaris', 'Parascaris', 'Baylisascaris', 'Toxocara']

# Candidate quartet-table column names per species, most specific first. The
# annotated table from branch_dnds_four.py uses as_gene / pu_gene / bsc_gene /
# tox_gene, but this list tolerates the other conventions in the project.
GENE_COL_CANDIDATES = {
    'Ascaris':       ['as_gene', 'ascaris_gene', 'gene_as', 'as_id', 'ascaris'],
    'Parascaris':    ['pu_gene', 'parascaris_gene', 'gene_pu', 'pu_id',
                      'parascaris'],
    'Baylisascaris': ['bsc_gene', 'baylisascaris_gene', 'gene_bsc', 'bsc_id',
                      'baylisascaris', 'bs_gene'],
    'Toxocara':      ['tox_gene', 'toxocara_gene', 'gene_tox', 'tox_id',
                      'toxocara', 'tc_gene'],
}

# Default FASTA names in the miniprot working directory.
DEFAULT_PROT = {'Ascaris': 'As_prot.fa', 'Parascaris': 'Pu_prot.fa',
                'Baylisascaris': 'Bsc_prot.fa', 'Toxocara': 'Tox_prot.fa'}
DEFAULT_CDS = {'Ascaris': 'As_cds.fa', 'Parascaris': 'Pu_cds.fa',
               'Baylisascaris': 'Bsc_cds.fa', 'Toxocara': 'Tox_cds.fa'}

STOPS = {'TAA', 'TAG', 'TGA'}

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

# Populated in main() before the worker pool forks. On Linux the children
# inherit these by copy-on-write, which avoids pickling four genomes per task.
PROT = {}
CDS = {}


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------

def read_fasta(path):
    """Read a FASTA into {id: sequence}, id truncated at first whitespace."""
    seqs = {}
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks).upper()
                name = line[1:].split()[0] if len(line) > 1 else ''
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(chunks).upper()
    return seqs


# ---------------------------------------------------------------------------
# Column detection
# ---------------------------------------------------------------------------

def detect_gene_columns(df, overrides=None):
    """
    Work out which table column holds each species' gene ID.

    Exact matches against GENE_COL_CANDIDATES come first; failing that, any
    column whose name contains the species name (or its short alias) and looks
    like an identifier column. Explicit --gene-col overrides always win.

    Returns {species: column_name}; raises SystemExit naming what is missing.
    """
    cols = {c.lower(): c for c in df.columns}
    found = {}

    for spec in overrides or []:
        if '=' not in spec:
            sys.exit(f"ERROR: --gene-col needs Species=column, got {spec!r}")
        sp, col = (x.strip() for x in spec.split('=', 1))
        if sp not in SPECIES_ORDER:
            sys.exit(f"ERROR: unknown species {sp!r} in --gene-col")
        if col not in df.columns:
            sys.exit(f"ERROR: column {col!r} not in the quartet table. "
                     f"Available: {list(df.columns)}")
        found[sp] = col

    for sp in SPECIES_ORDER:
        if sp in found:
            continue
        for cand in GENE_COL_CANDIDATES[sp]:
            if cand in cols:
                found[sp] = cols[cand]
                break

    missing = [s for s in SPECIES_ORDER if s not in found]
    if missing:
        sys.exit(f"ERROR: could not find gene-ID columns for {missing}.\n"
                 f"Columns present: {list(df.columns)}\n"
                 f"Set them explicitly, e.g. --gene-col 'Ascaris=as_gene'")
    return found


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

def align_proteins(records, workdir, mafft='mafft'):
    """
    Align protein sequences with MAFFT --auto.

    `records` is an OrderedDict of {label: protein string}. Returns the aligned
    dict in the same key order, or None if MAFFT failed or dropped a sequence.
    """
    infile = os.path.join(workdir, 'in.fa')
    with open(infile, 'w') as fh:
        for label, seq in records.items():
            fh.write(f'>{label}\n{seq}\n')

    try:
        proc = subprocess.run([mafft, '--auto', '--quiet', infile],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return None
    if proc.returncode != 0 or not proc.stdout:
        return None

    aligned, name, chunks = OrderedDict(), None, []
    for line in proc.stdout.splitlines():
        if line.startswith('>'):
            if name is not None:
                aligned[name] = ''.join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        aligned[name] = ''.join(chunks).upper()

    if len(aligned) != len(records):
        return None
    lengths = {len(v) for v in aligned.values()}
    if len(lengths) != 1:
        return None
    return OrderedDict((k, aligned[k]) for k in records if k in aligned)


def back_translate(aligned_prot, cds_by_label):
    """
    Map a protein alignment onto codons.

    Each non-gap residue consumes the next codon of that sequence's CDS; each
    gap becomes '---'. A CDS with fewer codons than the protein has residues
    means the protein and CDS FASTAs disagree about the gene model, so the
    quartet is rejected rather than silently truncated -- a frame-shifted tail
    would inflate both dN and dS for that gene.

    Returns {label: codon-aligned string} or None on any mismatch.
    """
    out = OrderedDict()
    for label, aln in aligned_prot.items():
        cds = cds_by_label[label]
        if len(cds) % 3:
            cds = cds[:len(cds) - (len(cds) % 3)]
        if cds[-3:] in STOPS:
            cds = cds[:-3]
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        if any(c in STOPS for c in codons):
            return None                       # internal stop: model mismatch

        n_res = len(aln.replace('-', ''))
        if len(codons) < n_res:
            return None                       # CDS too short for the protein

        pieces, pos = [], 0
        for aa in aln:
            if aa == '-':
                pieces.append('---')
            else:
                pieces.append(codons[pos])
                pos += 1
        out[label] = ''.join(pieces)

    lengths = {len(v) for v in out.values()}
    if len(lengths) != 1:
        return None
    return out


def process_quartet(task):
    """
    Build one quartet alignment. Runs in a worker process.

    Returns (index, status, n_codons). Status is 'OK' or a reason string; the
    caller tallies reasons so the log explains every dropped quartet.
    """
    idx, genes, outdir, mafft = task

    prots, cdss, labels = OrderedDict(), {}, {}
    for sp in SPECIES_ORDER:
        gid = genes[sp]
        p = PROT[sp].get(gid)
        c = CDS[sp].get(gid)
        if p is None:
            return idx, f'missing_protein_{sp}', 0
        if c is None:
            return idx, f'missing_cds_{sp}', 0
        label = f'{sp}|{gid}'
        labels[sp] = label
        prots[label] = p.rstrip('*').replace('*', 'X')
        cdss[label] = c

    if any(len(v) == 0 for v in prots.values()):
        return idx, 'empty_sequence', 0

    workdir = tempfile.mkdtemp(prefix='quartet_')
    try:
        aligned = align_proteins(prots, workdir, mafft)
        if aligned is None:
            return idx, 'alignment_failed', 0
        codon = back_translate(aligned, cdss)
        if codon is None:
            return idx, 'back_translation_failed', 0

        name = genes[SPECIES_ORDER[0]].replace('/', '_').replace('|', '_')
        path = os.path.join(outdir, f'{name}.fa')
        with open(path, 'w') as fh:
            for label, seq in codon.items():
                fh.write(f'>{label}\n{seq}\n')
        return idx, 'OK', len(next(iter(codon.values()))) // 3
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Rebuild four-way codon alignments from the v2 quartet '
                    'table, without rerunning codeml.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--quartets', required=True,
                   help='branch_dnds_four_results_v2/'
                        'branch_dnds_four_annotated.tsv')
    p.add_argument('--outdir', default='four_way_alignments_v2')
    p.add_argument('--threads', type=int, default=8,
                   help='Parallel MAFFT workers (default: 8)')

    for sp in SPECIES_ORDER:
        p.add_argument(f'--prot-{sp.lower()}', default=DEFAULT_PROT[sp],
                       help=f'{sp} protein FASTA (default: {DEFAULT_PROT[sp]})')
        p.add_argument(f'--cds-{sp.lower()}', default=DEFAULT_CDS[sp],
                       help=f'{sp} CDS FASTA (default: {DEFAULT_CDS[sp]})')

    p.add_argument('--gene-col', action='append', default=None,
                   metavar='SPECIES=COLUMN',
                   help="Override the quartet-table column for one species, "
                        "e.g. --gene-col 'Ascaris=as_gene'. Repeatable.")
    p.add_argument('--status-col', default='status',
                   help="Column marking successful quartets (default: "
                        "'status'). Pass 'none' to use every row.")
    p.add_argument('--require-status', default='OK',
                   help="Keep only rows whose status column equals this "
                        "(default: OK)")
    p.add_argument('--expect-min', type=int, default=5000,
                   help='Abort if fewer than this many quartets survive '
                        'filtering. Guards against silently running on a '
                        'pre-fix table (default: 5000; v2 has 5,954)')
    p.add_argument('--mafft', default='mafft')
    p.add_argument('--inspect', action='store_true',
                   help='Report detected columns, row counts and FASTA hit '
                        'rates, then exit')
    return p.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.quartets):
        sys.exit(f"ERROR: quartet table not found: {args.quartets}")
    if shutil.which(args.mafft) is None:
        sys.exit(f"ERROR: {args.mafft} not found in PATH. "
                 f"conda activate kaks first.")

    df = pd.read_csv(args.quartets, sep='\t')
    logger.info(f"Loaded {len(df):,} rows from {args.quartets}")

    gene_cols = detect_gene_columns(df, args.gene_col)
    for sp in SPECIES_ORDER:
        logger.info(f"  {sp:<15} -> column {gene_cols[sp]!r}")

    if args.status_col != 'none' and args.status_col in df.columns:
        before = len(df)
        df = df[df[args.status_col].astype(str) == args.require_status]
        logger.info(f"Status filter {args.status_col}=={args.require_status!r}: "
                    f"{before:,} -> {len(df):,}")
    elif args.status_col != 'none':
        logger.warning(f"No {args.status_col!r} column; using every row")

    df = df.dropna(subset=list(gene_cols.values()))
    logger.info(f"Quartets with all four gene IDs: {len(df):,}")

    prot_paths = {sp: getattr(args, f'prot_{sp.lower()}')
                  for sp in SPECIES_ORDER}
    cds_paths = {sp: getattr(args, f'cds_{sp.lower()}') for sp in SPECIES_ORDER}
    for sp in SPECIES_ORDER:
        for kind, path in (('protein', prot_paths[sp]), ('CDS', cds_paths[sp])):
            if not os.path.isfile(path):
                sys.exit(f"ERROR: {sp} {kind} FASTA not found: {path}")

    logger.info("Loading FASTAs")
    for sp in SPECIES_ORDER:
        PROT[sp] = read_fasta(prot_paths[sp])
        CDS[sp] = read_fasta(cds_paths[sp])
        ids = df[gene_cols[sp]].astype(str)
        hit_p = ids.isin(PROT[sp]).mean()
        hit_c = ids.isin(CDS[sp]).mean()
        logger.info(f"  {sp:<15} prot {len(PROT[sp]):>6,} seqs "
                    f"({100 * hit_p:5.1f}% of quartet IDs found)   "
                    f"cds {len(CDS[sp]):>6,} ({100 * hit_c:5.1f}%)")
        if hit_p < 0.5 or hit_c < 0.5:
            logger.warning(f"  {sp}: fewer than half the quartet IDs are in "
                           f"its FASTA. Wrong file, or the table uses a "
                           f"different ID convention.")

    if args.inspect:
        logger.info("\n--inspect: stopping before alignment.")
        return

    if len(df) < args.expect_min:
        sys.exit(f"ERROR: only {len(df):,} quartets, expected at least "
                 f"{args.expect_min:,}. The v2 table has 5,954; the pre-fix "
                 f"one had 2,771. Check --quartets points at the _v2 results. "
                 f"Override with --expect-min if this is deliberate.")

    os.makedirs(args.outdir, exist_ok=True)
    logging.getLogger().addHandler(logging.FileHandler(
        os.path.join(args.outdir, 'build_quartet_alignments.log')))

    tasks = []
    for i, row in enumerate(df.itertuples(index=False)):
        genes = {sp: str(getattr(row, gene_cols[sp].replace(' ', '_')))
                 for sp in SPECIES_ORDER}
        tasks.append((i, genes, args.outdir, args.mafft))

    logger.info(f"Aligning {len(tasks):,} quartets on {args.threads} workers")
    counts = Counter()
    codon_lens = []
    with Pool(args.threads) as pool:
        for n, (_, status, ncod) in enumerate(
                pool.imap_unordered(process_quartet, tasks, chunksize=16), 1):
            counts[status] += 1
            if status == 'OK':
                codon_lens.append(ncod)
            if n % 500 == 0:
                logger.info(f"  {n:,}/{len(tasks):,}  ok={counts['OK']:,}")

    logger.info(f"\nWrote {counts['OK']:,} alignments to {args.outdir}/")
    if codon_lens:
        logger.info(f"Median alignment length: "
                    f"{int(pd.Series(codon_lens).median()):,} codons")
    for status, n in counts.most_common():
        if status != 'OK':
            logger.info(f"  dropped {status}: {n:,}")

    if counts['OK'] < args.expect_min:
        logger.warning(f"Only {counts['OK']:,} alignments succeeded, below "
                       f"--expect-min {args.expect_min:,}. Check the dropped "
                       f"reasons above before using these downstream.")
    logger.info('Done.')


if __name__ == '__main__':
    main()
