#!/usr/bin/env python3
"""
reciprocal_best_hits.py — Identify 1:1 orthologs via reciprocal best BLAST hits.

Steps:
  1. Build BLAST databases for both species (blastp)
  2. Run bidirectional blastp
  3. Parse for reciprocal best hits (RBH)
  4. Output pairs TSV for kaks_pipeline.py

Usage:
  python reciprocal_best_hits.py \
    --prot1 As_prot.fa --prot2 Pu_prot.fa \
    --label1 Ascaris --label2 Parascaris \
    --threads 48 \
    --evalue 1e-10 \
    --outdir rbh_results

Output:
  rbh_results/
    orthologs.tsv          — pairs file (input for kaks_pipeline.py)
    blast_1vs2.tsv         — raw BLAST species1 vs species2
    blast_2vs1.tsv         — raw BLAST species2 vs species1
    rbh_stats.txt          — summary statistics
    rbh_identity_hist.svg  — alignment identity distribution

Dependencies:
  - NCBI BLAST+ (makeblastdb, blastp)
  - Python: pandas, matplotlib, numpy
"""

import argparse
import logging
import os
import subprocess
import sys
import shutil

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ===========================================================================
# BLAST database construction
# ===========================================================================

def make_blast_db(fasta_path, dbtype='prot', label=''):
    """
    Build a BLAST database from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to input FASTA.
    dbtype : str
        'prot' or 'nucl'.
    label : str
        Label for logging.

    Returns
    -------
    str : path to the database (same as fasta_path).
    """
    logging.info(f"Building BLAST database for {label}: {fasta_path}")

    cmd = [
        'makeblastdb',
        '-in', fasta_path,
        '-dbtype', dbtype,
        '-parse_seqids'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"makeblastdb failed: {result.stderr}")
        sys.exit(1)

    logging.info(f"  Database built successfully")
    return fasta_path


# ===========================================================================
# BLAST execution
# ===========================================================================

BLAST_OUTFMT = '6'
BLAST_COLUMNS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
]


def run_blastp(query, db, outfile, threads=8, evalue=1e-10, max_target_seqs=5):
    """
    Run blastp.

    Parameters
    ----------
    query : str
        Path to query FASTA.
    db : str
        Path to BLAST database.
    outfile : str
        Path to output file.
    threads : int
        Number of threads.
    evalue : float
        E-value threshold.
    max_target_seqs : int
        Maximum target sequences to report.

    Returns
    -------
    pd.DataFrame : BLAST results.
    """
    logging.info(f"Running blastp: {os.path.basename(query)} vs {os.path.basename(db)}")

    cmd = [
        'blastp',
        '-query', query,
        '-db', db,
        '-outfmt', BLAST_OUTFMT,
        '-evalue', str(evalue),
        '-max_target_seqs', str(max_target_seqs),
        '-num_threads', str(threads),
        '-out', outfile
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"blastp failed: {result.stderr}")
        sys.exit(1)

    # Parse results
    if os.path.getsize(outfile) == 0:
        logging.warning(f"  BLAST produced no hits")
        return pd.DataFrame(columns=BLAST_COLUMNS)

    df = pd.read_csv(outfile, sep='\t', header=None, names=BLAST_COLUMNS)
    logging.info(f"  {len(df)} total hits, {df['qseqid'].nunique()} query genes with hits")

    return df


def fasta_id_map(fasta_path):
    """
    Build {lowercased_id: exact_id} from a FASTA.

    BLAST rewrites sequence identifiers that look like NCBI accessions. An id
    such as 'g10000' (one letter + five digits) matches the accession pattern
    and comes back UPPERCASED as 'G10000', while 'g9999' (four digits) is
    treated as a local id and keeps its case. The forward and reverse searches
    can therefore report the same gene with different capitalisation, and the
    reciprocal merge -- an exact string join -- silently drops every such pair.

    Mapping BLAST's ids back onto the FASTA's own ids removes the whole class
    of problem and guarantees the output ids match the sequence files that the
    downstream Ka/Ks step reads.
    """
    id_map = {}
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                real = line[1:].split()[0]
                id_map[real.lower()] = real
    return id_map


def normalize_blast_ids(df, qmap, smap):
    """Rewrite qseqid/sseqid to the exact identifiers used in the FASTAs."""
    if len(df) == 0:
        return df
    df = df.copy()
    n_q = (~df['qseqid'].isin(qmap.values())).sum()
    n_s = (~df['sseqid'].isin(smap.values())).sum()
    df['qseqid'] = df['qseqid'].map(lambda x: qmap.get(str(x).lower(), x))
    df['sseqid'] = df['sseqid'].map(lambda x: smap.get(str(x).lower(), x))
    if n_q or n_s:
        logging.info(f"  Normalized ids to match FASTA: "
                     f"{n_q} query, {n_s} subject")
    return df


# ===========================================================================
# Best hit extraction
# ===========================================================================

def get_best_hits(blast_df):
    """
    Extract the single best hit per query gene (highest bitscore).

    If multiple hits share the same top bitscore, take the one with
    highest percent identity as tiebreaker.

    Parameters
    ----------
    blast_df : pd.DataFrame
        Raw BLAST output.

    Returns
    -------
    pd.DataFrame : one row per query gene.
    """
    if len(blast_df) == 0:
        return blast_df

    # Sort by bitscore (desc), then pident (desc) for tiebreaking
    df_sorted = blast_df.sort_values(
        ['qseqid', 'bitscore', 'pident'],
        ascending=[True, False, False]
    )

    # Keep first (best) hit per query
    best = df_sorted.groupby('qseqid').first().reset_index()
    logging.info(f"  {len(best)} best hits extracted")

    return best


# ===========================================================================
# Reciprocal best hit identification
# ===========================================================================

def find_reciprocal_best_hits(best_1vs2, best_2vs1):
    """
    Identify reciprocal best hits (RBH) between two BLAST searches.

    A pair (A, B) is an RBH if:
      - A's best hit in species2 is B
      - B's best hit in species1 is A

    Parameters
    ----------
    best_1vs2 : pd.DataFrame
        Best hits from species1 -> species2 BLAST.
    best_2vs1 : pd.DataFrame
        Best hits from species2 -> species1 BLAST.

    Returns
    -------
    pd.DataFrame : RBH pairs with columns:
        gene1, gene2, pident_fwd, pident_rev, bitscore_fwd, bitscore_rev,
        evalue_fwd, evalue_rev, align_len_fwd, align_len_rev
    """
    logging.info("Finding reciprocal best hits...")

    # Forward: gene1 -> gene2
    fwd = best_1vs2[['qseqid', 'sseqid', 'pident', 'bitscore', 'evalue', 'length']].copy()
    fwd.columns = ['gene1', 'gene2', 'pident_fwd', 'bitscore_fwd', 'evalue_fwd', 'align_len_fwd']

    # Reverse: gene2 -> gene1
    rev = best_2vs1[['qseqid', 'sseqid', 'pident', 'bitscore', 'evalue', 'length']].copy()
    rev.columns = ['gene2', 'gene1', 'pident_rev', 'bitscore_rev', 'evalue_rev', 'align_len_rev']

    # Merge on gene1 + gene2 — only pairs present in both directions survive
    rbh = fwd.merge(rev, on=['gene1', 'gene2'], how='inner')

    logging.info(f"  {len(fwd)} forward best hits")
    logging.info(f"  {len(rev)} reverse best hits")
    logging.info(f"  {len(rbh)} reciprocal best hits (1:1 orthologs)")

    # Sort by mean bitscore
    rbh['mean_bitscore'] = (rbh['bitscore_fwd'] + rbh['bitscore_rev']) / 2
    rbh = rbh.sort_values('mean_bitscore', ascending=False).reset_index(drop=True)

    return rbh


# ===========================================================================
# Plotting
# ===========================================================================

def plot_identity_distribution(rbh_df, outpath, label1='Species1', label2='Species2'):
    """Plot percent identity distribution of RBH pairs."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Percent identity histogram
    ax = axes[0]
    mean_pident = (rbh_df['pident_fwd'] + rbh_df['pident_rev']) / 2
    ax.hist(mean_pident, bins=50, color='#4C72B0', edgecolor='black', linewidth=0.3)
    ax.axvline(mean_pident.median(), color='red', linestyle='--', linewidth=1.5,
               label=f'Median = {mean_pident.median():.1f}%')
    ax.set_xlabel('Mean percent identity (%)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Number of ortholog pairs', fontsize=14, fontweight='bold')
    ax.set_title(f'RBH Identity: {label1} vs {label2}', fontsize=14, fontweight='bold')
    ax.tick_params(axis='both', labelsize=12)
    ax.legend(fontsize=11, frameon=False)
    ax.text(0.05, 0.95, f'n = {len(rbh_df):,}',
            transform=ax.transAxes, fontsize=12, ha='left', va='top')

    # Bitscore scatter: forward vs reverse
    ax = axes[1]
    ax.scatter(rbh_df['bitscore_fwd'], rbh_df['bitscore_rev'],
               s=6, alpha=0.3, color='#55A868', edgecolors='none')
    max_val = max(rbh_df['bitscore_fwd'].max(), rbh_df['bitscore_rev'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1, alpha=0.5)
    ax.set_xlabel(f'{label1} → {label2} bitscore', fontsize=14, fontweight='bold')
    ax.set_ylabel(f'{label2} → {label1} bitscore', fontsize=14, fontweight='bold')
    ax.set_title('Forward vs Reverse Bitscore', fontsize=14, fontweight='bold')
    ax.tick_params(axis='both', labelsize=12)

    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved identity plot: {outpath}")


# ===========================================================================
# Main
# ===========================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify 1:1 orthologs via reciprocal best BLAST hits.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python reciprocal_best_hits.py \\
    --prot1 As_prot.fa --prot2 Pu_prot.fa \\
    --label1 Ascaris --label2 Parascaris \\
    --threads 48 --outdir rbh_results

Then feed the output into kaks_pipeline.py:
  python kaks_pipeline.py \\
    --pairs rbh_results/orthologs.tsv \\
    --cds1 As_cds.fa --prot1 As_prot.fa \\
    --cds2 Pu_cds.fa --prot2 Pu_prot.fa ...
        """
    )

    parser.add_argument('--prot1', required=True,
                        help='Protein FASTA for species 1 (e.g., Ascaris)')
    parser.add_argument('--prot2', required=True,
                        help='Protein FASTA for species 2 (e.g., Parascaris)')
    parser.add_argument('--label1', default='Species1',
                        help='Label for species 1 (default: Species1)')
    parser.add_argument('--label2', default='Species2',
                        help='Label for species 2 (default: Species2)')
    parser.add_argument('--evalue', type=float, default=1e-10,
                        help='E-value threshold (default: 1e-10)')
    parser.add_argument('--max-target-seqs', type=int, default=5,
                        help='Max target seqs per query (default: 5)')
    parser.add_argument('--min-pident', type=float, default=30.0,
                        help='Minimum percent identity for RBH (default: 30)')
    parser.add_argument('--min-align-len', type=int, default=50,
                        help='Minimum alignment length in AA (default: 50)')
    parser.add_argument('--threads', type=int, default=8,
                        help='BLAST threads (default: 8)')
    parser.add_argument('--outdir', default='rbh_results',
                        help='Output directory (default: rbh_results)')
    parser.add_argument('--skip-blast', action='store_true',
                        help='Skip BLAST if output files already exist')

    return parser.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Logging
    log_path = os.path.join(args.outdir, 'rbh.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("=" * 60)
    logging.info("Reciprocal Best Hit Analysis")
    logging.info("=" * 60)
    logging.info(f"Arguments: {vars(args)}")

    # Check BLAST+ availability
    for tool in ['makeblastdb', 'blastp']:
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: {tool} not found in PATH. Load BLAST+ module or install.")

    blast_1vs2_path = os.path.join(args.outdir, 'blast_1vs2.tsv')
    blast_2vs1_path = os.path.join(args.outdir, 'blast_2vs1.tsv')

    if args.skip_blast and os.path.exists(blast_1vs2_path) and os.path.exists(blast_2vs1_path):
        logging.info("Skipping BLAST (--skip-blast, output files exist)")
        blast_1vs2 = pd.read_csv(blast_1vs2_path, sep='\t', header=None, names=BLAST_COLUMNS)
        blast_2vs1 = pd.read_csv(blast_2vs1_path, sep='\t', header=None, names=BLAST_COLUMNS)
        map1 = fasta_id_map(args.prot1)
        map2 = fasta_id_map(args.prot2)
        blast_1vs2 = normalize_blast_ids(blast_1vs2, map1, map2)
        blast_2vs1 = normalize_blast_ids(blast_2vs1, map2, map1)
    else:
        # Step 1: Build databases
        make_blast_db(args.prot1, label=args.label1)
        make_blast_db(args.prot2, label=args.label2)

        # Step 2: Bidirectional BLAST
        blast_1vs2 = run_blastp(
            args.prot1, args.prot2, blast_1vs2_path,
            threads=args.threads, evalue=args.evalue,
            max_target_seqs=args.max_target_seqs
        )
        blast_2vs1 = run_blastp(
            args.prot2, args.prot1, blast_2vs1_path,
            threads=args.threads, evalue=args.evalue,
            max_target_seqs=args.max_target_seqs
        )
        map1 = fasta_id_map(args.prot1)
        map2 = fasta_id_map(args.prot2)
        blast_1vs2 = normalize_blast_ids(blast_1vs2, map1, map2)
        blast_2vs1 = normalize_blast_ids(blast_2vs1, map2, map1)

    # Step 3: Extract best hits
    best_1vs2 = get_best_hits(blast_1vs2)
    best_2vs1 = get_best_hits(blast_2vs1)

    # Step 4: Find RBHs
    rbh = find_reciprocal_best_hits(best_1vs2, best_2vs1)

    if len(rbh) == 0:
        logging.error("No reciprocal best hits found. Check input files and parameters.")
        sys.exit(1)

    # Step 5: Filter by identity and alignment length
    n_pre = len(rbh)
    mean_pident = (rbh['pident_fwd'] + rbh['pident_rev']) / 2
    mean_alen = (rbh['align_len_fwd'] + rbh['align_len_rev']) / 2

    rbh = rbh[mean_pident >= args.min_pident]
    n_after_pident = len(rbh)
    logging.info(f"  Filtered by pident >= {args.min_pident}%: "
                 f"{n_pre} -> {n_after_pident}")

    mean_alen = (rbh['align_len_fwd'] + rbh['align_len_rev']) / 2
    rbh = rbh[mean_alen >= args.min_align_len]
    n_final = len(rbh)
    logging.info(f"  Filtered by align_len >= {args.min_align_len} AA: "
                 f"{n_after_pident} -> {n_final}")

    # Step 6: Save outputs
    # Full RBH table
    rbh_path = os.path.join(args.outdir, 'rbh_full.tsv')
    rbh.to_csv(rbh_path, sep='\t', index=False)
    logging.info(f"Full RBH table saved: {rbh_path}")

    # Pairs file for kaks_pipeline.py (just gene1 and gene2, no header)
    pairs_path = os.path.join(args.outdir, 'orthologs.tsv')
    rbh[['gene1', 'gene2']].to_csv(pairs_path, sep='\t', index=False, header=False)
    logging.info(f"Pairs file saved: {pairs_path}")

    # Step 7: Plot
    plot_identity_distribution(
        rbh,
        os.path.join(args.outdir, 'rbh_identity.svg'),
        label1=args.label1, label2=args.label2
    )

    # Step 8: Summary statistics
    mean_pident = (rbh['pident_fwd'] + rbh['pident_rev']) / 2
    summary_path = os.path.join(args.outdir, 'rbh_stats.txt')
    with open(summary_path, 'w') as fh:
        fh.write("Reciprocal Best Hit Summary\n")
        fh.write("=" * 40 + "\n")
        fh.write(f"{args.label1} proteins: {blast_1vs2['qseqid'].nunique()} with BLAST hits\n")
        fh.write(f"{args.label2} proteins: {blast_2vs1['qseqid'].nunique()} with BLAST hits\n")
        fh.write(f"RBH pairs (raw): {n_pre}\n")
        fh.write(f"RBH pairs (filtered): {n_final}\n")
        fh.write(f"Filters: pident >= {args.min_pident}%, "
                 f"align_len >= {args.min_align_len} AA\n")
        fh.write(f"\nPercent identity:\n")
        fh.write(f"  Median: {mean_pident.median():.1f}%\n")
        fh.write(f"  Mean:   {mean_pident.mean():.1f}%\n")
        fh.write(f"  Min:    {mean_pident.min():.1f}%\n")
        fh.write(f"  Max:    {mean_pident.max():.1f}%\n")
        fh.write(f"\nOutput files:\n")
        fh.write(f"  Pairs for Ka/Ks: {pairs_path}\n")
        fh.write(f"  Full RBH table:  {rbh_path}\n")

    logging.info(f"Summary saved: {summary_path}")
    logging.info(f"Done. Use {pairs_path} as --pairs input to kaks_pipeline.py")


if __name__ == '__main__':
    main()
