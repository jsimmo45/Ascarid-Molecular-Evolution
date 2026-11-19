#!/usr/bin/env python3
"""
kaks_pipeline.py — Ka/Ks calculation pipeline for ortholog pairs.

Automates:
  1. Protein alignment (MAFFT)
  2. Codon back-translation (built-in, no PAL2NAL needed)
  3. Ka/Ks calculation (KaKs_Calculator or PAML yn00)
  4. Filtering
  5. Plotting (Ks distribution, Ka/Ks histogram, genomic context)

Input:
  - Ortholog pairs file (TSV): species1_gene_id <tab> species2_gene_id
  - CDS FASTA for each species (nucleotide coding sequences)
  - Protein FASTA for each species

Output:
  - TSV of Ka, Ks, Ka/Ks per ortholog pair
  - Ks distribution histogram (SVG)
  - Ka/Ks distribution histogram (SVG)
  - Optional: Ka/Ks by genomic region if BED annotation provided

Dependencies:
  - Python: biopython, matplotlib, numpy, pandas
  - External: mafft, and one of [KaKs_Calculator, yn00 (PAML)]

Usage:
  python kaks_pipeline.py \\
    --pairs orthologs.tsv \\
    --cds1 ascaris_cds.fa --prot1 ascaris_prot.fa \\
    --cds2 parascaris_cds.fa --prot2 parascaris_prot.fa \\
    --method KaKs_Calculator \\
    --threads 8 \\
    --outdir kaks_results

Author: Generated for Ascaris/Parascaris comparative genomics
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Biopython imports
# ---------------------------------------------------------------------------
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit("ERROR: Biopython is required. Install with: pip install biopython")


# ===========================================================================
# FASTA I/O helpers
# ===========================================================================

def load_fasta(fasta_path):
    """Load FASTA file into dict {id: SeqRecord}."""
    logging.info(f"Loading FASTA: {fasta_path}")
    records = OrderedDict()
    for rec in SeqIO.parse(fasta_path, "fasta"):
        records[rec.id] = rec
    logging.info(f"  Loaded {len(records)} sequences")
    return records


def write_fasta(records, outpath):
    """Write list of SeqRecords to FASTA file."""
    with open(outpath, 'w') as fh:
        SeqIO.write(records, fh, "fasta")


# ===========================================================================
# Protein alignment with MAFFT
# ===========================================================================

def align_proteins(prot1_rec, prot2_rec, tmpdir):
    """
    Align two protein sequences using MAFFT.

    Parameters
    ----------
    prot1_rec, prot2_rec : SeqRecord
        Protein sequences for the ortholog pair.
    tmpdir : str
        Temporary directory for intermediate files.

    Returns
    -------
    dict : {gene_id: aligned_protein_string} or None on failure.
    """
    pair_id = f"{prot1_rec.id}__vs__{prot2_rec.id}"
    infile = os.path.join(tmpdir, f"{pair_id}.prot.fa")
    outfile = os.path.join(tmpdir, f"{pair_id}.prot.aln.fa")

    write_fasta([prot1_rec, prot2_rec], infile)

    cmd = ["mafft", "--auto", "--quiet", infile]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            logging.warning(f"  MAFFT failed for {pair_id}: {result.stderr.strip()}")
            return None
    except FileNotFoundError:
        sys.exit("ERROR: mafft not found in PATH. Install MAFFT first.")
    except subprocess.TimeoutExpired:
        logging.warning(f"  MAFFT timed out for {pair_id}")
        return None

    # Write and parse alignment
    with open(outfile, 'w') as fh:
        fh.write(result.stdout)

    aligned = OrderedDict()
    for rec in SeqIO.parse(outfile, "fasta"):
        aligned[rec.id] = str(rec.seq)

    if len(aligned) != 2:
        logging.warning(f"  Unexpected alignment output for {pair_id}")
        return None

    return aligned


# ===========================================================================
# Codon back-translation
# ===========================================================================

def back_translate(aligned_proteins, cds_records):
    """
    Convert a protein alignment to a codon alignment using CDS sequences.

    For each aligned protein, walk through the alignment: gaps map to '---',
    residues map to the next codon from the CDS.

    Parameters
    ----------
    aligned_proteins : dict
        {gene_id: aligned_protein_string} from MAFFT output.
    cds_records : dict
        {gene_id: SeqRecord} of nucleotide CDS sequences.

    Returns
    -------
    dict : {gene_id: codon_aligned_string} or None on failure.
    """
    codon_aln = OrderedDict()

    for gene_id, aln_prot in aligned_proteins.items():
        if gene_id not in cds_records:
            logging.warning(f"  CDS not found for {gene_id}")
            return None

        cds_seq = str(cds_records[gene_id].seq).upper()

        # Remove trailing stop codon from CDS if present
        last_codon = cds_seq[-3:]
        if last_codon in ("TAA", "TAG", "TGA"):
            cds_seq = cds_seq[:-3]

        # Remove stop codon character from protein alignment if present
        aln_prot_clean = aln_prot.rstrip('*')

        # Count non-gap residues in alignment
        n_residues = len(aln_prot_clean.replace('-', ''))
        n_codons = len(cds_seq) // 3

        if n_residues != n_codons:
            logging.warning(
                f"  Length mismatch for {gene_id}: "
                f"{n_residues} aligned residues vs {n_codons} CDS codons "
                f"(CDS len = {len(cds_seq)} nt)"
            )
            return None

        # Build codon alignment
        codon_str = []
        cds_pos = 0
        for aa in aln_prot_clean:
            if aa == '-':
                codon_str.append('---')
            else:
                codon = cds_seq[cds_pos * 3 : cds_pos * 3 + 3]
                if len(codon) != 3:
                    logging.warning(f"  Ran out of codons for {gene_id}")
                    return None
                codon_str.append(codon)
                cds_pos += 1

        codon_aln[gene_id] = ''.join(codon_str)

    # Verify all codon alignments are the same length
    lengths = [len(v) for v in codon_aln.values()]
    if len(set(lengths)) != 1:
        logging.warning(f"  Codon alignment length mismatch: {lengths}")
        return None

    return codon_aln


# ===========================================================================
# AXT format conversion (for KaKs_Calculator)
# ===========================================================================

def write_axt(codon_aln, pair_name, outpath):
    """
    Write codon alignment in AXT format for KaKs_Calculator.

    AXT format:
      >pair_name
      seq1_codons
      seq2_codons
    """
    ids = list(codon_aln.keys())
    with open(outpath, 'w') as fh:
        fh.write(f">{pair_name}\n")
        fh.write(f"{codon_aln[ids[0]]}\n")
        fh.write(f"{codon_aln[ids[1]]}\n")


# ===========================================================================
# PAML yn00 format conversion
# ===========================================================================

def write_yn00_input(codon_aln, outpath):
    """Write codon alignment in sequential PHYLIP format for yn00."""
    ids = list(codon_aln.keys())
    seq_len = len(codon_aln[ids[0]])
    with open(outpath, 'w') as fh:
        fh.write(f"  2  {seq_len}\n")
        for gid in ids:
            # yn00 needs names <=30 chars, no special characters
            safe_name = gid[:30].ljust(30)
            fh.write(f"{safe_name}  {codon_aln[gid]}\n")


def write_yn00_ctl(seqfile, outfile, ctlpath):
    """Write yn00 control file."""
    ctl = f"""\
      seqfile = {seqfile}
      outfile = {outfile}
      verbose = 0
      icode = 0
      weighting = 0
      commonf3x4 = 0
"""
    with open(ctlpath, 'w') as fh:
        fh.write(ctl)


# ===========================================================================
# Ka/Ks calculation — single pair
# ===========================================================================

def calc_kaks_pair(gene1_id, gene2_id, prot_db1, prot_db2, cds_db1, cds_db2,
                   method, tmpdir):
    """
    Calculate Ka/Ks for a single ortholog pair.

    Returns
    -------
    dict with keys: gene1, gene2, Ka, Ks, Ka_Ks, align_len, status
    """
    pair_name = f"{gene1_id}__vs__{gene2_id}"
    result = {
        'gene1': gene1_id, 'gene2': gene2_id,
        'Ka': np.nan, 'Ks': np.nan, 'Ka_Ks': np.nan,
        'align_len': 0, 'status': 'FAIL'
    }

    # --- Check sequences exist ---
    if gene1_id not in prot_db1:
        result['status'] = f'protein_not_found:{gene1_id}'
        return result
    if gene2_id not in prot_db2:
        result['status'] = f'protein_not_found:{gene2_id}'
        return result
    if gene1_id not in cds_db1:
        result['status'] = f'cds_not_found:{gene1_id}'
        return result
    if gene2_id not in cds_db2:
        result['status'] = f'cds_not_found:{gene2_id}'
        return result

    # --- Align proteins ---
    aligned = align_proteins(prot_db1[gene1_id], prot_db2[gene2_id], tmpdir)
    if aligned is None:
        result['status'] = 'alignment_failed'
        return result

    # --- Back-translate ---
    cds_combined = {gene1_id: cds_db1[gene1_id], gene2_id: cds_db2[gene2_id]}
    codon_aln = back_translate(aligned, cds_combined)
    if codon_aln is None:
        result['status'] = 'backtranslation_failed'
        return result

    # Alignment length in codons (excluding shared gap columns)
    ids = list(codon_aln.keys())
    seq1, seq2 = codon_aln[ids[0]], codon_aln[ids[1]]
    n_codons = 0
    for i in range(0, len(seq1), 3):
        c1, c2 = seq1[i:i+3], seq2[i:i+3]
        if c1 != '---' and c2 != '---':
            n_codons += 1
    result['align_len'] = n_codons

    if n_codons < 30:
        result['status'] = f'too_short:{n_codons}_codons'
        return result

    # --- Calculate Ka/Ks ---
    if method == 'KaKs_Calculator':
        result = _run_kaks_calculator(codon_aln, pair_name, tmpdir, result)
    elif method == 'yn00':
        result = _run_yn00(codon_aln, pair_name, tmpdir, result)
    else:
        result['status'] = f'unknown_method:{method}'

    return result


def _run_kaks_calculator(codon_aln, pair_name, tmpdir, result):
    """Run KaKs_Calculator on a codon alignment."""
    axt_file = os.path.join(tmpdir, f"{pair_name}.axt")
    out_file = os.path.join(tmpdir, f"{pair_name}.kaks")

    write_axt(codon_aln, pair_name, axt_file)

    cmd = ["KaKs_Calculator", "-i", axt_file, "-o", out_file, "-m", "YN"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if proc.returncode != 0:
            result['status'] = f'KaKs_Calculator_error:{proc.stderr.strip()[:80]}'
            return result
    except FileNotFoundError:
        sys.exit("ERROR: KaKs_Calculator not found in PATH.")
    except subprocess.TimeoutExpired:
        result['status'] = 'KaKs_Calculator_timeout'
        return result

    # Parse output (TSV with header)
    try:
        df = pd.read_csv(out_file, sep='\t')
        if len(df) == 0:
            result['status'] = 'KaKs_Calculator_empty_output'
            return result
        row = df.iloc[0]
        result['Ka'] = float(row['Ka']) if row['Ka'] != 'NA' else np.nan
        result['Ks'] = float(row['Ks']) if row['Ks'] != 'NA' else np.nan
        if not np.isnan(result['Ka']) and not np.isnan(result['Ks']) and result['Ks'] > 0:
            result['Ka_Ks'] = result['Ka'] / result['Ks']
        result['status'] = 'OK'
    except Exception as e:
        result['status'] = f'parse_error:{str(e)[:80]}'

    return result


def _run_yn00(codon_aln, pair_name, tmpdir, result):
    """Run PAML yn00 on a codon alignment."""
    # yn00 must be run from its own directory (writes temp files)
    workdir = os.path.join(tmpdir, f"yn00_{pair_name}")
    os.makedirs(workdir, exist_ok=True)

    seqfile = os.path.join(workdir, "input.phy")
    outfile = os.path.join(workdir, "yn00.out")
    ctlfile = os.path.join(workdir, "yn00.ctl")

    write_yn00_input(codon_aln, seqfile)
    write_yn00_ctl(seqfile, outfile, ctlfile)

    cmd = ["yn00", ctlfile]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=60, cwd=workdir)
        if proc.returncode != 0:
            result['status'] = f'yn00_error:{proc.stderr.strip()[:80]}'
            return result
    except FileNotFoundError:
        sys.exit("ERROR: yn00 (PAML) not found in PATH.")
    except subprocess.TimeoutExpired:
        result['status'] = 'yn00_timeout'
        return result

    # Parse yn00 output — look for Yang & Nielsen (2000) results
    # yn00 output format (after header "seq. seq.  S  N  t  kappa  omega  dN +- SE  dS +- SE"):
    #    2    1  5580.0 13404.0   0.1737  4.1171  0.0918 0.0148 +- 0.0011  0.1614 +- 0.0061
    # indices: 0=seq2 1=seq1 2=S 3=N 4=t 5=kappa 6=omega 7=dN 8=+- 9=SE_dN 10=dS 11=+- 12=SE_dS
    try:
        with open(outfile, 'r') as fh:
            lines = fh.readlines()

        ka_val, ks_val = np.nan, np.nan
        found_header = False
        for line in lines:
            # Look for the data header line in the YN section
            if "seq. seq." in line and "dN" in line and "dS" in line:
                found_header = True
                continue
            if found_header and line.strip():
                parts = line.split()
                if len(parts) >= 11:
                    try:
                        # Verify this is a data line (starts with seq indices)
                        int(parts[0])
                        int(parts[1])
                        ka_val = float(parts[7])   # dN
                        ks_val = float(parts[10])  # dS
                        break
                    except (ValueError, IndexError):
                        pass

        result['Ka'] = ka_val
        result['Ks'] = ks_val
        if not np.isnan(ka_val) and not np.isnan(ks_val) and ks_val > 0:
            result['Ka_Ks'] = ka_val / ks_val
        result['status'] = 'OK'

    except Exception as e:
        result['status'] = f'yn00_parse_error:{str(e)[:80]}'

    return result


# ===========================================================================
# Batch processing with parallelism
# ===========================================================================

def _worker(args):
    """Worker function for parallel processing (unpacks tuple)."""
    return calc_kaks_pair(*args)


def run_pipeline(pairs, prot_db1, prot_db2, cds_db1, cds_db2,
                 method, threads, tmpdir):
    """
    Run Ka/Ks calculation on all ortholog pairs.

    Parameters
    ----------
    pairs : list of (gene1_id, gene2_id) tuples
    method : str, 'KaKs_Calculator' or 'yn00'
    threads : int
    tmpdir : str

    Returns
    -------
    pd.DataFrame with columns: gene1, gene2, Ka, Ks, Ka_Ks, align_len, status
    """
    logging.info(f"Processing {len(pairs)} ortholog pairs with {threads} threads")
    logging.info(f"Ka/Ks method: {method}")

    results = []

    if threads == 1:
        # Sequential — easier to debug
        for i, (g1, g2) in enumerate(pairs):
            if (i + 1) % 200 == 0:
                logging.info(f"  Processed {i + 1}/{len(pairs)} pairs")
            res = calc_kaks_pair(g1, g2, prot_db1, prot_db2, cds_db1, cds_db2,
                                method, tmpdir)
            results.append(res)
    else:
        # Parallel — each worker gets its own temp subdirectory
        tasks = []
        for i, (g1, g2) in enumerate(pairs):
            worker_tmp = os.path.join(tmpdir, f"worker_{i}")
            os.makedirs(worker_tmp, exist_ok=True)
            tasks.append((g1, g2, prot_db1, prot_db2, cds_db1, cds_db2,
                          method, worker_tmp))

        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(_worker, t): i for i, t in enumerate(tasks)}
            done_count = 0
            for future in as_completed(futures):
                done_count += 1
                if done_count % 200 == 0:
                    logging.info(f"  Completed {done_count}/{len(pairs)} pairs")
                try:
                    results.append(future.result())
                except Exception as e:
                    idx = futures[future]
                    g1, g2 = pairs[idx]
                    logging.error(f"  Exception for {g1} vs {g2}: {e}")
                    results.append({
                        'gene1': g1, 'gene2': g2,
                        'Ka': np.nan, 'Ks': np.nan, 'Ka_Ks': np.nan,
                        'align_len': 0, 'status': f'exception:{str(e)[:80]}'
                    })

    df = pd.DataFrame(results)
    return df


# ===========================================================================
# Filtering
# ===========================================================================

def filter_results(df, min_align_len=100, max_ks=5.0, max_kaks=50.0):
    """
    Filter Ka/Ks results for reliability.

    Parameters
    ----------
    df : pd.DataFrame
        Raw results from run_pipeline().
    min_align_len : int
        Minimum alignment length in codons.
    max_ks : float
        Maximum Ks (saturated synonymous sites are unreliable).
    max_kaks : float
        Maximum Ka/Ks to exclude extreme outliers / calculation artifacts.

    Returns
    -------
    pd.DataFrame : filtered results.
    """
    n_total = len(df)
    n_ok = (df['status'] == 'OK').sum()
    logging.info(f"Filtering: {n_total} total pairs, {n_ok} successfully calculated")

    df_filt = df[df['status'] == 'OK'].copy()
    n_after_status = len(df_filt)

    df_filt = df_filt[df_filt['Ks'] > 0]
    n_after_ks_zero = len(df_filt)
    logging.info(f"  Removed {n_after_status - n_after_ks_zero} pairs with Ks = 0")

    df_filt = df_filt[df_filt['Ks'] <= max_ks]
    n_after_ks_max = len(df_filt)
    logging.info(f"  Removed {n_after_ks_zero - n_after_ks_max} pairs with Ks > {max_ks}")

    df_filt = df_filt[df_filt['Ka_Ks'] <= max_kaks]
    n_after_kaks = len(df_filt)
    logging.info(f"  Removed {n_after_ks_max - n_after_kaks} pairs with Ka/Ks > {max_kaks}")

    df_filt = df_filt[df_filt['align_len'] >= min_align_len]
    n_final = len(df_filt)
    logging.info(f"  Removed {n_after_kaks - n_final} pairs with align_len < {min_align_len}")

    logging.info(f"  Final: {n_final} pairs passed filtering "
                 f"({n_final/n_total*100:.1f}% of total)")

    # Summary statistics
    logging.info(f"  Median Ks = {df_filt['Ks'].median():.4f}")
    logging.info(f"  Median Ka/Ks = {df_filt['Ka_Ks'].median():.4f}")
    logging.info(f"  Mean Ks = {df_filt['Ks'].mean():.4f}")

    return df_filt


# ===========================================================================
# Plotting
# ===========================================================================

def plot_ks_distribution(df, outpath, species1="Ascaris", species2="Parascaris"):
    """
    Plot Ks distribution histogram.

    The modal peak corresponds to the speciation divergence event.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    ks_vals = df['Ks'].dropna()
    bins = np.arange(0, min(ks_vals.max() + 0.05, 3.05), 0.03)

    ax.hist(ks_vals, bins=bins, color='#4C72B0', edgecolor='black',
            linewidth=0.3, alpha=0.85)

    # Mark median
    median_ks = ks_vals.median()
    ax.axvline(median_ks, color='red', linestyle='--', linewidth=1.5,
               label=f'Median Ks = {median_ks:.3f}')

    ax.set_xlabel('Ks (synonymous substitutions/site)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Number of ortholog pairs', fontsize=14, fontweight='bold')
    ax.set_title(f'Ks Distribution: {species1} vs {species2}',
                 fontsize=16, fontweight='bold')
    ax.tick_params(axis='both', labelsize=12)
    ax.legend(fontsize=12, frameon=False)
    ax.set_xlim(0, 3)

    # Add n count
    ax.text(0.95, 0.95, f'n = {len(ks_vals):,}',
            transform=ax.transAxes, fontsize=12, ha='right', va='top')

    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved Ks distribution plot: {outpath}")


def plot_kaks_distribution(df, outpath, species1="Ascaris", species2="Parascaris"):
    """
    Plot Ka/Ks distribution histogram with selection regime annotations.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    kaks_vals = df['Ka_Ks'].dropna()
    # Cap display at 3 for visualization
    kaks_display = kaks_vals[kaks_vals <= 3.0]
    bins = np.arange(0, 3.05, 0.03)

    ax.hist(kaks_display, bins=bins, color='#55A868', edgecolor='black',
            linewidth=0.3, alpha=0.85)

    # Reference lines for selection regimes
    ax.axvline(1.0, color='red', linestyle='--', linewidth=1.5, label='Neutral (Ka/Ks = 1)')

    median_kaks = kaks_vals.median()
    ax.axvline(median_kaks, color='orange', linestyle='--', linewidth=1.5,
               label=f'Median Ka/Ks = {median_kaks:.3f}')

    ax.set_xlabel('Ka/Ks (dN/dS)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Number of ortholog pairs', fontsize=14, fontweight='bold')
    ax.set_title(f'Ka/Ks Distribution: {species1} vs {species2}',
                 fontsize=16, fontweight='bold')
    ax.tick_params(axis='both', labelsize=12)
    ax.legend(fontsize=12, frameon=False)
    ax.set_xlim(0, 3)

    # Annotate selection regimes
    ax.text(0.3, 0.85, 'Purifying\nselection', transform=ax.transAxes,
            fontsize=11, ha='center', va='center', color='#4C72B0', fontstyle='italic')
    ax.text(0.85, 0.85, 'Positive\nselection', transform=ax.transAxes,
            fontsize=11, ha='center', va='center', color='#C44E52', fontstyle='italic')

    n_positive = (kaks_vals > 1.0).sum()
    n_purifying = (kaks_vals < 1.0).sum()
    ax.text(0.95, 0.70,
            f'n = {len(kaks_vals):,}\n'
            f'Purifying: {n_purifying:,} ({n_purifying/len(kaks_vals)*100:.1f}%)\n'
            f'Positive: {n_positive:,} ({n_positive/len(kaks_vals)*100:.1f}%)',
            transform=ax.transAxes, fontsize=10, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved Ka/Ks distribution plot: {outpath}")


def plot_ka_vs_ks(df, outpath, species1="Ascaris", species2="Parascaris"):
    """
    Scatter plot of Ka vs Ks with neutral evolution line.
    """
    fig, ax = plt.subplots(figsize=(8, 7))

    ka = df['Ka'].dropna()
    ks = df['Ks'].dropna()
    mask = ka.index.intersection(ks.index)
    ka, ks = ka[mask], ks[mask]

    ax.scatter(ks, ka, s=8, alpha=0.3, color='#4C72B0', edgecolors='none')

    # Neutral line (Ka = Ks)
    max_val = max(ks.max(), ka.max()) * 1.05
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1.5, label='Ka = Ks (neutral)')

    ax.set_xlabel('Ks', fontsize=14, fontweight='bold')
    ax.set_ylabel('Ka', fontsize=14, fontweight='bold')
    ax.set_title(f'Ka vs Ks: {species1} vs {species2}',
                 fontsize=16, fontweight='bold')
    ax.tick_params(axis='both', labelsize=12)
    ax.legend(fontsize=12, frameon=False, loc='upper left')

    ax.text(0.95, 0.05, f'n = {len(ka):,}',
            transform=ax.transAxes, fontsize=12, ha='right', va='bottom')

    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved Ka vs Ks plot: {outpath}")


def plot_kaks_by_region(df, bed_path, outpath):
    """
    Compare Ka/Ks between genomic regions (e.g., CBR-adjacent vs genome-wide).

    The BED file should have columns: gene_id, region_label
    (e.g., 'CBR_adjacent', 'retained_interior', 'eliminated')
    """
    if bed_path is None or not os.path.exists(bed_path):
        logging.info("No region annotation provided, skipping regional Ka/Ks plot")
        return

    logging.info(f"Loading region annotations: {bed_path}")
    regions = pd.read_csv(bed_path, sep='\t', header=None,
                          names=['gene_id', 'region_label'])

    # Merge with Ka/Ks results (try both gene1 and gene2 columns)
    df_merged = df.merge(regions, left_on='gene1', right_on='gene_id', how='left')
    if df_merged['region_label'].isna().all():
        df_merged = df.merge(regions, left_on='gene2', right_on='gene_id', how='left')

    df_plot = df_merged.dropna(subset=['region_label', 'Ka_Ks'])
    if len(df_plot) == 0:
        logging.warning("No genes matched region annotations")
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    labels = sorted(df_plot['region_label'].unique())
    positions = range(len(labels))
    data = [df_plot[df_plot['region_label'] == l]['Ka_Ks'].values for l in labels]

    bp = ax.boxplot(data, positions=positions, widths=0.6, patch_artist=True,
                    showfliers=False)

    colors = ['#4C72B0', '#C44E52', '#55A868', '#8172B3', '#CCB974']
    for patch, color in zip(bp['boxes'], colors[:len(labels)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Overlay jittered points
    for i, d in enumerate(data):
        jitter = np.random.normal(0, 0.05, size=len(d))
        ax.scatter(np.full(len(d), i) + jitter, d,
                   alpha=0.3, s=10, color='black', zorder=3)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=12, fontweight='bold')
    ax.set_ylabel('Ka/Ks', fontsize=14, fontweight='bold')
    ax.set_title('Ka/Ks by Genomic Region', fontsize=16, fontweight='bold')
    ax.tick_params(axis='y', labelsize=12)
    ax.axhline(1.0, color='red', linestyle='--', linewidth=1, alpha=0.5)

    # Add n per group
    for i, (l, d) in enumerate(zip(labels, data)):
        ax.text(i, ax.get_ylim()[1] * 0.95, f'n={len(d)}',
                ha='center', fontsize=10)

    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved regional Ka/Ks plot: {outpath}")


# ===========================================================================
# Divergence time estimation
# ===========================================================================

def estimate_divergence_time(median_ks, rate=1.7e-8):
    """
    Estimate divergence time from median Ks.

    T = Ks / (2 * rate)

    Parameters
    ----------
    median_ks : float
    rate : float
        Synonymous substitution rate per site per year.
        Default 1.7e-8 is from C. elegans calibrations.

    Returns
    -------
    float : estimated divergence time in years.
    """
    t = median_ks / (2 * rate)
    return t


# ===========================================================================
# Main
# ===========================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Ka/Ks pipeline for ortholog pairs between two species.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python kaks_pipeline.py \\
    --pairs orthologs.tsv \\
    --cds1 ascaris_cds.fa --prot1 ascaris_prot.fa \\
    --cds2 parascaris_cds.fa --prot2 parascaris_prot.fa \\
    --method KaKs_Calculator \\
    --threads 8 \\
    --outdir kaks_results

Pairs file format (tab-separated, no header):
  ascaris_gene_001    parascaris_gene_001
  ascaris_gene_002    parascaris_gene_002
  ...
        """
    )

    # Required inputs
    parser.add_argument('--pairs', required=True,
                        help='TSV file: species1_gene_id <tab> species2_gene_id')
    parser.add_argument('--cds1', required=True,
                        help='CDS FASTA for species 1 (e.g., Ascaris)')
    parser.add_argument('--prot1', required=True,
                        help='Protein FASTA for species 1')
    parser.add_argument('--cds2', required=True,
                        help='CDS FASTA for species 2 (e.g., Parascaris)')
    parser.add_argument('--prot2', required=True,
                        help='Protein FASTA for species 2')

    # Method
    parser.add_argument('--method', default='KaKs_Calculator',
                        choices=['KaKs_Calculator', 'yn00'],
                        help='Ka/Ks calculation method (default: KaKs_Calculator)')

    # Filtering
    parser.add_argument('--min-align-len', type=int, default=100,
                        help='Minimum codon alignment length (default: 100)')
    parser.add_argument('--max-ks', type=float, default=5.0,
                        help='Maximum Ks to include (default: 5.0)')
    parser.add_argument('--max-kaks', type=float, default=50.0,
                        help='Maximum Ka/Ks to include (default: 50.0)')

    # Species names for plot labels
    parser.add_argument('--species1', default='Ascaris',
                        help='Species 1 name for plot labels')
    parser.add_argument('--species2', default='Parascaris',
                        help='Species 2 name for plot labels')

    # Optional region annotation
    parser.add_argument('--regions', default=None,
                        help='TSV: gene_id <tab> region_label (e.g., CBR_adjacent)')

    # Substitution rate for divergence time
    parser.add_argument('--sub-rate', type=float, default=1.7e-8,
                        help='Synonymous substitution rate/site/year (default: 1.7e-8)')

    # Output
    parser.add_argument('--outdir', default='kaks_results',
                        help='Output directory (default: kaks_results)')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of parallel threads (default: 1)')
    parser.add_argument('--keep-tmp', action='store_true',
                        help='Keep temporary alignment files')

    return parser.parse_args()


def main():
    args = parse_args()

    # --- Setup output directory ---
    os.makedirs(args.outdir, exist_ok=True)

    # --- Logging ---
    log_path = os.path.join(args.outdir, 'kaks_pipeline.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("=" * 60)
    logging.info("Ka/Ks Pipeline")
    logging.info("=" * 60)
    logging.info(f"Arguments: {vars(args)}")

    # --- Check external tool availability ---
    tool_name = args.method if args.method == 'yn00' else 'KaKs_Calculator'
    if shutil.which(tool_name) is None:
        sys.exit(f"ERROR: {tool_name} not found in PATH. "
                 f"Install it or add to PATH before running.")
    if shutil.which('mafft') is None:
        sys.exit("ERROR: mafft not found in PATH. Install MAFFT first.")
    logging.info(f"External tools verified: mafft, {tool_name}")

    # --- Load sequences ---
    prot_db1 = load_fasta(args.prot1)
    prot_db2 = load_fasta(args.prot2)
    cds_db1 = load_fasta(args.cds1)
    cds_db2 = load_fasta(args.cds2)

    # --- Load ortholog pairs ---
    logging.info(f"Loading ortholog pairs: {args.pairs}")
    pairs = []
    with open(args.pairs, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                parts = line.split()  # fallback to whitespace
            if len(parts) >= 2:
                pairs.append((parts[0], parts[1]))
    logging.info(f"  Loaded {len(pairs)} ortholog pairs")

    if len(pairs) == 0:
        sys.exit("ERROR: No ortholog pairs loaded. Check file format.")

    # --- Run pipeline ---
    tmpdir = tempfile.mkdtemp(prefix='kaks_tmp_')
    logging.info(f"Temporary directory: {tmpdir}")

    try:
        df_raw = run_pipeline(pairs, prot_db1, prot_db2, cds_db1, cds_db2,
                              args.method, args.threads, tmpdir)

        # --- Save raw results ---
        raw_path = os.path.join(args.outdir, 'kaks_raw.tsv')
        df_raw.to_csv(raw_path, sep='\t', index=False)
        logging.info(f"Raw results saved: {raw_path}")

        # --- Status summary ---
        status_counts = df_raw['status'].value_counts()
        logging.info("Status summary:")
        for status, count in status_counts.items():
            logging.info(f"  {status}: {count}")

        # --- Filter ---
        df_filt = filter_results(df_raw,
                                 min_align_len=args.min_align_len,
                                 max_ks=args.max_ks,
                                 max_kaks=args.max_kaks)

        filt_path = os.path.join(args.outdir, 'kaks_filtered.tsv')
        df_filt.to_csv(filt_path, sep='\t', index=False)
        logging.info(f"Filtered results saved: {filt_path}")

        if len(df_filt) == 0:
            logging.error("No pairs passed filtering. Check input data quality.")
            return

        # --- Plots ---
        # SVG for Illustrator compatibility
        plot_ks_distribution(
            df_filt,
            os.path.join(args.outdir, 'ks_distribution.svg'),
            species1=args.species1, species2=args.species2
        )
        plot_kaks_distribution(
            df_filt,
            os.path.join(args.outdir, 'kaks_distribution.svg'),
            species1=args.species1, species2=args.species2
        )
        plot_ka_vs_ks(
            df_filt,
            os.path.join(args.outdir, 'ka_vs_ks.svg'),
            species1=args.species1, species2=args.species2
        )

        # Optional: regional comparison
        if args.regions:
            plot_kaks_by_region(
                df_filt, args.regions,
                os.path.join(args.outdir, 'kaks_by_region.svg')
            )

        # --- Divergence time estimate ---
        median_ks = df_filt['Ks'].median()
        t_years = estimate_divergence_time(median_ks, rate=args.sub_rate)
        t_mya = t_years / 1e6
        logging.info(f"Divergence time estimate:")
        logging.info(f"  Median Ks = {median_ks:.4f}")
        logging.info(f"  Substitution rate = {args.sub_rate:.2e} /syn site/year")
        logging.info(f"  T = Ks / (2 * rate) = {t_years:.0f} years = {t_mya:.1f} MYA")
        logging.info(f"  NOTE: Rate from C. elegans; clade III nematodes may differ.")

        # Save summary
        summary_path = os.path.join(args.outdir, 'summary.txt')
        with open(summary_path, 'w') as fh:
            fh.write("Ka/Ks Pipeline Summary\n")
            fh.write("=" * 40 + "\n")
            fh.write(f"Total ortholog pairs: {len(pairs)}\n")
            fh.write(f"Successfully calculated: {(df_raw['status'] == 'OK').sum()}\n")
            fh.write(f"Passed filtering: {len(df_filt)}\n")
            fh.write(f"Method: {args.method}\n")
            fh.write(f"\nKs statistics:\n")
            fh.write(f"  Median: {df_filt['Ks'].median():.4f}\n")
            fh.write(f"  Mean: {df_filt['Ks'].mean():.4f}\n")
            fh.write(f"  Std: {df_filt['Ks'].std():.4f}\n")
            fh.write(f"\nKa/Ks statistics:\n")
            fh.write(f"  Median: {df_filt['Ka_Ks'].median():.4f}\n")
            fh.write(f"  Mean: {df_filt['Ka_Ks'].mean():.4f}\n")
            fh.write(f"  Genes under purifying selection (Ka/Ks < 1): "
                     f"{(df_filt['Ka_Ks'] < 1).sum()} "
                     f"({(df_filt['Ka_Ks'] < 1).mean()*100:.1f}%)\n")
            fh.write(f"  Genes under positive selection (Ka/Ks > 1): "
                     f"{(df_filt['Ka_Ks'] > 1).sum()} "
                     f"({(df_filt['Ka_Ks'] > 1).mean()*100:.1f}%)\n")
            fh.write(f"\nDivergence time estimate:\n")
            fh.write(f"  Rate: {args.sub_rate:.2e} /syn site/year\n")
            fh.write(f"  T = {t_mya:.1f} MYA\n")
            fh.write(f"  Caution: rate from C. elegans; actual rate may differ\n")

        logging.info(f"Summary saved: {summary_path}")
        logging.info("Pipeline complete.")

    finally:
        if not args.keep_tmp:
            shutil.rmtree(tmpdir, ignore_errors=True)
            logging.info("Temporary files cleaned up")
        else:
            logging.info(f"Temporary files kept at: {tmpdir}")


if __name__ == '__main__':
    main()
