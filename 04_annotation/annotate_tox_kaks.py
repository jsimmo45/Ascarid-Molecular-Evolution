#!/usr/bin/env python3
"""
annotate_tox_kaks.py — Annotate Toxocara Ka/Ks results with eliminated/retained labels.

Labels Toxocara genes by set membership against the eliminated transcript list
(no GFF or BED intersection needed — just ID matching). Optionally labels the
partner species (Ascaris/Parascaris) genes by GFF + retained/eliminated BED
intersection, reusing the same logic as annotate_kaks_regions.py.

Compares Ka/Ks, Ka, Ks, and dN/dS between retained and eliminated genes for
both species sides (where annotations are available).

Usage:
  python annotate_tox_kaks.py \
    --kaks kaks_Tox_As/kaks_filtered.tsv \
    --tox-eliminated /path/to/eliminated.transcripts.longest.bed \
    --species1 Toxocara --species2 Ascaris \
    --gff2 As_miniprot.gff \
    --retained2 AG.v50.retained_loose.bed \
    --eliminated2 AG.v50.eliminated_strict.bed \
    --outdir kaks_Tox_As_annotated

Dependencies: pandas, matplotlib, scipy
"""

import argparse
import logging
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
from scipy.stats import mannwhitneyu


# ===========================================================================
# Toxocara gene labeling (set membership)
# ===========================================================================

def load_tox_eliminated_ids(bed_path):
    """
    Extract eliminated transcript IDs from Toxocara eliminated BED.

    The BED12 file has MSTRG.xxx.y in column 4 (name field).
    Returns a set of transcript IDs (without .pN suffix).

    Parameters
    ----------
    bed_path : str
        Path to eliminated.transcripts.longest.bed

    Returns
    -------
    set : eliminated transcript IDs (e.g., {'MSTRG.524.3', 'MSTRG.909.5', ...})
    """
    logging.info(f"Loading Toxocara eliminated IDs from: {bed_path}")
    elim_ids = set()
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                elim_ids.add(parts[3])
    logging.info(f"  Loaded {len(elim_ids)} eliminated transcript IDs")
    return elim_ids


def strip_transdecoder_suffix(gene_id):
    """
    Strip TransDecoder .pN suffix from a gene ID.

    MSTRG.524.3.p1 -> MSTRG.524.3
    MP001234       -> MP001234  (no suffix, unchanged)

    Parameters
    ----------
    gene_id : str

    Returns
    -------
    str : transcript-level ID
    """
    return re.sub(r'\.p\d+$', '', gene_id)


def label_tox_genes(kaks_df, elim_ids, gene_col='gene1'):
    """
    Label Toxocara genes as 'retained' or 'eliminated'.

    Strips .pN suffix from TransDecoder IDs in the Ka/Ks results,
    then checks membership in the eliminated set.

    Parameters
    ----------
    kaks_df : DataFrame
        Ka/Ks filtered results.
    elim_ids : set
        Eliminated transcript IDs (no .pN suffix).
    gene_col : str
        Column containing Toxocara gene IDs.

    Returns
    -------
    Series : 'eliminated' or 'retained' for each row.
    """
    tx_ids = kaks_df[gene_col].apply(strip_transdecoder_suffix)
    labels = tx_ids.apply(lambda x: 'eliminated' if x in elim_ids else 'retained')
    n_elim = (labels == 'eliminated').sum()
    n_ret = (labels == 'retained').sum()
    logging.info(f"  Toxocara labels: {n_ret} retained, {n_elim} eliminated")
    return labels


# ===========================================================================
# Partner species labeling (GFF + BED intersection)
# ===========================================================================

def parse_gff_genes(gff_path):
    """
    Extract gene/mRNA positions from miniprot GFF.
    Returns dict: {gene_id: (chrom, start, end, midpoint)}.
    """
    logging.info(f"Parsing gene positions from: {gff_path}")
    genes = {}
    with open(gff_path, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'mRNA':
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]
            gene_id = None
            for field in attrs.split(';'):
                if field.startswith('ID='):
                    gene_id = field.split('=')[1]
                    break
            if gene_id is None:
                for field in attrs.split(';'):
                    if field.startswith('Target='):
                        gene_id = field.split('=')[1].split()[0]
                        break
            if gene_id and gene_id not in genes:
                genes[gene_id] = (chrom, start, end, (start + end) // 2)
    logging.info(f"  Extracted {len(genes)} gene positions")
    return genes


def load_bed_regions(bed_path):
    """Load BED file as list of (chrom, start, end) tuples."""
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions


def point_in_regions(chrom, midpoint, regions):
    """Check if a point falls within any region."""
    for r_chrom, r_start, r_end in regions:
        if chrom == r_chrom and r_start <= midpoint < r_end:
            return True
    return False


def label_partner_genes(kaks_df, gene_col, gff_path, retained_bed, eliminated_bed):
    """
    Label partner species genes as 'retained', 'eliminated', or 'ambiguous'
    using GFF positions + BED intersection.

    Parameters
    ----------
    kaks_df : DataFrame
    gene_col : str
        Column with partner gene IDs (e.g., 'gene2').
    gff_path : str
    retained_bed, eliminated_bed : str
        Paths to BED files.

    Returns
    -------
    Series : region label per row.
    """
    gene_pos = parse_gff_genes(gff_path)
    ret_regions = load_bed_regions(retained_bed)
    elim_regions = load_bed_regions(eliminated_bed)

    logging.info(f"  Retained regions: {len(ret_regions)}")
    logging.info(f"  Eliminated regions: {len(elim_regions)}")

    labels = []
    for gid in kaks_df[gene_col]:
        if gid not in gene_pos:
            labels.append('unannotated')
            continue
        chrom, start, end, mid = gene_pos[gid]
        in_ret = point_in_regions(chrom, mid, ret_regions)
        in_elim = point_in_regions(chrom, mid, elim_regions)
        if in_elim and not in_ret:
            labels.append('eliminated')
        elif in_ret and not in_elim:
            labels.append('retained')
        elif in_elim and in_ret:
            labels.append('ambiguous')
        else:
            labels.append('unannotated')

    labels = pd.Series(labels, index=kaks_df.index)
    for r in ['retained', 'eliminated', 'ambiguous', 'unannotated']:
        n = (labels == r).sum()
        if n > 0:
            logging.info(f"  Partner {r}: {n}")
    return labels


# ===========================================================================
# Comparison and plotting
# ===========================================================================

def compare_retained_vs_eliminated(df, region_col, metrics, label):
    """
    Mann-Whitney U comparison of retained vs eliminated genes.

    Parameters
    ----------
    df : DataFrame
        Annotated Ka/Ks results.
    region_col : str
        Column with 'retained'/'eliminated' labels.
    metrics : list of str
        Columns to compare (e.g., ['Ka_Ks', 'Ka', 'Ks']).
    label : str
        Label for logging (e.g., 'Toxocara' or 'Ascaris').

    Returns
    -------
    DataFrame : summary statistics.
    """
    ret = df[df[region_col] == 'retained']
    elim = df[df[region_col] == 'eliminated']

    logging.info(f"\n{'='*60}")
    logging.info(f"Comparison: {label} retained (n={len(ret)}) vs eliminated (n={len(elim)})")
    logging.info(f"{'='*60}")

    rows = []
    for metric in metrics:
        r_vals = ret[metric].dropna()
        e_vals = elim[metric].dropna()

        row = {
            'metric': metric,
            'ret_n': len(r_vals),
            'ret_median': r_vals.median() if len(r_vals) > 0 else np.nan,
            'ret_mean': r_vals.mean() if len(r_vals) > 0 else np.nan,
            'elim_n': len(e_vals),
            'elim_median': e_vals.median() if len(e_vals) > 0 else np.nan,
            'elim_mean': e_vals.mean() if len(e_vals) > 0 else np.nan,
            'fold_change': np.nan,
            'U_stat': np.nan,
            'p_value': np.nan
        }

        if len(r_vals) >= 5 and len(e_vals) >= 5:
            stat, pval = mannwhitneyu(r_vals, e_vals, alternative='two-sided')
            row['U_stat'] = stat
            row['p_value'] = pval
            if row['ret_median'] > 0:
                row['fold_change'] = row['elim_median'] / row['ret_median']

        rows.append(row)
        logging.info(f"  {metric}: ret median={row['ret_median']:.4f}, "
                     f"elim median={row['elim_median']:.4f}, "
                     f"fold={row['fold_change']:.2f}x, p={row['p_value']:.2e}")

    return pd.DataFrame(rows)


def plot_boxplots(df, region_col, outpath, species_label, partner_label):
    """
    Publication-quality boxplots of Ka/Ks, Ka, Ks for retained vs eliminated.
    Retained = blue/left, eliminated = red/right.
    """
    metrics = ['Ka_Ks', 'Ka', 'Ks']
    titles = ['Ka/Ks (dN/dS)', 'Ka (dN)', 'Ks (dS)']
    ylims = [1.0, 0.5, None]

    # Keep only metrics present in the data
    keep = [(m, t, y) for m, t, y in zip(metrics, titles, ylims) if m in df.columns]
    if not keep:
        logging.warning("No plottable metric columns; skipping boxplot.")
        return
    metrics, titles, ylims = zip(*keep)
    metrics, titles, ylims = list(metrics), list(titles), list(ylims)

    fig, axes = plt.subplots(1, len(metrics), figsize=(4.7 * len(metrics), 5))
    if len(metrics) == 1:
        axes = [axes]

    for ax, metric, title, ylim in zip(axes, metrics, titles, ylims):
        ret_vals = df[df[region_col] == 'retained'][metric].dropna()
        elim_vals = df[df[region_col] == 'eliminated'][metric].dropna()

        data = [ret_vals, elim_vals]
        colors = ['#4C72B0', '#C44E52']  # blue=retained, red=eliminated
        labels_tick = [
            f'Retained\n(n={len(ret_vals):,})',
            f'Eliminated\n(n={len(elim_vals):,})'
        ]

        bp = ax.boxplot(data, positions=[0, 1], widths=0.5,
                        patch_artist=True, showfliers=False, zorder=2)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        for median in bp['medians']:
            median.set_color('black')
            median.set_linewidth(2)

        # Jitter
        for i, (vals, color) in enumerate(zip(data, colors)):
            if len(vals) > 2000:
                idx = np.random.choice(len(vals), 2000, replace=False)
                v_plot = vals.values[idx]
            else:
                v_plot = vals.values
            jitter = np.random.normal(0, 0.06, size=len(v_plot))
            ax.scatter(np.full(len(v_plot), i) + jitter, v_plot,
                       s=4, alpha=0.15, color=color,
                       edgecolors='none', zorder=3)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(labels_tick, fontsize=12, fontweight='bold')
        ax.set_ylabel(title, fontsize=14, fontweight='bold')
        ax.tick_params(axis='y', labelsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if ylim is not None:
            ax.set_ylim(0, ylim)

        # Mann-Whitney p-value
        if len(ret_vals) >= 5 and len(elim_vals) >= 5:
            _, pval = mannwhitneyu(ret_vals, elim_vals, alternative='two-sided')
            pstr = f'p = {pval:.2e}'
            ax.text(0.5, 0.95, pstr, transform=ax.transAxes,
                    ha='center', va='top', fontsize=11, fontstyle='italic')

    fig.suptitle(f'{species_label} genes: {species_label}×{partner_label}',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    fig.savefig(outpath.replace('.pdf', '.svg'), bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved: {outpath}")


# ===========================================================================
# Main
# ===========================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate Toxocara Ka/Ks results with eliminated/retained labels.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example (Tox x Ascaris with partner annotation):
  python annotate_tox_kaks.py \\
    --kaks kaks_Tox_As/kaks_filtered.tsv \\
    --tox-eliminated /path/to/eliminated.transcripts.longest.bed \\
    --species1 Toxocara --species2 Ascaris \\
    --gff2 As_miniprot.gff \\
    --retained2 AG.v50.retained_loose.bed \\
    --eliminated2 AG.v50.eliminated_strict.bed \\
    --outdir kaks_Tox_As_annotated

Example (Tox x Baylisascaris, no partner annotation):
  python annotate_tox_kaks.py \\
    --kaks kaks_Tox_Bsc/kaks_filtered.tsv \\
    --tox-eliminated /path/to/eliminated.transcripts.longest.bed \\
    --species1 Toxocara --species2 Baylisascaris \\
    --outdir kaks_Tox_Bsc_annotated
        """
    )

    # Required
    parser.add_argument('--kaks', required=True,
                        help='Ka/Ks filtered TSV (gene1=Toxocara, gene2=partner)')
    parser.add_argument('--tox-eliminated', required=True,
                        help='Toxocara eliminated.transcripts.longest.bed')

    # Species labels
    parser.add_argument('--species1', default='Toxocara',
                        help='Species 1 name (default: Toxocara)')
    parser.add_argument('--species2', default='Partner',
                        help='Species 2 name (default: Partner)')

    # Optional: partner species annotation
    parser.add_argument('--gff2', default=None,
                        help='Miniprot GFF for partner species (optional)')
    parser.add_argument('--retained2', default=None,
                        help='Retained regions BED for partner species (optional)')
    parser.add_argument('--eliminated2', default=None,
                        help='Eliminated regions BED for partner species (optional)')

    parser.add_argument('--outdir', default='kaks_tox_annotated')
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Logging
    log_path = os.path.join(args.outdir, 'annotate_tox.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("=" * 60)
    logging.info("Toxocara Ka/Ks Annotation")
    logging.info("=" * 60)
    logging.info(f"Arguments: {vars(args)}")

    # --- Load Ka/Ks results ---
    kaks_df = pd.read_csv(args.kaks, sep='\t')
    logging.info(f"Loaded {len(kaks_df)} Ka/Ks pairs")
    logging.info(f"  gene1 examples (Tox): {list(kaks_df['gene1'].head(3))}")
    logging.info(f"  gene2 examples ({args.species2}): {list(kaks_df['gene2'].head(3))}")

    # --- Label Toxocara genes ---
    elim_ids = load_tox_eliminated_ids(args.tox_eliminated)
    kaks_df['tox_region'] = label_tox_genes(kaks_df, elim_ids, gene_col='gene1')
    kaks_df['tox_transcript'] = kaks_df['gene1'].apply(strip_transdecoder_suffix)

    # --- Optionally label partner species genes ---
    has_partner = (args.gff2 is not None and
                   args.retained2 is not None and
                   args.eliminated2 is not None)

    if has_partner:
        logging.info(f"\nAnnotating {args.species2} genes by GFF + BED")
        kaks_df['partner_region'] = label_partner_genes(
            kaks_df, 'gene2', args.gff2, args.retained2, args.eliminated2
        )

        # Consensus requires BOTH sides annotated. Toxocara is always
        # retained/eliminated (set membership), so the partner side is the
        # binding constraint: 'concordant' only when both agree, 'discordant'
        # when they disagree (lineage-specific elimination), 'partner_NA' when
        # the partner gene has no retained/eliminated assignment.
        def consensus(row):
            t = row['tox_region']
            p = row['partner_region']
            if p not in ('retained', 'eliminated'):
                return 'partner_NA'
            return t if t == p else 'discordant'
        kaks_df['consensus_region'] = kaks_df.apply(consensus, axis=1)

    # --- Save annotated table ---
    out_tsv = os.path.join(args.outdir, 'kaks_annotated.tsv')
    kaks_df.to_csv(out_tsv, sep='\t', index=False)
    logging.info(f"\nAnnotated table saved: {out_tsv}")

    # --- Summary counts ---
    logging.info(f"\nToxocara region counts:")
    for r in ['retained', 'eliminated']:
        logging.info(f"  {r}: {(kaks_df['tox_region'] == r).sum()}")

    if has_partner:
        logging.info(f"\n{args.species2} region counts:")
        for r in ['retained', 'eliminated', 'ambiguous', 'unannotated']:
            n = (kaks_df['partner_region'] == r).sum()
            if n > 0:
                logging.info(f"  {r}: {n}")

        logging.info(f"\nConsensus region counts:")
        for r in ['retained', 'eliminated', 'discordant', 'partner_NA']:
            n = (kaks_df['consensus_region'] == r).sum()
            if n > 0:
                logging.info(f"  {r}: {n}")
        n_disc = (kaks_df['consensus_region'] == 'discordant').sum()
        n_both = kaks_df['consensus_region'].isin(
            ['retained', 'eliminated', 'discordant']).sum()
        if n_both > 0:
            logging.info(f"  Lineage-specific (discordant) fraction: "
                         f"{n_disc}/{n_both} = {100*n_disc/n_both:.1f}% of "
                         f"genes with both sides annotated")

    # --- Comparisons ---
    # Use only metric columns that actually exist in the filtered TSV.
    candidate_metrics = ['Ka_Ks', 'Ka', 'Ks']
    metrics = [m for m in candidate_metrics if m in kaks_df.columns]
    missing = [m for m in candidate_metrics if m not in kaks_df.columns]
    if missing:
        logging.warning(f"Metric column(s) not found, skipping: {missing}. "
                        f"Available columns: {list(kaks_df.columns)}")
    if not metrics:
        sys.exit("ERROR: none of Ka_Ks/Ka/Ks present in the Ka/Ks table.")
    all_stats = []

    # 1. Toxocara-side comparison
    stats_tox = compare_retained_vs_eliminated(
        kaks_df, 'tox_region', metrics,
        f'{args.species1} (labeled by {args.species1} elimination status)'
    )
    stats_tox['comparison'] = f'{args.species1}_region'
    all_stats.append(stats_tox)

    plot_boxplots(kaks_df, 'tox_region',
                  os.path.join(args.outdir, f'kaks_tox_retained_vs_eliminated.pdf'),
                  species_label=args.species1, partner_label=args.species2)

    # 2. Partner-side comparison (if available)
    if has_partner:
        df_partner = kaks_df[kaks_df['partner_region'].isin(['retained', 'eliminated'])]
        if len(df_partner) > 0:
            stats_partner = compare_retained_vs_eliminated(
                df_partner, 'partner_region', metrics,
                f'{args.species2} (labeled by {args.species2} elimination status)'
            )
            stats_partner['comparison'] = f'{args.species2}_region'
            all_stats.append(stats_partner)

            plot_boxplots(df_partner, 'partner_region',
                          os.path.join(args.outdir,
                                       f'kaks_partner_retained_vs_eliminated.pdf'),
                          species_label=args.species2, partner_label=args.species1)

        # 3. Consensus comparison
        df_cons = kaks_df[kaks_df['consensus_region'].isin(['retained', 'eliminated'])]
        if len(df_cons) > 0:
            stats_cons = compare_retained_vs_eliminated(
                df_cons, 'consensus_region', metrics,
                'Consensus (both species agree)'
            )
            stats_cons['comparison'] = 'consensus'
            all_stats.append(stats_cons)

    # --- Save all stats ---
    if all_stats:
        stats_df = pd.concat(all_stats, ignore_index=True)
        stats_path = os.path.join(args.outdir, 'comparison_stats.tsv')
        stats_df.to_csv(stats_path, sep='\t', index=False)
        logging.info(f"\nComparison stats saved: {stats_path}")

    # --- Summary file ---
    summary_path = os.path.join(args.outdir, 'summary.txt')
    with open(summary_path, 'w') as fh:
        fh.write(f"Toxocara Ka/Ks Annotation Summary\n")
        fh.write(f"{'='*50}\n\n")
        fh.write(f"Comparison: {args.species1} x {args.species2}\n")
        fh.write(f"Total pairs: {len(kaks_df)}\n\n")

        fh.write(f"{args.species1} gene labels:\n")
        for r in ['retained', 'eliminated']:
            n = (kaks_df['tox_region'] == r).sum()
            fh.write(f"  {r}: {n}\n")
        fh.write(f"\n")

        if all_stats:
            stats_df = pd.concat(all_stats, ignore_index=True)
            fh.write(f"Comparisons:\n")
            fh.write(stats_df.to_string(index=False))
            fh.write(f"\n")

    logging.info(f"Summary saved: {summary_path}")
    logging.info("Done.")


if __name__ == '__main__':
    main()
