#!/usr/bin/env python3
"""
branch_dnds_four.py -- Four-taxon branch-specific dN/dS using codeml free-ratio.

Identifies 1:1:1:1 orthologs across Ascaris, Parascaris, Baylisascaris and
Toxocara from the six pairwise RBH files, runs codeml's free-ratio model on the
UNROOTED tree

    ((Ascaris,Parascaris),Baylisascaris,Toxocara);

and extracts a per-branch dN, dS and dN/dS for every gene. Compares retained
versus eliminated genes on each terminal branch.

Produces Figure 5A, Figure 5B and Table T2 (plus the S14 versions when
--tox-eliminated is given), and the annotated table that
`branch_dnds_decomposition.py` reads for 5C and 5D.

Terminal branches are identified by child tip number -- PAML assigns 1..4 from
the order of the sequences in the alignment, written here As, Pu, Bsc, Tox --
which is robust to how it numbers internal nodes.


2026-10-10 -- THE TREE WAS ROOTED; EVERY TOXOCARA BRANCH VALUE CHANGES
----------------------------------------------------------------------
Until this date codeml was given `(((Ascaris,Parascaris),Baylisascaris),
Toxocara);` with clock = 0. codeml does not unroot that: it warns and fits both
branches at the root, of which only the sum is estimable, so the Toxocara
branch and the root-to-ascarid branch split one length arbitrarily. On the
2026-10-09 table the "internal" branch subtended Toxocara alone in all 5,954
genes; where Toxocara dN was 0 the other half carried the whole path in 1,268
of 1,272 genes; and the two halves summed to a dS median of 1.277, against
1.267 for the Toxocara branch of the Ks tree. The Ascaris, Parascaris and
Baylisascaris branches are identifiable either way and should barely move.

  * `branch_lib.TREE_TOPOLOGY` is now unrooted, main() refuses a rooted tree
    before starting any codeml job, and a branch table that does not have
    exactly five branches becomes a per-gene status instead of being read.
  * The internal branch now is the {Ascaris, Parascaris} | {Baylisascaris,
    Toxocara} split, as the section further down always said it was.
  * Figure 5A: dN/dS at codeml's floor of 0.0001 means dN = 0. Those genes are
    drawn in a band at 0 below the log axis, with their count beside each box,
    instead of sitting on the bottom spine uncounted (the old test was `<= 0`,
    which never fires). Whiskers are the 5th and 95th percentiles -- 1.5 IQR
    whiskers on a log axis ran to the floor for every box -- and the stars and
    fold sit above each panel instead of on the jitter.
  * The cap now also drops a branch with dS = 0 under the cap, i.e. dN = dS =
    0. With no substitution on the branch, dN/dS does not affect the likelihood
    and codeml reports wherever its randomly started optimiser stopped (0.0001
    in one test run, 0.79 or 3.48 in the next), so whether such a gene passed
    the cap was chance. Logged per branch as "at 0/0".
  * codeml seeds its starting values from the clock and has no seed option, so
    a few genes with flat or multimodal likelihoods differ slightly between
    runs (in the end-to-end test, 1 of 64 by about 2% outside the 0/0 genes).
    Medians are not exactly reproducible in the third decimal.


2026-10-14 -- TICK LABELS AND ZERO COUNTS (figure only; no number changes)
--------------------------------------------------------------------------
  * 5A's log axis is labelled 0.001 ... 10 instead of 10^-3 ... 10^1
    (`branch_lib.decade_label`). The power form was mathtext, which reached
    the assembled figure in DejaVu Sans with 4.9 pt exponents.
  * The counts beside the zero band are 7 pt (FONT['annotation']); they were
    6 pt.
  Use --replot on the existing annotated table; the stats are unchanged.


2026-10-10, LATER -- --replot, AND MEASURED MARGINS
---------------------------------------------------
  * `--replot TABLE` redraws Figure 5A, 5B, Table T2 and the summary from an
    existing `branch_dnds_four_annotated.tsv`, with no ortholog search and no
    codeml. A full rerun refits every quartet from codeml's random starts (see
    above), so it would move the quoted medians and P values to fix a figure.
    From the annotated table onward the replot runs the same code as a full
    run (`report()`), so the same table, cap and seed give identical stats and
    the same jitter. It warns when --omega-cap or --seed differ from the run
    that wrote the table (read from that run's provenance.json).
  * 5A and 5B: the bottom margin is measured (`branch_lib.fit_bottom_margin`,
    1 mm under the lowest label) instead of a fixed fraction, which left
    5-6 mm of empty figure below the labels.
  * The 5A notes and docstrings said "dN = 0" for dN/dS at codeml's floor.
    On a long branch they differ: dN prints as (dN/dS) x dS, so a Toxocara
    gene at the floor shows dN = 0.0001 (21 floor genes, 1 with dN = 0 in the
    2026-10-10 run). The 5A band is now described as what it is.


WHAT CHANGED IN THE 2026-08-18 CONVERSION
-----------------------------------------
Six substantive changes beyond the style pass. Three of them alter a number.

1. FIGURE SIZE. Built at exactly the journal column width with sizes in mm.
   The old 20 x 6 inch canvas saved with bbox_inches='tight' is a 2.7x
   reduction when placed at 185 mm, so its 14 pt titles printed under 5.2 pt --
   below the 6 pt journal floor. Same defect that was found in the genome
   traces.

2. LOG AXIS ON 5A. The old panel autoscaled to the omega cap of 5 because the
   jitter layer draws every point up to the cap, which squashed all four boxes
   into the bottom ~8% of the axis. Branch omega spans two orders of magnitude
   and belongs on a log axis. The BOX STATISTICS ARE COMPUTED ON THE FULL SET
   so the drawn median matches the table. (How dN = 0 is drawn changed on
   2026-10-10; see above.)

3. SEM ON A MEDIAN, in the 5B bar panel. `Series.sem()` is the standard error
   of the MEAN. On a right-skewed distribution it is not even centred on the
   median it was drawn against. Replaced with a percentile bootstrap CI.

4. NEUTRALITY LINE REMOVED FROM 5B. Drawing omega = 1 on a panel whose tallest
   bar is ~0.4 spends 60% of the axis on empty space to make a point the data
   never approaches. The summary states the maximum branch median instead.
   The line stays on 5A, where a log axis makes it free.

5. REGION CALLS GO THROUGH `region_lib`. The old `annotate_gene_region` looped
   every gene against every BED interval with `iterrows()` AND tested
   `start <= mid <= end` on raw BED coordinates. BED is 0-based half-open, so
   that interval is one base too wide at its left edge. The error only bites a
   gene whose midpoint lands exactly on a BED start, but it meant this script
   and `annotate_regions.py` did not agree on what a region is.

6. BH CORRECTION ON T2. Eight tests are run (four branches x two region
   labellings) and the table reported raw P for all eight. q-values are now
   computed within each region labelling, matching the convention in
   `branch_dnds_decomposition.py` -- the two scripts previously disagreed about
   whether these tests needed correction at all.

Also: seeded RNG so jitter is reproducible, `make_outdir()` so a rerun cannot
overwrite the run a figure came from, and provenance.json.


2026-10-09
----------
Each worker task now carries only its quartet's eight sequence records
(`quartet_payload`). Passing the full FASTA dictionaries made the parent pickle
about 150 MB per task, serially; across 6,258 quartets that, not codeml, was
most of the v2 run's 1.5 h. No result changes.


THE INTERNAL BRANCH IS NOT AN ANCESTRAL BRANCH
----------------------------------------------
codeml runs with clock = 0 on the tree written UNROOTED (it does not unroot a
rooted string itself -- see 2026-10-10 above), and an unrooted four-taxon tree
has exactly one internal branch. It separates {Ascaris, Parascaris} from
{Baylisascaris, Toxocara}; it is not the branch leading to the ascarid clade,
and no rooting of this tree makes it so. The old column `internal_omega_min`
was a minimum over a list that always had one element, which invited that
reading. It is now `internal_omega` with an explicit `internal_subtends`
column, and the tip set is derived from codeml's own branch list rather than
assumed.

Usage
-----
  python branch_dnds_four.py \
    --rbh-as-pu rbh_results/orthologs.tsv \
    --rbh-as-bsc rbh_As_Bsc_v2/orthologs.tsv \
    --rbh-pu-bsc rbh_Pu_Bsc_v2/orthologs.tsv \
    --rbh-tox-as rbh_Tox_As/orthologs.tsv \
    --rbh-tox-pu rbh_Tox_Pu/orthologs.tsv \
    --rbh-tox-bsc rbh_Tox_Bsc_v2/orthologs.tsv \
    --cds-as As_cds.fa --prot-as As_prot.fa \
    ... \
    --retained1 AG.v50.retained_loose.bed \
    --eliminated1 AG.v50.eliminated_strict.bed \
    --gff1 As_miniprot.gff \
    --tox-eliminated eliminated.transcripts.longest.bed \
    --threads 48 --outdir results/branch_dnds_four_<stamp>

  # figures, T2 and summary again from an existing run -- no codeml
  python branch_dnds_four.py \
    --replot results/branch_dnds_four_<stamp>/branch_dnds_four_annotated.tsv \
    --outdir results/branch_dnds_four_replot_<stamp>

Dependencies: biopython, mafft, codeml (PAML), pandas, numpy, scipy, matplotlib
"""

import argparse
import json
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from io import StringIO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

from Bio import SeqIO

from paper_colors import (REGION, NEUTRAL, COL, apply_style, audit, save,
                          make_outdir, mm, SPECIES_LABEL, FONT)
from region_lib import (parse_gff_models, select_primary_alignment, parse_bed,
                        merge_intervals, in_intervals, benjamini_hochberg,
                        fmt_p, stars, fingerprint, git_commit, setup_logging)
from branch_lib import (TREE_TOPOLOGY, SPECIES_ORDER, SPECIES_KEY, METRICS,
                        METRIC_LABEL, parse_branch_table, label_branches,
                        check_unrooted, internal_branch_note, usable_mask,
                        split_and_compare, bootstrap_median_ci, at_zero,
                        as_plotted, zero_band_limits, zero_band_axis,
                        mark_axis_break, fit_bottom_margin, load_annotated)

_HERE = os.path.dirname(os.path.abspath(__file__))

CODEML_CTL_TEMPLATE = """\
      seqfile = {seqfile}
      treefile = {treefile}
      outfile = {outfile}
      noisy = 0
      verbose = 0
      runmode = 0
      seqtype = 1
      CodonFreq = 2
      model = 1
      NSsites = 0
      fix_kappa = 0
      kappa = 2
      fix_omega = 0
      omega = 0.4
      clock = 0
      getSE = 0
      RateAncestor = 0
      Small_Diff = .5e-6
"""


# ===========================================================================
# Four-way ortholog identification
# ===========================================================================

def load_pairs(path, invert=False):
    """Load an RBH orthologs.tsv (gene1<TAB>gene2) into a 1:1 dict.

    `invert` keys by gene2 instead, which normalises the Toxocara files (written
    Tox-first) into the same orientation as the anchor species.
    """
    d = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            d[b if invert else a] = a if invert else b
    return d


def find_four_way_orthologs(rbh_as_pu, rbh_as_bsc, rbh_pu_bsc,
                            rbh_tox_as, rbh_tox_pu, rbh_tox_bsc):
    """Quartets in which all six pairwise RBH relationships are consistent.

    Returns a list of (as_gene, pu_gene, bsc_gene, tox_gene) plus a dict of
    diagnostics recording how many quartets each consistency check removed --
    a large loss at one step points at an identifier mismatch in that RBH file
    rather than at biology.
    """
    logging.info('Finding four-way 1:1:1:1 orthologs...')

    d_as_pu = load_pairs(rbh_as_pu)                    # As  -> Pu
    d_as_bsc = load_pairs(rbh_as_bsc)                  # As  -> Bsc
    d_pu_bsc = load_pairs(rbh_pu_bsc)                  # Pu  -> Bsc
    d_as_tox = load_pairs(rbh_tox_as, invert=True)     # As  -> Tox
    d_pu_tox = load_pairs(rbh_tox_pu, invert=True)     # Pu  -> Tox
    d_bsc_tox = load_pairs(rbh_tox_bsc, invert=True)   # Bsc -> Tox

    logging.info('  pairwise sizes: As-Pu=%s As-Bsc=%s Pu-Bsc=%s '
                 'As-Tox=%s Pu-Tox=%s Bsc-Tox=%s',
                 f'{len(d_as_pu):,}', f'{len(d_as_bsc):,}', f'{len(d_pu_bsc):,}',
                 f'{len(d_as_tox):,}', f'{len(d_pu_tox):,}', f'{len(d_bsc_tox):,}')

    diag = {'anchor_as_pu': len(d_as_pu), 'no_bsc': 0, 'no_tox': 0,
            'pu_bsc_inconsistent': 0, 'pu_tox_inconsistent': 0,
            'bsc_tox_inconsistent': 0}

    quartets = []
    for as_gene, pu_gene in d_as_pu.items():
        bsc_gene = d_as_bsc.get(as_gene)
        if bsc_gene is None:
            diag['no_bsc'] += 1
            continue
        tox_gene = d_as_tox.get(as_gene)
        if tox_gene is None:
            diag['no_tox'] += 1
            continue
        if d_pu_bsc.get(pu_gene) != bsc_gene:
            diag['pu_bsc_inconsistent'] += 1
            continue
        if d_pu_tox.get(pu_gene) != tox_gene:
            diag['pu_tox_inconsistent'] += 1
            continue
        if d_bsc_tox.get(bsc_gene) != tox_gene:
            diag['bsc_tox_inconsistent'] += 1
            continue
        quartets.append((as_gene, pu_gene, bsc_gene, tox_gene))

    for k, v in diag.items():
        if k != 'anchor_as_pu':
            logging.info('    dropped at %-22s %s', k + ':', f'{v:,}')
    logging.info('  four-way orthologs: %s', f'{len(quartets):,}')
    return quartets, diag


# ===========================================================================
# Sequence helpers
# ===========================================================================

def load_fasta(path):
    """FASTA into {id: SeqRecord}."""
    return {rec.id: rec for rec in SeqIO.parse(path, 'fasta')}


def align_proteins(prot_recs, tmpdir, name):
    """Align N protein sequences with MAFFT. Returns {id: aligned_seq} or None."""
    infile = os.path.join(tmpdir, f'{name}.prot.fa')
    SeqIO.write(prot_recs, infile, 'fasta')
    try:
        result = subprocess.run(['mafft', '--auto', '--quiet', infile],
                                capture_output=True, text=True, timeout=180)
        if result.returncode != 0:
            return None
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None

    aligned = OrderedDict()
    for rec in SeqIO.parse(StringIO(result.stdout), 'fasta'):
        aligned[rec.id] = str(rec.seq)
    return aligned if len(aligned) == len(prot_recs) else None


def back_translate(aligned_proteins, cds_records):
    """Back-translate a protein alignment to a codon alignment."""
    codon_aln = OrderedDict()
    for gene_id, aln_prot in aligned_proteins.items():
        if gene_id not in cds_records:
            return None
        cds_seq = str(cds_records[gene_id].seq).upper()
        if cds_seq[-3:] in ('TAA', 'TAG', 'TGA'):
            cds_seq = cds_seq[:-3]

        aln_prot_clean = aln_prot.rstrip('*')
        if len(aln_prot_clean.replace('-', '')) != len(cds_seq) // 3:
            return None

        codons, pos = [], 0
        for aa in aln_prot_clean:
            if aa == '-':
                codons.append('---')
            else:
                codon = cds_seq[pos * 3: pos * 3 + 3]
                if len(codon) != 3:
                    return None
                codons.append(codon)
                pos += 1
        codon_aln[gene_id] = ''.join(codons)

    if len({len(v) for v in codon_aln.values()}) != 1:
        return None
    return codon_aln


# ===========================================================================
# codeml free-ratio
# ===========================================================================

def _empty_result():
    r = {f'{sp}_{m}': np.nan for sp in SPECIES_ORDER for m in METRICS}
    r['internal_omega'] = np.nan
    r['internal_dN'] = np.nan
    r['internal_dS'] = np.nan
    r['internal_subtends'] = ''
    r['status'] = 'FAIL'
    return r


def run_codeml_free_ratio(codon_aln, species_for_id, name, tmpdir,
                          min_codons=50):
    """Run codeml free-ratio on a four-way codon alignment."""
    result = _empty_result()

    workdir = os.path.join(tmpdir, f'codeml_{name}')
    os.makedirs(workdir, exist_ok=True)

    ids = list(codon_aln.keys())
    seq_len = len(codon_aln[ids[0]])

    n_informative = sum(1 for i in range(0, seq_len, 3)
                        if all(codon_aln[g][i:i + 3] != '---' for g in ids))
    if n_informative < min_codons:
        result['status'] = f'too_short:{n_informative}_codons'
        return result

    seqfile = os.path.join(workdir, 'input.phy')
    with open(seqfile, 'w') as fh:
        fh.write(f'  {len(ids)}  {seq_len}\n')
        for gid in ids:
            fh.write(f'{species_for_id[gid]:<14}  {codon_aln[gid]}\n')

    # PAML tree files need the "<ntax> <ntree>" header line.
    treefile = os.path.join(workdir, 'tree.nwk')
    with open(treefile, 'w') as fh:
        fh.write(f'  {len(ids)}  1\n{TREE_TOPOLOGY}\n')

    outfile = os.path.join(workdir, 'codeml.out')
    ctlfile = os.path.join(workdir, 'codeml.ctl')
    with open(ctlfile, 'w') as fh:
        fh.write(CODEML_CTL_TEMPLATE.format(seqfile=seqfile, treefile=treefile,
                                            outfile=outfile))

    # codeml commonly exits 1 on a successful run, so the output file's
    # existence is the success test, not the return code.
    try:
        proc = subprocess.run(['codeml', ctlfile], capture_output=True,
                              text=True, timeout=300, cwd=workdir)
        if proc.returncode != 0 and (not os.path.exists(outfile)
                                     or os.path.getsize(outfile) == 0):
            result['status'] = f'codeml_error:{proc.stderr.strip()[:80]}'
            return result
    except FileNotFoundError:
        result['status'] = 'codeml_not_found'
        return result
    except subprocess.TimeoutExpired:
        result['status'] = 'codeml_timeout'
        return result

    branch_data = parse_branch_table(outfile)
    if len(branch_data) < 4:
        result['status'] = f'parse_failed:{len(branch_data)}_branches'
        return result

    # label_branches() refuses any table without exactly 2n - 3 = 5 branches.
    # Six is what a rooted tree produces (the root branch split in two); that
    # must never be read as a Toxocara value again, so it becomes a status.
    try:
        terminal, internal = label_branches(branch_data)
    except ValueError:
        result['status'] = f'bad_branch_table:{len(branch_data)}_branches'
        return result
    for sp, data in terminal.items():
        for m in METRICS:
            result[f'{sp}_{m}'] = data[m]

    if internal:
        b = internal[0]
        result['internal_omega'] = b['omega']
        result['internal_dN'] = b['dN']
        result['internal_dS'] = b['dS']
        result['internal_subtends'] = '|'.join(b['subtends'])

    missing = [sp for sp in SPECIES_ORDER if np.isnan(result[f'{sp}_omega'])]
    if missing:
        result['status'] = f'missing_branches:{",".join(missing)}'
        return result

    result['status'] = 'OK'
    return result


# ===========================================================================
# Worker
# ===========================================================================

def quartet_payload(quartet, prot, cds, tmpdir, min_codons):
    """The argument tuple for one `process_quartet` call.

    Carries only this quartet's eight records. ProcessPoolExecutor pickles the
    arguments of every task it hands to a worker, and the earlier version passed
    the whole protein and CDS dictionaries of all four species -- about 150 MB
    per task, pickled serially in the parent, so roughly 900 GB and over an hour
    of pickling across 6,258 quartets. That, not codeml, is most of why the v2
    run took 1.5 h. Results are unchanged: a gene missing from a FASTA is simply
    absent from its one-entry dict, and process_quartet reports it the same way.
    """
    as_g, pu_g, bsc_g, tox_g = quartet

    def one(db, gid):
        return {gid: db[gid]} if gid in db else {}

    return (as_g, pu_g, bsc_g, tox_g,
            one(prot['as'], as_g), one(prot['pu'], pu_g),
            one(prot['bsc'], bsc_g), one(prot['tox'], tox_g),
            one(cds['as'], as_g), one(cds['pu'], pu_g),
            one(cds['bsc'], bsc_g), one(cds['tox'], tox_g),
            tmpdir, min_codons)


def process_quartet(args_tuple):
    """Align, back-translate and run codeml for one quartet."""
    (as_g, pu_g, bsc_g, tox_g,
     prot_as, prot_pu, prot_bsc, prot_tox,
     cds_as, cds_pu, cds_bsc, cds_tox, tmpdir, min_codons) = args_tuple

    name = re.sub(r'[^A-Za-z0-9_.]', '_', f'{as_g}__{pu_g}__{bsc_g}__{tox_g}')

    result = {'as_gene': as_g, 'pu_gene': pu_g, 'bsc_gene': bsc_g,
              'tox_gene': tox_g, 'align_len': 0}
    result.update(_empty_result())

    for gid, db, sp in ((as_g, prot_as, 'As'), (pu_g, prot_pu, 'Pu'),
                        (bsc_g, prot_bsc, 'Bsc'), (tox_g, prot_tox, 'Tox')):
        if gid not in db:
            result['status'] = f'prot_not_found:{sp}'
            return result
    for gid, db, sp in ((as_g, cds_as, 'As'), (pu_g, cds_pu, 'Pu'),
                        (bsc_g, cds_bsc, 'Bsc'), (tox_g, cds_tox, 'Tox')):
        if gid not in db:
            result['status'] = f'cds_not_found:{sp}'
            return result

    species_for_id = {as_g: 'Ascaris', pu_g: 'Parascaris',
                      bsc_g: 'Baylisascaris', tox_g: 'Toxocara'}

    aligned = align_proteins(
        [prot_as[as_g], prot_pu[pu_g], prot_bsc[bsc_g], prot_tox[tox_g]],
        tmpdir, name)
    if aligned is None:
        result['status'] = 'alignment_failed'
        return result

    codon_aln = back_translate(aligned, {as_g: cds_as[as_g], pu_g: cds_pu[pu_g],
                                         bsc_g: cds_bsc[bsc_g],
                                         tox_g: cds_tox[tox_g]})
    if codon_aln is None:
        result['status'] = 'backtranslation_failed'
        return result

    ids = list(codon_aln.keys())
    seq_len = len(codon_aln[ids[0]])
    result['align_len'] = sum(1 for i in range(0, seq_len, 3)
                              if all(codon_aln[g][i:i + 3] != '---'
                                     for g in ids))

    result.update(run_codeml_free_ratio(codon_aln, species_for_id, name,
                                        tmpdir, min_codons))
    return result


# ===========================================================================
# Region annotation -- through region_lib, so it agrees with annotate_regions
# ===========================================================================

def annotate_ascaris_region(genes, gff_path, retained_bed, eliminated_bed):
    """Region call per Ascaris gene, vectorised, on 1-based inclusive intervals.

    `region_lib.parse_bed` does the 0-based half-open to 1-based inclusive
    conversion once, and `merge_intervals` collapses overlaps before the
    membership test. Both matter for agreement with `annotate_regions.py`, which
    is what T1 and every region-labelled figure in the paper are built on.
    """
    models = parse_gff_models(gff_path, feature='mRNA')
    primary, diag = select_primary_alignment(models, id_col='gene_id')
    if diag.get('n_changed_vs_first_in_file'):
        logging.info('  primary-alignment rule changed %d of %d models '
                     'versus first-in-file', diag['n_changed_vs_first_in_file'],
                     len(primary))

    pos = primary.set_index('gene_id')
    mid = ((pos['start'] + pos['end']) // 2).astype(float)
    chrom = pos['chrom']

    ret_merged = merge_intervals(parse_bed(retained_bed))
    elim_merged = merge_intervals(parse_bed(eliminated_bed))

    calls = pd.Series('unannotated', index=pos.index, dtype=object)
    for ch, sub in chrom.groupby(chrom):
        idx = sub.index
        p = mid.loc[idx].to_numpy()
        is_e = (in_intervals(p, *elim_merged[ch]) if ch in elim_merged
                else np.zeros(len(p), bool))
        is_r = (in_intervals(p, *ret_merged[ch]) if ch in ret_merged
                else np.zeros(len(p), bool))
        lab = np.where(is_e & is_r, 'ambiguous',
                       np.where(is_e, 'eliminated',
                                np.where(is_r, 'retained', 'unannotated')))
        calls.loc[idx] = lab

    n_amb = int((calls == 'ambiguous').sum())
    if n_amb:
        logging.warning('  %d genes fall inside BOTH a retained and an '
                        'eliminated interval; reported as ambiguous rather '
                        'than silently resolved', n_amb)
    return genes.map(calls).fillna('unannotated')


def load_tox_eliminated_ids(bed_path):
    """Eliminated Toxocara transcript IDs from column 4 of the BED."""
    ids = set()
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                ids.add(parts[3])
    return ids


def strip_td_suffix(gene_id):
    """MSTRG.524.3.p1 -> MSTRG.524.3"""
    return re.sub(r'\.p\d+$', '', gene_id)


# ===========================================================================
# Stats
# ===========================================================================

def branch_stats_table(df, region_col, omega_cap, ds_max=np.inf, seed=0):
    """Per-branch retained-vs-eliminated stats for omega, with bootstrap CIs.

    BH q-values across the four branches within this region labelling. The old
    version reported four raw p-values with no correction while
    `branch_dnds_decomposition.py` corrected the same family, so T2 and Figure
    5C could not both be right.
    """
    rows = []
    for sp in SPECIES_ORDER:
        keep = usable_mask(df, sp, 'omega', omega_cap, apply_cap=True,
                           ds_max=ds_max)
        res, ret, elim = split_and_compare(df, keep, region_col, f'{sp}_omega')
        rlo, rhi = bootstrap_median_ci(ret, seed=seed)
        elo, ehi = bootstrap_median_ci(elim, seed=seed + 1)
        res.update({'branch': sp, 'grouping': region_col, 'metric': 'omega',
                    'ret_ci_lo': rlo, 'ret_ci_hi': rhi,
                    'elim_ci_lo': elo, 'elim_ci_hi': ehi,
                    'n_over_cap': int((pd.to_numeric(df[f'{sp}_omega'],
                                                     errors='coerce')
                                       > omega_cap).sum())})
        rows.append(res)

    stats = pd.DataFrame(rows)
    stats['q_value'] = benjamini_hochberg(stats['p_value'].to_numpy(float))
    for _, r in stats.iterrows():
        logging.info('  %-14s ret %.4f (n=%s)  elim %.4f (n=%s)  '
                     'fold %.2f  %s  q=%.3g %s',
                     r['branch'], r['ret_median'], f"{int(r['ret_n']):,}",
                     r['elim_median'], f"{int(r['elim_n']):,}",
                     r['fold_change'], fmt_p(r['p_value']).replace('$', ''),
                     r['q_value'], stars(r['q_value']))
    return stats


# ===========================================================================
# Figures
# ===========================================================================

def plot_branch_omega(df, region_col, outdir, stem, omega_cap, rng,
                      stats, yscale='log', width='double', height_mm=58.0,
                      panel_letters=False):
    """Figure 5A -- branch dN/dS by region, one panel per terminal branch.

    Box statistics come from the full masked set, so the median line is the
    median in the stats table. Whiskers are the 5th and 95th percentiles: on a
    log axis, 1.5 IQR whiskers reach the minimum of any right-skewed sample.

    dN/dS at codeml's floor of 0.0001 means no nonsynonymous change was
    detected on the branch. On the log axis those genes are drawn in a band at
    0 below the scale break, with their count beside each box -- they used to
    sit on the bottom spine, uncounted, because the old test for them was
    `<= 0`. (This is not the same count as dN = 0 in Figure 5C: on a long
    branch codeml prints such a gene's dN as 0.0001. See branch_lib.)

    The bottom margin is measured after layout (fit_bottom_margin), so the
    labels end 1 mm above the figure edge whatever their height.

    Returns {(branch, region): genes with dN/dS at the floor}.
    """
    apply_style()
    fig, axes = plt.subplots(1, 4, figsize=(COL[width], mm(height_mm)),
                             sharey=True, layout='none')

    series, zeros = [], {}
    for sp in SPECIES_ORDER:
        keep = usable_mask(df, sp, 'omega', omega_cap, apply_cap=True)
        _, ret, elim = split_and_compare(df, keep, region_col, f'{sp}_omega')
        series.append((ret, elim))
        zeros[(sp, 'retained')] = int(at_zero(ret, 'omega').sum())
        zeros[(sp, 'eliminated')] = int(at_zero(elim, 'omega').sum())

    log = yscale == 'log'
    if log:
        # hi_cap: every drawn value is <= omega_cap, so nothing is clipped.
        lo, hi, _ = zero_band_limits([s for pair in series for s in pair],
                                     'omega', hi_cap=omega_cap)

    for ax, sp, (ret, elim) in zip(axes, SPECIES_ORDER, series):
        if len(ret) < 5 or len(elim) < 5:
            ax.text(0.5, 0.5, 'insufficient data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=FONT['annotation'])
            ax.set_xlabel(SPECIES_LABEL[SPECIES_KEY[sp]])
            continue

        data = [as_plotted(ret, 'omega'), as_plotted(elim, 'omega')]
        bp = ax.boxplot(data, positions=[0, 1], widths=0.55,
                        patch_artist=True, showfliers=False, whis=(5, 95),
                        zorder=2)
        for patch, key in zip(bp['boxes'], ('retained', 'eliminated')):
            # Alpha in the FACE colour. set_alpha() would fade the black edge
            # to mid-grey and the box would lose its boundary at print size.
            patch.set_facecolor(to_rgba(REGION[key], 0.40))
            patch.set_edgecolor('black')
            patch.set_alpha(None)
        for med in bp['medians']:
            med.set_color('black')
            med.set_zorder(5)
        for art in bp['whiskers'] + bp['caps']:
            art.set_zorder(4)

        for i, (v, key) in enumerate(zip(data, ('retained', 'eliminated'))):
            if len(v) > 1500:
                v = v[rng.choice(len(v), 1500, replace=False)]
            ax.scatter(np.full(len(v), i) + rng.normal(0, 0.07, len(v)), v,
                       s=1.6, alpha=0.20, color=REGION[key], linewidths=0,
                       zorder=3, rasterized=True)

        row = stats[stats['branch'] == sp]
        q = float(row['q_value'].iloc[0]) if len(row) else np.nan
        fold = float(row['fold_change'].iloc[0]) if len(row) else np.nan

        if log:
            mark_axis_break(ax, zero_band_axis(ax, lo, hi))
            for i, key in enumerate(('retained', 'eliminated')):
                if zeros[(sp, key)]:
                    # House annotation size (7 pt); it was 1 pt smaller.
                    ax.text(i + 0.30, 0.0, f'{zeros[(sp, key)]:,}',
                            ha='left', va='center', zorder=6,
                            fontsize=FONT['annotation'],
                            color=NEUTRAL['ink_light'])
        else:
            ax.set_ylim(0, omega_cap)
        ax.axhline(1.0, color=NEUTRAL['rule'], ls='--', lw=0.9, zorder=1)
        ax.set_xlim(-0.6, 1.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'Retained\n(n = {len(ret):,})',
                            f'Eliminated\n(n = {len(elim):,})'])
        ax.set_xlabel(SPECIES_LABEL[SPECIES_KEY[sp]], labelpad=2)
        # Above the panel, not inside it: inside, the jitter of the retained
        # box reached the text on the Ascaris and Baylisascaris branches.
        ax.text(0.5, 1.03, f'{stars(q)}  {fold:.2f}×',
                transform=ax.transAxes, ha='center', va='bottom',
                fontsize=FONT['annotation'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].set_ylabel(METRIC_LABEL['omega'])
    # Direct-label the rule instead of a legend key. A one-entry legend inside
    # the axes has to sit somewhere, and every corner of these panels has
    # jitter in it -- audit() flagged the overlap. The label rides the line it
    # names, which is also what the house rule prefers for <= 4 series.
    axes[-1].text(1.58, 1.0, 'dN/dS = 1', ha='right', va='bottom',
                  fontsize=FONT['annotation'], color=NEUTRAL['rule'])
    if panel_letters:
        from paper_colors import panel_letter
        panel_letter(axes[0], 'A')

    # bottom is a starting value only; fit_bottom_margin() sets it from the
    # drawn labels. The fixed 0.30 left 5-6 mm of empty figure under them.
    fig.subplots_adjust(left=0.075, right=0.995, top=0.91, bottom=0.30,
                        wspace=0.10)
    logging.debug('  %s: bottom margin moved %.1f mm', stem,
                  fit_bottom_margin(fig))
    audit(fig)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return zeros


def plot_branch_medians(stats, outdir, stem, width='single', height_mm=52.0,
                        panel_letters=False):
    """Figure 5B -- median branch dN/dS per lineage, with bootstrap CIs.

    No neutrality line. The tallest bar is around 0.4, so drawing dN/dS = 1
    would spend most of the axis proving a point the data never approaches; the
    summary states the maximum branch median instead. The CIs are percentile
    bootstraps on the median -- the previous version drew `Series.sem()`, the
    standard error of the mean, against a median.
    """
    apply_style()
    fig, ax = plt.subplots(figsize=(COL[width], mm(height_mm)),
                           layout='none')

    x = np.arange(len(SPECIES_ORDER))
    w = 0.36
    for off, key, col in ((-w / 2, 'retained', 'ret'),
                          (+w / 2, 'eliminated', 'elim')):
        med, lo, hi = [], [], []
        for sp in SPECIES_ORDER:
            r = stats[stats['branch'] == sp]
            m = float(r[f'{col}_median'].iloc[0]) if len(r) else np.nan
            med.append(m)
            lo.append(m - float(r[f'{col}_ci_lo'].iloc[0]) if len(r) else 0)
            hi.append(float(r[f'{col}_ci_hi'].iloc[0]) - m if len(r) else 0)
        ax.bar(x + off, med, w, yerr=[lo, hi], capsize=2,
               facecolor=to_rgba(REGION[key], 0.55), edgecolor='black',
               linewidth=0.6, label=key.capitalize(),
               error_kw={'lw': 0.8, 'ecolor': 'black'})

    top = 0.0
    for sp in SPECIES_ORDER:
        r = stats[stats['branch'] == sp]
        if len(r):
            top = max(top, float(r['elim_ci_hi'].iloc[0]),
                      float(r['ret_ci_hi'].iloc[0]))
    ax.set_ylim(0, top * 1.28)

    for xi, sp in enumerate(SPECIES_ORDER):
        r = stats[stats['branch'] == sp]
        if not len(r):
            continue
        q = float(r['q_value'].iloc[0])
        y = max(float(r['elim_ci_hi'].iloc[0]), float(r['ret_ci_hi'].iloc[0]))
        ax.text(xi, y + top * 0.05, stars(q), ha='center', va='bottom',
                fontsize=FONT['annotation'])

    ax.set_xticks(x)
    ax.set_xticklabels([SPECIES_LABEL[SPECIES_KEY[sp]] for sp in SPECIES_ORDER])
    ax.set_ylabel('Median branch dN/dS')
    # Above the axes, not inside it. In the axes the key landed on the Ascaris
    # significance marker -- which audit() does not catch, because it tests
    # lines, collections and patches but not Text.
    ax.legend(frameon=False, fontsize=FONT['legend'], ncol=2, handlelength=1.2,
              loc='lower left', bbox_to_anchor=(0.0, 1.005), columnspacing=1.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if panel_letters:
        from paper_colors import panel_letter
        panel_letter(ax, 'B')

    # bottom: starting value only, as in plot_branch_omega().
    fig.subplots_adjust(left=0.19, right=0.99, top=0.89, bottom=0.20)
    logging.debug('  %s: bottom margin moved %.1f mm', stem,
                  fit_bottom_margin(fig))
    audit(fig)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))


# ===========================================================================
# Main
# ===========================================================================

#: The inputs a full run needs and a --replot does not take. argparse cannot
#: express "required unless --replot", so parse_args() checks them itself.
INPUT_ARGS = ('rbh_as_pu', 'rbh_as_bsc', 'rbh_pu_bsc', 'rbh_tox_as',
              'rbh_tox_pu', 'rbh_tox_bsc', 'cds_as', 'prot_as', 'cds_pu',
              'prot_pu', 'cds_bsc', 'prot_bsc', 'cds_tox', 'prot_tox',
              'retained1', 'eliminated1', 'gff1')

#: Defaults of the settings a replot must share with the run it redraws. One
#: definition, used by argparse and by source_run_settings(), which assumes
#: them for a source run that did not pass the flag. They have not changed
#: since the 2026-08-18 conversion.
DEFAULT_OMEGA_CAP = 5.0
DEFAULT_SEED = 0
DEFAULT_YSCALE = 'log'

def parse_args():
    p = argparse.ArgumentParser(
        description='Four-taxon branch-specific dN/dS via codeml free-ratio.',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required for a full run, refused with --replot (checked below).
    p.add_argument('--rbh-as-pu')
    p.add_argument('--rbh-as-bsc')
    p.add_argument('--rbh-pu-bsc')
    p.add_argument('--rbh-tox-as', help='Toxocara->Ascaris RBH')
    p.add_argument('--rbh-tox-pu', help='Toxocara->Parascaris RBH')
    p.add_argument('--rbh-tox-bsc', help='Toxocara->Baylisascaris RBH')

    for sp in ('as', 'pu', 'bsc', 'tox'):
        p.add_argument(f'--cds-{sp}')
        p.add_argument(f'--prot-{sp}')

    p.add_argument('--retained1', help='Ascaris retained BED')
    p.add_argument('--eliminated1', help='Ascaris eliminated BED')
    p.add_argument('--gff1', help='Ascaris miniprot GFF')
    p.add_argument('--tox-eliminated', default=None,
                   help='Toxocara eliminated transcripts BED (enables S14)')

    p.add_argument('--replot', metavar='ANNOTATED_TSV', default=None,
                   help='Redraw Figure 5A/5B (and S14A/B), Table T2 and the '
                        'summary from an existing branch_dnds_four_annotated.tsv. '
                        'No ortholog search and no codeml, so the numbers stay '
                        'those of the run that wrote the table; a full rerun '
                        'refits every quartet from random starts. The input '
                        'files above are not taken with --replot.')

    p.add_argument('--omega-cap', type=float, default=DEFAULT_OMEGA_CAP,
                   help='Mask a branch whose own omega exceeds this. Applied '
                        'per branch, so a gene is never dropped whole '
                        '(default: 5)')
    p.add_argument('--min-codons', type=int, default=50,
                   help='Minimum gap-free codons in the quartet alignment '
                        '(default: 50)')
    p.add_argument('--yscale', choices=['log', 'linear'], default=DEFAULT_YSCALE,
                   help='y-axis of Figure 5A. Linear to the omega cap squashes '
                        'every box into the bottom tenth of the panel '
                        '(default: log)')
    p.add_argument('--panel-letters', action='store_true',
                   help='Draw A/B on the panels. Off by default -- panel '
                        'letters are added at assembly.')
    p.add_argument('--seed', type=int, default=DEFAULT_SEED,
                   help='RNG seed for the bootstrap CIs and the 5A jitter '
                        f'(default: {DEFAULT_SEED})')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--keep-tmp', action='store_true')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()

    flag = lambda a: '--' + a.replace('_', '-')                    # noqa: E731
    if args.replot is None:
        missing = [flag(a) for a in INPUT_ARGS if getattr(args, a) is None]
        if missing:
            p.error('a full run needs ' + ', '.join(missing)
                    + ' (or redraw an existing run with --replot)')
    else:
        # Refused rather than ignored: a new --tox-eliminated passed here would
        # look as if it relabelled S14, and it would not.
        given = [flag(a) for a in INPUT_ARGS + ('tox_eliminated',)
                 if getattr(args, a) is not None]
        if given:
            p.error('--replot draws from the annotated table only; drop '
                    + ', '.join(given))
        if not os.path.isfile(args.replot):
            p.error(f'--replot: no such file: {args.replot}')

    if args.outdir is None:
        args.outdir = make_outdir('branch_dnds_four_replot' if args.replot
                                  else 'branch_dnds_four', base='results')
    return args


def source_run_settings(annotated_path):
    """Settings of the run that wrote `annotated_path`.

    Reads the provenance.json in the same directory. Returns (settings,
    provenance): settings maps 'omega_cap', 'seed' and 'yscale' to the value
    that run used -- from its command line, or the default when the flag was
    not given; provenance is the whole record. Returns ({}, None) when there
    is no readable provenance.json, or when it is not from this script.
    """
    prov_path = os.path.join(os.path.dirname(os.path.abspath(annotated_path)),
                             'provenance.json')
    if not os.path.isfile(prov_path):
        return {}, None
    try:
        with open(prov_path) as fh:
            prov = json.load(fh)
    except (OSError, ValueError) as exc:
        logging.warning('  could not read %s: %s', prov_path, exc)
        return {}, None
    if prov.get('script') != os.path.basename(__file__):
        logging.warning('  %s was written by %r, not %s; settings not checked',
                        prov_path, prov.get('script'),
                        os.path.basename(__file__))
        return {}, None

    argv = [str(a) for a in prov.get('argv', [])]
    settings = {'omega_cap': DEFAULT_OMEGA_CAP, 'seed': DEFAULT_SEED,
                'yscale': DEFAULT_YSCALE}
    for opt, key, cast in (('--omega-cap', 'omega_cap', float),
                           ('--seed', 'seed', int), ('--yscale', 'yscale', str)):
        for i, a in enumerate(argv):
            try:
                if a == opt and i + 1 < len(argv):
                    settings[key] = cast(argv[i + 1])
                elif a.startswith(opt + '='):
                    settings[key] = cast(a.split('=', 1)[1])
            except ValueError:
                logging.warning('  unreadable %s in %s: %r', opt, prov_path, a)
                settings.pop(key, None)
    return settings, prov


def run_codeml_stage(args):
    """Full run up to the annotated table: four-way orthologs, codeml on every
    quartet, region labels. Writes branch_dnds_four_raw.tsv and
    branch_dnds_four_annotated.tsv.

    Returns (df_ok, groupings, n_quartets, prov_extra, tmpdir).
    """
    # Before any codeml job: a rooted string would be fitted as rooted, with
    # two root branches of which only the sum is estimable (2026-10-10).
    try:
        check_unrooted(TREE_TOPOLOGY)
    except ValueError as exc:
        sys.exit(f'ERROR: {exc}')

    if shutil.which('codeml') is None:
        sys.exit('ERROR: codeml not found in PATH. Is the kaks env active?')
    if shutil.which('mafft') is None:
        sys.exit('ERROR: mafft not found in PATH.')

    quartets, ortho_diag = find_four_way_orthologs(
        args.rbh_as_pu, args.rbh_as_bsc, args.rbh_pu_bsc,
        args.rbh_tox_as, args.rbh_tox_pu, args.rbh_tox_bsc)
    if not quartets:
        sys.exit('ERROR: no four-way orthologs found. Check RBH orientations.')

    logging.info('Loading sequences...')
    prot = {k: load_fasta(getattr(args, f'prot_{k}'))
            for k in ('as', 'pu', 'bsc', 'tox')}
    cds = {k: load_fasta(getattr(args, f'cds_{k}'))
           for k in ('as', 'pu', 'bsc', 'tox')}

    tmpdir = tempfile.mkdtemp(prefix='branch4_')
    logging.info('Processing %s quartets on %d thread(s); tmp %s',
                 f'{len(quartets):,}', args.threads, tmpdir)

    results = []
    if args.threads == 1:
        for i, q in enumerate(quartets):
            if (i + 1) % 200 == 0:
                logging.info('  %d/%d', i + 1, len(quartets))
            results.append(process_quartet(
                quartet_payload(q, prot, cds, tmpdir, args.min_codons)))
    else:
        tasks = []
        for i, q in enumerate(quartets):
            wtmp = os.path.join(tmpdir, f'w_{i}')
            os.makedirs(wtmp, exist_ok=True)
            tasks.append(quartet_payload(q, prot, cds, wtmp, args.min_codons))
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            futures = {ex.submit(process_quartet, t): i
                       for i, t in enumerate(tasks)}
            done = 0
            for fut in as_completed(futures):
                done += 1
                if done % 200 == 0:
                    logging.info('  %d/%d', done, len(quartets))
                try:
                    results.append(fut.result())
                except Exception as exc:                      # noqa: BLE001
                    q = quartets[futures[fut]]
                    logging.error('  exception for %s: %s', q, exc)
                    r = {'as_gene': q[0], 'pu_gene': q[1], 'bsc_gene': q[2],
                         'tox_gene': q[3], 'align_len': 0}
                    r.update(_empty_result())
                    r['status'] = f'exception:{str(exc)[:60]}'
                    results.append(r)

    df = pd.DataFrame(results)
    df.to_csv(os.path.join(args.outdir, 'branch_dnds_four_raw.tsv'),
              sep='\t', index=False, float_format='%.10g')

    logging.info('Status:')
    for s, n in df['status'].value_counts().items():
        logging.info('  %-32s %s', s, f'{n:,}')

    df_ok = df[df['status'] == 'OK'].copy()
    logging.info('Usable quartets: %s of %s', f'{len(df_ok):,}', f'{len(df):,}')
    if not len(df_ok):
        sys.exit('ERROR: no quartet produced a complete branch table.')

    logging.info('Annotating by Ascaris region...')
    df_ok['as_region'] = annotate_ascaris_region(
        df_ok['as_gene'], args.gff1, args.retained1, args.eliminated1).values
    logging.info('  %s', df_ok['as_region'].value_counts().to_dict())

    groupings = ['as_region']
    if args.tox_eliminated:
        logging.info('Annotating by Toxocara elimination status...')
        tox_ids = load_tox_eliminated_ids(args.tox_eliminated)
        df_ok['tox_region'] = df_ok['tox_gene'].apply(
            lambda g: 'eliminated' if strip_td_suffix(g) in tox_ids
            else 'retained')
        logging.info('  %s', df_ok['tox_region'].value_counts().to_dict())
        groupings.append('tox_region')

    annotated = os.path.join(args.outdir, 'branch_dnds_four_annotated.tsv')
    df_ok.to_csv(annotated, sep='\t', index=False, float_format='%.10g')
    logging.info('Annotated table: %s', annotated)

    prov_extra = {
        'mode': 'full',
        'ortholog_diagnostics': ortho_diag,
        'inputs': {k: fingerprint(getattr(args, k))
                   for k in ('rbh_as_pu', 'rbh_as_bsc', 'rbh_pu_bsc',
                             'rbh_tox_as', 'rbh_tox_pu', 'rbh_tox_bsc', 'gff1',
                             'retained1', 'eliminated1', 'tox_eliminated')
                   if getattr(args, k, None)}}
    return df_ok, groupings, len(quartets), prov_extra, tmpdir


def load_for_replot(args):
    """The annotated table of an earlier full run, validated, plus what that
    run's provenance.json records.

    Returns (df_ok, groupings, n_quartets, prov_extra). n_quartets is None
    when the source run left no provenance.
    """
    settings, src = source_run_settings(args.replot)
    if src is None:
        logging.warning('  no readable provenance.json beside the table: '
                        'cannot check --omega-cap and --seed against the run '
                        'that wrote it')
    for key in ('omega_cap', 'seed'):
        if key in settings and settings[key] != getattr(args, key):
            logging.warning('  the source run used --%s %s, this replot %s: '
                            'Table T2 and the figures will NOT match that run',
                            key.replace('_', '-'), settings[key],
                            getattr(args, key))
    if 'yscale' in settings and settings['yscale'] != args.yscale:
        logging.info('  the source run drew 5A on a %s axis, this replot %s',
                     settings['yscale'], args.yscale)

    # load_annotated() refuses a table whose internal branch subtends anything
    # but Ascaris|Parascaris (made with the rooted tree). A table with no
    # internal_subtends column at all predates the 2026-10-09 conversion, and
    # every run before 2026-10-10 used the rooted tree, so refuse that too.
    df, groupings = load_annotated(args.replot)
    if 'internal_subtends' not in df.columns:
        sys.exit(f'ERROR: {args.replot} has no internal_subtends column, so '
                 'it was written before 2026-10-09, with the rooted tree. '
                 'Replot only a table from a 2026-10-10 or later run.')
    if 'status' in df.columns:
        bad = df['status'].astype(str) != 'OK'
        if bad.any():
            logging.warning('  %s rows with a status other than OK dropped',
                            f'{int(bad.sum()):,}')
            df = df[~bad].copy()
    if not len(df):
        sys.exit('ERROR: no usable rows in the annotated table.')

    keep = ('script', 'run_utc', 'git_commit', 'host', 'argv', 'n_quartets',
            'n_ok')
    prov_extra = {
        'mode': 'replot',
        'replot_of': fingerprint(args.replot),
        'source_run': ({k: src.get(k) for k in keep} if src else None),
        'ortholog_diagnostics': (src.get('ortholog_diagnostics')
                                 if src else None)}
    return df, groupings, (src or {}).get('n_quartets'), prov_extra


def report(df_ok, groupings, args, rng, n_quartets, prov_extra):
    """Everything from the annotated table onward: the internal-split check,
    masking counts, Table T2, Figure 5A/5B (and S14A/B), the summary and
    provenance.json.

    A full run and --replot both end here, so a replot draws exactly what a
    full run would have drawn from the same table with the same cap and seed.
    `rng` is used only for the 5A jitter, in the same order either way.
    """
    subtends = []
    if 'internal_subtends' in df_ok.columns:
        # dropna/astype: a table read back from disk has NaN where the full
        # run had '', and NaN is truthy.
        subtends = [s for s in
                    df_ok['internal_subtends'].dropna().astype(str).unique()
                    if s]
    logging.info('Internal branch: %s',
                 internal_branch_note([{'branch': '-',
                                        'subtends': tuple(subtends[0].split('|'))}]
                                      if subtends else []))
    if len(subtends) > 1:
        logging.warning('  %d different internal splits across genes: %s -- '
                        'the topology is not being recovered consistently',
                        len(subtends), subtends)
    elif subtends != ['Ascaris|Parascaris']:
        logging.warning('  internal branch subtends %s, expected '
                        'Ascaris|Parascaris for %s', subtends, TREE_TOPOLOGY)

    logging.info('Per-branch masking at dN/dS > %g, plus dS = 0 with dN/dS '
                 'under the cap (dN = dS = 0, dN/dS undefined); applied to that '
                 'branch only, no gene is dropped whole:', args.omega_cap)
    for sp in SPECIES_ORDER:
        om = pd.to_numeric(df_ok[f'{sp}_omega'], errors='coerce')
        ds = pd.to_numeric(df_ok[f'{sp}_dS'], errors='coerce')
        n_over = int((om > args.omega_cap).sum())
        n_00 = int(((om <= args.omega_cap) & (ds <= 0)).sum())
        logging.info('    %-14s %s over the cap, %s at 0/0', sp, f'{n_over:,}',
                     f'{n_00:,}')

    all_stats, n_zero_total = [], {}
    for g in groupings:
        logging.info('Branch dN/dS by %s:', g)
        st = branch_stats_table(df_ok, g, args.omega_cap, seed=args.seed)
        all_stats.append(st)

        sub = df_ok[df_ok[g].isin(['retained', 'eliminated'])]
        tag = 'fig5' if g == 'as_region' else 'figS14'
        n_zero_total[g] = plot_branch_omega(
            sub, g, args.outdir, f'{tag}A_branch_dnds_{g}', args.omega_cap,
            rng, st, yscale=args.yscale, panel_letters=args.panel_letters)
        plot_branch_medians(st, args.outdir, f'{tag}B_branch_medians_{g}',
                            panel_letters=args.panel_letters)

    stats_df = pd.concat(all_stats, ignore_index=True)
    stats_path = os.path.join(args.outdir, 'branch_stats.tsv')
    stats_df.to_csv(stats_path, sep='\t', index=False, float_format='%.6g')
    logging.info('Table T2: %s', stats_path)

    # --- summary ---------------------------------------------------------
    n_q = f'{n_quartets:,}' if n_quartets is not None else 'not recorded'
    L = ['Four-taxon branch-specific dN/dS', '']
    if args.replot:
        L += [f'REPLOT of         {os.path.abspath(args.replot)}',
              '                  (figures, T2 and summary only; codeml was '
              'not rerun)']
    L += [f'Tree              {TREE_TOPOLOGY}  (unrooted; runs before '
          '2026-10-10 used a rooted string, which codeml fits as rooted)',
          f'Four-way orthologs {n_q}',
          f'Usable quartets    {len(df_ok):,}',
          f'Per-branch cap     dN/dS <= {args.omega_cap:g}',
          f'Figure 5A y-axis   {args.yscale}', '',
          'Per-branch retained vs eliminated', '',
          stats_df.to_string(index=False), '', 'Notes', '']

    max_med = np.nanmax(stats_df[['ret_median', 'elim_median']].to_numpy(float))
    L += [f'  The largest branch median is {max_med:.3f}. No branch approaches',
          '  dN/dS = 1 in either compartment, which is why Figure 5B carries no',
          '  neutrality line -- drawing it would spend most of the axis on empty',
          '  space. The line stays on 5A, where the log scale makes it free.', '']

    if subtends:
        L += ['  ' + internal_branch_note(
            [{'branch': '-', 'subtends': tuple(subtends[0].split('|'))}]), '',
            '  Do NOT describe this as an ancestral or ascarid-stem branch. It',
            '  is the single internal edge of an UNROOTED four-taxon tree and',
            '  belongs to neither side of the split.', '']

    L += ['  Masking: a branch is dropped for a gene when its dN/dS exceeds the',
          f'  cap ({args.omega_cap:g}) or when dS = 0 under the cap (dN = dS = 0: no',
          '  substitution, dN/dS undefined, and codeml\'s value depends on its',
          '  random start). Counts per branch are in the log.', '']
    L += ['  dN/dS at codeml\'s lower bound, 0.0001: no nonsynonymous change was',
          '  detected on that branch. Figure 5A draws those genes in the band at',
          '  0 below the scale break; the grey number beside each box is this',
          '  count. Box statistics include them, so the drawn median matches the',
          '  table. This is not the dN = 0 count of Figure 5C: codeml prints dN',
          '  as (dN/dS) x dS to four decimals, so on a long branch (Toxocara) a',
          '  gene at the floor shows dN = 0.0001.', '']
    for g in groupings:
        zc = n_zero_total.get(g, {})
        L.append(f'  {g:<11} ' + '   '.join(
            f"{sp} {zc.get((sp, 'retained'), 0):,} / "
            f"{zc.get((sp, 'eliminated'), 0):,}" for sp in SPECIES_ORDER)
            + '   (retained / eliminated)')
    L.append('')

    L += ['  q-values are BH across the four branches within one region',
          '  labelling, matching branch_dnds_decomposition.py. Quoting raw P',
          '  from this table would disagree with Figure 5C.', '']

    with open(os.path.join(args.outdir, 'branch_dnds_four_summary.txt'),
              'w') as fh:
        fh.write('\n'.join(L) + '\n')

    prov = {'script': os.path.basename(__file__),
            'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(
                microsecond=0).isoformat(),
            'host': platform.node(),
            'python': platform.python_version(),
            'argv': sys.argv,
            'n_quartets': n_quartets, 'n_ok': int(len(df_ok))}
    prov.update(prov_extra)
    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump(prov, fh, indent=2)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'branch_dnds_four')
    apply_style()
    rng = np.random.default_rng(args.seed)

    # The PDF writer's font subsetting logs hundreds of INFO lines per figure
    # through fontTools, which buried the results in the 2026-10-09 job output.
    # Same silencing as branch_dnds_decomposition.py.
    for noisy in ('matplotlib.font_manager', 'matplotlib.texmanager',
                  'fontTools', 'fontTools.subset', 'fontTools.ttLib', 'PIL'):
        logging.getLogger(noisy).setLevel(logging.ERROR)

    tmpdir = None
    if args.replot:
        logging.info('Four-taxon branch-specific dN/dS -- REPLOT, no codeml')
        logging.info('  table:  %s', os.path.abspath(args.replot))
        logging.info('  outdir: %s', os.path.abspath(args.outdir))
        df_ok, groupings, n_quartets, prov_extra = load_for_replot(args)
    else:
        logging.info('Four-taxon branch-specific dN/dS')
        logging.info('  tree:   %s', TREE_TOPOLOGY)
        logging.info('  outdir: %s', os.path.abspath(args.outdir))
        df_ok, groupings, n_quartets, prov_extra, tmpdir = run_codeml_stage(args)

    report(df_ok, groupings, args, rng, n_quartets, prov_extra)

    if tmpdir and not args.keep_tmp:
        shutil.rmtree(tmpdir, ignore_errors=True)
    logging.info('Done.')


if __name__ == '__main__':
    main()
