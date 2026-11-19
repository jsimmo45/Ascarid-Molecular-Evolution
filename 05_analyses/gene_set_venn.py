#!/usr/bin/env python3
"""
gene_set_venn.py - which Ascaris genes enter which comparison (Figure 1).

Draws a three-set Venn diagram of the Ascaris genes in the three pairwise rate
tables that include Ascaris -- with P. univalens, B. schroederi and T. canis --
counting only the pairs that passed the rate filters (the kaks_filtered.tsv
tables, which every pairwise figure uses). Each region of the diagram gives its
retained genes (blue) above its eliminated genes (red), labelled by the
Ascaris annotation exactly as annotate_regions.py labels them: the same
region_lib calls on the same GFF and BEDs (primary alignment of each model,
region of its midpoint). The four-species quartets of the branch models are
counted against the centre and noted under the diagram.

Written in answer to Jianbin Wang's review of the figure legends (2026-10-26):
the gene counts change from figure to figure, and the diagram shows how the
pairwise sets relate before each figure states which subset it uses.
The note under the diagram uses American spelling ("center"), as the paper does.

Writes, in --outdir:
    gene_set_venn.{svg,pdf,png}     the panel, at a journal column width
    gene_set_venn_counts.tsv        every region: genes, retained, eliminated, other
    gene_set_venn_summary.txt       the counts and the checks below, for reading
    gene_set_venn.log               the run log
    provenance.json                 inputs (with checksums), arguments, versions

Checks (each stops the script with a message rather than drawing a wrong panel):
* Every Ascaris identifier in a rate table must be a model in the GFF.
* An identifier may appear only once per table (the pairs are 1:1 orthologs).
* With --annotated, the region call of every gene shared with that table must
  match its region column (for the paper: 9,085 retained and 322 eliminated
  Ascaris x P. univalens genes in results/annotate_regions_20260817).
Quartets outside the centre are not an error (a quartet comes from the ortholog
tables, before the pairwise rate filters); they are counted and reported.

Example
-------
    python gene_set_venn.py \\
      --gff As_miniprot.gff \\
      --retained-bed AG.v50.retained_loose.bed \\
      --eliminated-bed AG.v50.eliminated_strict.bed \\
      --chrom-sizes AG.v50.fa.fai \\
      --set parascaris=kaks_results/kaks_filtered.tsv:gene1 \\
      --set baylisascaris=kaks_As_Bsc_v2/kaks_filtered.tsv:gene1 \\
      --set toxocara=kaks_Tox_As/kaks_filtered.tsv:gene2 \\
      --quartets branch_dnds_four_results_v2/branch_dnds_four_annotated.tsv:as_gene \\
      --annotated results/annotate_regions_20260817/kaks_annotated.tsv \\
      --outdir results/gene_set_venn_20261026

Requires numpy, pandas and matplotlib, plus paper_colors.py and region_lib.py
from lib/, which are located automatically.
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import os
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Circle

# paper_colors.py and region_lib.py live in lib/ in the published layout and
# beside the scripts in a flat working tree; look in both. A missing module must
# stop the script: a fallback palette would quietly break the house colours.
_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, 'lib'), os.path.join(os.path.dirname(_HERE), 'lib')):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, _cand)
try:
    from paper_colors import (REGION, NEUTRAL, SPECIES_LABEL, FONT, LINEWIDTH,
                              apply_style, audit, figure, save)
    from region_lib import (parse_gff_models, select_primary_alignment, parse_bed,
                            merge_intervals, label_arrays, load_chrom_sizes,
                            derive_cbrs, annotate_positions, fingerprint,
                            git_commit, setup_logging)
except ImportError as exc:  # pragma: no cover - environment problem
    sys.exit(f'ERROR: cannot import the shared modules ({exc}). Run from the '
             'repository, or put its lib/ folder on PYTHONPATH.')

#: The three partners of Ascaris, in drawing order: top, lower left, lower right.
PARTNERS = ('parascaris', 'baylisascaris', 'toxocara')


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def parse_table_spec(spec, what):
    """'path:column' -> (path, column). The column is required: the Ascaris
    gene is gene1 in two tables and gene2 in the T. canis one."""
    if ':' not in spec:
        sys.exit(f'ERROR: {what} {spec!r} must be PATH:COLUMN, e.g. '
                 'kaks_results/kaks_filtered.tsv:gene1')
    path, col = spec.rsplit(':', 1)
    return path, col


def parse_set_spec(spec):
    """'key=path:column' -> (key, path, column), key one of PARTNERS."""
    if '=' not in spec:
        sys.exit(f'ERROR: --set {spec!r} must be KEY=PATH:COLUMN, with KEY one of '
                 f'{", ".join(PARTNERS)}')
    key, rest = spec.split('=', 1)
    if key not in PARTNERS:
        sys.exit(f'ERROR: --set key {key!r} is not one of {", ".join(PARTNERS)}')
    path, col = parse_table_spec(rest, '--set')
    return key, path, col


def read_ids(path, col, what):
    """The unique identifiers in one column of a tab-separated table.

    Read as strings: identifiers that look numeric must not lose zeros. A
    repeated identifier is an error, since every table here lists 1:1 pairs."""
    if not os.path.exists(path):
        sys.exit(f'ERROR: {what}: no such file {path}')
    df = pd.read_csv(path, sep='\t', dtype=str)
    if col not in df.columns:
        sys.exit(f'ERROR: {what}: column {col!r} not in {path} '
                 f'(columns: {", ".join(df.columns[:12])})')
    ids = df[col].dropna().str.strip()
    ids = ids[ids != '']
    dup = ids[ids.duplicated()]
    if len(dup):
        sys.exit(f'ERROR: {what}: {dup.nunique()} identifier(s) occur more than '
                 f'once in {col} of {path} (e.g. {dup.iloc[0]}); expected 1:1 pairs')
    logging.info('%s: %s rows, %s identifiers in %s (%s)', what, f'{len(df):,}',
                 f'{len(ids):,}', col, os.path.basename(path))
    return set(ids)


def region_calls(gff, retained_bed, eliminated_bed, chrom_sizes, terminal_tol,
                 cbr_window):
    """Region of every Ascaris model, as annotate_regions.annotate_species makes it.

    Same calls in the same order: primary alignment per model, midpoint, merged
    BED intervals, CBRs from the merged eliminated tracts, annotate_positions.
    Only the region column is used here. Returns (Series gene_id -> region, diag).
    """
    models = parse_gff_models(gff)
    genes, align_diag = select_primary_alignment(models)
    genes = genes.assign(midpoint=(genes['start'] + genes['end']) // 2)
    retained = parse_bed(retained_bed)
    eliminated = parse_bed(eliminated_bed, keep_label=True)
    ret_merged = merge_intervals(retained)
    elim_merged = merge_intervals(eliminated)
    elim_labels = label_arrays(eliminated, elim_merged)
    sizes, inferred = load_chrom_sizes(chrom_sizes, genes, eliminated)
    cbrs, _ = derive_cbrs(elim_merged, sizes, terminal_tol)
    annotated, counts = annotate_positions(genes, ret_merged, elim_merged,
                                           elim_labels, cbrs, cbr_window)
    diag = {'alignment': align_diag, 'region_counts_all_models': counts,
            'chrom_sizes_inferred': inferred}
    return annotated.set_index('gene_id')['region'], diag


# ---------------------------------------------------------------------------
# Counting
# ---------------------------------------------------------------------------

def collapse(region):
    """retained | eliminated | other (ambiguous or unannotated)."""
    return region if region in ('retained', 'eliminated') else 'other'


def venn_counts(sets, region):
    """One row per non-empty membership pattern of the three sets."""
    rows = []
    union = set().union(*sets.values())
    for pattern in itertools.product((True, False), repeat=3):
        if not any(pattern):
            continue
        members = [k for k, inside in zip(PARTNERS, pattern) if inside]
        ids = {g for g in union
               if all((g in sets[k]) == inside for k, inside in zip(PARTNERS, pattern))}
        calls = pd.Series([collapse(region.get(g, 'other')) for g in ids],
                          dtype=object)
        rows.append({'region': ' & '.join(members) + (' only' if len(members) == 1 else ''),
                     **{f'in_{k}': inside for k, inside in zip(PARTNERS, pattern)},
                     'genes': len(ids),
                     'retained': int((calls == 'retained').sum()),
                     'eliminated': int((calls == 'eliminated').sum()),
                     'other': int((calls == 'other').sum())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

#: Circle centres for radius 1, drawing order PARTNERS: top, lower left, lower
#: right. Centres 1.15 radii apart leave every region large enough for two lines
#: of 7 pt numbers at single-column width.
_D = 1.15
_CENTRES = {k: (_D / np.sqrt(3) * np.cos(a), _D / np.sqrt(3) * np.sin(a))
            for k, a in zip(PARTNERS, np.deg2rad([90.0, 210.0, 330.0]))}


def region_anchor(pattern, n_grid=401):
    """Where to print a region's numbers: the mean of a grid over the region.

    The mean of a fine grid is a stable label position for every region of a
    symmetric three-circle diagram, including the central one, without
    hand-placed offsets that break when the spacing changes."""
    xs = np.linspace(-2.0, 2.0, n_grid)
    gx, gy = np.meshgrid(xs, xs)
    mask = np.ones_like(gx, dtype=bool)
    for k, inside in zip(PARTNERS, pattern):
        cx, cy = _CENTRES[k]
        within = (gx - cx) ** 2 + (gy - cy) ** 2 <= 1.0
        mask &= within if inside else ~within
    return float(gx[mask].mean()), float(gy[mask].mean())


def draw(counts, set_sizes, quartet_note, outdir, width, height_mm):
    """The Venn panel. Numbers only, no areas: circles are not to scale."""
    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    ax.set_position([0.0, 0.10, 1.0, 0.86])
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(-2.05, 2.05)
    ax.set_ylim(-1.75, 2.15)

    # Neutral circles: the panel carries region colour in its numbers, and the
    # house rule keeps species colours out of panels that carry region colour.
    for k in PARTNERS:
        ax.add_patch(Circle(_CENTRES[k], 1.0, facecolor=NEUTRAL['ink'], alpha=0.04,
                            edgecolor='none', zorder=1))
        ax.add_patch(Circle(_CENTRES[k], 1.0, facecolor='none',
                            edgecolor=NEUTRAL['ink'], linewidth=LINEWIDTH['thin'],
                            zorder=2))

    for _, row in counts.iterrows():
        pattern = tuple(bool(row[f'in_{k}']) for k in PARTNERS)
        x, y = region_anchor(pattern)
        ax.text(x, y + 0.07, f"{row['retained']:,}", ha='center', va='bottom',
                fontsize=FONT['annotation'], color=REGION['retained'], zorder=3)
        ax.text(x, y - 0.07, f"{row['eliminated']:,}", ha='center', va='top',
                fontsize=FONT['annotation'], color=REGION['eliminated'], zorder=3)

    # Set labels outside each circle: the species name as a text of its own (so
    # the house pass sets it in italic) and the number of pairs beside it. The
    # top label sits above the circle with the pair count between name and
    # circle; the two lower labels sit below their circles' outer sides.
    top_y = _CENTRES['parascaris'][1] + 1.0          # top of the upper circle
    ax.text(0.0, top_y + 0.22, SPECIES_LABEL['parascaris'], ha='center',
            va='bottom', fontsize=FONT['axis_label'], color=NEUTRAL['ink'])
    ax.text(0.0, top_y + 0.05, f"{set_sizes['parascaris']:,} pairs", ha='center',
            va='bottom', fontsize=FONT['annotation'], color=NEUTRAL['ink'])
    for k, x in (('baylisascaris', -1.75), ('toxocara', 1.75)):
        ax.text(x, -1.36, SPECIES_LABEL[k], ha='center', va='bottom',
                fontsize=FONT['axis_label'], color=NEUTRAL['ink'])
        ax.text(x, -1.44, f'{set_sizes[k]:,} pairs', ha='center', va='top',
                fontsize=FONT['annotation'], color=NEUTRAL['ink'])

    # Key in the upper left corner, clear of every circle; the quartet note on
    # its own line under the diagram.
    fig.text(0.03, 0.975, 'Retained', ha='left', va='top',
             fontsize=FONT['legend'], color=REGION['retained'])
    fig.text(0.03, 0.935, 'Eliminated', ha='left', va='top',
             fontsize=FONT['legend'], color=REGION['eliminated'])
    fig.text(0.5, 0.03, quartet_note, ha='center', va='bottom',
             fontsize=FONT['legend'], color=NEUTRAL['ink'])

    problems = audit(fig)
    paths = save(fig, outdir, 'gene_set_venn')
    return paths, problems


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Venn diagram of the Ascaris genes in the pairwise rate tables '
                    '(Figure 1), with retained and eliminated counts per region.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--gff', required=True,
                   help='Ascaris gene models (miniprot GFF, as used by annotate_regions.py)')
    p.add_argument('--retained-bed', required=True, help='Ascaris retained regions')
    p.add_argument('--eliminated-bed', required=True, help='Ascaris eliminated regions')
    p.add_argument('--chrom-sizes', default=None,
                   help='.fai or two-column sizes; region calls do not depend on it')
    p.add_argument('--terminal-tol', type=int, default=1000,
                   help='as annotate_regions.py; affects CBR distances, not region calls')
    p.add_argument('--cbr-window', type=float, default=100000.0,
                   help='as annotate_regions.py; not used for region calls')
    p.add_argument('--set', dest='sets', action='append', required=True,
                   metavar='KEY=PATH:COLUMN',
                   help=f'one rate table per partner, KEY in {", ".join(PARTNERS)}, '
                        'COLUMN the Ascaris gene; give all three')
    p.add_argument('--quartets', default=None, metavar='PATH:COLUMN',
                   help='quartet table and its Ascaris gene column, counted against '
                        'the centre')
    p.add_argument('--annotated', default=None,
                   help='annotated rate table whose region calls must be reproduced')
    p.add_argument('--annotated-gene-col', default='gene1')
    p.add_argument('--annotated-region-col', default='region1')
    p.add_argument('--width', choices=['single', 'onehalf', 'double'], default='single')
    p.add_argument('--height-mm', type=float, default=70.0)
    p.add_argument('--outdir', required=True)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'gene_set_venn')
    apply_style()
    logging.info('gene_set_venn.py: %s', ' '.join(sys.argv[1:]))

    # --- inputs --------------------------------------------------------------
    specs = [parse_set_spec(s) for s in args.sets]
    keys = [k for k, _, _ in specs]
    if sorted(keys) != sorted(PARTNERS):
        sys.exit(f'ERROR: give exactly one --set for each of {", ".join(PARTNERS)} '
                 f'(got {", ".join(keys)})')
    sets = {k: read_ids(path, col, f'{k} table') for k, path, col in specs}
    set_sizes = {k: len(v) for k, v in sets.items()}

    region, region_diag = region_calls(args.gff, args.retained_bed,
                                       args.eliminated_bed, args.chrom_sizes,
                                       args.terminal_tol, args.cbr_window)
    region = region.to_dict()

    union = set().union(*sets.values())
    missing = sorted(g for g in union if g not in region)
    if missing:
        sys.exit(f'ERROR: {len(missing):,} Ascaris identifier(s) in the rate tables are '
                 f'not models in {args.gff} (e.g. {missing[0]}). The tables and the GFF '
                 'must come from the same gene set.')

    # --- cross-check against the published region calls -----------------------
    check = None
    if args.annotated:
        ann = pd.read_csv(args.annotated, sep='\t', dtype=str,
                          usecols=[args.annotated_gene_col, args.annotated_region_col])
        ann = ann.dropna().drop_duplicates(args.annotated_gene_col)
        shared = ann[ann[args.annotated_gene_col].isin(region.keys())]
        ours = shared[args.annotated_gene_col].map(region)
        bad = shared[ours != shared[args.annotated_region_col]]
        check = {'table': args.annotated, 'genes_compared': int(len(shared)),
                 'disagreements': int(len(bad)),
                 'retained': int((shared[args.annotated_region_col] == 'retained').sum()),
                 'eliminated': int((shared[args.annotated_region_col] == 'eliminated').sum())}
        logging.info('Region calls checked against %s: %s genes, %s disagreements',
                     os.path.basename(args.annotated), f'{len(shared):,}', len(bad))
        if len(bad):
            g = bad.iloc[0]
            sys.exit(f'ERROR: {len(bad):,} region calls differ from {args.annotated} '
                     f'(e.g. {g[args.annotated_gene_col]}: {region[g[args.annotated_gene_col]]} '
                     f'here, {g[args.annotated_region_col]} there). Use the GFF and BEDs '
                     'that annotated table was made from.')

    # --- counts ---------------------------------------------------------------
    counts = venn_counts(sets, region)
    centre_ids = set.intersection(*sets.values())

    quartet = None
    if args.quartets:
        qpath, qcol = parse_table_spec(args.quartets, '--quartets')
        q = read_ids(qpath, qcol, 'quartet table')
        q_missing = sorted(g for g in q if g not in region)
        if q_missing:
            sys.exit(f'ERROR: {len(q_missing):,} quartet identifier(s) are not models in '
                     f'{args.gff} (e.g. {q_missing[0]}).')
        calls = pd.Series([collapse(region[g]) for g in q], dtype=object)
        quartet = {'table': qpath, 'quartets': len(q),
                   'retained': int((calls == 'retained').sum()),
                   'eliminated': int((calls == 'eliminated').sum()),
                   'in_centre': len(q & centre_ids),
                   'outside_centre': len(q - centre_ids)}
        if quartet['outside_centre']:
            logging.warning('%s of %s quartets fall outside the centre: their pairwise '
                            'estimates did not pass the rate filters in every table.',
                            quartet['outside_centre'], quartet['quartets'])

    if quartet is None:
        quartet_note = ''
    elif quartet['outside_centre'] == 0:
        quartet_note = (f"Branch-model quartets: {quartet['quartets']:,}, "
                        'all in the center')
    else:
        quartet_note = (f"Branch-model quartets: {quartet['quartets']:,} "
                        f"({quartet['in_centre']:,} in the center)")

    # --- outputs --------------------------------------------------------------
    counts_path = os.path.join(args.outdir, 'gene_set_venn_counts.tsv')
    counts.to_csv(counts_path, sep='\t', index=False)
    paths, problems = draw(counts, set_sizes, quartet_note, args.outdir,
                           args.width, args.height_mm)

    lines = ['Ascaris genes in the pairwise rate tables (pairs that passed the rate filters)',
             '', 'Set sizes']
    for k in PARTNERS:
        n_ret = sum(collapse(region[g]) == 'retained' for g in sets[k])
        n_eli = sum(collapse(region[g]) == 'eliminated' for g in sets[k])
        lines.append(f'  {SPECIES_LABEL[k]:<14} {set_sizes[k]:>7,} genes '
                     f'({n_ret:,} retained, {n_eli:,} eliminated, '
                     f'{set_sizes[k] - n_ret - n_eli:,} other)')
    lines += ['', f'Union {len(union):,}; in all three {len(centre_ids):,}', '',
              'Regions (genes = retained + eliminated + other)']
    for _, r in counts.iterrows():
        lines.append(f"  {r['region']:<42} {r['genes']:>7,} = {r['retained']:,} + "
                     f"{r['eliminated']:,} + {r['other']:,}")
    if quartet:
        lines += ['', f"Quartets ({quartet['table']}): {quartet['quartets']:,} "
                      f"({quartet['retained']:,} retained, {quartet['eliminated']:,} "
                      f"eliminated); {quartet['in_centre']:,} in the centre, "
                      f"{quartet['outside_centre']:,} outside"]
    if check:
        lines += ['', f"Region calls reproduce {check['table']}: {check['genes_compared']:,} "
                      f"genes ({check['retained']:,} retained, {check['eliminated']:,} "
                      f"eliminated), {check['disagreements']} disagreements"]
    lines += ['', 'Figure audit: ' + ('clean' if not problems else '; '.join(problems))]
    summary_path = os.path.join(args.outdir, 'gene_set_venn_summary.txt')
    with open(summary_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

    prov = {'script': os.path.abspath(__file__), 'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(), 'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv, 'args': vars(args),
            'inputs': {**{'gff': fingerprint(args.gff),
                          'retained_bed': fingerprint(args.retained_bed),
                          'eliminated_bed': fingerprint(args.eliminated_bed),
                          'chrom_sizes': fingerprint(args.chrom_sizes)},
                       **{f'set_{k}': fingerprint(path) for k, path, _ in specs},
                       'quartets': fingerprint(parse_table_spec(args.quartets, '--quartets')[0])
                       if args.quartets else None,
                       'annotated': fingerprint(args.annotated)},
            'set_sizes': set_sizes, 'union': len(union), 'centre': len(centre_ids),
            'regions': counts.to_dict(orient='records'), 'quartets': quartet,
            'annotated_check': check, 'region_calls': region_diag,
            'outputs': paths + [counts_path, summary_path], 'audit': problems}
    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump(prov, fh, indent=2, default=str)
    logging.info('Done: %s', args.outdir)


if __name__ == '__main__':
    main()
