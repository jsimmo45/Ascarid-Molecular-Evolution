#!/usr/bin/env python3
"""
make_cai_reference.py - build a CAI reference gene set from expression data.

CAI measures how closely a gene's codon usage matches that of highly expressed
genes, so the reference set has to be defined by expression. Choosing it by
codon bias instead -- for example as the most biased retained genes -- makes the
reference high-CAI by construction and any retained-versus-eliminated comparison
partly circular. This script writes the gene ID list that
`codon_usage_bias.py --ref-genes` expects.

**Which tissues.** The default pools somatic tissues. Somatic cells have already
completed elimination, so a somatic reference cannot include an eliminated gene,
and the resulting w table describes the codon usage the retained genome is
translationally optimised for. Selecting on germline expression instead would
enrich the reference for the same expression class that eliminated genes belong
to, which trades one bias for a subtler one. Use --columns to override. With
--target-table, the reference IDs are counted by --target-region-col, so the
claim is checked on the analysed genes rather than assumed, and
--drop-target-region eliminated removes any that the analysis groups as
eliminated (a model whose workbook gene is retained can still sit where the v50
elimination call differs from the workbook's).

Selection by expression is independent of codon bias either way, which is the
property that removes the circularity.

**Translating workbook genes into model IDs.** The expression workbook names
genes `ascaris_gene_00002`; the rate tables are keyed on v50 miniprot model IDs
(`MP000003`). Use --bridge with `workbook_model_bridge.tsv`
(run_workbook_model_bridge.sh), which links the two by position: each workbook
gene is lifted from v37 to v50 by sequence and each model takes the lifted gene
it overlaps (results draft note 61). Every model assigned to a selected workbook
gene goes into the reference; a workbook gene can have several.

--id-map-gff is the superseded route. It matched the GFF `Target=` names
(`Ascaris_gene_00002_t01`) to the workbook names, but the query proteins come
from a later annotation release that numbers its genes differently, so equal
numbers are different genes (results draft notes 60 and 62). It is kept only to
reproduce the old reference, and warns when used.

Usage
-----
    python make_cai_reference.py \\
      --expression Ascaris_genes_15714_annotation_expression.xlsx \\
      --gene-col name --columns muscle intestine carcass \\
      --top-percent 5 \\
      --bridge results/workbook_model_bridge_<stamp>/workbook_model_bridge.tsv \\
      --target-table results/annotate_regions_20260817/kaks_annotated.tsv \\
      --out cai_reference_ascaris.txt --provenance cai_reference_provenance.json

run_make_cai_reference.sh runs this on Isaac with the inputs from config.sh.
--bridge needs workbook_model_bridge.py (05_analyses) and lib/ (region_lib,
rate_tables) on PYTHONPATH; the wrapper sets both.

Accepts .xlsx, .tsv or .csv. Requires pandas, and openpyxl for .xlsx.
"""

import argparse
import json
import os
import platform
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone

import numpy as np
import pandas as pd


def read_table(path, sheet=0):
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.xlsx', '.xls', '.xlsm'):
        return pd.read_excel(path, sheet_name=sheet)
    return pd.read_csv(path, sep='\t' if ext in ('.tsv', '.tab', '.txt') else ',')


def parse_args():
    p = argparse.ArgumentParser(
        description='Build a CAI reference gene list from expression data.')
    p.add_argument('--expression', required=True,
                   help='Expression table (.xlsx, .tsv or .csv)')
    p.add_argument('--sheet', default=0,
                   help='Worksheet name or index for .xlsx (default: first)')
    p.add_argument('--gene-col', default='name',
                   help='Gene ID column (default: name)')
    p.add_argument('--columns', nargs='+',
                   default=['muscle', 'intestine', 'carcass'],
                   help='Expression columns to pool (default: the somatic '
                        'tissues muscle intestine carcass)')
    p.add_argument('--statistic', choices=['mean', 'median', 'max'], default='mean',
                   help='How to combine the columns (default: mean)')
    p.add_argument('--top-percent', type=float, default=5.0,
                   help='Percentage of genes to keep, highest first (default: 5)')
    p.add_argument('--min-expression', type=float, default=None,
                   help='Optional absolute floor applied before the percentage')
    p.add_argument('--exclude-col', default=None,
                   help="Column marking genes to exclude, e.g. 'eliminated'. "
                        'Rows whose value is not --exclude-keep are dropped.')
    p.add_argument('--exclude-keep', default='.',
                   help="Value in --exclude-col that means 'keep' (default: .)")
    p.add_argument('--target-table', default=None,
                   help='Table the reference IDs have to match, normally the '
                        'annotated rate table. Overlap is checked before '
                        'writing, and if it is poor every other ID column in '
                        'the expression table is tried and the best reported. '
                        'An annotation workbook usually carries several naming '
                        'conventions and picking the wrong one fails silently '
                        'until the analysis is already running.')
    p.add_argument('--target-col', default='gene1',
                   help='ID column in --target-table (default: gene1)')
    p.add_argument('--target-region-col', default='region1',
                   help='Column of --target-table holding the retained/'
                        'eliminated call; the reference IDs are counted by it '
                        'and any eliminated one is listed (default: region1; '
                        'skipped if the table has no such column)')
    p.add_argument('--drop-target-region', default=None, metavar='VALUE',
                   help='Drop reference IDs whose --target-region-col value is '
                        'VALUE, normally eliminated: a model the analysis '
                        'groups as eliminated should not help define the codon '
                        'usage it is compared against. The dropped IDs are '
                        'listed in the report and the provenance.')
    p.add_argument('--probe', action='store_true',
                   help='Report the overlap of every candidate ID column '
                        'against --target-table and exit without writing.')
    m = p.add_mutually_exclusive_group()
    m.add_argument('--bridge', default=None,
                   help='workbook_model_bridge.tsv from '
                        'run_workbook_model_bridge.sh: translates the '
                        'expression gene names (--gene-col) into the v50 model '
                        'IDs the analysis uses, by position. Every model '
                        'assigned to a selected gene is kept.')
    m.add_argument('--id-map-gff', default=None,
                   help='SUPERSEDED (results draft notes 60, 62): translate '
                        'through the miniprot GFF Target= names. The query '
                        'proteins come from a later annotation release, so '
                        'equal numbers are different genes. Kept to reproduce '
                        'the old reference; use --bridge.')
    p.add_argument('--out', help='Output gene ID list (not needed with --probe)')
    p.add_argument('--provenance', default=None,
                   help='Optional JSON recording the inputs (with checksums), '
                        'settings and counts behind --out')
    return p.parse_args()


def normalise_id(text):
    """Fold an identifier to a comparable form.

    The GFF target is ``Ascaris_gene_00002_t01`` and the expression table says
    ``ascaris_gene_00002``: same gene, different capitalisation, and one carries
    a transcript suffix. Folding case and dropping a trailing ``_tNN`` or
    ``.N`` makes them comparable without inventing a mapping.
    """
    t = str(text).strip().lower()
    t = re.sub(r'[_.-]t\d+$', '', t)
    t = re.sub(r'\.\d+$', '', t)
    return t


def bridge_key(text):
    """The key read_bridge() gives a workbook gene name: stripped, lowercased.

    Nothing else is folded. The bridge records the workbook's own names, so a
    suffix rule could only make two different names collide.
    """
    return str(text).strip().lower()


def parse_miniprot_gff(path):
    """Map normalised query gene ID -> set of model IDs, from a miniprot GFF.

    SUPERSEDED: the Target names belong to a later annotation release than the
    expression workbook, so this map pairs different genes (results draft
    notes 60, 62). Kept for --id-map-gff only.

    One query can align more than once (miniprot reports Rank=1, 2, ...), so
    the mapping is one-to-many and every model is kept.
    """
    mapping = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'mRNA':
                continue
            attrs = parts[8]
            m_id = re.search(r'(?:^|;)ID=([^;]+)', attrs)
            m_tg = re.search(r'(?:^|;)Target=(\S+)', attrs)
            if m_id and m_tg:
                mapping[normalise_id(m_tg.group(1))].add(m_id.group(1).strip())
    return mapping


def bridge_map(path):
    """Map workbook gene name -> set of model IDs, from workbook_model_bridge.tsv.

    The file is read by workbook_model_bridge.read_bridge(), which defines the
    format for every script that uses the bridge and refuses an edited or
    malformed file. Models without a workbook gene are not in the map.

    Returns (mapping, info); `info` is read_bridge()'s fingerprint and counts.
    """
    try:
        from workbook_model_bridge import read_bridge
    except ImportError as exc:
        sys.exit('ERROR: --bridge needs workbook_model_bridge.py (05_analyses) '
                 'and lib/ (region_lib, rate_tables) on PYTHONPATH; '
                 f'run_make_cai_reference.sh sets both. ({exc})')
    model_to_name, info = read_bridge(path)
    mapping = defaultdict(set)
    for model, name in model_to_name.items():
        mapping[name].add(model)
    return mapping, info


def split_ids(series):
    """Expand cells holding several comma- or semicolon-separated IDs.

    Returns (set of ids, mapping from row index to its ids). Annotation
    workbooks often list every transcript for a gene in one cell.
    """
    out, per_row = set(), {}
    for idx, val in series.items():
        if pd.isna(val):
            continue
        ids = {t.strip() for t in str(val).replace(';', ',').split(',')
               if t.strip() and t.strip() != '.'}
        if ids:
            per_row[idx] = ids
            out |= ids
    return out, per_row


def candidate_columns(df):
    """Columns that plausibly hold identifiers: non-numeric, mostly distinct.

    Tested by what a column is not, rather than for ``dtype == object``. Pandas
    3 reports text columns as ``str`` and pandas 2 as ``object``, and an
    equality test against one of them silently matches nothing on the other --
    which turns this whole diagnostic into an empty table on exactly the machine
    it is needed on.
    """
    cands = []
    for c in df.columns:
        s = df[c]
        if pd.api.types.is_numeric_dtype(s) or pd.api.types.is_datetime64_any_dtype(s):
            continue
        vals = s.dropna().astype(str)
        vals = vals[vals.str.strip().ne('') & vals.str.strip().ne('.')]
        if len(vals) < 0.2 * len(df):
            continue
        if vals.nunique() < 0.5 * len(vals):
            continue
        # Judge token length after splitting: a cell listing twenty transcripts
        # is still an ID column, and a free-text description is not.
        tokens = vals.head(500).str.replace(';', ',').str.split(',').explode().str.strip()
        if tokens.empty or tokens.str.len().mean() > 60 or tokens.str.contains(' ').mean() > 0.3:
            continue
        cands.append(c)
    return cands


def report_overlap(df, target_ids, columns):
    """Overlap of each candidate column with the target ID space."""
    rows = []
    for c in columns:
        ids, _ = split_ids(df[c])
        hit = len(ids & target_ids)
        rows.append({'column': c, 'ids': len(ids), 'matching_target': hit,
                     'pct_of_target': 100.0 * hit / max(len(target_ids), 1)})
    if not rows:
        return pd.DataFrame(columns=['column', 'ids', 'matching_target',
                                     'pct_of_target'])
    return pd.DataFrame(rows).sort_values('matching_target', ascending=False)


def write_provenance(path, record, inputs):
    """Write `record`, plus checksums of `inputs` and of this script, as JSON.

    Uses region_lib.fingerprint (path, size, mtime, first 16 hex digits of the
    SHA-256), as every other script's provenance does.
    """
    try:
        from region_lib import fingerprint, git_commit
    except ImportError as exc:
        sys.exit('ERROR: --provenance needs lib/ (region_lib) on PYTHONPATH; '
                 f'run_make_cai_reference.sh sets it. ({exc})')
    record = dict(record)
    record['inputs'] = {k: fingerprint(v) for k, v in inputs.items() if v}
    record['script'] = {'path': os.path.abspath(__file__),
                        'fingerprint': fingerprint(os.path.abspath(__file__)),
                        'git_commit': git_commit(os.path.dirname(os.path.abspath(__file__)))}
    record['versions'] = {'python': platform.python_version(),
                          'pandas': pd.__version__, 'numpy': np.__version__}
    with open(path, 'w') as fh:
        json.dump(record, fh, indent=2, default=str)
        fh.write('\n')
    print(f'Provenance -> {path}')


def main():
    args = parse_args()
    if not os.path.isfile(args.expression):
        sys.exit(f'ERROR: not found: {args.expression}')

    prov = {'created_utc': datetime.now(timezone.utc).isoformat(),
            'command': sys.argv, 'settings': vars(args).copy()}

    sheet = args.sheet
    if isinstance(sheet, str) and sheet.isdigit():
        sheet = int(sheet)
    df = read_table(args.expression, sheet)
    print(f'Read {len(df):,} rows from {args.expression}')
    prov['expression_rows'] = int(len(df))

    # --- resolve the ID space before doing any work -------------------------
    target_ids, region_of = set(), {}
    if args.target_table:
        if not os.path.isfile(args.target_table):
            sys.exit(f'ERROR: not found: {args.target_table}')
        tgt = read_table(args.target_table)
        if args.target_col not in tgt.columns:
            sys.exit(f"ERROR: no '{args.target_col}' column in "
                     f'{args.target_table}. Columns: {tgt.columns.tolist()}')
        target_ids = {str(v) for v in tgt[args.target_col].dropna()}
        print(f'Target ID space: {len(target_ids):,} unique '
              f"'{args.target_col}' values from {args.target_table}")
        if args.target_region_col in tgt.columns:
            sub = tgt[[args.target_col, args.target_region_col]].dropna(
                subset=[args.target_col])
            region_of = dict(zip(sub[args.target_col].astype(str),
                                 sub[args.target_region_col].astype(str)))
        else:
            print(f"  (no '{args.target_region_col}' column: reference IDs "
                  'will not be counted by region)')
    if args.drop_target_region and not region_of:
        sys.exit('ERROR: --drop-target-region needs --target-table with a '
                 f"'{args.target_region_col}' column.")

    id_map, key_of = {}, normalise_id
    if args.bridge:
        id_map, binfo = bridge_map(args.bridge)
        if not id_map:
            sys.exit(f'ERROR: {args.bridge} gives no model a workbook gene.')
        key_of = bridge_key
        n_models = sum(len(v) for v in id_map.values())
        print(f'ID map (position bridge): {len(id_map):,} workbook genes -> '
              f'{n_models:,} models, of {binfo["models"]:,} models in '
              f'{args.bridge}')
        prov['bridge'] = binfo
    elif args.id_map_gff:
        if not os.path.isfile(args.id_map_gff):
            sys.exit(f'ERROR: not found: {args.id_map_gff}')
        print('WARNING: --id-map-gff is the superseded name join. The GFF '
              'Target= names come from a later annotation release than the '
              'expression table, so equal numbers are different genes and the '
              'reference written here is not the highly expressed genes '
              '(results draft notes 60, 62). Use --bridge.')
        id_map = parse_miniprot_gff(args.id_map_gff)
        n_models = sum(len(v) for v in id_map.values())
        print(f'ID map: {len(id_map):,} query genes -> {n_models:,} models '
              f'from {args.id_map_gff}')
    if id_map:
        prov['id_map'] = {'route': 'bridge' if args.bridge else 'gff_target_superseded',
                          'genes': int(len(id_map)), 'models': int(n_models)}
        if target_ids:
            mapped = {m for v in id_map.values() for m in v} & target_ids
            print(f'  {len(mapped):,} of {len(target_ids):,} target IDs are '
                  'reachable through the map')
            prov['id_map']['target_ids_reachable'] = int(len(mapped))
            if not mapped:
                sys.exit(f'ERROR: the map\'s model IDs do not match --target-col. '
                         f'The map gives e.g. '
                         f'{sorted(m for v in list(id_map.values())[:1] for m in v)[:2]}, '
                         f'target has e.g. {sorted(target_ids)[:2]}.')

    # With a map in hand the expression IDs are not expected to match the target
    # space directly, so the probe would fire on every run. Skip it.
    if args.probe or (target_ids and not id_map and args.gene_col in df.columns and
                      len(split_ids(df[args.gene_col])[0] & target_ids) <
                      0.2 * len(target_ids)):
        if not target_ids:
            sys.exit('ERROR: --probe needs --target-table.')
        table = report_overlap(df, target_ids, candidate_columns(df))
        if table.empty or table['matching_target'].max() == 0:
            sys.exit('ERROR: no column in the expression table shares an ID '
                     f'convention with {args.target_col} in {args.target_table}.\n'
                     f'  target IDs look like: '
                     f'{", ".join(sorted(target_ids)[:3])}\n'
                     '  The two files need a join key, or an ID mapping.')
        print('\nID column overlap with the target space:')
        print(table.to_string(index=False,
                              formatters={'pct_of_target': '{:.1f}%'.format}))
        best = table.iloc[0]
        if args.probe:
            print(f"\nBest match: --gene-col {best['column']} "
                  f"({best['pct_of_target']:.1f}% of the target IDs)")
            return
        sys.exit(
            f"\nERROR: '{args.gene_col}' matches almost none of the target IDs. "
            f"Rerun with --gene-col {best['column']}, which matches "
            f"{best['matching_target']:,} of {len(target_ids):,} "
            f"({best['pct_of_target']:.1f}%).")

    if args.gene_col not in df.columns:
        sys.exit(f"ERROR: no '{args.gene_col}' column. "
                 f'Columns present: {df.columns.tolist()}')
    if not args.out:
        sys.exit('ERROR: --out is required unless --probe is given.')
    missing = [c for c in args.columns if c not in df.columns]
    if missing:
        sys.exit(f'ERROR: expression columns not found: {missing}. '
                 f'Columns present: {df.columns.tolist()}')

    if args.exclude_col:
        if args.exclude_col not in df.columns:
            sys.exit(f"ERROR: no '{args.exclude_col}' column.")
        before = len(df)
        df = df[df[args.exclude_col].astype(str) == args.exclude_keep]
        print(f'Excluded {before - len(df):,} rows on {args.exclude_col}; '
              f'{len(df):,} remain')

    expr = df[args.columns].apply(pd.to_numeric, errors='coerce')
    combined = getattr(expr, args.statistic)(axis=1)
    df = df.assign(_expr=combined).dropna(subset=['_expr'])

    if args.min_expression is not None:
        before = len(df)
        df = df[df['_expr'] >= args.min_expression]
        print(f'Dropped {before - len(df):,} below --min-expression '
              f'{args.min_expression}; {len(df):,} remain')

    cutoff = np.percentile(df['_expr'], 100 - args.top_percent)
    ref = df[df['_expr'] >= cutoff]
    ids, _ = split_ids(ref[args.gene_col])

    print(f'\nCombined {args.columns} with {args.statistic}')
    print(f'Top {args.top_percent:g}% cutoff: {cutoff:.3f}')
    print(f'  expression range {ref["_expr"].min():.2f} - {ref["_expr"].max():.2f}, '
          f'median {ref["_expr"].median():.2f}')
    print(f'  {len(ids):,} selected gene IDs')
    prov['selection'] = {'genes_ranked': int(len(df)), 'cutoff': float(cutoff),
                         'rows_selected': int(len(ref)),
                         'gene_ids_selected': int(len(ids)),
                         'expression_min': float(ref['_expr'].min()),
                         'expression_median': float(ref['_expr'].median()),
                         'expression_max': float(ref['_expr'].max())}

    if id_map:
        translated, unmapped, several = set(), 0, 0
        for gid in ids:
            hit = id_map.get(key_of(gid))
            if hit:
                translated |= hit
                several += len(hit) > 1
            else:
                unmapped += 1
        print(f'  translated to {len(translated):,} model IDs '
              f'({unmapped:,} of {len(ids):,} selected genes had no model; '
              f'{several:,} had more than one)')
        prov['translation'] = {'models': int(len(translated)),
                               'genes_without_model': int(unmapped),
                               'genes_with_several_models': int(several)}
        ids = translated
    ids = sorted(ids)

    if args.drop_target_region:
        dropped = [i for i in ids if region_of.get(i) == args.drop_target_region]
        ids = [i for i in ids if region_of.get(i) != args.drop_target_region]
        more = f' (and {len(dropped) - 10} more)' if len(dropped) > 10 else ''
        print(f'Dropped {len(dropped)} reference ID(s) that '
              f'{args.target_region_col} calls {args.drop_target_region}'
              + (f": {', '.join(dropped[:10])}{more}" if dropped else ''))
        prov['dropped_by_region'] = {'column': args.target_region_col,
                                     'value': args.drop_target_region,
                                     'ids': dropped}

    if target_ids:
        hit = len(set(ids) & target_ids)
        print(f'Reference IDs found in the target space: {hit:,} of {len(ids):,}')
        prov['in_target'] = int(hit)
        if hit == 0:
            sys.exit('ERROR: none of the reference IDs appear in '
                     f'{args.target_table}. Rerun with --probe to see which '
                     'expression column shares its ID convention.')
        if hit < 30:
            print(f'WARNING: only {hit} reference genes will be usable. A CAI '
                  'w table built on that few genes is unstable; widen '
                  '--top-percent or check the ID column.')
        if region_of:
            by_region = Counter(region_of[i] for i in ids if i in region_of)
            print(f'  by {args.target_region_col}: ' + ', '.join(
                f'{k} {v:,}' for k, v in sorted(by_region.items())))
            prov['in_target_by_region'] = {str(k): int(v) for k, v in by_region.items()}
            elim = [i for i in ids if region_of.get(i) == 'eliminated']
            if elim:
                more = f' (and {len(elim) - 10} more)' if len(elim) > 10 else ''
                print(f'WARNING: {len(elim)} reference ID(s) lie in eliminated '
                      f'DNA by {args.target_region_col}. A somatic reference '
                      'should hold none; check them before using the list: '
                      f'{", ".join(elim[:10])}{more}')

    with open(args.out, 'w') as fh:
        fh.write('\n'.join(ids) + '\n')
    print(f'Reference set: {len(ids):,} IDs -> {args.out}')

    if args.provenance:
        prov['output'] = {'path': os.path.abspath(args.out), 'ids': int(len(ids))}
        write_provenance(args.provenance, prov,
                         {'expression': args.expression, 'bridge': args.bridge,
                          'id_map_gff': args.id_map_gff,
                          'target_table': args.target_table, 'output': args.out})


if __name__ == '__main__':
    main()
