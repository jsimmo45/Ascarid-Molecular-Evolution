#!/usr/bin/env python3
"""
Build the supplementary rate-heterogeneity table.

Two panels, one table:

  Panel A  Pairwise median Ks against the distance the four-way branch model
           predicts. Five of six pairs are reproduced exactly; the deviation
           falls on the pair --fit exact drops, in the direction saturation
           predicts. This is the additivity check.

  Panel B  Tajima relative-rate ratios at codon positions 1+2, position 3, and
           all sites. Three columns rather than three tables, because the point
           is that the ratio does not move between them -- which is what rules
           out saturation as the explanation.

The global-clock LRT is deliberately NOT given rows. At 8-12 M sites its P
values carry no information (any departure from the clock is significant), and
tabulating them implies otherwise. It goes in the caption as one line, which
this script also writes.

INPUTS
------
Panel A reads median Ks straight from the six kaks_*/kaks_filtered.tsv tables
rather than from any summary, so the numbers cannot drift from the source.
Branch lengths come from --branch, defaulting to the published --fit exact
values.

Panel B parses clock_test_sites*/clock_test_summary.txt. Those files carry n1,
n2 and the ratio; the accompanying results.tsv does not.

USAGE
-----
    python make_rate_table.py --outdir supp_rate_table

    # override a median if a table lives somewhere else
    python make_rate_table.py --ks 'As x Pu=0.1957'

OUTPUT
------
    supp_table_rates.tsv   tab-separated, for import into the supplement
    supp_table_rates.md    same content, readable / pasteable
    caption.txt            the LRT one-liner plus a draft legend
"""

import argparse
import os
import re
import sys

import pandas as pd

# Terminal branch lengths, four-way model, --fit exact.
DEFAULT_BRANCHES = {'Ascaris': 0.0461, 'Parascaris': 0.1496,
                    'Baylisascaris': 0.1461, 'Toxocara': 1.2666}
DEFAULT_INTERNAL = 0.0777

# Topology: {Ascaris, Parascaris} | {Baylisascaris, Toxocara}
CLADE_A = {'Ascaris', 'Parascaris'}

# (label, species pair, default kaks_filtered.tsv). Order is the table order.
PAIRS = [
    ('Ascaris x Parascaris',       ('Ascaris', 'Parascaris'),
     'kaks_results/kaks_filtered.tsv'),
    ('Ascaris x Baylisascaris',    ('Ascaris', 'Baylisascaris'),
     'kaks_As_Bsc_v2/kaks_filtered.tsv'),
    ('Parascaris x Baylisascaris', ('Parascaris', 'Baylisascaris'),
     'kaks_Pu_Bsc_v2/kaks_filtered.tsv'),
    ('Ascaris x Toxocara',         ('Ascaris', 'Toxocara'),
     'kaks_Tox_As/kaks_filtered.tsv'),
    ('Parascaris x Toxocara',      ('Parascaris', 'Toxocara'),
     'kaks_Tox_Pu/kaks_filtered.tsv'),
    ('Baylisascaris x Toxocara',   ('Baylisascaris', 'Toxocara'),
     'kaks_Tox_Bsc_v2/kaks_filtered.tsv'),
]

# Site subsets, in the column order Panel B uses.
SUBSETS = [('12', 'codon pos 1+2'), ('3', 'codon pos 3'), ('all', 'all sites')]


def predicted_distance(sp1, sp2, branches, internal):
    """
    Path length between two tips under the fixed topology.

    Same clade -> the two terminal branches. Different clades -> the terminals
    plus the internal branch, which is crossed exactly once.
    """
    d = branches[sp1] + branches[sp2]
    if (sp1 in CLADE_A) != (sp2 in CLADE_A):
        d += internal
    return d


def median_ks(path, col_candidates=('Ks', 'ks', 'dS', 'Ks_yn00')):
    """
    Median Ks from a kaks_filtered.tsv.

    Returns (median, n, column_used), or (nan, 0, None) if the file is absent so
    one missing pair does not abort the whole table.
    """
    if not os.path.isfile(path):
        return float('nan'), 0, None
    df = pd.read_csv(path, sep='\t')
    for c in col_candidates:
        if c in df.columns:
            v = pd.to_numeric(df[c], errors='coerce').dropna()
            v = v[v >= 0]
            return (float(v.median()), int(len(v)), c) if len(v) else \
                (float('nan'), 0, c)
    return float('nan'), 0, None


# Matches e.g.
#   Ascaris vs Parascaris (out: Baylisascaris)   64159  154962   0.41  <1e-300  <1e-300  ****
RRT_RE = re.compile(
    r'^\s*(?P<label>\S.*?\(out:\s*\w+\))\s+'
    r'(?P<n1>\d+)\s+(?P<n2>\d+)\s+(?P<ratio>[\d.]+)\s+'
    r'(?P<p>\S+)\s+(?P<q>\S+)')
LRT_RE = re.compile(r'2\*dlnL\s*=\s*(?P<stat>[\d.eE+-]+)\s+df\s*=\s*(?P<df>\d+)'
                    r'\s+P\s*=\s*(?P<p>\S+)')
SITES_RE = re.compile(r'Supermatrix:\s*([\d,]+)\s*sites')
GENES_RE = re.compile(r'Genes used:\s*([\d,]+)')


def parse_clock_summary(path):
    """
    Pull the LRT line and the Tajima rows out of a clock_test_summary.txt.

    Parsing the summary rather than clock_test_results.tsv is deliberate: only
    the summary carries n1, n2 and the ratio.
    """
    if not os.path.isfile(path):
        return None
    text = open(path).read()
    out = {'rrt': {}, 'lrt': None, 'sites': None, 'genes': None}

    m = LRT_RE.search(text)
    if m:
        out['lrt'] = {'stat': float(m.group('stat')), 'df': int(m.group('df')),
                      'p': m.group('p')}
    m = SITES_RE.search(text)
    if m:
        out['sites'] = int(m.group(1).replace(',', ''))
    m = GENES_RE.search(text)
    if m:
        out['genes'] = int(m.group(1).replace(',', ''))

    for line in text.splitlines():
        m = RRT_RE.match(line)
        if m:
            out['rrt'][m.group('label').strip()] = {
                'n1': int(m.group('n1')), 'n2': int(m.group('n2')),
                'ratio': float(m.group('ratio')), 'p': m.group('p'),
                'q': m.group('q')}
    return out


def parse_overrides(items, cast=float):
    """Parse repeated 'key=value' CLI overrides into a dict."""
    out = {}
    for spec in items or []:
        if '=' not in spec:
            sys.exit(f"ERROR: expected key=value, got {spec!r}")
        k, v = spec.split('=', 1)
        try:
            out[k.strip()] = cast(v)
        except ValueError:
            sys.exit(f"ERROR: {v!r} is not a number (in {spec!r})")
    return out


def main():
    p = argparse.ArgumentParser(
        description='Build the supplementary rate-heterogeneity table.')
    p.add_argument('--clock-dir-pattern', default='clock_test_sites{}',
                   help='Directory pattern per site subset '
                        '(default: clock_test_sites{})')
    p.add_argument('--kaks-dir', default='.',
                   help='Root the default kaks_*/kaks_filtered.tsv paths are '
                        'relative to (default: .)')
    p.add_argument('--ks', action='append', default=None, metavar='PAIR=VALUE',
                   help="Override a median Ks, e.g. --ks 'Ascaris x "
                        "Parascaris=0.1957'. Repeatable.")
    p.add_argument('--branch', action='append', default=None,
                   metavar='SPECIES=LENGTH',
                   help='Override a terminal branch length. Repeatable.')
    p.add_argument('--internal', type=float, default=DEFAULT_INTERNAL,
                   help=f'Internal branch length (default: {DEFAULT_INTERNAL})')
    p.add_argument('--outdir', default='supp_rate_table')
    args = p.parse_args()

    branches = dict(DEFAULT_BRANCHES)
    branches.update(parse_overrides(args.branch))
    ks_override = parse_overrides(args.ks)
    os.makedirs(args.outdir, exist_ok=True)

    # --- Panel A ----------------------------------------------------------
    rows_a = []
    for label, (sp1, sp2), rel in PAIRS:
        path = os.path.join(args.kaks_dir, rel)
        if label in ks_override:
            obs, n, col = ks_override[label], None, 'override'
        else:
            obs, n, col = median_ks(path)
            if pd.isna(obs):
                print(f"  WARNING: no median Ks for {label} ({path})",
                      file=sys.stderr)
        pred = predicted_distance(sp1, sp2, branches, args.internal)
        dev = 100 * (obs - pred) / pred if pred and not pd.isna(obs) else float('nan')
        rows_a.append({'pair': label, 'median_Ks': obs, 'n_orthologues': n,
                       'tree_predicted': pred, 'deviation_pct': dev,
                       'source_column': col})

    # --- Panel B ----------------------------------------------------------
    summaries, labels = {}, []
    for key, _ in SUBSETS:
        d = args.clock_dir_pattern.format(key)
        s = parse_clock_summary(os.path.join(d, 'clock_test_summary.txt'))
        if s is None:
            print(f"  WARNING: missing {d}/clock_test_summary.txt",
                  file=sys.stderr)
            continue
        summaries[key] = s
        for lab in s['rrt']:
            if lab not in labels:
                labels.append(lab)

    if not summaries:
        sys.exit("ERROR: no clock_test_summary.txt found. Check "
                 "--clock-dir-pattern.")

    rows_b = []
    for lab in labels:
        row = {'contrast': lab}
        for key, _ in SUBSETS:
            r = summaries.get(key, {}).get('rrt', {}).get(lab)
            row[f'ratio_{key}'] = r['ratio'] if r else float('nan')
            row[f'n1_{key}'] = r['n1'] if r else None
            row[f'n2_{key}'] = r['n2'] if r else None
        vals = [row[f'ratio_{k}'] for k, _ in SUBSETS
                if not pd.isna(row[f'ratio_{k}'])]
        # Fold difference, expressed so it reads as "N times faster".
        row['fold_faster'] = (1 / min(vals)) if vals and min(vals) > 0 else float('nan')
        rows_b.append(row)

    df_a, df_b = pd.DataFrame(rows_a), pd.DataFrame(rows_b)

    # --- write ------------------------------------------------------------
    tsv = os.path.join(args.outdir, 'supp_table_rates.tsv')
    with open(tsv, 'w') as fh:
        fh.write('# Panel A. Pairwise median Ks vs four-way tree prediction\n')
        df_a.to_csv(fh, sep='\t', index=False, float_format='%.4f')
        fh.write('\n# Panel B. Tajima relative-rate ratios by codon position\n')
        df_b.to_csv(fh, sep='\t', index=False, float_format='%.3f')

    md = os.path.join(args.outdir, 'supp_table_rates.md')
    with open(md, 'w') as fh:
        fh.write('**Panel A. Pairwise synonymous divergence against the '
                 'four-way branch model**\n\n')
        fh.write('| Pair | Median Ks | n | Tree-predicted | Deviation |\n')
        fh.write('|---|---|---|---|---|\n')
        for r in rows_a:
            n = f"{r['n_orthologues']:,}" if r['n_orthologues'] else '--'
            fh.write(f"| {r['pair']} | {r['median_Ks']:.4f} | {n} | "
                     f"{r['tree_predicted']:.4f} | "
                     f"{r['deviation_pct']:+.2f}% |\n")

        fh.write('\n**Panel B. Relative-rate ratios (Tajima 1D)**\n\n')
        fh.write('| Contrast | ' + ' | '.join(t for _, t in SUBSETS)
                 + ' | Fold difference |\n')
        fh.write('|---' * (len(SUBSETS) + 2) + '|\n')
        for r in rows_b:
            cells = ' | '.join(f"{r[f'ratio_{k}']:.2f}" for k, _ in SUBSETS)
            fh.write(f"| {r['contrast']} | {cells} | "
                     f"{r['fold_faster']:.1f}x |\n")

    cap = os.path.join(args.outdir, 'caption.txt')
    with open(cap, 'w') as fh:
        lrt_bits = []
        for key, title in SUBSETS:
            s = summaries.get(key)
            if s and s['lrt']:
                lrt_bits.append(f"{title}: 2dlnL = {s['lrt']['stat']:,.0f}")
        any_s = next(iter(summaries.values()))
        fh.write('Supplementary Table SX. Substitution-rate heterogeneity '
                 'among ascarid lineages.\n\n')
        fh.write(
            'Panel A. Median pairwise synonymous divergence (Ks, yn00) for each '
            'species pair,\nagainst the distance predicted by the four-way '
            'branch model. Five of six pairs are\nreproduced exactly; the '
            'single deviation falls on Parascaris x Toxocara, the pair\n'
            'excluded from the branch-length fit, and lies below prediction as '
            'expected for the\nlongest and most saturated path.\n\n')
        fh.write(
            'Panel B. Tajima 1D relative-rate ratios, computed at codon '
            'positions 1+2, position 3,\nand all sites. Each ratio is the '
            'number of substitutions inferred on the first lineage\ndivided by '
            'the number on the second, counted against the stated outgroup; '
            'equal rates\ngive 1. Ratios are stable across site subsets, so '
            'the asymmetry is not an artefact of\nsaturation at synonymous '
            'sites. Counts are conservative: no multiple-hit correction is\n'
            'applied, which compresses the ratio toward 1 for the faster '
            'lineage.\n\n')
        if lrt_bits:
            fh.write('A global molecular clock is rejected at every site '
                     'subset (likelihood ratio test,\nbaseml GTR+G, df = 2, '
                     'P < 1e-300; ' + '; '.join(lrt_bits) + ').\n'
                     'With ' + f"{any_s['sites']:,}"
                     + ' aligned sites the test has essentially unlimited '
                     'power, so the\nmagnitude of the rate difference in Panel '
                     'B, not the P value, is the informative result.\n\n')
        if any_s.get('genes'):
            fh.write(f"Based on {any_s['genes']:,} four-way orthologous gene "
                     f"alignments.\n")

    print(f"Wrote {tsv}")
    print(f"Wrote {md}")
    print(f"Wrote {cap}")
    print()
    print(open(md).read())


if __name__ == '__main__':
    main()
