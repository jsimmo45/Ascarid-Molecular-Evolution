#!/usr/bin/env python3
"""
assemble_main_figures.py -- copy every main-figure panel of the dN/dS paper
into one folder, named by figure and panel.

The companion of assemble_supplementary_figures.py. The analysis scripts write
their panels into their own run folders, under names that predate the final
figure order (fig5A_branch_dnds_as_region.svg, 6AB_expression_stage.svg,
Figure_polymorphism.svg, ...). This script leaves every analysis output where
it is and copies each panel, with its PDF and PNG siblings, to

    <outdir>/Fig<n>_<panel>.<ext>

from a manifest that records which run supplies each panel:

    main_figure_sources.tsv     figure<TAB>panel<TAB>path<TAB>note

`figure` is the figure number (1, 2, ...); `panel` a short name for the file
(letters, digits, '-', '_', '.'), which should start with the panel letter(s);
`path` the panel's SVG (or PDF where no SVG exists); `note` free text, which may
be empty. A path of '-' marks a panel made outside the pipeline, such as an
Illustrator schematic: it is listed in the summary and nothing is copied for
it. Relative paths resolve against --base-dir; ${VAR} expands from the
environment (the wrapper exports POPGEN_DIR). To refresh a panel, point its line
at the newer run and rerun.


Checks
------
Every problem is reported together, and nothing is copied unless all pass
(exit 1):
  * every source file exists, is not empty and is an .svg, .pdf or .png;
  * the figure numbers run 1..n with no gap;
  * no two rows produce the same output name, and no source is listed twice.
After copying, each copy's SHA-256 is compared with its source's.


Outputs (in --outdir)
---------------------
  Fig<n>_<panel>.<svg|pdf|png>     the copies (formats present at source)
  figure_manifest.tsv              one row per copied file: figure, panel,
                                   format, source, copy, bytes, SHA-256
  main_figures_summary.txt         files per figure, panels made outside the
                                   pipeline, formats missing at source
  assemble_main_figures.log
  provenance.json                  argv, manifest hash, versions


Usage
-----
  python assemble_main_figures.py \\
    --manifest 05_analyses/main_figure_sources.tsv \\
    --base-dir $PROJECT_DIR \\
    --outdir results/main_figures_<stamp>

  # check the manifest against the files without copying anything
  python assemble_main_figures.py ... --check-only

Dependencies: the standard library, and region_lib from lib/ (logging,
provenance).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import platform
import re
import shutil
import sys
from datetime import datetime, timezone

from region_lib import git_commit, setup_logging

SCRIPT_VERSION = '2026-10-14.1'
FORMATS = ('svg', 'pdf', 'png')
COLUMNS = ('figure', 'panel', 'path', 'note')
EXTERNAL = '-'                       # path of a panel made outside the pipeline
PANEL_OK = re.compile(r'^[A-Za-z0-9_.-]+$')
FIG_RE = re.compile(r'^(\d+)$')


def sha256(path):
    """Full SHA-256 of a file, read in 1 MB chunks."""
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_manifest(path, base_dir):
    """Rows of the manifest as dicts, with `src` the absolute path of the panel
    (None for a panel made outside the pipeline).

    Comment lines (#) and blank lines are skipped. The first non-comment line
    must be the header figure<TAB>panel<TAB>path<TAB>note. A row may omit the
    note; figure, panel and path must be non-empty.
    """
    rows, header, problems = [], None, []
    with open(path, encoding='utf-8') as fh:
        for ln, line in enumerate(fh, 1):
            s = line.rstrip('\n')
            if not s.strip() or s.lstrip().startswith('#'):
                continue
            parts = [p.strip() for p in s.split('\t')]
            if header is None:
                header = parts
                if tuple(header) != COLUMNS:
                    sys.exit(f'ERROR: {path}:{ln}: header must be '
                             f"{chr(9).join(COLUMNS)!r}, got {s!r}")
                continue
            if len(parts) == len(COLUMNS) - 1:
                parts.append('')                      # note omitted
            if len(parts) != len(COLUMNS) or not all(parts[:3]):
                problems.append(f'  {path}:{ln}: expected figure, panel and path '
                                f'(note optional), tab-separated; got {s!r}')
                continue
            row = dict(zip(COLUMNS, parts), line=ln, src=None)
            if row['path'] != EXTERNAL:
                e = os.path.expanduser(os.path.expandvars(row['path']))
                if '$' in e:
                    problems.append(f"  line {ln}: unexpanded variable in {row['path']!r} "
                                    '(export it, e.g. POPGEN_DIR)')
                    continue
                if not os.path.isabs(e):
                    e = os.path.join(base_dir, e)
                row['src'] = os.path.normpath(e)
            rows.append(row)
    if header is None:
        problems.append(f'  {path}: no header and no rows')
    if problems:
        sys.exit('ERROR: manifest problems:\n' + '\n'.join(problems))
    return rows


def out_stem(r):
    """Fig<n>_<panel>."""
    return f"Fig{int(r['figure'])}_{r['panel']}"


def check_rows(rows):
    """Every manifest-level problem, as a list of messages (empty if none)."""
    problems = []
    nums, names, srcs = set(), {}, {}
    for r in rows:
        ln = r['line']
        m = FIG_RE.match(r['figure'])
        if not m or int(m.group(1)) == 0:
            problems.append(f"  line {ln}: figure {r['figure']!r} is not a number from 1")
            continue
        nums.add(int(m.group(1)))
        if not PANEL_OK.match(r['panel']):
            problems.append(f"  line {ln}: panel {r['panel']!r} has characters other than "
                            "letters, digits, '-', '_' and '.'")
        name = out_stem(r)
        if name in names:
            problems.append(f'  line {ln}: output name {name} already used on line {names[name]}')
        names.setdefault(name, ln)
        if r['src'] is None:
            continue                                  # made outside the pipeline
        if r['src'] in srcs:
            problems.append(f"  line {ln}: source listed twice (also line {srcs[r['src']]}): "
                            f"{r['src']}")
        srcs.setdefault(r['src'], ln)
        if not os.path.isfile(r['src']):
            problems.append(f"  line {ln} (Fig {r['figure']} {r['panel']}): missing: {r['src']}")
        elif os.path.getsize(r['src']) == 0:
            problems.append(f"  line {ln} (Fig {r['figure']} {r['panel']}): empty: {r['src']}")
        elif os.path.splitext(r['src'])[1].lower().lstrip('.') not in FORMATS:
            problems.append(f"  line {ln}: path must be an .svg, .pdf or .png file: {r['src']}")
    if nums:
        gaps = sorted(set(range(1, max(nums) + 1)) - nums)
        if gaps:
            problems.append('  figure numbers have gaps: ' + ', '.join(str(g) for g in gaps))
    return problems


def siblings(src):
    """The panel's files in each format present beside it, keyed by format."""
    stem = os.path.splitext(src)[0]
    return {f: f'{stem}.{f}' for f in FORMATS
            if os.path.isfile(f'{stem}.{f}') and os.path.getsize(f'{stem}.{f}') > 0}


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description='Copy the main-figure panels into one folder, named by figure and panel.')
    ap.add_argument('--manifest', required=True,
                    help='main_figure_sources.tsv (figure, panel, path, note)')
    ap.add_argument('--base-dir', required=True,
                    help='directory that relative manifest paths resolve against '
                         '(PROJECT_DIR, the miniprot directory)')
    ap.add_argument('--outdir', required=True,
                    help='output folder; must not exist or be empty')
    ap.add_argument('--check-only', action='store_true',
                    help='check the manifest and the files, copy nothing')
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if os.path.isdir(args.outdir) and os.listdir(args.outdir):
        sys.exit(f'ERROR: --outdir exists and is not empty: {args.outdir} '
                 '(each run writes a fresh folder so that copies of two builds never mix)')
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'assemble_main_figures')
    logging.info('assemble_main_figures.py %s', SCRIPT_VERSION)
    base_dir = os.path.abspath(args.base_dir)
    logging.info('manifest: %s', os.path.abspath(args.manifest))
    logging.info('base dir: %s', base_dir)

    rows = read_manifest(args.manifest, base_dir)
    problems = check_rows(rows)
    figs = sorted({int(r['figure']) for r in rows if FIG_RE.match(r['figure'])})
    external = [r for r in rows if r['src'] is None]
    logging.info('%d panels in %d figures (1-%d), %d made outside the pipeline',
                 len(rows), len(figs), max(figs, default=0), len(external))
    if problems:
        for p in problems:
            logging.error(p.strip())
        sys.exit(f'ERROR: {len(problems)} problem(s) in the manifest or the files; '
                 'nothing copied.')
    logging.info('manifest checks passed: every source present, figures 1-%d with no gap, '
                 'unique output names', max(figs))
    if args.check_only:
        logging.info('--check-only: nothing copied')
        return 0

    out_rows, missing_fmt = [], []
    for r in sorted(rows, key=lambda r: (int(r['figure']), r['line'])):
        if r['src'] is None:
            logging.info('  Fig %-2s %-28s not copied (%s)', r['figure'], r['panel'],
                         r['note'] or 'made outside the pipeline')
            continue
        files = siblings(r['src'])
        lacking = [f for f in FORMATS if f not in files]
        if lacking:
            missing_fmt.append((r['figure'], r['panel'], lacking))
        for fmt, src in files.items():
            dest = os.path.join(args.outdir, f'{out_stem(r)}.{fmt}')
            shutil.copy2(src, dest)
            h_src, h_dst = sha256(src), sha256(dest)
            if h_src != h_dst:
                sys.exit(f'ERROR: copy differs from its source: {dest}')
            out_rows.append({'figure': r['figure'], 'panel': r['panel'], 'format': fmt,
                             'source': src, 'copy': os.path.basename(dest),
                             'bytes': os.path.getsize(dest), 'sha256': h_dst})
        logging.info('  Fig %-2s %-28s %s', r['figure'], r['panel'], ', '.join(files))

    with open(os.path.join(args.outdir, 'figure_manifest.tsv'), 'w', encoding='utf-8') as fh:
        cols = ['figure', 'panel', 'format', 'source', 'copy', 'bytes', 'sha256']
        fh.write('\t'.join(cols) + '\n')
        for o in out_rows:
            fh.write('\t'.join(str(o[c]) for c in cols) + '\n')

    # --- summary -------------------------------------------------------------------
    n_copied = len({(r['figure'], r['panel']) for r in rows if r['src'] is not None})
    L = [f'assemble_main_figures.py {SCRIPT_VERSION}   '
         f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
         f'Output: {os.path.abspath(args.outdir)}', '',
         f'{len(out_rows)} files copied ({n_copied} panels in {len(figs)} figures); every copy '
         'matches its source (SHA-256).', '']
    for f in figs:
        L.append(f'Fig {f}')
        for r in sorted((r for r in rows if int(r['figure']) == f), key=lambda r: r['line']):
            if r['src'] is None:
                L.append(f"  {out_stem(r):<34} not copied: "
                         f"{r['note'] or 'made outside the pipeline'}")
                continue
            fmts = [o['format'] for o in out_rows
                    if o['figure'] == r['figure'] and o['panel'] == r['panel']]
            L.append(f"  {out_stem(r):<34} {', '.join(fmts)}")
    if external:
        L += ['', 'Panels made outside the pipeline (add them in Illustrator):']
        L += [f"  Fig {r['figure']} {r['panel']}: {r['note'] or '-'}" for r in external]
    if missing_fmt:
        L += ['', 'Formats absent at source (copied what exists):']
        L += [f'  Fig {f} {p}: no {", ".join(lack)}' for f, p, lack in missing_fmt]
    txt = '\n'.join(L) + '\n'
    with open(os.path.join(args.outdir, 'main_figures_summary.txt'), 'w',
              encoding='utf-8') as fh:
        fh.write(txt)
    prov = {'script': os.path.abspath(__file__), 'version': SCRIPT_VERSION,
            'git_commit': git_commit(os.path.dirname(os.path.abspath(__file__))),
            'argv': sys.argv, 'manifest': os.path.abspath(args.manifest),
            'manifest_sha256': sha256(args.manifest), 'base_dir': base_dir,
            'python': platform.python_version(),
            'utc': datetime.now(timezone.utc).isoformat(),
            'panels': [{'figure': r['figure'], 'panel': r['panel'],
                        'source': r['src'] or EXTERNAL, 'note': r['note']}
                       for r in sorted(rows, key=lambda r: (int(r['figure']), r['line']))],
            'files': len(out_rows)}
    with open(os.path.join(args.outdir, 'provenance.json'), 'w', encoding='utf-8') as fh:
        json.dump(prov, fh, indent=2)
    logging.info('copied %d files into %s', len(out_rows), args.outdir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
