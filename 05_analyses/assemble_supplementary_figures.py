#!/usr/bin/env python3
"""
assemble_supplementary_figures.py -- copy every supplementary figure panel of
the dN/dS paper into one folder, named by its final figure number.

The analysis scripts name their outputs with the frozen figure-map numbers
(figS1_composition.svg, S8A_omega_above_one.svg, S18_min_s4_sweep.svg, ...).
Since 2026-10-12 the supplementary figures are numbered by first citation in
the text, so those names no longer match the manuscript; two of them (the
window sweep and the reference-bias figure) never did. This script leaves every
analysis output where it is and copies each panel, with its PDF and PNG
siblings, to

    <outdir>/FigS<nn>_<panel>.<ext>

from a manifest that holds the numbering in one place:

    supplementary_figure_sources.tsv     figure<TAB>was<TAB>panel<TAB>path

`figure` is the final number (S1, S2, ...), `was` the frozen-map number it
replaces, `panel` a short name for the file (letters, digits, '-', '_', '.'),
and `path` the panel's SVG (or PDF where no SVG exists). Relative paths resolve
against --base-dir; ${VAR} expands from the environment (the wrapper exports
POPGEN_DIR). Renumbering later means editing the `figure` column and rerunning.


Checks
------
Every problem is reported together, and nothing is copied unless all pass
(exit 1):
  * every source file exists and is not empty;
  * the figure numbers run S1..Sn with no gap (a vacated number must be closed
    up, not left empty);
  * every figure has one `was` number, and no `was` number belongs to two
    figures;
  * no two rows produce the same output name, and no source is listed twice.
After copying, each copy's SHA-256 is compared with its source's.


Outputs (in --outdir)
---------------------
  FigS<nn>_<panel>.<svg|pdf|png>       the copies (formats present at source)
  figure_manifest.tsv                  one row per copied file: figure, was,
                                       panel, format, source, copy, bytes, SHA-256
  supplementary_figures_summary.txt    old -> new numbers, files per figure,
                                       formats missing at source
  assemble_supplementary_figures.log
  provenance.json                      argv, manifest hash, versions


Usage
-----
  python assemble_supplementary_figures.py \\
    --manifest 05_analyses/supplementary_figure_sources.tsv \\
    --base-dir $PROJECT_DIR \\
    --outdir results/supplementary_figures_<stamp>

  # check the manifest against the files without copying anything
  python assemble_supplementary_figures.py ... --check-only

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

SCRIPT_VERSION = '2026-10-12.1'
FORMATS = ('svg', 'pdf', 'png')
COLUMNS = ('figure', 'was', 'panel', 'path')
PANEL_OK = re.compile(r'^[A-Za-z0-9_.-]+$')
FIG_RE = re.compile(r'^S(\d+)$')


def sha256(path):
    """Full SHA-256 of a file, read in 1 MB chunks."""
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_manifest(path, base_dir):
    """Rows of the manifest as dicts, with `src` the absolute path of the panel.

    Comment lines (#) and blank lines are skipped. The first non-comment line
    must be the header figure<TAB>was<TAB>panel<TAB>path. Variable expansion and
    path resolution follow make_supplementary_tables.py.
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
            if len(parts) != len(COLUMNS) or not all(parts):
                problems.append(f'  {path}:{ln}: expected four non-empty fields, got {s!r}')
                continue
            row = dict(zip(COLUMNS, parts), line=ln)
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


def check_rows(rows):
    """Every manifest-level problem, as a list of messages (empty if none)."""
    problems = []
    nums, was_of, fig_of_was = set(), {}, {}
    names, srcs = {}, {}
    for r in rows:
        ln = r['line']
        m, mw = FIG_RE.match(r['figure']), FIG_RE.match(r['was'])
        if not m:
            problems.append(f"  line {ln}: figure {r['figure']!r} is not S<number>")
            continue
        if not mw:
            problems.append(f"  line {ln}: was {r['was']!r} is not S<number>")
        if not PANEL_OK.match(r['panel']):
            problems.append(f"  line {ln}: panel {r['panel']!r} has characters other than "
                            "letters, digits, '-', '_' and '.'")
        nums.add(int(m.group(1)))
        if r['figure'] in was_of and was_of[r['figure']] != r['was']:
            problems.append(f"  line {ln}: {r['figure']} given as was {r['was']} here and "
                            f"{was_of[r['figure']]} earlier")
        was_of.setdefault(r['figure'], r['was'])
        if r['was'] in fig_of_was and fig_of_was[r['was']] != r['figure']:
            problems.append(f"  line {ln}: was {r['was']} assigned to {r['figure']} here and "
                            f"to {fig_of_was[r['was']]} earlier")
        fig_of_was.setdefault(r['was'], r['figure'])
        name = out_stem(r)
        if name in names:
            problems.append(f'  line {ln}: output name {name} already used on line {names[name]}')
        names.setdefault(name, ln)
        if r['src'] in srcs:
            problems.append(f"  line {ln}: source listed twice (also line {srcs[r['src']]}): "
                            f"{r['src']}")
        srcs.setdefault(r['src'], ln)
        if not os.path.isfile(r['src']):
            problems.append(f"  line {ln} ({r['figure']} {r['panel']}): missing: {r['src']}")
        elif os.path.getsize(r['src']) == 0:
            problems.append(f"  line {ln} ({r['figure']} {r['panel']}): empty: {r['src']}")
        elif os.path.splitext(r['src'])[1].lower().lstrip('.') not in FORMATS:
            problems.append(f"  line {ln}: path must be an .svg, .pdf or .png file: {r['src']}")
    if nums:
        gaps = sorted(set(range(1, max(nums) + 1)) - nums)
        if gaps:
            problems.append('  figure numbers have gaps (close them up): '
                            + ', '.join(f'S{g}' for g in gaps))
    return problems


def out_stem(r):
    """FigS<nn>_<panel>: zero-padded so that the folder sorts in figure order."""
    return f"FigS{int(FIG_RE.match(r['figure']).group(1)):02d}_{r['panel']}"


def siblings(src):
    """The panel's files in each format present beside it, keyed by format."""
    stem = os.path.splitext(src)[0]
    return {f: f'{stem}.{f}' for f in FORMATS
            if os.path.isfile(f'{stem}.{f}') and os.path.getsize(f'{stem}.{f}') > 0}


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description='Copy the supplementary figure panels into one folder under their '
                    'final numbers (first citation in the text).')
    ap.add_argument('--manifest', required=True,
                    help='supplementary_figure_sources.tsv (figure, was, panel, path)')
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
                 '(each run writes a fresh folder so that copies from two numberings '
                 'never mix)')
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'assemble_supplementary_figures')
    logging.info('assemble_supplementary_figures.py %s', SCRIPT_VERSION)
    base_dir = os.path.abspath(args.base_dir)
    logging.info('manifest: %s', os.path.abspath(args.manifest))
    logging.info('base dir: %s', base_dir)

    rows = read_manifest(args.manifest, base_dir)
    problems = check_rows(rows)
    figs = sorted({r['figure'] for r in rows}, key=lambda f: int(f[1:]))
    logging.info('%d panels in %d figures (S1-S%d)', len(rows), len(figs),
                 max((int(f[1:]) for f in figs), default=0))
    if problems:
        for p in problems:
            logging.error(p.strip())
        sys.exit(f'ERROR: {len(problems)} problem(s) in the manifest or the files; '
                 'nothing copied.')
    logging.info('manifest checks passed: every source present, numbers S1-S%d with no '
                 'gap, one old number per figure, unique output names', len(figs))
    if args.check_only:
        logging.info('--check-only: nothing copied')
        return 0

    out_rows, missing_fmt = [], []
    for r in sorted(rows, key=lambda r: (int(r['figure'][1:]), r['line'])):
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
            out_rows.append({'figure': r['figure'], 'was': r['was'], 'panel': r['panel'],
                             'format': fmt, 'source': src,
                             'copy': os.path.basename(dest),
                             'bytes': os.path.getsize(dest), 'sha256': h_dst})
        logging.info('  %-4s (was %-4s) %-28s %s', r['figure'], r['was'], r['panel'],
                     ', '.join(files))

    with open(os.path.join(args.outdir, 'figure_manifest.tsv'), 'w', encoding='utf-8') as fh:
        cols = ['figure', 'was', 'panel', 'format', 'source', 'copy', 'bytes', 'sha256']
        fh.write('\t'.join(cols) + '\n')
        for o in out_rows:
            fh.write('\t'.join(str(o[c]) for c in cols) + '\n')

    # --- summary -------------------------------------------------------------------
    was_of = {r['figure']: r['was'] for r in rows}
    n_panels = {f: sum(1 for r in rows if r['figure'] == f) for f in figs}
    L = [f'assemble_supplementary_figures.py {SCRIPT_VERSION}   '
         f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
         f'Output: {os.path.abspath(args.outdir)}', '',
         f'{len(out_rows)} files copied ({len(rows)} panels, {len(figs)} figures); every copy '
         'matches its source (SHA-256).', '',
         'Figure  was    panels  files']
    for f in figs:
        L.append(f'{f:<7} {was_of[f]:<6} {n_panels[f]:>6}  '
                 + ', '.join(o['copy'] for o in out_rows if o['figure'] == f and
                             o['format'] == 'svg'))
    vac = sorted(set(range(1, max(int(w[1:]) for w in was_of.values()) + 1))
                 - {int(w[1:]) for w in was_of.values()})
    if vac:
        L += ['', 'Old numbers not carried forward (vacated): '
              + ', '.join(f'S{v}' for v in vac)]
    if missing_fmt:
        L += ['', 'Formats absent at source (copied what exists):']
        L += [f'  {f} {p}: no {", ".join(lack)}' for f, p, lack in missing_fmt]
    txt = '\n'.join(L) + '\n'
    with open(os.path.join(args.outdir, 'supplementary_figures_summary.txt'), 'w',
              encoding='utf-8') as fh:
        fh.write(txt)
    prov = {'script': os.path.abspath(__file__), 'version': SCRIPT_VERSION,
            'git_commit': git_commit(os.path.dirname(os.path.abspath(__file__))),
            'argv': sys.argv, 'manifest': os.path.abspath(args.manifest),
            'manifest_sha256': sha256(args.manifest), 'base_dir': base_dir,
            'python': platform.python_version(),
            'utc': datetime.now(timezone.utc).isoformat(),
            'numbering': {f: was_of[f] for f in figs}, 'files': len(out_rows)}
    with open(os.path.join(args.outdir, 'provenance.json'), 'w', encoding='utf-8') as fh:
        json.dump(prov, fh, indent=2)
    logging.info('copied %d files into %s', len(out_rows), args.outdir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
