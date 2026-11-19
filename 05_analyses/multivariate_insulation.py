#!/usr/bin/env python3
"""
multivariate_insulation.py - do expression, base composition, compartment or
local gene spacing explain why eliminated genes sit in stronger domain
boundaries, or why they have higher dN/dS?

Figure 6D and 6E, Supplementary Figure S9B, and one row of Table S5 per run.
Replaces `kaks_multivariate.py`.

Writes:

    6D_insulation_by_region.{svg,pdf,png}  6D boxplot, 89 mm (Figure 6, middle row, beside 6C)
    6E_insulation_models.{svg,pdf,png}     6E coefficients of both models, 89 mm (bottom row, beside 6F)
    S9B_feature_correlations.*           Spearman matrix among the model variables, 89 mm
    multivariate_pergene.tsv             one row per gene: every variable, before filtering
    multivariate_models.tsv              every fitted model, one row per term
    multivariate_summary.txt             counts at each filter, models, partial correlations,
                                         EV1 orientation; ends with the Table S5 row
    provenance.json                      inputs, arguments, versions


Models
------
Ordinary least squares. The response and every continuous predictor are
z-scored within the fitted set; `is_eliminated` is 0/1, so its coefficient is
the eliminated-minus-retained difference in standard deviations of the
response. Intervals are beta +/- t(0.975, df) x SE.

    Model 1      dN/dS      ~ insulation + log_expr + gc + is_eliminated
    Model 2      insulation ~ log_expr + gc + is_eliminated
    + IG         either model plus log10 nearest intergenic distance (--gene-structure)
    + EV1        either model plus the compartment eigenvector (--eigenvector)

`log_expr` is log2(expression + 1) at the stage matching the insulation track
(testis with testis Hi-C); `gc` is G+C in a `--gc-window` window centred on the
gene. Figure 6E draws Model 1 and Model 2 side by side -- Results cites 6E for
both -- each without (filled) and with (open) intergenic distance.


What changed from kaks_multivariate.py, and why
-----------------------------------------------
**Positions come from the annotated table**, which `annotate_regions.py` took
from each model's best alignment (lowest Rank, then Identity, then span). The
old script re-parsed the GFF and kept the first record in the file, so a gene
whose secondary alignment was listed first was assigned insulation, GC and EV1
from the wrong locus.

**Expression is joined by position (2026-10-17).** `--bridge` reads
`workbook_model_bridge.tsv` (model ID -> workbook gene, from a sequence lift of
the workbook's v37 genes onto v50). The join it replaces, the model's GFF
`Target=` name, matched each model to a different gene: the query proteins come
from a later annotation release that numbers its genes differently (results
draft note 60). `--name-join` (with `--gff`) keeps it only to reproduce the
superseded runs; one of the two is required.

**6D prints P as plain text (hand edit of 2026-10-15, now in the code)**: no
asterisks beside it, as in 6A and 6C, and the exponent raised as its own 6 pt
text (`plain_text.draw_p`) instead of mathtext.

**Intergenic distance is a new covariate**, from `gene_structure_per_gene.tsv`
(`intergenic_nearest`, measured against every model in the GFF). It enters as
log10(d + 1); overlapping neighbours give d < 0 and are set to 0 first, and the
counts of overlapping, zero and missing distances are reported. **It is a
mediator, not a confounder**: eliminated regions are gene-poor as a property of
what was eliminated, so spacing lies downstream of the process that defines the
groups, and conditioning on it can remove real signal rather than control for
it. Both fits are therefore reported and drawn, and neither replaces the other.

**6E shows both models.** The old forest plot showed Model 1 only, while
Results cites 6E for both. Markers are neutral ink; the old plot coloured a
coefficient red when its interval excluded zero, which in this paper reads as
"eliminated".

**Figure 6 layout (2026-10-18).** 6D and 6E are separate single-column
panels, so that Figure 6 reads in order: the 6A-B strip, then 6C | 6D, then
6E | 6F, with the model schematic last, bottom right. They were one 185 mm
strip (6D narrow, 6E wide), which left 6F beside 6C, out of order. 6D is
`--d-height-mm` tall, which run_multivariate_insulation.sh sets to 6C's height
so the middle row lines up; 6E is `--e-height-mm` (default 48.2 mm, the
height of the 6F schematic beside it). In the narrower 6E the two models share one x label and the
key sits above the term labels in two lines.

**The expression term is named by its stage (2026-10-19).** 6E's row and the
S9B (now S13B) matrix said "Expression", but the term is one stage's TPM
(`--expr-col`, testis with testis Hi-C), while S12A plots the mean of the nine
stages: two variables under one name (results draft note 67). They now read
"Testis expression" and "Testis expr.", from `--expr-col`.

**The Table S5 row stays machine-readable.** `make_robustness_table.sh` takes
the last `is_eliminated` line of each summary. The summary now carries more
models, so it ends with an explicit Table S5 block -- Model 2, or Model 2 +EV1
when an eigenvector is given -- in the legacy line format, which keeps that
rule correct. `multivariate_models.tsv` holds every term of every model for
anything that wants more than one row.

**Partial correlation, both kinds.** The old script reported a Pearson partial
r on raw values; Results calls it rho. Both the Pearson partial and a rank-based
(Spearman) partial are reported, with degrees of freedom n - 2 - k.

**Removed.** `--coord-map`, which mapped germline positions onto somatic
eigenvectors: a somatic eigenvector is NaN over every eliminated gene, so it
cannot control for compartment in a comparison of eliminated with retained
genes. Pre-elimination eigenvectors are already in germline coordinates.
`--coord-check`, a one-off diagnostic of the v50 coordinate rebase.

Kept: per-chromosome EV1 sign correction (`--orient-ev-by-gc`); FAN-C's EV1
sign can invert between chromosomes, and four of 21 were inverted in the testis
eigenvector. Coordinates: table midpoints are 1-based; bedGraph bins are
0-based half-open, and the lookup converts.


Usage
-----
    python multivariate_insulation.py \\
      --kaks results/annotate_regions_20260817/kaks_annotated.tsv \\
      --bridge results/workbook_model_bridge_20261017_103009/workbook_model_bridge.tsv \\
      --insulation .../prepde/AG.v50.teste.is100001.ids80001.insulation.bedGraph \\
      --expression Ascaris_genes_15714_annotation_expression.xlsx --expr-col testis \\
      --genome AG.v50.fa --gc-window 20000 \\
      --gene-structure results/gene_structure_20260817_loci/gene_structure_per_gene.tsv \\
      [--eigenvector .../as_60hr_iced_100kb.matrix.eigenvector --orient-ev-by-gc] \\
      --stage-label 'Testis' --outdir results/multivariate_insulation_..._teste_noEV

Dependencies: numpy, pandas (+ openpyxl), scipy, matplotlib, paper_colors,
plain_text, region_lib, rate_tables, workbook_model_bridge
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import platform
import re
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
from scipy.stats import t as tdist, mannwhitneyu, rankdata, spearmanr
from scipy import __version__ as scipy_version

# A script that cannot import the shared style must stop, not fall back to
# defaults: a figure drawn without it would look plausible and be wrong.
from paper_colors import (REGION, NEUTRAL, FONT, LINEWIDTH, apply_style, audit,
                          save, figure, mm)
from plain_text import draw_p
from region_lib import (parse_gff_models, select_primary_alignment,
                        fingerprint, git_commit, setup_logging)
from rate_tables import read_rate_table
from workbook_model_bridge import read_bridge, shared_gene_counts

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Display names. Kept short: 6E is two coefficient panels sharing one column
#: of labels at journal width, and the caption defines each variable.
PRETTY = {'is_eliminated': 'Eliminated', 'insulation': 'Insulation score',
          'log_expr': 'Expression', 'gc': 'Local GC',
          'log_intergenic': 'Intergenic distance', 'ev1': 'Compartment (EV1)',
          'dN_dS': 'dN/dS'}

#: Shorter still for the S9B matrix, whose labels sit on two sides of a
#: single-column panel.
SHORT = {'dN_dS': 'dN/dS', 'insulation': 'Insulation', 'log_expr': 'Expression',
         'gc': 'Local GC', 'log_intergenic': 'Intergenic', 'ev1': 'EV1'}

#: Row order of 6E, top to bottom.
COEF_ROWS = ['is_eliminated', 'insulation', 'log_expr', 'gc', 'log_intergenic']


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def norm_chrom(c):
    """Strip zero-padding from numeric runs: chr06 -> chr6, chrX1 stays chrX1.

    The table, the bedGraphs, the eigenvectors and the FASTA come from different
    pipelines and do not agree on padding; every source passes through this.
    """
    c = str(c).strip()
    body = c[3:] if c.startswith('chr') else c
    return 'chr' + re.sub(r'(\d+)', lambda m: str(int(m.group(1))), body)


def load_track(path, value_name):
    """A 4-column bedGraph-style track: chrom, start, end, value.

    UCSC `track`/`browser`/`#` lines are skipped; NA, NaN and '.' values are
    dropped. Coordinates stay 0-based half-open, as in the file.
    """
    rows, skipped = [], 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith(('track', 'browser', '#')):
                skipped += 1
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 4:
                skipped += 1
                continue
            rows.append((norm_chrom(p[0]), p[1], p[2], p[3]))
    if not rows:
        sys.exit(f'ERROR: no data rows in {path}')
    df = pd.DataFrame(rows, columns=['chrom', 'start', 'end', value_name])
    # NA, NaN and '.' all become NaN under errors='coerce'. (Not replace() with
    # None: older pandas reads a None replacement as "forward-fill".)
    for c in ('start', 'end', value_name):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    n_all = len(df)
    df = df.dropna(subset=['start', 'end', value_name]).copy()
    df[['start', 'end']] = df[['start', 'end']].astype(np.int64)
    logging.info('  %s: %s bins with a value (%s without; %d header lines)',
                 os.path.basename(path), f'{len(df):,}', f'{n_all - len(df):,}',
                 skipped)
    return df


def assign_from_track(chrom, mid1, track, value_name):
    """Value of the bin containing each gene midpoint, or NaN.

    `mid1` is 1-based (GFF/table convention); the track is 0-based half-open,
    so position p is looked up as p - 1 in [start, end).
    Returns (array, counts dict).
    """
    by_chrom = {c: (g['start'].to_numpy(), g['end'].to_numpy(),
                    g[value_name].to_numpy())
                for c, g in track.sort_values(['chrom', 'start'])
                .groupby('chrom', sort=False)}
    out = np.full(len(chrom), np.nan)
    n_nochrom = n_gap = 0
    for i, (c, m) in enumerate(zip(chrom, mid1)):
        if c not in by_chrom or not np.isfinite(m):
            n_nochrom += 1
            continue
        s, e, v = by_chrom[c]
        p0 = m - 1
        j = int(np.searchsorted(s, p0, side='right')) - 1
        if 0 <= j < len(s) and s[j] <= p0 < e[j]:
            out[i] = v[j]
        else:
            n_gap += 1
    counts = {'assigned': int(np.isfinite(out).sum()),
              'sequence_absent': n_nochrom, 'no_bin': n_gap}
    logging.info('  %s: assigned %s (sequence absent %s, midpoint in no bin %s)',
                 value_name, f"{counts['assigned']:,}", f'{n_nochrom:,}',
                 f'{n_gap:,}')
    return out, counts


def load_genome(path, wanted):
    """FASTA sequences for the wanted (normalised) sequence names only."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None and name in wanted:
                    seqs[name] = ''.join(buf).upper()
                name = norm_chrom(line[1:].split()[0])
                buf = []
            elif name in wanted:
                buf.append(line.strip())
    if name is not None and name in wanted:
        seqs[name] = ''.join(buf).upper()
    missing = sorted(set(wanted) - set(seqs))
    logging.info('  genome: %d of %d needed sequences found%s', len(seqs),
                 len(wanted), f' (missing: {missing[:5]})' if missing else '')
    return seqs


def assign_gc(chrom, mid1, genome, window):
    """G+C fraction in a `window`-bp window centred on each midpoint; N and
    other non-ACGT characters are left out of the denominator."""
    half = window // 2
    out = np.full(len(chrom), np.nan)
    for i, (c, m) in enumerate(zip(chrom, mid1)):
        seq = genome.get(c)
        if seq is None or not np.isfinite(m):
            continue
        c0 = int(m) - 1
        sub = seq[max(0, c0 - half):min(len(seq), c0 + half)]
        gc = sub.count('G') + sub.count('C')
        at = sub.count('A') + sub.count('T')
        if gc + at:
            out[i] = gc / (gc + at)
    logging.info('  gc (window %s bp): assigned %s', f'{window:,}',
                 f'{int(np.isfinite(out).sum()):,}')
    return out


def target_names(gff_path):
    """THE SUPERSEDED NAME JOIN: model ID -> workbook-style name from the best
    alignment's `Target=` (query protein, suffix stripped, lowercased). Wrong:
    the query proteins come from a later annotation release, numbered
    differently (results draft note 60). Only for --name-join.
    Returns (dict, alignment diagnostics)."""
    models = parse_gff_models(gff_path)
    best, diag = select_primary_alignment(models)
    best = best.dropna(subset=['target'])
    names = (best['target'].astype(str)
             .str.replace(r'_t\d+$', '', regex=True).str.lower())
    return dict(zip(best['gene_id'].astype(str), names)), diag


def expression_lookup(name_of, xlsx_path, expr_col):
    """Model ID -> expression at `expr_col`, for the models `name_of` (model ID
    -> lowercased workbook name) maps. Returns (dict, workbook size)."""
    wb = pd.read_excel(xlsx_path)
    for need in ('name', expr_col):
        if need not in wb.columns:
            sys.exit(f"ERROR: workbook lacks column {need!r}; columns: "
                     f"{list(wb.columns)}")
    wb['name'] = wb['name'].astype(str).str.strip().str.lower()
    wb = wb.drop_duplicates('name', keep='first')
    lookup = dict(zip(wb['name'], pd.to_numeric(wb[expr_col], errors='coerce')))
    out = {gid: lookup.get(n, np.nan) for gid, n in name_of.items()}
    return out, len(wb)


def load_intergenic(path, id_col, value_col):
    """Nearest intergenic distance per gene from gene_structure_per_gene.tsv."""
    gs = pd.read_csv(path, sep='\t')
    if id_col not in gs.columns:
        guess = [c for c in ('gene1', 'gene_id', 'gene', 'id') if c in gs.columns]
        if not guess:
            sys.exit(f'ERROR: no gene column in {path}; columns: {list(gs.columns)}')
        id_col = guess[0]
    if value_col not in gs.columns:
        sys.exit(f'ERROR: column {value_col!r} not in {path}; columns: '
                 f'{list(gs.columns)}')
    v = pd.to_numeric(gs[value_col], errors='coerce')
    return dict(zip(gs[id_col].astype(str), v)), id_col


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def zscore(a):
    a = np.asarray(a, float)
    sd = a.std(ddof=0)
    return (a - a.mean()) / sd if sd > 0 else a * 0.0


def ols(X, y):
    """OLS with t-tests. X carries an intercept column. Raises on a
    rank-deficient design, which here means an empty group after filtering."""
    n, k = X.shape
    if np.linalg.matrix_rank(X) < k:
        raise np.linalg.LinAlgError(
            'design matrix is rank-deficient (a predictor is constant or '
            'collinear, usually one region empty after filtering)')
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - k
    mse = (resid @ resid) / dof
    se = np.sqrt(np.diag(mse * np.linalg.inv(X.T @ X)))
    tstat = beta / se
    p = 2 * tdist.sf(np.abs(tstat), dof)
    r2 = 1 - (resid @ resid) / ((y - y.mean()) ** 2).sum()
    return dict(beta=beta, se=se, t=tstat, p=p, r2=r2, dof=dof, n=n)


def fit_model(df, response, cont, binary, label):
    """Standardized OLS of `response` on `cont` (z-scored) and `binary` (0/1).

    A predictor constant within `df` is dropped with a warning rather than
    allowed to make the design singular. Returns a dict or None.
    """
    use_c = [c for c in cont if df[c].std(ddof=0) > 0]
    use_b = [c for c in binary if df[c].nunique() > 1]
    dropped = sorted((set(cont) - set(use_c)) | (set(binary) - set(use_b)))
    if dropped:
        logging.warning('  [%s] constant predictor(s) dropped: %s', label, dropped)
    cols = use_c + use_b
    if not cols:
        logging.warning('  [%s] no usable predictor; not fitted', label)
        return None
    X = np.column_stack([np.ones(len(df))]
                        + [zscore(df[c]) for c in use_c]
                        + [df[c].to_numpy(float) for c in use_b])
    try:
        m = ols(X, zscore(df[response]))
    except np.linalg.LinAlgError as exc:
        logging.warning('  [%s] %s; not fitted', label, exc)
        return None
    q = tdist.ppf(0.975, m['dof'])
    m.update(label=label, response=response, names=['intercept'] + cols,
             ci_lo=m['beta'] - q * m['se'], ci_hi=m['beta'] + q * m['se'],
             n_retained=int((df['is_eliminated'] == 0).sum()),
             n_eliminated=int((df['is_eliminated'] == 1).sum()))
    logging.info('[%s] %s ~ %s  (n = %s, R2 = %.4f)', label, response,
                 ' + '.join(cols), f"{m['n']:,}", m['r2'])
    for nm, b, se, p in zip(m['names'], m['beta'], m['se'], m['p']):
        logging.info('    %-14s beta=%+.4f  se=%.4f  p=%.3e', nm, b, se, p)
    return m


def partial_corr(x, y, Z, rank=False):
    """Partial correlation of x and y given the columns of Z.

    Residualises x and y on Z (with intercept) and correlates the residuals.
    With `rank`, every variable is rank-transformed first, giving the Spearman
    partial correlation. P uses n - 2 - k degrees of freedom.
    """
    x, y, Z = np.asarray(x, float), np.asarray(y, float), np.asarray(Z, float)
    if rank:
        x, y = rankdata(x), rankdata(y)
        Z = np.column_stack([rankdata(Z[:, j]) for j in range(Z.shape[1])])
    Zi = np.column_stack([np.ones(len(x)), Z])
    rx = x - Zi @ np.linalg.lstsq(Zi, x, rcond=None)[0]
    ry = y - Zi @ np.linalg.lstsq(Zi, y, rcond=None)[0]
    r = float(np.corrcoef(rx, ry)[0, 1])
    df = len(x) - 2 - Z.shape[1]
    tstat = r * np.sqrt(df / max(1e-300, 1 - r * r))
    return r, float(2 * tdist.sf(abs(tstat), df)), df


def orient_ev_by_gc(df, chrom_col='chrom', ev_col='ev1', gc_col='gc'):
    """Flip each sequence's EV1 so it correlates positively with GC (A = high GC).

    FAN-C orients by GC but can fail or flip on weak or fragmented sequences.
    Only the sign changes, and only where the sequence's own Spearman(EV1, GC)
    is negative. Returns the list of (sequence, rho before) flipped.
    """
    flipped = []
    sub = df.dropna(subset=[ev_col, gc_col])
    for c in sorted(sub[chrom_col].unique()):
        s = sub[sub[chrom_col] == c]
        if len(s) < 5 or s[ev_col].std(ddof=0) == 0:
            continue
        r = spearmanr(s[ev_col], s[gc_col])[0]
        if r < 0:
            df.loc[df[chrom_col] == c, ev_col] *= -1
            flipped.append((c, float(r)))
    return flipped


def ev_orientation_table(df, chrom_col='chrom', ev_col='ev1', gc_col='gc'):
    """Per-sequence Spearman(EV1, GC) and EV1 spread, after any correction."""
    rows = []
    sub = df.dropna(subset=[ev_col, gc_col])
    for c in sorted(sub[chrom_col].unique()):
        s = sub[sub[chrom_col] == c]
        if len(s) < 5 or s[ev_col].std(ddof=0) == 0:
            continue
        rows.append({'chrom': c, 'n': len(s),
                     'rho_ev_gc': float(spearmanr(s[ev_col], s[gc_col])[0]),
                     'ev1_sd': float(s[ev_col].std())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure helpers (same conventions as functional_enrichment.py)
# ---------------------------------------------------------------------------

def blank_figure(width, height_mm):
    """A journal-width canvas with no axes and no layout engine; panels are
    placed explicitly, because constrained_layout shrinks a panel to nothing
    to make room for text it cannot fit and only warns."""
    fig, ax = figure(width=width, height_mm=height_mm)
    fig.set_layout_engine('none')
    fig.delaxes(ax)
    return fig


def add_axes_mm(fig, left, bottom, right, top):
    w_in, h_in = fig.get_size_inches()
    return fig.add_axes([mm(left) / w_in, mm(bottom) / h_in,
                         (w_in - mm(left) - mm(right)) / w_in,
                         (h_in - mm(bottom) - mm(top)) / h_in])


def style_box(bp, faces, box_alpha):
    """Translucent face (alpha in the colour, not set_alpha, so the black edge
    stays black); box 2, jitter 3, whiskers and caps 4, median 5."""
    for patch, face in zip(bp['boxes'], faces):
        patch.set_facecolor(to_rgba(face, box_alpha))
        patch.set_edgecolor('black')
        patch.set_linewidth(LINEWIDTH['thin'])
        patch.set_zorder(2)
    for part in ('whiskers', 'caps'):
        for art in bp[part]:
            art.set_color('black')
            art.set_linewidth(LINEWIDTH['thin'])
            art.set_zorder(4)
    for med in bp['medians']:
        med.set_color('black')
        med.set_linewidth(LINEWIDTH['data'])
        med.set_zorder(5)


def drawn_texts(fig):
    """Every text object actually drawn (tick labels outside the view limits
    are returned by findobj but never drawn, so they are excluded)."""
    out = list(fig.texts)
    for leg in fig.legends:
        out += list(leg.get_texts())
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            tol = 1e-9 * max(abs(hi - lo), 1e-300)
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                loc = tick.get_loc()
                if loc is None or not (lo - tol <= loc <= hi + tol):
                    continue
                out += [lab for lab in (tick.label1, tick.label2)
                        if lab.get_visible()]
            out.append(axis.label)
            out.append(axis.get_offset_text())
        out += list(ax.texts)
        out += [t for t in (ax.title, getattr(ax, '_left_title', None),
                            getattr(ax, '_right_title', None)) if t is not None]
        leg = ax.get_legend()
        if leg is not None:
            out += list(leg.get_texts())
    return [t for t in out if t.get_visible() and (t.get_text() or '').strip()]


def offcanvas_text(fig, tol_px=1.0):
    """Drawn text running past the canvas edge, which audit() does not check."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    out = []
    for t in drawn_texts(fig):
        bb = t.get_window_extent(rend)
        if (bb.x0 < -tol_px or bb.y0 < -tol_px or bb.x1 > W + tol_px
                or bb.y1 > H + tol_px):
            out.append(f'text {t.get_text()[:28]!r} runs past the canvas edge')
    return out


def fit_to_canvas(fig, ax, pad_mm=0.8, iterations=4, min_frac=0.5):
    """Shrink one axes box until its labels, legend and side text lie on the
    canvas; refuses to go below `min_frac` of its starting size."""
    w_px, h_px = fig.bbox.width, fig.bbox.height
    pad = mm(pad_mm) * fig.dpi
    start = ax.get_position()
    for _ in range(iterations):
        fig.canvas.draw()
        bb = ax.get_tightbbox(fig.canvas.get_renderer())
        pos = ax.get_position()
        dx0 = max(0.0, pad - bb.x0) / w_px
        dx1 = max(0.0, bb.x1 - (w_px - pad)) / w_px
        dy0 = max(0.0, pad - bb.y0) / h_px
        dy1 = max(0.0, bb.y1 - (h_px - pad)) / h_px
        if not (dx0 or dx1 or dy0 or dy1):
            return
        new_w, new_h = pos.width - dx0 - dx1, pos.height - dy0 - dy1
        if new_w < min_frac * start.width or new_h < min_frac * start.height:
            logging.warning('  text around the panel needs more room than the '
                            'canvas allows; not shrinking further')
            return
        ax.set_position([pos.x0 + dx0, pos.y0 + dy0, new_w, new_h])


def finish(fig, outdir, stem, fit=True):
    """Fit, audit, save in every format, close. Returns (paths, problems)."""
    if fit:
        for ax in fig.axes:
            if ax.get_label() != '<colorbar>':
                fit_to_canvas(fig, ax)
    problems = list(audit(fig, verbose=False)) + offcanvas_text(fig)
    for prob in problems:
        logging.warning('  audit %s: %s', stem, prob)
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths, problems


def place_row(fig, axes, shares, gaps_mm, pad_mm=0.8, top_extra_mm=0.0,
              bottom_extra_mm=0.0, iterations=4):
    """Lay out panels side by side from their measured text.

    Each panel's own labels (left), side text (right) and the legend or title
    above are measured, the gaps between neighbours are kept, and the remaining
    width is split in `shares`. The row fills the canvas exactly whatever font
    the node renders with. `top_extra_mm` reserves room above everything for
    a figure-level legend, `bottom_extra_mm` below everything for a shared
    axis label.
    """
    W, H = fig.bbox.width, fig.bbox.height
    pad = mm(pad_mm) * fig.dpi
    extra = mm(top_extra_mm) * fig.dpi
    extra_bottom = mm(bottom_extra_mm) * fig.dpi
    gaps = [mm(g) * fig.dpi for g in gaps_mm]
    shares = np.asarray(shares, float) / np.sum(shares)
    for _ in range(iterations):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        tight = [ax.get_tightbbox(r) for ax in axes]
        boxes = [ax.get_window_extent(r) for ax in axes]
        lefts = [max(0.0, b.x0 - t.x0) for t, b in zip(tight, boxes)]
        rights = [max(0.0, t.x1 - b.x1) for t, b in zip(tight, boxes)]
        top = max(max(0.0, t.y1 - b.y1) for t, b in zip(tight, boxes)) + extra
        bottom = (max(max(0.0, b.y0 - t.y0) for t, b in zip(tight, boxes))
                  + extra_bottom)
        avail = W - 2 * pad - sum(lefts) - sum(rights) - sum(gaps)
        widths = avail * shares
        x = pad
        y0, hh = pad + bottom, H - 2 * pad - bottom - top
        for i, ax in enumerate(axes):
            x += lefts[i]
            ax.set_position([x / W, y0 / H, widths[i] / W, hh / H])
            x += widths[i] + rights[i] + (gaps[i] if i < len(gaps) else 0.0)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def draw_insulation_box(ax, df, rng, stage_label, jitter_max=1500,
                        box_alpha=0.40, cap_pct=(0.5, 99.5), box_width=0.55,
                        x_margin=0.55, jitter_sd=0.06):
    """6D: insulation score by region, retained left, eliminated right.

    `box_width`, `x_margin` (beyond each box centre, in box-position units) and
    `jitter_sd` set the drawing geometry only; a wider panel takes a smaller
    box width relative to its x range, so the boxes keep their printed size.

    Box = quartiles, whiskers = 5th-95th percentile; jitter on top. The y-range
    covers the central `cap_pct` of the pooled values, widened where needed so
    that every whisker ends inside it, plus 3% -- drawing only; the test uses
    every gene.
    Genes beyond the range are drawn at its limit as open triangles pointing
    the way the value lies (the S1 and S6 convention), so a clamped value is
    never read as a real one. Returns (test dict, counts of what was not drawn
    as is).
    """
    groups = ['retained', 'eliminated']
    data = [df.loc[df['region'] == g, 'insulation'].dropna().to_numpy()
            for g in groups]
    bp = ax.boxplot(data, positions=[0, 1], widths=box_width, patch_artist=True,
                    showfliers=False, whis=(5, 95), zorder=2)
    style_box(bp, [REGION[g] for g in groups], box_alpha)
    pooled = np.concatenate(data)
    lo, hi = np.percentile(pooled, cap_pct)
    # Widen the range to every whisker end. The pooled central 99% is set
    # almost entirely by the ~8,000 retained genes, and in the 2026-10-09 build
    # the eliminated 5th percentile lay below it: that whisker ran off the
    # bottom of the axes with no cap and looked as if it stopped at the floor.
    for d in data:
        if len(d):
            w_lo, w_hi = np.percentile(d, (5, 95))
            lo, hi = min(lo, w_lo), max(hi, w_hi)
    # A little room beyond that, so the triangles for clamped genes sit clear
    # of a whisker cap that sets the limit instead of on top of it.
    pad = 0.03 * float(hi - lo)
    lo, hi = float(lo) - pad, float(hi) + pad
    n_skip = 0
    beyond = {}
    for i, (g, d) in enumerate(zip(groups, data)):
        show = d if len(d) <= jitter_max else rng.choice(d, jitter_max,
                                                         replace=False)
        n_skip += len(d) - len(show)
        beyond[g] = {'below': int((d < lo).sum()), 'above': int((d > hi).sum())}
        x = i + rng.normal(0, jitter_sd, len(show))
        inside = (show >= lo) & (show <= hi)
        ax.scatter(x[inside], show[inside], s=1.4, linewidths=0,
                   color=REGION[g], alpha=0.35, zorder=3, rasterized=True)
        # Genes beyond the range: open triangles at the limit. Drawn as dots
        # (the 2026-10-09 build) they made a row on the axis floor and a flat
        # row at the top of the cloud, which sits mid-axis because the axis
        # extends above it for the P bracket, so it read as real values.
        for mask, y_at, marker in ((show < lo, lo, 'v'), (show > hi, hi, '^')):
            if mask.any():
                ax.scatter(x[mask], np.full(int(mask.sum()), y_at),
                           marker=marker, s=4.0, facecolors='none',
                           edgecolors=REGION[g], linewidths=0.4, zorder=4,
                           clip_on=False)
    u, p = mannwhitneyu(data[0], data[1], alternative='two-sided')
    span = hi - lo
    top = hi + 0.06 * span
    ax.plot([0, 0, 1, 1], [top, top + 0.03 * span, top + 0.03 * span, top],
            color='black', lw=0.7, zorder=6)
    # P alone, no asterisks (6A and 6C print P: one convention per figure),
    # and any exponent raised as its own 6 pt text rather than mathtext.
    draw_p(ax, 0.5, top + 0.045 * span, p, ha='center', va='bottom',
           fontsize=FONT['annotation'])
    ax.set_ylim(lo - 0.04 * span, top + 0.22 * span)
    ax.set_xlim(-x_margin, 1 + x_margin)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Retained\n(n = {len(data[0]):,})',
                        f'Eliminated\n(n = {len(data[1]):,})'])
    ax.set_ylabel('Insulation score\n(lower = stronger boundary)')
    if stage_label:
        ax.set_title(f'{stage_label} Hi-C', fontsize=FONT['axis_label'], pad=3)
    test = {'n_retained': len(data[0]), 'n_eliminated': len(data[1]),
            'median_retained': float(np.median(data[0])),
            'median_eliminated': float(np.median(data[1])),
            'U': float(u), 'p': float(p)}
    n_off = sum(b['below'] + b['above'] for b in beyond.values())
    return test, {'jitter_not_drawn': n_skip, 'beyond_y_range': n_off,
                  'beyond_by_region': beyond, 'y_range': [lo, hi],
                  'y_range_pct': list(cap_pct), 'jitter_max': jitter_max}


def nice_ticks(ax, max_ticks, reach=0.12):
    """At most `max_ticks` + 1 round x ticks on a narrow panel, labelled with
    one decimal (more if the step needs it) and a true minus sign.

    A tick that falls just outside the data range (within `reach` of the span)
    is brought in by widening the axis to it, so a panel does not stop a few
    per cent short of a round number and show only two ticks.
    """
    lo, hi = ax.get_xlim()
    ticks = MaxNLocator(nbins=max_ticks, steps=[1, 2, 5, 10]).tick_values(lo, hi)
    span = hi - lo
    above = [t for t in ticks if t > hi]
    below = [t for t in ticks if t < lo]
    if above and above[0] - hi <= reach * span:
        hi = above[0]
    if below and lo - below[-1] <= reach * span:
        lo = below[-1]
    ax.set_xlim(lo, hi)
    shown = [t for t in ticks if lo - 1e-9 * span <= t <= hi + 1e-9 * span]
    step = (ticks[1] - ticks[0]) if len(ticks) > 1 else 1.0
    decimals = max(1, int(-np.floor(np.log10(step))))

    def fmt(v, _pos):
        v = 0.0 if abs(v) < 1e-9 * span else v
        return f'{v:.{decimals}f}'.replace('-', '\u2212')

    ax.xaxis.set_major_locator(FixedLocator(shown))
    ax.xaxis.set_major_formatter(FuncFormatter(fmt))


def draw_coefficients(ax, fits, rows, title, show_labels, xlabel=True,
                      max_ticks=None):
    """6E: standardized coefficients with 95% CIs for one response.

    `fits` is [(model dict, style)] where style is 'filled' or 'open'; a second
    fit is offset vertically so both intervals stay visible. A row whose term
    is not in the model is left empty; the response row is marked as such.
    `xlabel=False` leaves the x label to the caller (one label shared by two
    panels); `max_ticks` caps the major ticks of a narrow panel.
    """
    ys = -np.arange(len(rows), dtype=float)
    offsets = {1: [0.0], 2: [0.17, -0.17]}[len(fits)]
    for (m, style), off in zip(fits, offsets):
        if m is None:
            continue
        for term, y in zip(rows, ys):
            if term not in m['names']:
                continue
            j = m['names'].index(term)
            ax.plot([m['ci_lo'][j], m['ci_hi'][j]], [y + off] * 2,
                    color=NEUTRAL['ink'], lw=LINEWIDTH['thin'],
                    solid_capstyle='butt', zorder=2)
            ax.plot([m['beta'][j]], [y + off], ls='none', marker='o', ms=3.4,
                    mfc=NEUTRAL['ink'] if style == 'filled' else 'white',
                    mec=NEUTRAL['ink'], mew=0.7, zorder=3)
    response = fits[0][0]['response'] if fits and fits[0][0] else None
    for term, y in zip(rows, ys):
        if term == response:
            ax.text(0.5, y, 'response', transform=ax.get_yaxis_transform(),
                    ha='center', va='center_baseline',
                    fontsize=FONT['annotation'], color=NEUTRAL['ink_light'],
                    style='italic')
    ax.axvline(0, color=NEUTRAL['rule'], lw=LINEWIDTH['rule'], ls=(0, (4, 2)),
               zorder=1)
    lo = min([0.0] + [m['ci_lo'][1:].min() for m, _ in fits if m is not None])
    hi = max([0.0] + [m['ci_hi'][1:].max() for m, _ in fits if m is not None])
    pad = 0.08 * (hi - lo if hi > lo else 1.0)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(ys[-1] - 0.6, ys[0] + 0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels([PRETTY[r] for r in rows] if show_labels else [])
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    if max_ticks:
        nice_ticks(ax, max_ticks)
    if xlabel:
        ax.set_xlabel('Standardized coefficient (95% CI)')
    ax.set_title(title, fontsize=FONT['axis_label'], pad=3)


#: Figure 6 geometry (2026-10-18): 6C | 6D in the middle row, 6E | 6F below.
#: The wrapper passes 6C's height for 6D, so these are only the defaults.
D_HEIGHT_MM = 60.0
E_HEIGHT_MM = 48.2


def mm_of(px, fig):
    """Display pixels to millimetres on this figure."""
    return px / fig.dpi * 25.4


def plot_box_panel(core, outdir, stem, rng, stage_label, height_mm):
    """6D on its own, one column wide: the right half of Figure 6's middle
    row, beside 6C. The boxes keep about the printed width they had in the
    strip (13 mm against 10 mm), and the extra width goes to the gap between
    them rather than to fatter boxes."""
    fig = blank_figure('single', height_mm)
    ax = add_axes_mm(fig, 18.0, 11.0, 3.0, 7.0)
    test, d_info = draw_insulation_box(ax, core, rng, stage_label,
                                       box_width=0.45, x_margin=0.65,
                                       jitter_sd=0.07)
    paths, problems = finish(fig, outdir, stem)
    return paths, problems, test, d_info


def plot_coefficient_panel(fits1, fits2, rows, outdir, stem, legend_labels,
                           height_mm, pad_mm=0.8, letter_room_mm=4.0):
    """6E on its own, one column wide: the left half of Figure 6's bottom row,
    beside the 6F schematic.

    The dN/dS model and the insulation model share the term labels on the left
    and one x label below, centred under both. The key (filled = without,
    open = with intergenic distance) sits above in two lines: on one line it
    is wider than the column. Space for both is measured, not guessed, before
    the panels are placed.
    """
    fig = blank_figure('single', height_mm)
    w_in, h_in = fig.get_size_inches()
    ax1 = add_axes_mm(fig, 26.0, 12.0, 34.0, 12.0)
    ax2 = add_axes_mm(fig, 58.0, 12.0, 2.0, 12.0)
    draw_coefficients(ax1, fits1, rows, 'dN/dS model', show_labels=True,
                      xlabel=False, max_ticks=4)
    draw_coefficients(ax2, fits2, rows, 'Insulation model', show_labels=False,
                      xlabel=False, max_ticks=4)
    key_mm = 0.0
    if legend_labels:
        handles = [Line2D([], [], ls='none', marker='o', ms=3.4,
                          mfc=NEUTRAL['ink'] if k == 0 else 'white',
                          mec=NEUTRAL['ink'], mew=0.7)
                   for k in range(len(legend_labels))]
        # Anchored 4 mm in from the left edge: the panel letter goes in the
        # top-left corner, above the term labels.
        leg = fig.legend(handles, legend_labels, loc='upper left',
                         bbox_to_anchor=(mm(pad_mm + letter_room_mm) / w_in,
                                         1 - mm(pad_mm) / h_in),
                         ncol=1, frameon=False, fontsize=FONT['legend'],
                         borderaxespad=0.0, handletextpad=0.3,
                         labelspacing=0.3)
        fig.canvas.draw()
        key_mm = mm_of(leg.get_window_extent(fig.canvas.get_renderer()).height,
                       fig) + 1.2
    xl = fig.text(0.5, mm(pad_mm) / h_in, 'Standardized coefficient (95% CI)',
                  ha='center', va='bottom', fontsize=FONT['axis_label'],
                  color=NEUTRAL['ink'])
    fig.canvas.draw()
    xl_mm = mm_of(xl.get_window_extent(fig.canvas.get_renderer()).height,
                  fig) + 1.0
    place_row(fig, [ax1, ax2], shares=[0.5, 0.5], gaps_mm=[4.0],
              pad_mm=pad_mm, top_extra_mm=key_mm, bottom_extra_mm=xl_mm)
    p1, p2 = ax1.get_position(), ax2.get_position()
    xl.set_x((p1.x0 + p2.x1) / 2)
    paths, problems = finish(fig, outdir, stem, fit=False)
    return paths, problems


def plot_correlations(df, cols, outdir, stem):
    """S9B: Spearman correlations among the model variables.

    Diverging purple-orange, not red-blue: in this paper red and blue mean
    eliminated and retained and nothing else.
    """
    sub = df[cols].dropna()
    rho = sub.corr(method='spearman').to_numpy()
    k = len(cols)
    fig = blank_figure('single', 74.0)
    ax = add_axes_mm(fig, 15, 15, 17, 6)
    im = ax.imshow(rho, cmap='PuOr', vmin=-1, vmax=1, aspect='auto')
    ax.set_xticks(range(k))
    ax.set_yticks(range(k))
    ax.set_xticklabels([SHORT[c] for c in cols], rotation=40, ha='right',
                       rotation_mode='anchor')
    ax.set_yticklabels([SHORT[c] for c in cols])
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    for i in range(k):
        for j in range(k):
            val = f'{rho[i, j]:.2f}'
            val = '0.00' if val == '-0.00' else val.replace('-', '\u2212')
            ax.text(j, i, val, ha='center', va='center',
                    fontsize=FONT['annotation'],
                    color='white' if abs(rho[i, j]) > 0.5 else NEUTRAL['ink'])
    ax.set_title(f'n = {len(sub):,} genes', fontsize=FONT['annotation'], pad=3)
    # Fit the matrix first, then hang the colour bar off its final right edge,
    # so the bar cannot drift away from a matrix that moved.
    fit_to_canvas(fig, ax)
    w_in, h_in = fig.get_size_inches()
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + mm(2.5) / w_in, pos.y0, mm(2.5) / w_in,
                        pos.height], label='<colorbar>')
    cb = fig.colorbar(im, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label('Spearman ρ')
    cb.ax.tick_params(labelsize=FONT['tick'])
    # If the bar's ticks and label run off the right edge, narrow the matrix
    # by the overflow and re-hang the bar.
    W = fig.bbox.width
    for _ in range(3):
        fig.canvas.draw()
        over = cax.get_tightbbox(fig.canvas.get_renderer()).x1 - (W - mm(0.8) * fig.dpi)
        if over <= 0:
            break
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width - over / W, pos.height])
        pos = ax.get_position()
        cax.set_position([pos.x1 + mm(2.5) / w_in, pos.y0, mm(2.5) / w_in,
                          pos.height])
    paths, problems = finish(fig, outdir, stem, fit=False)
    return paths, problems, pd.DataFrame(rho, index=cols, columns=cols), len(sub)


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Insulation, expression, composition, compartment and gene '
                    'spacing as predictors of dN/dS and of boundary strength; '
                    'Figure 6D-E, S9B, Table S5.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--kaks', required=True,
                   help='Annotated rate table (ANNOTATED_TABLE)')
    p.add_argument('--bridge', default=None,
                   help='workbook_model_bridge.tsv (workbook_model_bridge.py): '
                        'model ID -> workbook gene, by position. The join to use')
    p.add_argument('--name-join', action='store_true',
                   help='Join expression by the GFF Target= name instead, as '
                        'before 2026-10-17. Wrong (results draft note 60); only '
                        'to reproduce the superseded runs. Needs --gff')
    p.add_argument('--gff', default=None,
                   help='Ascaris miniprot GFF; used only by --name-join')
    p.add_argument('--insulation', required=True,
                   help='FAN-C insulation bedGraph, pre-elimination stage')
    p.add_argument('--expression', required=True,
                   help='Ascaris_genes_15714_annotation_expression.xlsx')
    p.add_argument('--expr-col', default='testis',
                   help='Workbook column matching the insulation stage '
                        '(testis | ovary | zygote); default testis')
    p.add_argument('--genome', required=True, help='Germline genome FASTA')
    p.add_argument('--gc-window', type=int, default=20000)
    p.add_argument('--gene-structure', default=None,
                   help='gene_structure_per_gene.tsv; adds the intergenic-'
                        'distance models')
    p.add_argument('--gs-gene-col', default='gene1')
    p.add_argument('--intergenic-col', default='intergenic_nearest')
    p.add_argument('--eigenvector', default=None,
                   help='Pre-elimination EV1 track (germline coordinates)')
    p.add_argument('--orient-ev-by-gc', action='store_true',
                   help='Flip each sequence\'s EV1 to correlate positively with '
                        'GC. Sign only.')
    p.add_argument('--ev-label', default=None,
                   help='Name of the eigenvector stage, for the summary and '
                        'Table S5 (e.g. "4-cell (60 h)")')
    p.add_argument('--gene-col', default='gene1')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--chrom-col', default='chrom1')
    p.add_argument('--pos-col', default='midpoint1',
                   help='1-based gene midpoint in the table (default midpoint1; '
                        'computed from start1/end1 if absent)')
    p.add_argument('--stage-label', default=None,
                   help='Insulation stage for the 6D title, e.g. Testis')
    p.add_argument('--d-height-mm', type=float, default=D_HEIGHT_MM,
                   help='Height of 6D in mm. run_multivariate_insulation.sh '
                        "passes 6C's height, so 6C and 6D line up in Figure "
                        f'6\'s middle row (default {D_HEIGHT_MM:g})')
    p.add_argument('--e-height-mm', type=float, default=E_HEIGHT_MM,
                   help='Height of 6E in mm; the default is the height of the '
                        f'6F schematic beside it ({E_HEIGHT_MM:g})')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    for flag, v in (('--d-height-mm', args.d_height_mm),
                    ('--e-height-mm', args.e_height_mm)):
        if not 35.0 <= v <= 120.0:
            sys.exit(f'ERROR: {flag} {v:g} is outside 35-120 mm.')
    if bool(args.bridge) == bool(args.name_join):
        sys.exit('ERROR: give --bridge workbook_model_bridge.tsv (the position '
                 'join), or --name-join to reproduce the superseded name join; '
                 'exactly one of the two.')
    if args.name_join and not args.gff:
        sys.exit('ERROR: --name-join needs --gff.')
    if args.outdir is None:
        args.outdir = os.path.join(
            'results', 'multivariate_insulation_'
            + datetime.now().strftime('%Y%m%d_%H%M%S'))
    return args


def write_block(lines, m, label=None):
    """A model in the legacy summary format, which make_robustness_table.sh
    parses. Do not change the line format."""
    label = label or (m['label'] if m else '?')
    if m is None:
        lines.append(f'[{label}] not fitted (constant predictor or empty group)')
        lines.append('')
        return
    lines.append(f"[{label}]  (standardized betas; n={m['n']}, R2={m['r2']:.4f}; "
                 f"retained {m['n_retained']}, eliminated {m['n_eliminated']})")
    for nm, b, se, p in zip(m['names'], m['beta'], m['se'], m['p']):
        lines.append(f'    {nm:<14} beta={b:+.4f}  se={se:.4f}  p={p:.3e}')
    lines.append('')


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'multivariate_insulation')
    for noisy in ('fontTools', 'fontTools.subset', 'fontTools.ttLib',
                  'matplotlib.font_manager', 'matplotlib.backends.backend_pdf'):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    logging.info('multivariate_insulation.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    apply_style()
    rng = np.random.default_rng(args.seed)
    # The expression term is one stage's TPM, not S12A's nine-stage mean, so
    # the figures name the stage ("Testis expression"); tables keep log_expr.
    stage_name = args.expr_col.replace('_', ' ').strip().capitalize()
    PRETTY['log_expr'] = f'{stage_name} expression'
    SHORT['log_expr'] = f'{stage_name} expr.'
    logging.info('  expression term labelled %r (6E) and %r (S9B)',
                 PRETTY['log_expr'], SHORT['log_expr'])

    # --- table and positions ------------------------------------------------
    tab = read_rate_table(args.kaks, require=('dN_dS',))
    for col in (args.gene_col, args.region_col, args.chrom_col):
        if col not in tab.columns:
            sys.exit(f'ERROR: column {col!r} not in {args.kaks}\n'
                     f'  columns: {list(tab.columns)}')
    if args.pos_col in tab.columns:
        mid = pd.to_numeric(tab[args.pos_col], errors='coerce')
    elif {'start1', 'end1'} <= set(tab.columns):
        mid = (pd.to_numeric(tab['start1'], errors='coerce')
               + pd.to_numeric(tab['end1'], errors='coerce')) // 2
        logging.info('  %s absent; midpoint computed from start1/end1',
                     args.pos_col)
    else:
        sys.exit(f'ERROR: neither {args.pos_col!r} nor start1/end1 in {args.kaks}')
    df = pd.DataFrame({
        'gene': tab[args.gene_col].astype(str),
        'region': tab[args.region_col],
        'chrom': tab[args.chrom_col].map(lambda c: norm_chrom(c)
                                         if pd.notna(c) else np.nan),
        'midpoint': mid.to_numpy(float),
        'dN': pd.to_numeric(tab['dN'], errors='coerce') if 'dN' in tab else np.nan,
        'dS': pd.to_numeric(tab['dS'], errors='coerce') if 'dS' in tab else np.nan,
        'dN_dS': pd.to_numeric(tab['dN_dS'], errors='coerce')})
    steps = {}

    def count(label, mask):
        steps[label] = (int((mask & (df['region'] == 'retained')).sum()),
                        int((mask & (df['region'] == 'eliminated')).sum()))
        logging.info('  %-38s retained %s, eliminated %s', label,
                     f'{steps[label][0]:,}', f'{steps[label][1]:,}')

    logging.info('Table: %s rows from %s', f'{len(df):,}',
                 os.path.basename(args.kaks))
    count('retained or eliminated', df['region'].isin(['retained',
                                                       'eliminated']))
    count('with a table position', df['chrom'].notna() & df['midpoint'].notna())

    # --- per-gene variables -----------------------------------------------
    logging.info('Insulation')
    ins = load_track(args.insulation, 'insulation')
    df['insulation'], ins_counts = assign_from_track(
        df['chrom'].to_numpy(), df['midpoint'].to_numpy(), ins, 'insulation')

    logging.info('Expression (%s)', args.expr_col)
    align_diag, bridge_info = None, None
    if args.bridge:
        name_of, bridge_info = read_bridge(args.bridge)
    else:
        logging.warning('--name-join: the superseded Target-name join, wrong '
                        'for most genes (results draft note 60)')
        name_of, align_diag = target_names(args.gff)
    expr, n_wb = expression_lookup(name_of, args.expression, args.expr_col)
    df['workbook_name'] = df['gene'].map(name_of)
    shared = shared_gene_counts(
        df['workbook_name'].where(df['region'].isin(['retained', 'eliminated'])),
        df['region'])
    df['expr'] = df['gene'].map(expr)
    df['log_expr'] = np.log2(df['expr'] + 1.0)
    logging.info('  expression matched for %s of %s genes',
                 f"{int(df['expr'].notna().sum()):,}", f'{len(df):,}')

    logging.info('GC')
    genome = load_genome(args.genome, set(df['chrom'].dropna()))
    df['gc'] = assign_gc(df['chrom'].to_numpy(), df['midpoint'].to_numpy(),
                         genome, args.gc_window)
    del genome

    ig_counts = None
    if args.gene_structure:
        logging.info('Intergenic distance')
        ig, ig_id = load_intergenic(args.gene_structure, args.gs_gene_col,
                                    args.intergenic_col)
        df['intergenic'] = df['gene'].map(ig)
        raw = df['intergenic']
        ig_counts = {'matched': int(raw.notna().sum()),
                     'missing': int(raw.isna().sum()),
                     'negative_set_to_0': int((raw < 0).sum()),
                     'zero': int((raw == 0).sum()), 'id_column': ig_id}
        df['log_intergenic'] = np.log10(raw.clip(lower=0) + 1.0)
        logging.info('  %s matched, %s missing; %s negative (overlapping '
                     'neighbours) set to 0, %s exactly 0', ig_counts['matched'],
                     ig_counts['missing'], ig_counts['negative_set_to_0'],
                     ig_counts['zero'])

    use_ev = bool(args.eigenvector)
    flipped, ev_table = [], pd.DataFrame()
    if use_ev:
        logging.info('Eigenvector')
        ev = load_track(args.eigenvector, 'ev1')
        df['ev1'], ev_counts = assign_from_track(
            df['chrom'].to_numpy(), df['midpoint'].to_numpy(), ev, 'ev1')
        if args.orient_ev_by_gc:
            flipped = orient_ev_by_gc(df)
            logging.info('  EV1 sign corrected on %d sequence(s)%s', len(flipped),
                         (': ' + ', '.join(f'{c} (rho was {r:+.2f})'
                                           for c, r in flipped)) if flipped else '')
        ev_table = ev_orientation_table(df)

    df['is_eliminated'] = (df['region'] == 'eliminated').astype(float)
    per_gene_cols = ['gene', 'chrom', 'midpoint', 'region', 'dN', 'dS', 'dN_dS',
                     'insulation', 'workbook_name', 'expr', 'log_expr', 'gc'] + \
        (['intergenic', 'log_intergenic'] if args.gene_structure else []) + \
        (['ev1'] if use_ev else [])
    df[per_gene_cols].to_csv(os.path.join(args.outdir, 'multivariate_pergene.tsv'),
                             sep='\t', index=False, float_format='%.6g')

    # --- analysis sets -------------------------------------------------------
    labelled = df['region'].isin(['retained', 'eliminated'])
    count('with insulation', labelled & df['insulation'].notna())
    count('with expression', labelled & df['log_expr'].notna())
    count('with GC', labelled & df['gc'].notna())
    core_mask = (labelled & np.isfinite(df['dN_dS']) & df['insulation'].notna()
                 & df['log_expr'].notna() & df['gc'].notna())
    count('core set (dN/dS, insulation, expression, GC)', core_mask)
    core = df[core_mask].copy()
    if len(core) < 30 or core['is_eliminated'].nunique() < 2:
        sys.exit('ERROR: fewer than 30 genes, or one region empty, after joining. '
                 'Check sequence names across the table, tracks and FASTA, and '
                 'the expression match count above.')
    ig_set = None
    if args.gene_structure:
        count('core set with intergenic distance',
              core_mask & df['log_intergenic'].notna())
        ig_set = core[core['log_intergenic'].notna()].copy()
    ev_set = None
    if use_ev:
        count('core set with EV1', core_mask & df['ev1'].notna())
        ev_set = core[core['ev1'].notna()].copy()
        if ev_set['is_eliminated'].sum() == 0:
            logging.warning('No eliminated gene has an EV1 value: this looks '
                            'like a somatic eigenvector, which cannot control '
                            'for compartment here. Use a pre-elimination one.')

    # --- models ------------------------------------------------------------
    base1, base2 = ['insulation', 'log_expr', 'gc'], ['log_expr', 'gc']
    models = {
        'Model 1': fit_model(core, 'dN_dS', base1, ['is_eliminated'], 'Model 1'),
        'Model 2': fit_model(core, 'insulation', base2, ['is_eliminated'],
                             'Model 2')}
    if ig_set is not None:
        models['Model 1, IG set'] = fit_model(ig_set, 'dN_dS', base1,
                                              ['is_eliminated'], 'Model 1, IG set')
        models['Model 1 +IG'] = fit_model(ig_set, 'dN_dS',
                                          base1 + ['log_intergenic'],
                                          ['is_eliminated'], 'Model 1 +IG')
        models['Model 2, IG set'] = fit_model(ig_set, 'insulation', base2,
                                              ['is_eliminated'], 'Model 2, IG set')
        models['Model 2 +IG'] = fit_model(ig_set, 'insulation',
                                          base2 + ['log_intergenic'],
                                          ['is_eliminated'], 'Model 2 +IG')
    if ev_set is not None:
        models['Model 1 +EV1'] = fit_model(ev_set, 'dN_dS', base1 + ['ev1'],
                                           ['is_eliminated'], 'Model 1 +EV1')
        models['Model 2 +EV1'] = fit_model(ev_set, 'insulation', base2 + ['ev1'],
                                           ['is_eliminated'], 'Model 2 +EV1')
    if models['Model 1'] is None or models['Model 2'] is None:
        sys.exit('ERROR: a primary model could not be fitted; see the log.')

    rows = []
    for key, m in models.items():
        if m is None:
            continue
        for j, nm in enumerate(m['names']):
            rows.append({'model': key, 'response': m['response'], 'term': nm,
                         'beta': m['beta'][j], 'se': m['se'][j], 't': m['t'][j],
                         'p': m['p'][j], 'ci_lo': m['ci_lo'][j],
                         'ci_hi': m['ci_hi'][j], 'n': m['n'],
                         'n_retained': m['n_retained'],
                         'n_eliminated': m['n_eliminated'], 'r2': m['r2'],
                         'df': m['dof']})
    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, 'multivariate_models.tsv'),
                              sep='\t', index=False, float_format='%.6g')

    # --- partial correlations -------------------------------------------------
    partials = {}
    raw_rho, raw_p = spearmanr(core['insulation'], core['dN_dS'])
    Z = core[['log_expr', 'gc']].to_numpy()
    partials['| expression, GC'] = (
        partial_corr(core['insulation'], core['dN_dS'], Z),
        partial_corr(core['insulation'], core['dN_dS'], Z, rank=True), len(core))
    if ig_set is not None:
        Zi = ig_set[['log_expr', 'gc', 'log_intergenic']].to_numpy()
        partials['| expression, GC, intergenic'] = (
            partial_corr(ig_set['insulation'], ig_set['dN_dS'], Zi),
            partial_corr(ig_set['insulation'], ig_set['dN_dS'], Zi, rank=True),
            len(ig_set))
    if ev_set is not None:
        Ze = ev_set[['log_expr', 'gc', 'ev1']].to_numpy()
        partials['| expression, GC, EV1'] = (
            partial_corr(ev_set['insulation'], ev_set['dN_dS'], Ze),
            partial_corr(ev_set['insulation'], ev_set['dN_dS'], Ze, rank=True),
            len(ev_set))

    # --- figures -------------------------------------------------------------
    audits = {}
    rows_e = [r for r in COEF_ROWS
              if r != 'log_intergenic' or ig_set is not None]
    if ig_set is not None:
        fits1 = [(models['Model 1'], 'filled'), (models['Model 1 +IG'], 'open')]
        fits2 = [(models['Model 2'], 'filled'), (models['Model 2 +IG'], 'open')]
        legend = [f"Without intergenic distance (n = {models['Model 1']['n']:,})",
                  f"With intergenic distance (n = {models['Model 1 +IG']['n']:,})"]
    else:
        fits1 = [(models['Model 1'], 'filled')]
        fits2 = [(models['Model 2'], 'filled')]
        legend = []
    logging.info('6D and 6E (single-column panels: 6C | 6D, then 6E | 6F)')
    _, audits['6D'], box_test, d_info = plot_box_panel(
        core, args.outdir, '6D_insulation_by_region', rng, args.stage_label,
        args.d_height_mm)
    _, audits['6E'] = plot_coefficient_panel(
        fits1, fits2, rows_e, args.outdir, '6E_insulation_models', legend,
        args.e_height_mm)

    corr_cols = ['dN_dS', 'insulation', 'log_expr', 'gc'] + \
        (['log_intergenic'] if ig_set is not None else []) + \
        (['ev1'] if use_ev else [])
    logging.info('S9B')
    _, audits['S9B'], rho_mat, n_corr = plot_correlations(
        core, corr_cols, args.outdir, 'S9B_feature_correlations')
    off = rho_mat.where(~np.eye(len(corr_cols), dtype=bool))
    max_abs = float(np.nanmax(np.abs(off.to_numpy())))

    # --- summary -------------------------------------------------------------
    L = ['Insulation, expression, composition, compartment and gene spacing', '',
         f'Table        {os.path.basename(args.kaks)}, {len(df):,} rows '
         f'(positions from the table: best alignment)',
         f'Insulation   {os.path.basename(args.insulation)}'
         + (f'  [{args.stage_label}]' if args.stage_label else ''),
         f'Expression   {os.path.basename(args.expression)} column '
         f"'{args.expr_col}', "
         + (f'via the position bridge {os.path.basename(args.bridge)} '
            f"({bridge_info['assigned']:,} of {bridge_info['models']:,} models "
            f'given a workbook gene)' if args.bridge else
            'via the SUPERSEDED GFF Target= name join (results draft note 60)'),
         f'GC window    {args.gc_window:,} bp, from {os.path.basename(args.genome)}',
         'Intergenic   ' + (os.path.basename(args.gene_structure)
                             if args.gene_structure else 'not supplied'),
         'Eigenvector  ' + (f'{os.path.basename(args.eigenvector)}'
                             + (f'  [{args.ev_label}]' if args.ev_label else '')
                             if use_ev else 'none (EV-free run)'), '',
         'Genes at each step                          retained  eliminated']
    for k, (r, e) in steps.items():
        L.append(f'  {k:<42} {r:>8,}  {e:>10,}')
    L += [f"  insulation lookup: {ins_counts['sequence_absent']:,} genes on a "
          f"sequence absent from the track, {ins_counts['no_bin']:,} in no bin "
          f"(NA bins are dropped on load)",
          f'  sharing their workbook gene with another analysed gene: retained '
          f'{shared[0]:,}, eliminated {shared[1]:,} (each keeps its value)']
    if ig_counts:
        L += [f"  intergenic distance ('{args.intergenic_col}', keyed on "
              f"'{ig_counts['id_column']}'): {ig_counts['matched']:,} matched, "
              f"{ig_counts['missing']:,} missing; "
              f"{ig_counts['negative_set_to_0']:,} negative (overlapping "
              f"neighbours) set to 0, {ig_counts['zero']:,} exactly 0; "
              f'entered as log10(d + 1)']
    L.append('')

    L += ['6D: insulation score by region (core set; Mann-Whitney, two-sided)',
          f"  retained median {box_test['median_retained']:+.4f} "
          f"(n = {box_test['n_retained']:,}); eliminated median "
          f"{box_test['median_eliminated']:+.4f} (n = {box_test['n_eliminated']:,}); "
          f"P = {box_test['p']:.3g}",
          f"  Drawn: y-range {d_info['y_range'][0]:+.3f} to "
          f"{d_info['y_range'][1]:+.3f} = central {d_info['y_range_pct'][0]:g}-"
          f"{d_info['y_range_pct'][1]:g}% of values, widened where needed to "
          f"every whisker end (5th-95th percentile of each region), plus 3%.",
          f"  Beyond it, drawn as open triangles at the limit (counts over all "
          f"genes; only those in the jitter sample are drawn): retained "
          f"{d_info['beyond_by_region']['retained']['below']:,} below / "
          f"{d_info['beyond_by_region']['retained']['above']:,} above, eliminated "
          f"{d_info['beyond_by_region']['eliminated']['below']:,} below / "
          f"{d_info['beyond_by_region']['eliminated']['above']:,} above.",
          f"  {d_info['jitter_not_drawn']:,} jitter points not drawn (at most "
          f"{d_info['jitter_max']:,} per region). Boxes and the test use every "
          f'gene.',
          f'  Panels: 6D 89 x {args.d_height_mm:.1f} mm, 6E 89 x '
          f'{args.e_height_mm:.1f} mm. Figure 6 reads 6A-B, then 6C | 6D, '
          f'then 6E | 6F.', '']

    L += ['PRIMARY MODELS (core set; no EV1, no intergenic distance)',
          '-' * 60]
    write_block(L, models['Model 1'])
    write_block(L, models['Model 2'])
    L += ['[Partial correlation] dN/dS vs insulation',
          f'    raw Spearman rho = {raw_rho:+.4f}, P = {raw_p:.3e} '
          f'(n = {len(core):,})']
    for k, ((r, p, dfree), (rr, pr, _), n) in partials.items():
        L.append(f'    {k:<30} Pearson partial r = {r:+.4f} (P = {p:.3e}); '
                 f'Spearman partial rho = {rr:+.4f} (P = {pr:.3e}); n = {n:,}, '
                 f'df = {dfree:,}')
    L.append('')

    if ig_set is not None:
        L += ['LOCAL GENE SPACING (genes with an intergenic distance)',
              '-' * 60,
              '  Intergenic distance is a mediator, not a confounder: eliminated '
              'regions are gene-poor',
              '  as a property of what was eliminated, so conditioning on it can '
              'remove real signal.',
              '  Both fits are reported; neither replaces the other.', '']
        for k in ('Model 1, IG set', 'Model 1 +IG', 'Model 2, IG set',
                  'Model 2 +IG'):
            write_block(L, models.get(k), k)
        for resp in ('1', '2'):
            a, b = models.get(f'Model {resp}, IG set'), models.get(f'Model {resp} +IG')
            if (a is not None and b is not None and 'is_eliminated' in a['names']
                    and 'is_eliminated' in b['names']):
                ja, jb = (a['names'].index('is_eliminated'),
                          b['names'].index('is_eliminated'))
                L.append(f"  Model {resp}: is_eliminated beta {a['beta'][ja]:+.4f} "
                         f"without -> {b['beta'][jb]:+.4f} with intergenic "
                         f"distance ({100 * (b['beta'][jb] / a['beta'][ja] - 1):+.1f}%)")
        L.append('')

    if use_ev:
        L += ['COMPARTMENT-CONTROLLED (genes with an EV1 value)', '-' * 60]
        write_block(L, models.get('Model 1 +EV1'), 'Model 1 +EV1')
        write_block(L, models.get('Model 2 +EV1'), 'Model 2 +EV1')
        ev_vals = core['ev1'].dropna()
        L += [f'  EV1 spread: sd {ev_vals.std():.4f}, range '
              f'[{ev_vals.min():+.3f}, {ev_vals.max():+.3f}], n = {len(ev_vals):,}. '
              f'A low spread means weak compartments and little control power.',
              f'  Sign corrected on {len(flipped)} sequence(s)'
              + (': ' + ', '.join(f'{c} (rho was {r:+.2f})' for c, r in flipped)
                 if flipped else '')]
        if len(ev_table):
            neg = ev_table[ev_table['rho_ev_gc'] < 0]
            L.append(f"  After correction: Spearman(EV1, GC) per sequence ranges "
                     f"{ev_table['rho_ev_gc'].min():+.2f} to "
                     f"{ev_table['rho_ev_gc'].max():+.2f}; {len(neg)} negative.")
            all_r = spearmanr(core['ev1'], core['gc'], nan_policy='omit')[0]
            L.append(f'  Genome-wide Spearman(EV1, GC) = {all_r:+.3f} (expect '
                     f'positive).')
        L.append('')

    L += [f'S9B: Spearman matrix over {n_corr:,} genes; largest off-diagonal '
          f'|rho| = {max_abs:.2f}.', '']
    L += ['Figure audit (paper_colors.audit, plus text past the canvas edge)']
    for k, probs in audits.items():
        L.append(f'  {k}: ' + ('clean' if not probs else '; '.join(probs)))
    L += ['', 'Notes', '',
          '  Standardized betas: response and continuous predictors z-scored '
          'within the fitted set;',
          '  is_eliminated is 0/1, so its beta is the eliminated-minus-retained '
          'difference in SDs.',
          '  6E intervals: beta +/- t(0.975, df) x SE. Filled = without, open = '
          'with intergenic distance.', '']

    # The Table S5 row. make_robustness_table.sh takes the LAST is_eliminated
    # line in this file, so this block must stay last.
    s5_key = 'Model 2 +EV1' if use_ev else 'Model 2'
    s5 = models.get(s5_key)
    L += ['TABLE S5 ROW: is_eliminated in the insulation model, compartment '
          'specification = '
          + (args.ev_label or os.path.basename(args.eigenvector) if use_ev
             else 'none'),
          f'[{s5_key}]']
    if s5 is not None and 'is_eliminated' in s5['names']:
        j = s5['names'].index('is_eliminated')
        L.append(f"    {'is_eliminated':<14} beta={s5['beta'][j]:+.4f}  "
                 f"se={s5['se'][j]:.4f}  p={s5['p'][j]:.3e}")
    else:
        # No is_eliminated line at all: make_robustness_table.sh would otherwise
        # take an earlier model's line and report it under this specification.
        L.append('    not estimable: no eliminated gene in the fitted set '
                 '(a somatic eigenvector?)')
    summary = '\n'.join(L) + '\n'
    with open(os.path.join(args.outdir, 'multivariate_summary.txt'), 'w') as fh:
        fh.write(summary)

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(), 'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'scipy': scipy_version,
                         'matplotlib': matplotlib.__version__},
            'argv': sys.argv, 'args': vars(args),
            'inputs': {k: fingerprint(v) for k, v in
                       (('kaks', args.kaks), ('bridge', args.bridge),
                        ('gff', args.gff),
                        ('insulation', args.insulation),
                        ('expression', args.expression),
                        ('genome', args.genome),
                        ('gene_structure', args.gene_structure),
                        ('eigenvector', args.eigenvector))},
            'gene_counts': {k: {'retained': r, 'eliminated': e}
                            for k, (r, e) in steps.items()},
            'join': {'mode': 'position bridge' if args.bridge
                     else 'name (superseded)', 'bridge': bridge_info,
                     'alignment': align_diag,
                     'shared_workbook_genes': {'retained': shared[0],
                                               'eliminated': shared[1]}},
            'workbook_genes': n_wb,
            'intergenic_counts': ig_counts,
            'ev1_sign_flipped': flipped,
            'figure_audit': audits,
        }, fh, indent=2, default=str)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
