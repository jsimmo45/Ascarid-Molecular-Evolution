#!/usr/bin/env python3
"""
workbook_model_bridge.py - link the expression workbook's genes to the v50
miniprot models by position, through a sequence lift from v37 to v50.

Why
---
Figure 6 joined the workbook to the models by name: each model's Target
(`Ascaris_gene_NNNNN_t01`) against the workbook's `ascaris_gene_NNNNN`. The two
come from different annotation releases, each numbering its genes serially along
its own assembly, so equal numbers are different genes (results draft note 60:
the two elimination calls agree barely better than chance, and gene spans are
uncorrelated). The workbook's coordinates are on the v37 germline assembly
(checked 2026-10-17: every workbook chromosome is in the v37 index, no gene runs
past a chromosome end, and 918 of 918 workbook-eliminated genes lie inside v37
eliminated regions against 0 of 14,796 retained genes). This script carries each
workbook gene onto v50 by aligning its v37 sequence, then gives each v50 model
the workbook gene it overlaps. No gene name is used to make the link.

Steps
-----
1. Cut each workbook gene's span (start-end, 1-based inclusive, as the workbook
   records it: len = end - start + 1) out of the v37 FASTA through its .fai.
2. Align the spans to the v50 FASTA with minimap2 (-x asm5 -c --secondary=no).
   A gene's alignment is its best piece plus any other pieces of a split
   alignment that lie beside it: same v50 sequence and strand, colinear, and
   spanning no more than 1.5 x the v37 span + 10 kb on v50 (an insertion in v50,
   or a v37 gap filled in v50, splits an alignment without moving the gene).
   The gene is lifted when that alignment covers at least --min-coverage of the
   span's non-N bases (query bases aligned, so v37 sequence missing from v50
   counts against it even inside one alignment), at --min-identity or more
   (1 - minimap2's gap-compressed divergence, de:f, so a long indel counts as
   one difference and does not drown the base-level identity) and mapping
   quality --min-mapq or more (best piece). Spans more than --max-n-frac N are
   not attempted. Every gene that fails is reported with its reason.
3. Each v50 model (its best alignment, region_lib.select_primary_alignment)
   takes the lifted workbook gene, on the same sequence and strand, with the
   highest reciprocal overlap: (overlap / model span) x (overlap / gene span).
   That favours a gene of the model's own size over a long host gene whose
   intron happens to contain it. The model is assigned if the overlap covers at
   least --min-model-overlap of the model and no other gene scores within
   --ambiguity of the best; otherwise it is left unassigned, with the reason.
   Models overlapped only by a gene on the other strand are counted separately:
   a large number would mean a strand error, not antisense genes.
4. Validation that does not use names. For the analysed genes (--kaks): the
   workbook's own elimination call against region1, and the model span against
   the workbook gene span. The name join of Figure 6 is scored the same way
   beside it, so the two can be compared directly.

Writes (into --outdir)
----------------------
    workbook_model_bridge.tsv     one row per v50 model: gene_id, workbook_name
                                  (empty if unassigned), the reason, and the
                                  evidence (overlap fractions, lift quality)
    workbook_lift.tsv             one row per workbook gene: v37 and v50
                                  coordinates, lift quality, reason
    workbook_genes_v50.bed        lifted workbook genes on v50 (BED6)
    bridge_validation.tsv         region1 x workbook call for both joins
    analysed_bridge_status.tsv    analysed genes: bridge outcome by region1
    workbook_genes_v37.fa         the extracted spans (minimap2 input)
    workbook_genes_v50.paf        minimap2 output; minimap2.log its stderr
    workbook_model_bridge_summary.txt
    provenance.json
    workbook_model_bridge.log

Downstream, functional_enrichment.py, multivariate_insulation.py and
germline_restricted_dnds.py take `--bridge workbook_model_bridge.tsv` and join
through `gene_id -> workbook_name` instead of the Target name. They read the
file with `read_bridge()` below, so the format is defined in one place.

Usage
-----
    python workbook_model_bridge.py \\
      --annotation Ascaris_genes_15714_annotation_expression.xlsx \\
      --old-genome $DB_DIR/ascaris/v37/AG.genome.fa \\
      --new-genome $DB_DIR/ascaris/v50/AG.v50.fa \\
      --gff1 As_miniprot.gff \\
      --kaks results/annotate_regions_20260817/kaks_annotated.tsv \\
      --threads 4 --outdir results/workbook_model_bridge_20261017_100000

Dependencies: numpy, pandas (+ openpyxl), scipy, minimap2 (binary),
region_lib, rate_tables, functional_enrichment (for load_workbook and the name
join it is compared against; imported in main(), so that the scripts importing
read_bridge() from here do not import each other in a circle)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from scipy import __version__ as scipy_version

from region_lib import (parse_gff_models, select_primary_alignment,
                        fingerprint, git_commit, setup_logging)
from rate_tables import read_rate_table

_HERE = os.path.dirname(os.path.abspath(__file__))
STRANDS = ('+', '-')
# Pieces of a split alignment are joined only while the joined alignment spans
# no more than this on v50: 1.5 x the v37 span + 10 kb.
CHAIN_MAX_STRETCH = 1.5
CHAIN_PAD_BP = 10_000


# ---------------------------------------------------------------------------
# FASTA access
# ---------------------------------------------------------------------------

def read_fai(path):
    """{name: (length, offset, line_bases, line_width)} from a samtools .fai."""
    if not os.path.exists(path):
        sys.exit(f'ERROR: no index {path}; run `samtools faidx` on the FASTA.')
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 5:
                out[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return out


def fasta_names(path):
    """Sequence names of a FASTA: from its .fai if there is one, else headers."""
    if os.path.exists(path + '.fai'):
        return set(read_fai(path + '.fai'))
    names = set()
    with open(path, 'rb') as fh:
        for line in fh:
            if line.startswith(b'>'):
                names.add(line[1:].split()[0].decode())
    return names


def fetch(fh, entry, start0, end0):
    """Sequence [start0, end0) of one record, read through its .fai entry.

    The byte position of base i is offset + (i // line_bases) * line_width +
    i % line_bases, so one seek and one read return the span with its line
    breaks, which are then removed (\\r as well, for CRLF files).
    """
    length, offset, lb, lw = entry
    start0, end0 = max(0, start0), min(length, end0)
    if end0 <= start0:
        return ''

    def pos(i):
        return offset + (i // lb) * lw + (i % lb)
    fh.seek(pos(start0))
    raw = fh.read(pos(end0 - 1) - pos(start0) + 1)
    return raw.replace(b'\n', b'').replace(b'\r', b'').decode('ascii').upper()


def extract_spans(wb, fasta, fai, out_fa):
    """Write every workbook gene's v37 span to `out_fa`; return per-gene rows."""
    rows = []
    with open(fasta, 'rb') as fh, open(out_fa, 'w') as out:
        for r in wb.itertuples(index=False):
            rec = {'workbook_name': r.name, 'v37_chr': r.chr,
                   'v37_start': r.start, 'v37_end': r.end,
                   'workbook_strand': r.strand, 'workbook_len': r.len,
                   'span_bp': np.nan, 'n_bases': np.nan, 'extracted': False,
                   'reason': ''}
            if r.chr not in fai:
                rec['reason'] = 'chromosome not in the v37 index'
            elif not (np.isfinite(r.start) and np.isfinite(r.end)) or r.end < r.start:
                rec['reason'] = 'no valid span in the workbook'
            else:
                seq = fetch(fh, fai[r.chr], int(r.start) - 1, int(r.end))
                if not seq:
                    rec['reason'] = 'span outside the v37 chromosome'
                else:
                    rec.update(span_bp=len(seq), n_bases=seq.count('N'),
                               extracted=True)
                    out.write(f'>{r.name}\n')
                    for i in range(0, len(seq), 80):
                        out.write(seq[i:i + 80] + '\n')
            rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# minimap2
# ---------------------------------------------------------------------------

def minimap2_version(exe):
    try:
        return subprocess.run([exe, '--version'], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as exc:
        sys.exit(f'ERROR: cannot run minimap2 ({exe!r}): {exc}\n'
                 f'  Install it into the environment (conda install -c bioconda '
                 f'minimap2) or pass --minimap2 /path/to/minimap2.')


def run_minimap2(exe, preset, threads, target, query, paf_path, log_path):
    """Align `query` to `target`; PAF to `paf_path`, minimap2's stderr to
    `log_path`. Exits on a non-zero return code, showing the end of the log."""
    cmd = [exe, '-x', preset, '-c', '--secondary=no', '-t', str(threads),
           target, query]
    logging.info('  %s', ' '.join(cmd))
    with open(paf_path, 'w') as out, open(log_path, 'w') as err:
        res = subprocess.run(cmd, stdout=out, stderr=err)
    if res.returncode != 0:
        with open(log_path) as fh:
            tail = fh.read()[-2000:]
        sys.exit(f'ERROR: minimap2 exited with {res.returncode}\n{tail}')
    return cmd


PAF_COLS = ['qname', 'qlen', 'qstart', 'qend', 'paf_strand', 'tname', 'tlen',
            'tstart', 'tend', 'nmatch', 'alnlen', 'mapq', 'de', 'ins_bp']
_CIGAR_INS = re.compile(r'(\d+)I')


def read_paf(path):
    """Primary records (tp:A:P, which includes the pieces of a split
    alignment) from a PAF, as a DataFrame; empty if none.

    de   gap-compressed divergence (de:f, written with -c); NaN if absent
    ins_bp  query bases inside the alignment that are absent from the target
            (I operations of the cg:Z CIGAR): v37 sequence missing from v50
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 12:
                continue
            tags = {t[:2]: t[5:] for t in f[12:] if len(t) > 5 and t[2] == ':'}
            if tags.get('tp', 'P') != 'P':
                continue
            de = float(tags['de']) if 'de' in tags else np.nan
            ins = sum(int(n) for n in _CIGAR_INS.findall(tags.get('cg', '')))
            rows.append((f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5],
                         int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10]),
                         int(f[11]), de, ins))
    return pd.DataFrame(rows, columns=PAF_COLS)


def union_length(intervals):
    """Bases covered by half-open intervals [(start, end), ...]."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def chain_pieces(recs, span_bp):
    """The alignment of one gene: its best piece and the colinear pieces beside it.

    `recs` are the gene's primary PAF records as dicts. Starting from the piece
    with most matching bases, pieces on the same v50 sequence and strand are
    added, in order of matching bases, when the chain stays colinear (query
    order = target order, reversed on the minus strand) and its v50 extent stays
    within CHAIN_MAX_STRETCH x span + CHAIN_PAD_BP.
    """
    recs = sorted(recs, key=lambda r: -r['nmatch'])
    anchor = recs[0]
    max_extent = CHAIN_MAX_STRETCH * span_bp + CHAIN_PAD_BP
    chain = [anchor]
    for r in recs[1:]:
        if r['tname'] != anchor['tname'] or r['paf_strand'] != anchor['paf_strand']:
            continue
        trial = chain + [r]
        if (max(x['tend'] for x in trial) - min(x['tstart'] for x in trial)
                > max_extent):
            continue
        t = [x['tstart'] for x in sorted(trial, key=lambda x: x['qstart'])]
        steps = np.diff(t)
        if (anchor['paf_strand'] == '+' and (steps <= 0).any()) or \
           (anchor['paf_strand'] == '-' and (steps >= 0).any()):
            continue
        chain = trial
    return chain


def lift_genes(spans, paf, min_identity, min_coverage, min_mapq, max_n_frac):
    """One row per workbook gene: where it lies on v50, and whether it lifted.

    coverage = query bases aligned by the chain (the union of its query
    intervals less the bases inserted relative to v50) / non-N bases of the
    span, at most 1, so a v37 gap does not count against a gene but v37
    sequence that v50 lacks does. identity = 1 - de, averaged over the chain's
    pieces weighted by matching bases (matching bases / alignment columns if
    the PAF has no de tag). mapq = the best piece's. len_ratio = v50 extent /
    v37 span, > 1 where v50 carries sequence v37 lacked.
    """
    by_gene = {q: g.to_dict('records') for q, g in paf.groupby('qname', sort=False)}
    flip = {'+': '-', '-': '+'}
    rows = []
    for s in spans.itertuples(index=False):
        rec = s._asdict()
        rec.update(v50_chr='', v50_start=np.nan, v50_end=np.nan, paf_strand='',
                   v50_strand='.', identity=np.nan, coverage=np.nan,
                   mapq=np.nan, n_pieces=0, n_pieces_chained=0, lift_mode='',
                   len_ratio=np.nan, lifted=False, reason_lift='')
        if not s.extracted:
            rec['reason_lift'] = s.reason
            rows.append(rec)
            continue
        n_frac = s.n_bases / s.span_bp
        if n_frac > max_n_frac:
            rec['reason_lift'] = f'span more than {max_n_frac:g} N'
            rows.append(rec)
            continue
        recs = by_gene.get(s.workbook_name)
        if not recs:
            rec['reason_lift'] = 'no alignment to v50'
            rows.append(rec)
            continue
        chain = chain_pieces(recs, s.span_bp)
        anchor = chain[0]
        covered = max(0, union_length([(x['qstart'], x['qend']) for x in chain])
                      - sum(x['ins_bp'] for x in chain))
        non_n = max(1.0, s.span_bp - s.n_bases)
        w = np.array([x['nmatch'] for x in chain], float)
        de = np.array([x['de'] for x in chain], float)
        if np.isfinite(de).all() and w.sum() > 0:
            identity = 1.0 - float((de * w).sum() / w.sum())
        else:
            identity = (sum(x['nmatch'] for x in chain)
                        / max(1, sum(x['alnlen'] for x in chain)))
        lo = min(x['tstart'] for x in chain)
        hi = max(x['tend'] for x in chain)
        ps = anchor['paf_strand']
        rec.update(
            v50_chr=anchor['tname'], v50_start=lo + 1, v50_end=hi,
            paf_strand=ps,
            # The gene's v50 strand: its v37 strand, flipped if the span
            # aligned in reverse.
            v50_strand=(s.workbook_strand if ps == '+' else
                        flip.get(str(s.workbook_strand), '.')),
            identity=identity, coverage=min(1.0, covered / non_n),
            mapq=anchor['mapq'],
            n_pieces=len(recs), n_pieces_chained=len(chain),
            lift_mode='single' if len(chain) == 1 else 'chained',
            len_ratio=(hi - lo) / s.span_bp)
        if rec['coverage'] < min_coverage:
            rec['reason_lift'] = f'coverage below {min_coverage:g}'
        elif rec['identity'] < min_identity:
            rec['reason_lift'] = f'identity below {min_identity:g}'
        elif rec['mapq'] < min_mapq:
            rec['reason_lift'] = f'mapping quality below {min_mapq} (repeat)'
        else:
            rec['lifted'] = True
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Models <- lifted workbook genes
# ---------------------------------------------------------------------------

def assign_models(models, lifted, min_model_overlap, ambiguity):
    """One row per model: the lifted workbook gene it takes, or the reason not.

    Within each sequence every model is compared with every lifted gene (at most
    a couple of thousand of each), in an overlap matrix zeroed where both
    strands are known and differ. Score = (overlap / model span) x (overlap /
    gene span). The strand-blind matrix is kept only to count models whose
    sole overlapping gene is on the other strand.
    """
    rows = []
    g_by_chr = {c: g.reset_index(drop=True) for c, g in lifted.groupby('v50_chr')}
    for c, m in models.groupby('chrom', sort=False):
        m = m.reset_index(drop=True)
        ms, me = m['start'].to_numpy(np.int64), m['end'].to_numpy(np.int64)
        mlen = (me - ms + 1).astype(float)
        mst = m['strand'].astype(str).to_numpy()
        base = [{'gene_id': gid, 'model_chrom': c, 'model_start': int(ms[i]),
                 'model_end': int(me[i]), 'model_strand': mst[i],
                 'workbook_name': np.nan, 'candidate': np.nan,
                 'n_candidates': 0, 'overlap_bp': 0.0, 'frac_model': 0.0,
                 'frac_gene': 0.0, 'score': 0.0, 'second_score': 0.0,
                 'reason': ''} for i, gid in enumerate(m['gene_id'])]
        g = g_by_chr.get(c)
        if g is None or g.empty:
            for rec in base:
                rec['reason'] = 'no lifted workbook gene on this sequence'
            rows.extend(base)
            continue
        gs = g['v50_start'].to_numpy(np.int64)
        ge = g['v50_end'].to_numpy(np.int64)
        glen = (ge - gs + 1).astype(float)
        gst = g['v50_strand'].astype(str).to_numpy()
        ov_any = (np.minimum(me[:, None], ge[None, :])
                  - np.maximum(ms[:, None], gs[None, :]) + 1).clip(min=0)
        ov_any = ov_any.astype(float)
        known = np.isin(mst, STRANDS)[:, None] & np.isin(gst, STRANDS)[None, :]
        ov = np.where(known & (mst[:, None] != gst[None, :]), 0.0, ov_any)
        score = (ov / mlen[:, None]) * (ov / glen[None, :])
        order = np.argsort(-score, axis=1, kind='stable')
        idx = np.arange(len(m))
        best = order[:, 0]
        s1 = score[idx, best]
        s2 = score[idx, order[:, 1]] if score.shape[1] > 1 else np.zeros(len(m))
        n_cand = (ov > 0).sum(axis=1)
        any_ov = (ov_any > 0).any(axis=1)
        for i, rec in enumerate(base):
            j = best[i]
            if ov[i, j] <= 0:
                rec['reason'] = ('only an opposite-strand workbook gene overlaps it'
                                 if any_ov[i] else
                                 'no lifted workbook gene overlaps it')
                continue
            rec.update(candidate=g['workbook_name'].iat[j],
                       n_candidates=int(n_cand[i]), overlap_bp=float(ov[i, j]),
                       frac_model=float(ov[i, j] / mlen[i]),
                       frac_gene=float(ov[i, j] / glen[j]),
                       score=float(s1[i]), second_score=float(s2[i]))
            if rec['frac_model'] < min_model_overlap:
                rec['reason'] = f'overlap below {min_model_overlap:g} of the model'
            elif s2[i] >= ambiguity * s1[i]:
                rec['reason'] = 'ambiguous: two workbook genes overlap it alike'
            else:
                rec['workbook_name'] = rec['candidate']
                rec['reason'] = 'assigned'
        rows.extend(base)
    out = pd.DataFrame(rows)
    shared = out.dropna(subset=['workbook_name']).groupby('workbook_name').size()
    out['models_sharing_gene'] = out['workbook_name'].map(shared)
    return out


# ---------------------------------------------------------------------------
# Reading the bridge (for the scripts that join expression through it)
# ---------------------------------------------------------------------------

def read_bridge(path):
    """Model ID -> workbook gene name, from workbook_model_bridge.tsv.

    Only models given a workbook gene are in the map; any other model has no
    entry and joins to nothing, as an unmatched gene always has. Names are
    stripped and lowercased, the key functional_enrichment.load_workbook()
    gives the workbook. Read as text throughout, so no identifier is ever
    turned into a number or a missing value.

    Returns (map, info). `info` holds the file's fingerprint and counts, for a
    summary and provenance.
    """
    if not path or not os.path.exists(path):
        sys.exit(f'ERROR: bridge file {path!r} not found; run '
                 f'run_workbook_model_bridge.sh first.')
    b = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False,
                    na_values=[''])
    missing = [c for c in ('gene_id', 'workbook_name', 'reason')
               if c not in b.columns]
    if missing:
        sys.exit(f'ERROR: {path} lacks column(s) {missing}; is it '
                 f'workbook_model_bridge.tsv?')
    b['gene_id'] = b['gene_id'].str.strip()
    dup = b['gene_id'].duplicated()
    if dup.any():
        sys.exit(f'ERROR: {path} lists {int(dup.sum())} model ID(s) twice, '
                 f"e.g. {b.loc[dup, 'gene_id'].head(3).tolist()}")
    asg = b[b['workbook_name'].notna()].copy()
    odd = asg['reason'].ne('assigned')
    if odd.any():
        sys.exit(f'ERROR: {path}: {int(odd.sum())} row(s) carry a workbook name '
                 f'without reason "assigned"; the file has been edited.')
    asg['workbook_name'] = asg['workbook_name'].str.strip().str.lower()
    mapping = dict(zip(asg['gene_id'], asg['workbook_name']))
    info = {'file': fingerprint(path), 'models': int(len(b)),
            'assigned': int(len(asg)),
            'distinct_workbook_genes': int(asg['workbook_name'].nunique()),
            'reasons': {str(k): int(v) for k, v in
                        b['reason'].value_counts().items()}}
    logging.info('Bridge %s: %s models, %s given a workbook gene (%s distinct)',
                 os.path.basename(path), f"{info['models']:,}",
                 f"{info['assigned']:,}",
                 f"{info['distinct_workbook_genes']:,}")
    return mapping, info


def shared_gene_counts(names, regions):
    """(retained, eliminated) counts of analysed genes whose workbook gene is
    also another analysed gene's. `names` and `regions` are aligned Series."""
    shared = names.notna() & names.duplicated(keep=False)
    return (int((shared & regions.eq('retained')).sum()),
            int((shared & regions.eq('eliminated')).sum()))


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def score_join(analysed, name_of, wb, label):
    """region1 x workbook elimination call, and span agreement, for one join."""
    d = analysed.copy()
    d['workbook_name'] = d['gene1'].map(name_of)
    d = d.merge(wb[['name', 'wb_eliminated', 'len']], left_on='workbook_name',
                right_on='name', how='left')
    j = d[d['wb_eliminated'].notna()]
    tab = pd.crosstab(j['region1'], j['wb_eliminated'].astype(bool)).reindex(
        index=['retained', 'eliminated'], columns=[False, True], fill_value=0)
    a, b = int(tab.loc['eliminated', True]), int(tab.loc['eliminated', False])
    c, dd = int(tab.loc['retained', True]), int(tab.loc['retained', False])
    if (a + b) and (c + dd):
        odds, p = fisher_exact([[a, b], [c, dd]])
    else:
        odds, p = np.nan, np.nan
    k = j.dropna(subset=['model_span', 'len'])
    rho = (spearmanr(k['model_span'], k['len']).correlation
           if len(k) > 2 else np.nan)
    return {'join': label, 'n_analysed': len(d), 'n_joined': len(j),
            'elim_elim': a, 'elim_ret': b, 'ret_elim': c, 'ret_ret': dd,
            'pct_eliminated_agree': 100.0 * a / (a + b) if a + b else np.nan,
            'pct_retained_agree': 100.0 * dd / (c + dd) if c + dd else np.nan,
            'odds_ratio': odds, 'fisher_p': p, 'span_spearman_rho': rho}


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Link workbook genes (v37 coordinates) to v50 miniprot '
                    'models by sequence lift and overlap.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotation', required=True,
                   help='Ascaris_genes_15714_annotation_expression.xlsx')
    p.add_argument('--old-genome', required=True,
                   help='v37 germline FASTA the workbook coordinates are on '
                        '(with its .fai)')
    p.add_argument('--new-genome', required=True,
                   help='v50 germline FASTA the miniprot GFF is on')
    p.add_argument('--gff1', required=True, help='v50 miniprot GFF')
    p.add_argument('--kaks', default=None,
                   help='Annotated rate table, for validation on the analysed '
                        'genes (recommended)')
    p.add_argument('--gene-col', default='gene1')
    p.add_argument('--region-col', default='region1')
    p.add_argument('--wb-elim-col', default='eliminated')
    p.add_argument('--wb-elim-pattern', default=r'^eliminated')
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--preset', default='asm5',
                   help='minimap2 preset (default asm5: two assemblies of one '
                        'genome)')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--min-identity', type=float, default=0.95)
    p.add_argument('--min-coverage', type=float, default=0.90,
                   help='Fraction of the span\'s non-N bases the alignment must '
                        'cover (default 0.90)')
    p.add_argument('--min-mapq', type=int, default=5)
    p.add_argument('--max-n-frac', type=float, default=0.5,
                   help='Spans with a larger N fraction are not lifted '
                        '(default 0.5)')
    p.add_argument('--min-model-overlap', type=float, default=0.5,
                   help='Fraction of a model the workbook gene must cover '
                        '(default 0.5)')
    p.add_argument('--ambiguity', type=float, default=0.8,
                   help='A model is ambiguous when the second gene scores at '
                        'least this fraction of the best (default 0.8)')
    p.add_argument('--outdir', default=None,
                   help='Default: results/workbook_model_bridge_<YYYYmmdd_HHMMSS>')
    args = p.parse_args()
    for name in ('min_identity', 'min_coverage', 'max_n_frac',
                 'min_model_overlap', 'ambiguity'):
        v = getattr(args, name)
        if not 0 < v <= 1:
            sys.exit(f"ERROR: --{name.replace('_', '-')} must lie in (0, 1]; got {v}")
    if args.threads < 1:
        sys.exit('ERROR: --threads must be at least 1')
    if args.outdir is None:
        args.outdir = os.path.join('results', 'workbook_model_bridge_' +
                                   datetime.now().strftime('%Y%m%d_%H%M%S'))
    return args


def main():
    # The workbook is read, and the old name join reproduced, exactly as
    # Figure 6 did. Imported here, not at the top: functional_enrichment.py
    # imports read_bridge() from this module.
    from functional_enrichment import load_workbook, target_bridge
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'workbook_model_bridge')
    logging.info('workbook_model_bridge.py')
    logging.info('  outdir: %s', os.path.abspath(args.outdir))
    for f in (args.old_genome, args.new_genome, args.gff1, args.annotation):
        if not os.path.exists(f):
            sys.exit(f'ERROR: missing input {f}')
    mm2_version = minimap2_version(args.minimap2)
    logging.info('  minimap2 %s (%s)', mm2_version, args.minimap2)

    # --- models, and a check that they are on --new-genome (before aligning) ------
    models, align_diag = select_primary_alignment(parse_gff_models(args.gff1))
    models['gene_id'] = models['gene_id'].astype(str)
    gff_seqs = set(models['chrom'])
    v50_seqs = fasta_names(args.new_genome)
    absent = sorted(gff_seqs - v50_seqs)
    logging.info('  GFF sequences: %d, of which %d are in %s', len(gff_seqs),
                 len(gff_seqs) - len(absent), os.path.basename(args.new_genome))
    if len(absent) > 0.5 * len(gff_seqs):
        sys.exit(f'ERROR: {len(absent)} of {len(gff_seqs)} GFF sequences are not '
                 f'in {args.new_genome} (e.g. {absent[:5]}); the GFF is on another '
                 f'assembly or its names differ.')
    if absent:
        logging.warning('  GFF sequences missing from the v50 FASTA: %s', absent[:10])

    # --- rate table for validation: read now, so a bad table fails before aligning --
    analysed = None
    if args.kaks:
        rates = read_rate_table(args.kaks, require=())
        for col in (args.gene_col, args.region_col):
            if col not in rates.columns:
                sys.exit(f'ERROR: column {col!r} not in {args.kaks}')
        analysed = (rates.rename(columns={args.gene_col: 'gene1',
                                          args.region_col: 'region1'})
                    [['gene1', 'region1']].copy())
        analysed['gene1'] = analysed['gene1'].astype(str)
        analysed = analysed[analysed['region1'].isin(['retained', 'eliminated'])]
        n_in = int(analysed['gene1'].isin(set(models['gene_id'])).sum())
        logging.info('  rate table: %s retained/eliminated rows, %s of them '
                     'models in the GFF', f'{len(analysed):,}', f'{n_in:,}')
        if n_in < 0.5 * max(1, len(analysed)):
            sys.exit(f'ERROR: fewer than half the {args.gene_col} identifiers of '
                     f'{args.kaks} are mRNA IDs in {args.gff1}.')

    # --- workbook ---------------------------------------------------------------
    wb = load_workbook(args.annotation)
    need = ['chr', 'start', 'end', 'strand', 'len', args.wb_elim_col]
    missing = [c for c in need if c not in wb.columns]
    if missing:
        sys.exit(f'ERROR: workbook lacks {missing}; columns: {list(wb.columns)}')
    for c in ('start', 'end', 'len'):
        wb[c] = pd.to_numeric(wb[c], errors='coerce')
    wb['chr'] = wb['chr'].astype(str).str.strip()
    wb['strand'] = wb['strand'].astype(str).str.strip()
    wb['wb_eliminated'] = (wb[args.wb_elim_col].astype(str).str.strip()
                           .str.contains(args.wb_elim_pattern, regex=True))
    off = (wb['len'] - (wb['end'] - wb['start'] + 1)).abs()
    logging.info('  workbook len vs end - start + 1: %d of %d genes differ '
                 '(coordinates are taken as 1-based inclusive)',
                 int((off > 0).sum()), len(wb))

    # --- 1. extract v37 spans ------------------------------------------------------
    fai = read_fai(args.old_genome + '.fai')
    fa_path = os.path.join(args.outdir, 'workbook_genes_v37.fa')
    spans = extract_spans(wb, args.old_genome, fai, fa_path)
    n_ext = int(spans['extracted'].sum())
    logging.info('Extracted %s of %s workbook spans from v37 (median %.0f bp)',
                 f'{n_ext:,}', f'{len(spans):,}', spans['span_bp'].median())
    if n_ext == 0:
        sys.exit('ERROR: no span extracted; check that the workbook coordinates '
                 'are on --old-genome.')

    # --- 2. align to v50 ----------------------------------------------------------
    paf_path = os.path.join(args.outdir, 'workbook_genes_v50.paf')
    mm2_cmd = run_minimap2(args.minimap2, args.preset, args.threads,
                           args.new_genome, fa_path, paf_path,
                           os.path.join(args.outdir, 'minimap2.log'))
    paf = read_paf(paf_path)
    logging.info('  %s primary PAF records for %s genes', f'{len(paf):,}',
                 f'{paf["qname"].nunique():,}')
    lift = lift_genes(spans, paf, args.min_identity, args.min_coverage,
                      args.min_mapq, args.max_n_frac)
    lift['wb_eliminated'] = lift['workbook_name'].map(
        dict(zip(wb['name'], wb['wb_eliminated'])))
    lift = lift.drop(columns=['reason'])
    lift.to_csv(os.path.join(args.outdir, 'workbook_lift.tsv'), sep='\t',
                index=False, float_format='%.6g')
    lifted = lift[lift['lifted']].copy()
    lifted[['v50_start', 'v50_end']] = lifted[['v50_start', 'v50_end']].astype(np.int64)
    with open(os.path.join(args.outdir, 'workbook_genes_v50.bed'), 'w') as fh:
        for r in lifted.sort_values(['v50_chr', 'v50_start']).itertuples():
            fh.write(f'{r.v50_chr}\t{r.v50_start - 1}\t{r.v50_end}\t'
                     f'{r.workbook_name}\t{int(round(1000 * r.identity))}\t'
                     f'{r.v50_strand}\n')
    lift_tab = (pd.crosstab(lift['reason_lift'].replace('', 'lifted'),
                            lift['wb_eliminated'].map({True: 'eliminated',
                                                       False: 'retained'}))
                .reindex(columns=['retained', 'eliminated'], fill_value=0))
    lift_tab = lift_tab.loc[lift_tab.sum(axis=1).sort_values(ascending=False).index]
    logging.info('Lifted %s of %s workbook genes to v50', f'{len(lifted):,}',
                 f'{len(lift):,}')
    no_model_seq = sorted(set(lifted['v50_chr']) - gff_seqs)
    if no_model_seq:
        logging.info('  %d v50 sequences carry lifted genes but no model',
                     len(no_model_seq))

    # --- 3. models <- lifted genes ---------------------------------------------------
    bridge = assign_models(models, lifted, args.min_model_overlap, args.ambiguity)
    keep = ['workbook_name', 'v37_chr', 'v37_start', 'v37_end', 'workbook_strand',
            'v50_chr', 'v50_start', 'v50_end', 'v50_strand', 'identity',
            'coverage', 'mapq', 'lift_mode', 'wb_eliminated']
    bridge = bridge.merge(lifted[keep], on='workbook_name', how='left')
    bridge.to_csv(os.path.join(args.outdir, 'workbook_model_bridge.tsv'),
                  sep='\t', index=False, float_format='%.6g')
    b_reasons = bridge['reason'].value_counts()
    n_assigned = int((bridge['reason'] == 'assigned').sum())
    n_shared = int((bridge['models_sharing_gene'] > 1).sum())
    n_genes_used = int(bridge['workbook_name'].nunique())
    logging.info('Assigned %s of %s v50 models a workbook gene (%s distinct genes)',
                 f'{n_assigned:,}', f'{len(bridge):,}', f'{n_genes_used:,}')
    n_opp = int(b_reasons.get('only an opposite-strand workbook gene overlaps it', 0))
    if n_opp > 0.2 * max(1, n_assigned):
        logging.warning('  %s models overlap only opposite-strand genes, against '
                        '%s assigned: check the strand handling before using '
                        'the bridge', f'{n_opp:,}', f'{n_assigned:,}')

    # --- 4. validation on the analysed genes ------------------------------------------
    val = pd.DataFrame()
    status = pd.DataFrame()
    same_gene = np.nan
    if analysed is not None:
        span = dict(zip(models['gene_id'], models['end'] - models['start'] + 1))
        analysed['model_span'] = analysed['gene1'].map(span)
        bridge_map = dict(zip(bridge['gene_id'], bridge['workbook_name']))
        name_of, _, _ = target_bridge(args.gff1)
        val = pd.DataFrame([score_join(analysed, bridge_map, wb, 'position bridge'),
                            score_join(analysed, name_of, wb,
                                       'name join (Figure 6)')])
        by_bridge = analysed['gene1'].map(bridge_map)
        by_name = analysed['gene1'].map(name_of)
        ok = by_bridge.notna() & by_name.notna()
        same_gene = (100.0 * float((by_bridge[ok] == by_name[ok]).mean())
                     if ok.any() else np.nan)
        val['pct_same_gene_both_joins'] = same_gene
        val.to_csv(os.path.join(args.outdir, 'bridge_validation.tsv'), sep='\t',
                   index=False, float_format='%.6g')
        reason_of = dict(zip(bridge['gene_id'], bridge['reason']))
        status = (pd.crosstab(analysed['gene1'].map(reason_of)
                              .fillna('model not in the GFF'),
                              analysed['region1'])
                  .reindex(columns=['retained', 'eliminated'], fill_value=0))
        status = status.loc[status.sum(axis=1).sort_values(ascending=False).index]
        status.to_csv(os.path.join(args.outdir, 'analysed_bridge_status.tsv'),
                      sep='\t')

    # --- summary -----------------------------------------------------------------------
    L = ['Workbook genes (v37) -> v50 miniprot models, by sequence lift and overlap',
         '',
         f'Workbook   {os.path.basename(args.annotation)}, {len(wb):,} genes',
         f'v37        {args.old_genome}',
         f'v50        {args.new_genome}',
         f'Models     {os.path.basename(args.gff1)}, {len(models):,} '
         f'on {len(gff_seqs)} sequences',
         f'minimap2   {mm2_version}: -x {args.preset} -c --secondary=no',
         '',
         f'1. Extraction: {n_ext:,} of {len(spans):,} spans (median '
         f'{spans["span_bp"].median():,.0f} bp; '
         f'{int((spans["n_bases"] > 0).sum()):,} contain N)',
         '',
         f'2. Lift (coverage >= {args.min_coverage:g} of non-N bases, identity >= '
         f'{args.min_identity:g}, MAPQ >= {args.min_mapq}): {len(lifted):,} of '
         f'{len(lift):,} genes',
         f'   {"":<45} {"retained":>9} {"eliminated":>11}']
    for k, r in lift_tab.iterrows():
        L.append(f'   {k:<45} {r["retained"]:>9,} {r["eliminated"]:>11,}')
    if len(lifted):
        q = lifted['identity'].quantile([0.01, 0.5]).to_numpy()
        lr = lifted['len_ratio'].quantile([0.01, 0.5, 0.99]).to_numpy()
        L += [f'   lifted: identity median {q[1]:.4f} (1st percentile {q[0]:.4f}); '
              f'{int((lifted["lift_mode"] == "chained").sum()):,} from a split '
              f'alignment rejoined',
              f'   v50 length / v37 span: median {lr[1]:.3f} (1st-99th percentile '
              f'{lr[0]:.3f}-{lr[2]:.3f})']
    L += ['',
          f'3. Models (same strand, overlap >= {args.min_model_overlap:g} of the '
          f'model, ambiguity {args.ambiguity:g}): {n_assigned:,} of '
          f'{len(bridge):,} assigned, {n_genes_used:,} distinct workbook genes']
    for k, n in b_reasons.items():
        L.append(f'     {k:<52} {n:>7,}')
    L.append(f'   models sharing their workbook gene with another model: {n_shared:,}')
    asg = bridge[bridge['reason'] == 'assigned']
    if len(asg):
        L.append(f'   assigned: median overlap {asg["frac_model"].median():.2f} of the '
                 f'model, {asg["frac_gene"].median():.2f} of the workbook gene')
    if len(val):
        L += ['', '4. Validation on the analysed genes (no names used): workbook '
              'elimination call vs region1, and model span vs workbook gene span',
              f'   {"join":<22} {"joined":>7} {"elim/elim":>10} {"elim/ret":>9} '
              f'{"ret/elim":>9} {"ret/ret":>8} {"elim agree":>11} '
              f'{"ret agree":>10} {"OR":>9} {"span rho":>9}']
        for r in val.itertuples():
            L.append(f'   {r.join:<22} {r.n_joined:>7,} {r.elim_elim:>10,} '
                     f'{r.elim_ret:>9,} {r.ret_elim:>9,} {r.ret_ret:>8,} '
                     f'{r.pct_eliminated_agree:>10.1f}% '
                     f'{r.pct_retained_agree:>9.1f}% {r.odds_ratio:>9.3g} '
                     f'{r.span_spearman_rho:>+9.2f}')
        L += [f'   Analysed genes given the same workbook gene by both joins: '
              f'{same_gene:.1f}%',
              '   A correct join agrees on both calls and has a strongly positive '
              'span rho.',
              '', '   Analysed genes by bridge outcome:',
              f'   {"":<52} {"retained":>9} {"eliminated":>11}']
        for k, r in status.iterrows():
            L.append(f'   {k:<52} {r["retained"]:>9,} {r["eliminated"]:>11,}')
    L += ['', 'Use: pass --bridge workbook_model_bridge.tsv to the expression '
          'scripts, which then join through gene_id -> workbook_name.', '']
    summary = '\n'.join(L)
    with open(os.path.join(args.outdir, 'workbook_model_bridge_summary.txt'),
              'w') as fh:
        fh.write(summary + '\n')

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({
            'script': os.path.basename(__file__),
            'git_commit': git_commit(_HERE),
            'run_utc': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            'host': platform.node(), 'python': platform.python_version(),
            'packages': {'numpy': np.__version__, 'pandas': pd.__version__,
                         'scipy': scipy_version},
            'minimap2': {'version': mm2_version, 'command': mm2_cmd},
            'chaining': {'max_stretch': CHAIN_MAX_STRETCH, 'pad_bp': CHAIN_PAD_BP},
            'argv': sys.argv, 'args': vars(args),
            'inputs': {k: fingerprint(v) for k, v in
                       (('annotation', args.annotation), ('gff1', args.gff1),
                        ('kaks', args.kaks),
                        ('old_genome_fai', args.old_genome + '.fai'),
                        ('new_genome_fai', args.new_genome + '.fai'))},
            'counts': {'workbook_genes': len(wb), 'extracted': n_ext,
                       'lifted': len(lifted), 'models': len(bridge),
                       'models_assigned': n_assigned,
                       'distinct_workbook_genes_assigned': n_genes_used,
                       'models_sharing_gene': n_shared},
            'lift_reasons': {k: v.to_dict() for k, v in lift_tab.iterrows()},
            'model_reasons': b_reasons.to_dict(),
            'gff_sequences_missing_from_v50': absent,
            'alignment_selection': align_diag,
            'validation': val.to_dict(orient='records'),
            'analysed_status': {k: v.to_dict() for k, v in status.iterrows()},
        }, fh, indent=2, default=str)

    print()
    print(summary)
    logging.info('Done.')


if __name__ == '__main__':
    main()
