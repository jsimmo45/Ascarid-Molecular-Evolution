#!/usr/bin/env python3
"""
toxocara_annotation_diagnostic.py -- Is the Toxocara dN elevation an artefact
of TransDecoder gene models? Supplementary figure S9 (S13 before the
2026-10-12 renumbering; the file stem figS13_ is kept, and the assembly copies
it as FigS9) and Table T7.

BACKGROUND
----------
The four-way branch analysis (branch_dnds_four.py, unrooted tree, 2026-10-10)
gives eliminated/retained dN fold changes of Ascaris 1.96, Parascaris 2.07,
Baylisascaris 3.17 and Toxocara 2.00 (as_region). Toxocara CDS are TransDecoder
ORFs on 2017 StringTie transcript models (82% complete), so the standing worry
is that its eliminated genes -- reproductive-tissue biased -- carry poorer gene
models than its retained genes, and that frameshift and truncation noise
inflates their rates. Since the rerun the same branch also shows a small dS
excess in eliminated genes (1.10-fold, q 3.0e-5), which gene-model noise would
produce too, because it inflates dN and dS together.

This script tests that WITHOUT re-annotating anything, using files already on
disk. Three tests, in increasing order of how much they assume:

  TEST 1  PREMISE.  Do Toxocara eliminated genes actually have worse gene
          models than retained ones? Compares ORF completeness rate, ORF
          length, alignment coverage and percent identity between the two
          classes. If these are equal, class-specific annotation noise does not
          exist and tests 2-3 are moot.

  TEST 2  RESTRICTION.  Recompute the fold changes in progressively cleaner
          subsets (complete ORFs only; high-coverage alignments only; all
          filters). Restricting rather than splitting keeps power: the
          eliminated class is 112 genes (as_region) or 165 (tox_region), so
          splitting it into complete-vs-partial arms would leave a handful in
          one arm. dN and dN/dS use the capped rows and dS every gene, exactly
          as Figure 5 does.

          dS IS THE NEW PART (2026-10-10). If the Toxocara dS excess comes from
          poor gene models, it should shrink in complete ORFs while the dN
          excess holds. That decides how Results words the dS result.

  TEST 3  REGRESSION.  Fit log(dN) ~ eliminated (+ gene-model covariates) and
          ask whether the eliminated coefficient survives. CONTROL: the same
          Toxocara-derived covariates are applied to all four species. If they
          attenuate the Toxocara coefficient specifically, that is annotation
          noise; if they attenuate all four roughly equally, the covariates are
          proxying gene class and the attenuation means nothing. Do not report
          test 3 without the control.

INTERPRETATION
--------------
  Equal ORF quality by class (test 1)          -> concern is dead, stop here
  Fold holds in clean subset (test 2)          -> concern is dead
  Coefficient survives covariates (test 3)     -> concern is dead
  Toxocara coefficient attenuates and the      -> real annotation problem;
    other three do not (test 3 control)           re-annotate before publishing


CONVERSION, 2026-10-10
----------------------
Moved into the repo (05_analyses/) from the root of the miniprot directory,
where no repo script produced S13 or T7. What changes:

  * Input. The 2026-10-10 annotated table, read through branch_lib's
    load_annotated(), which refuses a table made with the rooted tree (every
    run before that date; its Toxocara values were half a branch). The v2
    table the last run used was rooted AND came from the pre-fix quartet set.
  * Masking. Genes usable on a branch come from branch_lib.usable_mask(), the
    definition Figure 5 uses: dN/dS <= the cap and dS > 0 for dN and dN/dS, no
    cap for dS. The local copy lacked dS > 0, so a gene with no substitution
    on a branch passed or failed the cap on codeml's random start.
  * TEST 2 reports dS and dN/dS beside dN, with BH q within each labelling and
    metric (across branches x strata). Stars come from q; a stratum with fewer
    than --min-stratum-n eliminated genes is hatched and carries no star.
  * Figure in the house style (paper_colors): one figure per labelling, dN on
    top and dS below, built at 185 mm with 8/7 pt type, bars in ink because a
    fold change is a metric, the rule at 1 behind the bars, n in the tick
    labels, no suptitle or footer (legend notes go to the summary). Saved as
    figS13_fold_by_stratum_<grouping>.{svg,pdf,png}.
  * The summary ends with the numbers Results and Methods quote (ORF
    completeness by class, the complete-ORF dN range, dS in complete ORFs, the
    test 3 control), computed rather than left to be read off tables.
  * provenance.json; logging via region_lib.

FIGURE FIXES, 2026-10-12 (figure only; every table is unchanged)
-----------------------------------------------------------------
  * Tick labels "Complete ORF" and "n = 112" (were "ORF" and "n=112") when they
    fit, which they do with the two strata the TransDecoder evidence supports;
    with more strata the short names come back (tick_labels() measures, and
    the summary says which form was drawn).
  * Each row's y-axis is shared by the figures of both labellings (the dN axes
    were 0-3.5 and 0-4), since the two are read species by species.
  * The summary's legend notes say what each file drew: tick-label form,
    y-limits, and that nothing on the figure names its labelling.

INPUTS
------
The annotated table is required. Every evidence source is optional, and each
contributes different covariates:

  --transdecoder  TransDecoder .pep or .cds FASTA for the focal species.
                  Headers carry "type:complete", "type:5prime_partial",
                  "type:3prime_partial" or "type:internal", and "len:".
                  Contributes: orf_type, orf_len
  --aln-dir       Directory of per-gene codeml input alignments (FASTA or
                  PHYLIP). Contributes: aln_len, focal_cov, focal_ident.
                  branch_dnds_four.py does not keep these; align_len from the
                  annotated table is used instead (--table-covariates auto).
  --rbh           Reciprocal-best-hit TSV(s). Column names are auto-detected
                  from common BLAST outfmt-6 headers.
                  Contributes: rbh_pident, rbh_cov, rbh_bitscore

ID MATCHING IS THE MAIN FAILURE MODE. Run --inspect FIRST: it reports how many
IDs matched from each source and exits. The Baylisascaris BLAST identifier bug
(g10000 uppercased to G10000, silently dropping every gene above g9999) is why.

USAGE
-----
  python toxocara_annotation_diagnostic.py \\
    --annotated results/branch_dnds_four_20261010_113902/branch_dnds_four_annotated.tsv \\
    --transdecoder TG.transcripts.longest.20065.fasta.transdecoder.pep \\
    --rbh rbh_As_Tox_v2/rbh.tsv rbh_Pu_Tox_v2/rbh.tsv rbh_Tox_Bsc_v2/rbh.tsv \\
    --outdir results/toxocara_annotation_diagnostic_<stamp> --inspect
  # then the same without --inspect

OUTPUTS
-------
  annotation_covariates.tsv          per-gene covariates (for reuse)
  test1_premise.tsv                  T7: gene-model quality by class
  test2_strata.tsv                   T7: fold changes, species x labelling x
                                     metric x stratum, with q
  test3_regression.tsv               T7: OLS coefficients, all species, +/- covariates
  figS13_fold_by_stratum_<grouping>.{svg,pdf,png}
  diagnostic_summary.txt             verdicts, numbers for the text, legend notes
  diagnostic.log, provenance.json

Statistics: Mann-Whitney U for class comparisons, Fisher's exact test for ORF
completeness, BH-FDR q for TEST 2; OLS on log(dN + pseudocount) with HC3
heteroscedasticity-robust standard errors (dN is skewed and its variance grows
with the mean, so classical SEs would be too small).

Dependencies: pandas, numpy, scipy, matplotlib; paper_colors, region_lib and
branch_lib from lib/.
"""

import argparse
import glob
import gzip
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
from scipy.stats import mannwhitneyu, fisher_exact
from scipy.stats import t as t_dist

from paper_colors import (NEUTRAL, COL, FONT, apply_style, audit, save,
                          make_outdir, mm, SPECIES_LABEL)
from region_lib import (benjamini_hochberg, stars, fingerprint, git_commit,
                        setup_logging)
from branch_lib import (SPECIES_ORDER, SPECIES_KEY, usable_mask, compare,
                        load_annotated, fit_bottom_margin)

_HERE = os.path.dirname(os.path.abspath(__file__))

ORF_TYPES = ['complete', '5prime_partial', '3prime_partial', 'internal']

# The pipeline names per-species ID columns by abbreviation (as_gene, pu_gene,
# bsc_gene, tox_gene) rather than by full genus, so both forms are searched.
SPECIES_ALIASES = {
    'Ascaris': ['ascaris', 'as', 'asu'],
    'Parascaris': ['parascaris', 'pu', 'pun', 'puni'],
    'Baylisascaris': ['baylisascaris', 'bsc', 'bs'],
    'Toxocara': ['toxocara', 'tox', 'tc'],
}

# Columns that may already exist in the annotated table and are usable as
# alignment-quality covariates without any external file.
TABLE_COVARIATE_CANDIDATES = ['align_len', 'aln_len', 'alignment_len',
                              'alignment_length', 'n_sites', 'ncodons']

# Covariates that are DOWNSTREAM OF SEQUENCE DIVERGENCE rather than measures of
# gene-model quality. Eliminated genes really are faster-evolving at the protein
# level - that is the finding under test - so their BLAST identity and bitscore
# are lower for biological reasons, not annotation reasons. Putting these into
# the regression would be regressing dN on a proxy for itself: the eliminated
# coefficient would collapse and be misread as "gene models explain it".
# They stay in TEST 1 as description, flagged, but are kept out of TEST 3 unless
# --include-identity is passed.
DIVERGENCE_CONFOUNDED = ('rbh_pident', 'rbh_bitscore', 'focal_ident')

# Structural covariates - properties of the gene model itself, not of how
# similar it is to another species. Only these support a clean causal reading.
# ORF completeness is the one true member; the others are partial (a diverged
# gene also aligns over less of its length), which is why losing the
# TransDecoder output degrades the regression rather than merely shrinking it.
STRUCTURAL_CLEAN = ('orf_type', 'orf_len')

# Never treated as an ID column when falling back to a substring search
ID_EXCLUDE_SUFFIXES = ('_dn', '_ds', '_omega', '_region', '_t', '_n', '_s',
                       '_len', '_length', '_status')

# Candidate RBH column names, in preference order, mapped to our covariate name
RBH_CANDIDATES = {
    'rbh_pident': ['pident', 'percent_identity', 'identity', 'pct_identity',
                   'perc_ident'],
    'rbh_cov': ['qcovs', 'qcovhsp', 'query_coverage', 'coverage', 'qcov',
                'scov', 'scovs'],
    'rbh_bitscore': ['bitscore', 'bit_score', 'score'],
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def smart_open(path):
    """Open plain or gzipped text."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def normalise_id(raw, strip_suffix=True):
    """
    Reduce an identifier to the form most likely to match across files.

    TransDecoder appends '.p1', '.p2' ... to the parent transcript ID, and
    aligners often carry a '|' or whitespace-delimited description. This strips
    the description, then optionally the TransDecoder ORF suffix.

    NOTE: case is preserved deliberately. BLAST's accession normalisation
    uppercases IDs matching its accession pattern (this is what silently
    dropped every Baylisascaris gene above g9999), so a case-insensitive match
    here would hide exactly the class of bug we have already been bitten by.
    Use --lower-ids only if --inspect shows a case mismatch you have confirmed
    is benign.
    """
    if raw is None:
        return None
    s = str(raw).strip().split()[0]
    s = s.split('|')[0]
    if strip_suffix:
        s = re.sub(r'\.p\d+$', '', s)
    return s


def flag_true(series):
    """
    Plain boolean Series from a 1.0/0.0/NaN flag column, NaN treated as False.

    Flags are stored as float rather than bool because a bool column cannot
    hold the "no evidence for this gene" state, and pandas raises rather than
    silently upcasting when NaN is assigned into one.
    """
    return (pd.to_numeric(series, errors='coerce') == 1).fillna(False).astype(bool)


FLAG_COLS = ('orf_complete', 'high_cov', 'long_aln')


def continuous_covariates(cov):
    """
    Numeric covariate columns available for testing.

    Derived boolean flags and the categorical ORF type are excluded; everything
    else that has at least one finite value counts. Detecting these dynamically
    means covariates read straight out of the annotated table (align_len and
    friends) are picked up without touching the test functions.
    """
    out = []
    for c in cov.columns:
        if c == 'orf_type' or c in FLAG_COLS:
            continue
        if pd.to_numeric(cov[c], errors='coerce').notna().any():
            out.append(c)
    return out


def as_flag(condition, missing):
    """
    Encode a boolean condition as a 1.0/0.0/NaN float flag.

    Parameters
    ----------
    condition : boolean Series
    missing : boolean Series
        True where the underlying evidence was absent, so the flag is unknown.
    """
    out = condition.astype(float)
    out[missing.to_numpy()] = np.nan
    return out


def detect_column(df, candidates, contains=None):
    """
    Find the first column matching any candidate name (case-insensitive), or
    failing that the first column whose name contains all of `contains`.

    Returns None if nothing matches.
    """
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if contains:
        for c in df.columns:
            if all(tok.lower() in c.lower() for tok in contains):
                return c
    return None


# ---------------------------------------------------------------------------
# Evidence source parsers
# ---------------------------------------------------------------------------

def parse_transdecoder(path, strip_suffix=True):
    """
    Extract ORF completeness class and ORF length from TransDecoder headers.

    A TransDecoder header looks like:
      >TG.1234.p1 GENE.TG.1234~~TG.1234.p1  ORF type:complete len:412 (+),...

    Returns
    -------
    DataFrame with columns [focal_id, orf_type, orf_len]
    """
    rows = []
    n_headers = 0
    with smart_open(path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            n_headers += 1
            head = line[1:].rstrip('\n')
            fid = normalise_id(head, strip_suffix=strip_suffix)

            m = re.search(r'type:(\w+)', head)
            otype = m.group(1) if m else None
            if otype is not None and otype not in ORF_TYPES:
                # TransDecoder occasionally writes other tokens; keep them but
                # they will not be used as the 'complete' reference level
                pass

            m = re.search(r'len:(\d+)', head)
            olen = int(m.group(1)) if m else np.nan

            rows.append({'focal_id': fid, 'orf_type': otype, 'orf_len': olen})

    df = pd.DataFrame(rows)
    if len(df):
        # A transcript can yield several ORFs; keep the longest, which is what
        # a longest-ORF CDS set would have used.
        df = (df.sort_values('orf_len', ascending=False)
                .drop_duplicates('focal_id', keep='first')
                .reset_index(drop=True))
    logging.info(f"TransDecoder: {n_headers} headers -> {len(df)} unique IDs; "
                 f"type parsed for {int(df['orf_type'].notna().sum()) if len(df) else 0}")
    if len(df) and df['orf_type'].notna().any():
        counts = df['orf_type'].value_counts().to_dict()
        logging.info(f"  ORF types: {counts}")
    return df


def read_alignment(path):
    """
    Read a FASTA or PHYLIP alignment into {name: sequence}.

    Detects PHYLIP by a first line of exactly two integers (the codeml
    convention). Supports sequential PHYLIP, which is what pal2nal writes.

    Returns
    -------
    dict, or {} if the file could not be parsed
    """
    try:
        with smart_open(path) as fh:
            lines = [ln.rstrip('\n') for ln in fh]
    except OSError:
        return {}

    lines = [ln for ln in lines if ln.strip() != '']
    if not lines:
        return {}

    first = lines[0].split()
    is_phylip = len(first) == 2 and all(tok.isdigit() for tok in first)

    seqs = {}
    if is_phylip:
        ntax, nchar = int(first[0]), int(first[1])
        name, buf = None, []
        for ln in lines[1:]:
            parts = ln.split(None, 1)
            starts_new = (name is None
                          or (len(''.join(buf)) >= nchar and len(seqs) < ntax))
            if starts_new and parts:
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = parts[0]
                buf = [parts[1].replace(' ', '')] if len(parts) > 1 else []
            else:
                buf.append(ln.replace(' ', ''))
        if name is not None:
            seqs[name] = ''.join(buf)
    else:
        name, buf = None, []
        for ln in lines:
            if ln.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = ln[1:].split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            seqs[name] = ''.join(buf)

    return seqs


def alignment_metrics(seqs, focal_pattern):
    """
    Per-alignment quality metrics for the focal sequence.

    focal_cov   fraction of alignment columns where the focal sequence is not a
                gap. A truncated or frameshifted gene model shows up here.
    focal_ident mean percent identity of the focal sequence to each other
                sequence, over columns where both are ungapped.
    aln_len     alignment length in columns.

    Returns
    -------
    dict or None if the focal sequence is absent
    """
    if not seqs:
        return None

    focal_name = None
    pat = focal_pattern.lower()
    for name in seqs:
        if pat in name.lower():
            focal_name = name
            break
    if focal_name is None:
        return None

    f = seqs[focal_name].upper()
    aln_len = len(f)
    if aln_len == 0:
        return None

    gapchars = set('-?.')
    f_arr = np.frombuffer(f.encode(), dtype='S1')
    f_gap = np.isin(f_arr, [c.encode() for c in gapchars])
    focal_cov = float((~f_gap).sum()) / aln_len

    idents = []
    for name, s in seqs.items():
        if name == focal_name:
            continue
        s = s.upper()
        if len(s) != aln_len:
            continue
        s_arr = np.frombuffer(s.encode(), dtype='S1')
        s_gap = np.isin(s_arr, [c.encode() for c in gapchars])
        both = (~f_gap) & (~s_gap)
        if both.sum() == 0:
            continue
        idents.append(float((f_arr[both] == s_arr[both]).sum()) / both.sum())

    return {'aln_len': aln_len, 'focal_cov': focal_cov,
            'focal_ident': float(np.mean(idents)) if idents else np.nan}


def parse_alignments(aln_dir, patterns, focal_pattern, strip_suffix=True,
                     id_strip=None):
    """
    Walk an alignment directory and compute per-gene metrics.

    The gene ID is taken from the filename stem, with known alignment suffixes
    removed. --aln-id-strip appends further suffixes to remove.

    Returns
    -------
    DataFrame with columns [focal_id, aln_len, focal_cov, focal_ident]
    """
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(aln_dir, '**', pat),
                               recursive=True))
    files = sorted(set(files))
    logging.info(f"Alignments: {len(files)} files matched in {aln_dir}")

    strip_tokens = ['.aln', '.fa', '.fasta', '.phy', '.phylip', '.nuc',
                    '.codon', '.pal2nal', '.mafft', '.gz', '.cds', '.best']
    if id_strip:
        strip_tokens.extend(id_strip)

    rows, n_ok, n_nofocal, n_bad = [], 0, 0, 0
    for fp in files:
        stem = os.path.basename(fp)
        changed = True
        while changed:
            changed = False
            for tok in strip_tokens:
                if stem.lower().endswith(tok):
                    stem = stem[: -len(tok)]
                    changed = True
        gid = normalise_id(stem, strip_suffix=strip_suffix)

        seqs = read_alignment(fp)
        if not seqs:
            n_bad += 1
            continue
        m = alignment_metrics(seqs, focal_pattern)
        if m is None:
            n_nofocal += 1
            continue
        m['focal_id'] = gid
        rows.append(m)
        n_ok += 1

    logging.info(f"  parsed OK: {n_ok}; focal sequence not found: {n_nofocal}; "
                 f"unreadable: {n_bad}")
    if n_nofocal and n_ok == 0:
        logging.warning(f"  no alignment contained a sequence matching "
                        f"'{focal_pattern}' - set --aln-focal-pattern to the "
                        f"name used in your alignment files")

    df = pd.DataFrame(rows)
    if len(df):
        df = df.drop_duplicates('focal_id', keep='first').reset_index(drop=True)
    return df


def parse_rbh(paths, strip_suffix=True):
    """
    Pull identity/coverage/bitscore out of reciprocal-best-hit tables.

    Column names vary between BLAST output formats, so they are auto-detected.
    The focal ID column is taken to be the first column unless a column named
    like a query ID is found.

    Returns
    -------
    DataFrame with columns [focal_id, rbh_pident, rbh_cov, rbh_bitscore]
    """
    frames = []
    for path in paths:
        try:
            df = pd.read_csv(path, sep='\t')
        except Exception as exc:
            logging.warning(f"RBH: could not read {path}: {exc}")
            continue

        id_col = detect_column(df, ['qseqid', 'query', 'query_id', 'qid',
                                    'gene', 'gene_id'])
        if id_col is None:
            id_col = df.columns[0]
            logging.warning(f"RBH {os.path.basename(path)}: no obvious ID "
                            f"column, using first column '{id_col}'")

        out = pd.DataFrame(
            {'focal_id': df[id_col].map(
                lambda v: normalise_id(v, strip_suffix=strip_suffix))})
        found = []
        for our, cands in RBH_CANDIDATES.items():
            col = detect_column(df, cands)
            if col is not None:
                out[our] = pd.to_numeric(df[col], errors='coerce')
                found.append(f'{our}<-{col}')
        logging.info(f"RBH {os.path.basename(path)}: {len(df)} rows, "
                     f"id='{id_col}', mapped [{', '.join(found) or 'nothing'}]")
        if len(out.columns) > 1:
            frames.append(out)

    if not frames:
        return pd.DataFrame(columns=['focal_id'])

    allf = pd.concat(frames, ignore_index=True)
    num = [c for c in allf.columns if c != 'focal_id']
    # A gene appears in several pairwise RBH files; take the median so one bad
    # pairing does not dominate.
    return allf.groupby('focal_id', as_index=False)[num].median()


# ---------------------------------------------------------------------------
# Assembling the covariate table
# ---------------------------------------------------------------------------

def find_focal_id_col(df, focal, override=None):
    """Locate the column holding focal-species gene IDs in the annotated table."""
    if override:
        if override not in df.columns:
            sys.exit(f"ERROR: --focal-id-col '{override}' not in the annotated "
                     f"table. Available: {list(df.columns)}")
        return override

    lower = {c.lower(): c for c in df.columns}
    suffixes = ('_gene', '_id', '_gene_id', '_transcript', '_protein', '_seq',
                '_prot', '')

    # Full genus first, then the pipeline's abbreviations (tox_gene etc.)
    for name in [focal.lower()] + SPECIES_ALIASES.get(focal, []):
        for suffix in suffixes:
            cand = f'{name}{suffix}'
            if cand in lower:
                return lower[cand]

    # Last resort: any column mentioning the species that is not a metric
    for name in [focal.lower()] + SPECIES_ALIASES.get(focal, []):
        for c in df.columns:
            cl = c.lower()
            if name in cl and not cl.endswith(ID_EXCLUDE_SUFFIXES):
                return c
    return None


def make_key(df, col, strip_suffix=True, lower=False):
    """Normalised join key from one column of the annotated table."""
    key = df[col].map(lambda v: normalise_id(v, strip_suffix=strip_suffix))
    return key.str.lower() if lower else key


def choose_join_column(df, src_ids, preferred=None, strip_suffix=True,
                       lower=False):
    """
    Pick the annotated-table column whose IDs best match a source's IDs.

    Different evidence sources key on different things: TransDecoder and RBH
    files carry focal-species gene IDs, but per-gene alignment files are
    normally named by quartet/ortholog ID instead. Rather than make the user
    work that out, every candidate column is scored by match rate and the best
    one wins. `preferred` is tried first and kept on a tie.

    Returns
    -------
    (column name or None, match rate in [0, 1])
    """
    src = set(src_ids.dropna().astype(str))
    if lower:
        src = {s.lower() for s in src}
    if not src:
        return None, 0.0

    cands = []
    if preferred and preferred in df.columns:
        cands.append(preferred)
    cands += [c for c in df.columns if c not in cands]

    best, best_rate = None, -1.0
    for c in cands:
        try:
            key = make_key(df, c, strip_suffix, lower)
        except Exception:
            continue
        rate = float(key.isin(src).mean())
        if rate > best_rate:
            best, best_rate = c, rate
    return best, max(best_rate, 0.0)


def build_covariates(df, focal_id_col, sources, min_coverage,
                     strip_suffix=True, lower_ids=False, overrides=None):
    """
    Attach every available evidence source to the annotated table's genes.

    Each source is joined on whichever annotated-table column matches its IDs
    best (see choose_join_column), so alignments keyed by quartet ID and
    FASTA/RBH files keyed by focal gene ID both work without configuration.

    Parameters
    ----------
    sources : dict of name -> DataFrame with a 'focal_id' column
    overrides : dict of name -> forced join column

    Returns
    -------
    (DataFrame of covariates, dict of match/join statistics)
    """
    cov = pd.DataFrame(index=df.index)
    stats = {'n_genes': len(df), 'joins': {}}
    overrides = overrides or {}

    for name, src in sources.items():
        if src is None or not len(src):
            stats[f'{name}_matched'] = 0
            continue

        forced = overrides.get(name)
        if forced and forced not in df.columns:
            sys.exit(f"ERROR: --{name}-id-col '{forced}' is not a column of "
                     f"the annotated table. Available: {list(df.columns)}")
        col, rate = choose_join_column(
            df, src['focal_id'], preferred=forced or focal_id_col,
            strip_suffix=strip_suffix, lower=lower_ids)
        if forced:
            col = forced
        stats['joins'][name] = (col, rate)

        key = make_key(df, col, strip_suffix, lower_ids)
        s = src.drop_duplicates('focal_id', keep='first').copy()
        if lower_ids:
            s['focal_id'] = s['focal_id'].astype(str).str.lower()
        s = s.set_index('focal_id')

        merged = s.reindex(key.to_numpy())
        merged.index = df.index
        newcols = list(s.columns)
        for c in newcols:
            cov[c] = merged[c]

        matched = int(merged[newcols[0]].notna().sum()) if newcols else 0
        stats[f'{name}_matched'] = matched
        pct = 100 * matched / max(len(df), 1)
        logging.info(f"Merged {name}: joined on '{col}', "
                     f"{matched}/{len(df)} genes matched ({pct:.1f}%)")
        if matched == 0:
            logging.warning(f"  ZERO matches from {name} - this is an ID "
                            f"mismatch, not a data problem. Run --inspect.")
        elif pct < 50:
            logging.warning(f"  low match rate from {name} ({pct:.1f}%) - "
                            f"check --inspect before trusting the result")

    # Derived flags, stored as 1.0 / 0.0 / NaN (see flag_true)
    if 'orf_type' in cov.columns:
        cov['orf_complete'] = as_flag(cov['orf_type'].eq('complete'),
                                      cov['orf_type'].isna())
    if 'focal_cov' in cov.columns:
        v = pd.to_numeric(cov['focal_cov'], errors='coerce')
        cov['high_cov'] = as_flag(v >= min_coverage, v.isna())
    elif 'rbh_cov' in cov.columns:
        v = pd.to_numeric(cov['rbh_cov'], errors='coerce')
        # BLAST reports coverage as a percentage; ours is a fraction
        scale = 100.0 if v.max(skipna=True) > 1.5 else 1.0
        cov['high_cov'] = as_flag((v / scale) >= min_coverage, v.isna())

    return cov, stats


# ---------------------------------------------------------------------------
# TEST 1 - premise
# ---------------------------------------------------------------------------

def test1_premise(df, cov, groupings):
    """
    Do eliminated genes have worse focal-species gene models than retained?

    Continuous covariates: Mann-Whitney. Binary covariates (ORF completeness):
    Fisher exact on the 2x2 of class by complete/not.

    A uniformly null table here means class-specific annotation noise does not
    exist, and the whole concern is answered.

    Returns
    -------
    DataFrame
    """
    rows = []
    continuous = continuous_covariates(cov)

    for g in groupings:
        cls = df[g]
        for metric in continuous:
            v = pd.to_numeric(cov[metric], errors='coerce')
            r = v[(cls == 'retained') & v.notna()]
            e = v[(cls == 'eliminated') & v.notna()]
            if len(r) < 5 or len(e) < 5:
                continue
            _, p = mannwhitneyu(r, e, alternative='two-sided')
            rm = r.median()
            rows.append({
                'grouping': g, 'metric': metric, 'test': 'mannwhitney',
                'ret_n': len(r), 'elim_n': len(e),
                'ret_value': rm, 'elim_value': e.median(),
                'ratio': (e.median() / rm) if rm else np.nan,
                'p_value': p,
                # a difference here is expected on biological grounds and is
                # NOT evidence of a gene-model problem
                'divergence_confounded': metric in DIVERGENCE_CONFOUNDED,
            })

        if 'orf_complete' in cov.columns:
            oc = pd.to_numeric(cov['orf_complete'], errors='coerce')
            ok = oc.notna()
            comp = flag_true(oc)
            r_ok = int(((cls == 'retained') & ok & comp).sum())
            r_no = int(((cls == 'retained') & ok & ~comp).sum())
            e_ok = int(((cls == 'eliminated') & ok & comp).sum())
            e_no = int(((cls == 'eliminated') & ok & ~comp).sum())
            if (r_ok + r_no) >= 5 and (e_ok + e_no) >= 5:
                _, p = fisher_exact([[e_ok, e_no], [r_ok, r_no]])
                r_rate = r_ok / (r_ok + r_no)
                e_rate = e_ok / (e_ok + e_no)
                rows.append({
                    'grouping': g, 'metric': 'orf_complete_rate',
                    'test': 'fisher',
                    'ret_n': r_ok + r_no, 'elim_n': e_ok + e_no,
                    'ret_value': r_rate, 'elim_value': e_rate,
                    'ratio': (e_rate / r_rate) if r_rate else np.nan,
                    'p_value': p,
                    'divergence_confounded': False,
                })

    return pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# TEST 2 - restriction
# ---------------------------------------------------------------------------

def branch_keep(df, sp, omega_cap):
    """Genes usable for dN on one branch.

    branch_lib.usable_mask, the set Figure 5 uses: dN/dS <= cap and dS > 0.
    Until 2026-10-10 this was a local copy without dS > 0, so a gene with no
    substitution on the branch (dN = dS = 0, dN/dS wherever codeml's random
    start left it) passed or failed the cap by chance. TEST 3 calls this.
    """
    return usable_mask(df, sp, 'dN', omega_cap, apply_cap=True)


#: Metrics for TEST 2, with whether the dN/dS cap applies -- Figure 5's
#: convention: dN and dN/dS on the capped rows, dS on every gene, because the
#: cap is a filter on dS and testing dS after it is circular.
TEST2_METRICS = (('dN', True), ('dS', False), ('omega', True))


def test2_strata(df, cov, groupings, omega_cap, min_stratum_n=20):
    """
    Fold change per labelling x metric x species x stratum.

    Every species is reported, not just the focal one. The focal-derived
    filters should shrink the focal fold change if annotation noise is
    inflating it; if they shrink every species equally, the filter is
    selecting a different gene class rather than cleaning up annotation.

    q_value is BH within each (labelling, metric), across species x strata.
    The strata are nested subsets, so their tests are positively dependent,
    which BH tolerates.

    Returns
    -------
    DataFrame
    """
    strata = build_strata(cov)
    rows = []
    for g in groupings:
        cls = df[g]
        for metric, capped in TEST2_METRICS:
            for sp in SPECIES_ORDER:
                base = usable_mask(df, sp, metric, omega_cap, apply_cap=capped)
                vals = pd.to_numeric(df[f'{sp}_{metric}'], errors='coerce')
                for sname, sfilt in strata.items():
                    keep = base if sfilt is None else (base & sfilt.fillna(False))
                    r = vals[keep & (cls == 'retained')].dropna()
                    e = vals[keep & (cls == 'eliminated')].dropna()
                    c = compare(r, e)
                    rows.append({
                        'grouping': g, 'metric': metric, 'capped': capped,
                        'branch': sp, 'stratum': sname,
                        'ret_n': c['ret_n'], 'elim_n': c['elim_n'],
                        'ret_median': c['ret_median'],
                        'elim_median': c['elim_median'],
                        'fold_change': c['fold_change'],
                        'p_value': c['p_value'],
                        # A stratum can be cleaned down to a handful of
                        # eliminated genes. Its fold change is then noise, so
                        # it is flagged and drawn without a star.
                        'underpowered': c['elim_n'] < min_stratum_n})
    out = pd.DataFrame(rows)
    out['q_value'] = np.nan
    for _, block in out.groupby(['grouping', 'metric']):
        ok = block['p_value'].notna()
        if ok.any():
            out.loc[block.index[ok], 'q_value'] = benjamini_hochberg(
                block.loc[ok, 'p_value'].to_numpy(float))
    return out


def build_strata(cov):
    """
    Named boolean filters, from loosest to strictest.

    Only strata whose covariates are actually available are returned.
    """
    strata = {'all': None}
    parts = []
    for flag, label in (('orf_complete', 'complete_orf'),
                        ('high_cov', 'high_cov'),
                        ('long_aln', 'long_aln')):
        if flag in cov.columns and pd.to_numeric(
                cov[flag], errors='coerce').notna().any():
            f = flag_true(cov[flag])
            strata[label] = f
            parts.append(f)

    # 'clean' is only meaningful when it is stricter than any single filter
    if len(parts) >= 2:
        combined = parts[0]
        for f in parts[1:]:
            combined = combined & f
        strata['clean'] = combined
    return strata


# ---------------------------------------------------------------------------
# TEST 3 - regression
# ---------------------------------------------------------------------------

def ols_hc3(X, y):
    """
    OLS with HC3 heteroscedasticity-consistent standard errors.

    dN is right-skewed with variance growing in the mean even after logging, so
    classical SEs are anti-conservative. HC3 is the recommended small-sample
    choice (MacKinnon & White 1985).

    Parameters
    ----------
    X : ndarray (n, k), must already include an intercept column
    y : ndarray (n,)

    Returns
    -------
    dict with beta, se, t, p, n, k, r2
    """
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta

    # leverage, clipped so a near-singular point cannot blow up HC3
    h = np.einsum('ij,jk,ik->i', X, XtX_inv, X)
    h = np.clip(h, 0.0, 1 - 1e-8)

    omega = (resid / (1.0 - h)) ** 2
    cov = XtX_inv @ (X.T * omega) @ X @ XtX_inv
    se = np.sqrt(np.clip(np.diag(cov), 0, np.inf))

    with np.errstate(divide='ignore', invalid='ignore'):
        tstat = np.where(se > 0, beta / se, np.nan)
    dof = max(n - k, 1)
    p = 2 * t_dist.sf(np.abs(tstat), df=dof)

    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {'beta': beta, 'se': se, 't': tstat, 'p': p,
            'n': n, 'k': k, 'r2': r2}


def design_matrix(elim, cov_df, covariates, use_covariates):
    """
    Build [intercept, eliminated, (z-scored continuous), (ORF-type dummies)].

    Continuous covariates are z-scored so their coefficients are comparable and
    the intercept stays interpretable. ORF type is dummy-coded with 'complete'
    as the reference level, so each dummy reads as "penalty relative to a
    complete ORF".

    Returns
    -------
    (X ndarray, list of column names, boolean row mask of usable rows)
    """
    cols = {'intercept': np.ones(len(elim)), 'eliminated': elim.astype(float)}
    usable = np.ones(len(elim), dtype=bool)

    if use_covariates:
        for c in covariates:
            if c == 'orf_type':
                continue
            v = pd.to_numeric(cov_df[c], errors='coerce')
            usable &= v.notna().to_numpy()
            sd = v.std(ddof=0)
            cols[c] = ((v - v.mean()) / sd if sd and np.isfinite(sd) and sd > 0
                       else v * 0.0).to_numpy()

        if 'orf_type' in covariates and 'orf_type' in cov_df.columns:
            ot = cov_df['orf_type']
            usable &= ot.notna().to_numpy()
            present = [t for t in ORF_TYPES[1:] if (ot == t).any()]
            for t in present:
                cols[f'orf_{t}'] = ot.eq(t).astype(float).to_numpy()

    names = list(cols.keys())
    X = np.column_stack([np.asarray(cols[n], dtype=float) for n in names])
    return X, names, usable


def test3_regression(df, cov, groupings, omega_cap, pseudocount=None,
                     include_identity=False):
    """
    log(dN) ~ eliminated, with and without gene-model covariates, per species.

    Fitted for EVERY species using the SAME focal-derived covariates. The
    comparison of interest is not whether the eliminated coefficient shrinks in
    the focal species, but whether it shrinks MORE than in the other three. Only
    a focal-specific attenuation implicates focal annotation quality; a uniform
    attenuation means the covariates track gene class.

    Returns
    -------
    DataFrame of coefficient rows
    """
    available = continuous_covariates(cov)
    if 'orf_type' in cov.columns and cov['orf_type'].notna().any():
        available = ['orf_type'] + available

    if not include_identity:
        dropped = [c for c in available if c in DIVERGENCE_CONFOUNDED]
        available = [c for c in available if c not in DIVERGENCE_CONFOUNDED]
        if dropped:
            logging.info(f"test3: excluding divergence-confounded covariates "
                         f"{dropped} - they are downstream of dN, so adjusting "
                         f"for them would absorb the effect under test "
                         f"(--include-identity to override)")

    clean = [c for c in available if c in STRUCTURAL_CLEAN]
    if not clean:
        logging.warning("test3: NO structural gene-model covariate is "
                        "available (ORF completeness/length missing).")
        logging.warning("  The remaining covariates are only partly "
                        "independent of divergence, so attenuation here is a "
                        "weak upper bound on any annotation effect, not a "
                        "clean estimate. Do not lead a rebuttal with it.")
    if not available:
        logging.warning("test3: no covariates available; fitting the "
                        "unadjusted model only")

    rows = []
    for g in groupings:
        cls = df[g]
        for sp in SPECIES_ORDER:
            base = branch_keep(df, sp, omega_cap)
            dn = pd.to_numeric(df[f'{sp}_dN'], errors='coerce')

            for use_cov in (False, True):
                if use_cov and not available:
                    continue
                X, names, usable = design_matrix(
                    (cls == 'eliminated').to_numpy(), cov, available, use_cov)

                mask = base.to_numpy() & usable & dn.notna().to_numpy()
                # need both classes represented
                if mask.sum() < 30:
                    continue
                y_raw = dn[mask].to_numpy(dtype=float)
                elim_col = X[mask][:, names.index('eliminated')]
                if elim_col.sum() < 5 or (1 - elim_col).sum() < 5:
                    continue

                pc = pseudocount
                if pc is None:
                    nz = y_raw[y_raw > 0]
                    pc = (nz.min() / 2.0) if len(nz) else 1e-6
                y = np.log(y_raw + pc)

                Xm = X[mask]
                # drop any all-constant column (e.g. a dummy with no members
                # left after masking) apart from the intercept
                keepc = [0] + [i for i in range(1, Xm.shape[1])
                               if np.ptp(Xm[:, i]) > 0]
                Xm, names_m = Xm[:, keepc], [names[i] for i in keepc]

                fit = ols_hc3(Xm, y)
                idx = names_m.index('eliminated')
                rows.append({
                    'grouping': g, 'branch': sp,
                    'model': 'covariate_adjusted' if use_cov else 'unadjusted',
                    'n': fit['n'], 'r2': fit['r2'],
                    'pseudocount': pc,
                    'covariates': ','.join(
                        [c for c in names_m
                         if c not in ('intercept', 'eliminated')]) or '-',
                    'beta_eliminated': fit['beta'][idx],
                    'se_eliminated': fit['se'][idx],
                    't_eliminated': fit['t'][idx],
                    'p_eliminated': fit['p'][idx],
                    'fold_from_beta': float(np.exp(fit['beta'][idx])),
                })

    res = pd.DataFrame(rows)
    if len(res):
        # attenuation: how much of the unadjusted effect the covariates absorb
        piv = res.pivot_table(index=['grouping', 'branch'], columns='model',
                              values='beta_eliminated')
        if {'unadjusted', 'covariate_adjusted'}.issubset(piv.columns):
            att = (1 - piv['covariate_adjusted'] / piv['unadjusted']) * 100
            att = att.rename('pct_attenuation').reset_index()
            res = res.merge(att, on=['grouping', 'branch'], how='left')
    return res


# ---------------------------------------------------------------------------
# Figure S9 (S13 before the 2026-10-12 renumbering; file stem kept)
# ---------------------------------------------------------------------------

#: Strata in drawing order, and their tick labels. The full names are used
#: when they fit: with the two strata the TransDecoder evidence supports, each
#: bar has about 19 mm. With all five, a bar has about 8 mm and the full names
#: at 7 pt would run into their neighbours, so the five-character forms are
#: used and the legend notes spell them out. tick_labels() measures, and the
#: summary records which form was drawn.
STRATUM_ORDER = ('all', 'complete_orf', 'high_cov', 'long_aln', 'clean')
STRATUM_LABEL_LONG = {'all': 'All', 'complete_orf': 'Complete ORF',
                      'high_cov': 'High coverage', 'long_aln': 'Long alignment',
                      'clean': 'Clean'}
STRATUM_LABEL = {'all': 'All', 'complete_orf': 'ORF', 'high_cov': 'Cov',
                 'long_aln': 'Aln', 'clean': 'Clean'}

#: Tick-label forms, tried in order until every panel's labels clear their
#: neighbours: (stratum names, n format, description for the summary).
TICK_FORMS = ((STRATUM_LABEL_LONG, 'n = {:,}', 'full stratum names, "n = 112"'),
              (STRATUM_LABEL, 'n = {:,}', 'short stratum names, "n = 112"'),
              (STRATUM_LABEL, 'n={:,}', 'short stratum names, compact "n=112" '
                                        '(the spaced form did not fit)'))

#: Bar tones, composited onto white so the rule at 1 does not show through.
BAR_STRENGTH = 0.80
WEAK_STRENGTH = 0.25


def on_white(color, strength):
    """`color` at `strength` composited onto white, as an opaque RGB tuple."""
    r, g, b, _ = to_rgba(color)
    return tuple(strength * c + (1.0 - strength) for c in (r, g, b))


def tick_labels(fig, axes_grid, order, ns_grid, gap_pt=2.0):
    """Set the stratum tick labels in the first form of TICK_FORMS whose
    labels clear their neighbours (by `gap_pt`) in every panel; the last form
    is kept if none does, with a warning. `ns_grid[ri][ci]` lists the
    eliminated n per stratum for that panel. Returns the form's description.
    """
    gap_px = gap_pt / 72 * fig.dpi
    for k, (names, nfmt, desc) in enumerate(TICK_FORMS):
        for ri, row in enumerate(axes_grid):
            for ci, ax in enumerate(row):
                ax.set_xticklabels([f'{names[st]}\n{nfmt.format(n)}'
                                    for st, n in zip(order, ns_grid[ri][ci])])
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        clash = False
        for row in axes_grid:
            for ax in row:
                boxes = sorted((t.get_window_extent(r) for t in
                                ax.get_xticklabels() if t.get_text()),
                               key=lambda b: b.x0)
                if any(b.x0 - a.x1 < gap_px for a, b in zip(boxes, boxes[1:])):
                    clash = True
        if not clash:
            return desc
        if k == len(TICK_FORMS) - 1:
            logging.warning('  tick labels still touch in the most compact form;'
                            ' check the figure')
            return desc + ' -- STILL TOUCHING'
    return desc


def plot_strata(strata_df, grouping, outdir, stem, width='double',
                height_mm=104.0):
    """
    Supplementary Figure S9 (S13 before renumbering) for one region
    labelling: eliminated/retained fold change by gene-model stratum, dN on
    top and dS below, one column per branch.

    The strata are defined by FOCAL-species gene models, so a focal column whose
    bars fall while the other three hold is the signature of annotation-driven
    inflation; bars that move together in every column mean the filter selects
    a different gene class. Bars are in ink -- a fold change is a metric, not a
    region. A stratum with too few eliminated genes is hatched light grey and
    carries no star. n = eliminated genes, in the tick labels ("n = 112";
    tick_labels() picks the widest form that fits).

    Each row's y-axis is shared by the four columns AND by the figures of
    every labelling in `strata_df`, because the two figures are read species
    by species; the limit comes from the tallest bar of the metric in any
    labelling.

    Returns (saved paths, facts for the summary), or ([], {}) if the
    labelling has no rows.
    """
    sub = strata_df[strata_df['grouping'] == grouping]
    order = [s for s in STRATUM_ORDER if s in set(sub['stratum'])]
    if not len(sub) or not order:
        return [], {}

    apply_style()
    rows = (('dN', 'dN fold change\n(eliminated / retained)'),
            ('dS', 'dS fold change\n(eliminated / retained)'))
    fig, axes = plt.subplots(len(rows), len(SPECIES_ORDER),
                             figsize=(COL[width], mm(height_mm)),
                             sharey='row', layout='none')
    x = np.arange(len(order))
    facts = {'ylim': {}, 'n_eliminated': {}}
    ns_grid = []
    for ri, (metric, ylabel) in enumerate(rows):
        msub = sub[sub['metric'] == metric]
        # The y-limit from every labelling, so the figures share their axes.
        finite = strata_df.loc[strata_df['metric'] == metric,
                               'fold_change'].to_numpy(float)
        finite = finite[np.isfinite(finite)]
        # Room above the tallest bar for its value and star.
        ymax = max(float(finite.max()) * 1.32, 1.35) if len(finite) else 2.0
        facts['ylim'][metric] = round(ymax, 3)
        ns_grid.append([])
        for ci, sp in enumerate(SPECIES_ORDER):
            ax = axes[ri, ci]
            s = msub[msub['branch'] == sp].set_index('stratum')

            def get(st, col, default):
                return s.loc[st, col] if st in s.index else default

            folds = [float(get(st, 'fold_change', np.nan)) for st in order]
            weak = [bool(get(st, 'underpowered', True)) for st in order]
            ns = [int(get(st, 'elim_n', 0)) for st in order]
            qs = [float(get(st, 'q_value', np.nan)) for st in order]
            ns_grid[ri].append(ns)
            facts['n_eliminated'][f'{metric} {sp}'] = dict(zip(order, ns))

            face = [on_white(NEUTRAL['ink'], WEAK_STRENGTH if w else BAR_STRENGTH)
                    for w in weak]
            bars = ax.bar(x, np.nan_to_num(folds, nan=0.0), 0.68, color=face,
                          edgecolor='black', linewidth=0.6, zorder=2)
            for b, w in zip(bars, weak):
                if w:
                    b.set_hatch('////')
            ax.axhline(1.0, color=NEUTRAL['rule'], ls='--', lw=0.9, zorder=1)

            for xi in range(len(order)):
                if not np.isfinite(folds[xi]):
                    continue
                # White box: a value just under 1 would otherwise sit on, and
                # be struck through by, the rule at 1.
                ax.text(xi, folds[xi] + 0.02 * ymax, f'{folds[xi]:.2f}',
                        ha='center', va='bottom', fontsize=FONT['annotation'],
                        zorder=4, bbox=dict(facecolor='white', edgecolor='none',
                                            pad=0.4))
                if not weak[xi]:
                    ax.text(xi, folds[xi] + 0.11 * ymax, stars(qs[xi]),
                            ha='center', va='bottom',
                            fontsize=FONT['annotation'])

            ax.set_xticks(x)                 # labels: tick_labels(), below
            ax.set_xlim(-0.6, len(order) - 0.4)
            ax.set_ylim(0, ymax)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            if ri == 0:
                ax.set_title(SPECIES_LABEL[SPECIES_KEY[sp]],
                             fontsize=FONT['panel_title'], pad=4)
            if ci == 0:
                ax.set_ylabel(ylabel)

    fig.subplots_adjust(left=0.085, right=0.995, top=0.95, bottom=0.12,
                        wspace=0.12, hspace=0.45)
    facts['tick_labels'] = tick_labels(fig, axes, order, ns_grid)
    logging.info('  tick labels: %s', facts['tick_labels'])
    fit_bottom_margin(fig)
    facts['audit'] = list(audit(fig))
    paths = save(fig, outdir, stem)
    plt.close(fig)
    logging.info('  wrote %s', ', '.join(os.path.basename(p) for p in paths))
    return paths, facts


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def p_text(p):
    """Plain-text P or q."""
    if p is None or not np.isfinite(p):
        return 'NA'
    return f'{p:.2g}' if p < 0.001 else f'{p:.3f}'


def text_numbers(t1, t2, t3, focal, groupings):
    """The numbers Results and Methods quote, as summary lines."""
    L = []
    oc = t1[t1['metric'] == 'orf_complete_rate'] if len(t1) else t1
    for _, r in oc.iterrows():
        L.append(f"  ORF completeness [{r['grouping']}]: eliminated "
                 f"{100 * r['elim_value']:.1f}% (n = {int(r['elim_n']):,}) vs "
                 f"retained {100 * r['ret_value']:.1f}% (n = {int(r['ret_n']):,}),"
                 f" Fisher P = {p_text(r['p_value'])}")

    f2 = t2[(t2['branch'] == focal)] if len(t2) else t2
    for metric in ('dN', 'dS'):
        for g in groupings:
            s = f2[(f2['metric'] == metric) & (f2['grouping'] == g)]
            if not len(s):
                continue
            parts = []
            for st in STRATUM_ORDER:
                r = s[s['stratum'] == st]
                if not len(r):
                    continue
                r = r.iloc[0]
                flag = ' [underpowered]' if r['underpowered'] else ''
                parts.append(f"{st} {r['fold_change']:.2f} (n = "
                             f"{int(r['elim_n'])}, q = {p_text(r['q_value'])})"
                             f"{flag}")
            L.append(f'  {focal} {metric} fold [{g}]: ' + '; '.join(parts))

    co = f2[(f2['metric'] == 'dN') & (f2['stratum'] == 'complete_orf')
            & ~f2['underpowered']] if len(f2) else f2
    if len(co):
        lo, hi = co['fold_change'].min(), co['fold_change'].max()
        L.append(f'  -> complete-ORF {focal} dN fold across labellings: '
                 f'{lo:.2f} to {hi:.2f} (Results: "persists at ... fold")')
    ds_all = f2[(f2['metric'] == 'dS') & (f2['stratum'] == 'all')]
    ds_co = f2[(f2['metric'] == 'dS') & (f2['stratum'] == 'complete_orf')]
    if len(ds_all) and len(ds_co):
        L.append(f'  -> {focal} dS excess, all genes vs complete ORFs: ' + '; '.join(
            f"[{g}] {ds_all[ds_all['grouping'] == g]['fold_change'].iloc[0]:.2f}"
            f" (q {p_text(ds_all[ds_all['grouping'] == g]['q_value'].iloc[0])})"
            f" -> {ds_co[ds_co['grouping'] == g]['fold_change'].iloc[0]:.2f}"
            f" (q {p_text(ds_co[ds_co['grouping'] == g]['q_value'].iloc[0])})"
            for g in groupings
            if (ds_all['grouping'] == g).any() and (ds_co['grouping'] == g).any()))
        L.append('     A dS excess that shrinks to ~1 in complete ORFs while dN '
                 'holds supports the')
        L.append('     gene-model-noise sentence; one that holds does not.')

    if len(t3) and 'pct_attenuation' in t3.columns:
        att = (t3.dropna(subset=['pct_attenuation'])
                 .drop_duplicates(['grouping', 'branch']))
        for g in groupings:
            a = att[att['grouping'] == g]
            if not len(a) or not (a['branch'] == focal).any():
                continue
            fa = float(a.loc[a['branch'] == focal, 'pct_attenuation'].iloc[0])
            others = a.loc[a['branch'] != focal, 'pct_attenuation']
            L.append(f'  TEST 3 control [{g}]: {focal} attenuation {fa:+.1f}% vs '
                     f'{others.min():+.1f}% to {others.max():+.1f}% in the other '
                     'three')
    return L


def legend_notes(args, focal, fig_facts=None):
    """Lines for the S9 legend (S13 before renumbering), with what each
    figure file actually drew."""
    L = [
        f'Bars: eliminated / retained median per stratum; strata defined by '
        f'{focal} gene models.',
        f'dN: branch dN/dS <= {args.omega_cap:g} and dS > 0 on that branch '
        '(Figure 5 masking); dS: every gene.',
        'n under each stratum = ELIMINATED genes in the stratum (retained n '
        'per cell: test2_strata.tsv).',
        f'Hatched light grey: fewer than {args.min_stratum_n} eliminated genes; '
        'not interpreted, no star.',
        'Asterisks: BH-FDR q within each labelling and metric, across branches '
        'x strata.',
        'Dashed line: no difference (fold change 1).',
        'Strata: All = every usable gene; Complete ORF (short form "ORF") = '
        'complete TransDecoder ORFs only;',
        f'High coverage ("Cov") = alignment or BLAST coverage >= '
        f'{args.min_coverage:g}; Long alignment ("Aln") = alignment length >= '
        '--min-align-len (if set);',
        'Clean = every available filter at once. Only strata with data are '
        'drawn.']
    for g, f in (fig_facts or {}).items():
        if not f:
            L.append(f'figS13_fold_by_stratum_{g}: no rows, not drawn.')
            continue
        lim = ', '.join(f'{m} 0-{v:g}' for m, v in f['ylim'].items())
        L.append(f'figS13_fold_by_stratum_{g} ({g} labelling): tick labels '
                 f'{f["tick_labels"]}; y-axes {lim}, shared with the other '
                 'labelling\'s figure. Nothing on the figure names its '
                 'labelling, so the caption must.')
        if f.get('audit'):
            L.append('  AUDIT: ' + '; '.join(f['audit']))
    return L


def write_summary(path, t1, t2, t3, focal, groupings, match_stats, args,
                  fig_facts=None):
    """Human-readable verdict for each test, the numbers for the text, and the
    S9 legend notes (S13 before renumbering)."""
    with open(path, 'w') as fh:
        fh.write("Gene-model diagnostic for the elevated rates on the "
                 f"{focal} branch\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Annotated table : {args.annotated}\n")
        fh.write(f"Focal species   : {focal}\n")
        fh.write(f"dN/dS cap       : {args.omega_cap:g} (dS > 0 required on "
                 "the capped rows)\n")
        fh.write(f"Coverage cutoff : {args.min_coverage:g}\n")
        fh.write(f"Genes           : {match_stats.get('n_genes', 0):,}\n")
        for src in ('transdecoder', 'alignments', 'rbh'):
            fh.write(f"  matched from {src:<14}: "
                     f"{match_stats.get(src + '_matched', 0):,}\n")
        fh.write("\n")

        fh.write("-" * 70 + "\nNumbers for the text\n" + "-" * 70 + "\n")
        fh.write('\n'.join(text_numbers(t1, t2, t3, focal, groupings)) + "\n\n")

        fh.write("-" * 70 + "\n")
        fh.write("TEST 1  PREMISE - do eliminated genes have worse gene "
                 "models?\n")
        fh.write("-" * 70 + "\n")
        if len(t1):
            fh.write(t1.to_string(index=False))
            fh.write("\n\n")
            if t1.get('divergence_confounded', pd.Series(dtype=bool)).any():
                fh.write("NOTE: metrics flagged divergence_confounded "
                         "(percent identity, bitscore) are\nexpected to be "
                         "lower in eliminated genes because those genes really "
                         "are\nfaster-evolving - that is the finding under "
                         "test, not a gene-model defect.\nJudge annotation "
                         "quality from the structural metrics (ORF "
                         "completeness,\nORF length) and from coverage and "
                         "alignment length.\n\n")
            sig = t1[(t1['p_value'] < 0.05)
                     & ~t1.get('divergence_confounded',
                               pd.Series(False, index=t1.index))]
            if not len(sig):
                fh.write("VERDICT: no gene-model quality metric differs between "
                         "eliminated and\nretained genes, so class-specific "
                         "annotation noise cannot be inflating\ntheir rates. "
                         "Tests 2 and 3 are confirmatory.\n\n")
            else:
                fh.write("VERDICT: gene-model quality DOES differ by class for "
                         f"{len(sig)} metric(s)\nbelow. The premise of the "
                         "concern holds, so tests 2 and 3 decide whether\nit "
                         "actually drives the result.\n")
                fh.write(sig[['grouping', 'metric', 'ret_value', 'elim_value',
                              'ratio', 'p_value']].to_string(index=False))
                fh.write("\n\n")
        else:
            fh.write("No covariates available - supply --transdecoder, "
                     "--aln-dir or --rbh.\n\n")

        fh.write("-" * 70 + "\n")
        fh.write("TEST 2  RESTRICTION - does the fold change hold in clean "
                 "genes?\n")
        fh.write("-" * 70 + "\n")
        if len(t2):
            for metric, _ in TEST2_METRICS:
                for g in groupings:
                    piv = (t2[(t2['grouping'] == g) & (t2['metric'] == metric)]
                           .pivot_table(index='branch', columns='stratum',
                                        values='fold_change')
                           .reindex(SPECIES_ORDER))
                    piv = piv[[s for s in STRATUM_ORDER if s in piv.columns]]
                    label = 'dN/dS' if metric == 'omega' else metric
                    fh.write(f"\n[{g}] {label} fold change by stratum\n")
                    fh.write(piv.to_string(float_format=lambda v: f'{v:.3f}'))
                    fh.write("\n")
            fh.write("\nRead DOWN the focal row for the effect of cleaning, and "
                     "ACROSS the other\nspecies for whether the same filter "
                     "moves them too. A focal-only drop is\nannotation; a "
                     "uniform drop means the filter selected a different gene "
                     "class.\nn, medians, P and q per cell: test2_strata.tsv.\n\n")

        fh.write("-" * 70 + "\n")
        fh.write("TEST 3  REGRESSION - log(dN) ~ eliminated (+ covariates)\n")
        fh.write("-" * 70 + "\n")
        used = set()
        if len(t3) and 'covariates' in t3.columns:
            for s in t3['covariates'].dropna():
                used.update(x for x in str(s).split(',') if x and x != '-')
        if not (used & set(STRUCTURAL_CLEAN)) and not any(
                u.startswith('orf_') for u in used):
            fh.write("WARNING: no structural gene-model covariate (ORF "
                     "completeness or ORF\nlength) was available, so this "
                     "regression adjusts only for covariates that\nare "
                     "themselves partly a function of divergence. Any "
                     "attenuation below is an\nUPPER BOUND on a possible "
                     "annotation effect, not an estimate of one.\n\n")
        if len(t3):
            show = ['grouping', 'branch', 'model', 'n', 'beta_eliminated',
                    'se_eliminated', 'p_eliminated', 'fold_from_beta']
            if 'pct_attenuation' in t3.columns:
                show.append('pct_attenuation')
            fh.write(t3[show].to_string(index=False))
            fh.write("\n\nCONTROL: compare the focal attenuation with the "
                     "other three species.\nSimilar attenuation across all "
                     "four means the covariates are proxying gene\nclass, "
                     "not focal annotation quality, and the adjustment "
                     "carries no\nevidence about gene models.\n")
        else:
            fh.write("Not enough data to fit.\n")
        fh.write("\n")

        fh.write("-" * 70 + "\nSupplementary figure S9 legend notes (S13 "
                 "before renumbering)\n" + "-" * 70 + "\n")
        fh.write('\n'.join(legend_notes(args, focal, fig_facts)) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Test whether elevated branch rates on the focal species "
                    "are driven by gene-model quality (supplementary fig. S9, "
                    "S13 before renumbering; Table T7).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotated', required=True,
                   help='branch_dnds_four_annotated.tsv from a 2026-10-10 or '
                        'later branch_dnds_four.py run')
    p.add_argument('--focal-species', default='Toxocara',
                   choices=SPECIES_ORDER,
                   help='Species whose gene models are under suspicion '
                        '(default: Toxocara)')
    p.add_argument('--focal-id-col', default=None,
                   help='Column in the annotated table holding focal-species '
                        'gene IDs. Auto-detected if omitted.')
    p.add_argument('--transdecoder', default=None,
                   help='TransDecoder .pep or .cds FASTA for the focal species')
    p.add_argument('--aln-dir', default=None,
                   help='Directory of per-gene codeml input alignments')
    p.add_argument('--aln-pattern', default='*.fa,*.fasta,*.aln,*.phy,*.nuc',
                   help='Comma-separated glob patterns for alignment files')
    p.add_argument('--aln-focal-pattern', default=None,
                   help='Substring identifying the focal sequence inside an '
                        'alignment (default: the focal species name)')
    p.add_argument('--aln-id-strip', default=None,
                   help='Comma-separated extra filename suffixes to strip when '
                        'deriving a gene ID from an alignment filename')
    p.add_argument('--rbh', nargs='*', default=None,
                   help='One or more reciprocal-best-hit TSVs')
    p.add_argument('--transdecoder-id-col', default=None,
                   help='Force the annotated-table column to join the '
                        'TransDecoder FASTA on (default: best match)')
    p.add_argument('--alignments-id-col', default=None,
                   help='Force the annotated-table column to join alignments '
                        'on (default: best match)')
    p.add_argument('--rbh-id-col', default=None,
                   help='Force the annotated-table column to join RBH tables '
                        'on (default: best match)')
    p.add_argument('--no-strip-orf-suffix', action='store_true',
                   help="Do not strip TransDecoder '.p1' suffixes when "
                        "matching IDs")
    p.add_argument('--lower-ids', action='store_true',
                   help='Lowercase all IDs before matching. Use only after '
                        '--inspect confirms a benign case mismatch; BLAST '
                        'accession normalisation has caused silent ID loss in '
                        'this pipeline before.')
    p.add_argument('--omega-cap', type=float, default=5.0,
                   help='Per-branch dN/dS cap, matching branch_dnds_four.py '
                        '(default: 5)')
    p.add_argument('--min-coverage', type=float, default=0.7,
                   help='Alignment coverage threshold for the high_cov stratum '
                        '(default: 0.7)')
    p.add_argument('--include-identity', action='store_true',
                   help='Allow BLAST percent identity and bitscore into the '
                        'TEST 3 regression. OFF by default: these are '
                        'downstream of the divergence being measured, so '
                        'adjusting for them absorbs the very effect under test.')
    p.add_argument('--table-covariates', default='auto',
                   help="Numeric columns already in the annotated table to use "
                        "as covariates, comma-separated. 'auto' (default) picks "
                        "up known alignment-quality columns such as align_len; "
                        "'none' disables.")
    p.add_argument('--min-align-len', type=int, default=0,
                   help='If >0, add a long_aln stratum of genes whose alignment '
                        'is at least this many columns (default: 0, off)')
    p.add_argument('--min-stratum-n', type=int, default=20,
                   help='Eliminated-gene count below which a stratum is marked '
                        'underpowered: reported, but hatched and unstarred '
                        '(default: 20)')
    p.add_argument('--pseudocount', type=float, default=None,
                   help='Added to dN before logging. Default: half the smallest '
                        'nonzero dN in each fit.')
    p.add_argument('--groupings', default=None,
                   help='Comma-separated region columns. Default: auto-detect '
                        'as_region and tox_region.')
    p.add_argument('--inspect', action='store_true',
                   help='Report ID match rates from every source and exit. RUN '
                        'THIS FIRST.')
    p.add_argument('--outdir', default=None)
    args = p.parse_args()
    if not 0 < args.min_coverage <= 1:
        p.error('--min-coverage must be in (0, 1].')
    if args.omega_cap <= 0:
        p.error('--omega-cap must be positive.')
    if args.outdir is None:
        args.outdir = make_outdir('toxocara_annotation_diagnostic',
                                  base='results')
    return args


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir, 'diagnostic')
    for noisy in ('matplotlib.font_manager', 'matplotlib.texmanager',
                  'fontTools', 'fontTools.subset', 'fontTools.ttLib', 'PIL'):
        logging.getLogger(noisy).setLevel(logging.ERROR)

    logging.info("=" * 60)
    logging.info(f"Gene-model diagnostic - focal species {args.focal_species}")
    logging.info("=" * 60)
    logging.info(f"outdir: {os.path.abspath(args.outdir)}")

    if not (args.transdecoder or args.aln_dir or args.rbh):
        # Not fatal: the annotated table may itself carry usable covariates
        # (align_len). Whether anything usable exists is checked below.
        logging.warning("No external evidence source given (--transdecoder, "
                        "--aln-dir, --rbh).")
        logging.warning("Falling back to covariates already in the annotated "
                        "table. TEST 1 will cover alignment length but NOT ORF "
                        "completeness, which is the most informative input.")

    # Validates the rate columns and refuses a rooted-tree table.
    groupings = ([g.strip() for g in args.groupings.split(',')]
                 if args.groupings else None)
    df, groupings = load_annotated(args.annotated, groupings)
    if 'status' in df.columns:
        bad = df['status'].astype(str) != 'OK'
        if bad.any():
            logging.warning(f"  {int(bad.sum()):,} rows with a status other "
                            "than OK dropped")
            df = df[~bad].reset_index(drop=True)
    logging.info(f"Groupings: {groupings}")

    focal_id_col = find_focal_id_col(df, args.focal_species, args.focal_id_col)
    if focal_id_col is None:
        sys.exit(f"ERROR: could not find a {args.focal_species} gene-ID column. "
                 f"Pass --focal-id-col. Available columns:\n"
                 f"  {list(df.columns)}")
    logging.info(f"Focal ID column: '{focal_id_col}'  e.g. "
                 f"{list(df[focal_id_col].dropna().astype(str)[:3])}")

    strip = not args.no_strip_orf_suffix

    td = aln = rbh = None
    if args.transdecoder:
        if not os.path.isfile(args.transdecoder):
            sys.exit(f"ERROR: --transdecoder not found: {args.transdecoder}")
        td = parse_transdecoder(args.transdecoder, strip_suffix=strip)
    if args.aln_dir:
        if not os.path.isdir(args.aln_dir):
            sys.exit(f"ERROR: --aln-dir not found: {args.aln_dir}")
        aln = parse_alignments(
            args.aln_dir,
            [p.strip() for p in args.aln_pattern.split(',') if p.strip()],
            args.aln_focal_pattern or args.focal_species,
            strip_suffix=strip,
            id_strip=([s.strip() for s in args.aln_id_strip.split(',')]
                      if args.aln_id_strip else None))
    if args.rbh:
        missing = [f for f in args.rbh if not os.path.isfile(f)]
        if missing:
            sys.exit(f"ERROR: --rbh file(s) not found: {missing}")
        rbh = parse_rbh(args.rbh, strip_suffix=strip)

    sources = {'transdecoder': td, 'alignments': aln, 'rbh': rbh}
    overrides = {'transdecoder': args.transdecoder_id_col,
                 'alignments': args.alignments_id_col,
                 'rbh': args.rbh_id_col}
    cov, match_stats = build_covariates(df, focal_id_col, sources,
                                        args.min_coverage,
                                        strip_suffix=strip,
                                        lower_ids=args.lower_ids,
                                        overrides=overrides)

    # Covariates that are already columns of the annotated table.
    if args.table_covariates.lower() != 'none':
        if args.table_covariates.lower() == 'auto':
            tcols = [c for c in TABLE_COVARIATE_CANDIDATES if c in df.columns]
        else:
            tcols = [c.strip() for c in args.table_covariates.split(',')
                     if c.strip()]
            bad = [c for c in tcols if c not in df.columns]
            if bad:
                sys.exit(f"ERROR: --table-covariates not in the annotated "
                         f"table: {bad}. Available: {list(df.columns)}")
        for c in tcols:
            v = pd.to_numeric(df[c], errors='coerce')
            if not v.notna().any():
                logging.warning(f"table covariate '{c}' has no numeric values "
                                f"- skipping")
                continue
            name = c if c not in cov.columns else f'tbl_{c}'
            cov[name] = v.to_numpy()
            logging.info(f"Table covariate '{c}' -> {name} "
                         f"({int(v.notna().sum())}/{len(df)} finite)")

    if args.min_align_len > 0:
        alncol = next((c for c in ('align_len', 'aln_len', 'tbl_align_len',
                                   'alignment_len', 'alignment_length')
                       if c in cov.columns), None)
        if alncol is None:
            logging.warning("--min-align-len set but no alignment-length "
                            "column is available - ignoring")
        else:
            v = pd.to_numeric(cov[alncol], errors='coerce')
            cov['long_aln'] = as_flag(v >= args.min_align_len, v.isna())
            logging.info(f"long_aln stratum: {alncol} >= {args.min_align_len} "
                         f"-> {int(flag_true(cov['long_aln']).sum())} genes")

    if not continuous_covariates(cov) and 'orf_type' not in cov.columns:
        sys.exit(
            "ERROR: no usable covariates were found.\n"
            "  Nothing matched from --transdecoder / --aln-dir / --rbh, and the\n"
            "  annotated table has no numeric alignment-quality column either.\n"
            f"  Annotated table columns: {list(df.columns)}\n"
            "  Supply an evidence source, or name a column explicitly with\n"
            "  --table-covariates <col>.")

    if args.inspect:
        logging.info("")
        logging.info("--- INSPECT ---")
        logging.info(f"focal ID column '{focal_id_col}' e.g. "
                     f"{list(df[focal_id_col].dropna().astype(str)[:4])}")
        for name, src in sources.items():
            if src is None or not len(src):
                logging.info(f"{name:<14} not supplied")
                continue
            col, rate = match_stats['joins'].get(name, (None, 0.0))
            m = match_stats.get(f'{name}_matched', 0)
            logging.info(f"{name:<14} IDs e.g. "
                         f"{list(src['focal_id'].dropna().astype(str)[:4])}")
            logging.info(f"{'':<14} best join column: '{col}' "
                         f"(match rate {100 * rate:.1f}%) -> "
                         f"{m}/{len(df)} genes")
        logging.info(f"Covariate columns produced: {list(cov.columns)}")
        logging.info("A match rate near zero means the IDs differ in form. "
                     "Compare the examples above,")
        logging.info("then fix with --focal-id-col, --alignments-id-col, "
                     "--aln-id-strip,")
        logging.info("--no-strip-orf-suffix or --lower-ids before running for "
                     "real.")
        logging.info("--- inspect only, exiting ---")
        return

    for name, (col, rate) in match_stats.get('joins', {}).items():
        logging.info(f"join: {name} <- '{col}' ({100 * rate:.1f}%)")

    cov_path = os.path.join(args.outdir, 'annotation_covariates.tsv')
    pd.concat([df[[focal_id_col] + groupings].reset_index(drop=True),
               cov.reset_index(drop=True)],
              axis=1).to_csv(cov_path, sep='\t', index=False)
    logging.info(f"Covariates: {cov_path}")

    logging.info("\n=== TEST 1: premise ===")
    t1 = test1_premise(df, cov, groupings)
    if len(t1):
        for _, r in t1.iterrows():
            logging.info(f"  [{r['grouping']}] {r['metric']:<18} "
                         f"ret={r['ret_value']:.4g} elim={r['elim_value']:.4g} "
                         f"ratio={r['ratio']:.3f} p={r['p_value']:.3g} "
                         f"{stars(r['p_value'])}")
    t1.to_csv(os.path.join(args.outdir, 'test1_premise.tsv'),
              sep='\t', index=False)

    logging.info("\n=== TEST 2: restriction ===")
    t2 = test2_strata(df, cov, groupings, args.omega_cap,
                      min_stratum_n=args.min_stratum_n)
    weak = t2[t2['underpowered'] & t2['stratum'].ne('all')
              & t2['metric'].eq('dN')]
    if len(weak):
        logging.warning(f"  {len(weak)} dN stratum/species cells have fewer "
                        f"than {args.min_stratum_n} eliminated genes and are "
                        "not interpretable:")
        for _, r in weak.iterrows():
            logging.warning(f"    [{r['grouping']}] {r['branch']} / "
                            f"{r['stratum']}: n_elim={int(r['elim_n'])}")
    for metric in ('dN', 'dS'):
        for g in groupings:
            piv = (t2[(t2['grouping'] == g) & (t2['metric'] == metric)]
                   .pivot_table(index='branch', columns='stratum',
                                values='fold_change').reindex(SPECIES_ORDER))
            logging.info(f"  [{g}] {metric} fold change by stratum\n"
                         f"{piv.to_string()}")
    t2.to_csv(os.path.join(args.outdir, 'test2_strata.tsv'),
              sep='\t', index=False, float_format='%.6g')

    logging.info("\n=== TEST 3: regression ===")
    t3 = test3_regression(df, cov, groupings, args.omega_cap, args.pseudocount,
                          include_identity=args.include_identity)
    if len(t3):
        for _, r in t3.iterrows():
            logging.info(f"  [{r['grouping']}] {r['branch']:<14} "
                         f"{r['model']:<19} beta={r['beta_eliminated']:+.4f} "
                         f"(SE {r['se_eliminated']:.4f}) "
                         f"p={r['p_eliminated']:.3g} n={int(r['n'])}")
        if 'pct_attenuation' in t3.columns:
            att = (t3.dropna(subset=['pct_attenuation'])
                     .drop_duplicates(['grouping', 'branch']))
            for _, r in att.iterrows():
                logging.info(f"  attenuation [{r['grouping']}] "
                             f"{r['branch']:<14} {r['pct_attenuation']:+.1f}%")
    t3.to_csv(os.path.join(args.outdir, 'test3_regression.tsv'),
              sep='\t', index=False, float_format='%.6g')

    logging.info("\n=== Figure S9 (S13 before renumbering) ===")
    fig_facts = {}
    for g in groupings:
        _, fig_facts[g] = plot_strata(t2, g, args.outdir,
                                      f'figS13_fold_by_stratum_{g}')

    summary_path = os.path.join(args.outdir, 'diagnostic_summary.txt')
    write_summary(summary_path, t1, t2, t3, args.focal_species, groupings,
                  match_stats, args, fig_facts)
    for line in text_numbers(t1, t2, t3, args.focal_species, groupings):
        logging.info(line)

    with open(os.path.join(args.outdir, 'provenance.json'), 'w') as fh:
        json.dump({'script': os.path.basename(__file__),
                   'git_commit': git_commit(_HERE),
                   'run_utc': datetime.now(timezone.utc).replace(
                       microsecond=0).isoformat(),
                   'host': platform.node(),
                   'python': platform.python_version(),
                   'argv': sys.argv, 'n_genes': int(len(df)),
                   'groupings': groupings,
                   'match_stats': {k: v for k, v in match_stats.items()
                                   if k != 'joins'},
                   'joins': {k: list(v) for k, v in
                             match_stats.get('joins', {}).items()},
                   'inputs': {'annotated': fingerprint(args.annotated),
                              'transdecoder': fingerprint(args.transdecoder),
                              'rbh': [fingerprint(f) for f in (args.rbh or [])]},
                   'figures': fig_facts},
                  fh, indent=2, default=str)
    logging.info(f"\nSummary: {summary_path}")
    logging.info("Done.")


if __name__ == '__main__':
    main()
