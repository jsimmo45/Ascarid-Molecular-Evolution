#!/usr/bin/env python3
"""
Molecular-clock tests for the four-way ascarid quartets.

WHY THIS EXISTS
---------------
The four-way branch model gives terminal branch lengths that differ by a factor
of 3.25 between Ascaris (0.0461) and Parascaris (0.1496), even though those two
branches subtend exactly the same elapsed time. If that rate heterogeneity is
real, then no single substitution rate can be applied across these lineages,
pairwise Ks cannot be converted to divergence times, and the sensible thing is
to report rates rather than dates.

This script turns that observation into a formal result, two independent ways:

  TEST 1 - Global-clock likelihood ratio test (PAML baseml).
      Fit the concatenated alignment twice, once with clock = 0 (every branch
      free) and once with clock = 1 (a single rate for the whole tree), then
      compare with an LRT. Degrees of freedom are taken from PAML's own
      parameter counts rather than assumed, and cross-checked against the
      theoretical n - 2.

  TEST 2 - Tajima's 1D relative-rate test (Tajima 1993, Genetics 135:599).
      For a sister pair and an outgroup, count sites where each sister differs
      from the outgroup while the other matches it. Under equal rates the two
      counts are equal in expectation. This is model-free: no substitution
      model, no branch-length estimation, nothing to misspecify. It is the
      right companion to the LRT precisely because it can fail for entirely
      different reasons.

Both tests run off one concatenated supermatrix built from the per-gene codon
alignments already produced by branch_dnds_four.py.

USAGE
-----
Always run --inspect first. Sequence headers in the alignment files are gene
IDs, not species names, and the script has to learn which is which:

    python molecular_clock_test.py --aln-dir four_way_alignments --inspect

Then the real run:

    python molecular_clock_test.py \\
      --aln-dir four_way_alignments \\
      --outdir clock_test \\
      --sites 12

SITE SUBSETS (--sites)
----------------------
  all  every aligned nucleotide position
  12   codon positions 1 and 2 only  -- RECOMMENDED
  3    codon position 3 only

Position 3 carries most of the synonymous signal and is where the Toxocara
branch is saturated (median Ks 1.39-1.43). Running the clock test on positions
1 and 2 keeps the test away from the saturated sites, so a rejection cannot be
dismissed as a saturation artifact. Run 'all' as well if you want to show the
answer does not depend on the choice.

OUTPUT
------
  supermatrix.phy          concatenated alignment, PHYLIP sequential
  tree.trees               topology with PAML's "<ntax> <ntrees>" header
  clock0/, clock1/         PAML working directories (baseml writes several
                           files into its cwd, so each fit gets its own)
  clock_test_results.tsv   machine-readable results for both tests
  clock_test_summary.txt   human-readable summary with the interpretation
  clock_test.log

GOTCHAS BAKED IN
----------------
  * PAML returns exit code 1 on a successful run. Success is decided by whether
    the output file exists and contains a parsable lnL, never by the exit code.
  * PAML tree files need a "<ntaxa> <ntrees>" first line (the "3 1" header in
    the lab notes -- it is "4 1" here because there are four taxa).
  * clock = 1 requires a ROOTED tree and method = 0. PAML's method = 1
    (one branch at a time) silently does not support clock models.
  * baseml scatters rst / rub / lnf / mlb into its working directory, which is
    why each fit runs in its own subdirectory.

Author: generated for the Ka/Ks ascarid PDE project.
"""

import argparse
import glob
import logging
import os
import re
import subprocess
import sys
from collections import Counter, OrderedDict

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Terminal-branch order used everywhere else in this project.
SPECIES_ORDER = ['Ascaris', 'Parascaris', 'Baylisascaris', 'Toxocara']

# Default topology. Matches four_way_kaks.py: the internal split is
# {Ascaris, Parascaris} | {Baylisascaris, Toxocara}. Written rooted, which
# clock = 1 requires; baseml unroots it internally for clock = 0.
DEFAULT_TOPOLOGY = '((Ascaris,Parascaris),(Baylisascaris,Toxocara));'

# Sister pairs implied by DEFAULT_TOPOLOGY. Tajima's test requires the outgroup
# to be equally distant from BOTH test taxa, so only sister pairs are valid --
# comparing Ascaris with Baylisascaris using Toxocara as outgroup would be
# wrong, because Toxocara is closer to Baylisascaris than to Ascaris.
DEFAULT_SISTERS = [('Ascaris', 'Parascaris'), ('Baylisascaris', 'Toxocara')]

# Header patterns observed in this project's FASTA files. Ascaris is
# deliberately absent: its identifier format has not been confirmed, so it is
# resolved by elimination (see assign_species) rather than guessed.
DEFAULT_PATTERNS = OrderedDict([
    ('Toxocara',      r'^MSTRG\.'),        # StringTie models + TransDecoder ORFs
    ('Parascaris',    r'^MP\d'),           # miniprot Pu models
    ('Baylisascaris', r'^[gG]\d+'),        # Augustus; note the g/G case bug
])

DNA = set('ACGT')
STOPS = {'TAA', 'TAG', 'TGA'}

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# FASTA / alignment handling
# ---------------------------------------------------------------------------

def read_fasta(path):
    """
    Minimal FASTA reader.

    Biopython is avoided here so the script has no dependency beyond numpy;
    these files are simple aligned FASTA and need no format tolerance.

    Returns
    -------
    OrderedDict of {header: sequence}, sequence uppercased with whitespace
    stripped. Header is truncated at the first whitespace, matching how PAML
    and most aligners treat it.
    """
    seqs = OrderedDict()
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks).upper()
                name = line[1:].split()[0] if len(line) > 1 else ''
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(chunks).upper()
    return seqs


def discover_alignments(aln_dir, pattern):
    """Sorted list of alignment files under aln_dir matching pattern."""
    if not os.path.isdir(aln_dir):
        sys.exit(f"ERROR: --aln-dir not found: {aln_dir}")
    files = sorted(glob.glob(os.path.join(aln_dir, pattern)))
    if not files:
        # Give the user something actionable rather than an empty result.
        present = sorted({os.path.splitext(f)[1]
                          for f in os.listdir(aln_dir)})[:10]
        sys.exit(f"ERROR: no files match {pattern!r} in {aln_dir}. "
                 f"Extensions present: {present or '(directory is empty)'}. "
                 f"Set --aln-glob accordingly.")
    return files


def compile_patterns(pattern_args):
    """
    Build {species: compiled regex} from --species-pattern Name=regex.

    User-supplied patterns override the defaults for the same species name.
    """
    pats = OrderedDict((k, re.compile(v)) for k, v in DEFAULT_PATTERNS.items())
    for spec in pattern_args or []:
        if '=' not in spec:
            sys.exit(f"ERROR: --species-pattern needs Name=regex, got {spec!r}")
        name, rx = spec.split('=', 1)
        name = name.strip()
        if name not in SPECIES_ORDER:
            sys.exit(f"ERROR: unknown species {name!r} in --species-pattern; "
                     f"expected one of {SPECIES_ORDER}")
        try:
            pats[name] = re.compile(rx)
        except re.error as exc:
            sys.exit(f"ERROR: bad regex for {name}: {exc}")
    return pats


def assign_species(headers, patterns):
    """
    Map each FASTA header in one alignment to a species.

    Resolution order, most reliable first:

      1. The header IS a species name exactly.
      2. A delimiter-separated token of the header is a species name, e.g.
         'Ascaris|AgB00001' or 'Baylisascaris_g9000'.
      3. The header contains a species name as a bare substring.
      4. The header matches that species' identifier regex.
      5. Exactly one species and exactly one header remain unmatched, so the
         leftover header must belong to the leftover species. This is what
         resolves Ascaris without hard-coding a guess at its ID format.

    NAME NESTING: 'Ascaris' is a substring of BOTH 'Parascaris' and
    'Baylisascaris'. A naive substring scan in SPECIES_ORDER order will happily
    assign 'Baylisascaris|g9000' to Ascaris and then fail, or worse succeed and
    silently swap two species. Rules 2 and 3 therefore walk species in
    DESCENDING NAME LENGTH, so Baylisascaris and Parascaris claim their headers
    before Ascaris is ever considered. Do not "simplify" this back to
    SPECIES_ORDER.

    Returns
    -------
    dict of {species: header}, containing only species that could be resolved.
    """
    assigned = {}
    used = set()

    # Longest first: Baylisascaris, Parascaris, Toxocara, Ascaris.
    by_length = sorted(SPECIES_ORDER, key=len, reverse=True)
    token_split = re.compile(r'[|_:;,/\s.\-]+')

    for sp in SPECIES_ORDER:                                   # rule 1
        for h in headers:
            if h not in used and h == sp:
                assigned[sp] = h
                used.add(h)
                break

    for sp in by_length:                                       # rule 2
        if sp in assigned:
            continue
        for h in headers:
            if h in used:
                continue
            if sp.lower() in [t.lower() for t in token_split.split(h) if t]:
                assigned[sp] = h
                used.add(h)
                break

    for sp in by_length:                                       # rule 3
        if sp in assigned:
            continue
        for h in headers:
            if h not in used and sp.lower() in h.lower():
                assigned[sp] = h
                used.add(h)
                break

    for sp, rx in patterns.items():                            # rule 3
        if sp in assigned:
            continue
        for h in headers:
            if h not in used and rx.search(h):
                assigned[sp] = h
                used.add(h)
                break

    missing_sp = [s for s in SPECIES_ORDER if s not in assigned]
    free_h = [h for h in headers if h not in used]
    if len(missing_sp) == 1 and len(free_h) == 1:              # rule 4
        assigned[missing_sp[0]] = free_h[0]

    return assigned


def inspect(files, patterns, n_show=5):
    """
    Report what the headers look like and how well assignment works, then exit.

    Run this before anything else. Silent mis-assignment of a species would
    invalidate every number the script produces, so the mapping is made visible
    rather than assumed.
    """
    logger.info(f"Inspecting {len(files)} alignment files")
    logger.info(f"\nFirst {n_show} files, headers as found:")
    for path in files[:n_show]:
        seqs = read_fasta(path)
        logger.info(f"  {os.path.basename(path)}: {list(seqs.keys())}")

    resolved = Counter()
    complete = 0
    sample = files[:2000]
    for path in sample:
        got = assign_species(list(read_fasta(path).keys()), patterns)
        for sp in got:
            resolved[sp] += 1
        if len(got) == len(SPECIES_ORDER):
            complete += 1

    logger.info(f"\nAssignment over {len(sample)} files:")
    for sp in SPECIES_ORDER:
        n = resolved[sp]
        logger.info(f"  {sp:<15} {n:>6} / {len(sample)} "
                    f"({100 * n / len(sample):5.1f}%)")
    logger.info(f"  all four present: {complete} / {len(sample)} "
                f"({100 * complete / len(sample):5.1f}%)")

    if complete < len(sample):
        logger.warning(
            "\nSome files did not resolve to all four species. If a whole "
            "species is at 0%, pass its pattern explicitly, e.g.\n"
            "  --species-pattern 'Ascaris=^AgB'\n"
            "Files missing any species are skipped by the real run.")
    else:
        logger.info("\nAll files resolved. Rerun without --inspect.")


# ---------------------------------------------------------------------------
# Supermatrix
# ---------------------------------------------------------------------------

def subset_sites(seq, sites):
    """
    Keep the requested codon positions.

    'all' returns the sequence unchanged; '12' drops every third base; '3'
    keeps only every third base. The sequence must already be a codon
    alignment, which is what branch_dnds_four.py produces.
    """
    if sites == 'all':
        return seq
    if sites == '12':
        return ''.join(seq[i] for i in range(len(seq)) if i % 3 != 2)
    if sites == '3':
        return seq[2::3]
    raise ValueError(f"unknown --sites {sites!r}")


def build_supermatrix(files, patterns, sites, max_genes=None,
                      drop_stop_codons=True):
    """
    Concatenate per-gene alignments into one alignment per species.

    A gene is used only if all four species are present AND all four sequences
    have the same length AND that length is a multiple of three. Anything else
    is counted and skipped -- a ragged gene would shift every downstream codon
    position and silently corrupt the --sites subsetting.

    Returns
    -------
    (seqs, stats) where seqs is {species: concatenated string} and stats is a
    dict of counts explaining what was dropped and why.
    """
    parts = {sp: [] for sp in SPECIES_ORDER}
    stats = Counter()
    used_lengths = []

    for path in files:
        if max_genes and stats['used'] >= max_genes:
            break
        stats['seen'] += 1
        seqs = read_fasta(path)
        if not seqs:
            stats['empty'] += 1
            continue

        got = assign_species(list(seqs.keys()), patterns)
        if len(got) != len(SPECIES_ORDER):
            stats['incomplete_taxa'] += 1
            continue

        block = {sp: seqs[got[sp]] for sp in SPECIES_ORDER}
        lengths = {len(v) for v in block.values()}
        if len(lengths) != 1:
            stats['ragged'] += 1
            continue
        L = lengths.pop()
        if L == 0:
            stats['empty'] += 1
            continue
        if L % 3 != 0:
            stats['not_codon_multiple'] += 1
            continue

        # Trailing stop codons are not homologous positions across species and
        # PAML will complain about them in a codon context. Drop the last codon
        # only if it is a stop in at least one species.
        if drop_stop_codons and any(v[-3:] in STOPS for v in block.values()):
            block = {sp: v[:-3] for sp, v in block.items()}
            L -= 3
            stats['stop_trimmed'] += 1
            if L == 0:
                stats['empty'] += 1
                continue

        for sp in SPECIES_ORDER:
            parts[sp].append(subset_sites(block[sp], sites))
        used_lengths.append(len(subset_sites(block[SPECIES_ORDER[0]], sites)))
        stats['used'] += 1

    if not stats['used']:
        sys.exit("ERROR: no usable alignments. Run --inspect to check that "
                 "headers map to all four species.")

    out = {sp: ''.join(parts[sp]) for sp in SPECIES_ORDER}
    final = {len(v) for v in out.values()}
    assert len(final) == 1, 'supermatrix rows differ in length'
    stats['sites'] = final.pop()
    stats['mean_gene_len'] = int(np.mean(used_lengths)) if used_lengths else 0
    return out, stats


def write_phylip(path, seqs):
    """
    PHYLIP sequential, the format baseml expects.

    Names are species names, padded to a fixed width. PAML truncates long names
    and chokes on some punctuation, so using clean species labels here avoids a
    class of silent failures that gene IDs would invite.
    """
    n = len(seqs)
    L = len(next(iter(seqs.values())))
    with open(path, 'w') as fh:
        fh.write(f"  {n}  {L}\n")
        for name, seq in seqs.items():
            fh.write(f"{name[:30]:<32}{seq}\n")
    return L


def top_level_children(nwk):
    """
    Split the outermost Newick clade into its immediate children.

    Returns a list of substrings, or None if `nwk` is not a parenthesised
    clade. Commas inside nested clades are ignored.
    """
    s = nwk.strip().rstrip(';').strip()
    if not (s.startswith('(') and s.endswith(')')):
        return None
    inner = s[1:-1]
    kids, depth, start = [], 0, 0
    for i, ch in enumerate(inner):
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        elif ch == ',' and depth == 0:
            kids.append(inner[start:i].strip())
            start = i + 1
    kids.append(inner[start:].strip())
    return kids


def is_rooted(nwk):
    """A rooted bifurcating Newick has exactly two top-level children."""
    kids = top_level_children(nwk)
    return kids is not None and len(kids) == 2


def unroot_newick(nwk):
    """
    Convert a rooted bifurcating topology to its unrooted (trifurcating) form.

    WHY THIS MATTERS: given a rooted tree with clock = 0, baseml estimates a
    length for BOTH root branches even though only their sum is identifiable.
    That leaves the free model over-parameterised by one, and the LRT then
    takes df = 3 instead of the correct n - 2 = 2 for four taxa. Feeding the
    unrooted form to the clock = 0 fit and the rooted form to the clock = 1 fit
    is what makes the degrees of freedom come out right.

        ((A,B),(C,D));  ->  (A,B,(C,D));

    Returns the input unchanged if it is already unrooted or cannot be parsed.
    """
    kids = top_level_children(nwk)
    if kids is None or len(kids) != 2:
        return nwk.strip()
    # Splice whichever child is itself a clade up into the root.
    for i, j in ((0, 1), (1, 0)):
        k = kids[i]
        if k.startswith('(') and k.endswith(')'):
            return f"({k[1:-1]},{kids[j]});"
    return nwk.strip()                       # both children are leaves


def write_tree(path, topology, ntaxa):
    """
    PAML tree file: a "<ntaxa> <ntrees>" header line, then the Newick string.

    Omitting that header is the single most common PAML input error in this
    project's notes.
    """
    with open(path, 'w') as fh:
        fh.write(f"{ntaxa} 1\n")
        fh.write(topology.strip().rstrip(';') + ';\n')


# ---------------------------------------------------------------------------
# PAML
# ---------------------------------------------------------------------------

BASEML_CTL = """\
      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}

        noisy = 3
      verbose = 0
      runmode = 0

        model = {model}
        Mgene = 0

        clock = {clock}
    fix_kappa = 0
        kappa = 5
    fix_alpha = 0
        alpha = 0.5
       Malpha = 0
        ncatG = {ncatG}

        nhomo = 0
        getSE = 0
 RateAncestor = 0
   Small_Diff = 7e-6
    cleandata = {cleandata}
       method = 0
"""

# lnL(ntime:  5  np:  8):  -123456.789012      +0.000000
LNL_RE = re.compile(r'lnL\s*\(\s*ntime:\s*(\d+)\s+np:\s*(\d+)\s*\)\s*:\s*'
                    r'(-?\d+\.?\d*)')


def run_baseml(workdir, phy, tree, clock, model, ncatG, cleandata, exe):
    """
    Run one baseml fit in its own directory and parse the result.

    baseml writes rst, rub, lnf and mlb into its working directory with fixed
    names, so two fits in one directory would overwrite each other's
    intermediates. Each fit therefore gets a private subdirectory.

    IMPORTANT: baseml returns a non-zero exit code on perfectly successful
    runs. Success is judged solely by whether the output file exists and
    contains a parsable lnL line.

    Returns
    -------
    dict with lnL (float), np (int, number of free parameters), ntime (int),
    and outfile (path).
    """
    os.makedirs(workdir, exist_ok=True)
    ctl_path = os.path.join(workdir, 'baseml.ctl')
    outfile = 'mlb.out'

    # Relative paths: baseml is invoked with cwd=workdir.
    with open(ctl_path, 'w') as fh:
        fh.write(BASEML_CTL.format(
            seqfile=os.path.relpath(phy, workdir),
            treefile=os.path.relpath(tree, workdir),
            outfile=outfile, clock=clock, model=model, ncatG=ncatG,
            cleandata=cleandata))

    logger.info(f"  baseml clock={clock} in {workdir}/ ...")
    proc = subprocess.run([exe, 'baseml.ctl'], cwd=workdir,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True)

    out_path = os.path.join(workdir, outfile)
    if not os.path.isfile(out_path):
        tail = '\n'.join(proc.stdout.splitlines()[-25:])
        sys.exit(f"ERROR: baseml produced no output for clock={clock}.\n"
                 f"(exit code {proc.returncode} is NOT diagnostic for PAML)\n"
                 f"--- last lines of baseml output ---\n{tail}")

    with open(out_path) as fh:
        text = fh.read()
    m = LNL_RE.search(text)
    if not m:
        sys.exit(f"ERROR: no lnL line found in {out_path}. baseml may have "
                 f"failed after creating the file; inspect it directly.")

    ntime, npar, lnl = int(m.group(1)), int(m.group(2)), float(m.group(3))
    logger.info(f"    lnL = {lnl:.4f}   np = {npar}   ntime = {ntime}")
    return {'lnL': lnl, 'np': npar, 'ntime': ntime, 'outfile': out_path}


def chi2_sf(x, df):
    """
    Upper-tail probability of a chi-square variate, without scipy.

    scipy is available in the kaks environment, but keeping this dependency-free
    means the script also runs anywhere Python and numpy exist. Uses the
    regularised upper incomplete gamma via a continued fraction (Lentz), with a
    series expansion in the region where the fraction converges slowly.
    """
    if not np.isfinite(x) or x < 0:
        return float('nan')
    if x == 0:
        return 1.0
    a, xx = df / 2.0, x / 2.0

    # Log gamma via Lanczos, adequate for the precision needed here.
    def lgamma(z):
        g = [676.5203681218851, -1259.1392167224028, 771.32342877765313,
             -176.61502916214059, 12.507343278686905, -0.13857109526572012,
             9.9843695780195716e-6, 1.5056327351493116e-7]
        if z < 0.5:
            return (np.log(np.pi / np.sin(np.pi * z))
                    - lgamma(1 - z))
        z -= 1
        a0 = 0.99999999999980993
        t = z + 7.5
        for i, c in enumerate(g):
            a0 += c / (z + i + 1)
        return 0.5 * np.log(2 * np.pi) + (z + 0.5) * np.log(t) - t + np.log(a0)

    lg = lgamma(a)
    if xx < a + 1:                                   # series for P(a, x)
        term = 1.0 / a
        total = term
        n = 1
        while n < 10000:
            term *= xx / (a + n)
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
            n += 1
        p_lower = total * np.exp(-xx + a * np.log(xx) - lg)
        return max(0.0, min(1.0, 1.0 - p_lower))

    # continued fraction for Q(a, x)
    tiny = 1e-300
    b = xx + 1 - a
    c = 1 / tiny
    d = 1 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-15:
            break
    q = np.exp(-xx + a * np.log(xx) - lg) * h
    return max(0.0, min(1.0, q))


# ---------------------------------------------------------------------------
# Tajima's relative-rate test
# ---------------------------------------------------------------------------

def tajima_rrt(s1, s2, out):
    """
    Tajima's 1D relative-rate test (Tajima 1993, Genetics 135:599-607).

    Consider only sites where all three sequences carry an unambiguous A/C/G/T.
    Let

        n1 = sites where s2 == outgroup and s1 differs
             -> a substitution occurred on lineage 1
        n2 = sites where s1 == outgroup and s2 differs
             -> a substitution occurred on lineage 2

    Under equal rates in the two lineages since their common ancestor,
    E(n1) = E(n2), and

        chi2 = (n1 - n2)^2 / (n1 + n2),  df = 1.

    Sites where all three differ are uninformative for this contrast and are
    excluded, as are sites where all three agree.

    The test is model-free: it never estimates a branch length or assumes a
    substitution model, so it cannot be dismissed as model misspecification.
    Its assumption is only that the outgroup is genuinely outside the pair --
    which is why the caller must pass sister taxa.

    Returns
    -------
    dict with n1, n2, n_informative, n_compared, chi2, p_value, and the
    lineage-1/lineage-2 rate ratio (n1/n2).
    """
    if not (len(s1) == len(s2) == len(out)):
        raise ValueError('sequences differ in length')

    a1 = np.frombuffer(s1.encode(), dtype='S1')
    a2 = np.frombuffer(s2.encode(), dtype='S1')
    ao = np.frombuffer(out.encode(), dtype='S1')

    valid = np.ones(len(a1), dtype=bool)
    for arr in (a1, a2, ao):
        ok = np.zeros(len(arr), dtype=bool)
        for base in DNA:
            ok |= (arr == base.encode())
        valid &= ok

    v1, v2, vo = a1[valid], a2[valid], ao[valid]
    n1 = int(np.sum((v2 == vo) & (v1 != vo)))
    n2 = int(np.sum((v1 == vo) & (v2 != vo)))

    tot = n1 + n2
    chi2 = ((n1 - n2) ** 2) / tot if tot else float('nan')
    p = chi2_sf(chi2, 1) if tot else float('nan')
    ratio = (n1 / n2) if n2 else float('inf') if n1 else float('nan')

    return {'n1': n1, 'n2': n2, 'n_informative': tot,
            'n_compared': int(valid.sum()), 'chi2': chi2,
            'p_value': p, 'rate_ratio': ratio}


def bh_qvalues(pvals):
    """Benjamini-Hochberg adjusted p-values; NaN passes through as NaN."""
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    adj = np.empty(n)
    adj[order] = ranked
    q[ok] = adj
    return q


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def stars(p):
    """Conventional significance markers."""
    if p is None or not np.isfinite(p):
        return 'n/a'
    for thr, s in ((1e-4, '****'), (1e-3, '***'), (1e-2, '**'), (0.05, '*')):
        if p < thr:
            return s
    return 'ns'


def fmt_p(p):
    """
    Format a p-value for display.

    Chi-square tail probabilities underflow to exactly 0.0 in double precision
    once the statistic passes roughly 1400. Printing '0' invites a reviewer to
    read it as a missing value, so an exact zero is reported as a bound.
    """
    if p is None or not np.isfinite(p):
        return 'n/a'
    if p == 0.0:
        return '<1e-300'
    return f'{p:.3g}'


def write_outputs(outdir, args, stats, lrt, rrt_rows, phy_len):
    """Write the TSV and the human-readable summary."""
    tsv = os.path.join(outdir, 'clock_test_results.tsv')
    with open(tsv, 'w') as fh:
        fh.write('test\tcomparison\tstatistic\tvalue\tdf\tp_value\tq_value\n')
        fh.write(f"global_clock_LRT\twhole_tree\t2*dlnL\t{lrt['stat']:.6f}\t"
                 f"{lrt['df']}\t{lrt['p_value']:.6g}\tNA\n")
        for r in rrt_rows:
            fh.write(f"tajima_RRT\t{r['label']}\tchi2\t{r['chi2']:.6f}\t1\t"
                     f"{r['p_value']:.6g}\t{r['q_value']:.6g}\n")
    logger.info(f"Results: {tsv}")

    path = os.path.join(outdir, 'clock_test_summary.txt')
    with open(path, 'w') as fh:
        fh.write('Molecular-clock tests, four-way ascarid quartets\n')
        fh.write('=' * 62 + '\n\n')
        fh.write(f"Alignments:      {args.aln_dir}/{args.aln_glob}\n")
        fh.write(f"Codon positions: {args.sites}\n")
        fh.write(f"Topology:        {args.topology}\n")
        fh.write(f"Genes used:      {stats['used']} of {stats['seen']} "
                 f"files seen\n")
        fh.write(f"Supermatrix:     {phy_len:,} sites x "
                 f"{len(SPECIES_ORDER)} taxa\n")
        skipped = [(k, stats[k]) for k in
                   ('incomplete_taxa', 'ragged', 'not_codon_multiple', 'empty')
                   if stats[k]]
        if skipped:
            fh.write("Skipped:         "
                     + ', '.join(f"{k}={v}" for k, v in skipped) + '\n')
        fh.write('\n')

        fh.write('-' * 62 + '\n')
        fh.write('TEST 1  Global-clock likelihood ratio test (baseml)\n')
        fh.write('-' * 62 + '\n')
        fh.write(f"  clock = 0 (free branches)   lnL = {lrt['lnL0']:15.4f}   "
                 f"np = {lrt['np0']}\n")
        fh.write(f"  clock = 1 (global clock)    lnL = {lrt['lnL1']:15.4f}   "
                 f"np = {lrt['np1']}\n\n")
        fh.write(f"  2*dlnL = {lrt['stat']:.4f}   df = {lrt['df']}   "
                 f"P = {fmt_p(lrt['p_value'])}   {stars(lrt['p_value'])}\n")
        if lrt['df'] != lrt['df_theory']:
            fh.write(f"  NOTE: PAML's parameter counts give df = {lrt['df']}, "
                     f"but n-2 = {lrt['df_theory']}.\n        The PAML count is "
                     f"used. Check the tree if these disagree.\n")
        fh.write('\n')
        if lrt['p_value'] < 0.05:
            fh.write("  => The strict molecular clock is REJECTED. Substitution\n"
                     "     rates differ among these lineages, so a single rate\n"
                     "     cannot be used to convert Ks into divergence times.\n")
        else:
            fh.write("  => The clock is not rejected at P < 0.05.\n")
        fh.write('\n')

        fh.write('-' * 62 + '\n')
        fh.write("TEST 2  Tajima's 1D relative-rate test (model-free)\n")
        fh.write('-' * 62 + '\n')
        width = max(34, max(len(r['label']) for r in rrt_rows) + 2)
        fh.write(f"  {'comparison':<{width}}{'n1':>8}{'n2':>8}{'ratio':>8}"
                 f"{'P':>10}{'q':>10}\n")
        for r in rrt_rows:
            fh.write(f"  {r['label']:<{width}}{r['n1']:>8}{r['n2']:>8}"
                     f"{r['rate_ratio']:>8.2f}{fmt_p(r['p_value']):>10}"
                     f"{fmt_p(r['q_value']):>10}  {stars(r['q_value'])}\n")
        fh.write('\n  n1 = substitutions on the first lineage, n2 on the '
                 'second,\n  counted against the outgroup. Equal under a '
                 'clock.\n  q = Benjamini-Hochberg across the tests shown.\n\n')

        sig = [r for r in rrt_rows if np.isfinite(r['q_value'])
               and r['q_value'] < 0.05]
        if sig:
            fh.write("  => Rate equality is rejected for: "
                     + ', '.join(r['label'] for r in sig) + '\n')
            fh.write("     This is independent of any substitution model.\n")
        else:
            fh.write("  => No pair rejects rate equality.\n")
        fh.write('\n')

        fh.write('-' * 62 + '\n')
        fh.write('INTERPRETATION\n')
        fh.write('-' * 62 + '\n')
        fh.write(
            "Both tests address one question: do these lineages accumulate\n"
            "substitutions at a common rate? If they do not, then pairwise Ks\n"
            "measures path length rather than elapsed time, the same node\n"
            "dates differently depending on which pair is used, and reporting\n"
            "divergence times from a single substitution rate is not\n"
            "defensible. Report median Ks and branch lengths instead, and\n"
            "leave absolute dating to a relaxed-clock analysis with more taxa\n"
            "and more than one calibration.\n")
    logger.info(f"Summary: {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Global-clock LRT and Tajima relative-rate tests for the '
                    'four-way ascarid quartets.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Run --inspect first to confirm species assignment.')

    p.add_argument('--aln-dir', required=True,
                   help='Directory of per-gene codon alignments (FASTA), e.g. '
                        'four_way_alignments')
    p.add_argument('--aln-glob', default='*.fa',
                   help="Filename pattern inside --aln-dir (default: '*.fa'; "
                        "try '*.fasta' or '*.aln' if nothing matches)")
    p.add_argument('--outdir', default='clock_test')

    p.add_argument('--sites', default='12', choices=['all', '12', '3'],
                   help="Codon positions to use. '12' (default) avoids the "
                        "saturated third positions, so a clock rejection "
                        "cannot be blamed on Toxocara saturation. '3' isolates "
                        "the synonymous signal. 'all' uses everything.")
    p.add_argument('--topology', default=DEFAULT_TOPOLOGY,
                   help='Rooted Newick topology with species names as tips. '
                        'clock=1 requires a rooted tree. Default: '
                        + DEFAULT_TOPOLOGY)
    p.add_argument('--sisters', default=None,
                   help="Sister pairs for Tajima's test, as "
                        "'A+B,C+D'. Default matches the topology: "
                        "'Ascaris+Parascaris,Baylisascaris+Toxocara'. Only "
                        "sister pairs are valid -- the outgroup must be "
                        "equidistant from both test taxa.")

    p.add_argument('--species-pattern', action='append', default=None,
                   metavar='NAME=REGEX',
                   help='Override the header regex for one species, e.g. '
                        "--species-pattern 'Ascaris=^AgB'. Repeatable.")
    p.add_argument('--inspect', action='store_true',
                   help='Report header formats and species-assignment rates, '
                        'then exit without running anything.')

    p.add_argument('--model', type=int, default=7,
                   help='baseml nucleotide model (default: 7 = GTR/REV)')
    p.add_argument('--ncatG', type=int, default=5,
                   help='Discrete gamma rate categories (default: 5)')
    p.add_argument('--cleandata', type=int, default=0, choices=[0, 1],
                   help='PAML cleandata: 0 keeps sites with gaps/ambiguity '
                        '(default), 1 removes them')
    p.add_argument('--baseml', default='baseml',
                   help='Path to the baseml executable (default: baseml)')

    p.add_argument('--max-genes', type=int, default=None,
                   help='Use at most this many genes. For a fast smoke test; '
                        'omit for the real run.')
    p.add_argument('--keep-stop-codons', action='store_true',
                   help='Do not trim a trailing stop codon')
    return p.parse_args()


def parse_sisters(arg):
    """Parse 'A+B,C+D' into [(A, B), (C, D)] with validation."""
    if not arg:
        return DEFAULT_SISTERS
    pairs = []
    for chunk in arg.split(','):
        bits = [b.strip() for b in chunk.split('+')]
        if len(bits) != 2:
            sys.exit(f"ERROR: --sisters needs 'A+B' per pair, got {chunk!r}")
        for b in bits:
            if b not in SPECIES_ORDER:
                sys.exit(f"ERROR: unknown species {b!r} in --sisters")
        pairs.append((bits[0], bits[1]))
    return pairs


def main():
    args = parse_args()
    patterns = compile_patterns(args.species_pattern)
    files = discover_alignments(args.aln_dir, args.aln_glob)

    if args.inspect:
        inspect(files, patterns)
        return

    os.makedirs(args.outdir, exist_ok=True)
    logging.getLogger().addHandler(
        logging.FileHandler(os.path.join(args.outdir, 'clock_test.log')))

    sisters = parse_sisters(args.sisters)
    logger.info(f"Found {len(files)} alignment files in {args.aln_dir}")
    logger.info(f"Codon positions: {args.sites}")

    # --- supermatrix ------------------------------------------------------
    seqs, stats = build_supermatrix(files, patterns, args.sites,
                                    args.max_genes,
                                    drop_stop_codons=not args.keep_stop_codons)
    logger.info(f"Concatenated {stats['used']} genes -> {stats['sites']:,} "
                f"sites x {len(SPECIES_ORDER)} taxa")
    for key in ('incomplete_taxa', 'ragged', 'not_codon_multiple', 'empty',
                'stop_trimmed'):
        if stats[key]:
            logger.info(f"  {key}: {stats[key]}")

    phy = os.path.join(args.outdir, 'supermatrix.phy')
    phy_len = write_phylip(phy, seqs)

    # Two tree files, not one. clock = 1 needs the rooted topology; clock = 0
    # needs the unrooted one, or baseml estimates an unidentifiable extra root
    # branch and the LRT df comes out one too large. See unroot_newick().
    if not is_rooted(args.topology):
        sys.exit(f"ERROR: --topology must be rooted (two top-level clades) "
                 f"because clock=1 requires a rooted tree. Got: "
                 f"{args.topology}")
    tree_rooted = os.path.join(args.outdir, 'tree_rooted.trees')
    tree_unrooted = os.path.join(args.outdir, 'tree_unrooted.trees')
    write_tree(tree_rooted, args.topology, len(SPECIES_ORDER))
    write_tree(tree_unrooted, unroot_newick(args.topology),
               len(SPECIES_ORDER))
    logger.info(f"Wrote {phy}")
    logger.info(f"  rooted   (clock=1): {args.topology.strip()}")
    logger.info(f"  unrooted (clock=0): {unroot_newick(args.topology)}")

    # --- TEST 1: global-clock LRT ----------------------------------------
    logger.info("\nTEST 1: global-clock likelihood ratio test")
    fit0 = run_baseml(os.path.join(args.outdir, 'clock0'), phy, tree_unrooted,
                      0, args.model, args.ncatG, args.cleandata, args.baseml)
    fit1 = run_baseml(os.path.join(args.outdir, 'clock1'), phy, tree_rooted,
                      1, args.model, args.ncatG, args.cleandata, args.baseml)

    # df from PAML's own parameter counts rather than an assumed n-2: if the
    # topology or model is changed, the assumed value would silently go wrong.
    df = fit0['np'] - fit1['np']
    df_theory = len(SPECIES_ORDER) - 2
    if df <= 0:
        sys.exit(f"ERROR: non-positive df ({fit0['np']} - {fit1['np']}). "
                 f"The clock model did not constrain anything; check that "
                 f"--topology is rooted.")
    stat = 2.0 * (fit0['lnL'] - fit1['lnL'])
    if stat < 0:
        logger.warning(f"2*dlnL is negative ({stat:.4f}). The constrained fit "
                       f"beat the free one, which means an optimisation "
                       f"failure, not evidence for the clock. Treating as 0.")
        stat = 0.0
    p_lrt = chi2_sf(stat, df)
    logger.info(f"  2*dlnL = {stat:.4f}, df = {df}, P = {fmt_p(p_lrt)} "
                f"{stars(p_lrt)}")
    if df != df_theory:
        logger.warning(f"  df from PAML ({df}) != n-2 ({df_theory})")

    lrt = {'lnL0': fit0['lnL'], 'lnL1': fit1['lnL'], 'np0': fit0['np'],
           'np1': fit1['np'], 'stat': stat, 'df': df, 'df_theory': df_theory,
           'p_value': p_lrt}

    # --- TEST 2: Tajima relative-rate ------------------------------------
    logger.info("\nTEST 2: Tajima relative-rate tests")
    rrt_rows = []
    for a, b in sisters:
        outgroups = [s for s in SPECIES_ORDER if s not in (a, b)]
        for og in outgroups:
            r = tajima_rrt(seqs[a], seqs[b], seqs[og])
            r['label'] = f'{a} vs {b} (out: {og})'
            rrt_rows.append(r)
            logger.info(f"  {r['label']:<44} n1={r['n1']:>7} "
                        f"n2={r['n2']:>7} ratio={r['rate_ratio']:.2f} "
                        f"P={fmt_p(r['p_value'])}")

    qs = bh_qvalues([r['p_value'] for r in rrt_rows])
    for r, q in zip(rrt_rows, qs):
        r['q_value'] = float(q)

    write_outputs(args.outdir, args, stats, lrt, rrt_rows, phy_len)
    logger.info('Done.')


if __name__ == '__main__':
    main()
