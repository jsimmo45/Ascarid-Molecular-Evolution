#!/bin/bash
#SBATCH -J s1_composition
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=s1_composition_%j.out
#SBATCH --error=s1_composition_%j.err
#
# Supplementary Figure S2 (file figS1_composition, named before the
# renumbering): GC3s and coding GC of eliminated vs retained Ascaris genes,
# dN/dS against GC3s, and a GC3s-matched dN/dS contrast that tests the Results
# claim that composition biases the main contrast downward.
#
# 2026-10-19 (results draft note 67): the C key shows a blue and a red open
# triangle for genes above the axis, as drawn, and every P is plain text with
# a raised exponent (lib/plain_text.py). Statistics unchanged.
#
#   sbatch run_compositional_analysis.sh
#   sbatch run_compositional_analysis.sh <path/to/kaks_annotated.tsv>
#
# Inputs (config.sh; relative paths are taken from PROJECT_DIR):
#   ANNOTATED_TABLE        argument 1 overrides it
#   AS_CDS                 default As_cds.fa
#   CODON_PER_GENE_TABLE   optional; default
#                          ${RESULTS_DIR}/codon_usage_ascaris_20260817_v3/codon_usage_per_gene.tsv
#                          (the Figure 2 run). Used only for a per-gene check
#                          that S1 and Figure 2 compute the same GC3s; skipped
#                          with a note if the file is not there.
#
# Replaces the May root script of the same name, which counted ATG and TGG in
# GC3 and drew the superseded gene-density profile. The old copy still sits in
# the miniprot directory, so the wrapper refuses any compositional_analysis.py
# that does not import paper_colors.
#
# Single-threaded: ~9,400 sequences and 500 matched replicates, about a
# minute. 1 node / 1 task / 1 CPU.

set -euo pipefail

# Under sbatch, SLURM copies the batch script to a spool directory on the compute
# node, so BASH_SOURCE points at /var/spool/slurm/... and neither config.sh nor
# the .py file is anywhere near it. SLURM_SUBMIT_DIR is the directory sbatch was
# run from, which is on the shared filesystem, so search from there as well.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

search_up() {                       # search_up <start dir> <filename>
  local p="$1"
  [[ -n "$p" && -d "$p" ]] || return 1
  p="$(cd "$p" && pwd)"
  for _ in 1 2 3 4 5 6; do
    [[ -f "${p}/$2" ]] && { printf '%s' "${p}/$2"; return 0; }
    [[ "${p}" == "/" ]] && break
    p="$(dirname "${p}")"
  done
  return 1
}

# The repo may sit *below* the directory sbatch was run from, so probe one level
# down as well, accepting only a directory that also looks like the repo.
search_down() {                     # search_down <start dir>
  local c
  for c in "$1"/*/config.sh; do
    [[ -f "$c" ]] || continue
    local r; r="$(dirname "$c")"
    [[ -d "${r}/05_analyses" || -f "${r}/config.example.sh" ]] || continue
    printf '%s' "$c"; return 0
  done
  return 1
}

CONFIG="${KAKS_CONFIG:-}"
if [[ -z "${CONFIG}" ]]; then
  for d in "${HERE}" "${SLURM_SUBMIT_DIR:-}" "${PWD}"; do
    CONFIG="$(search_up "$d" config.sh || true)"
    [[ -n "${CONFIG}" ]] && break
    CONFIG="$(search_down "$d" || true)"
    [[ -n "${CONFIG}" ]] && break
  done
fi
if [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then
  echo "ERROR: no config.sh found." >&2
  echo "  searched upward from: ${HERE}" >&2
  echo "                        ${SLURM_SUBMIT_DIR:-<not under sbatch>}" >&2
  echo "                        ${PWD}" >&2
  echo "  Set KAKS_CONFIG=/full/path/to/config.sh and resubmit." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" \
           "${REPO_ROOT}/06_figures" "${REPO_ROOT}" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script compositional_analysis.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: compositional_analysis.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# The May script of the same name is still in the miniprot directory. Only the
# converted one imports the shared figure style, so refuse anything that does not.
grep -q '^from paper_colors import' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the old May script (no paper_colors import)." >&2
  echo "  Use Ascarid-Molecular-Evolution/05_analyses/compositional_analysis.py." >&2
  exit 1; }
grep -q '^from plain_text import draw_p' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the version from before 2026-10-19 (black key triangle; mathtext P)." >&2
  exit 1; }

PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHONPATH

# The kaks environment, with its bin/ first on PATH: on Isaac compute nodes the
# Spack python otherwise shadows the conda one. config.sh can override all three.
CONDA_SH="${CONDA_SH:-}"
CONDA_ENV="${CONDA_ENV:-kaks}"
CONDA_BIN_PREPEND="${CONDA_BIN_PREPEND:-}"
if [[ -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
fi
PY="$(command -v python || command -v python3)"

# Probe the imports: paper_colors, region_lib, rate_tables and plain_text
# live in lib/.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("paper_colors", "region_lib", "rate_tables", "plain_text",
                       "scipy", "matplotlib")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, plain_text, scipy, matplotlib all importable")
PYEOF

cd "${PROJECT_DIR}"

TABLE="${1:-${ANNOTATED_TABLE:-}}"
AS_CDS="${AS_CDS:-As_cds.fa}"
CODON_PER_GENE_TABLE="${CODON_PER_GENE_TABLE:-${RESULTS_DIR}/codon_usage_ascaris_20260817_v3/codon_usage_per_gene.tsv}"

[[ -n "${TABLE}" ]] || { echo "ERROR: no table. Pass one as argument 1 or set ANNOTATED_TABLE in config.sh." >&2; exit 1; }
[[ -s "${TABLE}" ]] || { echo "ERROR: missing or empty: ${TABLE}" >&2; exit 1; }
[[ -s "${AS_CDS}" ]] || { echo "ERROR: missing or empty CDS FASTA: ${AS_CDS} (set AS_CDS in config.sh)" >&2; exit 1; }

# The gene and region columns must exist, or nothing joins.
for col in gene1 region1; do
  head -1 "${TABLE}" | tr '\t' '\n' | grep -q -x -F "${col}" || {
    echo "ERROR: no ${col} column in ${TABLE}; is this annotate_regions.py output?" >&2
    exit 1; }
done

CC_ARG=()
if [[ -s "${CODON_PER_GENE_TABLE}" ]]; then
  CC_ARG=(--codon-per-gene "${CODON_PER_GENE_TABLE}")
else
  echo "note: no codon_usage_bias per-gene table at ${CODON_PER_GENE_TABLE};"
  echo "      the GC3s cross-check against Figure 2 is skipped."
  d="$(dirname "${CODON_PER_GENE_TABLE}")"
  if [[ -d "${d}" ]]; then
    echo "      .tsv files in ${d}:"
    find "${d}" -maxdepth 1 -name '*.tsv' -printf '        %f\n' 2>/dev/null || true
  fi
fi

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/compositional_analysis_${STAMP}"
printf '%-8s %s\n' "Table:" "${TABLE}" "CDS:" "${AS_CDS}" "Output:" "${OUTDIR}"

# S1 is emitted whole by this one script, so it carries its own panel letters.
"${PY}" "${SCRIPT}" --annotated "${TABLE}" --cds "${AS_CDS}" \
  "${CC_ARG[@]}" --panel-letters --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/composition_summary.txt"
