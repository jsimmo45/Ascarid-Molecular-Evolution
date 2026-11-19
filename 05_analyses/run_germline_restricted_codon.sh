#!/bin/bash
#SBATCH -J germ_codon
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=germ_codon_%j.out
#SBATCH --error=germ_codon_%j.err
#
# Codon bias (dENC) of eliminated vs retained germline-restricted genes
# (germline_restricted_codon.py). Jianbin Wang, 2026-10-28: would codon usage
# help distinguish eliminated genes from genes that are only repressed in the
# soma? Reads the classified genes of the S12B-C run and the per-gene codon
# table of the fig. 2B-C run; nothing upstream is rerun.
#
#   sbatch run_germline_restricted_codon.sh
#   sbatch run_germline_restricted_codon.sh <germline_restricted_dnds run> <codon_usage run>
#
# Defaults: the runs in the manuscript, results/germline_restricted_dnds_20261019_153637
# (S12B-C) and results/codon_usage_ascaris_20261027 (fig. 2B-C); if either is
# missing, the newest run of that analysis under RESULTS_DIR, named in the log.
#
# Optional settings (environment): METRIC (dENC, default | ENC | GC3s),
# CALIPER (0.25, as S12C), GC3S_CALIPER (0.02), RUN_TAG.
#
# Single-threaded pandas/scipy on two tables: 1 node / 1 task / 1 CPU. The
# matching (500 + 1,000 replicates per row, three rows) takes a few minutes.
#
# WHAT TO READ: germline_codon_summary.txt (printed at the end), "What to say":
#   the elimination contrast (eliminated - retained, both germline-restricted),
#   raw and matched, and the restriction contrast (germline-restricted -
#   somatic peak, retained). The figure germline_codon_contrasts.svg is the
#   S12C layout on dENC, a candidate panel if the result goes in.

set -euo pipefail

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
  echo "ERROR: no config.sh found. Set KAKS_CONFIG=/full/path/to/config.sh." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" \
           "${REPO_ROOT}" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script germline_restricted_codon.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: germline_restricted_codon.py not found near ${REPO_ROOT}." >&2; exit 1; }
SDIR="$(dirname "${SCRIPT}")"
echo "script: ${SCRIPT}"
# It imports S12's groups, contrasts and matching from germline_restricted_dnds.py,
# which imports functional_enrichment, multivariate_insulation, workbook_model_bridge.
for m in germline_restricted_dnds.py functional_enrichment.py \
         multivariate_insulation.py workbook_model_bridge.py; do
  [[ -f "${SDIR}/${m}" ]] || { echo "ERROR: ${m} is not beside ${SCRIPT}." >&2; exit 1; }
done
export PYTHONPATH="${SDIR}:${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
fi
PY="$(command -v python || command -v python3)"

cd "${PROJECT_DIR}"

# resolve <run dir or file> <file name> <default run> <glob of runs>
resolve() {
  local given="$1" name="$2" default="$3" glob="$4" d
  if [[ -n "${given}" ]]; then
    [[ -f "${given}" ]] && { printf '%s' "${given}"; return 0; }
    [[ -f "${given%/}/${name}" ]] && { printf '%s' "${given%/}/${name}"; return 0; }
    echo "ERROR: no ${name} at ${given}" >&2; return 1
  fi
  [[ -f "${RESULTS_DIR}/${default}/${name}" ]] && {
    printf '%s' "${RESULTS_DIR}/${default}/${name}"; return 0; }
  for d in $(ls -td "${RESULTS_DIR}"/${glob} 2>/dev/null); do
    [[ -f "${d}/${name}" ]] && {
      echo "WARNING: ${default} not found; using the newest run, ${d}" >&2
      printf '%s' "${d}/${name}"; return 0; }
  done
  echo "ERROR: no ${name} under ${RESULTS_DIR}/${glob}" >&2; return 1
}
CLASSIFIED="$(resolve "${1:-}" genes_classified.tsv \
              germline_restricted_dnds_20261019_153637 'germline_restricted_dnds_*')"
CODON="$(resolve "${2:-}" codon_usage_per_gene.tsv \
         codon_usage_ascaris_20261027 'codon_usage_ascaris_*')"

METRIC="${METRIC:-dENC}"
CALIPER="${CALIPER:-0.25}"
GC3S_CALIPER="${GC3S_CALIPER:-0.02}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="${RESULTS_DIR}/germline_restricted_codon_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-12s %s\n' "Classified:" "${CLASSIFIED}" "Codon:" "${CODON}" \
  "Metric:" "${METRIC}" "Caliper:" "${CALIPER} (log2 expr), GC3s ${GC3S_CALIPER}" \
  "Output:" "${OUTDIR}"
echo ""
for col in gene1 region1 log_expr group_default; do
  head -1 "${CLASSIFIED}" | tr '\t' '\n' | grep -qx "${col}" || {
    echo "ERROR: ${CLASSIFIED} has no ${col} column." >&2; exit 1; }
done
for col in gene_id region dENC GC3s ENC ENC_at_ceiling; do
  head -1 "${CODON}" | tr '\t' '\n' | grep -qx "${col}" || {
    echo "ERROR: ${CODON} has no ${col} column." >&2; exit 1; }
done

mkdir -p "${OUTDIR}"
"${PY}" "${SCRIPT}" \
  --classified "${CLASSIFIED}" \
  --codon "${CODON}" \
  --metric "${METRIC}" \
  --caliper "${CALIPER}" \
  --gc3s-caliper "${GC3S_CALIPER}" \
  --n-boot 2000 \
  --n-match 500 \
  --n-match-boot 1000 \
  --seed 0 \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/germline_codon_summary.txt"
