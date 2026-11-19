#!/bin/bash
#SBATCH -J germ_gc3s
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=germ_gc3s_%j.out
#SBATCH --error=germ_gc3s_%j.err
#
# Codon-bias contrasts of germline-restricted genes with GC3s controlled
# (germline_codon_gc3s_adjusted.py): regression on every expressed gene with a
# GC3s spline, plus GC3s-matched restriction and elimination contrasts. Follows
# germline_restricted_codon.py (job germ_codon, 2026-10-29), whose
# elimination contrast was confounded by composition (eliminated
# germline-restricted genes GC3s 0.475, retained 0.397) and whose restriction
# contrast was not matched on GC3s. Reads that run's joined per-gene table;
# nothing upstream is rerun.
#
#   sbatch run_germline_codon_gc3s_adjusted.sh
#   sbatch run_germline_codon_gc3s_adjusted.sh <germline_restricted_codon run dir or genes table>
#
# Default input: results/germline_restricted_codon_20261029_143308/germline_codon_genes.tsv;
# if missing, the newest germline_restricted_codon run under RESULTS_DIR, named
# in the log.
#
# Optional settings (environment): KNOTS (GC3s spline knots, 5), EXPR_KNOTS (4),
# GC3S_CALIPER (0.01, as codon_usage_bias.py), CALIPER (0.25, log2 expression),
# RUN_TAG.
#
# Single-threaded numpy/pandas: 1 node / 1 task / 1 CPU. About 15 minutes,
# most of it the median-regression bootstraps.
#
# WHAT TO READ: gc3s_adjusted_summary.txt (printed at the end), "What to say":
#   the elimination contrast at equal GC3s and expression, with its reading
#   (CHANGES / HOLDS, firmer / UNRESOLVED), and the restriction contrast at
#   equal GC3s.

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
SCRIPT="$(find_script germline_codon_gc3s_adjusted.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: germline_codon_gc3s_adjusted.py not found near ${REPO_ROOT}." >&2; exit 1; }
SDIR="$(dirname "${SCRIPT}")"
echo "script: ${SCRIPT}"
# It imports germline_restricted_codon.py, which imports germline_restricted_dnds.py
# (and through it functional_enrichment, multivariate_insulation, workbook_model_bridge).
for m in germline_restricted_codon.py germline_restricted_dnds.py \
         functional_enrichment.py multivariate_insulation.py workbook_model_bridge.py; do
  [[ -f "${SDIR}/${m}" ]] || { echo "ERROR: ${m} is not beside ${SCRIPT}." >&2; exit 1; }
done
export PYTHONPATH="${SDIR}:${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS}" MKL_NUM_THREADS="${OMP_NUM_THREADS}"

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
GENES="$(resolve "${1:-}" germline_codon_genes.tsv \
         germline_restricted_codon_20261029_143308 'germline_restricted_codon_*')"

KNOTS="${KNOTS:-5}"
EXPR_KNOTS="${EXPR_KNOTS:-4}"
GC3S_CALIPER="${GC3S_CALIPER:-0.01}"
CALIPER="${CALIPER:-0.25}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="${RESULTS_DIR}/germline_codon_gc3s_adjusted_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-12s %s\n' "Genes:" "${GENES}" "Knots:" "GC3s ${KNOTS}, expression ${EXPR_KNOTS}" \
  "Calipers:" "GC3s ${GC3S_CALIPER}, log2 expression ${CALIPER}" "Output:" "${OUTDIR}"
echo ""
for col in gene1 region1 expression log_expr group_default dENC GC3s ENC_at_ceiling; do
  head -1 "${GENES}" | tr '\t' '\n' | grep -qx "${col}" || {
    echo "ERROR: ${GENES} has no ${col} column." >&2; exit 1; }
done

mkdir -p "${OUTDIR}"
"${PY}" "${SCRIPT}" \
  --genes "${GENES}" \
  --metric dENC \
  --knots "${KNOTS}" \
  --expr-knots "${EXPR_KNOTS}" \
  --gc3s-caliper "${GC3S_CALIPER}" \
  --caliper "${CALIPER}" \
  --n-boot 2000 \
  --n-boot-median 1000 \
  --n-match 500 \
  --n-match-boot 1000 \
  --seed 0 \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/gc3s_adjusted_summary.txt"
