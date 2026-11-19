#!/bin/bash
#SBATCH -J pos_sel
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=pos_sel_%j.out
#SBATCH --error=pos_sel_%j.err
#
# Supplementary Figure S8A and Supplementary Table T3 -- genes with dN/dS > 1.
#
#   sbatch run_positive_selection.sh
#   bash   run_positive_selection.sh
#
# S8 pairs this with S8B_relaxation_k.* from run_cbr_changepoint_relax.sh; the
# two panels are assembled in Illustrator.
#
# NOTE ON NAMING. These are genes whose pairwise dN/dS point estimate exceeds 1.
# That is a threshold, not a test -- pairwise yn00 offers no LRT against
# omega = 1. The script reports a dS sensitivity table precisely because the
# omega > 1 tail is where estimation noise concentrates. Read that table before
# writing the sentence.
#
# Trivial runtime; the whole job is a Fisher test and two small figures.

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
SCRIPT="$(find_script positive_selection.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: positive_selection.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

# region_lib.py must sit beside it -- positive_selection.py and gene_density_cbr.py
# both import it, and that shared import is what keeps the CBR definition
# identical across the panels of Figure 4.
[[ -f "$(dirname "${SCRIPT}")/region_lib.py" ]] || {
  echo "ERROR: region_lib.py is not beside ${SCRIPT}." >&2; exit 1; }

export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"




if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
fi
PY="$(command -v python || command -v python3)"

cd "${PROJECT_DIR}"

TABLE="${ANNOTATED_TABLE:-}"
if [[ -z "${TABLE}" ]]; then
  echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2
  exit 1
fi
SPECIES="${SPECIES1:-Ascaris}"
THRESHOLD="${OMEGA_THRESHOLD:-1.0}"
MIN_DS="${MIN_DS:-0}"

STAMP="$(date +%Y%m%d)"
OUTDIR="${RESULTS_DIR}/positive_selection_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-12s %s\n' "Table:" "${TABLE}" "Species:" "${SPECIES}" \
  "Threshold:" "dN/dS > ${THRESHOLD}" "Min dS:" "${MIN_DS}" "Output:" "${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
[[ -s "${TABLE}" ]] || { echo "ERROR: missing or empty: ${TABLE}" >&2; exit 1; }
for col in Ka_Ks Ka Ks region1; do
  head -1 "${TABLE}" | tr '\t' '\n' | grep -qx "${col}" || {
    echo "ERROR: ${TABLE} has no ${col} column." >&2; exit 1; }
done
printf '  %-44s %8d rows\n' "$(basename "${TABLE}")" "$(( $(wc -l < "${TABLE}") - 1 ))"
if head -1 "${TABLE}" | tr '\t' '\n' | grep -qx 'align_len'; then
  echo "  align_len present -- the dS/length diagnostic will be drawn"
else
  echo "  NOTE: no align_len column; the diagnostic will show dS only"
fi

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --annotated "${TABLE}" \
  --threshold "${THRESHOLD}" \
  --min-ds "${MIN_DS}" \
  --ds-cuts 0 0.01 0.02 0.05 0.10 \
  --species "${SPECIES}" \
  --width single \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/positive_selection_summary.txt"
