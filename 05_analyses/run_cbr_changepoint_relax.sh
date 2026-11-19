#!/bin/bash
#SBATCH -J cbr_relax
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=cbr_relax_%j.out
#SBATCH --error=cbr_relax_%j.err
#
# Figure 3D -- per-arm slopes and the step in constraint at the CBR.
# Supplementary fig. S3B (file S8B_relaxation_k) -- the relaxation exponent k.
#
#   sbatch run_cbr_changepoint_relax.sh
#   bash   run_cbr_changepoint_relax.sh
#
# Replaces run_cbr_changepoint_relax.sh's old form. Three changes worth knowing:
#   * It reads signed_cbr_dist1 from the annotated table instead of re-signing
#     min_cbr_dist1 by category, so 3A, 3B and 3D share one CBR set and one sign.
#   * Panel C (the typeset equations) is gone -- the models belong in Methods.
#   * Panel B is now its own figure, S3B (file stem S8B_relaxation_k).
#
# --slope-window doubles as the drawn x-limit, so the axis stops where the model
# was fitted. --terminal-tol must match annotate_regions.py or the cross-check
# fails, which is what it is for.
#
# Single-threaded scipy: 1 node / 1 task / 1 CPU. The bootstrap refits two
# regressions --n-boot times; 2,000 takes well under a minute.

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
SCRIPT="$(find_script cbr_changepoint_relax.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: cbr_changepoint_relax.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

# region_lib.py must sit beside it -- cbr_changepoint_relax.py and gene_density_cbr.py
# both import it, and that shared import is what keeps the CBR definition
# identical across the panels of Figure 3.
[[ -f "$(dirname "${SCRIPT}")/region_lib.py" ]] || {
  echo "ERROR: region_lib.py is not beside ${SCRIPT}." >&2; exit 1; }

export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

# The 2026-10-30 version draws the 10 of log10(dN/dS) as its own Arial text
# (Fig. 3D y label, S3B x label) through lib/plain_text.py; an older script
# writes Unicode subscript digits, which Illustrator sets in Myriad Pro.
grep -q '^def _set_log_omega_label' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} predates 2026-10-30 (no _set_log_omega_label); copy" >&2
  echo "  in the new cbr_changepoint_relax.py first." >&2; exit 1; }
# plain_text.py is imported from the first of these on PYTHONPATH.
PLAIN_TEXT=""
for c in "$(dirname "${SCRIPT}")/plain_text.py" "${REPO_ROOT}/lib/plain_text.py"; do
  [[ -f "$c" ]] && { PLAIN_TEXT="$c"; break; }
done
grep -q '^def replace_ylabel' "${PLAIN_TEXT:-/dev/null}" || {
  echo "ERROR: ${PLAIN_TEXT:-plain_text.py (not found)} has no replace_ylabel;" >&2
  echo "  copy in the 2026-10-30 lib/plain_text.py first." >&2; exit 1; }

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

TABLE="${ANNOTATED_TABLE:-}"
if [[ -z "${TABLE}" ]]; then
  echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2
  echo "  Point it at the annotate_regions.py run, e.g." >&2
  echo "    ANNOTATED_TABLE=\"results/annotate_regions_YYYYMMDD/kaks_annotated.tsv\"" >&2
  exit 1
fi
ELIM="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
SIZES="${AS_CHROM_SIZES:-}"
SPECIES="${SPECIES1:-Ascaris}"
TERMINAL_TOL="${TERMINAL_TOL:-1000}"
SLOPE_WINDOW="${SLOPE_WINDOW:-150}"
OMEGA_MAX="${OMEGA_MAX:-10}"
N_BOOT="${N_BOOT:-2000}"

# Time in the stamp, so a second run on one day cannot overwrite the first.
# run_cbr_edge_check.sh finds the newest cbr_changepoint_* by time either way.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/cbr_changepoint_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-12s %s\n' "Table:" "${TABLE}" "Elim:" "${ELIM}" \
  "Sizes:" "${SIZES:-<none>}" "Species:" "${SPECIES}" \
  "Window:" "+/- ${SLOPE_WINDOW} kb" "omega max:" "${OMEGA_MAX}" \
  "Bootstrap:" "${N_BOOT}" "TermTol:" "${TERMINAL_TOL} bp" "Output:" "${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${TABLE}" "${ELIM}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done
if ! head -1 "${TABLE}" | tr '\t' '\n' | grep -qx 'signed_cbr_dist1'; then
  echo "ERROR: ${TABLE} has no signed_cbr_dist1 column." >&2
  echo "  It was written by the old annotate_kaks_regions.py. Run" >&2
  echo "  run_annotate_regions.sh and repoint ANNOTATED_TABLE at its output." >&2
  exit 1
fi
for col in Ka_Ks region1 chrom1 midpoint1; do
  head -1 "${TABLE}" | tr '\t' '\n' | grep -qx "${col}" || {
    echo "ERROR: ${TABLE} has no ${col} column." >&2; exit 1; }
done
printf '  %-44s %8d rows\n' "$(basename "${TABLE}")" "$(( $(wc -l < "${TABLE}") - 1 ))"

SIZE_ARG=()
[[ -n "${SIZES}" ]] && SIZE_ARG=(--chrom-sizes "${SIZES}")

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --annotated "${TABLE}" \
  --eliminated-bed "${ELIM}" \
  "${SIZE_ARG[@]}" \
  --terminal-tol "${TERMINAL_TOL}" \
  --slope-window "${SLOPE_WINDOW}" \
  --omega-max "${OMEGA_MAX}" \
  --n-boot "${N_BOOT}" \
  --n-bins 24 \
  --seed 0 \
  --species "${SPECIES}" \
  --width single \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/cbr_changepoint_summary.txt"
