#!/bin/bash
#SBATCH -J cbr_edge
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=cbr_edge_%j.out
#SBATCH --error=cbr_edge_%j.err
#
# Fig. 3D edge check (cbr_edge_check.py): is the retained bin next to the CBR,
# whose mean sits near the eliminated level, different from the rest of the
# retained arm? Also tests the genes whose span reaches across their CBR, and
# lists the genes of both edge bins. No figure. Jianbin Wang's question,
# 2026-10-26 and 2026-10-28.
#
#   sbatch run_cbr_edge_check.sh                      # stats of the newest figure run
#   sbatch run_cbr_edge_check.sh <cbr_changepoint run dir>
#
# The arm fits are checked against that run's cbr_changepoint_stats.tsv (the
# figure in the manuscript is results/cbr_changepoint_20261027); the job stops
# if they differ, because then it is not testing the genes the figure drew. It
# uses config.sh's ANNOTATED_TABLE, AS_ELIMINATED_BED, AS_CHROM_SIZES and
# TERMINAL_TOL exactly as run_cbr_changepoint_relax.sh does.
#
# Optional: FLAG_BEDS="a.bed b.bed" flags genes overlapping those intervals
# (for example a BED of the CBR zones themselves, or a strict retained set);
# flagged genes leave the "clear" gene set.
#
# Since 2026-10-30 it also counts the genes beyond the fig. 3D y-axis (drawn
# nowhere, used everywhere), with the figure's own rule and cap; the cap is
# read from the figure run's provenance.json, or set CAP_PCT to override.
#
# Single-threaded numpy/scipy on one table: 1 node / 1 task / 1 CPU.
#
# WHAT TO READ: cbr_edge_summary.txt (printed at the end), "What to say".
#   all       the drawn retained bin against the rest of the retained arm
#   clear     the same without genes that reach across their CBR
#   crossing  retained genes reaching into eliminated DNA against the others
#   3D legend the genes beyond the 3D y-axis, for the figure legend
# The arm-fit lines must read "same" against the figure run.

set -euo pipefail

# Under sbatch, SLURM copies the batch script to a spool directory, so search
# for config.sh from SLURM_SUBMIT_DIR and PWD as well as from here.
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
SCRIPT="$(find_script cbr_edge_check.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: cbr_edge_check.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# It imports the figure's fit and filter from cbr_changepoint_relax.py, and the
# CBR definition from region_lib.py; both must sit beside it.
for dep in cbr_changepoint_relax.py region_lib.py; do
  [[ -f "$(dirname "${SCRIPT}")/${dep}" || -f "${REPO_ROOT}/lib/${dep}" ]] || {
    echo "ERROR: ${dep} is neither beside ${SCRIPT} nor in ${REPO_ROOT}/lib." >&2
    exit 1; }
done
export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"
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
[[ -n "${TABLE}" ]] || { echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2; exit 1; }
ELIM="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
SIZES="${AS_CHROM_SIZES:-}"
TERMINAL_TOL="${TERMINAL_TOL:-1000}"
SLOPE_WINDOW="${SLOPE_WINDOW:-150}"
OMEGA_MAX="${OMEGA_MAX:-10}"
N_BOOT="${N_BOOT:-2000}"

# The figure run to check against: argument 1, else the newest under RESULTS_DIR.
FIGRUN="${1:-}"
if [[ -z "${FIGRUN}" ]]; then
  FIGRUN="$(ls -td "${RESULTS_DIR}"/cbr_changepoint_* 2>/dev/null | head -1 || true)"
fi
STATS="${FIGRUN%/}/cbr_changepoint_stats.tsv"
[[ -s "${STATS}" ]] || { echo "ERROR: no cbr_changepoint_stats.tsv in '${FIGRUN}'." >&2
                         echo "  Pass the figure run directory as the first argument." >&2; exit 1; }

STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="${RESULTS_DIR}/cbr_edge_check_${STAMP}"

printf '%-12s %s\n' "Table:" "${TABLE}" "Elim:" "${ELIM}" \
  "Sizes:" "${SIZES:-<none>}" "Figure run:" "${FIGRUN}" \
  "Window:" "+/- ${SLOPE_WINDOW} kb" "omega max:" "${OMEGA_MAX}" \
  "Flag BEDs:" "${FLAG_BEDS:-<none>}" "Output:" "${OUTDIR}"
echo ""

for f in "${TABLE}" "${ELIM}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done
for col in signed_cbr_dist1 Ka_Ks region1 chrom1 midpoint1; do
  head -1 "${TABLE}" | tr '\t' '\n' | grep -qx "${col}" || {
    echo "ERROR: ${TABLE} has no ${col} column." >&2; exit 1; }
done
for col in start1 end1; do
  head -1 "${TABLE}" | tr '\t' '\n' | grep -qx "${col}" || \
    echo "WARNING: ${TABLE} has no ${col} column: spans cannot be placed against the CBR." >&2
done

SIZE_ARG=()
[[ -n "${SIZES}" ]] && SIZE_ARG=(--chrom-sizes "${SIZES}")
CAP_ARG=()
[[ -n "${CAP_PCT:-}" ]] && CAP_ARG=(--cap-pct "${CAP_PCT}")
grep -q '^def off_axis' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} has no off_axis(): it is the version before 2026-10-30;" \
       "copy the new cbr_edge_check.py beside this wrapper." >&2; exit 1; }
FLAG_ARGS=()
for b in ${FLAG_BEDS:-}; do
  [[ -s "$b" ]] || { echo "ERROR: FLAG_BEDS entry missing or empty: $b" >&2; exit 1; }
  FLAG_ARGS+=(--flag-bed "$b")
done

mkdir -p "${OUTDIR}"

"${PY}" "${SCRIPT}" \
  --annotated "${TABLE}" \
  --eliminated-bed "${ELIM}" \
  "${SIZE_ARG[@]}" \
  --terminal-tol "${TERMINAL_TOL}" \
  --stats "${STATS}" \
  --slope-window "${SLOPE_WINDOW}" \
  --omega-max "${OMEGA_MAX}" \
  --n-bins 24 \
  --min-bin 3 \
  --n-boot "${N_BOOT}" \
  --seed 0 \
  "${CAP_ARG[@]}" \
  "${FLAG_ARGS[@]}" \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/cbr_edge_summary.txt"
