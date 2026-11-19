#!/bin/bash
#SBATCH -J dnds_decomp
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=dnds_decomposition_%j.out
#SBATCH --error=dnds_decomposition_%j.err
#
# Figure 2A (dN / dS / dN-dS decomposition), Figure S3A (consensus-region
# decomposition) and Figure S3B (dN vs dS scatter).
#
# Runs two ways, unchanged:
#     sbatch run_dnds_decomposition.sh          on the cluster
#     bash   run_dnds_decomposition.sh          anywhere else
#
# Everything machine-specific is in config.sh (see config.example.sh). The SLURM
# directives above are static because sbatch reads them before any shell code
# runs, so they cannot be filled from a variable; they are ignored entirely when
# the script is run with bash.
#
# 2026-10-12: config.sh and the .py file are found from SLURM_SUBMIT_DIR and
# PWD as well as from the script's own location, which under sbatch is
# SLURM's spool copy (see below). Run it from the repo or anywhere under it.

set -euo pipefail

# Under sbatch, SLURM copies the batch script to a spool directory on the compute
# node, so BASH_SOURCE points at /var/spool/slurm/... and neither config.sh nor
# the .py file is anywhere near it (jobs 6383411-6383414 and 6383612 failed on
# exactly that, 2026-10-12). SLURM_SUBMIT_DIR is the directory sbatch was run
# from, which is on the shared filesystem, so search from there as well, as
# run_genome_trace.sh does. $KAKS_CONFIG still overrides.
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
  echo "  Copy config.example.sh to config.sh and edit it, or set" >&2
  echo "  KAKS_CONFIG=/full/path/to/config.sh and resubmit." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

# The .py file, by the same reasoning: never assume it sits beside BASH_SOURCE.
REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" \
           "${REPO_ROOT}/06_figures" "${REPO_ROOT}" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script dnds_decomposition.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: dnds_decomposition.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

# --- environment ------------------------------------------------------------
if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
else
  echo "No conda config; using python on PATH: $(command -v python || command -v python3)"
fi
PY="$(command -v python || command -v python3)"

cd "${PROJECT_DIR}"

# --- dated output directory -------------------------------------------------
# A rerun lands in a new directory rather than overwriting the run a figure was
# built from. Pass RUN_TAG=something to distinguish runs made on the same day.
STAMP="$(date +%Y%m%d)"
OUTDIR="${RESULTS_DIR}/dnds_decomposition_${STAMP}${RUN_TAG:+_${RUN_TAG}}"
mkdir -p "${OUTDIR}"
echo "Output: ${OUTDIR}"

# --- preflight --------------------------------------------------------------
ANNOTATED="${ANNOTATED_TABLE:-kaks_regions/kaks_annotated.tsv}"
if [[ ! -f "${ANNOTATED}" ]]; then
  echo "ERROR: annotated table not found: ${PROJECT_DIR}/${ANNOTATED}" >&2
  exit 1
fi
printf '  %-44s %8d rows\n' "${ANNOTATED}" "$(( $(wc -l < "${ANNOTATED}") - 1 ))"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --annotated "${ANNOTATED}" \
  --width double \
  --scatter-width single \
  --outdir "${OUTDIR}"

echo "===== DONE ====="
cat "${OUTDIR}/decomposition_summary.txt"
