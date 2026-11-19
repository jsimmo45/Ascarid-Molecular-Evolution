#!/bin/bash
#SBATCH -J as_mapfull
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/mapfull_%A_%a.out
#SBATCH -e logs/mapfull_%A_%a.err

# ---------------------------------------------------------------------------
# Full-depth mapping for the germline-positive cohort only.
#
# Driven by cohort_*/selected_runs.txt, ONE RUN PER LINE, indexed directly by
# the array task id — not by manifest row number. The row-index route already
# cost one wrong submission; a run accession read from the selection file is
# self-evidently correct in the log.
#
# Usage:
#   sbatch --array=1-$(wc -l < cohort_20260820/selected_runs.txt)%6 run_map_full.sh cohort_20260820
#
# 10 h and 16 CPUs: a full 25-45M-pair library is ~10-15x the pilot subsample,
# so budget ~1-3 h each; the ceiling is for the deepest. %6 keeps the condo
# usable — these write ~15-25 GB per BAM.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Configuration. Every machine-specific path and the conda setup come from
# config.sh at the top of the repository (copy config.example.sh). Under sbatch
# the batch script runs from a spool copy, so config.sh is searched for upward
# from this file, from the submission directory and from $PWD;
# KAKS_CONFIG=/full/path/to/config.sh overrides the search.
# ---------------------------------------------------------------------------
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
  echo "ERROR: no config.sh found (copy config.example.sh to config.sh)." >&2
  echo "  searched upward from: ${HERE}" >&2
  echo "                        ${SLURM_SUBMIT_DIR:-<not under sbatch>}" >&2
  echo "                        ${PWD}" >&2
  echo "  Set KAKS_CONFIG=/full/path/to/config.sh and resubmit." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=../config.example.sh
source "${CONFIG}"
REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
REQUIRED=(POPGEN_DIR POPGEN_WORK_DIR AS_GENOME_FA)
for v in "${REQUIRED[@]}"; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: ${v} is not set in ${CONFIG}" >&2; exit 1; }
done

find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/07_polymorphism" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script map_reads.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: map_reads.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

COHORT="${1:?usage: sbatch --array=1-N run_map_full.sh <cohort_dir>}"
POPGEN="${POPGEN_DIR}"
GENOME="${AS_GENOME_FA}"

RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${COHORT}/selected_runs.txt")
if [[ -z "$RUN" ]]; then
  echo "No run on line ${SLURM_ARRAY_TASK_ID} of ${COHORT}/selected_runs.txt" >&2
  exit 1
fi
echo "Task ${SLURM_ARRAY_TASK_ID} -> ${RUN}"

export TMPDIR="${POPGEN}/tmp/full_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

python "${SCRIPT}" \
  --manifest "${POPGEN}/reads/ena_manifest.tsv" \
  --reads-dir "${POPGEN}/reads" \
  --genome "$GENOME" \
  --outdir "${POPGEN}/bam_full" \
  --run "$RUN" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --subsample-pairs 0

rm -rf "$TMPDIR"
