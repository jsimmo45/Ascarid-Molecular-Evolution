#!/bin/bash
#SBATCH -J depthmatch
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/dmatch_%A_%a.out
#SBATCH -e logs/dmatch_%A_%a.err

# ---------------------------------------------------------------------------
# Equalise depth across region classes and worms before genotyping.
#
#   sbatch --array=1-10 run_depth_match.sh cohort_20260820 8
#                                          ^cohort dir     ^target depth
#
# SET THE TARGET DEPTH DELIBERATELY. Run the audit on bam_full first and pick a
# value the whole cohort can reach in BOTH classes; the lowest-gf worm in the
# selection is at 0.40, so its eliminated depth is ~40% of its retained depth.
# Worms that cannot reach the target are SKIPPED (a .SKIPPED file is written),
# never left un-thinned — mixing matched and unmatched samples reintroduces the
# exact bias this step exists to remove.
#
# 8 CPUs: samtools view/sort and mosdepth all thread modestly; the job is I/O
# bound on a 15-25 GB BAM.
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
REQUIRED=(POPGEN_DIR POPGEN_WORK_DIR AS_V50_DIR)
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
SCRIPT="$(find_script depth_match.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: depth_match.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

COHORT="${1:?usage: sbatch --array=1-N run_depth_match.sh <cohort_dir> <target_depth>}"
TARGET="${2:?give a target depth, e.g. 8}"
POPGEN="${POPGEN_DIR}"
V50="${AS_V50_DIR}"

RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${COHORT}/selected_runs.txt")
[[ -z "$RUN" ]] && { echo "No run on line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
echo "Task ${SLURM_ARRAY_TASK_ID} -> ${RUN} (target ${TARGET}x)"

python "${SCRIPT}" \
  --bam "${POPGEN}/bam_full/${RUN}.bam" \
  --targets-bed "${POPGEN}/targets/callable_cds_targets.bed" \
  --eliminated-bed "${V50}/AG.v50.eliminated_strict.bed" \
  --target-depth "${TARGET}" \
  --min-mq 30 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --outdir "${POPGEN}/bam_matched"
