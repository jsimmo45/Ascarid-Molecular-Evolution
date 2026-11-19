#!/bin/bash
#SBATCH -J as_map
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --array=1-3
#SBATCH -o logs/map_%A_%a.out
#SBATCH -e logs/map_%A_%a.err

# ---------------------------------------------------------------------------
# Map PRJNA1195031 runs to the AG.v50 GERMLINE assembly.
#
#   Pilot (default):  --array=1-3  with SUBSAMPLE=3000000  -> ~20 min/run
#   Full set:         sbatch --array=1-64 --time=08:00:00 run_map_reads.sh
#                     and set SUBSAMPLE=0 below.
#
# 16 CPUs: bwa-mem2 scales well to ~16; samtools sort gets half of them at
# 1G/thread, so peak RSS is ~18-20 GB (8 GB sort + ~10 GB bwa-mem2 index).
# Do not raise --cpus-per-task without lowering `-m` in map_reads.py.
#
# Build the bwa-mem2 index ONCE before the first array (~10 GB RAM, ~25 min):
#   sbatch --wrap="bwa-mem2 index $AS_GENOME_FA" \
#     -A ACF-UTK0146 --partition=condo-jwang166 --qos=condo \
#     --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=24G --time=02:00:00
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

POPGEN="${POPGEN_DIR}"
GENOME="${AS_GENOME_FA}"
SUBSAMPLE=3000000            # 0 = map all reads

export TMPDIR="${POPGEN}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

python "${SCRIPT}" \
  --manifest "${POPGEN}/reads/ena_manifest.tsv" \
  --reads-dir "${POPGEN}/reads" \
  --genome "$GENOME" \
  --outdir "${POPGEN}/bam" \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --subsample-pairs "${SUBSAMPLE}"

rm -rf "$TMPDIR"
