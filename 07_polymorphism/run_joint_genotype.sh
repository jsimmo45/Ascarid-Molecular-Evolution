#!/bin/bash
#SBATCH -J jointgt
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/jointgt_%A_%a.out
#SBATCH -e logs/jointgt_%A_%a.err

# ---------------------------------------------------------------------------
# All-sites joint genotyping over the callable CDS targets, one autosome per
# array task.
#
#   sbatch --array=1-19%8 run_joint_genotype.sh          # chr01..chr19
#   sbatch --dependency=afterok:<jobid> --wrap="..."     # then merge, below
#
# Merge step (single task, after the array succeeds):
#   python joint_genotype.py --merge --outdir $POPGEN/vcf --threads 8
#
# Only chr01-chr19: the analysis is autosomes-only (X has a different Ne and
# this cohort is mixed-sex).
#
# mpileup and call each get half the CPUs since they run as a pipe. Restricting
# to ~14 Mb of callable CDS rather than the whole 276 Mb genome is what makes
# an all-sites callset affordable here.
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
SCRIPT="$(find_script joint_genotype.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: joint_genotype.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
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
V50="${AS_V50_DIR}"
CHROM=$(printf "chr%02d" "${SLURM_ARRAY_TASK_ID}")
echo "Task ${SLURM_ARRAY_TASK_ID} -> ${CHROM}"

# Depth-matched BAMs only. If bam_matched is empty, depth_match.py has not run
# and calling on bam_full would bake in the eliminated-region depth deficit.
shopt -s nullglob
BAMS=("${POPGEN}"/bam_matched/*.matched.bam)
if (( ${#BAMS[@]} == 0 )); then
  echo "No depth-matched BAMs in ${POPGEN}/bam_matched — run run_depth_match.sh first" >&2
  exit 1
fi
echo "${#BAMS[@]} depth-matched sample(s)"

python "${SCRIPT}" \
  --bams "${BAMS[@]}" \
  --genome "${V50}/AG.v50.fa" \
  --targets-bed "${POPGEN}/targets/callable_cds_targets.bed" \
  --chrom "${CHROM}" \
  --min-mq 30 \
  --min-bq 20 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --outdir "${POPGEN}/vcf"
