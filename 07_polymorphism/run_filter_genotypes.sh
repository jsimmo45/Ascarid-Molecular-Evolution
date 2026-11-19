#!/bin/bash
#SBATCH -J gtfilter
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/gtfilter_%j.out
#SBATCH -e logs/gtfilter_%j.err

# ---------------------------------------------------------------------------
# Filter the all-sites callset and split SNPs+invariant from indels.
#
# Depth thresholds assume BAMs matched to 13x by depth_match.py:
#   --min-dp 6   ~0.45x the matched mean
#   --max-dp 30  ~2.3x  the matched mean; excess depth AFTER matching is the
#                signature of a collapsed duplication, not of good data
#
# --max-missing 4 allows up to four of the ten worms to be missing at a site
# (genotypes already masked on DP/GQ count as missing). These are the settings
# of the published set, vcf/filtered_miss4 (12,434,792 sites), confirmed by
# rerunning them on 2026-10-26; an earlier run with --max-missing 1 wrote
# vcf/filtered.
#
# Read the per-class allele-balance drop rates in the log. If eliminated sites
# are dropped far more often than retained ones, residual collapsed paralogy is
# real and every downstream rate must be per callable site.
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
SCRIPT="$(find_script filter_genotypes.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: filter_genotypes.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
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

python "${SCRIPT}" \
  --vcf "${POPGEN}/vcf/allsites.bcf" \
  --eliminated-bed "${V50}/AG.v50.eliminated_strict.bed" \
  --outdir "${POPGEN}/vcf/filtered_miss4" \
  --min-dp 6 --max-dp 30 --min-gq 20 \
  --max-missing 4 \
  --min-mq 40 --max-mq0f 0.10 \
  --indel-buffer 5 \
  --min-ab 0.25 --max-ab 0.75 --min-het 3 \
  --threads "${SLURM_CPUS_PER_TASK}"
