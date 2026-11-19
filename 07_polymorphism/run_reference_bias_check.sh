#!/bin/bash
#SBATCH -J refbias
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/refbias_%j.out
#SBATCH -e logs/refbias_%j.err

# ---------------------------------------------------------------------------
# Supplementary fig. S17 — reference bias, four panels.
#
#   sbatch run_reference_bias_check.sh filtered_miss4
#
# MUST be given the PRE-allele-balance callset (step3_gapfilter.bcf). The
# [0.25, 0.75] screen trims exactly the tail panel A measures, so running this
# on the post-screen file reports a cleaner distribution than the data had.
#
# Memory and time match run_pi_by_region.sh: it builds the same codon index.
#
# WHAT TO READ IN THE LOG:
#   1. TEST 2's eliminated/retained hom-alt ratio. Below 1 means eliminated
#      regions are LESS diverged from the reference, which bounds how much
#      class-asymmetric bias is even possible. That single number does more
#      work than either of the other tests.
#   2. TEST 3's per-stratum folds, now printed on panel D as well. piN/piS
#      rising with local variant density is expected — dense clusters are the
#      misalignment signature — and does not matter provided the CLASS FOLD
#      holds within every stratum.
#   3. "FIGURE AUDIT -- clean".
#   4. Any "omitted from the panel" warning: a class with no stratum above the
#      site floor is dropped rather than drawn as an invisible line.
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
REQUIRED=(PROJECT_DIR POPGEN_DIR POPGEN_WORK_DIR AS_V50_DIR)
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
SCRIPT="$(find_script reference_bias_check.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: reference_bias_check.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

export PAPER_COLORS_DIR="${REPO_ROOT}/lib"

SET="${1:-filtered_miss4}"
POPGEN="${POPGEN_DIR}"
V50="${AS_V50_DIR}"
MINIPROT="${PROJECT_DIR}"

python "${SCRIPT}" \
  --vcf "${POPGEN}/vcf/${SET}/step3_gapfilter.bcf" \
  --genome "${V50}/AG.v50.fa" \
  --gff "${MINIPROT}/As_miniprot.gff" \
  --gene-callable "${POPGEN}/targets/gene_callable.tsv" \
  --window 150 \
  --tag "${SET}" \
  --outdir "${POPGEN}/results/refbias"
