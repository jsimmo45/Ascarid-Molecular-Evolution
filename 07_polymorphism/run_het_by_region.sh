#!/bin/bash
#SBATCH -J hetregion
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/hetregion_%j.out
#SBATCH -e logs/hetregion_%j.err

# ---------------------------------------------------------------------------
# Per-sample heterozygosity split by region class.
#
#   sbatch run_het_by_region.sh filtered_miss4
#
# Run this BEFORE run_jackknife_individuals.sh: it tells you which libraries
# (if any) belong in that script's --drop-set. The jackknife cannot find a
# shared artifact on its own; this is what identifies the group.
#
# Light compared with pi_by_region.py — a flat CDS interval map, not the codon
# index, and by default only segregating sites — hence 12G and one hour rather
# than 48G and four.
#
# WHAT TO READ IN THE LOG: the eli/ret column, NOT the absolute het rates. The
# ratio is internally controlled. A library that is globally het-inflated has a
# high het rate in both classes and a normal RATIO, and cannot bias a contrast
# between the classes — it moves numerator and denominator together. A library
# with an off-median ratio is inflating one class only, and that is the one
# that can manufacture the result.
#
# Then read the allele-balance block. Heterozygous calls centring well below
# 0.5 in eliminated regions are collapsed paralogs, not heterozygotes: the
# "alternate allele" is a minority read stack from a second copy.
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
REQUIRED=(PROJECT_DIR POPGEN_DIR POPGEN_WORK_DIR)
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
SCRIPT="$(find_script het_by_region.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: het_by_region.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

SET="${1:-filtered_miss4}"
POPGEN="${POPGEN_DIR}"
MINIPROT="${PROJECT_DIR}"

python "${SCRIPT}" \
  --vcf "${POPGEN}/vcf/${SET}/snps_invariant.bcf" \
  --gff "${MINIPROT}/As_miniprot.gff" \
  --gene-callable "${POPGEN}/targets/gene_callable.tsv" \
  --ratio-tolerance 0.25 \
  --outdir "${POPGEN}/results/het_by_region_${SET}"
