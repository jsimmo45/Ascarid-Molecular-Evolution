#!/bin/bash
#SBATCH -J dupcheck
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/dupcheck_%j.out
#SBATCH -e logs/dupcheck_%j.err

# ---------------------------------------------------------------------------
# Are any two of the ten libraries the same worm?
#
#   sbatch run_check_duplicate_individuals.sh filtered_miss4
#
# Landeryou et al. claim 54 worms for 54 runs, but the submission carries ONE
# BioSample for all of them and Supplementary Data 2 is column-scrambled, so
# nothing in the metadata confirms one worm per run. If two of our ten are the
# same individual, n is really nine, that worm is double-counted in every
# allele frequency, and every bootstrap CI is slightly too narrow.
#
# Cheap: 200k SNPs, 45 pairs, one core, minutes. Small --mem and one CPU are
# correct here — the genotype matrix is 200k x 10 int8.
#
# WHAT TO READ IN THE LOG: the VERDICT line. Expect a clean separation —
# unrelated Ascaris pairs disagree at roughly half of biallelic SNP genotypes,
# a same-worm pair at one or two percent. There is no grey zone in between, so
# a marginal call would itself be evidence something is wrong upstream.
#
# Also check the "relatedness" range. It is rescaled so the median pair reads
# 0; a full sib would read near 0.5. Worms from one host are not independent
# draws even when they are distinct individuals, and that is worth a methods
# sentence if it shows up.
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
REQUIRED=(POPGEN_DIR POPGEN_WORK_DIR)
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
SCRIPT="$(find_script check_duplicate_individuals.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: check_duplicate_individuals.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
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

python "${SCRIPT}" \
  --vcf "${POPGEN}/vcf/${SET}/snps_invariant.bcf" \
  --max-sites 200000 \
  --outdir "${POPGEN}/results/duplicates_${SET}"
