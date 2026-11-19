#!/bin/bash
#SBATCH -J pi_region
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/pi_%j.out
#SBATCH -e logs/pi_%j.err

# ---------------------------------------------------------------------------
# piN/piS and the MK test, by region class.
#
#   sbatch run_pi_by_region.sh filtered_miss4     # primary
#   sbatch run_pi_by_region.sh filtered           # sensitivity (max-missing 1)
#
# Single-threaded (numpy + a bcftools query pipe), but MEMORY HUNGRY: the
# genome is held as strings (~276 MB) and the codon index is built from a
# Python list of ~14 M tuples before conversion to numpy, which peaks several
# GB. Hence --mem=48G and 1 CPU. If it OOMs, that list is the place to look.
#
# READ THE "REF/codon agreement" LINE IN THE LOG. It validates coordinates,
# strand and reading frame simultaneously against data the script did not
# generate. The script exits non-zero above 0.1% disagreement; anything above
# ~0.01% deserves investigation before the numbers are believed.
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
REQUIRED=(PROJECT_DIR POPGEN_DIR POPGEN_WORK_DIR AS_V50_DIR ANNOTATED_TABLE)
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
SCRIPT="$(find_script pi_by_region.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: pi_by_region.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
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
V50="${AS_V50_DIR}"
MINIPROT="${PROJECT_DIR}"
# The annotated rate table the other analyses read: ANNOTATED_TABLE in config.sh,
# relative to PROJECT_DIR. The published run read the May table,
# kaks_regions/kaks_annotated.tsv; its gene1, Ka and Ks columns, all this script
# takes from it, are identical to those of the 2026-08-17 table for all 9,407 genes.
case "${ANNOTATED_TABLE}" in
  /*) KAKS_TABLE="${ANNOTATED_TABLE}" ;;
  *)  KAKS_TABLE="${PROJECT_DIR}/${ANNOTATED_TABLE}" ;;
esac

python "${SCRIPT}" \
  --vcf "${POPGEN}/vcf/${SET}/snps_invariant.bcf" \
  --genome "${V50}/AG.v50.fa" \
  --gff "${MINIPROT}/As_miniprot.gff" \
  --gene-callable "${POPGEN}/targets/gene_callable.tsv" \
  --kaks "${KAKS_TABLE}" \
  --min-sites 100 \
  --n-boot 2000 \
  --outdir "${POPGEN}/results/pi_${SET}"
