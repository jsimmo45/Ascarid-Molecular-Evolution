#!/bin/bash

#SBATCH -J build_quartet_aln
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=build_quartet_aln_%j.out
#SBATCH --error=build_quartet_aln_%j.err

# ---------------------------------------------------------------------------
# run_make_rate_table.sh
#
# Build the supplementary rate-heterogeneity table from results already on
# disk. Pure formatting -- no alignment, no PAML, seconds to run -- so this is
# a plain bash script, NOT an sbatch job. Run it on the login node:
#
#   bash run_make_rate_table.sh
#
# Reads:
#   kaks_results/kaks_filtered.tsv           As x Pu
#   kaks_As_Bsc_v2/kaks_filtered.tsv         As x Bsc
#   kaks_Pu_Bsc_v2/kaks_filtered.tsv         Pu x Bsc
#   kaks_Tox_As/kaks_filtered.tsv            As x Tox
#   kaks_Tox_Pu/kaks_filtered.tsv            Pu x Tox
#   kaks_Tox_Bsc_v2/kaks_filtered.tsv        Bsc x Tox
#   clock_test_sites{12,3,all}/clock_test_summary.txt
#
# Note the _v2 suffixes on the three Baylisascaris tables. A path without _v2
# there is the pre-fix data.
# ---------------------------------------------------------------------------

set -euo pipefail

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
REQUIRED=(PROJECT_DIR)
for v in "${REQUIRED[@]}"; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: ${v} is not set in ${CONFIG}" >&2; exit 1; }
done

find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/06_outputs" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script make_rate_table.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: make_rate_table.py not found near ${REPO_ROOT}/06_outputs." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV:?set CONDA_ENV in config.sh}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
fi

cd "${PROJECT_DIR}" || exit 1

# Branch lengths from four_way_kaks.py --fit exact. If the four-way analysis is
# ever refit, change them here -- Panel A's predicted column depends on them.
python "${SCRIPT}" \
  --branch 'Ascaris=0.0461' \
  --branch 'Parascaris=0.1496' \
  --branch 'Baylisascaris=0.1461' \
  --branch 'Toxocara=1.2666' \
  --internal 0.0777 \
  --outdir supp_rate_table

echo ""
echo "===== DONE ====="
echo "Check that Panel A's deviation column is ~0.00% on five of six pairs."
echo "Anything else means the branch lengths above no longer match the tree."
