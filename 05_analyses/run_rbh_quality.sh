#!/bin/bash
#SBATCH -J s18_rbhqual
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=s18_rbhqual_%j.out
#SBATCH --error=s18_rbhqual_%j.err
#
# Supplementary fig. S18 (S10 before the 2026-10-12 renumbering; the output
# stem figS10_rbh_quality is kept, and the assembly copies it as FigS18):
# quality of the six reciprocal-best-hit ortholog sets (protein identity;
# forward/reverse bitscore agreement), drawn from the existing
# reciprocal_best_hits.py outputs -- no BLAST is rerun -- and restricted to the
# pairs that entered each rate table, so the n match Methods.
#
# 2026-10-12: panel letters A and B; the cluster names in A placed clear of
# every curve; log2 set as log₂; the summary gives panel B's quartiles and
# whiskers for the legend. Figure only -- every number is unchanged.
#
#   sbatch run_rbh_quality.sh
#
# Inputs (config.sh; relative paths are taken from PROJECT_DIR):
#   RBH_AS_PU ... RBH_TOX_BSC   orthologs.tsv of each pair, as for
#                               run_branch_dnds_four.sh; rbh_full.tsv and
#                               rbh_stats.txt are read from the same folders
#   KAKS_AS_PU_TABLE ...        each pair's rate table; default
#   KAKS_TOX_BSC_TABLE          <kaks folder>/kaks_filtered.tsv (below). A
#                               missing one is skipped with a note; the pairs
#                               are then not ordered by median dS.
#
# Single-threaded, six tables of ~10,000 rows: seconds. 1 node / 1 task / 1 CPU.

set -euo pipefail

# Under sbatch, SLURM copies the batch script to a spool directory on the compute
# node, so BASH_SOURCE points at /var/spool/slurm/... and neither config.sh nor
# the .py file is anywhere near it. SLURM_SUBMIT_DIR is the directory sbatch was
# run from, which is on the shared filesystem, so search from there as well.
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
  echo "ERROR: no config.sh found." >&2
  echo "  searched upward from: ${HERE}" >&2
  echo "                        ${SLURM_SUBMIT_DIR:-<not under sbatch>}" >&2
  echo "                        ${PWD}" >&2
  echo "  Set KAKS_CONFIG=/full/path/to/config.sh and resubmit." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" \
           "${REPO_ROOT}/06_figures" "${REPO_ROOT}" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script rbh_quality.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: rbh_quality.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHONPATH

# The kaks environment, with its bin/ first on PATH: on Isaac compute nodes the
# Spack python otherwise shadows the conda one. config.sh can override all three.
CONDA_SH="${CONDA_SH:-}"
CONDA_ENV="${CONDA_ENV:-kaks}"
CONDA_BIN_PREPEND="${CONDA_BIN_PREPEND:-}"
if [[ -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
fi
PY="$(command -v python || command -v python3)"

# Probe the imports: paper_colors, region_lib and rate_tables live in lib/.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("paper_colors", "region_lib", "rate_tables", "scipy",
                       "matplotlib")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, scipy, matplotlib all importable")
PYEOF

cd "${PROJECT_DIR}"

RBH_AS_PU="${RBH_AS_PU:-rbh_results/orthologs.tsv}"
RBH_AS_BSC="${RBH_AS_BSC:-rbh_As_Bsc_v2/orthologs.tsv}"
RBH_PU_BSC="${RBH_PU_BSC:-rbh_Pu_Bsc_v2/orthologs.tsv}"
RBH_TOX_AS="${RBH_TOX_AS:-rbh_Tox_As/orthologs.tsv}"
RBH_TOX_PU="${RBH_TOX_PU:-rbh_Tox_Pu/orthologs.tsv}"
RBH_TOX_BSC="${RBH_TOX_BSC:-rbh_Tox_Bsc_v2/orthologs.tsv}"
KAKS_AS_PU_TABLE="${KAKS_AS_PU_TABLE:-kaks_results/kaks_filtered.tsv}"
KAKS_AS_BSC_TABLE="${KAKS_AS_BSC_TABLE:-kaks_As_Bsc_v2/kaks_filtered.tsv}"
KAKS_PU_BSC_TABLE="${KAKS_PU_BSC_TABLE:-kaks_Pu_Bsc_v2/kaks_filtered.tsv}"
KAKS_TOX_AS_TABLE="${KAKS_TOX_AS_TABLE:-kaks_Tox_As/kaks_filtered.tsv}"
KAKS_TOX_PU_TABLE="${KAKS_TOX_PU_TABLE:-kaks_Tox_Pu/kaks_filtered.tsv}"
KAKS_TOX_BSC_TABLE="${KAKS_TOX_BSC_TABLE:-kaks_Tox_Bsc_v2/kaks_filtered.tsv}"

ARGS=()
missing=0
add_pair() {                        # add_pair <label> <orthologs.tsv> <rate table>
  local rbh; rbh="$(dirname "$2")/rbh_full.tsv"
  if [[ ! -s "${rbh}" ]]; then
    echo "ERROR: missing or empty: ${rbh} ($1)" >&2; missing=1; return
  fi
  if [[ -s "$3" ]]; then
    ARGS+=(--pair "$1" "${rbh}" "$3")
    printf '  %-30s %s  +  %s\n' "$1" "${rbh}" "$3"
  else
    ARGS+=(--pair "$1" "${rbh}")
    printf '  %-30s %s  (no rate table at %s)\n' "$1" "${rbh}" "$3"
  fi
}
echo "pairs:"
add_pair "Ascaris × P. univalens"       "${RBH_AS_PU}"   "${KAKS_AS_PU_TABLE}"
add_pair "Ascaris × B. schroederi"      "${RBH_AS_BSC}"  "${KAKS_AS_BSC_TABLE}"
add_pair "P. univalens × B. schroederi" "${RBH_PU_BSC}"  "${KAKS_PU_BSC_TABLE}"
add_pair "T. canis × Ascaris"           "${RBH_TOX_AS}"  "${KAKS_TOX_AS_TABLE}"
add_pair "T. canis × P. univalens"      "${RBH_TOX_PU}"  "${KAKS_TOX_PU_TABLE}"
add_pair "T. canis × B. schroederi"     "${RBH_TOX_BSC}" "${KAKS_TOX_BSC_TABLE}"
[[ "${missing}" -eq 0 ]] || { echo "  (set RBH_* in config.sh)" >&2; exit 1; }

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/rbh_quality_${STAMP}"
echo "output: ${OUTDIR}"

"${PY}" "${SCRIPT}" "${ARGS[@]}" --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/rbh_quality_summary.txt"
