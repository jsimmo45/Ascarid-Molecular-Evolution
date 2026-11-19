#!/bin/bash
#SBATCH -J supp_figs
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=supp_figs_%j.out
#SBATCH --error=supp_figs_%j.err
#
# Copies every supplementary figure panel into one folder under its final number
# (first citation in the text): results/supplementary_figures_<stamp>/FigS<nn>_<panel>.*
# The analysis outputs are not touched.
#
#   sbatch run_assemble_supplementary_figures.sh
#   sbatch run_assemble_supplementary_figures.sh --check-only   # verify, copy nothing
#
# Inputs: supplementary_figure_sources.tsv beside assemble_supplementary_figures.py
# (figure<TAB>was<TAB>panel<TAB>path; relative paths resolve against PROJECT_DIR
# from config.sh, or SUPP_FIG_BASE if set). Arguments are passed to the script.
#
# File copies only: about 34 panels x 3 formats, seconds. 1 node / 1 task / 1 CPU.
#
# Exit status: 0 copied (or checked); 1 a source is missing, a number is skipped
# or repeated, or two rows collide; nothing is copied in that case.

set -euo pipefail

# Under sbatch, SLURM copies the batch script to a spool directory on the compute
# node, so BASH_SOURCE points at /var/spool/slurm/... SLURM_SUBMIT_DIR is the
# directory sbatch was run from, on the shared filesystem, so search from there too.
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
SCRIPT="$(find_script assemble_supplementary_figures.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: assemble_supplementary_figures.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
MANIFEST="${SUPP_FIG_SOURCES:-$(dirname "${SCRIPT}")/supplementary_figure_sources.tsv}"
[[ -s "${MANIFEST}" ]] || { echo "ERROR: missing or empty manifest: ${MANIFEST}" >&2; exit 1; }
echo "manifest: ${MANIFEST}"

PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHONPATH

# The manifest expands ${POPGEN_DIR} (S14-S17 live in the population-genetics
# tree). config.sh must set it.
export POPGEN_DIR="${POPGEN_DIR:?set POPGEN_DIR in config.sh}"
echo "popgen: ${POPGEN_DIR}"

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

"${PY}" -c "import region_lib" 2>/dev/null \
  || { echo "ERROR: region_lib not importable; lib/ must be on PYTHONPATH." >&2; exit 1; }

cd "${PROJECT_DIR}"

BASE="${SUPP_FIG_BASE:-${PROJECT_DIR}}"
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/supplementary_figures_${STAMP}"
printf '%-9s %s\n' "Base:" "${BASE}" "Output:" "${OUTDIR}"

rc=0
"${PY}" "${SCRIPT}" --manifest "${MANIFEST}" --base-dir "${BASE}" \
  --outdir "${OUTDIR}" "$@" || rc=$?

echo ""
if [[ -s "${OUTDIR}/supplementary_figures_summary.txt" ]]; then
  echo "===== SUMMARY (exit ${rc}) ====="
  cat "${OUTDIR}/supplementary_figures_summary.txt"
elif [[ ${rc} -eq 0 ]]; then
  echo "===== CHECK ONLY: manifest and files OK (exit 0) ====="
else
  echo "===== FAILED (exit ${rc}); nothing copied. Every problem is listed above and in the .err file. ====="
fi
exit "${rc}"
