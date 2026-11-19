#!/bin/bash
#SBATCH -J x_auto_chk
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=x_auto_chk_%j.out
#SBATCH --error=x_auto_chk_%j.err
#
# Ascaris X chromosomes vs autosomes (dS, dN, dN/dS on retained genes, tested
# on genes and on chromosomes), the male-to-female mutation ratio implied by
# the dS gap, and the eliminated-versus-retained contrast within autosomes and
# within the X. Decides the optional sex-chromosome sentence (results draft,
# change note 6) and rules out X/autosome composition as a confounder of the
# main result.
#
#   sbatch run_x_autosome_check.sh
#   sbatch run_x_autosome_check.sh <path/to/kaks_annotated.tsv>
#
# Input: argument 1, else ANNOTATED_TABLE from config.sh.
#
# Single-threaded pandas/scipy on ~9,400 rows; the exact chromosome-level
# permutation is 42,504 assignments, a few seconds. 1 node / 1 task / 1 CPU.

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
SCRIPT="$(find_script x_autosome_check.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: x_autosome_check.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

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

# Probe the imports: region_lib and rate_tables live in lib/.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("region_lib", "rate_tables", "scipy")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("modules: region_lib, rate_tables, scipy all importable")
PYEOF

cd "${PROJECT_DIR}"

TABLE="${1:-${ANNOTATED_TABLE:-}}"
[[ -n "${TABLE}" ]] || { echo "ERROR: no table. Pass one as argument 1 or set ANNOTATED_TABLE in config.sh." >&2; exit 1; }
[[ -s "${TABLE}" ]] || { echo "ERROR: missing or empty: ${TABLE}" >&2; exit 1; }

# The chromosome column must exist, or every gene is 'other' and nothing runs.
head -1 "${TABLE}" | tr '\t' '\n' | grep -q -x -F chrom1 || {
  echo "ERROR: no chrom1 column in ${TABLE}; is this annotate_regions.py output?" >&2
  exit 1; }

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/x_autosome_check_${STAMP}"
printf '%-8s %s\n' "Table:" "${TABLE}" "Output:" "${OUTDIR}"

"${PY}" "${SCRIPT}" --annotated "${TABLE}" --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/x_autosome_summary.txt"
