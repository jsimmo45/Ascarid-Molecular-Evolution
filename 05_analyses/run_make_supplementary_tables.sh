#!/bin/bash
#SBATCH -J supp_tables
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=supp_tables_%j.out
#SBATCH --error=supp_tables_%j.err
#
# Supplementary Tables S1-S11: one workbook (Supplementary_Tables.xlsx), one TSV
# per table part, the captions, and a summary of every cross-table check.
#
#   sbatch run_make_supplementary_tables.sh
#   sbatch run_make_supplementary_tables.sh --only S6,S11
#   sbatch run_make_supplementary_tables.sh --source x_dosage=results/<newer run>/per_chromosome.tsv
#
# Inputs: supplementary_table_sources.tsv beside make_supplementary_tables.py
# (key<TAB>path; relative paths resolve against PROJECT_DIR from config.sh, or
# SUPP_TABLE_BASE if set). Arguments are passed straight to the script.
#
# Single-threaded pandas/scipy/openpyxl on ~40 small tables and three of ~10,000
# rows; under a minute. 1 node / 1 task / 1 CPU. 4 GB covers openpyxl holding
# the whole workbook (~600,000 cells) in memory.
#
# Exit status: 0 all checks passed; 2 tables written but a check failed (read
# the summary printed below); 1 inputs missing, nothing written.
#
# 2026-11-01: table S9 has a part C, the insulation result at five Hi-C stages
# (testis, ovary, zygote, 48 h, 4-cell (60 h)), read from the stage_* lines of
# the manifest; the testis lines now point at job 6436588 (2026-10-19), whose
# models equal the 2026-10-17 run's. This wrapper refuses an older script,
# which would ignore those lines and build S9 without part C.
#
# 2026-11-05: table S11 has a row for positive selection, built from the
# McDonald-Kreitman summary of pi_by_region.py (manifest key mk_summary:
# ${POPGEN_DIR}/results/pi_filtered_miss4/mk_summary.tsv), with DoS and alpha at
# MAF floors of 0.10, 0.15 and 0.20. The wrapper refuses a script without it.

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
SCRIPT="$(find_script make_supplementary_tables.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: make_supplementary_tables.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
grep -q '^STAGE_RUNS = \[' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} predates 2026-11-01 (no S9 part C); copy in the new script first." >&2
  exit 1; }
grep -q '^MK_FLOORS = ' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} predates 2026-11-05 (no McDonald-Kreitman row in S11); copy in the new script first." >&2
  exit 1; }
MANIFEST="${SUPP_TABLE_SOURCES:-$(dirname "${SCRIPT}")/supplementary_table_sources.tsv}"
[[ -s "${MANIFEST}" ]] || { echo "ERROR: missing or empty manifest: ${MANIFEST}" >&2; exit 1; }
echo "manifest: ${MANIFEST}"

export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

# The manifest expands ${REPO_ROOT} (the curated S11 file) and ${POPGEN_DIR}
# (the population-sample audit and the McDonald-Kreitman summary of S11, outside
# the project tree). config.sh must set POPGEN_DIR.
export REPO_ROOT
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

# Probe the imports: region_lib and rate_tables live in lib/.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) or the environment, and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("region_lib", "rate_tables", "scipy", "openpyxl")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("modules: region_lib, rate_tables, scipy, openpyxl all importable")
PYEOF

cd "${PROJECT_DIR}"

BASE="${SUPP_TABLE_BASE:-${PROJECT_DIR}}"
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/supplementary_tables_${STAMP}"
printf '%-9s %s\n' "Base:" "${BASE}" "Output:" "${OUTDIR}"

rc=0
"${PY}" "${SCRIPT}" --manifest "${MANIFEST}" --base-dir "${BASE}" \
  --outdir "${OUTDIR}" "$@" || rc=$?

echo ""
if [[ -s "${OUTDIR}/supplementary_tables_summary.txt" ]]; then
  echo "===== SUMMARY (exit ${rc}) ====="
  cat "${OUTDIR}/supplementary_tables_summary.txt"
else
  echo "===== FAILED before writing (exit ${rc}); see the .err file ====="
fi
exit "${rc}"
