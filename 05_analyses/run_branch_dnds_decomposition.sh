#!/bin/bash
#SBATCH -J fig5_grid
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=fig5_grid_%j.out
#SBATCH --error=fig5_grid_%j.err
#
# Figure 5 and supplementary fig. S8, each WHOLE, from branch_dnds_decomposition.py:
# three rows of four panels (A branch dN/dS, B branch dN, C branch dS; one
# column per terminal branch), built at 185 mm with every axes placed in mm.
#   fig5_branch_grid_as_region.{svg,pdf,png}     Figure 5 (Ascaris labels)
#   figS14_branch_grid_tox_region.{svg,pdf,png}  S8 (T. canis labels; S14 is
#                                                S8's frozen-map number)
# Replaces the hand assembly of 5A (branch_dnds_four.py) over 5C, which needed
# 5C scaled to 5A's width (2026-10-22; results draft note 78). Statistics and
# dnds_decomposition_stats.tsv are unchanged.
#
# Reads the existing four-way codeml table -- no recomputation, so single-threaded
# pandas/scipy/matplotlib on 1 node / 1 task / 1 CPU, about a minute.
#
#   sbatch run_branch_dnds_decomposition.sh
#   sbatch run_branch_dnds_decomposition.sh <path/to/branch_dnds_four_annotated.tsv>
#   FIGURE=both sbatch run_branch_dnds_decomposition.sh     # + the old 5C/5D files
#   DECADE_MM=6.5 sbatch run_branch_dnds_decomposition.sh   # shorter figure
#
# Input: argument 1, else ANNOTATED_FOUR from config.sh, else the newest
# branch_dnds_four_*/branch_dnds_four_annotated.tsv under RESULTS_DIR. The table
# behind the current Figure 5 is
#   results/branch_dnds_four_20261010_113902/branch_dnds_four_annotated.tsv
# (job 6364785, unrooted tree); pass it explicitly to be certain.
#
# 2026-10-22: job name fig5_grid (dnds_decomp was shared with
# run_dnds_decomposition.sh, Figure 2A); config.sh and the .py are found from
# SLURM_SUBMIT_DIR and PWD as well as the script's own location, which under
# sbatch is SLURM's spool copy; the conda environment defaults to kaks with its
# bin/ first on PATH; refuses a branch_dnds_decomposition.py without the grid.
# Later on 2026-10-22: refuses the first grid version as well (no
# _jitter_sample), whose sample of points could leave a band that holds one
# gene without its point, beside a printed "1".

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
SCRIPT="$(find_script branch_dnds_decomposition.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: branch_dnds_decomposition.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# Only the latest 2026-10-22 version draws the grid with every band shown.
grep -q '^def _jitter_sample' "${SCRIPT}" || {
  if grep -q '^def plot_branch_grid' "${SCRIPT}"; then
    echo "ERROR: ${SCRIPT} is the first 2026-10-22 version (no _jitter_sample)." >&2
  else
    echo "ERROR: ${SCRIPT} is the version from before 2026-10-22 (no plot_branch_grid)." >&2
  fi
  echo "  Copy the new branch_dnds_decomposition.py into 05_analyses/." >&2
  exit 1; }

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

# Probe the imports: branch_lib, region_lib and paper_colors live in lib/.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("branch_lib", "region_lib", "paper_colors", "scipy",
                       "matplotlib")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("modules: branch_lib, region_lib, paper_colors, scipy, matplotlib all importable")
PYEOF

cd "${PROJECT_DIR}"

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${OUTDIR:-${RESULTS_DIR}/dnds_decomposition_${STAMP}}"
OMEGA_CAP="${OMEGA_CAP:-5}"
STAR_STAT="${STAR_STAT:-q}"
FDR_FAMILY="${FDR_FAMILY:-metric}"
FIGURE="${FIGURE:-grid}"
DECADE_MM="${DECADE_MM:-7}"
SEED="${SEED:-0}"

ANNOTATED="${1:-${ANNOTATED_FOUR:-}}"
if [[ -z "${ANNOTATED}" ]]; then
  ANNOTATED="$(ls -1dt "${RESULTS_DIR}"/branch_dnds_four_*/branch_dnds_four_annotated.tsv \
               2>/dev/null | head -n1 || true)"
fi
if [[ -z "${ANNOTATED}" || ! -f "${ANNOTATED}" ]]; then
  echo "ERROR: no annotated table. Pass one as argument 1, set ANNOTATED_FOUR" >&2
  echo "       in config.sh, or run run_branch_dnds_four.sh first." >&2
  exit 1
fi

N_GENES=$(( $(wc -l < "${ANNOTATED}") - 1 ))
echo "input : ${ANNOTATED}"
echo "genes : ${N_GENES}"
PROV="$(dirname "${ANNOTATED}")/provenance.json"
if [[ -f "${PROV}" ]]; then
  "${PY}" - "${PROV}" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
print("run   : %s  commit %s  quartets %s -> ok %s"
      % (d.get("run_utc", "?"), d.get("git_commit", "?"),
         d.get("n_quartets", "?"), d.get("n_ok", "?")))
EOF
else
  echo "run   : no provenance.json beside this table -- check that it is the"
  echo "        unrooted-tree table (load_annotated refuses a rooted one)."
fi
echo "output: ${OUTDIR}"

# --mask-ds is deliberately NOT passed: the dN/dS cap is a filter on dS, so
# testing dS after applying it is circular. --star-stat q --fdr-family metric
# are the defaults, written out so the log shows them: BH within each metric,
# across the four branches and both labellings (8 tests per row).
"${PY}" "${SCRIPT}" \
  --annotated "${ANNOTATED}" \
  --omega-cap "${OMEGA_CAP}" \
  --star-stat "${STAR_STAT}" \
  --fdr-family "${FDR_FAMILY}" \
  --figure "${FIGURE}" \
  --decade-mm "${DECADE_MM}" \
  --seed "${SEED}" \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
ls -1 "${OUTDIR}"
echo ""
sed -n '/Figure legend notes: the grid/,/^Grid layout/p' "${OUTDIR}/decomposition_summary.txt"
echo ""
echo "Manifests: point Figure 5 and S8 at the grid files, one row each, in"
echo "  main_figure_sources.tsv:          5   A-C_branch_grid   ${OUTDIR}/fig5_branch_grid_as_region.svg"
echo "  supplementary_figure_sources.tsv: S8  S14  A-C_branch_grid  ${OUTDIR}/figS14_branch_grid_tox_region.svg"
echo "replacing the four 5A-5D and four S8 A-D rows."
