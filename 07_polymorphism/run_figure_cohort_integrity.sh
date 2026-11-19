#!/bin/bash
#SBATCH -J s19_cohort
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/s19_%j.out
#SBATCH -e logs/s19_%j.err

# ---------------------------------------------------------------------------
# Supplementary fig. S17 — cohort integrity, assembled from three runs. (S19
# before the 2026-10-12 renumbering; the output keeps the S19 stem, so the
# folder is results/figures/s19_cohort_integrity_<date>_<set>.)
#
#   sbatch run_figure_cohort_integrity.sh filtered_miss4
#
# Since 2026-10-31 the figure is built at 185 mm, the width of every other
# figure, with the left column sized so panel A's matrix fills it.
#
# Reads three summary TSVs and DRAWS. It recomputes nothing, so it cannot
# disagree with the scripts that produced them. Those three must have run
# first:
#
#   run_check_duplicate_individuals.sh   -> pairwise.tsv
#   run_het_by_region.sh                 -> het_by_region.tsv
#   run_jackknife_individuals.sh         -> jackknife_summary.tsv
#
# Seconds, 4G, one CPU. Fine on a login node.
#
# BEFORE DRAWING it checks that all three tables name the same worms and exits
# if they do not. Mixing a stale run into the assembly is the obvious failure
# mode here and it would be completely invisible in the finished figure.
#
# WHAT TO READ IN THE LOG:
#   1. "Cohort agrees across all three inputs".
#   2. Panel C's count of worms outside the +/-25% band. C is the panel that
#      carries B: a 2.5x spread in absolute heterozygosity looks alarming on
#      its own, and what decides whether it can bias a BETWEEN-class contrast
#      is whether the RATIO moves.
#   3. Whether panel D reports its x-axis as zoomed. When every jackknife value
#      sits far above 1, the no-difference line is dropped rather than
#      compressing all the points into a tenth of the panel — and the caption
#      then has to say the axis does not include 1.
#   4. The "Layout:" line (A's matrix in mm, residual under 0.2 mm) and no
#      "Off the canvas" warning.
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
SCRIPT="$(find_script figure_cohort_integrity.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: figure_cohort_integrity.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

# An older figure_cohort_integrity.py draws the equal-column layout that left
# a blank strip down the right side. Refuse it rather than redraw that.
# (SCRIPT is set only where a config.sh header locates the script; here it is
# the copy in this directory.)
grep -q '^def fit_left_column' "${SCRIPT:-figure_cohort_integrity.py}" || {
  echo "ERROR: ${SCRIPT:-figure_cohort_integrity.py} predates 2026-10-31 (no fit_left_column);" >&2
  echo "  copy in the new script first." >&2; exit 1; }

export PAPER_COLORS_DIR="${REPO_ROOT}/lib"

SET="${1:-filtered_miss4}"
POPGEN="${POPGEN_DIR}"

python "${SCRIPT}" \
  --pairwise  "${POPGEN}/results/duplicates_${SET}/pairwise.tsv" \
  --het       "${POPGEN}/results/het_by_region_${SET}/het_by_region.tsv" \
  --jackknife "${POPGEN}/results/jackknife_${SET}/jackknife_summary.tsv" \
  --dup-discordance 0.05 \
  --ratio-tolerance 0.25 \
  --tag "${SET}" \
  --outdir "${POPGEN}/results/figures"
