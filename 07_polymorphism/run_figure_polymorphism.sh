#!/bin/bash
#SBATCH -J figpoly
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/figpoly_%j.out
#SBATCH -e logs/figpoly_%j.err

# ---------------------------------------------------------------------------
# Figure 6 (Figure 7 until the 2026-10-26 renumbering; figure_polymorphism.py),
# rebuilt from the inputs a previous build
# recorded in its provenance.json (--inputs-from). It draws and recomputes
# only panel A's bootstrap, which is seeded, so the same inputs reproduce the
# panel A folds and intervals exactly; the log says whether they did.
# Single-threaded pandas/matplotlib on summary tables: one CPU.
#
# No wrapper existed before 2026-10-14 (the August-October builds ran through
# sbatch --wrap), which is why the inputs are taken from the provenance.json
# of the figure in the manuscript rather than retyped here.
#
# Usage (submit from 07_polymorphism):
#   sbatch run_figure_polymorphism.sh
#   sbatch run_figure_polymorphism.sh <previous run dir> [<pi4_matched dir>]
#
# Arg 1: the previous figure_polymorphism run directory; default
#   ${POPGEN}/results/figures/figure_polymorphism_20261015_pi_filtered_miss4,
#   the build in the manuscript until the 2026-10-29 bar-order rerun;
#   --inputs-from reuses exactly the inputs that build recorded.
# Arg 2 (optional): a pi4_matched_control output directory, under
#   ${POPGEN}/results or absolute, whose pi4_matched_summary.tsv and
#   windows.tsv replace the recorded ones, e.g.
#   pi4_matched_5k_20261014_111442_pi_filtered_miss4 (same windows and values
#   as the 08-21 tables the 10-08 build used; results draft note 51).
#
# Output: ${POPGEN}/results/figures/figure_polymorphism_<YYYYMMDD>_pi_filtered_miss4/
# Date only: a second run on the same day overwrites the first. The main
# manifest (Ascarid-Molecular-Evolution/05_analyses/main_figure_sources.tsv,
# the five Fig 6 lines) names this folder.
#
# WHAT TO READ IN THE LOG:
#   1. The "(from the previous run)" lines: the inputs actually used.
#   2. "Panel A reproduces the previous run" -- must say so. "DIFFERS" means
#      an input changed since 10-08: stop and find out which.
#   3. Panel B values, in the drawn order since 2026-10-29: retained genes
#      0.2859, eliminated genes 0.5196, retained pi4-matched windows 0.3402,
#      eliminated windows 0.4892; folds 1.817 (genes) and 1.438 (windows),
#      the 1.82 and 1.44 of Results.
#   4. "FIGURE AUDIT -- clean" for each of the five figures.
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
SCRIPT="$(find_script figure_polymorphism.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: figure_polymorphism.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

POPGEN="${POPGEN_DIR}"
export PAPER_COLORS_DIR="${REPO_ROOT}/lib"

PREV="${1:-${POPGEN}/results/figures/figure_polymorphism_20261015_pi_filtered_miss4}"
PI4="${2:-}"
TAG=pi_filtered_miss4
OUTBASE="${POPGEN}/results/figures"

if [ ! -s "${PREV}/provenance.json" ]; then
  echo "ERROR: no provenance.json in ${PREV}." >&2
  echo "  Name the inputs instead: python figure_polymorphism.py --help" >&2
  exit 1
fi

EXTRA=()
if [ -n "${PI4}" ]; then
  case "${PI4}" in
    /*) ;;
    *) PI4="${POPGEN}/results/${PI4}" ;;
  esac
  for f in pi4_matched_summary.tsv windows.tsv; do
    if [ ! -s "${PI4}/${f}" ]; then
      echo "ERROR: input missing or empty: ${PI4}/${f}" >&2
      exit 1
    fi
  done
  EXTRA=(--matched-summary "${PI4}/pi4_matched_summary.tsv"
         --windows "${PI4}/windows.tsv")
fi

# Record the inputs in the log.
echo "[$(date '+%F %T')] inputs from ${PREV}"
[ -n "${PI4}" ] && echo "  pi4-matched tables from ${PI4}"
echo "  outdir ${OUTBASE}/figure_polymorphism_$(date '+%Y%m%d')_${TAG}"

# Refuse the pre-2026-10-29 script, which draws panel B in the nested order.
if ! grep -q "^PAIR_GAP = " figure_polymorphism.py; then
  echo "ERROR: figure_polymorphism.py in $(pwd) predates the 2026-10-29" \
       "panel B order (no PAIR_GAP). Copy the new version here first." >&2
  exit 1
fi

python "${SCRIPT}" \
  --inputs-from "${PREV}" \
  "${EXTRA[@]}" \
  --outdir "${OUTBASE}" \
  --tag "${TAG}" \
  --log-level INFO
