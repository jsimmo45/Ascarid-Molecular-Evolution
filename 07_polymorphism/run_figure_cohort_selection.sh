#!/bin/bash
#SBATCH -J s14_cohort
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/s14_cohort_%j.out
#SBATCH -e logs/s14_cohort_%j.err

# ---------------------------------------------------------------------------
# Supplementary fig. S14 (stem S16_cohort_selection), panels A-D, assembled
# by figure_cohort_selection.py from three summary tables. It recomputes
# nothing, so it is quick; one CPU.
#
# Until 2026-10-14 this script ran only through sbatch --wrap, so the
# arguments of the August-October builds were not recorded. These are the
# arguments documented in the script and behind the recorded inputs
# (cohort_20260820/sample_metadata.tsv, bam_matched depth tables, the
# allele-balance table under results/figures).
#
# Usage (submit from 07_polymorphism):
#   sbatch run_figure_cohort_selection.sh             # newest table
#   sbatch run_figure_cohort_selection.sh <path>      # a given table or dir
#   sbatch --dependency=afterok:<ab_panel job> run_figure_cohort_selection.sh
#
# Arg: allele_balance_screen.tsv, or a directory searched for it; default
# ${POPGEN}/results/figures, where the newest copy wins. The log line
# "Using .../allele_balance_screen.tsv" says which one was taken.
#
# Output: ${POPGEN}/results/figures/s16_cohort_selection_<YYYYMMDD>/
# Date only: a second run on the same day overwrites the first. The
# supplementary manifest (05_analyses/supplementary_figure_sources.tsv, row
# S14) names this folder, so point it at the new date before reassembling.
#
# WHAT TO READ IN THE LOG:
#   1. "Using ..." for the allele-balance table: the run you meant.
#   2. "D  the screen removes ..." : must equal the per-class percentages in
#      the allele_balance_panel.py log, and is what Results and the legend
#      quote.
#   3. FIGURE AUDIT -- clean.
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
SCRIPT="$(find_script figure_cohort_selection.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: figure_cohort_selection.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

# The 2026-10-31 script labels the worm in panel B's band "Unassigned", as the
# legend does; an older copy says "Unknown". Refuse it rather than redraw that.
# (SCRIPT is set only where a config.sh header locates the script.)
grep -q '"unknown": "Unassigned"' "${SCRIPT:-figure_cohort_selection.py}" || {
  echo "ERROR: ${SCRIPT:-figure_cohort_selection.py} predates 2026-10-31 (key still \"Unknown\");" >&2
  echo "  copy in the new script first." >&2; exit 1; }

POPGEN="${POPGEN_DIR}"
export PAPER_COLORS_DIR="${REPO_ROOT}/lib"

AB="${1:-${POPGEN}/results/figures}"
META=cohort_20260820
DEPTH="${POPGEN}/bam_matched"
OUTBASE="${POPGEN}/results/figures"

for d in "${META}" "${DEPTH}"; do
  if [ ! -d "${d}" ]; then
    echo "ERROR: directory not found: ${d}" >&2
    exit 1
  fi
done
if [ ! -e "${AB}" ]; then
  echo "ERROR: allele-balance table or directory not found: ${AB}" >&2
  exit 1
fi

# Record the inputs in the log.
echo "[$(date '+%F %T')] metadata=${META} depth=${DEPTH} allele_balance=${AB}"
echo "  outdir ${OUTBASE}/s16_cohort_selection_$(date '+%Y%m%d')"

python "${SCRIPT}" \
  --metadata "${META}" \
  --depth-match-glob "${DEPTH}" \
  --allele-balance "${AB}" \
  --outdir "${OUTBASE}" \
  --log-level INFO
