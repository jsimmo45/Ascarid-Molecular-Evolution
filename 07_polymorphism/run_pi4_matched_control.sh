#!/bin/bash
#SBATCH -J pi4_matched
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/pi4_matched_%j.out
#SBATCH -e logs/pi4_matched_%j.err

# ---------------------------------------------------------------------------
# pi4-matched control: is the eliminated-class pi0/pi4 elevation explained by
# lower neutral diversity (reduced Ne)? Genes are grouped into positionally
# contiguous windows with a floor on four-fold sites, and eliminated windows
# are compared with retained windows of matched pi4. Feeds Fig 7B and 7C.
#
# Usage:
#   sbatch run_pi4_matched_control.sh pi_filtered_miss4 5000  # headline
#   sbatch run_pi4_matched_control.sh pi_clustfilt 5000       # robustness
#
# Args: the pi_by_region output directory under ${POPGEN}/results, then
# --min-s4 (four-fold sites per window; default 5000, the Fig 7 setting).
# Use the same pi_by_region directory as the figure and the sweep.
#
# Output: ${POPGEN}/results/pi4_matched_<min-s4>_<YYYYMMDD_HHMMSS>_<PI_SET>
# Timestamped, so a rerun cannot overwrite the run a figure came from. Pass
# that directory name as the second argument to run_figure_polymorphism.sh.
#
# No wrapper existed before 2026-10-08, so the commands behind the August
# outputs (pi4_matched, pi4_matched_5k, pi4_matched_5k_miss4) are not
# recorded. The bootstrap is seeded (default_rng(0)), so a rerun on the same
# pi_by_region directory reproduces the windows, the matched fold and its
# intervals exactly.
#
# Why rerun pi4_matched_5k_miss4 (2026-08-21): it was written by an earlier
# version of pi4_matched_control.py whose trend statistic divided the
# eliminated excess by the scatter of single retained windows (z = 1.54 in
# Results). The current script tests the pooled eliminated ratio against its
# own standard error plus that of the prediction. Its summary has the columns
# se_eliminated, se_prediction, p_vs_retained_trend and z_vs_window_scatter;
# the 08-21 file lacks them.
#
# WHAT TO READ IN THE LOG:
#   1. TEST 1: the matched fold, its interval, and how many retained windows
#      fall inside the eliminated pi4 range. Must equal the 08-21 run (13
#      retained windows, 0.3402 against 0.4892 on 9 eliminated windows, fold
#      1.438). Quote the window count wherever the fold is quoted.
#   2. TEST 2: "Difference = ..., z SE above the trend, P = ..." replaces
#      "z = 1.54; P = 0.12" in Results. The predicted value must equal the
#      08-21 run (0.3728). If P < 0.05 the "Decomposition" line gives the
#      Ne share against the trend at the retained pooled pi4; it replaces
#      the ~43% in Results. If P >= 0.05 no share is printed.
#   3. Agreement with the sweep: the fold must equal the naive row at the
#      same --min-s4 in pi4_sweep_<date>_<PI_SET> (1.438 at 5,000 in
#      pi4_sweep_20260822_pi_filtered_miss4). A mismatch means the two
#      scripts' windowing has diverged.
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
SCRIPT="$(find_script pi4_matched_control.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: pi4_matched_control.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

PI_SET="${1:-pi_filtered_miss4}"
MIN_S4="${2:-5000}"
POPGEN="${POPGEN_DIR}"
PER_GENE="${POPGEN}/results/${PI_SET}/per_gene.tsv"
CALLABLE="${POPGEN}/targets/gene_callable.tsv"

# 5000 -> "5k", matching the existing pi4_matched_5k naming; a value that is
# not a whole number of thousands is written as given.
S4_LABEL=$(awk -v s="${MIN_S4}" 'BEGIN { if (s % 1000 == 0) printf "%dk", s / 1000; else printf "%s", s }')
STAMP=$(date '+%Y%m%d_%H%M%S')
OUTDIR="${POPGEN}/results/pi4_matched_${S4_LABEL}_${STAMP}_${PI_SET}"

for f in "${PER_GENE}" "${CALLABLE}"; do
  if [ ! -s "${f}" ]; then
    echo "ERROR: input missing or empty: ${f}" >&2
    exit 1
  fi
done
if [ -e "${OUTDIR}" ]; then
  echo "ERROR: ${OUTDIR} already exists." >&2
  exit 1
fi

# Record the inputs in the log.
echo "[$(date '+%F %T')] PI_SET=${PI_SET} MIN_S4=${MIN_S4} OUTDIR=${OUTDIR}"

python "${SCRIPT}" \
  --per-gene "${PER_GENE}" \
  --gene-callable "${CALLABLE}" \
  --min-s4 "${MIN_S4}" \
  --n-boot 2000 \
  --outdir "${OUTDIR}" \
  --log-level INFO
