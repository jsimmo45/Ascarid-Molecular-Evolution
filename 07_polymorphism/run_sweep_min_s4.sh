#!/bin/bash
#SBATCH -J s4sweep
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/s4sweep_%j.out
#SBATCH -e logs/s4sweep_%j.err

# ---------------------------------------------------------------------------
# How small can the pi4-matched windows get before the match stops meaning
# anything?
#
#   sbatch run_sweep_min_s4.sh pi_filtered_miss4
#
# The current --min-s4 5000 leaves only 8 retained windows in the eliminated
# pi4 range — the weakest number in the figure. Lowering --min-s4 buys more
# windows trivially, which is exactly why this needs a guard rather than a
# plain sweep.
#
# Reads TSVs only, no VCF and no genome, so it is fast and small: 8G and one
# CPU. Runs fine on a login node too.
#
# WHAT TO READ IN THE LOG, in this order:
#
#   1. "Self-test passed" — the window assignment here must reproduce
#      pi4_matched_control.make_windows exactly, or the sweep is describing a
#      different estimator than the paper uses.
#
#   2. The rho vs SH rho columns. rho is the naive Spearman, which shares SNPs
#      between pi4 and pi0/pi4 and so is negative from sampling noise alone.
#      SH rho estimates pi4 and pi0/pi4 from disjoint halves of each window's
#      genes, so it cannot be. Where the two track each other the windows hold
#      a real regional signal; where SH rho collapses toward zero while rho
#      stays strong, the windows are too small and the naive number is noise.
#
#   3. The FOLD DRIFT block. Shrinking windows also widens the matched pi4
#      range, which inflates the fold by a route the split-half check does NOT
#      cover. If the fold rises as windows shrink, quote a large-window value.
#
# Do NOT pick the setting with the most matched windows. Smaller windows always
# produce more overlap; the extra windows are the artifact.
#
# OUTPUT MOVED: results now land in
#   results/pi4_sweep_<YYYYMMDD>_<PI_SET>/
# and the figure is S18_min_s4_sweep.{svg,pdf,png}, built at 185 x 62 mm.
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
SCRIPT="$(find_script sweep_min_s4.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: sweep_min_s4.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

export PAPER_COLORS_DIR="${REPO_ROOT}/lib"

PI_SET="${1:-pi_filtered_miss4}"
POPGEN="${POPGEN_DIR}"

python "${SCRIPT}" \
  --per-gene "${POPGEN}/results/${PI_SET}/per_gene.tsv" \
  --gene-callable "${POPGEN}/targets/gene_callable.tsv" \
  --values 1500,2000,3000,5000,7500,10000,15000,20000 \
  --n-boot 2000 \
  --n-splits 25 \
  --gap-tolerance 0.15 \
  --tag "${PI_SET}" \
  --outdir "${POPGEN}/results"
