#!/bin/bash
#SBATCH -J pis_position
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/pis_position_%j.out
#SBATCH -e logs/pis_position_%j.err

# ---------------------------------------------------------------------------
# Neighbourhood and linked-selection controls (pis_vs_position.py):
#   - piS and piN/piS for retained genes within vs beyond --proximal-kb of a
#     CBR, and for eliminated genes;
#   - eliminated genes against neighbourhood-matched (CBR-proximal) retained;
#   - Spearman(piS, CBR distance) per class, Spearman(piS, tract size).
# Feeds the "Neighbourhood" and "Linked selection" sentences in Results.
#
# Usage (submit from 07_polymorphism, so logs/ resolves here):
#   sbatch run_pis_vs_position.sh pi_filtered_miss4   # primary (headline)
#   sbatch run_pis_vs_position.sh pi_clustfilt        # robustness only
#
# Arg: the pi_by_region output directory under ${POPGEN}/results. Use the
# same one as Fig 7 (pi_filtered_miss4 since 2026-10-08).
#
# The other inputs are fixed to the ones every Figure 4 and Fig 7 number
# uses:
#   --kaks  ANNOTATED_TABLE in config.sh; for the paper,
#           results/annotate_regions_20260817/kaks_annotated.tsv (corrected
#           CBR set: abutting tracts merged, terminal boundaries dropped)
#   --gene-callable   ${POPGEN}/targets/gene_callable.tsv
#   --eliminated-bed  ${V50}/AG.v50.eliminated_strict.bed
#
# Output: ${POPGEN}/results/pis_position_<YYYYMMDD_HHMMSS>_<PI_SET>
# Timestamped, so a rerun cannot overwrite an earlier one.
#
# Why this wrapper exists: results/pis_position was written 2026-08-20 17:10
# without one. That is before both current per_gene.tsv files (08-21 10:04
# and 10:05) existed, and the run recorded no inputs, so the callset and the
# CBR table behind it cannot be recovered.
#
# WHAT TO READ:
#   1. pis_position_summary.tsv: piS and piN/piS for "retained, <= 150 kb",
#      "retained, > 150 kb" and "eliminated (all)", and the Spearman rows.
#      Results quotes the near/far piS ratio and the tract-size rho and P.
#   2. The two 'contrast' rows at the end of the same table (also in the .err
#      log under "Contrasts"): the near-minus-far piN/piS difference and the
#      eliminated / neighbourhood-matched fold, each with its gene-bootstrap
#      95% interval (estimate, ci_lo, ci_hi) and P(boot > null). These are
#      the two intervals in the Results "Neighbourhood" sentence; the
#      script computes them since 2026-10-09.
#   3. Gene counts. The August run used 6,882 retained / 264 eliminated (the
#      genes carrying both a polymorphism estimate and a Ka/Ks row). A
#      different count means a different gene set, not only a different
#      callset.
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
REQUIRED=(PROJECT_DIR POPGEN_DIR POPGEN_WORK_DIR AS_V50_DIR ANNOTATED_TABLE)
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
SCRIPT="$(find_script pis_vs_position.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: pis_vs_position.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
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
POPGEN="${POPGEN_DIR}"
V50="${AS_V50_DIR}"

PER_GENE="${POPGEN}/results/${PI_SET}/per_gene.tsv"
# The annotated rate table (ANNOTATED_TABLE in config.sh, relative to PROJECT_DIR).
case "${ANNOTATED_TABLE}" in
  /*) KAKS_TABLE="${ANNOTATED_TABLE}" ;;
  *)  KAKS_TABLE="${PROJECT_DIR}/${ANNOTATED_TABLE}" ;;
esac
KAKS="${KAKS_TABLE}"
CALLABLE="${POPGEN}/targets/gene_callable.tsv"
ELIM_BED="${V50}/AG.v50.eliminated_strict.bed"

STAMP=$(date '+%Y%m%d_%H%M%S')
OUTDIR="${POPGEN}/results/pis_position_${STAMP}_${PI_SET}"

for f in "${PER_GENE}" "${KAKS}" "${CALLABLE}" "${ELIM_BED}"; do
  if [ ! -s "${f}" ]; then
    echo "ERROR: input missing or empty: ${f}" >&2
    exit 1
  fi
done
if [ -e "${OUTDIR}" ]; then
  echo "ERROR: ${OUTDIR} already exists." >&2
  exit 1
fi

# Record every input in the log; the August run did not, which is why this
# rerun is needed.
echo "[$(date '+%F %T')] PI_SET=${PI_SET}"
echo "  per-gene        ${PER_GENE}"
echo "  kaks            ${KAKS}"
echo "  gene-callable   ${CALLABLE}"
echo "  eliminated-bed  ${ELIM_BED}"
echo "  outdir          ${OUTDIR}"

python "${SCRIPT}" \
  --per-gene "${PER_GENE}" \
  --kaks "${KAKS}" \
  --gene-callable "${CALLABLE}" \
  --eliminated-bed "${ELIM_BED}" \
  --proximal-kb 150 \
  --n-boot 2000 \
  --seed 0 \
  --outdir "${OUTDIR}" \
  --log-level INFO
