#!/bin/bash
#SBATCH -J ab_panel
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/ab_panel_%j.out
#SBATCH -e logs/ab_panel_%j.err

# ---------------------------------------------------------------------------
# Supplementary fig. S14D (stem S16D): what the allele-balance screen removed,
# by region class. allele_balance_panel.py re-applies filter_genotypes.py's
# screen to every biallelic SNP in classified CDS of the PRE-screen callset
# (step3_gapfilter.bcf) and checks each verdict against the screen's own
# excluded_sites.txt.
#
# Usage (submit from 07_polymorphism; the -o/-e paths resolve against the
# submission directory):
#   sbatch run_allele_balance_panel.sh                  # filtered_miss4
#   sbatch run_allele_balance_panel.sh filtered         # another callset
#
# Arg: the callset directory under ${POPGEN}/vcf (default filtered_miss4, the
# callset of the headline).
#
# Output: ${POPGEN}/results/figures/s16d_allele_balance_<YYYYMMDD>_<SET>/
#   allele_balance_screen.tsv  the table figure_cohort_selection.py draws as
#                              S14D (it takes the newest one under
#                              results/figures)
#   S16D_allele_balance.{svg,pdf,png}, provenance.json
#
# One CPU: the work is one bcftools query pipe over ~12.4 M records, with
# Python parsing the ~150,000 variant lines it passes through.
#
# WHAT TO READ IN THE LOG:
#   1. "Check against .../excluded_sites.txt": must end "AGREE". Any
#      disagreement exits 3 (so an afterok job does not run) and prints
#      examples; it means the screen parameters or the files do not match the
#      filter_genotypes.py run.
#   2. The per-class "removed N of M variant sites" lines and the fold. These
#      replace 7.23% / 24.92% (3.4x) in Results, the S14D legend and Methods.
#   3. The not-judged share: sites with fewer than three usable heterozygotes
#      are kept by the screen and stay in the callset. There is no "floor".
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
REQUIRED=(PROJECT_DIR POPGEN_DIR POPGEN_WORK_DIR)
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
SCRIPT="$(find_script allele_balance_panel.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: allele_balance_panel.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

SET="${1:-filtered_miss4}"
POPGEN="${POPGEN_DIR}"
MINIPROT="${PROJECT_DIR}"
export PAPER_COLORS_DIR="${REPO_ROOT}/lib"

VCF="${POPGEN}/vcf/${SET}/step3_gapfilter.bcf"
EXCL="${POPGEN}/vcf/${SET}/excluded_sites.txt"
GFF="${MINIPROT}/As_miniprot.gff"
CALLABLE="${POPGEN}/targets/gene_callable.tsv"
OUTBASE="${POPGEN}/results/figures"

for f in "${VCF}" "${GFF}" "${CALLABLE}" "${PAPER_COLORS_DIR}/paper_colors.py"; do
  if [ ! -s "${f}" ]; then
    echo "ERROR: input missing or empty: ${f}" >&2
    exit 1
  fi
done
# excluded_sites.txt is empty only if the screen dropped nothing, so it is
# required to exist, not to be non-empty.
if [ ! -f "${EXCL}" ]; then
  echo "ERROR: ${EXCL} not found; the panel cannot be checked against the" \
       "screen that was applied." >&2
  exit 1
fi
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH" >&2; exit 1; }

# Record the inputs in the log.
echo "[$(date '+%F %T')] SET=${SET}"
echo "  vcf      ${VCF}"
echo "  excluded ${EXCL} ($(wc -l < "${EXCL}") lines)"
echo "  outdir   ${OUTBASE}/s16d_allele_balance_$(date '+%Y%m%d')_${SET}"

python "${SCRIPT}" \
  --vcf "${VCF}" \
  --excluded-sites "${EXCL}" \
  --gff "${GFF}" \
  --gene-callable "${CALLABLE}" \
  --outdir "${OUTBASE}" \
  --tag "${SET}" \
  --log-level INFO
