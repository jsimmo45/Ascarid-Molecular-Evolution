#!/bin/bash
#SBATCH -J clock_test
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=clock_test_%j.out
#SBATCH --error=clock_test_%j.err
#
# ---------------------------------------------------------------------------
# run_molecular_clock_test.sh
#
# Global-clock LRT (PAML baseml) + Tajima relative-rate tests on the four-way
# quartet alignments. Formalises the rate heterogeneity already visible in the
# branch lengths: Ascaris 0.0461 vs Parascaris 0.1496 on branches that subtend
# the same elapsed time.
#
#   sbatch run_molecular_clock_test.sh                 # sites 12, 3 and all
#   sbatch run_molecular_clock_test.sh "12"            # one subset only
#   bash   run_molecular_clock_test.sh inspect         # header check, no job
#
# baseml is single-threaded, hence 1 node / 1 task / 1 CPU. On ~6,000 genes it
# runs in minutes; the 4 h wall time is slack, not an estimate.
# ---------------------------------------------------------------------------

set -euo pipefail

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
REQUIRED=(PROJECT_DIR)
for v in "${REQUIRED[@]}"; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: ${v} is not set in ${CONFIG}" >&2; exit 1; }
done

find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/03_rates" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script molecular_clock_test.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: molecular_clock_test.py not found near ${REPO_ROOT}/03_rates." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV:?set CONDA_ENV in config.sh}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
fi

cd "${PROJECT_DIR}" || exit 1

MODE="${1:-12 3 all}"

# --- locate the alignments --------------------------------------------------
# CRITICAL: four_way_alignments/ is dated 07-31, BEFORE the Baylisascaris
# identifier fix landed on 08-08, and its provenance was never resolved (see
# claude/kaks_analysis_status.md, "Unresolved provenance"). If it was built
# pre-fix it holds ~2,771 quartets instead of ~5,954, and every number below
# would be computed on half the Baylisascaris genes -- exactly the failure that
# cost two days on run_branch_dnds_decomposition.sh.
#
# Prefer any alignment directory written by the v2 run, and refuse to proceed
# on a suspiciously small set.
MIN_ALIGNMENTS=5000

ALN_DIR=""
for cand in \
    branch_dnds_four_results_v2/alignments \
    branch_dnds_four_results_v2/aln \
    four_way_alignments_v2 \
    four_way_alignments ; do
  if [[ -d "$cand" ]]; then ALN_DIR="$cand"; break; fi
done

if [[ -z "$ALN_DIR" ]]; then
  echo "ERROR: no alignment directory found. Looked for:" >&2
  echo "  branch_dnds_four_results_v2/{alignments,aln}" >&2
  echo "  four_way_alignments_v2 , four_way_alignments" >&2
  echo "Find them with:  find . -name '*.fa' -path '*align*' | head" >&2
  exit 1
fi

# Alignments may be .fa, .fasta or .aln depending on which script wrote them.
GLOB='*.fa'
for g in '*.fa' '*.fasta' '*.aln' '*.fna'; do
  if compgen -G "$ALN_DIR/$g" > /dev/null; then GLOB="$g"; break; fi
done

N_ALN=$(find "$ALN_DIR" -maxdepth 1 -name "$GLOB" | wc -l)
echo "Alignment directory: $ALN_DIR"
echo "Pattern:             $GLOB"
echo "Files:               $N_ALN"

if [[ "$N_ALN" -lt "$MIN_ALIGNMENTS" ]]; then
  echo "" >&2
  echo "ERROR: only $N_ALN alignments, expected >= $MIN_ALIGNMENTS." >&2
  echo "This is almost certainly a PRE-FIX set built before the" >&2
  echo "Baylisascaris identifier bug was corrected on 2026-08-08" >&2
  echo "(pre-fix quartet count was 2771; v2 is 5954)." >&2
  echo "" >&2
  echo "Rebuild the alignments from rbh_*_v2/orthologs.tsv, or point" >&2
  echo "ALN_DIR at the directory the v2 branch_dnds run wrote." >&2
  exit 1
fi

# --- header check only ------------------------------------------------------
# Sequence headers are gene IDs, not species names. The script resolves
# Toxocara / Parascaris / Baylisascaris by ID pattern and Ascaris by
# elimination. Confirm that before trusting any result: a silently
# mis-assigned species would invalidate everything.
if [[ "$MODE" == "inspect" ]]; then
  python "${SCRIPT}" \
    --aln-dir "$ALN_DIR" --aln-glob "$GLOB" --inspect
  exit 0
fi

# --- run --------------------------------------------------------------------
# Three site subsets, because the obvious objection to a clock rejection is
# that Toxocara's third positions are saturated (median Ks 1.39-1.43):
#   12   codon positions 1+2 -- the headline result, saturation-free
#   3    third positions only -- where the synonymous signal lives
#   all  everything, for completeness
# If all three reject, saturation is not the explanation.
for SITES in $MODE; do
  echo ""
  echo "===== codon positions: $SITES ====="
  python "${SCRIPT}" \
    --aln-dir "$ALN_DIR" \
    --aln-glob "$GLOB" \
    --sites "$SITES" \
    --topology '((Ascaris,Parascaris),(Baylisascaris,Toxocara));' \
    --model 7 \
    --ncatG 5 \
    --outdir "clock_test_sites${SITES}"
done

echo ""
echo "===== DONE ====="
for SITES in $MODE; do
  echo ""
  echo "--- clock_test_sites${SITES} ---"
  sed -n '/TEST 1/,/^$/p;/TEST 2/,/^$/p' \
    "clock_test_sites${SITES}/clock_test_summary.txt" 2>/dev/null
done
