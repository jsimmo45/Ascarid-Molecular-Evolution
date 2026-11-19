#!/bin/bash
#SBATCH -J build_quartet_aln
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=build_quartet_aln_%j.out
#SBATCH --error=build_quartet_aln_%j.err
#
# ---------------------------------------------------------------------------
# run_build_quartet_alignments.sh
#
# Regenerate the four-way codon alignments that branch_dnds_four.py discarded.
#
#   bash   run_build_quartet_alignments.sh inspect    # column/ID check, no job
#   sbatch run_build_quartet_alignments.sh            # build
#
# WHY: branch_dnds_four.py aligns into tempfile.mkdtemp(prefix='branch4_') and
# deletes it unless --keep-tmp, which the v2 run (SLURM 6089193) did not pass.
# The only alignments on disk are four_way_alignments/ -- 2,772 files dated
# 07-31, built BEFORE the Baylisascaris identifier fix, holding roughly half
# the Baylisascaris genes.
#
# This rebuilds them from the v2 quartet table WITHOUT rerunning codeml, which
# is where the original 1h32m went. Alignment alone is a few minutes: 200
# quartets took 28 s on 8 workers in testing, so ~6,000 on 24 workers is
# roughly 5 minutes. The 2 h wall time is slack.
#
# 24 CPUs because the work is one MAFFT process per quartet and the pool is
# sized to --threads; nothing here benefits from more than one node.
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
SCRIPT="$(find_script build_quartet_alignments.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: build_quartet_alignments.py not found near ${REPO_ROOT}/03_rates." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV:?set CONDA_ENV in config.sh}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
fi

cd "${PROJECT_DIR}" || exit 1

MODE="${1:-build}"

QUARTETS=branch_dnds_four_results_v2/branch_dnds_four_annotated.tsv
OUTDIR=four_way_alignments_v2
THREADS=24

# --- refuse the pre-fix table ----------------------------------------------
# Same failure mode that cost two days on run_branch_dnds_decomposition.sh:
# the pre-fix branch_dnds_four_results/ table has 2,771 quartets, the v2 one
# has 5,954. Checked here as well as inside the script so a mistyped path is
# caught before conda even matters.
if [[ ! -s "$QUARTETS" ]]; then
  echo "ERROR: quartet table not found: $QUARTETS" >&2
  echo "Do NOT substitute branch_dnds_four_results/ -- that is the pre-fix" >&2
  echo "table with 2,771 quartets." >&2
  exit 1
fi

EXPECTED=5954
N_ROWS=$(( $(wc -l < "$QUARTETS") - 1 ))
echo "Quartet table: $QUARTETS"
echo "Rows:          $N_ROWS (expected $EXPECTED)"
if [[ "$N_ROWS" -ne "$EXPECTED" ]]; then
  echo "" >&2
  echo "ERROR: expected $EXPECTED quartets, found $N_ROWS." >&2
  echo "2771 means this is the PRE-FIX table. Point QUARTETS at the _v2" >&2
  echo "results directory." >&2
  exit 1
fi

# --- column and ID check only ----------------------------------------------
# Confirms which table column holds each species' gene ID and what fraction of
# those IDs are actually present in the corresponding FASTA. A low hit rate
# means the wrong FASTA or a different ID convention, and would otherwise show
# up only as a pile of silently dropped quartets.
if [[ "$MODE" == "inspect" ]]; then
  python "${SCRIPT}" --quartets "$QUARTETS" --inspect
  exit 0
fi

# --- build ------------------------------------------------------------------
python "${SCRIPT}" \
  --quartets "$QUARTETS" \
  --outdir "$OUTDIR" \
  --threads "$THREADS" \
  --expect-min 5000

echo ""
echo "===== DONE ====="
echo "Alignments written: $(find "$OUTDIR" -maxdepth 1 -name '*.fa' | wc -l)"
echo ""
echo "Next:"
echo "  bash   run_molecular_clock_test.sh inspect"
echo "  sbatch run_molecular_clock_test.sh"
echo ""
echo "run_molecular_clock_test.sh already prefers $OUTDIR over the pre-fix"
echo "four_way_alignments/, so no edit is needed there."
