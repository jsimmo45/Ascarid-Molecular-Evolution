#!/bin/bash
#SBATCH -J rbh_bsc_fix
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=rbh_bsc_fix_%j.out
#SBATCH --error=rbh_bsc_fix_%j.err
#
# Re-run the three Baylisascaris RBH comparisons with the BLAST id
# normalization fix. BLAST uppercases accession-like ids (one letter + five
# digits), so 'g10000' came back as 'G10000' in one direction only and the
# reciprocal merge dropped every Bsc gene above g9999.
#
# --skip-blast reuses the existing blast_1vs2.tsv / blast_2vs1.tsv, so this is
# seconds of pandas work rather than re-running blastp. The original output
# directories are left untouched; results go to *_v2.

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
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/02_orthologs" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script reciprocal_best_hits.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: reciprocal_best_hits.py not found near ${REPO_ROOT}/02_orthologs." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV:?set CONDA_ENV in config.sh}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
fi

cd "${PROJECT_DIR}" || exit 1

# comparison : prot1 : prot2 : label1 : label2
COMPARISONS=(
  "As_Bsc:As_prot.fa:Bsc_prot.fa:Ascaris:Baylisascaris"
  "Pu_Bsc:Pu_prot.fa:Bsc_prot.fa:Parascaris:Baylisascaris"
  "Tox_Bsc:Tox_prot.fa:Bsc_prot.fa:Toxocara:Baylisascaris"
)

for entry in "${COMPARISONS[@]}"; do
  IFS=':' read -r NAME P1 P2 L1 L2 <<< "$entry"
  OLD="rbh_${NAME}"
  NEW="rbh_${NAME}_v2"

  echo "=============================================="
  echo "$NAME   ($L1 vs $L2)"
  echo "=============================================="

  if [[ ! -s "$OLD/blast_1vs2.tsv" || ! -s "$OLD/blast_2vs1.tsv" ]]; then
    echo "  SKIP: missing BLAST output in $OLD" >&2
    continue
  fi

  # Reuse the existing BLAST tables; keep the originals intact.
  mkdir -p "$NEW"
  cp -f "$OLD/blast_1vs2.tsv" "$OLD/blast_2vs1.tsv" "$NEW/"

  python "${SCRIPT}" \
    --prot1 "$P1" --prot2 "$P2" \
    --label1 "$L1" --label2 "$L2" \
    --evalue 1e-10 --min-pident 30 --min-align-len 50 \
    --skip-blast \
    --outdir "$NEW"

  echo "--- before vs after ---"
  echo "  old: $(grep 'RBH pairs (filtered)' $OLD/rbh_stats.txt)"
  echo "  new: $(grep 'RBH pairs (filtered)' $NEW/rbh_stats.txt)"
  echo "  old max Bsc id: $(cut -f2 $OLD/orthologs.tsv | sed 's/[^0-9]//g' | sort -n | tail -1)"
  echo "  new max Bsc id: $(cut -f2 $NEW/orthologs.tsv | sed 's/[^0-9]//g' | sort -n | tail -1)"
  echo ""
done

echo "===== DONE ====="
echo "If the new counts are much larger and max Bsc id now exceeds 9999,"
echo "regenerate downstream: Ka/Ks for these three pairs, then"
echo "branch_dnds_four.py, then branch_dnds_decomposition.py."
