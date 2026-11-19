#!/bin/bash
#SBATCH -J rbh_bsc_blast
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# All-against-all blastp and reciprocal best hits for Ascaris x B. schroederi
# and P. univalens x B. schroederi.
#
# Steps 2-3 of the original run_baylisascaris.sh (May 2026): same commands,
# inputs and output folders, and the arguments recorded on the "Arguments:"
# line of rbh_As_Bsc/rbh.log and rbh_Pu_Bsc/rbh.log. The other recorded
# settings are reciprocal_best_hits.py defaults: e-value 1e-10,
# max_target_seqs 5, minimum identity 30%, minimum alignment 50 aa.
#
# Their BLAST tables (blast_1vs2.tsv, blast_2vs1.tsv in rbh_As_Bsc/ and
# rbh_Pu_Bsc/) are the ones run_rbh_bsc_fix.sh re-merged into the reported
# ortholog lists (rbh_*_Bsc_v2). The lists these runs wrote at the time lost
# every B. schroederi gene above g9999 and were not used;
# reciprocal_best_hits.py now carries the fix, so a rerun writes correct lists
# here as well.
#
# Order: after 01_gene_sets/run_translate_cds.sh; before run_rbh_bsc_fix.sh.
# A pair whose orthologs.tsv already exists is skipped, so the reported runs are
# never overwritten and a timed-out job resumes where it stopped.
# blastp runs on 48 threads (--threads 48): 1 node, 1 task, 48 CPUs.
# Usage: cd 02_orthologs && sbatch run_rbh_bsc.sh

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


# name : protein set 1 : protein set 2 : label 1 : label 2   (output rbh_<name>/)
COMPARISONS=(
  "As_Bsc:As_prot.fa:Bsc_prot.fa:Ascaris:Baylisascaris"
  "Pu_Bsc:Pu_prot.fa:Bsc_prot.fa:Parascaris:Baylisascaris"
)

command -v blastp >/dev/null || { echo "FATAL: blastp not on PATH (conda env kaks)" >&2; exit 1; }
cd "${PROJECT_DIR}" || exit 1
echo "python : $(command -v python)"
echo "blastp : $(command -v blastp)"
echo "workdir: $PWD"

# Check every input before the first long run.
for entry in "${COMPARISONS[@]}"; do
  IFS=':' read -r NAME P1 P2 L1 L2 <<< "$entry"
  for f in "$P1" "$P2"; do
    [[ -s "$f" ]] || { echo "FATAL: missing or empty input for $NAME: ${PROJECT_DIR}/$f" >&2; exit 1; }
  done
done

ran=(); skipped=()
for entry in "${COMPARISONS[@]}"; do
  IFS=':' read -r NAME P1 P2 L1 L2 <<< "$entry"
  OUT="rbh_${NAME}"
  echo "=============================================="
  echo "$NAME   ($L1 x $L2)   -> $OUT"
  echo "=============================================="
  if [[ -s "$OUT/orthologs.tsv" ]]; then
    echo "  SKIP: $OUT/orthologs.tsv exists; move $OUT aside to rerun"
    skipped+=("$NAME"); continue
  fi
  python "$SCRIPT" \
    --prot1 "$P1" \
    --prot2 "$P2" \
    --label1 "$L1" --label2 "$L2" \
    --threads 48 \
    --outdir "$OUT"
  [[ -s "$OUT/orthologs.tsv" ]] || { echo "FATAL: $OUT/orthologs.tsv not written" >&2; exit 1; }
  grep -E 'with BLAST hits|RBH pairs|Filters' "$OUT/rbh_stats.txt" | sed 's/^/  /' || true
  ran+=("$NAME")
done

echo "===== DONE =====   ran: ${ran[*]:-none}   skipped: ${skipped[*]:-none}"
echo "Next: run_rbh_bsc_fix.sh (re-merges these BLAST tables into rbh_*_Bsc_v2)."
