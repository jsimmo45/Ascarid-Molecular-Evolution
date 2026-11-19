#!/bin/bash
#SBATCH -J kaks_bsc_v2
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
# Pairwise dN and dS (MAFFT protein alignment, codon back-translation, yn00)
# for Ascaris, P. univalens and T. canis, each x B. schroederi, on the ortholog
# lists corrected by 02_orthologs/run_rbh_bsc_fix.sh.
#
# REBUILT FROM THE RUN LOGS. The wrapper of the reported runs (2026-08-08) was
# not kept. The arguments below are those on the "Arguments:" line of
# kaks_As_Bsc_v2/, kaks_Pu_Bsc_v2/ and kaks_Tox_Bsc_v2/kaks_pipeline.log. The
# other recorded settings are kaks_pipeline.py defaults: min_align_len 100,
# max_ks 5.0, max_kaks 50.0, sub_rate 1.7e-08, no regions file, keep_tmp off.
#
# Inputs (in PROJECT_DIR): rbh_As_Bsc_v2/, rbh_Pu_Bsc_v2/ and rbh_Tox_Bsc_v2/ from
# 02_orthologs/run_rbh_bsc_fix.sh, and <stem>_cds.fa / <stem>_prot.fa for
# Ascaris (As), P. univalens (Pu), T. canis (Tox) and B. schroederi (Bsc, from
# 01_gene_sets/run_translate_cds.sh).
# A pair whose kaks_filtered.tsv already exists is skipped, so the reported
# tables are never overwritten and a timed-out job resumes where it stopped.
# kaks_pipeline.py runs 48 workers (--threads 48): 1 node, 1 task, 48 CPUs.
# Usage: cd 03_rates && sbatch run_kaks_bsc_v2.sh

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
SCRIPT="$(find_script kaks_pipeline.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: kaks_pipeline.py not found near ${REPO_ROOT}/03_rates." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV:?set CONDA_ENV in config.sh}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
fi


# name : file stem 1 : file stem 2 : species 1 : species 2
# (orthologs from rbh_<name>/, results to kaks_<name>/)
PAIRS=(
  "As_Bsc_v2:As:Bsc:Ascaris:Baylisascaris"
  "Pu_Bsc_v2:Pu:Bsc:Parascaris:Baylisascaris"
  "Tox_Bsc_v2:Tox:Bsc:Toxocara:Baylisascaris"
)

for tool in mafft yn00; do
  command -v "$tool" >/dev/null || { echo "FATAL: $tool not on PATH (conda env kaks)" >&2; exit 1; }
done
cd "${PROJECT_DIR}" || exit 1
echo "python : $(command -v python)"
echo "workdir: $PWD"

# Check every input before the first long run.
for entry in "${PAIRS[@]}"; do
  IFS=':' read -r NAME S1 S2 SP1 SP2 <<< "$entry"
  for f in "rbh_${NAME}/orthologs.tsv" "${S1}_cds.fa" "${S1}_prot.fa" "${S2}_cds.fa" "${S2}_prot.fa"; do
    [[ -s "$f" ]] || { echo "FATAL: missing or empty input for $NAME: ${PROJECT_DIR}/$f" >&2; exit 1; }
  done
done

ran=(); skipped=()
for entry in "${PAIRS[@]}"; do
  IFS=':' read -r NAME S1 S2 SP1 SP2 <<< "$entry"
  OUT="kaks_${NAME}"
  echo "=============================================="
  echo "$NAME   ($SP1 x $SP2)   -> $OUT"
  echo "=============================================="
  if [[ -s "$OUT/kaks_filtered.tsv" ]]; then
    echo "  SKIP: $OUT/kaks_filtered.tsv exists; move $OUT aside to rerun"
    skipped+=("$NAME"); continue
  fi
  python "$SCRIPT" \
    --pairs "rbh_${NAME}/orthologs.tsv" \
    --cds1 "${S1}_cds.fa" --prot1 "${S1}_prot.fa" \
    --cds2 "${S2}_cds.fa" --prot2 "${S2}_prot.fa" \
    --method yn00 \
    --threads 48 \
    --species1 "$SP1" --species2 "$SP2" \
    --outdir "$OUT"
  [[ -s "$OUT/kaks_filtered.tsv" ]] || { echo "FATAL: $OUT/kaks_filtered.tsv not written" >&2; exit 1; }
  echo "  ortholog pairs in: $(wc -l < "rbh_${NAME}/orthologs.tsv")"
  if [[ -s "$OUT/summary.txt" ]]; then sed -n '1,25s/^/  /p' "$OUT/summary.txt"; fi
  ran+=("$NAME")
done

echo "===== DONE =====   ran: ${ran[*]:-none}   skipped: ${skipped[*]:-none}"
