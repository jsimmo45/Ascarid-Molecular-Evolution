#!/bin/bash
#SBATCH -J tox_annotate
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=tox_annotate_%j.out
#SBATCH --error=tox_annotate_%j.err

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
REQUIRED=(PROJECT_DIR TOX_ELIMINATED_BED AS_RETAINED_BED AS_ELIMINATED_BED PU_RETAINED_BED PU_ELIMINATED_BED)
for v in "${REQUIRED[@]}"; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: ${v} is not set in ${CONFIG}" >&2; exit 1; }
done

find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/04_annotation" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script annotate_tox_kaks.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: annotate_tox_kaks.py not found near ${REPO_ROOT}/04_annotation." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV:?set CONDA_ENV in config.sh}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
fi

WORKDIR="${PROJECT_DIR}"
cd "$WORKDIR"

TOX_ELIM="${TOX_ELIMINATED_BED}"

# --- 1. Tox x Ascaris (with Ascaris partner annotation) ---
echo "===== Tox x Ascaris ====="
python "${SCRIPT}" \
  --kaks kaks_Tox_As/kaks_filtered.tsv \
  --tox-eliminated "$TOX_ELIM" \
  --species1 Toxocara --species2 Ascaris \
  --gff2 As_miniprot.gff \
  --retained2 "${AS_RETAINED_BED}" \
  --eliminated2 "${AS_ELIMINATED_BED}" \
  --outdir kaks_Tox_As_annotated

# --- 2. Tox x Parascaris (with Parascaris partner annotation) ---
echo "===== Tox x Parascaris ====="
python "${SCRIPT}" \
  --kaks kaks_Tox_Pu/kaks_filtered.tsv \
  --tox-eliminated "$TOX_ELIM" \
  --species1 Toxocara --species2 Parascaris \
  --gff2 Pu_miniprot.gff \
  --retained2 "${PU_RETAINED_BED}" \
  --eliminated2 "${PU_ELIMINATED_BED}" \
  --outdir kaks_Tox_Pu_annotated

# --- 3. Tox x Baylisascaris (no partner annotation available) ---
echo "===== Tox x Baylisascaris ====="
python "${SCRIPT}" \
  --kaks kaks_Tox_Bsc_v2/kaks_filtered.tsv \
  --tox-eliminated "$TOX_ELIM" \
  --species1 Toxocara --species2 Baylisascaris \
  --outdir kaks_Tox_Bsc_annotated

echo "===== DONE ====="
echo ""
echo "--- Quick summary ---"
for d in kaks_Tox_As_annotated kaks_Tox_Pu_annotated kaks_Tox_Bsc_annotated; do
  echo "==> $d/summary.txt <=="
  cat "$d/summary.txt" 2>/dev/null || echo "  (not found)"
  echo ""
done
