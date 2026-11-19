#!/bin/bash
#SBATCH -J sat_diag
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=saturation_%j.out
#SBATCH --error=saturation_%j.err
#
# Supplementary Figure S12 -- dS saturation diagnostic.
#
# Answers the objection the paper invites by using Toxocara at median dS
# 1.39-1.43 against a ~1.5 reliability ceiling for pairwise yn00: how many pairs
# are past the threshold, does dN/dS depend on dS, and does the eliminated-gene
# contrast survive dropping the saturated tail.
#
#   sbatch run_saturation_diagnostic.sh
#   bash   run_saturation_diagnostic.sh
#
# Single-threaded pandas/scipy, hence 1 node / 1 task / 1 CPU. The bootstrap in
# panel C is the slow part; a couple of minutes at the default 2,000 replicates.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

find_config() {
  local d="${HERE}"
  for _ in 1 2 3 4; do
    [[ -f "${d}/config.sh" ]] && { printf '%s' "${d}/config.sh"; return 0; }
    d="$(dirname "${d}")"
  done
  return 1
}
CONFIG="${KAKS_CONFIG:-$(find_config || true)}"
if [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then
  echo "ERROR: no config.sh found in ${HERE} or up to four levels above it." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
else
  echo "No conda config; using python on PATH: $(command -v python || command -v python3)"
fi
PY="$(command -v python || command -v python3)"

cd "${PROJECT_DIR}"

STAMP="$(date +%Y%m%d)"
OUTDIR="${RESULTS_DIR}/saturation_diagnostic_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

# The Toxocara annotated table carries tox_region -- Toxocara's OWN elimination
# calls, not projected from Ascaris. Panel C is skipped without it, and panel C
# is the one that actually tests whether the region contrast depends on the
# saturated tail, so this is worth having.
ANNOTATED="${TOX_ANNOTATED_TABLE:-kaks_Tox_As_annotated/kaks_annotated.tsv}"

echo "Output: ${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
echo "Inputs:"
MISSING=0
for f in "${KAKS_AS_PU}" "${KAKS_AS_BSC}" "${KAKS_PU_BSC}" \
         "${KAKS_TOX_AS}" "${KAKS_TOX_PU}" "${KAKS_TOX_BSC}"; do
  if [[ ! -s "$f" ]]; then
    echo "  MISSING  $f" >&2
    MISSING=1
  else
    printf '  %-44s %8d rows\n' "$f" "$(( $(wc -l < "$f") - 1 ))"
  fi
done
[[ "${MISSING}" -eq 0 ]] || { echo "ERROR: missing rate tables (see above)." >&2; exit 1; }

REGION_ARG=()
if [[ -s "${ANNOTATED}" ]]; then
  printf '  %-44s %8d rows  (panel C)\n' "${ANNOTATED}" \
    "$(( $(wc -l < "${ANNOTATED}") - 1 ))"
  REGION_ARG=(--annotated "${ANNOTATED}" --region-col "${TOX_REGION_COL:-tox_region}")
else
  echo "  NOTE: ${ANNOTATED} not found -- panel C will be skipped." >&2
fi

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT_DIR:-${HERE}}/saturation_diagnostic.py" \
  --as-pu   "${KAKS_AS_PU}" \
  --as-bsc  "${KAKS_AS_BSC}" \
  --pu-bsc  "${KAKS_PU_BSC}" \
  --tox-as  "${KAKS_TOX_AS}" \
  --tox-pu  "${KAKS_TOX_PU}" \
  --tox-bsc "${KAKS_TOX_BSC}" \
  "${REGION_ARG[@]}" \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/saturation_summary.txt"
