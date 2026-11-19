#!/bin/bash
#SBATCH -J four_way_dnds
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=four_way_dnds_%j.out
#SBATCH --error=four_way_dnds_%j.err
#
# Figure 1C (six-comparison Ks overlay) and Figure 1D (Ks distance tree).
#
# Runs two ways, unchanged:
#     sbatch run_four_way_dnds.sh          on the cluster
#     bash   run_four_way_dnds.sh          anywhere else
#
# Everything machine-specific is in config.sh (see config.example.sh). The SLURM
# directives above are static because sbatch reads them before any shell code
# runs, so they cannot be filled from a variable; they are ignored entirely when
# the script is run with bash.

set -euo pipefail

# --- locate and load config -------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Search upward for config.sh so one config serves every wrapper, whether they
# sit together in one directory or are filed per figure. $KAKS_CONFIG overrides.
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
  echo "       Copy config.example.sh to config.sh and edit it," >&2
  echo "       or point KAKS_CONFIG at one." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

# --- environment ------------------------------------------------------------
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

# --- dated output directory -------------------------------------------------
# A rerun lands in a new directory rather than overwriting the run a figure was
# built from. Pass RUN_TAG=something to distinguish runs made on the same day.
STAMP="$(date +%Y%m%d)"
OUTDIR="${RESULTS_DIR}/four_way_dnds_${STAMP}${RUN_TAG:+_${RUN_TAG}}"
mkdir -p "${OUTDIR}"
echo "Output: ${OUTDIR}"

# --- preflight --------------------------------------------------------------
# Fail before burning a job on a missing or superseded input, and print the row
# count of each table so a pre-fix Baylisascaris file is obvious in the log.
for f in "${KAKS_AS_PU}" "${KAKS_AS_BSC}" "${KAKS_PU_BSC}" \
         "${KAKS_TOX_AS}" "${KAKS_TOX_PU}" "${KAKS_TOX_BSC}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: input not found: ${PROJECT_DIR}/${f}" >&2
    exit 1
  fi
  printf '  %-40s %8d rows\n' "${f}" "$(( $(wc -l < "${f}") - 1 ))"
done

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT_DIR:-${HERE}}/four_way_dnds.py" \
  --as-pu   "${KAKS_AS_PU}" \
  --as-bsc  "${KAKS_AS_BSC}" \
  --pu-bsc  "${KAKS_PU_BSC}" \
  --tox-as  "${KAKS_TOX_AS}" \
  --tox-pu  "${KAKS_TOX_PU}" \
  --tox-bsc "${KAKS_TOX_BSC}" \
  --fit exact \
  --bootstrap 1000 \
  --tree-width single \
  --overlay-width single \
  --outdir "${OUTDIR}"

echo "===== DONE ====="
cat "${OUTDIR}/four_way_summary.txt"
