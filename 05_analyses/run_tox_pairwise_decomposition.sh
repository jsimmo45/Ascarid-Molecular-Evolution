#!/bin/bash
#SBATCH -J s11_toxpair
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=s11_toxpair_%j.out
#SBATCH --error=s11_toxpair_%j.err
#
# Supplementary fig. S11: eliminated vs retained dN/dS, dN and dS in the
# pairwise comparisons of T. canis with each ascarid, labelled by T. canis and,
# where the partner has an elimination annotation, by the partner. Reads the
# tables annotate_tox_kaks.py wrote; nothing is re-annotated.
#
#   sbatch run_tox_pairwise_decomposition.sh
#
# Inputs (config.sh; relative paths are taken from PROJECT_DIR):
#   TOX_AS_ANNOTATED    default kaks_Tox_As_annotated/kaks_annotated.tsv
#   TOX_PU_ANNOTATED    default kaks_Tox_Pu_annotated/kaks_annotated.tsv
#   TOX_BSC_ANNOTATED   default kaks_Tox_Bsc_v2_annotated/kaks_annotated.tsv if
#                       that exists, else kaks_Tox_Bsc_annotated/kaks_annotated.tsv
# The partner-labelled rows need a partner_region column; a table without one
# contributes only its T. canis row.
#
# Single-threaded; 2,000 bootstrap replicates per comparison and metric, about
# 30 s. 1 node / 1 task / 1 CPU.

set -euo pipefail

# Under sbatch, SLURM copies the batch script to a spool directory on the compute
# node, so BASH_SOURCE points at /var/spool/slurm/... and neither config.sh nor
# the .py file is anywhere near it. SLURM_SUBMIT_DIR is the directory sbatch was
# run from, which is on the shared filesystem, so search from there as well.
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
  echo "ERROR: no config.sh found." >&2
  echo "  searched upward from: ${HERE}" >&2
  echo "                        ${SLURM_SUBMIT_DIR:-<not under sbatch>}" >&2
  echo "                        ${PWD}" >&2
  echo "  Set KAKS_CONFIG=/full/path/to/config.sh and resubmit." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" \
           "${REPO_ROOT}/06_figures" "${REPO_ROOT}" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script tox_pairwise_decomposition.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: tox_pairwise_decomposition.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHONPATH

# The kaks environment, with its bin/ first on PATH: on Isaac compute nodes the
# Spack python otherwise shadows the conda one. config.sh can override all three.
CONDA_SH="${CONDA_SH:-}"
CONDA_ENV="${CONDA_ENV:-kaks}"
CONDA_BIN_PREPEND="${CONDA_BIN_PREPEND:-}"
if [[ -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
fi
PY="$(command -v python || command -v python3)"

# Probe the imports: paper_colors, region_lib and rate_tables live in lib/.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("paper_colors", "region_lib", "rate_tables", "scipy",
                       "matplotlib")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, scipy, matplotlib all importable")
PYEOF

cd "${PROJECT_DIR}"

TOX_AS="${TOX_AS_ANNOTATED:-kaks_Tox_As_annotated/kaks_annotated.tsv}"
TOX_PU="${TOX_PU_ANNOTATED:-kaks_Tox_Pu_annotated/kaks_annotated.tsv}"
if [[ -n "${TOX_BSC_ANNOTATED:-}" ]]; then
  TOX_BSC="${TOX_BSC_ANNOTATED}"
elif [[ -s kaks_Tox_Bsc_v2_annotated/kaks_annotated.tsv ]]; then
  TOX_BSC=kaks_Tox_Bsc_v2_annotated/kaks_annotated.tsv
else
  TOX_BSC=kaks_Tox_Bsc_annotated/kaks_annotated.tsv
fi

has_col() { head -1 "$1" | tr '\t' '\n' | grep -q -x -F "$2"; }

ARGS=()
missing=0
add_sets() {                        # add_sets <label> <table> <partner name or ->
  if [[ ! -s "$2" ]]; then
    echo "ERROR: missing or empty: $2 ($1)" >&2; missing=1; return
  fi
  has_col "$2" tox_region || { echo "ERROR: no tox_region column in $2" >&2; missing=1; return; }
  ARGS+=(--set "$1" "$2" tox_region "T. canis")
  printf '  %-26s %s  (T. canis labels)\n' "$1" "$2"
  if [[ "$3" != "-" ]]; then
    if has_col "$2" partner_region; then
      ARGS+=(--set "$1" "$2" partner_region "$3")
      printf '  %-26s %s  (%s labels)\n' "$1" "$2" "$3"
    else
      echo "  note: no partner_region in $2; the $3-labelled row is skipped"
    fi
  fi
}
echo "comparisons:"
add_sets "T. canis × Ascaris"       "${TOX_AS}"  "Ascaris"
add_sets "T. canis × P. univalens"  "${TOX_PU}"  "P. univalens"
add_sets "T. canis × B. schroederi" "${TOX_BSC}" "-"
[[ "${missing}" -eq 0 ]] || { echo "  (set TOX_*_ANNOTATED in config.sh)" >&2; exit 1; }

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/tox_pairwise_${STAMP}"
echo "output: ${OUTDIR}"

"${PY}" "${SCRIPT}" "${ARGS[@]}" --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/tox_pairwise_summary.txt"
