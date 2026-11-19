#!/bin/bash
#SBATCH -J s7_microsyn
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=s7_microsyn_%j.out
#SBATCH --error=s7_microsyn_%j.err
#
# Supplementary fig. S7 and table T8: conservation of local gene order between
# Ascaris and P. univalens for retained, eliminated and CBR-spanning adjacent
# gene pairs, with the Fisher test the May script never wrote out, an odds
# ratio with CI, and a sweep over the rank-distance threshold.
#
#   sbatch run_microsynteny.sh
#   sbatch run_microsynteny.sh <path/to/kaks_annotated.tsv>
#
# Inputs (config.sh; relative paths are taken from PROJECT_DIR):
#   ANNOTATED_TABLE   argument 1 overrides it
#   AS_MINIPROT_GFF   default As_miniprot.gff
#   PU_MINIPROT_GFF   default Pu_miniprot.gff
#
# Replaces the May root script of the same name, which still sits in the
# miniprot directory, so the wrapper refuses any microsynteny.py that does not
# import paper_colors.
#
# Single-threaded, two GFF parses and one pass over ~20,000 adjacent pairs:
# under a minute. 1 node / 1 task / 1 CPU.

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
SCRIPT="$(find_script microsynteny.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: microsynteny.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# The May script of the same name is still in the miniprot directory. Only the
# converted one imports the shared figure style, so refuse anything that does not.
grep -q '^from paper_colors import' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the old May script (no paper_colors import)." >&2
  echo "  Use Ascarid-Molecular-Evolution/05_analyses/microsynteny.py." >&2
  exit 1; }

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

TABLE="${1:-${ANNOTATED_TABLE:-}}"
GFF1="${AS_MINIPROT_GFF:-As_miniprot.gff}"
GFF2="${PU_MINIPROT_GFF:-Pu_miniprot.gff}"

[[ -n "${TABLE}" ]] || { echo "ERROR: no table. Pass one as argument 1 or set ANNOTATED_TABLE in config.sh." >&2; exit 1; }
missing=0
for f in "${TABLE}" "${GFF1}" "${GFF2}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; missing=1; }
done
[[ "${missing}" -eq 0 ]] || { echo "  (set AS_MINIPROT_GFF / PU_MINIPROT_GFF in config.sh)" >&2; exit 1; }

# The gene and region columns must exist, or nothing joins.
for col in gene1 gene2 region1; do
  head -1 "${TABLE}" | tr '\t' '\n' | grep -q -x -F "${col}" || {
    echo "ERROR: no ${col} column in ${TABLE}; is this annotate_regions.py output?" >&2
    exit 1; }
done

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/microsynteny_${STAMP}"
printf '%-8s %s\n' "Table:" "${TABLE}" "GFF 1:" "${GFF1}" "GFF 2:" "${GFF2}" "Output:" "${OUTDIR}"

# S7 is emitted whole by this one script, so it carries its own panel letters.
"${PY}" "${SCRIPT}" --annotated "${TABLE}" --gff1 "${GFF1}" --gff2 "${GFF2}" \
  --max-distance 3 --panel-letters --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/microsynteny_summary.txt"
