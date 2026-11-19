#!/bin/bash
#SBATCH -J s9a_compart
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=s9a_compart_%j.out
#SBATCH --error=s9a_compart_%j.err
#
# Supplementary Figure S13A (file figS9A_compartments, named before the
# renumbering): dN, dS and dN/dS of retained genes in A vs B compartments of
# the somatic (L1) eigenvector, after converting germline gene positions to
# somatic coordinates and orienting EV1 per chromosome by GC.
#
# 2026-10-19 (results draft note 67): genes above the y-axis limit are open
# triangles at the limit, not a row of dots; no hatch on B (both boxes retained
# blue, as the hand fix); P as plain text with a raised exponent
# (lib/plain_text.py); the axis starts 2% below 0. Statistics unchanged.
#
#   sbatch run_kaks_compartments.sh
#   sbatch run_kaks_compartments.sh <path/to/kaks_annotated.tsv>
#
# Inputs (config.sh; relative paths are taken from PROJECT_DIR):
#   ANNOTATED_TABLE       argument 1 overrides it
#   SOMATIC_EIGENVECTOR   default ${HIC_COMPARTMENT_DIR}/postpde/eigenvectors/
#                         as_10day_iced_100kb.matrix.eigenvector
#   GERM_TO_SOMA_MAP      default ${DB_DIR}/ascaris/v50_postPDE/prepde_to_postpde_coordinates.bed
#   AS_GENOME_FA          default ${DB_DIR}/ascaris/v50/AG.v50.fa (GC for orienting EV1)
#   SOMA_GENOME_FA        the post-elimination assembly; default: the one FASTA in
#                         the coordinate map's folder, if there is exactly one.
#                         With it, each somatic chromosome's strand and offset are
#                         read from sequence instead of assumed (2026-10-12: the
#                         first run flagged seven chromosomes that may be reversed).
#   SOMATIC_EV_LABEL      default "L1 (as_10day)"
#
# Replaces the May root script of the same name (still in the miniprot
# directory), so the wrapper refuses any kaks_compartments.py that does not
# import paper_colors.
#
# Single-threaded; reading the two assemblies is most of the run, a few
# minutes and ~3 GB. 1 node / 1 task / 1 CPU.

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
SCRIPT="$(find_script kaks_compartments.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: kaks_compartments.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# The May script of the same name is still in the miniprot directory. Only the
# converted one imports the shared figure style, so refuse anything that does not.
grep -q '^from paper_colors import' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the old May script (no paper_colors import)." >&2
  echo "  Use Ascarid-Molecular-Evolution/05_analyses/kaks_compartments.py." >&2
  exit 1; }
grep -q '^from plain_text import draw_p' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the version from before 2026-10-19 (capped genes drawn as dots; B hatched)." >&2
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

# Probe the imports: paper_colors, region_lib, rate_tables and plain_text
# live in lib/.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("paper_colors", "region_lib", "rate_tables", "plain_text",
                       "scipy", "matplotlib")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, plain_text, scipy, matplotlib all importable")
PYEOF

cd "${PROJECT_DIR}"

TABLE="${1:-${ANNOTATED_TABLE:-}}"
SOMA_EV="${SOMATIC_EIGENVECTOR:-${HIC_COMPARTMENT_DIR}/postpde/eigenvectors/as_10day_iced_100kb.matrix.eigenvector}"
CMAP="${GERM_TO_SOMA_MAP:-${DB_DIR}/ascaris/v50_postPDE/prepde_to_postpde_coordinates.bed}"
GENOME="${AS_GENOME_FA:-${DB_DIR}/ascaris/v50/AG.v50.fa}"
EV_LABEL="${SOMATIC_EV_LABEL:-L1 (as_10day)}"
SOMA_FA="${SOMA_GENOME_FA:-}"
if [[ -z "${SOMA_FA}" && -d "$(dirname "${CMAP}")" ]]; then
  mapfile -t cands < <(find "$(dirname "${CMAP}")" -maxdepth 1 -type f \
    \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' -o -name '*.fa.gz' \
       -o -name '*.fasta.gz' \) | sort)
  if [[ ${#cands[@]} -eq 1 ]]; then
    SOMA_FA="${cands[0]}"
  else
    echo "note: ${#cands[@]} FASTA files beside the coordinate map; set SOMA_GENOME_FA"
    echo "      in config.sh to the post-elimination assembly to check strands from sequence."
    printf '        %s\n' "${cands[@]}"
  fi
fi
SOMA_ARG=()
[[ -n "${SOMA_FA}" ]] && SOMA_ARG=(--soma-genome "${SOMA_FA}")

[[ -n "${TABLE}" ]] || { echo "ERROR: no table. Pass one as argument 1 or set ANNOTATED_TABLE in config.sh." >&2; exit 1; }
missing=0
for f in "${TABLE}" "${SOMA_EV}" "${CMAP}" "${GENOME}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; missing=1; }
done
[[ "${missing}" -eq 0 ]] || { echo "  (set SOMATIC_EIGENVECTOR / GERM_TO_SOMA_MAP / AS_GENOME_FA in config.sh)" >&2; exit 1; }

# The gene, region and chromosome columns must exist, or nothing joins.
for col in gene1 region1 chrom1; do
  head -1 "${TABLE}" | tr '\t' '\n' | grep -q -x -F "${col}" || {
    echo "ERROR: no ${col} column in ${TABLE}; is this annotate_regions.py output?" >&2
    exit 1; }
done
echo "coordinate map, first lines:"
head -3 "${CMAP}" | sed 's/^/  /'

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/kaks_compartments_${STAMP}"
printf '%-12s %s\n' "Table:" "${TABLE}" "Eigenvector:" "${SOMA_EV}" \
  "Coord map:" "${CMAP}" "Genome:" "${GENOME}" \
  "Somatic:" "${SOMA_FA:-none (strands assumed from the map)}" "Output:" "${OUTDIR}"

"${PY}" "${SCRIPT}" --annotated "${TABLE}" --eigenvector "${SOMA_EV}" \
  --coord-map "${CMAP}" --genome "${GENOME}" --ev-label "${EV_LABEL}" \
  "${SOMA_ARG[@]}" --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/compartment_summary.txt"
