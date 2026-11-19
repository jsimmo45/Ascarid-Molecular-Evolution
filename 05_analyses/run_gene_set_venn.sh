#!/bin/bash
#SBATCH -J gene_set_venn
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=gene_set_venn_%j.out
#SBATCH --error=gene_set_venn_%j.err
#
# Figure 1 -- Venn diagram of the Ascaris genes in the three pairwise rate tables
# that include Ascaris (with P. univalens, B. schroederi and T. canis), with the
# retained and eliminated genes of every region and the branch-model quartets
# counted against the centre. Added 2026-10-26 for Jianbin Wang's review.
#
#   sbatch run_gene_set_venn.sh
#   bash   run_gene_set_venn.sh
#
# Region calls are made as annotate_regions.py makes them and are checked
# against ANNOTATED_TABLE: the run stops if a single call differs.
#
# Single-threaded pandas: 1 node / 1 task / 1 CPU; the GFF and three tables fit
# in well under 1 GB.

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

# The repo may sit *below* the directory sbatch was run from -- submitting from
# the data directory with `sbatch Ascarid-Molecular-Evolution/05_analyses/...`
# is the obvious thing to do, and an upward-only search never sees it. So probe
# one level down as well, accepting only a directory that also looks like the
# repo, so a stray config.sh elsewhere cannot be picked up by accident.
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
  echo "  Set KAKS_CONFIG=/full/path/to/config.sh, or submit with" >&2
  echo "    sbatch --export=ALL,KAKS_CONFIG=/full/path/config.sh $0" >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

# The .py cannot be located from BASH_SOURCE under sbatch either. Anchor on the
# repo root, which is wherever config.sh was found.
REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
find_script() {                     # find_script <basename>
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" \
           "${REPO_ROOT}/06_figures" "${REPO_ROOT}" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script gene_set_venn.py || true)"
if [[ -z "${SCRIPT}" ]]; then
  echo "ERROR: gene_set_venn.py not found near ${REPO_ROOT}." >&2
  exit 1
fi
echo "script: ${SCRIPT}"

# paper_colors.py lives in lib/ in the published layout but beside the scripts in
# a flat working tree. Put both on PYTHONPATH rather than depending on which.
export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

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

GFF="${AS_MINIPROT_GFF:-As_miniprot.gff}"
RET="${AS_RETAINED_BED:-${DB_DIR}/ascaris/v50/AG.v50.retained_loose.bed}"
ELIM="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
SIZES="${AS_CHROM_SIZES:-}"
# The rate tables every pairwise figure uses (pairs that passed the rate
# filters). The Ascaris gene is gene1 in the first two and gene2 in the T. canis
# table, which was run with T. canis as species 1.
T_PU="${KAKS_AS_PU:-kaks_results/kaks_filtered.tsv}"
T_BSC="${KAKS_AS_BSC:-kaks_As_Bsc_v2/kaks_filtered.tsv}"
T_TOX="${KAKS_TOX_AS:-kaks_Tox_As/kaks_filtered.tsv}"
QUARTETS="${QUARTET_TABLE:-branch_dnds_four_results_v2/branch_dnds_four_annotated.tsv}"
ANNOT="${ANNOTATED_TABLE:-results/annotate_regions_20260817/kaks_annotated.tsv}"
TERMINAL_TOL="${TERMINAL_TOL:-1000}"

STAMP="$(date +%Y%m%d)"
OUTDIR="${RESULTS_DIR}/gene_set_venn_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

echo "Models:     ${GFF}"
echo "Retained:   ${RET}"
echo "Eliminated: ${ELIM}"
echo "Tables:     ${T_PU} | ${T_BSC} | ${T_TOX}"
echo "Quartets:   ${QUARTETS}"
echo "Check:      ${ANNOT}"
echo "Output:     ${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${GFF}" "${RET}" "${ELIM}" "${T_PU}" "${T_BSC}" "${T_TOX}" "${QUARTETS}" "${ANNOT}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done
SIZES_ARG=()
if [[ -n "${SIZES}" ]]; then
  [[ -s "${SIZES}" ]] || { echo "ERROR: missing or empty: ${SIZES}" >&2; exit 1; }
  SIZES_ARG=(--chrom-sizes "${SIZES}")
fi

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --gff "${GFF}" \
  --retained-bed "${RET}" \
  --eliminated-bed "${ELIM}" \
  "${SIZES_ARG[@]}" \
  --terminal-tol "${TERMINAL_TOL}" \
  --set "parascaris=${T_PU}:gene1" \
  --set "baylisascaris=${T_BSC}:gene1" \
  --set "toxocara=${T_TOX}:gene2" \
  --quartets "${QUARTETS}:as_gene" \
  --annotated "${ANNOT}" \
  --width single \
  --height-mm 70 \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
