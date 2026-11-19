#!/bin/bash
#SBATCH -J mv_insul
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=mv_insul_%j.out
#SBATCH --error=mv_insul_%j.err
#
# Figure 6D and 6E (two single-column panels), Supplementary Figure S9B,
# Table S5.
#
#   sbatch run_multivariate_insulation.sh               # teste, all three EV runs
#   sbatch run_multivariate_insulation.sh teste 60hr    # one run
#
#   STAGE  insulation + expression stage:
#          teste | ovary | 0hr | 48hr | 60hr                    (default teste)
#   EV     eigenvector for the compartment control:
#          none | teste | ovary | 0hr | 48hr | 60hr | all       (default all)
#          'all' = none, teste and 60hr in turn: the three rows of Table S5.
#
# 2026-11-01: STAGE 48hr (2-4-cell embryo, before elimination) and 60hr (4-cell,
# as elimination begins) added for the stage series of the insulation result,
# each with the workbook's early_embryo expression, the pre-elimination embryo
# sample; testis, ovary and zygote keep their own columns. Run them with EV
# none:  sbatch run_multivariate_insulation.sh 48hr none  (and 60hr). There is
# no post-elimination stage: after elimination the eliminated DNA is gone from
# the somatic cells, so eliminated genes have no insulation to score.
#
# Figure 6D and 6E come from the noEV run (identical in all three runs, since
# neither panel uses EV1); S9B from the 60hr run, whose EV1 has the clearest
# compartment signal.
#
# 2026-10-19 (results draft note 67): the expression term is labelled by its
# stage, "Testis expression" on 6E and "Testis expr." on S9B (now S13B), since
# it is testis TPM, not S12A's nine-stage mean. Labels only; run after
# run_functional_enrichment.sh so 6D takes the new 6C's height.
#
# 2026-10-18: Figure 6 is laid out 6A-B (strip), then 6C | 6D, then 6E | 6F, so
# 6D and 6E are now separate 89 mm panels (6D_insulation_by_region,
# 6E_insulation_models) instead of the 185 mm 6DE strip. 6D is drawn at 6C's
# height, read from the newest results/functional_enrichment_*/
# 6C_functional_class_forest.svg (D_HEIGHT_MM overrides it; 60 mm if none is
# found). 6E is E_HEIGHT_MM tall, default 48.2 mm, the height of the 6F
# schematic beside it. Statistics are unchanged: the run should reproduce the 2026-10-17
# numbers printed at the end.
#
# 2026-10-17: expression is joined by position, through
# workbook_model_bridge.tsv (run_workbook_model_bridge.sh), not by the GFF
# Target= name, which matched each model to a different gene (results draft
# note 60). The core set and every model with expression in it change. The
# bridge is WB_BRIDGE in config.sh, or else the newest
# results/workbook_model_bridge_*/ folder. 6D now prints P without asterisks
# and with its exponent as raised 6 pt text (lib/plain_text.py), the hand edit
# of 2026-10-15. After the run, point supplementary_table_sources.tsv at the
# three new folders.
#
# Replaces run_kaks_multivariate.sh and kaks_multivariate.py, loose in the
# miniprot directory. What changes:
#   * Gene positions come from ANNOTATED_TABLE (best alignment), not the first
#     GFF record.
#   * Nearest intergenic distance is added as a covariate (log10), reported
#     with and without, because it is a mediator. GENE_STRUCTURE_TABLE in
#     config.sh, or the 08-17 gene_structure run by default. NO_IG=1 skips it.
#   * 6E shows the dN/dS model and the insulation model.
#   * Each summary ends with an explicit Table S5 block, so the last
#     is_eliminated line is still the right one for make_robustness_table.sh.
#
# Single-threaded numpy/scipy: 1 node / 1 task / 1 CPU. 8G covers the genome
# FASTA (only the sequences the table uses are kept).

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
SCRIPT="$(find_script multivariate_insulation.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: multivariate_insulation.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
grep -q -- "--bridge" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} predates the position bridge (no --bridge option)." >&2
  exit 1; }
grep -q -- "'--d-height-mm'" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} still draws the 6DE strip (the version from before 2026-10-18)." >&2
  exit 1; }
grep -q "PRETTY\['log_expr'\] = " "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the version from before 2026-10-19 (the testis term is labelled 'Expression')." >&2
  exit 1; }
[[ -f "$(dirname "${SCRIPT}")/workbook_model_bridge.py" ]] || {
  echo "ERROR: workbook_model_bridge.py is not beside ${SCRIPT}." >&2; exit 1; }

SDIR="$(dirname "${SCRIPT}")"
export PYTHONPATH="${SDIR}:${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

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

"${PY}" - <<'PYEOF' || { echo "ERROR: fix the environment (lib/ on PYTHONPATH, lib/plain_text.py present, openpyxl installed, workbook_model_bridge.py beside the script) and resubmit." >&2; exit 1; }
import importlib, importlib.util, sys
missing = [m for m in ("paper_colors", "region_lib", "rate_tables", "plain_text",
                       "workbook_model_bridge", "openpyxl")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
wmb = importlib.import_module("workbook_model_bridge")
if not hasattr(wmb, "read_bridge"):
    print("ERROR: workbook_model_bridge.py has no read_bridge(); copy the "
          "2026-10-17 version", file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, plain_text, "
      "workbook_model_bridge, openpyxl all importable")
PYEOF

cd "${PROJECT_DIR}"

STAGE="${1:-teste}"
EV="${2:-all}"
case "${STAGE}" in
  teste) EXPR_COL=testis; INS_TOK=teste; STAGE_LABEL="Testis" ;;
  ovary) EXPR_COL=ovary;  INS_TOK=ovary; STAGE_LABEL="Ovary" ;;
  0hr)   EXPR_COL=zygote; INS_TOK=0hr;   STAGE_LABEL="1-cell (0 h)" ;;
  48hr)  EXPR_COL=early_embryo; INS_TOK=48hr; STAGE_LABEL="48 h embryo" ;;
  60hr)  EXPR_COL=early_embryo; INS_TOK=60hr; STAGE_LABEL="4-cell (60 h)" ;;
  *) echo "ERROR: unknown STAGE '${STAGE}' (teste | ovary | 0hr | 48hr | 60hr)" >&2; exit 1 ;;
esac
if [[ "${EV}" == "all" ]]; then EV_LIST="none teste 60hr"; else EV_LIST="${EV}"; fi
ev_label() {
  case "$1" in
    teste) echo "testis" ;;          ovary) echo "ovary" ;;
    0hr)   echo "1-cell (0 h)" ;;    48hr)  echo "48 h embryo" ;;
    60hr)  echo "4-cell (60 h)" ;;   *)     echo "$1" ;;
  esac
}

# Every path gets a default, so a variable missing from config.sh is reported
# by the preflight rather than killing the script under `set -u`.
TABLE="${ANNOTATED_TABLE:-}"
[[ -n "${TABLE}" ]] || { echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2; exit 1; }
XLSX="${EXPRESSION_XLSX:-Ascaris_genes_15714_annotation_expression.xlsx}"
# The position bridge: WB_BRIDGE in config.sh, else the newest bridge run.
BRIDGE="${WB_BRIDGE:-}"
if [[ -z "${BRIDGE}" ]]; then
  LATEST="$(ls -1d "${RESULTS_DIR}"/workbook_model_bridge_*/ 2>/dev/null | sort | tail -1 || true)"
  [[ -n "${LATEST}" ]] && BRIDGE="${LATEST%/}/workbook_model_bridge.tsv"
fi
[[ -n "${BRIDGE}" ]] || {
  echo "ERROR: no workbook_model_bridge.tsv (WB_BRIDGE unset and none under ${RESULTS_DIR}/)." >&2
  echo "  Run run_workbook_model_bridge.sh first, or set WB_BRIDGE in config.sh." >&2
  exit 1; }
GENOME="${AS_GENOME_FA:-${DB_DIR}/ascaris/v50/AG.v50.fa}"
INS_DIR="${INSULATION_DIR:-${ASCARIS_HIC_DIR}/insulation_score_updated/insulation_scores/prepde}"
EV_DIR="${EIGENVECTOR_DIR:-${HIC_COMPARTMENT_DIR}/prepde/eigenvectors}"
GS_TABLE="${GENE_STRUCTURE_TABLE:-${RESULTS_DIR}/gene_structure_20260817_loci/gene_structure_per_gene.tsv}"
GC_WINDOW="${GC_WINDOW:-20000}"
NO_IG="${NO_IG:-0}"
INSULATION="${INS_DIR}/AG.v50.${INS_TOK}.is100001.ids80001.insulation.bedGraph"

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"

# Panel heights. 6D shares Figure 6's middle row with 6C, so it is drawn at the
# height of the newest 6C (matplotlib writes the SVG height in points).
D_SRC="D_HEIGHT_MM"
D_HEIGHT="${D_HEIGHT_MM:-}"
if [[ -z "${D_HEIGHT}" ]]; then
  FE_LATEST="$(ls -1d "${RESULTS_DIR}"/functional_enrichment_*/ 2>/dev/null | sort | tail -1 || true)"
  C_SVG="${FE_LATEST%/}/6C_functional_class_forest.svg"
  if [[ -n "${FE_LATEST}" && -s "${C_SVG}" ]]; then
    D_HEIGHT="$("${PY}" - "${C_SVG}" <<'PYEOF'
import re, sys
head = open(sys.argv[1], encoding='utf-8', errors='replace').read(4000)
m = re.search(r'<svg[^>]*?\sheight="([0-9.]+)pt"', head, re.S)
print(f'{float(m.group(1)) * 25.4 / 72:.1f}' if m else '')
PYEOF
)"
    [[ -n "${D_HEIGHT}" ]] && D_SRC="6C, ${C_SVG}"
  fi
fi
if [[ -z "${D_HEIGHT}" ]]; then D_HEIGHT=60; D_SRC="default (no 6C found)"; fi
E_HEIGHT="${E_HEIGHT_MM:-48.2}"

printf '%-12s %s\n' "Table:" "${TABLE}" "Bridge:" "${BRIDGE}" "Workbook:" "${XLSX}" \
  "Genome:" "${GENOME}" "Insulation:" "${INSULATION}" "Stage:" "${STAGE} (${EXPR_COL})" \
  "EV runs:" "${EV_LIST}" "Intergenic:" "$([[ "${NO_IG}" == 1 ]] && echo skipped || echo "${GS_TABLE}")" \
  "6D height:" "${D_HEIGHT} mm (${D_SRC})" "6E height:" "${E_HEIGHT} mm"
echo ""

# --- preflight --------------------------------------------------------------
missing=0
for f in "${TABLE}" "${BRIDGE}" "${XLSX}" "${GENOME}" "${INSULATION}"; do
  [[ -s "$f" ]] || { echo "  MISSING: $f" >&2; missing=1; }
done
IG_ARG=()
if [[ "${NO_IG}" != "1" ]]; then
  if [[ -s "${GS_TABLE}" ]]; then
    IG_ARG=(--gene-structure "${GS_TABLE}")
  else
    echo "  MISSING: ${GS_TABLE}" >&2
    echo "  (set GENE_STRUCTURE_TABLE in config.sh, or NO_IG=1 to run without it)" >&2
    missing=1
  fi
fi
for ev in ${EV_LIST}; do
  case "${ev}" in
    none) ;;
    teste|ovary|0hr|48hr|60hr)
      f="${EV_DIR}/as_${ev}_iced_100kb.matrix.eigenvector"
      [[ -s "$f" ]] || { echo "  MISSING: $f" >&2; missing=1; } ;;
    *) echo "ERROR: unknown EV '${ev}' (none | teste | ovary | 0hr | 48hr | 60hr | all)" >&2
       exit 1 ;;
  esac
done
[[ "${missing}" -eq 0 ]] || { echo "ERROR: inputs missing, nothing run." >&2; exit 1; }

# --- runs -------------------------------------------------------------------
OUTDIRS=()
for ev in ${EV_LIST}; do
  if [[ "${ev}" == "none" ]]; then
    TAG="${STAGE}_noEV"; EV_ARGS=()
  else
    TAG="${STAGE}_EV${ev}"
    # Pre-elimination eigenvectors are in germline coordinates: no coordinate
    # map. --orient-ev-by-gc fixes FAN-C's per-chromosome sign (A = high GC).
    EV_ARGS=(--eigenvector "${EV_DIR}/as_${ev}_iced_100kb.matrix.eigenvector"
             --orient-ev-by-gc --ev-label "$(ev_label "${ev}")")
  fi
  OUTDIR="${RESULTS_DIR}/multivariate_insulation_${STAMP}_${TAG}${RUN_TAG:+_${RUN_TAG}}"
  OUTDIRS+=("${OUTDIR}")
  mkdir -p "${OUTDIR}"
  echo ""
  echo "===== ${TAG} ====="
  "${PY}" "${SCRIPT}" \
    --kaks "${TABLE}" \
    --bridge "${BRIDGE}" \
    --insulation "${INSULATION}" \
    --expression "${XLSX}" --expr-col "${EXPR_COL}" \
    --genome "${GENOME}" --gc-window "${GC_WINDOW}" \
    "${IG_ARG[@]}" "${EV_ARGS[@]}" \
    --stage-label "${STAGE_LABEL}" \
    --d-height-mm "${D_HEIGHT}" --e-height-mm "${E_HEIGHT}" \
    --seed 0 \
    --outdir "${OUTDIR}" > "${OUTDIR}/console.txt" 2>&1 \
    || { echo "  FAILED: see ${OUTDIR}/console.txt"; tail -5 "${OUTDIR}/console.txt"; continue; }
  echo "  done: ${OUTDIR}"
done

echo ""
echo "===== DONE ====="
first=1
for d in "${OUTDIRS[@]}"; do
  echo "--- $(basename "$d") ---"
  if [[ ! -s "$d/multivariate_summary.txt" ]]; then echo "  no summary (run failed)"; continue; fi
  if [[ "${first}" == 1 ]]; then
    # Identical in every run except the EV1 line, so shown once.
    sed -n '/^Genes at each step/,/^$/p' "$d/multivariate_summary.txt" | sed 's/^/  /'
    grep -A1 '^6D:' "$d/multivariate_summary.txt" | sed 's/^/  /'
    first=0
  fi
  grep 'core set with EV1\|Sign corrected\|Genome-wide Spearman(EV1' "$d/multivariate_summary.txt" | sed 's/^ */  /' || true
  sed -n '/^Figure audit/,/^$/p' "$d/multivariate_summary.txt" | sed 's/^/  /'
done

echo "===== 2026-10-17 RUN (job 6427846, position bridge) vs THIS RUN: must match ====="
echo "2026-10-17 (noEV run): n = 9,216 (8,907 retained, 309 eliminated);"
echo "  Model 1 is_eliminated beta = +0.7963, P = 2.632e-40; insulation -0.0243 (P = 2.082e-02);"
echo "          log_expr -0.1176; gc +0.0499; R2 = 0.0349"
echo "  Model 2 is_eliminated beta = -0.9385, P = 1.242e-57; log_expr +0.1238; gc -0.1164; R2 = 0.0532"
echo "  partial correlation dN/dS ~ insulation | expression, GC: r = -0.0473, P = 5.628e-06"
echo "  6D medians: retained +0.0494, eliminated -0.1262; P = 1.66e-24"
for d in "${OUTDIRS[@]}"; do
  [[ "$d" == *_noEV* && -s "$d/multivariate_summary.txt" ]] || continue
  echo "this run ($(basename "$d")):"
  sed -n '/^PRIMARY MODELS/,/^LOCAL GENE SPACING\|^COMPARTMENT\|^S9B/p' "$d/multivariate_summary.txt" \
    | grep -v '^LOCAL GENE SPACING\|^COMPARTMENT\|^S9B' | sed 's/^/  /'
  grep -A2 '^LOCAL GENE SPACING' "$d/multivariate_summary.txt" >/dev/null && \
    grep '^  Model [12]: is_eliminated beta' "$d/multivariate_summary.txt" | sed 's/^/  /'
done
echo ""
echo "Table S5 (is_eliminated, insulation model). 2026-10-17: none -0.9385 (se 0.0582, P 1.242e-57);"
echo "  testis EV1 -0.9330 (0.0583, 6.902e-57); 4-cell EV1 -0.9387 (0.0583, 1.300e-57)."
for d in "${OUTDIRS[@]}"; do
  [[ -s "$d/multivariate_summary.txt" ]] || continue
  printf '  %-58s %s\n' "$(basename "$d")" \
    "$(grep 'is_eliminated' "$d/multivariate_summary.txt" | tail -1 | sed 's/^ *//')"
done
echo ""
echo "Figure 6D: the noEV run's 6D_insulation_by_region.svg (beside 6C)."
echo "Figure 6E: the noEV run's 6E_insulation_models.svg (beside 6F)."
echo "S9B: the EV60hr run's S9B_feature_correlations.svg."
