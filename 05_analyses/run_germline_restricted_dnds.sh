#!/bin/bash
#SBATCH -J germ_restr_dnds
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=germ_restr_dnds_%j.out
#SBATCH --error=germ_restr_dnds_%j.err
#
# Reviewer control for Figure 6 (M. Labrador, 2026-10-16): is the dN/dS excess
# of eliminated genes a germline-restriction effect? Compares eliminated with
# retained germline-restricted genes (somatic expression at the noise floor of
# genes absent from the soma), retained germline-restricted with somatic-peak
# genes, a stage- and level-matched version, the Figure 6E model with tissue
# specificity added, and a check of region1 against the workbook's own
# elimination call.
#
#   sbatch run_germline_restricted_dnds.sh
#   bash   run_germline_restricted_dnds.sh
#
# Optional settings (environment or config.sh):
#   FLOOR_SOURCE      workbook (default) | analysed | concordant
#   FLOOR_QUANTILES   default 0.95,0.90,0.75 (first = default definition)
#   TAU_THRESHOLDS    default 0.7,0.8,0.9
#   GERMLINE_STAGES   default testis,ovary,zygote,early_embryo (as Figure 6)
#   CALIPER           default 0.25 (log2 peak TPM + 1)
#   B_HEIGHT_MM       S12B height; default: the height of the newest S12A
#                     (results/functional_enrichment_*/S6_dnds_vs_expression.svg),
#                     else 72 mm
#   C_WIDTH           S12C width: single (default) | onehalf | double
#   C_HEIGHT_MM       S12C height; default S12_BOTTOM_MM
#   S12_BOTTOM_MM     height of the S12 bottom row, C | D (default 80); D is
#                     drawn at C's height by run_dnds_within_class.sh
#   RUN_TAG           appended to the output folder name
#
# 2026-10-19: the figure is Supplementary Figure S12 B and C, written as two
# files (S12B_germline_groups, S12C_germline_contrasts) with no panel letters.
# S12B sits beside S12A (functional_enrichment.py's S6_dnds_vs_expression) and
# is drawn at its height, so run functional_enrichment first (afterok).
# Later the same day: S12C is 89 mm wide, beside S12D (dnds_within_class.py)
# on the bottom row, both S12_BOTTOM_MM tall; C's n column is a table axes, as
# 6C's. Run this before run_dnds_within_class.sh, which reads C's height.
#
# 2026-10-17: expression is joined by position, through
# workbook_model_bridge.tsv, as Figure 6 now is (results draft note 60; the
# 2026-10-16 run used the superseded name join and is void). The bridge is
# WB_BRIDGE in config.sh, or else the newest results/workbook_model_bridge_*/.
#
# Stage definitions are imported from functional_enrichment.py, the model from
# multivariate_insulation.py and the bridge reader from workbook_model_bridge.py,
# so all three must be the 2026-10-19 versions in 05_analyses/.
# Single-threaded pandas/scipy: 1 node / 1 task / 1 CPU. About a minute on a table of the real size; 8G covers the workbook,
# the GFF and the annotated table.

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
SCRIPT="$(find_script germline_restricted_dnds.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: germline_restricted_dnds.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
grep -q -- "'--b-height-mm'" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} still draws the two-panel figure (the version from before 2026-10-19)." >&2
  exit 1; }
grep -q -- "ax_t=ax_t" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} draws S12C with its n column off the plot's edge (before the 89 mm S12C of 2026-10-19)." >&2
  exit 1; }
SDIR="$(dirname "${SCRIPT}")"
# The imported modules must sit beside the script and be the converted ones:
# the May functional_enrichment.py in the miniprot directory has no
# target_bridge and would fail on import, or worse, be found first.
for m in functional_enrichment.py multivariate_insulation.py workbook_model_bridge.py; do
  [[ -f "${SDIR}/${m}" ]] || { echo "ERROR: ${m} is not beside ${SCRIPT}." >&2; exit 1; }
  [[ "${m}" == workbook_model_bridge.py ]] && continue
  grep -q '^from paper_colors import' "${SDIR}/${m}" || {
    echo "ERROR: ${SDIR}/${m} is an unconverted script (no paper_colors import)." >&2
    exit 1; }
done

export PYTHONPATH="${SDIR}:${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

# The kaks environment, with its bin/ first on PATH: on Isaac compute nodes the
# Spack python otherwise shadows the conda one. The defaults are the three
# lines every Ka/Ks job carries; config.sh can override them.
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

# Import, not just find: a module that exists but fails to import, or an old
# functional_enrichment.py without target_bridge, is what a file test misses.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix the environment (lib/ on PYTHONPATH, openpyxl installed, converted 05_analyses scripts) and resubmit." >&2; exit 1; }
import importlib, importlib.util, sys
missing = [m for m in ("paper_colors", "region_lib", "rate_tables", "plain_text",
                       "openpyxl")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
fe = importlib.import_module("functional_enrichment")
mi = importlib.import_module("multivariate_insulation")
wmb = importlib.import_module("workbook_model_bridge")
need = [("functional_enrichment", fe, n) for n in ("target_bridge", "load_workbook", "fit_to_canvas",
                                                   "blank_figure", "add_axes_mm", "offcanvas_text")] + \
       [("multivariate_insulation", mi, "fit_model")] + \
       [("workbook_model_bridge", wmb, n) for n in ("read_bridge", "shared_gene_counts")]
lacking = [f"{a}.{n}" for a, mod, n in need if not hasattr(mod, n)]
if lacking:
    print("ERROR: missing in the imported scripts: " + ", ".join(lacking), file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, plain_text, openpyxl, "
      "functional_enrichment, multivariate_insulation, workbook_model_bridge")
PYEOF

# Working directory, explicitly: relative paths in config.sh are relative to it.
cd "${PROJECT_DIR}"

TABLE="${ANNOTATED_TABLE:-}"
[[ -n "${TABLE}" ]] || { echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2; exit 1; }
XLSX="${EXPRESSION_XLSX:-Ascaris_genes_15714_annotation_expression.xlsx}"
GFF1="${AS_MINIPROT_GFF:-As_miniprot.gff}"
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
FLOOR_SOURCE="${FLOOR_SOURCE:-workbook}"
FLOOR_QUANTILES="${FLOOR_QUANTILES:-0.95,0.90,0.75}"
TAU_THRESHOLDS="${TAU_THRESHOLDS:-0.7,0.8,0.9}"
GERMLINE_STAGES="${GERMLINE_STAGES:-testis,ovary,zygote,early_embryo}"
CALIPER="${CALIPER:-0.25}"
C_WIDTH="${C_WIDTH:-single}"
C_HEIGHT="${C_HEIGHT_MM:-${S12_BOTTOM_MM:-80}}"

# S12B shares the supplement's top row with S12A, so it is drawn at the height
# of the newest S12A (matplotlib writes the SVG height in points).
B_SRC="B_HEIGHT_MM"
B_HEIGHT="${B_HEIGHT_MM:-}"
if [[ -z "${B_HEIGHT}" ]]; then
  FE_LATEST="$(ls -1d "${RESULTS_DIR}"/functional_enrichment_*/ 2>/dev/null | sort | tail -1 || true)"
  A_SVG="${FE_LATEST%/}/S6_dnds_vs_expression.svg"
  if [[ -n "${FE_LATEST}" && -s "${A_SVG}" ]]; then
    B_HEIGHT="$("${PY}" - "${A_SVG}" <<'PYEOF'
import re, sys
head = open(sys.argv[1], encoding='utf-8', errors='replace').read(4000)
m = re.search(r'<svg[^>]*?\sheight="([0-9.]+)pt"', head, re.S)
print(f'{float(m.group(1)) * 25.4 / 72:.1f}' if m else '')
PYEOF
)"
    [[ -n "${B_HEIGHT}" ]] && B_SRC="S12A, ${A_SVG}"
  fi
fi
if [[ -z "${B_HEIGHT}" ]]; then B_HEIGHT=72; B_SRC="default (no S12A found)"; fi

# Seconds in the stamp: a date-only stamp sends two runs on the same day into
# one folder, the second overwriting the first.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/germline_restricted_dnds_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-17s %s\n' "Table:" "${TABLE}" "Workbook:" "${XLSX}" "Bridge:" "${BRIDGE}" "GFF:" "${GFF1}" \
  "Floor source:" "${FLOOR_SOURCE}" "Floor quantiles:" "${FLOOR_QUANTILES}" \
  "tau thresholds:" "${TAU_THRESHOLDS}" "Germline stages:" "${GERMLINE_STAGES}" \
  "Caliper:" "${CALIPER}" "S12B height:" "${B_HEIGHT} mm (${B_SRC})" \
  "S12C size:" "${C_WIDTH}, $([[ -n "${C_HEIGHT}" ]] && echo "${C_HEIGHT} mm" || echo "height from rows")" \
  "Output:" "${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${TABLE}" "${XLSX}" "${BRIDGE}" "${GFF1}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done

# The rate table's gene1 IDs must be the GFF's mRNA IDs, or nothing joins and
# every count below is zero.
COL="$(head -1 "${TABLE}" | tr '\t' '\n' | grep -n -x -F gene1 | cut -d: -f1 | head -1)"
[[ -n "${COL}" ]] || { echo "ERROR: no gene1 column in ${TABLE}" >&2; exit 1; }
MATCH="$(awk -v col="${COL}" -F'\t' '
  NR==FNR { if ($3=="mRNA") { split($9, a, /[;=]/); for (i=1;i<length(a);i++) if (a[i]=="ID") ids[a[i+1]]=1 } next }
  FNR==1 { next }
  { if ($col in ids) n++ }
  END { print n+0 }' "${GFF1}" "${TABLE}")"
ROWS=$(( $(wc -l < "${TABLE}") - 1 ))
printf '  %-28s %6d / %6d table rows found in the GFF\n' "gene1 identifiers" "${MATCH}" "${ROWS}"
[[ "${MATCH}" -gt $(( ROWS / 2 )) ]] || {
  echo "ERROR: fewer than half the gene1 identifiers appear in ${GFF1}." >&2; exit 1; }

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --kaks "${TABLE}" \
  --annotation "${XLSX}" \
  --bridge "${BRIDGE}" \
  --gff1 "${GFF1}" \
  --floor-source "${FLOOR_SOURCE}" \
  --floor-quantiles "${FLOOR_QUANTILES}" \
  --tau-thresholds "${TAU_THRESHOLDS}" \
  --germline-stages "${GERMLINE_STAGES}" \
  --caliper "${CALIPER}" \
  --b-height-mm "${B_HEIGHT}" \
  --c-width "${C_WIDTH}" \
  ${C_HEIGHT:+--c-height-mm "${C_HEIGHT}"} \
  --seed 0 \
  --outdir "${OUTDIR}"

SUM="${OUTDIR}/germline_restricted_summary.txt"
echo ""
echo "===== DONE ====="
echo "Summary: ${SUM}"
echo ""

echo "===== CHECK AGAINST FIGURE 6A ====="
echo "functional_enrichment.py with the same bridge must give the same peak-stage counts"
echo "and shares (its 'with a stage of peak expression' line and 6A category rows)."
echo "Superseded name join, for reference: 8,059 retained and 310 eliminated;"
echo "  germline 35.2% / 41.3%, somatic 35.2% / 33.5%, low 29.7% / 25.2%."
grep -m1 '^  this run:' "${SUM}" | sed 's/^  /  /'
echo ""
echo "===== ELIMINATION CALLS: region1 x workbook ====="
sed -n '/^Elimination calls/,/^  Percentages/p' "${SUM}"
echo ""
echo "===== DEFAULT DEFINITION (first floor quantile) ====="
# The first block of the contrast table is the default definition.
sed -n '/^Contrasts (/,/all_expressed/p' "${SUM}"
echo ""
sed -n '/^Matched elimination contrast/,/^$/p' "${SUM}"
sed -n '/^Outcome at the default definition/,/reading:/p' "${SUM}"
echo ""

echo "===== FIGURE (S12 B and C; S12A is functional_enrichment's S6_dnds_vs_expression) ====="
for s in S12B_germline_groups S12C_germline_contrasts; do
  if [[ -s "${OUTDIR}/${s}.svg" ]]; then
    echo "  OK       ${OUTDIR}/${s}.svg"
  else
    echo "  MISSING  ${OUTDIR}/${s}.svg"
  fi
done
sed -n '/^Figure audit/,/^$/p' "${SUM}"
