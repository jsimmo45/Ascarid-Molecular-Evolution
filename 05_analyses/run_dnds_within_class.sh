#!/bin/bash
#SBATCH -J dnds_class
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=dnds_class_%j.out
#SBATCH --error=dnds_class_%j.err
#
# Reviewer control for Figure 6C (2026-10-19, results draft note 70): does
# functional class account for the dN/dS excess of eliminated genes? Two of the
# 73 classes of 6C (Unknown, Proteasome) are enriched among eliminated genes
# after Bonferroni correction. This tests the excess within classes:
#   - eliminated / retained fold of medians, with and without the Bonferroni
#     classes and the unknown-type classes;
#   - each eliminated gene matched to a retained gene of its own class;
#   - van Elteren rank-sum test stratified by class;
#   - the Figure 6E model form with class fixed effects;
#   - per-class folds for the classes drawn in 6C;
#   - retained genes only: does each class evolve faster on its own?
# Draws Supplementary Figure S12D (S12D_function_classes: the folds by gene set
# and by class), 89 mm wide and as tall as the newest S12C, so C | D make the
# bottom row of S12. Run run_germline_restricted_dnds.sh first. Writes
# results/dnds_within_class_<stamp>/ (figure, summary, tables, provenance).
#
#   jid=$(sbatch --parsable --exclude=clr0808 run_germline_restricted_dnds.sh)
#   sbatch --exclude=clr0808 --dependency=afterok:${jid} run_dnds_within_class.sh
#   bash   run_dnds_within_class.sh
#
# Optional settings (environment or config.sh):
#   UNKNOWN_PATTERN   regex for the unknown-type classes, any case
#                     (default 'unknown|uncharacteri[sz]ed')
#   MIN_ELIM_CLASS    per-class contrasts for classes with at least this many
#                     eliminated genes (default 10, the 6C forest set)
#   D_WIDTH           S12D width: single (default) | onehalf | double
#   D_HEIGHT_MM       S12D height; default: the height of the newest S12C
#                     (results/germline_restricted_dnds_*/S12C_germline_contrasts.svg)
#                     when that S12C is 89 mm wide, else S12_BOTTOM_MM
#   S12_BOTTOM_MM     height of the S12 bottom row (default 80, as
#                     run_germline_restricted_dnds.sh)
#   RUN_TAG           appended to the output folder name
#
# The join, the class column and the class tests are imported from
# functional_enrichment.py, the fold and matching from
# germline_restricted_dnds.py, the model from multivariate_insulation.py and
# the bridge reader from workbook_model_bridge.py, so all four must sit beside
# this script (the 2026-10-19 versions). The expression workbook is joined by
# position through workbook_model_bridge.tsv, as in Figure 6: WB_BRIDGE in
# config.sh, or else the newest results/workbook_model_bridge_*/.
#
# Single-threaded pandas/scipy: 1 node / 1 task / 1 CPU. About half a minute on
# a table of the real size, figure included; 4G covers the workbook and the
# annotated table.

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
SCRIPT="$(find_script dnds_within_class.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: dnds_within_class.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
grep -q -- "'--unknown-pattern'" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is not the 2026-10-19 version (no --unknown-pattern)." >&2
  exit 1; }
grep -q -- "'--d-width'" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} draws no figure (the version from before S12D, 2026-10-19)." >&2
  exit 1; }
SDIR="$(dirname "${SCRIPT}")"
# The imported modules must sit beside the script and be the converted ones:
# the May functional_enrichment.py in the miniprot directory has no
# category_tests and would fail on import, or worse, be found first.
for m in functional_enrichment.py germline_restricted_dnds.py \
         multivariate_insulation.py workbook_model_bridge.py; do
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
# script without the functions used here, is what a file test misses.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix the environment (lib/ on PYTHONPATH, openpyxl installed, converted 05_analyses scripts) and resubmit." >&2; exit 1; }
import importlib, importlib.util, sys
missing = [m for m in ("paper_colors", "region_lib", "rate_tables", "plain_text",
                       "openpyxl")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
need = {"functional_enrichment": ("load_workbook", "category_tests", "CLASS_COL"),
        "germline_restricted_dnds": ("ratio_contrast", "match_once"),
        "multivariate_insulation": ("fit_model",),
        "workbook_model_bridge": ("read_bridge", "shared_gene_counts")}
lacking = []
for mod, names in need.items():
    m = importlib.import_module(mod)
    lacking += [f"{mod}.{n}" for n in names if not hasattr(m, n)]
if lacking:
    print("ERROR: missing in the imported scripts: " + ", ".join(lacking), file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, plain_text, openpyxl, "
      "functional_enrichment, germline_restricted_dnds, multivariate_insulation, "
      "workbook_model_bridge")
PYEOF

# Working directory, explicitly: relative paths in config.sh are relative to it.
cd "${PROJECT_DIR}"

TABLE="${ANNOTATED_TABLE:-}"
[[ -n "${TABLE}" ]] || { echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2; exit 1; }
XLSX="${EXPRESSION_XLSX:-Ascaris_genes_15714_annotation_expression.xlsx}"
# The position bridge: WB_BRIDGE in config.sh, else the newest bridge run.
BRIDGE="${WB_BRIDGE:-}"
if [[ -z "${BRIDGE}" ]]; then
  LATEST="$(ls -1d "${RESULTS_DIR}"/workbook_model_bridge_*/ 2>/dev/null | sort | tail -1 || true)"
  [[ -n "${LATEST}" ]] && BRIDGE="${LATEST%/}/workbook_model_bridge.tsv"
fi
[[ -n "${BRIDGE}" && -s "${BRIDGE}" ]] || {
  echo "ERROR: no workbook_model_bridge.tsv (WB_BRIDGE unset and none under ${RESULTS_DIR}/)." >&2
  echo "  Run run_workbook_model_bridge.sh first, or set WB_BRIDGE in config.sh." >&2
  exit 1; }
UNKNOWN_PATTERN="${UNKNOWN_PATTERN:-unknown|uncharacteri[sz]ed}"
MIN_ELIM_CLASS="${MIN_ELIM_CLASS:-10}"
D_WIDTH="${D_WIDTH:-single}"

# S12D shares the bottom row with S12C, so it is drawn at the height of the
# newest S12C (matplotlib writes the SVG size in points). Only an 89 mm S12C
# counts: a wider one is the layout from before C | D shared the row, and its
# height says nothing about the row.
D_SRC="D_HEIGHT_MM"
D_HEIGHT="${D_HEIGHT_MM:-}"
if [[ -z "${D_HEIGHT}" ]]; then
  GR_LATEST="$(ls -1d "${RESULTS_DIR}"/germline_restricted_dnds_*/ 2>/dev/null | sort | tail -1 || true)"
  C_SVG="${GR_LATEST%/}/S12C_germline_contrasts.svg"
  if [[ -n "${GR_LATEST}" && -s "${C_SVG}" ]]; then
    C_SIZE="$("${PY}" - "${C_SVG}" <<'PYEOF'
import re, sys
head = open(sys.argv[1], encoding='utf-8', errors='replace').read(4000)
m = re.search(r'<svg\b([^>]*)>', head, re.S)
attrs = m.group(1) if m else ''
w = re.search(r'\swidth="([0-9.]+)pt"', attrs)
h = re.search(r'\sheight="([0-9.]+)pt"', attrs)
print(f'{float(w.group(1)) * 25.4 / 72:.1f} {float(h.group(1)) * 25.4 / 72:.1f}'
      if w and h else '')
PYEOF
)" || C_SIZE=""
    C_W="${C_SIZE% *}"
    C_H="${C_SIZE#* }"
    if [[ -z "${C_SIZE}" ]]; then
      echo "NOTE: could not read the size of ${C_SVG}; S12D uses S12_BOTTOM_MM."
    elif awk -v w="${C_W}" 'BEGIN { exit !(w > 88.5 && w < 89.5) }'; then
      D_HEIGHT="${C_H}"
      D_SRC="S12C, ${C_SVG}"
    else
      echo "NOTE: the newest S12C is ${C_W} mm wide, not 89 (${C_SVG}): rerun"
      echo "      run_germline_restricted_dnds.sh for the C | D row. S12D uses S12_BOTTOM_MM."
    fi
  fi
fi
if [[ -z "${D_HEIGHT}" && -n "${S12_BOTTOM_MM:-}" ]]; then
  D_HEIGHT="${S12_BOTTOM_MM}"; D_SRC="S12_BOTTOM_MM"
elif [[ -z "${D_HEIGHT}" ]]; then
  D_HEIGHT=80; D_SRC="S12_BOTTOM_MM default"
fi

# Seconds in the stamp: a date-only stamp sends two runs on the same day into
# one folder, the second overwriting the first.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/dnds_within_class_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-17s %s\n' "Table:" "${TABLE}" "Workbook:" "${XLSX}" "Bridge:" "${BRIDGE}" \
  "Unknown pattern:" "${UNKNOWN_PATTERN}" "Min elim/class:" "${MIN_ELIM_CLASS}" \
  "S12D size:" "${D_WIDTH}, ${D_HEIGHT} mm (${D_SRC})" \
  "Output:" "${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${TABLE}" "${XLSX}" "${BRIDGE}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done

# The rate table's gene1 IDs must be the bridge's model IDs (the GFF's mRNA
# IDs), or nothing joins and every count below is zero.
COL="$(head -1 "${TABLE}" | tr '\t' '\n' | grep -n -x -F gene1 | cut -d: -f1 | head -1)"
[[ -n "${COL}" ]] || { echo "ERROR: no gene1 column in ${TABLE}" >&2; exit 1; }
MATCH="$(awk -v col="${COL}" -F'\t' '
  NR==FNR { if (FNR == 1) { for (i = 1; i <= NF; i++) if ($i == "gene_id") g = i; next }
            ids[$g] = 1; next }
  FNR==1 { next }
  { if ($col in ids) n++ }
  END { print n+0 }' "${BRIDGE}" "${TABLE}")"
ROWS=$(( $(wc -l < "${TABLE}") - 1 ))
printf '  %-28s %6d / %6d table rows found in the bridge\n' "gene1 identifiers" "${MATCH}" "${ROWS}"
[[ "${MATCH}" -gt $(( ROWS / 2 )) ]] || {
  echo "ERROR: fewer than half the gene1 identifiers appear in ${BRIDGE}." >&2; exit 1; }

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --kaks "${TABLE}" \
  --annotation "${XLSX}" \
  --bridge "${BRIDGE}" \
  --unknown-pattern "${UNKNOWN_PATTERN}" \
  --min-elim-class "${MIN_ELIM_CLASS}" \
  --d-width "${D_WIDTH}" \
  --d-height-mm "${D_HEIGHT}" \
  --seed 0 \
  --outdir "${OUTDIR}"

SUM="${OUTDIR}/dnds_within_class_summary.txt"
echo ""
echo "===== DONE ====="
echo "Summary: ${SUM}"
echo ""

echo "===== CHECK AGAINST FIGURE 6C ====="
echo "The class tests must be Figure 6C's (73 classes, Unknown 0.0020, Proteasome 0.017 on 2026-10-19)."
FE_LATEST="$(ls -1d "${RESULTS_DIR}"/functional_enrichment_*/ 2>/dev/null | sort | tail -1 || true)"
if [[ -n "${FE_LATEST}" && -s "${FE_LATEST%/}/functional_enrichment_summary.txt" ]]; then
  echo "  6C run (${FE_LATEST%/}):"
  grep -m1 'classes tested (>=' "${FE_LATEST%/}/functional_enrichment_summary.txt" | sed 's/^ */    /'
else
  echo "  no functional_enrichment run found under ${RESULTS_DIR}/ to compare"
fi
echo "  this run:"
sed -n '/^Class tests (as Figure 6C)/,/^Unknown-type classes/p' "${SUM}" | sed '$d' | sed 's/^/  /'
echo ""

echo "===== CONTRASTS ====="
sed -n '/^Eliminated \/ retained dN\/dS/,/^$/p' "${SUM}"
sed -n '/^Class-matched/,/^$/p' "${SUM}"
grep '^van Elteren' "${SUM}" || true
echo ""
sed -n '/^Models (standardized OLS/,/^$/p' "${SUM}"
sed -n '/^Per class, classes with/,/^$/p' "${SUM}"
grep '^Outcome' "${SUM}" || true
echo ""

echo "===== FIGURE (S12D, beside S12C on the bottom row) ====="
if [[ -s "${OUTDIR}/S12D_function_classes.svg" ]]; then
  echo "  OK       ${OUTDIR}/S12D_function_classes.svg"
else
  echo "  MISSING  ${OUTDIR}/S12D_function_classes.svg"
fi
echo "  height ${D_HEIGHT} mm, from ${D_SRC}"
sed -n '/^Figure audit/,/^$/p' "${SUM}"
sed -n '/^Legend notes (S12D/,/^$/p' "${SUM}"
