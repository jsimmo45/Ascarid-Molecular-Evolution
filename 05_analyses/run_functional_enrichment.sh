#!/bin/bash
#SBATCH -J func_enrich
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=func_enrich_%j.out
#SBATCH --error=func_enrich_%j.err
#
# Figure 6A + 6B (one 185 mm strip), Figure 6C, Supplementary Figure S12A
# (written as S6_dnds_vs_expression; run_germline_restricted_dnds.sh draws
# S12B at its height, so run that one afterok this).
#
# 2026-10-19 (results draft note 67): S12A draws the genes with zero mean
# expression at log2(0.1) instead of dropping them, so its n matches 6A's; the
# 2 of its log2 is a 6 pt subscript drawn by lib/plain_text.py; 6C gains a
# column of Bonferroni-adjusted P per class. Tests and tables are unchanged.
#
#   sbatch run_functional_enrichment.sh
#   bash   run_functional_enrichment.sh
#
# 2026-10-17: the expression workbook is joined by position, through
# workbook_model_bridge.tsv (run_workbook_model_bridge.sh), not by the GFF
# Target= name, which matched each model to a different gene (results draft
# note 60). Every number in Figure 6A-C and S6 changes. The bridge is WB_BRIDGE
# in config.sh, or else the newest results/workbook_model_bridge_*/ folder.
# The 2026-10-15 hand edits are now in the code: 6C class names in sentence
# case, S6 axis "log2", P exponents as raised 6 pt text (lib/plain_text.py).
#
# Earlier changes from the May script: denominators are genes the workbook can
# classify; early embryo counts as germline; 6C plots odds ratios with exact
# 95% CIs for classes with >= 10 eliminated genes.
#
# The DONE block prints the values the Results draft quotes (computed with the
# superseded name join) beside this run's, so the text can be updated.
#
# Single-threaded pandas/scipy: 1 node / 1 task / 1 CPU. 8G covers the
# workbook, the bridge and the annotated table.

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
SCRIPT="$(find_script functional_enrichment.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: functional_enrichment.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# The May script of the same name still sits in the miniprot directory, which
# is where SLURM_SUBMIT_DIR points if this is submitted from there. Only the
# converted one imports paper_colors, so refuse anything that does not.
grep -q '^from paper_colors import' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the old, unconverted script (no paper_colors import)." >&2
  echo "  Use Ascarid-Molecular-Evolution/05_analyses/functional_enrichment.py." >&2
  exit 1; }
grep -q -- "--bridge" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} predates the position bridge (no --bridge option)." >&2
  exit 1; }
grep -q "replace_xlabel" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the version from before 2026-10-19 (S12A drops zero-expression genes; 6C has no P column)." >&2
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

# Probe the imports, not for sibling files: paper_colors, region_lib,
# rate_tables and plain_text live in lib/, workbook_model_bridge beside the
# script, and openpyxl is what pandas needs to read the workbook. A module that
# exists but fails to import is what a file test misses.
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
pt = importlib.import_module("plain_text")
if not all(hasattr(pt, n) for n in ("replace_xlabel", "subscript_runs", "draw_runs")):
    print("ERROR: lib/plain_text.py has no replace_xlabel(); copy the "
          "2026-10-19 version", file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, plain_text, "
      "workbook_model_bridge, openpyxl all importable")
PYEOF

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

# Seconds in the stamp: a date-only stamp sends two runs on the same day into
# one folder, the second overwriting the first.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/functional_enrichment_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-10s %s\n' "Table:" "${TABLE}" "Workbook:" "${XLSX}" "Bridge:" "${BRIDGE}" \
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
  --min-class-genes 5 \
  --forest-min-elim 10 \
  --seed 0 \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
echo "Summary: ${OUTDIR}/functional_enrichment_summary.txt"
echo ""

# The 2026-10-17 bridge run (results draft note 63, job 6427845) beside this
# one. Only S12A's n and rho (zero-expression genes now drawn) may move; every
# test and table is expected to be identical.
cols() {                            # cols <tsv> <col,col,...>: named columns
  awk -F'\t' -v want="$2" '
    NR == 1 { n = split(want, w, ","); for (i = 1; i <= NF; i++) h[$i] = i
              for (j = 1; j <= n; j++) printf "%s%s", w[j], (j < n ? "\t" : "\n"); next }
    { for (j = 1; j <= n; j++) printf "%s%s", $(h[w[j]]), (j < n ? "\t" : "\n") }' "$1"
}
SUM="${OUTDIR}/functional_enrichment_summary.txt"
echo "===== JOIN ====="
grep -m1 '^Join ' "${SUM}" | sed 's/^/  /'
sed -n '/^Genes at each step/,/^$/p' "${SUM}" | sed 's/^/  /'
grep -m1 'Sharing their workbook gene' "${SUM}" | sed 's/^ */  /'
echo ""
echo "===== 2026-10-17 RUN (NOTE 63) vs THIS RUN ====="
echo "10-17: 8,930 retained and 315 eliminated genes with a peak stage."
grep -m1 'with a stage of peak expression' "${SUM}" | sed 's/^/  this run: /'
echo ""
echo "10-17: germline 74.9% vs 41.6%, OR 4.19, P = 4.9e-32; somatic 2.5% vs 44.7%,"
echo "       P = 3.1e-65; low 22.5% vs 13.7%, P = 3.1e-05; chi-square 221.4, df 2."
cols "${OUTDIR}/expression_category_tests.tsv" \
  category,pct_eliminated,pct_retained,odds_ratio,p,q_bh,chi2_homogeneity,chi2_p | sed 's/^/  /'
echo ""
echo "10-17: testis 139 / 975, 4.04-fold, P = 4.6e-48, Bonferroni 4.6e-47."
cols "${OUTDIR}/expression_stage_tests.tsv" \
  category,n_eliminated_in,n_retained_in,fold_enrichment,p,p_bonferroni | sed 's/^/  /'
echo ""
echo "10-17: 73 classes; Bonferroni P Unknown 0.0020, Proteasome 0.017; eight drawn."
echo "       This run's 6C P column is listed under 'P column of 6C'."
sed -n '/^Functional class (6C)/,/^$/p' "${SUM}" | sed 's/^/  /'
echo "10-17: S12 rho = -0.176, P = 5.0e-65, n = 8,925 retained / 315 eliminated"
echo "       (zero-expression genes left out); this run draws them at log2(0.1)."
sed -n '/^dN\/dS vs mean expression/,/^  6A counts/p' "${SUM}" | sed 's/^/  /'
echo ""

echo "===== FIGURES ====="
for f in 6AB_expression_stage 6C_functional_class_forest S6_dnds_vs_expression; do  # S6_ = S12A
  if [[ -s "${OUTDIR}/${f}.svg" ]]; then echo "  OK       ${OUTDIR}/${f}.svg"
  else echo "  MISSING  ${OUTDIR}/${f}.svg"; fi
done
sed -n '/^Figure audit/,/^$/p' "${OUTDIR}/functional_enrichment_summary.txt"
