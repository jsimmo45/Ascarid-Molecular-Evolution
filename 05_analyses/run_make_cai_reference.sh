#!/bin/bash
#SBATCH -J cai_ref
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=cai_ref_%j.out
#SBATCH --error=cai_ref_%j.err
#
# Builds the Ascaris CAI reference (the most highly expressed workbook genes in
# somatic tissue) and translates it into v50 model IDs through the position
# bridge, replacing the reference built through the GFF Target= names (results
# draft note 62). Writes results/cai_reference_<stamp>/:
#   cai_reference_ascaris.txt     one model ID per line (codon_usage_bias.py --ref-genes)
#   cai_reference_summary.txt     the script's report
#   provenance.json               inputs with checksums, settings, counts
#
#   sbatch --exclude=clr0808 run_make_cai_reference.sh
#
# Then point CAI_REFERENCE in config.sh at the new list (the job prints the
# line) and rerun:  sbatch --exclude=clr0808 run_codon_usage_bias.sh ascaris
#
# Settings (environment or config.sh):
#   EXPRESSION_XLSX   workbook (default Ascaris_genes_15714_annotation_expression.xlsx)
#   WB_BRIDGE         workbook_model_bridge.tsv (default: the newest
#                     results/workbook_model_bridge_*/)
#   ANNOTATED_TABLE   rate table whose gene1 IDs the reference must match (required)
#   CAI_COLUMNS       somatic columns to pool (default "muscle intestine carcass")
#   CAI_TOP_PERCENT   share of genes kept, highest first (default 5)
#   CAI_EXCLUDE_ELIMINATED  1 (default): rank only the genes the workbook calls
#                     retained (its 'eliminated' column is '.'), and drop any
#                     reference model that region1 calls eliminated. Four
#                     workbook-eliminated genes (an MSP and three hypothetical
#                     proteins) reach the somatic top 5% although their DNA is
#                     absent from somatic cells, most likely through reads from
#                     retained paralogues; and one retained gene's model
#                     (MP014343, 2026-10-17) sits where the v50 call is
#                     eliminated. 0 ranks every gene and drops nothing, as the
#                     first reference did.
#   CAI_REFERENCE     the reference in use now; only read, to compare with the new one
#   RUN_TAG           appended to the output folder name
#
# Reads one 15,714-row workbook and two ~10,000-row tables; seconds.
# Single-threaded pandas: 1 node / 1 task / 1 CPU.

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
SCRIPT="$(find_script make_cai_reference.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: make_cai_reference.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# An older copy without --bridge would fall back to the name join.
grep -q -- "'--bridge'" "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} has no --bridge option (the version from before 2026-10-17)." >&2; exit 1; }
SDIR="$(dirname "${SCRIPT}")"
# The bridge reader is imported from workbook_model_bridge.py, which must sit
# beside the script.
[[ -f "${SDIR}/workbook_model_bridge.py" ]] || {
  echo "ERROR: workbook_model_bridge.py is not beside ${SCRIPT}." >&2; exit 1; }

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
# workbook_model_bridge.py without read_bridge, is what a file test misses.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix the environment (lib/ on PYTHONPATH, openpyxl installed, the 2026-10-17 workbook_model_bridge.py) and resubmit." >&2; exit 1; }
import importlib, importlib.util, sys
missing = [m for m in ("region_lib", "rate_tables", "openpyxl")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
wmb = importlib.import_module("workbook_model_bridge")
if not hasattr(wmb, "read_bridge"):
    print("ERROR: workbook_model_bridge lacks read_bridge", file=sys.stderr)
    sys.exit(1)
print("modules: region_lib, rate_tables, openpyxl, workbook_model_bridge all importable")
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
if [[ -z "${BRIDGE}" ]]; then
  echo "ERROR: no workbook_model_bridge.tsv (WB_BRIDGE unset and none under ${RESULTS_DIR}/)." >&2
  echo "  Run run_workbook_model_bridge.sh first, or set WB_BRIDGE in config.sh." >&2
  exit 1
fi
read -r -a COLUMNS <<< "${CAI_COLUMNS:-muscle intestine carcass}"
TOP="${CAI_TOP_PERCENT:-5}"
EXCLUDE_ARG=()
case "${CAI_EXCLUDE_ELIMINATED:-1}" in
  1) EXCLUDE_ARG=(--exclude-col eliminated --exclude-keep . --drop-target-region eliminated) ;;
  0) ;;
  *) echo "ERROR: CAI_EXCLUDE_ELIMINATED must be 1 or 0, not '${CAI_EXCLUDE_ELIMINATED}'." >&2; exit 1 ;;
esac
OLD_REF="${CAI_REFERENCE:-cai_reference_ascaris.txt}"

# Seconds in the stamp: a date-only stamp sends two runs on the same day into
# one folder, the second overwriting the first.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/cai_reference_${STAMP}${RUN_TAG:+_${RUN_TAG}}"
OUT="${OUTDIR}/cai_reference_ascaris.txt"

printf '%-13s %s\n' "Workbook:" "${XLSX}" "Bridge:" "${BRIDGE}" "Table:" "${TABLE}" \
  "Columns:" "${COLUMNS[*]} (mean)" "Top percent:" "${TOP}" \
  "Ranked:" "$([[ ${#EXCLUDE_ARG[@]} -gt 0 ]] && echo "workbook-retained genes only; models region1 calls eliminated dropped" || echo "every gene; nothing dropped")" \
  "In use now:" "${OLD_REF}" "Output:" "${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${XLSX}" "${BRIDGE}" "${TABLE}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done
# The rate table's gene1 IDs must be the bridge's model IDs, or nothing joins.
COL="$(head -1 "${TABLE}" | tr '\t' '\n' | grep -n -x -F gene1 | cut -d: -f1 | head -1)"
[[ -n "${COL}" ]] || { echo "ERROR: no gene1 column in ${TABLE}" >&2; exit 1; }
BCOL="$(head -1 "${BRIDGE}" | tr '\t' '\n' | grep -n -x -F gene_id | cut -d: -f1 | head -1)"
[[ -n "${BCOL}" ]] || { echo "ERROR: no gene_id column in ${BRIDGE}; is it workbook_model_bridge.tsv?" >&2; exit 1; }
MATCH="$(awk -v bcol="${BCOL}" -v col="${COL}" -F'\t' '
  NR==FNR { if (FNR > 1) ids[$bcol]=1; next }
  FNR==1 { next }
  { if ($col in ids) n++ }
  END { print n+0 }' "${BRIDGE}" "${TABLE}")"
ROWS=$(( $(wc -l < "${TABLE}") - 1 ))
printf '  %-28s %6d / %6d table rows found in the bridge\n' "gene1 identifiers" "${MATCH}" "${ROWS}"
[[ "${MATCH}" -gt $(( ROWS / 2 )) ]] || {
  echo "ERROR: fewer than half the gene1 identifiers are bridge model IDs; wrong bridge or table?" >&2; exit 1; }
echo ""

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --expression "${XLSX}" \
  --gene-col name \
  --columns "${COLUMNS[@]}" \
  --statistic mean \
  --top-percent "${TOP}" \
  "${EXCLUDE_ARG[@]}" \
  --bridge "${BRIDGE}" \
  --target-table "${TABLE}" \
  --target-col gene1 \
  --target-region-col region1 \
  --out "${OUT}" \
  --provenance "${OUTDIR}/provenance.json" \
  | tee "${OUTDIR}/cai_reference_summary.txt"

[[ -s "${OUT}" ]] || { echo "ERROR: ${OUT} was not written." >&2; exit 1; }

echo ""
echo "===== DONE ====="
echo "Reference: ${OUT} ($(wc -l < "${OUT}") model IDs)"
if [[ -s "${OLD_REF}" ]]; then
  SHARED="$(comm -12 <(sort -u "${OLD_REF}") <(sort -u "${OUT}") | wc -l)"
  echo "In use now: ${OLD_REF} ($(wc -l < "${OLD_REF}") IDs); ${SHARED} shared with the new list."
else
  echo "In use now: ${OLD_REF} not found (nothing to compare)."
fi
echo ""
echo "Next: set this line in config.sh, then rerun the Ascaris codon-usage job."
echo "  CAI_REFERENCE=\"\${RESULTS_DIR}/$(basename "${OUTDIR}")/cai_reference_ascaris.txt\""
echo "  sbatch --exclude=clr0808 run_codon_usage_bias.sh ascaris"
echo "The codon-usage job prints 'CAI ref:' at the top; check it names this file."
