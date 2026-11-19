#!/bin/bash
#SBATCH -J wb_bridge
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=wb_bridge_%j.out
#SBATCH --error=wb_bridge_%j.err
#
# Links the expression workbook's genes (v37 coordinates) to the v50 miniprot
# models by sequence lift (minimap2) and overlap, to replace the name join that
# results draft note 60 found to be wrong. Writes workbook_model_bridge.tsv
# (model gene_id -> workbook_name) and validates it, beside the old name join,
# against region1 on the analysed genes.
#
#   sbatch --exclude=clr0808 run_workbook_model_bridge.sh
#   bash   run_workbook_model_bridge.sh
#
# minimap2 is taken from MINIMAP2 (full path to the binary), else from PATH
# once the kaks environment is active, else from `module load minimap2`.
# It is not in the kaks environment yet; to add it:
#   conda install -n kaks -c bioconda minimap2
#
# Optional settings (environment or config.sh):
#   AS_V37_GENOME_FA  v37 FASTA the workbook is on (needs its .fai)
#                     default ${DB_DIR}/ascaris/v37/AG.genome.fa
#   AS_GENOME_FA      v50 FASTA the miniprot GFF is on
#                     default ${DB_DIR}/ascaris/v50/AG.v50.fa
#   MINIMAP2          path to the minimap2 binary
#   RUN_TAG           appended to the output folder name
#
# minimap2 runs on the 4 CPUs; the Python steps are single-threaded. Indexing
# v50 and aligning ~15,700 gene spans should take a few minutes (a 100 Mb test
# took 9 s and 0.5 GB in minimap2); 8G leaves room for a ~300 Mb genome.

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
SCRIPT="$(find_script workbook_model_bridge.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: workbook_model_bridge.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
SDIR="$(dirname "${SCRIPT}")"
# The workbook reader and the old name join are imported from the converted
# functional_enrichment.py, which must sit beside the script: the May version
# in the miniprot directory has no target_bridge.
[[ -f "${SDIR}/functional_enrichment.py" ]] || {
  echo "ERROR: functional_enrichment.py is not beside ${SCRIPT}." >&2; exit 1; }
grep -q '^from paper_colors import' "${SDIR}/functional_enrichment.py" || {
  echo "ERROR: ${SDIR}/functional_enrichment.py is an unconverted script (no paper_colors import)." >&2
  exit 1; }

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
missing = [m for m in ("paper_colors", "region_lib", "rate_tables", "openpyxl")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
fe = importlib.import_module("functional_enrichment")
lacking = [n for n in ("target_bridge", "load_workbook") if not hasattr(fe, n)]
if lacking:
    print("ERROR: functional_enrichment lacks " + ", ".join(lacking), file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, rate_tables, openpyxl, functional_enrichment")
PYEOF

# minimap2: MINIMAP2, then PATH (the kaks environment is active now), then the
# module system. `module` is a shell function that batch shells may not have.
MM2="${MINIMAP2:-}"
if [[ -z "${MM2}" ]]; then
  MM2="$(command -v minimap2 || true)"
fi
if [[ -z "${MM2}" ]] && type module >/dev/null 2>&1; then
  module load minimap2 >/dev/null 2>&1 || true
  MM2="$(command -v minimap2 || true)"
fi
if [[ -z "${MM2}" || ! -x "${MM2}" ]]; then
  echo "ERROR: minimap2 not found (MINIMAP2 unset, not on PATH, no module)." >&2
  echo "  Install it into the kaks environment and resubmit:" >&2
  echo "    conda install -n ${CONDA_ENV} -c bioconda minimap2" >&2
  echo "  or point MINIMAP2 at a binary: MINIMAP2=/path/to/minimap2 sbatch run_workbook_model_bridge.sh" >&2
  exit 1
fi
echo "minimap2: ${MM2} ($("${MM2}" --version))"

# Working directory, explicitly: relative paths in config.sh are relative to it.
cd "${PROJECT_DIR}"

TABLE="${ANNOTATED_TABLE:-}"
[[ -n "${TABLE}" ]] || { echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2; exit 1; }
XLSX="${EXPRESSION_XLSX:-Ascaris_genes_15714_annotation_expression.xlsx}"
GFF1="${AS_MINIPROT_GFF:-As_miniprot.gff}"
V37="${AS_V37_GENOME_FA:-${DB_DIR}/ascaris/v37/AG.genome.fa}"
V50="${AS_GENOME_FA:-${DB_DIR}/ascaris/v50/AG.v50.fa}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Seconds in the stamp: a date-only stamp sends two runs on the same day into
# one folder, the second overwriting the first.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/workbook_model_bridge_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-10s %s\n' "Table:" "${TABLE}" "Workbook:" "${XLSX}" "GFF:" "${GFF1}" \
  "v37:" "${V37}" "v50:" "${V50}" "Threads:" "${THREADS}" "Output:" "${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${TABLE}" "${XLSX}" "${GFF1}" "${V37}" "${V37}.fai" "${V50}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done

# The workbook's chromosomes must be v37 sequence names, or nothing extracts.
# A few absent names only cost their genes (reported as not lifted); most
# absent means the wrong assembly.
WB_CHR="$("${PY}" -c "import pandas as pd, sys; print(' '.join(sorted(pd.read_excel(sys.argv[1], usecols=['chr'])['chr'].astype(str).str.strip().unique())))" "${XLSX}")"
N_WB=0; N_IN=0; ABSENT=""
for c in ${WB_CHR}; do
  N_WB=$((N_WB + 1))
  if cut -f1 "${V37}.fai" | grep -qx -F "$c"; then
    N_IN=$((N_IN + 1))
  else
    ABSENT="${ABSENT} $c"
  fi
done
printf '  %-28s %6d / %6d in %s\n' "workbook chromosomes" "${N_IN}" "${N_WB}" "$(basename "${V37}.fai")"
[[ -z "${ABSENT}" ]] || echo "  not in the v37 index (their genes will not lift):${ABSENT}"
[[ "${N_IN}" -gt $(( N_WB / 2 )) ]] || {
  echo "ERROR: most workbook chromosomes are not in ${V37}.fai; is the workbook on this assembly?" >&2; exit 1; }

# The rate table's gene1 IDs must be the GFF's mRNA IDs, or validation joins
# nothing.
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
echo ""

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --annotation "${XLSX}" \
  --old-genome "${V37}" \
  --new-genome "${V50}" \
  --gff1 "${GFF1}" \
  --kaks "${TABLE}" \
  --minimap2 "${MM2}" \
  --threads "${THREADS}" \
  --outdir "${OUTDIR}"

SUM="${OUTDIR}/workbook_model_bridge_summary.txt"
echo ""
echo "===== DONE ====="
echo "Bridge:  ${OUTDIR}/workbook_model_bridge.tsv"
echo "Summary: ${SUM}"
for f in workbook_model_bridge.tsv workbook_lift.tsv bridge_validation.tsv; do
  if [[ -s "${OUTDIR}/${f}" ]]; then
    echo "  OK       ${f}"
  else
    echo "  MISSING  ${f}"
  fi
done
