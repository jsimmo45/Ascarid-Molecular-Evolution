#!/bin/bash
#SBATCH -J tox_annot_diag
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=tox_annot_diag_%j.out
#SBATCH --error=tox_annot_diag_%j.err
#
# Supplementary Figure S13 and Table T7: is the Toxocara rate elevation an
# artefact of TransDecoder gene models? Three tests (premise, restriction to
# clean gene models, covariate regression with a four-species control); see the
# script's docstring. Since 2026-10-10 the restriction also reports dS, which
# decides how Results words the small Toxocara dS excess.
#
#   sbatch run_toxocara_annotation_diagnostic.sh
#   sbatch run_toxocara_annotation_diagnostic.sh <path/to/branch_dnds_four_annotated.tsv>
#
# The job runs --inspect first and then the real analysis, so the ID match rates
# are in the job log as a record of what was joined to what. A match rate near
# zero is an ID mismatch, not a result: fix it before reading anything else.
#
# INPUT TABLE: argument 1, else ANNOTATED_FOUR from config.sh, else the newest
# timestamped branch_dnds_four_<date>_<time>/ under RESULTS_DIR. It must come
# from a 2026-10-10 or later run (unrooted tree); the script refuses older
# tables. Use results/branch_dnds_four_20261010_113902.
#
# EVIDENCE FILES (relative to PROJECT_DIR, the miniprot directory), overridable
# in config.sh:
#   TOX_TRANSDECODER_PEP  TransDecoder peptides; carries the ORF types, the most
#                         informative covariate. Default below.
#   TOX_DIAG_RBH          space-separated RBH tables involving Toxocara.
# Per-gene alignments are not kept by branch_dnds_four.py; align_len from the
# annotated table is used instead (--table-covariates auto).
#
# Replaces the run_tox_annotation_diagnostic.sh in the miniprot directory,
# which hardcoded the working directory and the pre-2026-10-10 v2 table.
#
# Single-threaded pandas/scipy: 1 node / 1 task / 1 CPU. 8G covers the
# TransDecoder FASTA and the RBH tables.

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
SCRIPT="$(find_script toxocara_annotation_diagnostic.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: toxocara_annotation_diagnostic.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# The unconverted script of the same name sits in the miniprot directory, which
# is SLURM_SUBMIT_DIR if this is submitted from there. Only the converted one
# imports branch_lib, so refuse anything that does not.
grep -q '^from branch_lib import' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the old, unconverted script (no branch_lib import)." >&2
  echo "  Use Ascarid-Molecular-Evolution/05_analyses/toxocara_annotation_diagnostic.py." >&2
  exit 1; }

export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

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

# Probe the imports, and that lib/branch_lib.py is the 2026-10-10 (afternoon)
# version, which has the margin helper this script imports.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) or update lib/branch_lib.py, and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("paper_colors", "region_lib", "branch_lib")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
import branch_lib
need = [f for f in ("usable_mask", "load_annotated", "fit_bottom_margin", "compare")
        if not hasattr(branch_lib, f)]
if need:
    print("ERROR: lib/branch_lib.py lacks " + ", ".join(need) + " (older version)",
          file=sys.stderr)
    sys.exit(1)
print("modules: paper_colors, region_lib, branch_lib (current) all importable")
PYEOF

cd "${PROJECT_DIR}"

ANNOTATED="${1:-${ANNOTATED_FOUR:-}}"
if [[ -z "${ANNOTATED}" ]]; then
  ANNOTATED="$(ls -1dt "${RESULTS_DIR}"/branch_dnds_four_[0-9]*/branch_dnds_four_annotated.tsv \
               2>/dev/null | head -n1 || true)"
fi
if [[ -z "${ANNOTATED}" || ! -f "${ANNOTATED}" ]]; then
  echo "ERROR: no annotated table. Pass one as argument 1, set ANNOTATED_FOUR" >&2
  echo "       in config.sh, or run run_branch_dnds_four.sh first." >&2
  exit 1
fi
echo "input : ${ANNOTATED} ($(( $(wc -l < "${ANNOTATED}") - 1 )) genes)"
PROV="$(dirname "${ANNOTATED}")/provenance.json"
if [[ -f "${PROV}" ]]; then
  "${PY}" - "${PROV}" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
print("run   : %s  commit %s  quartets %s -> ok %s"
      % (d.get("run_utc", "?"), d.get("git_commit", "?"),
         d.get("n_quartets", "?"), d.get("n_ok", "?")))
EOF
fi

TRANSDECODER="${TOX_TRANSDECODER_PEP:-TG.transcripts.longest.20065.fasta.transdecoder.pep}"
read -r -a RBH_FILES <<< "${TOX_DIAG_RBH:-rbh_As_Tox_v2/rbh.tsv rbh_Pu_Tox_v2/rbh.tsv rbh_Tox_Bsc_v2/rbh.tsv}"

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/toxocara_annotation_diagnostic_${STAMP}"

ARGS=(--annotated "${ANNOTATED}" --focal-species Toxocara --outdir "${OUTDIR}")
# The TransDecoder FASTA is what the job is for: without it there is no ORF
# completeness, no complete-ORF stratum and nothing Results quotes. The first
# 2026-10-10 run printed a NOTE to the .err file and carried on with alignment
# length alone, so it now stops unless ALLOW_NO_TRANSDECODER=1.
if [[ -f "${TRANSDECODER}" ]]; then
  ARGS+=(--transdecoder "${TRANSDECODER}")
  echo "transdecoder: ${TRANSDECODER}"
elif [[ "${ALLOW_NO_TRANSDECODER:-0}" == "1" ]]; then
  echo "WARNING: no TransDecoder FASTA ('${TRANSDECODER}'); running on alignment" >&2
  echo "         length alone because ALLOW_NO_TRANSDECODER=1." >&2
else
  echo "ERROR: TransDecoder FASTA not found: ${PROJECT_DIR}/${TRANSDECODER}" >&2
  echo "  Find it with:" >&2
  echo "    find ${PROJECT_DIR} -maxdepth 4 -name '*transdecoder.pep'" >&2
  echo "  then set TOX_TRANSDECODER_PEP (full path) in config.sh or on the" >&2
  echo "  sbatch command line and resubmit." >&2
  exit 1
fi
EXISTING_RBH=()
for f in "${RBH_FILES[@]}"; do
  if [[ -f "$f" ]]; then EXISTING_RBH+=("$f"); else echo "WARNING: RBH table not found: $f (coverage stratum needs it)"; fi
done
if [[ ${#EXISTING_RBH[@]} -gt 0 ]]; then
  ARGS+=(--rbh "${EXISTING_RBH[@]}")
  echo "rbh   : ${EXISTING_RBH[*]}"
fi
echo "outdir: ${OUTDIR}"

echo "=============================================================="
echo "STEP 1/2  --inspect (ID matching only, no computation)"
echo "=============================================================="
"${PY}" "${SCRIPT}" "${ARGS[@]}" --inspect

echo
echo "=============================================================="
echo "STEP 2/2  full diagnostic"
echo "=============================================================="
"${PY}" "${SCRIPT}" "${ARGS[@]}" \
  --omega-cap 5 \
  --min-coverage 0.7 \
  --min-stratum-n 20 \
  --table-covariates auto

echo
echo "===== DONE ====="
sed -n '/^Numbers for the text/,/^TEST 1/p' "${OUTDIR}/diagnostic_summary.txt" | sed '$d'
echo "Full summary: ${OUTDIR}/diagnostic_summary.txt (TEST 1 first)."
