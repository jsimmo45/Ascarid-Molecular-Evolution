#!/bin/bash
#SBATCH -J branch_dnds4
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=branch_dnds4_%j.out
#SBATCH --error=branch_dnds4_%j.err
#
# Figure 5A, Figure 5B, Table T2, and the annotated table that
# branch_dnds_decomposition.py reads for 5C and 5D. With --tox-eliminated it
# also emits the S14 versions.
#
#   sbatch run_branch_dnds_four.sh
#   bash   run_branch_dnds_four.sh
#
# RESOURCES: the old header asked for --ntasks-per-node=48, which requests 48
# SLURM TASKS. This is ONE python process using a ProcessPoolExecutor, so what
# it needs is 48 CPUs inside one task: --nodes=1 --ntasks=1 --cpus-per-task=48.
# The old form happened to allocate roughly the right core count, but it is the
# wrong idiom and would launch 48 copies of the script under srun.
#
# codeml is the cost: one free-ratio fit per quartet, a few thousand quartets.
# THREADS defaults to the allocated CPU count rather than a hardcoded 48, so a
# smaller allocation does not leave the pool oversubscribed.
#
# 2026-10-09: the script now lives in 03_rates/, which find_script did not
# search. Submitted from 05_analyses it stopped with "not found"; submitted from
# the miniprot directory it would have run the OLD, unconverted
# branch_dnds_four.py sitting there. 03_rates is now searched first, and a
# script without the branch_lib import is refused. --time cut from 24 h to 3 h:
# the v2 run took 1.5 h, most of it pickling whole FASTA dictionaries into every
# task, which the script no longer does.
#
# 2026-10-10: the job refuses to start unless lib/branch_lib.py writes the tree
# unrooted (the earlier rooted tree split the Toxocara branch in two). The
# 2026-10-09 run took 29 min on 48 CPUs.

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
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/03_rates" \
           "${REPO_ROOT}/05_analyses" "${REPO_ROOT}/06_figures" "${REPO_ROOT}" \
           "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script branch_dnds_four.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: branch_dnds_four.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"
# The unconverted script of the same name still sits in the miniprot directory.
# Only the converted one imports branch_lib, so refuse anything that does not.
grep -q '^from branch_lib import' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} is the old, unconverted script (no branch_lib import)." >&2
  echo "  Use Ascarid-Molecular-Evolution/03_rates/branch_dnds_four.py." >&2
  exit 1; }


export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"


if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
fi
PY="$(command -v python || command -v python3)"

# Probe the IMPORTS, not for sibling files. branch_lib, region_lib and
# paper_colors live in lib/, which the PYTHONPATH above already carries -- and
# branch_dnds_four.py belongs in 03_rates/ while branch_dnds_decomposition.py
# belongs in 05_analyses/, so "beside the script" is not a place a shared module
# can be. An import probe is also the stronger check: a file that exists but
# fails to import is exactly the failure a file test misses.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("branch_lib", "region_lib", "paper_colors")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("modules: branch_lib, region_lib, paper_colors all importable")
PYEOF

# 2026-10-10: codeml must get the UNROOTED tree. The branch_lib.py from before
# that date writes it rooted, which codeml fits with two unidentifiable root
# branches (the Toxocara branch split in two). Refuse to start on it.
"${PY}" - <<'PYEOF' || { echo "ERROR: lib/branch_lib.py is the pre-2026-10-10 version (rooted tree). Replace it and resubmit." >&2; exit 1; }
import branch_lib
branch_lib.check_unrooted(branch_lib.TREE_TOPOLOGY)
print("tree:", branch_lib.TREE_TOPOLOGY, "(unrooted)")
PYEOF

cd "${PROJECT_DIR}"

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/branch_dnds_four_${STAMP}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
OMEGA_CAP="${OMEGA_CAP:-5}"
YSCALE="${YSCALE:-log}"
SEED="${SEED:-0}"

# Every path gets a `:-` default. Under `set -u` a bare reference to a config
# variable the user has not defined aborts the whole script instead of letting
# `need` report it -- that is how a missing PU_ELIMINATED_BED once killed two
# species' runs silently.
RBH_AS_PU="${RBH_AS_PU:-rbh_results/orthologs.tsv}"
RBH_AS_BSC="${RBH_AS_BSC:-rbh_As_Bsc_v2/orthologs.tsv}"
RBH_PU_BSC="${RBH_PU_BSC:-rbh_Pu_Bsc_v2/orthologs.tsv}"
RBH_TOX_AS="${RBH_TOX_AS:-rbh_Tox_As/orthologs.tsv}"
RBH_TOX_PU="${RBH_TOX_PU:-rbh_Tox_Pu/orthologs.tsv}"
RBH_TOX_BSC="${RBH_TOX_BSC:-rbh_Tox_Bsc_v2/orthologs.tsv}"
AS_CDS="${AS_CDS:-As_cds.fa}";   AS_PROT="${AS_PROT:-As_prot.fa}"
PU_CDS="${PU_CDS:-Pu_cds.fa}";   PU_PROT="${PU_PROT:-Pu_prot.fa}"
BSC_CDS="${BSC_CDS:-Bsc_cds.fa}"; BSC_PROT="${BSC_PROT:-Bsc_prot.fa}"
TOX_CDS="${TOX_CDS:-Tox_cds.fa}"; TOX_PROT="${TOX_PROT:-Tox_prot.fa}"
AS_GFF="${AS_MINIPROT_GFF:-As_miniprot.gff}"
AS_RET="${AS_RETAINED_BED:-${DB_DIR}/ascaris/v50/AG.v50.retained_loose.bed}"
AS_ELIM="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
TOX_ELIM="${TOX_ELIMINATED_BED:-${DB_DIR}/toxocara/final_release/eliminated.transcripts.longest.bed}"

missing=0
need () { if [[ -f "$2" ]]; then printf '  %-14s %s\n' "$1" "$2";
           else printf '  %-14s MISSING  %s\n' "$1" "$2"; missing=1; fi; }
echo "inputs:"
need rbh_as_pu   "${RBH_AS_PU}"
need rbh_as_bsc  "${RBH_AS_BSC}"
need rbh_pu_bsc  "${RBH_PU_BSC}"
need rbh_tox_as  "${RBH_TOX_AS}"
need rbh_tox_pu  "${RBH_TOX_PU}"
need rbh_tox_bsc "${RBH_TOX_BSC}"
for f in "${AS_CDS}" "${AS_PROT}" "${PU_CDS}" "${PU_PROT}" \
         "${BSC_CDS}" "${BSC_PROT}" "${TOX_CDS}" "${TOX_PROT}"; do
  need seq "$f"
done
need gff1        "${AS_GFF}"
need retained1   "${AS_RET}"
need eliminated1 "${AS_ELIM}"
[[ "${missing}" -eq 0 ]] || { echo "ERROR: inputs missing, nothing run." >&2; exit 1; }

TOX_ARG=()
if [[ -f "${TOX_ELIM}" ]]; then
  TOX_ARG=(--tox-eliminated "${TOX_ELIM}")
  echo "  tox_elim       ${TOX_ELIM}  (S14 panels enabled)"
else
  echo "  tox_elim       MISSING -- S14 panels will be skipped"
fi

command -v codeml >/dev/null || { echo "ERROR: codeml not in PATH." >&2; exit 1; }
command -v mafft  >/dev/null || { echo "ERROR: mafft not in PATH." >&2; exit 1; }

echo "threads: ${THREADS}   outdir: ${OUTDIR}"
"${PY}" "${SCRIPT}" \
  --rbh-as-pu "${RBH_AS_PU}" --rbh-as-bsc "${RBH_AS_BSC}" \
  --rbh-pu-bsc "${RBH_PU_BSC}" --rbh-tox-as "${RBH_TOX_AS}" \
  --rbh-tox-pu "${RBH_TOX_PU}" --rbh-tox-bsc "${RBH_TOX_BSC}" \
  --cds-as "${AS_CDS}" --prot-as "${AS_PROT}" \
  --cds-pu "${PU_CDS}" --prot-pu "${PU_PROT}" \
  --cds-bsc "${BSC_CDS}" --prot-bsc "${BSC_PROT}" \
  --cds-tox "${TOX_CDS}" --prot-tox "${TOX_PROT}" \
  --retained1 "${AS_RET}" --eliminated1 "${AS_ELIM}" --gff1 "${AS_GFF}" \
  "${TOX_ARG[@]}" \
  --omega-cap "${OMEGA_CAP}" --yscale "${YSCALE}" --seed "${SEED}" \
  --threads "${THREADS}" --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
# Print the summary from the directory this run actually wrote. The previous
# version cat'd branch_dnds_four_results/ while --outdir was
# branch_dnds_four_results_v2/, so the console output described a different run.
cat "${OUTDIR}/branch_dnds_four_summary.txt"
echo ""
echo "next: sbatch run_branch_dnds_decomposition.sh ${OUTDIR}/branch_dnds_four_annotated.tsv"
