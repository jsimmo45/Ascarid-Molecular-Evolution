#!/bin/bash
#SBATCH -J cbr_profile
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=cbr_profile_%j.out
#SBATCH --error=cbr_profile_%j.err
#
# Figure 3A -- dN/dS, dN and dS against distance to the nearest CBR.
# Also draws the single-chromosome track, Figure 3C: chromosome 2 by default
# since 2026-10-30 (TRACK_CHROMOSOME; was chr01). Both of its terminal
# eliminated tracts are assembled, it carries the internal break that makes
# 2a and 2b, and it is the chromosome of Fig. 1A.
#
#   sbatch run_cbr_meta_profile.sh
#   bash   run_cbr_meta_profile.sh
#
# Replaces run_kaks_landscape.sh.
#
# Reads signed_cbr_dist1 from the annotated table rather than recomputing it,
# and cross-checks the table against a fresh computation before drawing. A table
# from the old annotate_kaks_regions.py has no signed_cbr_dist1 at all and will
# stop the script with a message saying so.
#
# --terminal-tol MUST match the value annotate_regions.py was run with, or the
# cross-check fails -- which is the point of it.
#
# Single-threaded numpy/pandas: 1 node / 1 task / 1 CPU. The bootstrap is the
# only real work and it is vectorised.

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
SCRIPT="$(find_script cbr_meta_profile.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: cbr_meta_profile.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

# The 2026-10-30 version sets the 3A fold labels as plain Arial text
# ("**** 1.69x"); an older one writes them in mathtext, which opens in
# Illustrator in DejaVu Sans. Refuse it rather than redraw the old labels.
grep -q '^def fold_label' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} predates 2026-10-30 (no fold_label); copy in the" >&2
  echo "  new cbr_meta_profile.py first." >&2; exit 1; }

# region_lib.py must sit beside it -- cbr_meta_profile.py and gene_density_cbr.py
# both import it, and that shared import is what keeps the CBR definition
# identical across the panels of Figure 3.
[[ -f "$(dirname "${SCRIPT}")/region_lib.py" ]] || {
  echo "ERROR: region_lib.py is not beside ${SCRIPT}." >&2; exit 1; }

export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"


if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
fi
PY="$(command -v python || command -v python3)"

cd "${PROJECT_DIR}"

TABLE="${ANNOTATED_TABLE:-}"
if [[ -z "${TABLE}" ]]; then
  echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2
  echo "  Point it at the annotate_regions.py run, e.g." >&2
  echo "    ANNOTATED_TABLE=\"results/annotate_regions_YYYYMMDD/kaks_annotated.tsv\"" >&2
  exit 1
fi
ELIM="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
SIZES="${AS_CHROM_SIZES:-}"
SPECIES="${SPECIES1:-Ascaris}"
TERMINAL_TOL="${TERMINAL_TOL:-1000}"
WINDOW_KB="${WINDOW_KB:-50}"
MIN_GENES="${MIN_GENES:-10}"
# The first chromosome is the panel (chromosome_track.*); any others are drawn
# the same way as chromosome_track_<name>.* for comparison.
CHROMOSOME="${TRACK_CHROMOSOME:-chr02,chr01}"

# Time in the stamp, so a second run on one day cannot overwrite the first.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/cbr_meta_profile_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-10s %s\n' "Table:" "${TABLE}" "Elim:" "${ELIM}" \
  "Sizes:" "${SIZES:-<none>}" "Species:" "${SPECIES}" \
  "Track:" "${CHROMOSOME:-<most genes>}" "Bin:" "${WINDOW_KB} kb" "Min/bin:" "${MIN_GENES} genes" \
  "TermTol:" "${TERMINAL_TOL} bp" "Output:" "${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${TABLE}" "${ELIM}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done

# The old annotate_kaks_regions.py wrote no signed distance. Catch that here with
# a clear message rather than letting the script report a missing column.
if ! head -1 "${TABLE}" | tr '\t' '\n' | grep -qx 'signed_cbr_dist1'; then
  echo "ERROR: ${TABLE} has no signed_cbr_dist1 column." >&2
  echo "  It was written by the old annotate_kaks_regions.py. Run" >&2
  echo "  run_annotate_regions.sh and repoint ANNOTATED_TABLE at its output." >&2
  exit 1
fi
for col in Ka_Ks Ka Ks chrom1 midpoint1; do
  head -1 "${TABLE}" | tr '\t' '\n' | grep -qx "${col}" || {
    echo "ERROR: ${TABLE} has no ${col} column." >&2; exit 1; }
done
printf '  %-44s %8d rows\n' "$(basename "${TABLE}")" "$(( $(wc -l < "${TABLE}") - 1 ))"

SIZE_ARG=()
if [[ -n "${SIZES}" ]]; then
  [[ -s "${SIZES}" ]] || { echo "ERROR: missing or empty: ${SIZES}" >&2; exit 1; }
  SIZE_ARG=(--chrom-sizes "${SIZES}")
else
  echo "  NOTE: no AS_CHROM_SIZES. The cross-check still runs, but without"
  echo "        sequence lengths a telomeric boundary is treated as a CBR in the"
  echo "        recomputation and the check may report a false disagreement."
fi
CHROM_ARG=()
[[ -n "${CHROMOSOME}" ]] && CHROM_ARG=(--chromosome "${CHROMOSOME}")

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
# --max-dist-kb 500 (was 1000 until 2026-10-26): both arms run to 500 kb, as the
# gene-density panel does. The retained arm used to reach 1000 kb while the
# eliminated arm stops near 500 kb, where bins fall below --min-genes.
"${PY}" "${SCRIPT}" \
  --annotated "${TABLE}" \
  --eliminated-bed "${ELIM}" \
  "${SIZE_ARG[@]}" \
  "${CHROM_ARG[@]}" \
  --terminal-tol "${TERMINAL_TOL}" \
  --window-kb "${WINDOW_KB}" \
  --max-dist-kb 500 \
  --min-genes "${MIN_GENES}" \
  --step-window-kb 150 \
  --n-boot 2000 \
  --species "${SPECIES}" \
  --width single \
  --height-mm 130 \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/cbr_meta_profile_summary.txt"
