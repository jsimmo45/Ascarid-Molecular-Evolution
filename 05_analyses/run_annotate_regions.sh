#!/bin/bash
#SBATCH -J annot_regions
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=annot_regions_%j.out
#SBATCH --error=annot_regions_%j.err
#
# Upstream: region calls, tract types and CBR distances for every ortholog pair.
# Writes kaks_annotated.tsv, plus supplementary fig. S7 (S5 before the
# 2026-10-12 renumbering; the files keep their S5A_/S5B_ stems).
#
#   sbatch run_annotate_regions.sh
#   bash   run_annotate_regions.sh
#
# Replaces run_annotate_kaks.sh.
#
# Regenerates kaks_annotated.tsv, into a new timestamped folder. Every
# downstream script reads that table through ANNOTATED_TABLE in config.sh. The
# 2026-08-17 conversion changed it in three ways that alter results (sign
# convention on CBR distance, the CBR set, locus choice for multiply aligned
# genes); every figure was rerun against that table.
#
# 2026-10-09: the stats table gains the pooled CBR-proximal vs interior test
# (family proximity_pooled) that Results quotes. The annotation code did not
# change, so on the same inputs the table should come out byte-identical to the
# one in use. The DONE block checks that against ANNOTATED_TABLE, and checks that
# every earlier test is reproduced, and says whether anything downstream needs
# rerunning. Do not repoint ANNOTATED_TABLE unless it reports a difference you
# have explained.
#
# 2026-10-12: figure changes only (S7 review). Brackets sit above the drawn
# data, points above the cap are counted in the summary, S7A is not hatched,
# S7B has no title. The table and every test are unchanged, so the DONE block
# should again report kaks_annotated.tsv byte-identical to ANNOTATED_TABLE.
#
# Single-threaded numpy/pandas: 1 node / 1 task / 1 CPU. 16G because both GFFs
# and the full rate table are held at once.

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
SCRIPT="$(find_script annotate_regions.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: annotate_regions.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

# region_lib.py must sit beside it -- annotate_regions.py and gene_density_cbr.py
# both import it, and that shared import is what keeps the CBR definition
# identical across the panels of Figure 4.
[[ -f "$(dirname "${SCRIPT}")/region_lib.py" ]] || {
  echo "ERROR: region_lib.py is not beside ${SCRIPT}." >&2; exit 1; }

PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHONPATH

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
fi
PY="$(command -v python || command -v python3)"

cd "${PROJECT_DIR}"

KAKS="${KAKS_AS_PU:-kaks_results/kaks_filtered.tsv}"
GFF1="${AS_MINIPROT_GFF:-As_miniprot.gff}"
GFF2="${PU_MINIPROT_GFF:-Pu_miniprot.gff}"
RET1="${AS_RETAINED_BED:-${DB_DIR}/ascaris/v50/AG.v50.retained_loose.bed}"
ELIM1="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
RET2="${PU_RETAINED_BED:-${DB_DIR}/parascaris/v1.13/somatic_genome/germ_to_soma_mapping.bed}"
ELIM2="${PU_ELIMINATED_BED:-${DB_DIR}/parascaris/v1.13/pu_v1.13_eliminated.bed}"
SIZES1="${AS_CHROM_SIZES:-}"
SIZES2="${PU_CHROM_SIZES:-}"

# Must match gene_density_cbr.py, or Figure 4's panels use different CBR sets.
CBR_WINDOW="${CBR_WINDOW:-100000}"
TERMINAL_TOL="${TERMINAL_TOL:-1000}"

# Seconds in the stamp: a date-only stamp sends two runs on the same day into
# one folder, the second overwriting the first.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/annotate_regions_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

printf '%-10s %s\n' "Rates:" "${KAKS}" "GFF1:" "${GFF1}" "GFF2:" "${GFF2}" \
  "Ret1:" "${RET1}" "Elim1:" "${ELIM1}" "Ret2:" "${RET2}" "Elim2:" "${ELIM2}" \
  "Sizes1:" "${SIZES1:-<none>}" "Sizes2:" "${SIZES2:-<none>}" \
  "CBR win:" "${CBR_WINDOW} bp" "TermTol:" "${TERMINAL_TOL} bp" \
  "Output:" "${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${KAKS}" "${GFF1}" "${GFF2}" "${RET1}" "${ELIM1}" "${RET2}" "${ELIM2}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done

# Sequence names must agree between each GFF and its BEDs. They come from
# different pipelines, and a mismatch is not an error in either file: every gene
# would simply come back unannotated. The script fails on a total mismatch, but
# a *partial* one is silent, so report the fraction here.
check_names() {                     # check_names <gff> <bed> <label>
  local shared bed_seq
  shared="$(awk -F'\t' '
    NR==FNR { if ($1 !~ /^#/ && $3=="mRNA") g[$1]=1; next }
    $1 !~ /^#/ { if ($1 in g) s[$1]=1 }
    END { print length(s) }' "$1" "$2")"
  bed_seq="$(awk -F'\t' '$1 !~ /^#/ {s[$1]=1} END {print length(s)}' "$2")"
  printf '  %-42s %3d / %3d BED sequences match the GFF\n' "$3" "${shared}" "${bed_seq}"
  if [[ "${shared}" -eq 0 ]]; then
    echo "ERROR: no sequence name shared between $1 and $2" >&2
    echo "  BED: $(awk -F'\t' '$1 !~ /^#/ {print $1}' "$2" | head -3 | tr '\n' ' ')" >&2
    echo "  GFF: $(awk -F'\t' '$3=="mRNA" {print $1}' "$1" | head -3 | tr '\n' ' ')" >&2
    exit 1
  fi
}
check_names "${GFF1}" "${RET1}"  "retained ${SPECIES1:-Ascaris}"
check_names "${GFF1}" "${ELIM1}" "eliminated ${SPECIES1:-Ascaris}"
check_names "${GFF2}" "${RET2}"  "retained ${SPECIES2:-Parascaris}"
check_names "${GFF2}" "${ELIM2}" "eliminated ${SPECIES2:-Parascaris}"

# Gene-identifier namespace: the rate table's gene columns must exist in the GFFs.
for pair in "1:${GFF1}" "2:${GFF2}"; do
  N="${pair%%:*}"; G="${pair#*:}"
  COL="$(head -1 "${KAKS}" | tr '\t' '\n' | grep -n -x -F "gene${N}" | cut -d: -f1 | head -1)"
  if [[ -z "${COL}" ]]; then
    echo "ERROR: column 'gene${N}' not found in ${KAKS}" >&2
    head -1 "${KAKS}" | tr '\t' '\n' | sed 's/^/    /' >&2
    exit 1
  fi
  MATCH="$(awk -v col="${COL}" -F'\t' '
    NR==FNR { if ($3=="mRNA") { split($9, a, /[;=]/); for (i=1;i<length(a);i++) if (a[i]=="ID") ids[a[i+1]]=1 } next }
    FNR==1 { next }
    { if ($col in ids) n++ }
    END { print n+0 }' "${G}" "${KAKS}")"
  ROWS=$(( $(wc -l < "${KAKS}") - 1 ))
  printf '  %-42s %6d / %6d rate rows found in the GFF\n' "gene${N} identifiers" "${MATCH}" "${ROWS}"
  [[ "${MATCH}" -gt $(( ROWS / 2 )) ]] || {
    echo "ERROR: fewer than half the gene${N} identifiers appear in ${G}." >&2; exit 1; }
done

SIZE_ARGS=()
[[ -n "${SIZES1}" ]] && SIZE_ARGS+=(--chrom-sizes1 "${SIZES1}")
[[ -n "${SIZES2}" ]] && SIZE_ARGS+=(--chrom-sizes2 "${SIZES2}")
if [[ -z "${SIZES1}" || -z "${SIZES2}" ]]; then
  echo "  NOTE: a .fai is missing for at least one species. Without it a terminal"
  echo "        eliminated tract's telomeric edge is counted as a CBR, which"
  echo "        changes which genes are called CBR-proximal."
fi

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT}" \
  --kaks "${KAKS}" \
  --gff1 "${GFF1}" --gff2 "${GFF2}" \
  --retained1 "${RET1}" --eliminated1 "${ELIM1}" \
  --retained2 "${RET2}" --eliminated2 "${ELIM2}" \
  "${SIZE_ARGS[@]}" \
  --cbr-window "${CBR_WINDOW}" \
  --terminal-tol "${TERMINAL_TOL}" \
  --species1 "${SPECIES1:-Ascaris}" --species2 "${SPECIES2:-Parascaris}" \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/annotate_regions_summary.txt"
echo ""

NEW_TABLE="${OUTDIR}/kaks_annotated.tsv"
NEW_STATS="${OUTDIR}/annotate_regions_stats.tsv"

echo "===== POOLED CBR-PROXIMAL vs INTERIOR (all regions) ====="
awk -F'\t' 'NR == 1 || $1 == "proximity_pooled"' "${NEW_STATS}"
echo ""

# Is this run a drop-in replacement for the table downstream scripts read?
echo "===== COMPARISON WITH THE TABLE IN USE ====="
CUR_TABLE="${ANNOTATED_TABLE:-}"
if [[ -n "${CUR_TABLE}" && -s "${CUR_TABLE}" ]]; then
  echo "  ANNOTATED_TABLE: ${CUR_TABLE}"
  if cmp -s "${NEW_TABLE}" "${CUR_TABLE}"; then
    echo "  kaks_annotated.tsv: byte-identical. Downstream results stay valid;"
    echo "  leave ANNOTATED_TABLE where it is."
  else
    echo "  kaks_annotated.tsv: DIFFERS from the table in use."
    CMP_OUT="$(cmp "${NEW_TABLE}" "${CUR_TABLE}" 2>&1 || true)"
    echo "    ${CMP_OUT}"
    # Show the first differing line from both files, so a difference in number
    # formatting alone can be told apart from a difference in values.
    L="$(printf '%s\n' "${CMP_OUT}" | sed -n 's/.*line \([0-9][0-9]*\).*/\1/p')"
    if [[ -n "${L}" ]]; then
      echo "    line ${L}, this run then the table in use (first 300 characters):"
      sed -n "${L}p" "${NEW_TABLE}" | cut -c1-300 | sed 's/^/      new:    /'
      sed -n "${L}p" "${CUR_TABLE}" | cut -c1-300 | sed 's/^/      in use: /'
    fi
    echo "  Downstream results were computed on the table in use, not this one."
    echo "  Find the cause before repointing ANNOTATED_TABLE."
  fi
  CUR_STATS="$(dirname "${CUR_TABLE}")/annotate_regions_stats.tsv"
  if [[ -s "${CUR_STATS}" ]]; then
    # The new rows are inserted by family, so drop them and the remainder
    # must match the earlier stats table line for line.
    if cmp -s <(awk -F'\t' '$1 != "proximity_pooled"' "${NEW_STATS}") "${CUR_STATS}"; then
      echo "  annotate_regions_stats.tsv: every earlier test reproduced exactly;"
      echo "  the only new rows are proximity_pooled."
    else
      echo "  annotate_regions_stats.tsv: earlier tests DIFFER from ${CUR_STATS}"
      echo "  (< this run, > the table in use):"
      diff <(awk -F'\t' '$1 != "proximity_pooled"' "${NEW_STATS}") "${CUR_STATS}" \
        | head -8 | cut -c1-200 | sed 's/^/    /' || true
    fi
  else
    echo "  (no annotate_regions_stats.tsv beside ANNOTATED_TABLE to compare with)"
  fi
else
  echo "  ANNOTATED_TABLE is unset or missing in config.sh. To use this run:"
  echo "      ANNOTATED_TABLE=\"${OUTDIR#${PROJECT_DIR}/}/kaks_annotated.tsv\""
fi
