#!/bin/bash
#SBATCH -J codon_usage
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=codon_usage_%j.out
#SBATCH --error=codon_usage_%j.err
#
# Codon usage bias (ENC, GC3s, dENC, CAI) for retained versus eliminated genes.
# Figure 2B-C for Ascaris; supplementary S2 for the other species.
#
#   sbatch run_codon_usage_bias.sh ascaris
#   sbatch run_codon_usage_bias.sh parascaris
#   sbatch run_codon_usage_bias.sh toxocara
#   bash   run_codon_usage_bias.sh ascaris --check    # preflight only
#
# Replaces run_codon_usage_bias.sh (Ascaris hardcoded), run_codon_usage_tox.sh
# (one-off) and run_codon_usage_bias_multi.sh. The Python was always
# species-agnostic; only the old wrapper pinned it to Ascaris, which is the sole
# reason Figure 2B-C is Ascaris-only.
#
# Baylisascaris is deliberately absent. It has CDS and a native Augustus
# annotation but no eliminated-regions BED, so its genes can only be labelled by
# projection from Ascaris orthologs. Codon usage on projected labels is a
# materially weaker claim and should not be run as though it were equivalent.
#
# Single-threaded pandas/numpy throughout, hence 1 node / 1 task / 1 CPU.
#
# 2026-10-12: config.sh and the .py file are found from SLURM_SUBMIT_DIR and
# PWD as well as from the script's own location, which under sbatch is
# SLURM's spool copy (see below). Run it from the repo or anywhere under it.
# Everything machine-specific lives in config.sh; see config.example.sh.
#
# PANEL_B in config.sh chooses the boxplot metric in the Figure 2 strip. It
# defaults to denc: dENC is the only metric here that separates translational
# selection from GC3s composition, and the two groups differ in GC3s. Set
# PANEL_B="cai" to go back to the original panel.

set -euo pipefail

# Under sbatch, SLURM copies the batch script to a spool directory on the compute
# node, so BASH_SOURCE points at /var/spool/slurm/... and neither config.sh nor
# the .py file is anywhere near it (jobs 6383411-6383414 and 6383612 failed on
# exactly that, 2026-10-12). SLURM_SUBMIT_DIR is the directory sbatch was run
# from, which is on the shared filesystem, so search from there as well, as
# run_genome_trace.sh does. $KAKS_CONFIG still overrides.
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
  echo "  Copy config.example.sh to config.sh and edit it, or set" >&2
  echo "  KAKS_CONFIG=/full/path/to/config.sh and resubmit." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

# The .py file, by the same reasoning: never assume it sits beside BASH_SOURCE.
REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" \
           "${REPO_ROOT}/06_figures" "${REPO_ROOT}" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script codon_usage_bias.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: codon_usage_bias.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

SPECIES="${1:-}"
CHECK_ONLY="${2:-}"
if [[ -z "${SPECIES}" ]]; then
  echo "ERROR: give a species: ascaris | parascaris | toxocara" >&2
  exit 1
fi

# --- environment ------------------------------------------------------------
if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
  echo "python: $(command -v python)"
else
  echo "No conda config; using python on PATH: $(command -v python || command -v python3)"
fi
PY="$(command -v python || command -v python3)"

cd "${PROJECT_DIR}"

##############################################################################
# Per-species inputs. Each species needs a CDS FASTA, a table carrying its gene
# IDs, and a column in that table holding ITS OWN retained/eliminated calls.
##############################################################################
REF_ARG=()
case "${SPECIES}" in

  ascaris)
    # Figure 2B-C. region1 = Ascaris' own calls.
    CDS="${AS_CDS:-As_cds.fa}"
    TABLE="${ANNOTATED_TABLE:-kaks_regions/kaks_annotated.tsv}"
    GENE_COL=gene1
    REGION_COL=region1
    # Expression-based CAI reference: see make_cai_reference.py. Without it CAI
    # falls back to the most codon-biased retained genes, which makes the
    # retained-vs-eliminated CAI contrast partly circular.
    [[ -s "${CAI_REFERENCE:-cai_reference_ascaris.txt}" ]] && \
      REF_ARG=(--ref-genes "${CAI_REFERENCE:-cai_reference_ascaris.txt}")
    ;;

  parascaris)
    # region2 = Parascaris' own calls, from pu_v1.13_eliminated.bed (1,259
    # eliminated genes), not projected from Ascaris -- an independent
    # replication with better gene models than Toxocara and no saturation
    # caveat. Already present in the annotated table; no new annotation needed.
    CDS="${PU_CDS:-Pu_cds.fa}"
    TABLE="${ANNOTATED_TABLE:-kaks_regions/kaks_annotated.tsv}"
    GENE_COL=gene2
    REGION_COL=region2
    ;;

  toxocara)
    # Column names confirmed 2026-08-12 against the actual header. The Toxocara
    # annotation script writes 'tox_region', not 'region1'.
    #
    # gene1 is Toxocara here (species1 in the Tox x Ascaris pairing) and holds
    # TransDecoder ORF IDs like MSTRG.17238.1.p1, which is what Tox_cds.fa is
    # keyed on. Do NOT use tox_transcript -- that is the StringTie transcript
    # one level up, and it will not match the CDS FASTA.
    #
    # partner_region holds the Ascaris calls for the same pairs; swapping to it
    # gives the reciprocal-labelling control rather than Toxocara's own
    # annotation.
    #
    # No expression-based reference exists: the Toxocara TPM matrices cover
    # eliminated genes only. CAI therefore uses the fallback here and is
    # reported as supporting evidence. Lead with dENC, which needs no reference.
    CDS="${TOX_CDS:-Tox_cds.fa}"
    TABLE="${TOX_ANNOTATED_TABLE:-kaks_Tox_As_annotated/kaks_annotated.tsv}"
    GENE_COL=gene1
    REGION_COL=tox_region
    ;;

  *)
    echo "ERROR: unknown species '${SPECIES}'" >&2
    echo "Choose: ascaris | parascaris | toxocara" >&2
    exit 1
    ;;
esac

STAMP="$(date +%Y%m%d)"
OUTDIR="${RESULTS_DIR}/codon_usage_${SPECIES}_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

echo "Species:     ${SPECIES}"
echo "CDS:         ${CDS}"
echo "Table:       ${TABLE}"
echo "Gene column: ${GENE_COL}"
echo "Region col:  ${REGION_COL}"
echo "CAI ref:     ${REF_ARG[1]:-<fallback, partly circular>}"
echo "Panel 2B:    ${PANEL_B:-denc}"
echo "Output:      ${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
# The Python fails clearly on a bad column, but only after parsing a whole
# genome's CDS. Checking here is instant. Gene-ID mismatch between FASTA headers
# and the table column is the likeliest failure and is otherwise silent until
# the very end.
for f in "${CDS}" "${TABLE}"; do
  [[ -s "${f}" ]] || { echo "ERROR: missing or empty: ${PROJECT_DIR}/${f}" >&2; exit 1; }
done

HEADER="$(head -1 "${TABLE}")"
col_index () { echo "${HEADER}" | tr '\t' '\n' | grep -n -x -F "$1" | cut -d: -f1 | head -1; }
GENE_IDX="$(col_index "${GENE_COL}" || true)"
REGION_IDX="$(col_index "${REGION_COL}" || true)"

if [[ -z "${GENE_IDX}" || -z "${REGION_IDX}" ]]; then
  echo "ERROR: column not found in ${TABLE}" >&2
  echo "  gene column   '${GENE_COL}'   -> ${GENE_IDX:-NOT FOUND}" >&2
  echo "  region column '${REGION_COL}' -> ${REGION_IDX:-NOT FOUND}" >&2
  echo "  columns present:" >&2
  echo "${HEADER}" | tr '\t' '\n' | sed 's/^/    /' >&2
  exit 1
fi

echo "Region labels in ${REGION_COL}:"
awk -v c="${REGION_IDX}" -F'\t' 'NR>1 {n[$c]++} END {for (k in n) printf "  %-14s %6d\n", k, n[k]}' \
  "${TABLE}" | sort -k2 -rn

N_MATCH="$(awk -v col="${GENE_IDX}" '
  NR==FNR { if (substr($0,1,1) == ">") { h=substr($0,2); split(h,a,/[ \t]/); ids[a[1]]=1 } next }
  FNR == 1 { next }
  { split($0, f, "\t"); if (f[col] in ids) n++ }
  END { print n+0 }' "${CDS}" "${TABLE}")"
N_ROWS=$(( $(wc -l < "${TABLE}") - 1 ))
echo ""
echo "Gene-ID overlap: ${N_MATCH} / ${N_ROWS} table rows found in ${CDS}"

if [[ "${N_MATCH}" -lt $(( N_ROWS / 2 )) ]]; then
  echo "" >&2
  echo "ERROR: fewer than half the IDs in '${GENE_COL}' appear in ${CDS}." >&2
  echo "Wrong FASTA for this species, or a different ID convention." >&2
  echo "  first three table IDs:" >&2
  awk -v c="${GENE_IDX}" -F'\t' 'NR>1 && NR<=4 {print "    " $c}' "${TABLE}" >&2
  echo "  first three FASTA headers:" >&2
  grep '^>' "${CDS}" | head -3 | sed 's/^/    /' >&2
  exit 1
fi

if [[ "${CHECK_ONLY}" == "--check" ]]; then
  echo ""
  echo "Preflight passed. Rerun without --check to analyse."
  exit 0
fi

# --- run --------------------------------------------------------------------
mkdir -p "${OUTDIR}"
"${PY}" "${SCRIPT}" \
  --cds-fasta "${CDS}" \
  --rates-table "${TABLE}" \
  --gene-col "${GENE_COL}" \
  --region-col "${REGION_COL}" \
  --min-codons 100 \
  --panel-b "${PANEL_B:-denc}" \
  "${REF_ARG[@]}" \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/codon_usage_summary.txt"
