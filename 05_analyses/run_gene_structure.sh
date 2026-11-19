#!/bin/bash
#SBATCH -J gene_struct
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=gene_structure_%j.out
#SBATCH --error=gene_structure_%j.err
#
# Figure 2D -- genomic span, intron fraction and CDS length by region, plus
# per-gene intergenic distance written to the output table.
#
#   sbatch run_gene_structure.sh
#   bash   run_gene_structure.sh
#
# Ascaris only. Parascaris models are miniprot projections and Toxocara models
# are TransDecoder ORFs on assembled transcripts, so intron structure and
# genomic span are not measured comparably in either -- running this for them
# would produce numbers that look comparable and are not.
#
# Single-threaded pandas/scipy, hence 1 node / 1 task / 1 CPU.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

find_config() {
  local d="${HERE}"
  for _ in 1 2 3 4; do
    [[ -f "${d}/config.sh" ]] && { printf '%s' "${d}/config.sh"; return 0; }
    d="$(dirname "${d}")"
  done
  return 1
}
CONFIG="${KAKS_CONFIG:-$(find_config || true)}"
if [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then
  echo "ERROR: no config.sh found in ${HERE} or up to four levels above it." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

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

CDS="${AS_CDS:-As_cds.fa}"
GFF="${AS_MINIPROT_GFF:-As_miniprot.gff}"
TABLE="${ANNOTATED_TABLE:-kaks_regions/kaks_annotated.tsv}"

STAMP="$(date +%Y%m%d)"
OUTDIR="${RESULTS_DIR}/gene_structure_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

echo "CDS:     ${CDS}"
echo "Models:  ${GFF}"
echo "Table:   ${TABLE}"
echo "Output:  ${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${CDS}" "${GFF}" "${TABLE}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: ${PROJECT_DIR}/$f" >&2; exit 1; }
done
printf '  %-44s %8d rows\n' "${TABLE}" "$(( $(wc -l < "${TABLE}") - 1 ))"
printf '  %-44s %8d mRNA features\n' "${GFF}" "$(awk -F'\t' '$3=="mRNA"' "${GFF}" | wc -l)"

# Gene-ID overlap between the table and the GFF. Intergenic distance is computed
# from every model in the GFF, but the per-gene join still needs the analysed
# genes to be present there, and a silent mismatch would leave every span NaN.
HEADER="$(head -1 "${TABLE}")"
GENE_IDX="$(echo "${HEADER}" | tr '\t' '\n' | grep -n -x -F "${GENE_COL:-gene1}" | cut -d: -f1 | head -1)"
if [[ -z "${GENE_IDX}" ]]; then
  echo "ERROR: column '${GENE_COL:-gene1}' not found in ${TABLE}" >&2
  echo "${HEADER}" | tr '\t' '\n' | sed 's/^/    /' >&2
  exit 1
fi
N_MATCH="$(awk -v col="${GENE_IDX}" -F'\t' '
  NR==FNR { if ($3=="mRNA") { split($9, a, /[;=]/); for (i=1;i<length(a);i++) if (a[i]=="ID") ids[a[i+1]]=1 } next }
  FNR==1 { next }
  { if ($col in ids) n++ }
  END { print n+0 }' "${GFF}" "${TABLE}")"
N_ROWS=$(( $(wc -l < "${TABLE}") - 1 ))
echo "  Gene-ID overlap: ${N_MATCH} / ${N_ROWS} table rows found in ${GFF}"
if [[ "${N_MATCH}" -lt $(( N_ROWS / 2 )) ]]; then
  echo "ERROR: fewer than half the gene IDs appear in the GFF." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
"${PY}" "${SCRIPT_DIR:-${HERE}}/gene_structure.py" \
  --cds-fasta "${CDS}" \
  --gff "${GFF}" \
  --rates-table "${TABLE}" \
  --gene-col "${GENE_COL:-gene1}" \
  --region-col "${REGION_COL:-region1}" \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/gene_structure_summary.txt"
