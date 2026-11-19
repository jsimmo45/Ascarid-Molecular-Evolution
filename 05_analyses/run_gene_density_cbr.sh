#!/bin/bash
#SBATCH -J gene_density
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=gene_density_%j.out
#SBATCH --error=gene_density_%j.err
#
# Figure 3B -- gene density against signed distance to the nearest CBR.
#
#   sbatch run_gene_density_cbr.sh
#   bash   run_gene_density_cbr.sh
#
# Replaces run_replot_gene_density.sh. The old script reported "genes per kb"
# about 600x too high because it pooled genes across every CBR while dividing by
# a single bin width; this one divides by the base pairs of genome actually at
# that distance. See the docstring in gene_density_cbr.py.
#
# Single-threaded numpy/pandas, hence 1 node / 1 task / 1 CPU. The base-pair
# denominator tiles the genome at 1 kb, so peak memory is a few hundred MB on a
# 300 Mb assembly; 8G is generous.

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

# The repo may sit *below* the directory sbatch was run from -- submitting from
# the data directory with `sbatch Ascarid-Molecular-Evolution/05_analyses/...`
# is the obvious thing to do, and an upward-only search never sees it. So probe
# one level down as well, accepting only a directory that also looks like the
# repo, so a stray config.sh elsewhere cannot be picked up by accident.
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
  echo "  Set KAKS_CONFIG=/full/path/to/config.sh, or submit with" >&2
  echo "    sbatch --export=ALL,KAKS_CONFIG=/full/path/config.sh $0" >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=config.example.sh
source "${CONFIG}"

# The .py cannot be located from BASH_SOURCE under sbatch either. Anchor on the
# repo root, which is wherever config.sh was found.
REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
find_script() {                     # find_script <basename>
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" \
           "${REPO_ROOT}/06_figures" "${REPO_ROOT}" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script gene_density_cbr.py || true)"
if [[ -z "${SCRIPT}" ]]; then
  echo "ERROR: gene_density_cbr.py not found near ${REPO_ROOT}." >&2
  exit 1
fi
echo "script: ${SCRIPT}"

# The 2026-10-30 version carries the 2026-10-15 hand edits (plain y label,
# "mean" in the lower corner, top margin level with 3A, upright region labels).
# Refuse an older one rather than redraw what was fixed by hand.
grep -q '^def _label_mean_corner' "${SCRIPT}" || {
  echo "ERROR: ${SCRIPT} predates 2026-10-30 (no _label_mean_corner); copy" >&2
  echo "  in the new gene_density_cbr.py first." >&2; exit 1; }

# paper_colors.py lives in lib/ in the published layout but beside the scripts in
# a flat working tree. Put both on PYTHONPATH rather than depending on which.
export PYTHONPATH="$(dirname "${SCRIPT}"):${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

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

SPECIES="${SPECIES:-Ascaris}"
GFF="${AS_MINIPROT_GFF:-As_miniprot.gff}"
ELIM="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
SIZES="${AS_CHROM_SIZES:-}"
# Override per run:  sbatch --export=ALL,TERMINAL_TOL=0,RUN_TAG=tol0 ...
TERMINAL_TOL="${TERMINAL_TOL:-1000}"

# Time in the stamp, so a second run on one day cannot overwrite the first.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/gene_density_cbr_${STAMP}${RUN_TAG:+_${RUN_TAG}}"

echo "Species: ${SPECIES}"
echo "Models:  ${GFF}"
echo "Elim:    ${ELIM}"
echo "Sizes:   ${SIZES:-<none: lengths inferred from the last feature>}"
echo "TermTol: ${TERMINAL_TOL} bp"
echo "Output:  ${OUTDIR}"
echo ""

# --- preflight --------------------------------------------------------------
for f in "${GFF}" "${ELIM}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done
printf '  %-44s %8d mRNA features\n' "$(basename "${GFF}")" \
  "$(awk -F'\t' '$3=="mRNA"' "${GFF}" | wc -l)"
printf '  %-44s %8d intervals\n' "$(basename "${ELIM}")" \
  "$(grep -vc '^#' "${ELIM}")"

# Sequence names must agree between the GFF and the BED. They come from
# different pipelines, and a mismatch is not an error in either file -- the
# script would simply find no CBR on any sequence carrying genes, exclude every
# gene, and report an empty profile rather than a wrong one. Catch it here.
N_SHARED="$(awk -F'\t' '
  NR==FNR { if ($1 !~ /^#/ && $3=="mRNA") g[$1]=1; next }
  $1 !~ /^#/ { if ($1 in g) s[$1]=1 }
  END { print length(s) }' "${GFF}" "${ELIM}")"
N_BED_SEQ="$(awk -F'\t' '$1 !~ /^#/ {s[$1]=1} END {print length(s)}' "${ELIM}")"
echo "  Sequence names shared with the GFF: ${N_SHARED} / ${N_BED_SEQ} in the BED"
if [[ "${N_SHARED}" -eq 0 ]]; then
  echo "ERROR: no sequence name in ${ELIM} appears in ${GFF}." >&2
  echo "  BED examples: $(awk -F'\t' '$1 !~ /^#/ {print $1}' "${ELIM}" | head -3 | tr '\n' ' ')" >&2
  echo "  GFF examples: $(awk -F'\t' '$3=="mRNA" {print $1}' "${GFF}" | head -3 | tr '\n' ' ')" >&2
  exit 1
fi

SIZES_ARG=()
if [[ -n "${SIZES}" ]]; then
  [[ -s "${SIZES}" ]] || { echo "ERROR: missing or empty: ${SIZES}" >&2; exit 1; }
  SIZES_ARG=(--chrom-sizes "${SIZES}")
else
  echo "  NOTE: no AS_CHROM_SIZES in config.sh. Sequence lengths will be taken"
  echo "        from the last annotated feature, which underestimates every one"
  echo "        of them and makes the outermost retained bins read slightly"
  echo "        dense. Point AS_CHROM_SIZES at the .fai to remove this."
fi

mkdir -p "${OUTDIR}"

# --- run --------------------------------------------------------------------
# 25 kb bins is the resolution the figure outline specifies. --extra-windows
# draws 10 and 50 kb alongside it in a separate control panel, which answers the
# obvious referee question about binning without cluttering the main panel.
"${PY}" "${SCRIPT}" \
  --gff "${GFF}" \
  --eliminated-bed "${ELIM}" \
  "${SIZES_ARG[@]}" \
  --species "${SPECIES}" \
  --window-kb 25 \
  --extra-windows 10 50 \
  --max-dist-kb 500 \
  --tile-bp 1000 \
  --terminal-tol "${TERMINAL_TOL:-1000}" \
  --width single \
  --height-mm 55 \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/gene_density_summary.txt"
