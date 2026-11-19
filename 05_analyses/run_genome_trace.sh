#!/bin/bash
#SBATCH -J genome_trace
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=genome_trace_%j.out
#SBATCH --error=genome_trace_%j.err
#
# Figure 3  -- dN/dS along whole chromosomes, one run per species.
# Figure S6 -- the dN and dS runs, and the other two species (S4 before the
#              2026-10-12 renumbering; the STAMP below names the run, not the
#              figure number).
#
#   sbatch run_genome_trace.sh ascaris
#   sbatch run_genome_trace.sh all
#   sbatch run_genome_trace.sh ascaris Ka_Ks          # one metric only
#
# Replaces run_plot_kaks_genome_trace.sh. Changes worth knowing:
#   * The continuous 240 Mb Parascaris panel is gone. It was unreadable at
#     column width and only the wrapped version was ever usable.
#   * --min-genes 10, matching the CBR meta-profile bin floor. A window holding
#     one gene was being drawn as if its
#     "median" were comparable to the rest of the line.
#   * The figure is built at the journal column width with type in points. The
#     old 20-inch canvas printed its 13 pt labels at under 5 pt once placed.
#
# 2026-10-09:
#   * A window whose neighbours on both sides fall below --min-genes is drawn
#     as a dot. Before, it was not drawn at all (13 of 2,114 Baylisascaris
#     windows in S4v2). The count is in each summary's Truncation block.
#   * Baylisascaris panels are labelled 'Scaffold N' instead of the full
#     GWHBAVG accession, which overlapped in four of the five S4v2 rows. The
#     summary reports any adjacent labels that still collide.
#   * Parascaris: carets for windows above the y-cap are drawn and counted on
#     their own slice only, not once per slice of the wrapped chromosome.
#
# 2026-10-12 (figure review of S6):
#   * y labels are plain text (dN/dS, dN, dS, upright), and every grid carries
#     one "Position (Mb)" title under it.
#   * Parascaris: the single germline chromosome, chrX in the assembly, is
#     labelled 'Germline chromosome' (--label-sub '^chrX$' ...).
#   * Each summary lists the sequences the filters remove, with lengths and
#     gene counts, under "Sequences not drawn" (Baylisascaris: the scaffolds
#     below --min-chrom-mb 5). The legend n counts only genes on the sequences
#     drawn.
#
# 2026-10-12 (later): Baylisascaris draws every scaffold of 5 Mb or more (26).
#   The 24-panel cap (--max-chroms 24) had dropped scaffolds 20 (5.85 Mb) and
#   27 (6.33 Mb), which pass the size floor. Rows are now packed to BSC_ROW_MB
#   (60) like the Ascaris rows, so the 26 panels fill five rows of nearly equal
#   width at the same canvas height, instead of a sixth row holding one panel.
#
# Figure 3:  sbatch run_genome_trace.sh ascaris Ka_Ks
# S6:        STAMP=S6v4 ROW_HEIGHT=12 NO_LEGEND=1 sbatch run_genome_trace.sh all
#            (Baylisascaris alone after the cap change:
#             STAMP=S6v5 ROW_HEIGHT=12 NO_LEGEND=1 sbatch run_genome_trace.sh baylisascaris)
#
# METRICS defaults to all three, which is what S6 needs alongside Figure 3.
# Set YMAX to put two species on one scale -- by default each run picks its own
# percentile and separately-drawn species are NOT comparable.

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
SCRIPT="$(find_script genome_trace.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: genome_trace.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

# region_lib.py must sit beside it -- genome_trace.py and gene_density_cbr.py
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

WHICH="${1:-all}"
METRICS="${2:-Ka_Ks Ka Ks}"
# Seconds matter. A date-only stamp meant every invocation on the same day wrote
# into the SAME results directory, which defeats the whole point of make_outdir()
# -- and if two invocations overlapped, both wrote provenance.json at once and
# left it unparseable. Override STAMP to deliberately add to an existing run.
STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"

TABLE="${ANNOTATED_TABLE:-}"
[[ -n "${TABLE}" ]] || { echo "ERROR: ANNOTATED_TABLE is not set in config.sh." >&2; exit 1; }

# Per-species tuning, overridable from the environment.
ASC_ROW_MB="${ASC_ROW_MB:-75}"     # row packing target; the summary prints the
                                   # actual Mb per row, so tune from that
PU_WRAP_MB="${PU_WRAP_MB:-60}"     # Parascaris wrap length
BSC_ROW_MB="${BSC_ROW_MB:-60}"     # Baylisascaris row packing target (26
                                   # scaffolds >= 5 Mb -> five rows of 52-59 Mb)
YMAX="${YMAX:-}"                   # set to compare species on one scale
ROW_HEIGHT="${ROW_HEIGHT:-26}"     # mm per chromosome row; drop to ~14-16 for a
                                   # supplementary panel that has to share a page
MIN_GENES="${MIN_GENES:-10}"   # windows below this are dropped, not drawn.
                                   # 10 matches the CBR meta-profile bin floor
NO_LEGEND="${NO_LEGEND:-0}"        # 1 puts the legend in the caption instead,
                                   # buying back 7 mm per panel

# Every path gets a default. With `set -u`, referencing a config variable the
# user has not defined aborts the whole script -- which is how a missing
# PU_ELIMINATED_BED killed the parascaris and baylisascaris runs instead of
# skipping them. Defaults let `need` report the missing file and carry on.
AS_GFF="${AS_MINIPROT_GFF:-As_miniprot.gff}"
AS_ELIM="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
AS_SIZES="${AS_CHROM_SIZES:-${DB_DIR}/ascaris/v50/AG.v50.fa.fai}"

PU_GFF="${PU_MINIPROT_GFF:-Pu_miniprot.gff}"
PU_ELIM="${PU_ELIMINATED_BED:-${DB_DIR}/parascaris/v1.13/pu_v1.13_eliminated.bed}"
PU_SIZES="${PU_CHROM_SIZES:-${DB_DIR}/parascaris/v1.13/pu_v1.13.masked.fa.fai}"

BSC_TABLE="${BSC_KAKS_TABLE:-kaks_As_Bsc_v2/kaks_filtered.tsv}"
BSC_GFF_F="${BSC_GFF:-../gff3/Bsc_Augustus_V1_July2023_b/Bsc_Augustus_V1.sorted.gff3}"
BSC_SIZES="${BSC_CHROM_SIZES:-${PROJECT_DIR}/../genomes/GWHBAVG00000000.chrom.sizes}"
BSC_ELIM="${BSC_ELIMINATED_BED:-}"

need () {
  local missing=0 f
  for f in "$@"; do
    [[ -s "$f" ]] || { echo "  MISSING: $f" >&2; missing=1; }
  done
  return $missing
}

run_one () {
  local sp="$1" metric="$2"
  local ymax_arg=() common_arg=(--row-height-mm "${ROW_HEIGHT}")
  [[ -n "${YMAX}" ]] && ymax_arg=(--ymax "${YMAX}")
  [[ "${NO_LEGEND}" == "1" ]] && common_arg+=(--no-legend)

  case "$sp" in
    ascaris)
      # 19 autosomes + chrX1..chrX5. Rows are packed to ASC_ROW_MB rather than a
      # fixed count: panel width is proportional to sequence length and Ascaris
      # chromosomes run 6.5-17.5 Mb, so a fixed count gives visibly ragged rows.
      # --break-before '^chrX' keeps the sex chromosomes off the autosome rows.
      need "${TABLE}" "${AS_GFF}" "${AS_ELIM}" "${AS_SIZES}" || return 1
      "${PY}" "${SCRIPT}" \
        --kaks-table "${TABLE}" \
        --gff "${AS_GFF}" \
        --elim-bed "${AS_ELIM}" \
        --chrom-sizes "${AS_SIZES}" \
        --species Ascaris --metric "${metric}" \
        --gene-col gene1 --region-col region1 \
        --window 500 --step 100 --min-genes "${MIN_GENES}" \
        --exclude-chroms '^chrUn' \
        --row-mb "${ASC_ROW_MB}" --break-before '^chrX' \
        "${ymax_arg[@]}" "${common_arg[@]}" \
        --outdir "${RESULTS_DIR}/genome_trace_${STAMP}_ascaris"
      ;;

    parascaris)
      # The germline is ONE 240 Mb chromosome, so it is wrapped across rows.
      # --include-chroms '^chrX$' drops the 50 unplaced fragments (listed, with
      # their genes, under "Sequences not drawn" in the summary).
      # --no-sex-highlight: this chrX is the whole germline genome, not a sex
      # chromosome, and accenting it as one would be wrong. For the same reason
      # the panels say 'Germline chromosome 0-60' and so on, not 'chrX 0-60';
      # the windows table keeps the name chrX.
      need "${TABLE}" "${PU_GFF}" "${PU_ELIM}" "${PU_SIZES}" || return 1
      "${PY}" "${SCRIPT}" \
        --kaks-table "${TABLE}" \
        --gff "${PU_GFF}" \
        --elim-bed "${PU_ELIM}" \
        --chrom-sizes "${PU_SIZES}" \
        --species Parascaris --metric "${metric}" \
        --gene-col gene2 --region-col region2 \
        --window 500 --step 100 --min-genes "${MIN_GENES}" \
        --include-chroms '^chrX$' --no-sex-highlight \
        --label-sub '^chrX$' 'Germline chromosome' \
        --wrap-mb "${PU_WRAP_MB}" --ncols 1 \
        "${ymax_arg[@]}" "${common_arg[@]}" \
        --outdir "${RESULTS_DIR}/genome_trace_${STAMP}_parascaris"
      ;;

    baylisascaris)
      # No eliminated-region calls exist for Baylisascaris. Colouring its genes
      # by their Ascaris ortholog's status would assert a call this genome does
      # not support, so --region-col none: one neutral colour, window line only.
      # That still answers whether dN/dS rises toward chromosome ends.
      #
      # Annotation is Augustus ('transcript' features), not miniprot; the script
      # auto-detects and reports the match count per feature type.
      #
      # Sequence names are 15-character accessions (GWHBAVG00000005), wider
      # than the narrower panels; --label-sub draws them as 'Scaffold 5'. Only
      # the label changes -- the windows table keeps the accession.
      # --min-chrom-mb 5 decides which scaffolds are drawn (no panel cap since
      # 2026-10-12: all 26 of 5 Mb or more); the summary names every sequence
      # left out and why, for the legend. Rows are packed to BSC_ROW_MB.
      need "${BSC_TABLE}" "${BSC_GFF_F}" "${BSC_SIZES}" || return 1
      local elim_arg=() region_arg=(--region-col none)
      if [[ -n "${BSC_ELIM}" && -s "${BSC_ELIM}" ]]; then
        echo "  (Baylisascaris eliminated BED found -- shading and colouring)"
        elim_arg=(--elim-bed "${BSC_ELIM}")
        region_arg=(--region-col region1)
      else
        echo "  (no Baylisascaris eliminated BED -- unshaded, one colour)"
      fi
      "${PY}" "${SCRIPT}" \
        --kaks-table "${BSC_TABLE}" \
        --gff "${BSC_GFF_F}" \
        --chrom-sizes "${BSC_SIZES}" \
        "${elim_arg[@]}" "${region_arg[@]}" \
        --species Baylisascaris --metric "${metric}" \
        --gene-col gene2 \
        --window 500 --step 100 --min-genes "${MIN_GENES}" \
        --min-chrom-mb 5 --row-mb "${BSC_ROW_MB}" \
        --label-sub '^GWHBAVG0*' 'Scaffold ' \
        "${ymax_arg[@]}" "${common_arg[@]}" \
        --outdir "${RESULTS_DIR}/genome_trace_${STAMP}_baylisascaris"
      ;;

    *)
      echo "ERROR: unknown species '$sp' (ascaris|parascaris|baylisascaris|all)" >&2
      return 1
      ;;
  esac
}

[[ "${WHICH}" == "all" ]] && LIST="ascaris parascaris baylisascaris" || LIST="${WHICH}"

for sp in ${LIST}; do
  for m in ${METRICS}; do
    echo ""
    echo "===== ${sp} / ${m} ====="
    run_one "${sp}" "${m}" || echo "  (skipped ${sp} / ${m})"
  done
done

echo ""
echo "===== DONE ====="
for d in "${RESULTS_DIR}"/genome_trace_"${STAMP}"_*; do
  [[ -d "$d" ]] || continue
  echo "--- $(basename "$d") ---"
  for f in "$d"/*_genome_trace_summary.txt; do
    [[ -e "$f" ]] || continue
    echo "  [$(basename "$f" _genome_trace_summary.txt)]"
    sed -n '1,14p' "$f"
    # Dropped, clipped and isolated windows, and label collisions.
    awk '/^Truncation$/{p=1} /^Notes$/{p=0} p' "$f"
    echo ""
  done
done

# Completeness check. The Ka_Ks/Ka/Ks loop is easy to leave half-finished -- a
# walltime kill or a bad path takes out the tail of the loop and the script
# still exits 0, because run_one failures are caught and reported as skips.
echo "===== FIGURE COUNT ====="
for d in "${RESULTS_DIR}"/genome_trace_"${STAMP}"_*; do
  [[ -d "$d" ]] || continue
  n=$(ls "$d"/*.svg 2>/dev/null | wc -l)
  want=$(echo ${METRICS} | wc -w)
  if [[ "${n}" -eq "${want}" ]]; then
    printf '  OK       %-44s %s/%s svg\n' "$(basename "$d")" "${n}" "${want}"
  else
    printf '  INCOMPLETE %-42s %s/%s svg  <-- rerun\n' "$(basename "$d")" "${n}" "${want}"
  fi
done
