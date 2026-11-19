#!/bin/bash
#SBATCH -J jackknife
#SBATCH -A ACF-UTK0146
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=48G
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH -o logs/jackknife_%j.out
#SBATCH -e logs/jackknife_%j.err

# ---------------------------------------------------------------------------
# Leave-one-worm-out jackknife of the piN/piS class contrast.
#
# Usage:
#   sbatch run_jackknife_individuals.sh filtered_miss4        # primary (headline)
#   sbatch run_jackknife_individuals.sh filtered_miss4 3      # clustered robustness
#   sbatch run_jackknife_individuals.sh filtered_miss4 0 100  # main-estimate gene set
#
# Args: callset directory, then --max-local-variants (default 0 = off). MATCH
# THE SECOND ARG TO THE pi_by_region RUN YOU ARE CHECKING:
#   0  primary callset — the headline estimate in Fig 7 and the text
#   3  clustered-variant filter, as in pi_clustfilt — robustness only
# A mismatch jackknifes a different analysis than the one being reported.
#
# Third arg: --min-sites (default 0). The jackknife in the paper (S17, Methods
# "Dependence on individual worms") was run with every gene pooled, which is
# what 0 does: full-cohort fold 1.887. 100 drops genes with fewer than 100
# callable coding sites before pooling, exactly as pi_by_region.py does for the
# main estimate (1.884). Until 2026-10-13 the script accepted --min-sites and
# ignored it, so every earlier run was a 0 run whatever it was passed.
#
# Output: ${POPGEN}/results/jackknife_<SET>            at 0
#         ${POPGEN}/results/jackknife_<SET>_clust<N>   at N > 0
#         ... with _ms<M> appended when --min-sites M is not 0
# so a filtered run cannot overwrite the headline result, and the job refuses
# to start if its output directory already exists. The first run (2026-08-21)
# used 3, when Fig 7 was built from pi_clustfilt; that output was renamed
# jackknife_filtered_miss4_clust3 on 2026-10-08.
#
# ONE pass, not eleven: every subset's pi is recovered from the full-set
# allele counts by subtraction, so the bcftools stream and the codon index are
# paid for once. Time and memory therefore match run_pi_by_region.sh rather
# than eleven times it — same 48G for the same reason (the codon index peaks
# several GB while it is still a Python list), and one CPU because it is
# numpy plus a single query pipe.
#
# WHAT TO READ IN THE LOG:
#   1. REF/codon agreement — same self-check as pi_by_region; if this is not
#      ~0% the coordinates are wrong and nothing below matters.
#   2. The full-cohort fold. It should land near the matching pi_by_region
#      fold (1.88 primary; 2.16 pi_clustfilt, where the clust3 jackknife gave
#      2.135). If it does not, the second arg does not match the analysis
#      being reported. With a third arg of 100 it should EQUAL pi_by_region's
#      (1.884 primary); with 0 it differs in the third decimal (1.887), because
#      the genes under 100 callable coding sites are pooled too.
#   3. The leave-one-out fold range. If its bottom stays above 1.0, no single
#      worm makes the result.
#   4. The LIMIT OF THIS TEST note. A flat jackknife rules out a one-library
#      artifact, NOT a shared one. Run het_by_region.py first; if it names a
#      group, re-run this with --drop-set to remove them together.
#
# The jackknife SE it reports is over INDIVIDUALS, unlike the gene bootstrap
# in pi_by_region.py. That is the number to quote at a reviewer asking whether
# ten worms is enough.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Configuration. Every machine-specific path and the conda setup come from
# config.sh at the top of the repository (copy config.example.sh). Under sbatch
# the batch script runs from a spool copy, so config.sh is searched for upward
# from this file, from the submission directory and from $PWD;
# KAKS_CONFIG=/full/path/to/config.sh overrides the search.
# ---------------------------------------------------------------------------
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
  echo "ERROR: no config.sh found (copy config.example.sh to config.sh)." >&2
  echo "  searched upward from: ${HERE}" >&2
  echo "                        ${SLURM_SUBMIT_DIR:-<not under sbatch>}" >&2
  echo "                        ${PWD}" >&2
  echo "  Set KAKS_CONFIG=/full/path/to/config.sh and resubmit." >&2
  exit 1
fi
echo "config: ${CONFIG}"
# shellcheck source=../config.example.sh
source "${CONFIG}"
REPO_ROOT="$(cd "$(dirname "${CONFIG}")" && pwd)"
REQUIRED=(PROJECT_DIR POPGEN_DIR POPGEN_WORK_DIR AS_V50_DIR)
for v in "${REQUIRED[@]}"; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: ${v} is not set in ${CONFIG}" >&2; exit 1; }
done

find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/07_polymorphism" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script jackknife_individuals.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: jackknife_individuals.py not found near ${REPO_ROOT}/07_polymorphism." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${POPGEN_CONDA_ENV:?set POPGEN_CONDA_ENV in config.sh}"
  [[ -n "${POPGEN_BIN_PREPEND:-}" ]] && export PATH="${POPGEN_BIN_PREPEND}:${PATH}"
fi

cd "${POPGEN_WORK_DIR}" || exit 1
mkdir -p logs

SET="${1:-filtered_miss4}"
CLUST="${2:-0}"
MIN_SITES="${3:-0}"
if [ "${CLUST}" = "0" ]; then SUFFIX=""; else SUFFIX="_clust${CLUST}"; fi
if [ "${MIN_SITES}" != "0" ]; then SUFFIX="${SUFFIX}_ms${MIN_SITES}"; fi
POPGEN="${POPGEN_DIR}"
V50="${AS_V50_DIR}"
MINIPROT="${PROJECT_DIR}"
OUTDIR="${POPGEN}/results/jackknife_${SET}${SUFFIX}"

# jackknife_individuals.py has no make_outdir(), so a rerun into an existing
# directory would overwrite the run it came from. Refuse instead.
if [ -e "${OUTDIR}" ]; then
  echo "ERROR: ${OUTDIR} already exists. Move it aside before rerunning." >&2
  exit 1
fi

# Record the inputs in the log; the August logs did not, which is why the
# callset behind each result had to be reconstructed from wrapper defaults.
echo "[$(date '+%F %T')] SET=${SET} CLUST=${CLUST} MIN_SITES=${MIN_SITES} OUTDIR=${OUTDIR}"

# Three explicit drop-sets, chosen from the het_by_region.py result rather
# than guessed. Only the FIRST tests the threat to the headline.
#
#   1. SRR31757523 + SRR31757531 — the only two libraries whose eliminated
#      heterozygosity EXCEEDS their retained (eli/ret 1.18 and 1.16 against a
#      cohort median of 0.85). That is the anti-conservative direction: excess
#      het in the eliminated class inflates eliminated pi and pushes the fold
#      up on its own. If the fold survives dropping both, the collapsed-
#      paralogy route to a false positive is closed.
#
#   2. The whole SRR31757xxx block. 523 and 531 hold the cohort's two lowest
#      RETAINED het rates, which is what makes their ratios exceed 1 — so this
#      may be a batch property rather than two independent libraries. A batch
#      artifact is invisible to leave-one-out by construction, so it has to be
#      removed as a group to be tested at all.
#
#   3. SRR31692059 + SRR31675198 — the two libraries with the global IBS0
#      deficit and the highest absolute het (0.52, 0.54). Their eli/ret ratios
#      are 0.66 and 0.85, i.e. at or BELOW the median, so their excess het sits
#      in retained and works AGAINST the result. Dropped anyway to show the
#      fold does not depend on them either.
python "${SCRIPT}" \
  --vcf "${POPGEN}/vcf/${SET}/snps_invariant.bcf" \
  --genome "${V50}/AG.v50.fa" \
  --gff "${MINIPROT}/As_miniprot.gff" \
  --gene-callable "${POPGEN}/targets/gene_callable.tsv" \
  --max-local-variants "${CLUST}" \
  --local-window 150 \
  --min-sites "${MIN_SITES}" \
  --drop-set SRR31757523,SRR31757531 \
  --drop-set SRR31757523,SRR31757530,SRR31757531 \
  --drop-set SRR31692059,SRR31675198 \
  --outdir "${OUTDIR}"
