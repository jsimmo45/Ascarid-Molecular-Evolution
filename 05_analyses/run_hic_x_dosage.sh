#!/bin/bash
#SBATCH -J hic_x_dosage
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=hic_x_dosage_%j.out
#SBATCH --error=hic_x_dosage_%j.err
#
# Do Ascaris males carry one copy of each X? Read depth of the testis Hi-C
# library against the ovary library, per chromosome. Every X should sit at
# ~0.5x the ovary depth and every autosome at ~1x. Backs the male-biased-
# mutation sentence in the Discussion (results draft, note 24) with libraries
# whose sex is known from dissection.
#
#   sbatch run_hic_x_dosage.sh
#   sbatch run_hic_x_dosage.sh <testis pairs file> <ovary pairs file>
#
# INPUTS, first match wins:
#   1. arguments 1 and 2 (one testis file, one ovary file);
#   2. HIC_MALE_PAIRS and HIC_FEMALE_PAIRS from config.sh, space-separated
#      (replicates of one sex are pooled);
#   3. auto-discovery: HiC-Pro allValidPairs (or .pairs), plain or gzipped,
#      under ASCARIS_HIC_DIR (default below) whose names contain teste/testis
#      and ovary/ovari. It needs EXACTLY ONE file per sex, and otherwise stops
#      and lists what it found, because a repeat-masked or downsampled copy
#      beside the primary file must not be picked by accident.
# Assembly: AS_CHROM_SIZES and AS_ELIMINATED_BED from config.sh, defaulting to
# the v50 germline files genome_trace uses. All Hi-C is in germline coordinates.
#
# CORES: 2, one per library. The script reads each file in one pass and caps
# its workers at the number of libraries, so a third core would sit idle.
# TIME: about 1-2 M pairs per second per worker on plain text (slower if
# gzipped), so a few minutes per 100 M pairs; 1 h is generous.

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
SCRIPT="$(find_script hic_x_dosage.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: hic_x_dosage.py not found near ${REPO_ROOT}." >&2; exit 1; }
echo "script: ${SCRIPT}"

SCRIPT_PARENT="$(dirname "${SCRIPT}")"
export PYTHONPATH="${SCRIPT_PARENT}:${REPO_ROOT}/lib${PYTHONPATH:+:${PYTHONPATH}}"

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

# Probe the imports: region_lib lives in lib/.
"${PY}" - <<'PYEOF' || { echo "ERROR: fix PYTHONPATH (lib/ must be on it) and resubmit." >&2; exit 1; }
import importlib.util, sys
missing = [m for m in ("region_lib", "numpy", "pandas")
           if importlib.util.find_spec(m) is None]
if missing:
    print("ERROR: not importable: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
import pandas
print("modules: region_lib, numpy, pandas %s all importable" % pandas.__version__)
PYEOF

cd "${PROJECT_DIR}"

# --- inputs -------------------------------------------------------------------
SIZES="${AS_CHROM_SIZES:-${DB_DIR}/ascaris/v50/AG.v50.fa.fai}"
ELIM="${AS_ELIMINATED_BED:-${DB_DIR}/ascaris/v50/AG.v50.eliminated_strict.bed}"
HIC_DIR="${ASCARIS_HIC_DIR:?set ASCARIS_HIC_DIR in config.sh}"

discover() {                        # discover <root> <name fragment>...
  local root="$1"; shift
  local expr=() p
  for p in "$@"; do
    [[ ${#expr[@]} -gt 0 ]] && expr+=(-o)
    expr+=(-iname "*${p}*")
  done
  find "${root}" -maxdepth "${HIC_FIND_DEPTH:-6}" \( -type f -o -type l \) \
       \( -name '*allValidPairs' -o -name '*allValidPairs.gz' \
          -o -name '*.pairs' -o -name '*.pairs.gz' \) \
       \( "${expr[@]}" \) 2>/dev/null | sort
}

MALE=(); FEMALE=()
if [[ -n "${1:-}" || -n "${2:-}" ]]; then
  [[ -n "${1:-}" && -n "${2:-}" ]] || {
    echo "ERROR: give both files: sbatch run_hic_x_dosage.sh <testis> <ovary>" >&2; exit 1; }
  MALE=("$1"); FEMALE=("$2"); SOURCE="arguments"
elif [[ -n "${HIC_MALE_PAIRS:-}" && -n "${HIC_FEMALE_PAIRS:-}" ]]; then
  read -r -a MALE <<< "${HIC_MALE_PAIRS}"
  read -r -a FEMALE <<< "${HIC_FEMALE_PAIRS}"
  SOURCE="config.sh"
else
  [[ -d "${HIC_DIR}" ]] || {
    echo "ERROR: no inputs given and ASCARIS_HIC_DIR does not exist: ${HIC_DIR}" >&2; exit 1; }
  echo "Discovering pairs files under ${HIC_DIR} (depth ${HIC_FIND_DEPTH:-6}) ..."
  mapfile -t MC < <(discover "${HIC_DIR}" teste testis)
  mapfile -t FC < <(discover "${HIC_DIR}" ovary ovari)
  if (( ${#MC[@]} != 1 || ${#FC[@]} != 1 )); then
    echo "ERROR: found ${#MC[@]} testis and ${#FC[@]} ovary pairs files; need exactly one of each." >&2
    for f in "${MC[@]}"; do echo "  testis candidate: $(ls -lhL "$f" | awk '{print $5}')  $f" >&2; done
    for f in "${FC[@]}"; do echo "  ovary candidate : $(ls -lhL "$f" | awk '{print $5}')  $f" >&2; done
    echo "  Pick the primary (unmasked, full-depth) files and pass them:" >&2
    echo "    sbatch run_hic_x_dosage.sh <testis.allValidPairs> <ovary.allValidPairs>" >&2
    echo "  or set HIC_MALE_PAIRS and HIC_FEMALE_PAIRS in config.sh (replicates are pooled)." >&2
    echo "  If nothing was found, set ASCARIS_HIC_DIR in config.sh to the HiC-Pro output." >&2
    exit 1
  fi
  MALE=("${MC[0]}"); FEMALE=("${FC[0]}"); SOURCE="auto-discovered"
fi

for f in "${MALE[@]}" "${FEMALE[@]}" "${SIZES}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done
ELIM_ARGS=()
if [[ -s "${ELIM}" ]]; then
  ELIM_ARGS=(--eliminated-bed "${ELIM}")
else
  echo "ERROR: eliminated BED missing or empty: ${ELIM}" >&2
  echo "  Without it, somatic carryover in either gonad could read as a dosage" >&2
  echo "  difference. Set AS_ELIMINATED_BED in config.sh." >&2
  exit 1
fi

echo "inputs (${SOURCE}):"
for f in "${MALE[@]}"; do echo "  testis: $(ls -lhL "$f" | awk '{print $5}')  $f"; done
for f in "${FEMALE[@]}"; do echo "  ovary : $(ls -lhL "$f" | awk '{print $5}')  $f"; done
echo "  sizes : ${SIZES}"
echo "  elim  : ${ELIM}"

STAMP="${STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUTDIR="${RESULTS_DIR}/hic_x_dosage_${STAMP}"
echo "outdir: ${OUTDIR}"

"${PY}" "${SCRIPT}" \
  --male "${MALE[@]}" \
  --female "${FEMALE[@]}" \
  --chrom-sizes "${SIZES}" \
  "${ELIM_ARGS[@]}" \
  --binsize 100000 \
  --min-mapq 30 \
  --min-frac 0.2 \
  --max-fold 5 \
  --min-run 5 \
  --workers "${SLURM_CPUS_PER_TASK:-1}" \
  --outdir "${OUTDIR}"

echo ""
echo "===== DONE ====="
cat "${OUTDIR}/hic_x_dosage_summary.txt"
