#!/bin/bash
#SBATCH -J gs_rerun
#SBATCH -A ACF-UTK0146
#SBATCH --partition=condo-jwang166
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=gs_rerun_%j.out
#SBATCH --error=gs_rerun_%j.err
#
# Rerun gene_structure.py (Figure 2D) on exactly the inputs of an earlier run.
#
#   sbatch run_gene_structure_rerun.sh [SOURCE_RUN] [NEW_RUN]
#
#   SOURCE_RUN  folder under results/ whose provenance.json supplies the inputs
#               and arguments (default gene_structure_20260817_loci, the
#               current Fig 2D and the source of table S1's gene columns)
#   NEW_RUN     new folder under results/ (default gene_structure_<today>_loci);
#               the job refuses to write into a folder that already exists
#
# Written for the 2026-10-13 changes to gene_structure.py, which touch the
# figure and the summary text only. After the run the job checks, and prints:
#   - gene_structure_per_gene.tsv and gene_structure_summary.tsv against the
#     source run's: byte-identical, or failing that equal to 1e-12 (a different
#     pandas can print a mean's last digit differently), or DIFFERENT;
#   - every input's checksum against the one recorded in the source run;
#   - the new legend-notes block, for the Fig 2D legend.
# Repoint Fig 2D only if the tables match and the inputs are unchanged.
#
# Single-threaded pandas/matplotlib: 1 node / 1 task / 1 CPU; a few minutes.

set -euo pipefail

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
REQUIRED=(PROJECT_DIR)
for v in "${REQUIRED[@]}"; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: ${v} is not set in ${CONFIG}" >&2; exit 1; }
done

find_script() {
  local c
  for c in "${SCRIPT_DIR:-}" "${HERE}" "${REPO_ROOT}/05_analyses" "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$c" && -f "${c}/$1" ]] && { printf '%s' "${c}/$1"; return 0; }
  done
  return 1
}
SCRIPT="$(find_script gene_structure.py || true)"
[[ -n "${SCRIPT}" ]] || { echo "ERROR: gene_structure.py not found near ${REPO_ROOT}/05_analyses." >&2; exit 1; }
echo "script: ${SCRIPT}"

if [[ -n "${CONDA_SH:-}" && -f "${CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV:?set CONDA_ENV in config.sh}"
  [[ -n "${CONDA_BIN_PREPEND:-}" ]] && export PATH="${CONDA_BIN_PREPEND}:${PATH}"
fi

MINIPROT="${MINIPROT:-${PROJECT_DIR}}"
RESULTS="${MINIPROT}/results"
SRC="${1:-gene_structure_20260817_loci}"
NEW="${2:-gene_structure_$(date +%Y%m%d)_loci}"
SRC_DIR="${RESULTS}/${SRC}"
NEW_DIR="${RESULTS}/${NEW}"

[[ -f "${SCRIPT}" ]] || { echo "ERROR: ${SCRIPT} not found." >&2; exit 1; }
[[ -f "${SRC_DIR}/provenance.json" ]] || {
  echo "ERROR: no provenance.json in ${SRC_DIR}." >&2; exit 1; }
[[ -e "${NEW_DIR}" ]] && {
  echo "ERROR: ${NEW_DIR} already exists; pass another NEW_RUN." >&2; exit 1; }

# The earlier run wrote relative paths from here (outdir results/...), so a
# relative path in its provenance resolves against this directory.
cd "${MINIPROT}" || exit 1
echo "[$(date '+%F %T')] ${SRC} -> ${NEW}"
echo "script: ${SCRIPT}"
echo "python: $(command -v python)"

# --- the earlier run's inputs and arguments ---------------------------------
ARGFILE="$(mktemp)"
trap 'rm -f "${ARGFILE}"' EXIT
python - "${SRC_DIR}/provenance.json" > "${ARGFILE}" <<'EOF'
import json, os, sys
prov = json.load(open(sys.argv[1]))
args, inputs = prov.get('args', {}), prov.get('inputs', {})
out = []
for key, flag in (('cds_fasta', '--cds-fasta'), ('gff', '--gff'),
                  ('rates_table', '--rates-table')):
    path = (inputs.get(key) or {}).get('path') or args.get(key)
    if not path or not os.path.isfile(path):
        sys.exit(f'ERROR: {key} recorded in the provenance not found: {path}')
    out += [flag, path]
for key, flag in (('gene_col', '--gene-col'), ('region_col', '--region-col'),
                  ('box_alpha', '--box-alpha'), ('cap_pct', '--cap-pct'),
                  ('seed', '--seed'), ('width', '--width'),
                  ('height_mm', '--height-mm')):
    if args.get(key) is not None:
        out += [flag, str(args[key])]
if args.get('panel_letters'):
    out.append('--panel-letters')
print('\n'.join(out))
EOF
mapfile -t ARGS < "${ARGFILE}"
echo "arguments from ${SRC}/provenance.json:"
printf '  %s\n' "${ARGS[@]}"

python "${SCRIPT}" "${ARGS[@]}" --outdir "${NEW_DIR}"

# --- checks -----------------------------------------------------------------
echo
echo "=== checks against ${SRC} ==="
python - "${SRC_DIR}" "${NEW_DIR}" <<'EOF'
import filecmp, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

src, new = Path(sys.argv[1]), Path(sys.argv[2])
ok = True
for name in ('gene_structure_per_gene.tsv', 'gene_structure_summary.tsv'):
    a, b = src / name, new / name
    if filecmp.cmp(a, b, shallow=False):
        print(f'IDENTICAL       {name}')
        continue
    x, y = pd.read_csv(a, sep='\t'), pd.read_csv(b, sep='\t')
    same = list(x.columns) == list(y.columns) and x.shape == y.shape
    if same:
        for c in x.columns:
            if pd.api.types.is_numeric_dtype(x[c]) and pd.api.types.is_numeric_dtype(y[c]):
                same &= bool(np.allclose(x[c], y[c], rtol=1e-12, atol=0,
                                         equal_nan=True))
            else:
                same &= bool((x[c].fillna('NA').astype(str)
                              == y[c].fillna('NA').astype(str)).all())
    print(f"{'EQUAL TO 1e-12' if same else 'DIFFERENT':<15} {name}")
    ok &= same
old = json.loads((src / 'provenance.json').read_text())
cur = json.loads((new / 'provenance.json').read_text())
unrecorded = False
for k in ('cds_fasta', 'gff', 'rates_table'):
    o = (old.get('inputs', {}).get(k) or {}).get('sha256_16')
    n = (cur.get('inputs', {}).get(k) or {}).get('sha256_16')
    if o is None:
        # An older provenance without checksums: the tables are the check.
        print(f"{'NOT RECORDED':<15} {k:<12} -> {n}")
        unrecorded = True
        continue
    print(f"{'SAME INPUT' if o == n else 'INPUT CHANGED':<15} {k:<12} {o} -> {n}")
    ok &= (o == n)
if unrecorded:
    print('(the source run recorded no checksum for some inputs, so only the '
          'tables confirm them)')
print('packages then', old.get('packages'), '/ now', cur.get('packages'))
print()
print('VERDICT: ' + ('figure-only change confirmed; repoint Fig 2D to '
                     f'{new.name}/gene_structure_panel.svg (table S1 can keep '
                     f'{src.name})' if ok else
                     'NOT a figure-only rerun; do not repoint until the '
                     'difference is explained'))
EOF

echo
sed -n '/^Legend notes/,/^Reporting notes/p' "${NEW_DIR}/gene_structure_summary.txt" | sed '$d'
echo "[$(date '+%F %T')] done: ${NEW_DIR}"
