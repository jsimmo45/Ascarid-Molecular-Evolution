# config.example.sh — copy to config.sh and edit for your machine.
#
#   cp config.example.sh config.sh
#
# config.sh is gitignored. Every run_*.sh finds it by searching upward from the
# folder it is run or submitted from (set KAKS_CONFIG=/full/path/to/config.sh
# to use one elsewhere), so the paths and environment below are the only
# machine-specific settings in the repository. A wrapper that needs a setting
# you have not filled in stops and names it.
#
# File names below are the ones the published runs used; README.md says where
# each input comes from.

# ---------------------------------------------------------------------------
# Environments (envs/kaks.yml, envs/popgen.yml)
# ---------------------------------------------------------------------------
# conda's activation script. Leave empty to skip conda and use the python and
# tools already on PATH.
CONDA_SH="/path/to/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="kaks"              # 01_gene_sets to 06_outputs
POPGEN_CONDA_ENV="popgen"     # 07_polymorphism

# Optional bin directories put first on PATH after activation. Needed only where
# a system or module python shadows the conda one (as on the cluster the paper
# used, where these were the two environments' bin/ folders). Leave empty
# otherwise.
CONDA_BIN_PREPEND=""
POPGEN_BIN_PREPEND=""

# ---------------------------------------------------------------------------
# Working directories
# ---------------------------------------------------------------------------
# 01-06: gene sets (*_prot.fa, *_cds.fa, *_miniprot.gff), ortholog and rate
# tables (rbh_*/, kaks_*/) and the dated run directories under results/.
PROJECT_DIR="/path/to/work"

# Root for dated run directories: each wrapper writes
# ${RESULTS_DIR}/<analysis>_<YYYYMMDD_HHMMSS>/, so a rerun never overwrites the
# run a figure was built from.
RESULTS_DIR="${PROJECT_DIR}/results"

# Directory holding the .py scripts. Leave empty: each wrapper finds its .py
# beside itself, which is what the repository layout gives you.
SCRIPT_DIR=""

# 07: reads, BAMs, VCFs, callable targets and results/ of the population sample.
POPGEN_DIR="/path/to/popgen"
# The folder the 07 wrappers run in (sample-selection folders such as
# cohort_20260820/ and supp_audit/, logs/).
POPGEN_WORK_DIR="${PROJECT_DIR}/07_polymorphism"

# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------
DB_DIR="/path/to/db"

# Ascaris germline assembly AG.v50 and its retained and eliminated regions
AS_V50_DIR="${DB_DIR}/ascaris/v50"
AS_GENOME_FA="${AS_V50_DIR}/AG.v50.fa"
AS_CHROM_SIZES="${AS_V50_DIR}/AG.v50.fa.fai"
AS_ELIMINATED_BED="${AS_V50_DIR}/AG.v50.eliminated_strict.bed"
AS_RETAINED_BED="${AS_V50_DIR}/AG.v50.retained_loose.bed"
# The previous assembly (v37), on which the expression table's gene models sit
AS_V37_GENOME_FA="${DB_DIR}/ascaris/v37/AG.genome.fa"
# Germline-to-somatic coordinates, to place genes on the somatic Hi-C map
GERM_TO_SOMA_MAP="${DB_DIR}/ascaris/v50_postPDE/prepde_to_postpde_coordinates.bed"

# P. univalens germline assembly (pu_v3; file names carry pu_v1.13)
PU_ELIMINATED_BED="${DB_DIR}/parascaris/v1.13/pu_v1.13_eliminated.bed"
PU_RETAINED_BED="${DB_DIR}/parascaris/v1.13/somatic_genome/germ_to_soma_mapping.bed"
PU_CHROM_SIZES="${DB_DIR}/parascaris/v1.13/pu_v1.13.masked.fa.fai"

# T. canis eliminated transcripts (Wang et al. 2017)
TOX_ELIMINATED_BED="${DB_DIR}/toxocara/final_release/eliminated.transcripts.longest.bed"

# B. schroederi Augustus CDS, translated by 01_gene_sets/run_translate_cds.sh
BSC_CDS_FA="/path/to/Bsc_Augustus_V1.sorted.fasta"
# B. schroederi GFF and chromosome sizes for 05_analyses/run_genome_trace.sh.
# Unset, that wrapper looks for them beside PROJECT_DIR
# (../gff3/Bsc_Augustus_V1_July2023_b/ and ../genomes/).
# BSC_GFF="/path/to/Bsc_Augustus_V1.sorted.gff3"
# BSC_CHROM_SIZES="/path/to/GWHBAVG00000000.chrom.sizes"

# ---------------------------------------------------------------------------
# Hi-C products of the companion study (Simmons et al. 2026; GEO GSE314626)
# ---------------------------------------------------------------------------
# FAN-C compartment eigenvectors at 100 kb
HIC_COMPARTMENT_DIR="/path/to/compartment_analysis/ascaris/fanc_hic_files_100kb"
EIGENVECTOR_DIR="${HIC_COMPARTMENT_DIR}/prepde/eigenvectors"
SOMATIC_EIGENVECTOR="${HIC_COMPARTMENT_DIR}/postpde/eigenvectors/as_10day_iced_100kb.matrix.eigenvector"
# Hi-C pairs and insulation scores
ASCARIS_HIC_DIR="/path/to/ascaris_hic_reads"
INSULATION_DIR="${ASCARIS_HIC_DIR}/insulation_score_updated/insulation_scores/prepde"

# ---------------------------------------------------------------------------
# Pairwise Ka/Ks tables — paths relative to PROJECT_DIR
#
# The _v2 suffix marks tables regenerated after the BLAST identifier fix. Any
# Baylisascaris table without it is pre-fix and is roughly half the true size.
# ---------------------------------------------------------------------------
KAKS_AS_PU="kaks_results/kaks_filtered.tsv"
KAKS_AS_BSC="kaks_As_Bsc_v2/kaks_filtered.tsv"
KAKS_PU_BSC="kaks_Pu_Bsc_v2/kaks_filtered.tsv"
KAKS_TOX_AS="kaks_Tox_As/kaks_filtered.tsv"
KAKS_TOX_PU="kaks_Tox_Pu/kaks_filtered.tsv"
KAKS_TOX_BSC="kaks_Tox_Bsc_v2/kaks_filtered.tsv"

# Annotated rate table written by 05_analyses/run_annotate_regions.sh: per-gene
# dN, dS and dN/dS with region calls and CBR distances, relative to PROJECT_DIR.
# Most downstream analyses read it; that wrapper prints the line to paste here.
ANNOTATED_TABLE=""

# ---------------------------------------------------------------------------
# SLURM
# ---------------------------------------------------------------------------
# The wrappers' #SBATCH lines name the account, partition and QOS of the cluster
# the paper used. They are not read from here; override them when you submit:
#   sbatch -A <account> -p <partition> --qos=<qos> run_x.sh
# Run without SLURM with bash, setting what the wrapper would get from SLURM:
#   SLURM_CPUS_PER_TASK=8 bash run_x.sh
#   SLURM_ARRAY_TASK_ID=3 SLURM_CPUS_PER_TASK=8 bash run_map_reads.sh   # array jobs
