# Ascarid-Molecular-Evolution

Code for the molecular evolution of genes removed by programmed DNA elimination
in four ascarid nematodes (*Ascaris*, *Parascaris univalens*, *Baylisascaris
schroederi* and *Toxocara canis*): pairwise and branch-specific dN/dS, constraint
around chromosomal breakage regions (CBRs), codon usage, gene structure,
chromatin context and polymorphism in an *Ascaris* population sample.

Every analysis is a Python script with a command-line interface
(`python script.py --help`) and a `run_*.sh` wrapper that holds the arguments
used for the paper. Each run writes a dated folder with a `provenance.json`
recording its inputs and arguments.

## Layout

| Folder | Contents |
|---|---|
| `lib/` | shared modules: figure style (`paper_colors.py`), region and CBR definitions (`region_lib.py`), table readers |
| `01_gene_sets/` | translation of the *B. schroederi* gene models |
| `02_orthologs/` | 1:1 orthologs by reciprocal best BLAST hits, six species pairs |
| `03_rates/` | pairwise dN and dS (MAFFT, yn00); four-species alignments, codeml branch models, clock tests |
| `04_annotation/` | region labels on the *T. canis* pairwise tables |
| `05_analyses/` | region annotation, the analyses behind each figure and table, figure and table assembly |
| `06_outputs/` | rate and additivity table |
| `07_polymorphism/` | population sample: download, mapping, germline screen, genotyping, diversity, figures |
| `envs/` | conda environments |
| `docs/code_manifest.tsv` | every file, where it came from on the analysis server, and what it produces |

## Setup

1. Create the environments (built for Linux, linux-64):

   ```bash
   conda env create -f envs/kaks.yml      # 01_gene_sets to 06_outputs
   conda env create -f envs/popgen.yml    # 07_polymorphism
   ```

   `envs/kaks.exact.yml` and `envs/popgen.exact.yml` list every package of the
   two environments as used for the paper.

2. Copy the settings template and fill in the paths:

   ```bash
   cp config.example.sh config.sh
   ```

   `config.sh` holds every machine-specific setting: the conda activation script,
   the working directories (`PROJECT_DIR`, `POPGEN_DIR`), the reference data
   (`DB_DIR` and the files under it) and the Hi-C products. Each wrapper finds
   it by searching upward from where it is run; a setting it needs but cannot
   find stops it with the setting's name.

3. Put the inputs in place (below) and run the wrappers in order.

### SLURM, or not

The `#SBATCH` lines name the account, partition and QOS of the cluster the paper
used. Override them when you submit:

```bash
cd 03_rates
sbatch -A <account> -p <partition> --qos=<qos> run_kaks_toxocara.sh
```

Without SLURM, run a wrapper with `bash`, setting what SLURM would provide:

```bash
SLURM_CPUS_PER_TASK=8 bash run_pi_by_region.sh
SLURM_ARRAY_TASK_ID=3 SLURM_CPUS_PER_TASK=8 bash run_map_reads.sh     # array jobs
```

The `07_polymorphism` wrappers write their SLURM logs to `logs/` in the folder
they are submitted from, so submit them from `07_polymorphism/` (which has one).

## Inputs

| Input | Setting in `config.sh` | Source |
|---|---|---|
| *Ascaris* germline assembly AG.v50 | `AS_GENOME_FA`, `AS_CHROM_SIZES` | https://dnaelimination.utk.edu/protocols-data/ (Simmons et al. 2026) |
| *Ascaris* retained and eliminated regions | `AS_RETAINED_BED`, `AS_ELIMINATED_BED` | Wang et al. 2020, transferred to AG.v50 by Simmons et al. 2026 |
| *Ascaris* germline-to-somatic coordinates | `GERM_TO_SOMA_MAP` | Simmons et al. 2026 |
| *Ascaris* previous assembly (v37) | `AS_V37_GENOME_FA` | the assembly the expression table's gene models are on (Wang et al. 2020) |
| *P. univalens* germline assembly pu_v3 and its germline-to-somatic mapping | `PU_RETAINED_BED`, `PU_ELIMINATED_BED`, `PU_CHROM_SIZES` | https://dnaelimination.utk.edu/protocols-data/; pu_v3 differs from the published assembly (BioProject PRJNA1067589; Simmons et al. 2024) only in the orientation of somatic chromosome 1. File names carry pu_v1.13. |
| *T. canis* eliminated transcripts | `TOX_ELIMINATED_BED` | Wang et al. 2017 (GEO GSE99523) |
| *B. schroederi* gene models (Augustus v3.4.0) | `BSC_CDS_FA`, `BSC_GFF` | predicted for the comparative analyses of Simmons et al. 2024 on assembly GWHBAVG00000000 (Xie et al. 2022), before that study merged them with RNA-seq-based models; not the annotation published with the assembly; to be deposited with the publication |
| *Ascaris* gene annotation and expression table | `EXPRESSION_XLSX` (default `Ascaris_genes_15714_annotation_expression.xlsx` in `PROJECT_DIR`) | gene models and nine-stage expression summary of Wang et al. 2020; to be deposited with the publication |
| Hi-C compartment eigenvectors, insulation scores, pairs | `HIC_COMPARTMENT_DIR`, `ASCARIS_HIC_DIR` and the files under them | GEO GSE314626, processed as in Simmons et al. 2026 (HiC-Pro v3.1.0, FAN-C v0.9.28) |
| Population sample, *A. lumbricoides* | `POPGEN_DIR` (reads/) | BioProject PRJNA1195031 (Landeryou et al. 2025); `07_polymorphism/run_fetch_ena_reads.sh` |
| Runs of the population sample to download | `POPGEN_WORK_DIR` (copy `07_polymorphism/supp_audit/` there) | `analysed_runs.txt`: the 54 runs analysed by Landeryou et al. 2025, matched to SRA runs by read count because the accessions printed in their supplementary table do not all agree (`supplementary_audit.tsv`); `runs_not_in_paper.txt`: the BioProject's other eight runs |
| Germline and somatic tissue controls | `POPGEN_DIR` (controls/) | BioProject PRJNA62057 (SRR502751, SRR502756, SRR5573034); `07_polymorphism/fetch_ena_reads_wang.py` |

## Gene sets

Made in `PROJECT_DIR` before step 2.

* *Ascaris* and *P. univalens*: previously annotated proteins, one per locus
  (the longest isoform), mapped onto AG.v50 and pu_v3 with miniprot
  v0.18-r281, with coding sequences and proteins extracted by gffread v0.12.7:
  `As_miniprot.gff`, `As_cds.fa`, `As_prot.fa` (18,588 models) and
  `Pu_miniprot.gff`, `Pu_cds.fa`, `Pu_prot.fa` (19,143 models). The
  *P. univalens* proteins (`Pu_largest_tx.16603.faa`) are the 16,603 gene
  models of Simmons et al. 2024; the *Ascaris* proteins
  (`As_largest_tx.17658.faa`) are an unpublished 2023 laboratory annotation of
  17,658 genes on the assembly of Wang et al. 2020. Both will be deposited
  with the publication.

  ```bash
  miniprot -t 48 --gff AG.v50.fa As_largest_tx.17658.faa > As_miniprot_raw.gff
  miniprot -t 48 --gff pu_v1.13.masked.fa Pu_largest_tx.16603.faa > Pu_miniprot_raw.gff
  # gffread cannot parse the ##PAF lines miniprot writes into its GFF
  grep -v '^##PAF' As_miniprot_raw.gff > As_miniprot.gff
  grep -v '^##PAF' Pu_miniprot_raw.gff > Pu_miniprot.gff
  gffread As_miniprot.gff -g AG.v50.fa -x As_cds.fa -y As_prot.fa
  gffread Pu_miniprot.gff -g pu_v1.13.masked.fa -x Pu_cds.fa -y Pu_prot.fa
  ```

* *T. canis*: open reading frames called with TransDecoder v5.5.0 on the
  StringTie transcripts of Wang et al. 2017, saved as `Tox_cds.fa` and
  `Tox_prot.fa` (13,835 ORFs):

  ```bash
  TransDecoder.LongOrfs -t TG.transcripts.longest.20065.fasta
  TransDecoder.Predict  -t TG.transcripts.longest.20065.fasta --single_best_only
  ```

* *B. schroederi*: `01_gene_sets/run_translate_cds.sh` writes `Bsc_prot.fa` and
  `Bsc_cds.fa`.

## Run order

Each wrapper's header lists its inputs and outputs; `docs/code_manifest.tsv`
maps every script to the figures and tables it produces.

1. **Gene sets**: `01_gene_sets/run_translate_cds.sh`.
2. **Orthologs** (`02_orthologs/`): `run_rbh.sh` (*Ascaris* × *P. univalens*),
   `run_rbh_bsc.sh` (*Ascaris* and *P. univalens* × *B. schroederi*),
   `run_rbh_toxocara.sh` (*T. canis* × the other three), then
   `run_rbh_bsc_fix.sh`, which re-merges the three *B. schroederi* pairs into
   the reported `rbh_*_Bsc_v2/` sets (an identifier bug in an earlier version
   dropped *B. schroederi* genes above g9999; the current
   `reciprocal_best_hits.py` contains the fix).
3. **Pairwise rates** (`03_rates/`): `run_kaks_pipeline_2.sh` (*Ascaris* ×
   *P. univalens*), `run_kaks_toxocara.sh`, `run_kaks_bsc_v2.sh`.
4. **Region labels**: `05_analyses/run_annotate_regions.sh`, then set
   `ANNOTATED_TABLE` in `config.sh` to the table it writes (it prints the line);
   `04_annotation/run_annotate_tox_kaks.sh`.
5. **Four-species analyses**: `03_rates/run_branch_dnds_four.sh`,
   `run_build_quartet_alignments.sh`, `run_molecular_clock_test.sh`;
   `06_outputs/run_make_rate_table.sh`.
6. **Figures and tables** (`05_analyses/`): the analysis wrappers below, then
   `run_make_supplementary_tables.sh`, `run_assemble_main_figures.sh` and
   `run_assemble_supplementary_figures.sh`.
7. **Polymorphism** (`07_polymorphism/`): `run_fetch_ena_reads.sh`;
   `run_map_reads.sh` and `run_map_controls.sh` (screening subsample);
   `run_germline_fraction_audit.sh`; `select_germline_samples.py` (writes the
   cohort folder with `selected_runs.txt`); `run_map_full.sh <cohort>`; `run_mappability_mask.sh`;
   `make_callable_targets.py`; `run_depth_match.sh <cohort> <depth>`;
   `run_joint_genotype.sh` and its merge step (in its header);
   `run_filter_genotypes.sh`; `run_pi_by_region.sh`; then the remaining
   wrappers below.

Some wrappers name the run folders the published figures were built from (for
example the default input of `run_figure_polymorphism.sh`). After rerunning
the steps upstream, pass your own run folder as the argument.

### Figures and tables

Most scripts have a wrapper `run_<script>.sh` beside them; `select_germline_samples.py`,
`make_callable_targets.py` and `fetch_ena_reads_wang.py` are run directly.
Figure numbers written inside the scripts' comments follow earlier drafts of
the paper; this table gives the published numbering.

| Figure or table | Script |
|---|---|
| Fig. 1B | `05_analyses/gene_set_venn.py` |
| Fig. 1C–D | `05_analyses/four_way_dnds.py` |
| Fig. 2A, S1 | `05_analyses/dnds_decomposition.py` |
| Fig. 2B–C, S4, table S4 | `05_analyses/codon_usage_bias.py` (CAI reference: `make_cai_reference.py`) |
| Fig. 2D, table S1 gene-structure columns | `05_analyses/gene_structure.py` |
| S6 | `05_analyses/genome_trace.py` |
| Fig. 3A, 3C | `05_analyses/cbr_meta_profile.py` |
| Fig. 3B | `05_analyses/gene_density_cbr.py` |
| Fig. 3D, S3B | `05_analyses/cbr_changepoint_relax.py` |
| Fig. 4, S8, table S7 | `03_rates/branch_dnds_four.py`, `05_analyses/branch_dnds_decomposition.py` |
| Fig. 5A–C, S12A | `05_analyses/functional_enrichment.py` |
| Fig. 5D–E, S13B, table S9 | `05_analyses/multivariate_insulation.py` |
| Fig. 5, S12, S13B, table S9 (expression joined to the gene models) | `05_analyses/workbook_model_bridge.py` |
| Fig. 6A, 6D | `07_polymorphism/pi_by_region.py` |
| Fig. 6B–C | `07_polymorphism/pi4_matched_control.py` |
| Fig. 6 (assembly) | `07_polymorphism/figure_polymorphism.py` |
| Results text (neighborhood control) | `07_polymorphism/pis_vs_position.py` |
| S2 | `05_analyses/compositional_analysis.py` |
| S3A, table S3 | `05_analyses/positive_selection.py` |
| S5, table S5 | `05_analyses/microsynteny.py` |
| S7, tables S1 and S7 (region labels) | `05_analyses/annotate_regions.py` |
| S9, table S8 | `05_analyses/toxocara_annotation_diagnostic.py` |
| S10 | `05_analyses/tox_pairwise_decomposition.py` (labels: `04_annotation/annotate_tox_kaks.py`) |
| S11 | `05_analyses/saturation_diagnostic.py` |
| S12B–C | `05_analyses/germline_restricted_dnds.py` |
| S12D | `05_analyses/dnds_within_class.py` |
| S13A | `05_analyses/kaks_compartments.py` |
| S14, table S10 | `07_polymorphism/figure_cohort_selection.py`, `germline_fraction_audit.py`, `select_germline_samples.py`, `depth_match.py`, `allele_balance_panel.py` |
| S15 | `07_polymorphism/sweep_min_s4.py` |
| S16 | `07_polymorphism/reference_bias_check.py` |
| S17 | `07_polymorphism/figure_cohort_integrity.py`, `check_duplicate_individuals.py`, `het_by_region.py`, `jackknife_individuals.py` |
| S18 | `05_analyses/rbh_quality.py` |
| Table S2 | `03_rates/molecular_clock_test.py`, `06_outputs/make_rate_table.py` |
| Table S6 | `05_analyses/x_autosome_check.py`, `hic_x_dosage.py` |
| Table S11 | `05_analyses/polymorphism_alternatives.tsv` (curated rows) and the McDonald–Kreitman summary written by `07_polymorphism/pi_by_region.py` |
| All supplementary tables | `05_analyses/make_supplementary_tables.py` (inputs in `supplementary_table_sources.tsv`) |

## Software

Python 3.10 (divergence) and 3.11 (polymorphism) with NumPy, pandas, SciPy,
matplotlib, Biopython and openpyxl; BLAST+, MAFFT, PAML (yn00, codeml, baseml),
minimap2 and gffread; bwa-mem2, SAMtools, BCFtools, GenMap and mosdepth. Exact
versions are in `envs/`. The gene sets were made with miniprot v0.18-r281 and
TransDecoder v5.5.0, which are not in the environments.

## Citation

Cite the release you used by its Zenodo DOI. `CITATION.cff` holds the citation
details, which GitHub shows under "Cite this repository". The paper is in
preparation.

## License

MIT (see `LICENSE`).
