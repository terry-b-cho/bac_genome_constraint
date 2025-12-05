# Constraints on Bacterial Genome Size

HST.508 Final Project (Fall 2025)

**Authors:** Aidan Pavao, Terry Cho

## Project Overview

Bacterial genomes vary nearly an order of magnitude in size. Interestingly, bacteria from similar environments across diverse taxa tend to have similar genome sizes (e.g., human skin commensals ~2.5 Mb). This project investigates how environmental factors constrain bacterial genome size through statistical analysis of gene content scaling laws, stratified by environment.

**Key Questions:**
- How does gene content scale with genome size across different functional categories?
- Do scaling relationships differ between environments?
- What functional categories show environment-specific scaling patterns?

## Data

**Source:** NCBI Reference Genome Dataset  
**Total Genomes:** 5,986 complete reference genomes (downloaded and organized)  
**Location:** `data/ncbi/`

**Selection Criteria:**
- Reference genomes only
- Complete assembly level
- Released since 2010
- Annotated (RefSeq and GenBank)
- Excluded atypical, MAGs, and multi-isolate assemblies

**For detailed data documentation:** See `scripts/0_data_download/README.md`

## Repository Structure

```
bac_genome_constraint/
├── data/                              # Input datasets
│   └── ncbi/                          # NCBI genome assemblies (126 GB)
│       ├── assemblies/                # 5,986 genome directories
│       └── metadata/                  # Summary files and manifests
│
├── scripts/                           # Analysis pipelines
│   ├── 0_data_download/               # NCBI data download and organization
│   ├── 1_exploratory_analyses/        # Preliminary data exploration
│   ├── 2_JGIgold_KEGG_anayses/       # KEGG pathway annotation
│   ├── 3_GO_analyses/                 # Gene Ontology (GO) term extraction
│   └── 4_statistical_analyses/        # Main statistical pipeline
│       ├── 01_build_master_table_env.py      # Build master table + QC
│       ├── 02_define_env_cohorts.py           # Environment filtering
│       ├── 03_fit_global_scaling.py           # Global scaling laws
│       ├── 04_fit_env_scaling_and_Z.py       # Environment-specific scaling + Z-scores
│       ├── 04.5_plot_intermediate_figures.py # Intermediate diagnostics
│       ├── 05_map_go_labels.py                # GO term annotation
│       ├── 06_make_scaling_figures.py         # Publication figures
│       ├── 07_scaling_extensions_tf_mobile_nutrient.py  # TF/mobile elements
│       ├── extract_metabolic_go_terms.py       # Metabolism-focused analysis
│       ├── create_prevalence_filtered_terms.py # GO term prevalence filtering
│       └── prevalence_utils.py                 # Shared prevalence utilities
│
└── results/                           # Output files
    ├── 1_exploratory_analyses_out/    # Preliminary analysis outputs
    ├── 3_GO_analyses/                 # GO term counts and tables
    └── 4_statistical_analyses/        # Statistical analysis outputs
        ├── 01_master_table/           # Master data tables
        ├── 02_env_cohorts/             # Environment-filtered data
        ├── 03_global_scaling/          # Global scaling parameters
        ├── 04_env_scaling/             # Environment-specific parameters + Z-scores
        ├── 04.5_intermediate_figures/  # Diagnostic figures
        ├── 05_go_labels/               # GO term labels and metabolic terms
        ├── 06_figures/                 # Publication-quality figures
        └── 07_extensions/              # TF/mobile element analyses
```

## Analysis Pipeline

### Main Statistical Pipeline (`scripts/4_statistical_analyses/`)

The core analysis follows a sequential pipeline:

1. **Script 01:** Build master table with QC filtering (CheckM completeness/contamination, gene counts)
2. **Script 02:** Define environment cohorts (GOLD metadata) and filter to environments with ≥20 genomes
3. **Script 03:** Fit global power-law scaling laws: `n_c(g) = β × n(g)^α` for each GO category
4. **Script 04:** Fit environment-specific scaling laws and compute Z-scores comparing env vs. global parameters
5. **Script 05:** Map GO term IDs to human-readable labels
6. **Script 06:** Generate publication-quality figures (scaling plots, Z-statistics, environment comparisons)
7. **Script 07:** Extend analyses to transcription factors and mobile elements

**Optional Features:**
- **Prevalence filtering:** Filter GO terms by prevalence threshold (e.g., `--prevalence-threshold 95` for terms present in ≥95% of genomes)
- **Metabolic focus:** Generate metabolism-specific plots using `extract_metabolic_go_terms.py`

**For detailed methods and figure descriptions:** See `scripts/4_statistical_analyses/00_Figure_descriptions_n_methods.md`

### Key Outputs

- **Scaling parameters:** Exponents (α) and offsets (β) for each GO category, globally and per-environment
- **Z-statistics:** Quantify environment-specific deviations from global scaling
- **Figures:** Publication-quality plots showing scaling relationships, Z-scores, and environment stratification

## Usage

### Running the Full Pipeline

```bash
# Basic run (no prevalence filtering)
python3 scripts/4_statistical_analyses/01_build_master_table_env.py
python3 scripts/4_statistical_analyses/02_define_env_cohorts.py
python3 scripts/4_statistical_analyses/03_fit_global_scaling.py
python3 scripts/4_statistical_analyses/04_fit_env_scaling_and_Z.py
python3 scripts/4_statistical_analyses/05_map_go_labels.py
python3 scripts/4_statistical_analyses/06_make_scaling_figures.py

# With prevalence filtering (e.g., 95% threshold)
python3 scripts/4_statistical_analyses/01_build_master_table_env.py --prevalence-threshold 95
python3 scripts/4_statistical_analyses/02_define_env_cohorts.py --prevalence-threshold 95
# ... (continue with --prevalence-threshold 95 for all scripts)
```

### Test Mode

Most scripts support `--test-mode` for quick testing on subsets:
```bash
python3 scripts/4_statistical_analyses/03_fit_global_scaling.py --test-mode
```

## Dependencies

- Python 3.x
- pandas, numpy
- matplotlib, seaborn (for plotting)
- scipy (optional, for statistical functions; falls back to manual implementation)
- pyarrow (for Parquet file I/O)

## Contributions

- **Aidan Pavao:** Primary data analysis, visualization, modeling
- **Terry Cho:** Data setup, repository organization, documentation

## Timeline

- **Data criteria due:** Nov 7, 2025
- **Preliminary analysis due:** Nov 17, 2025
- **Presentation:** Dec 8–10, 2025
- **Final paper due:** Dec 12, 2025
