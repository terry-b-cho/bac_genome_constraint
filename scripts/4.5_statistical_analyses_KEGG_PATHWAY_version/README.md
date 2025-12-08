# KEGG Pathway Statistical Analyses Pipeline

## Overview

This directory contains a complete statistical analysis pipeline for KEGG data (reactions, KOs, and pathways), adapted from the GO term analysis pipeline. The pipeline performs scaling law analyses and computes environment-specific variations.

## Status: ✅ CORE PIPELINE COMPLETE

**All essential statistical analyses (Scripts 01-04) are functional and tested.**

---

## Quick Start

### Running the Complete Pipeline

```bash
cd scripts/4.5_statistical_analyses_KEGG_PATHWAY_version

# Run all core scripts
python 01_build_master_table_env.py
python 02_define_env_cohorts.py
python 03_fit_global_scaling.py
bash 04_fit_env_scaling_and_Z.sh
bash 05_map_kegg_labels.sh  # Optional, for better labels
```

### With Test Mode (Fast)

```bash
python 01_build_master_table_env.py --test-mode
python 02_define_env_cohorts.py --test-mode
python 03_fit_global_scaling.py --test-mode
bash 04_fit_env_scaling_and_Z.sh --test-mode
```

### With Prevalence Filtering

```bash
# Use 95% or 99% prevalence threshold
python 01_build_master_table_env.py --prevalence-threshold 95
python 02_define_env_cohorts.py --prevalence-threshold 95
python 03_fit_global_scaling.py --prevalence-threshold 95
bash 04_fit_env_scaling_and_Z.sh --prevalence-threshold 95
```

---

## Pipeline Scripts

### ✅ Script 01: Build Master Table
**File**: `01_build_master_table_env.py`  
**Purpose**: Integrate KEGG counts with genome metadata and apply QC filters  
**Input**: KEGG count tables from `results/3.5_KEGG_n_reaction_analyses/`  
**Output**: Master tables for all three data types  
**Status**: Complete & Tested

### ✅ Script 02: Define Environment Cohorts
**File**: `02_define_env_cohorts.py`  
**Purpose**: Select environments with ≥20 genomes  
**Output**: Filtered master tables and environment lists  
**Status**: Complete & Tested

### ✅ Script 03: Fit Global Scaling
**File**: `03_fit_global_scaling.py`  
**Purpose**: Fit power-law scaling: nc(g) = β × n(g)^α  
**Output**: Global scaling parameters (α, β, R², SE, CI)  
**Status**: Complete & Tested

**Results**:
- Reactions: 448/485 categories fitted
- KOs: 684/739 categories fitted
- Pathways: 198/202 categories fitted

### ✅ Script 04: Environment-Specific Fits & Z-Scores
**Files**: `04_fit_env_scaling_and_Z_single.py`, `04_fit_env_scaling_and_Z.sh`  
**Purpose**: Fit per-environment scaling and compute Z-scores  
**Output**: Environment-specific parameters and Z-scores  
**Status**: Complete & Tested

**Results**:
- Reactions: 2,668 env×category fits, 361 categories with Z-scores
- KOs: 4,617 env×category fits, 574 categories with Z-scores
- Pathways: 1,420 env×category fits, 190 categories with Z-scores

### 🔄 Script 05: Map KEGG Labels
**Files**: `05_map_kegg_labels_single.py`, `05_map_kegg_labels.sh`  
**Purpose**: Fetch KEGG names and definitions from API  
**Status**: Running in background (uses caching)

### 📋 Scripts 04.5, 06, 07: Plotting & Extensions
**Status**: Templates and minimal versions provided  
**See**: `PLOTTING_SCRIPTS_STATUS.md` and `04.5_plot_intermediate_figures_TEMPLATE.md`

---

## Data Types Processed

The pipeline processes three types of KEGG data:

| Type | Column Pattern | Example | Count |
|------|----------------|---------|-------|
| **Reactions** | R + 5 digits | R00130 | 485 |
| **KOs** | K + 5 digits | K00005 | 739 |
| **Pathways** | map/rn + 5 digits | map00010 | 202 |

All outputs include a suffix indicating the data type:
- `*_reactions.tsv` - Reaction data
- `*_ko.tsv` - KO data
- `*_pathway.tsv` - Pathway data

---

## Input Data

### KEGG Count Tables
**Location**: `results/3.5_KEGG_n_reaction_analyses/`

**Base files**:
- `all_reactions_counts_table_kegg_reaction.txt`
- `all_ko_counts_table_kegg_reaction.txt`
- `all_pathway_counts_table_kegg_reaction.txt`

**Prevalence-filtered files** (optional):
- `*_prev_95.txt` - 95% prevalence threshold
- `*_prev_99.txt` - 99% prevalence threshold

### Genome Metadata
**File**: `results/2_JGIgold_KEGG_anayses_out/05_genome_feature_matrix.tsv`  
**Contains**: Genome accessions, environments, gene counts, QC metrics

---

## Output Structure

```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/
│
├── 01_master_table/
│   ├── master_table_raw_reactions.tsv           (3088 genomes)
│   ├── master_table_raw_ko.tsv                  (3088 genomes)
│   ├── master_table_raw_pathway.tsv             (3088 genomes)
│   ├── master_table_high_quality_reactions.tsv  (2209 genomes)
│   ├── master_table_high_quality_ko.tsv         (2209 genomes)
│   ├── master_table_high_quality_pathway.tsv    (2209 genomes)
│   └── qc_01_master_table_*.log
│
├── 02_env_cohorts/
│   ├── valid_environments_min20_reactions.tsv   (8 environments)
│   ├── valid_environments_min20_ko.tsv          (8 environments)
│   ├── valid_environments_min20_pathway.tsv     (8 environments)
│   ├── master_table_env_filtered_reactions.tsv  (2164 genomes × 485 reactions)
│   ├── master_table_env_filtered_ko.tsv         (2164 genomes × 739 KOs)
│   ├── master_table_env_filtered_pathway.tsv    (2164 genomes × 202 pathways)
│   └── qc_02_env_cohorts_*.log
│
├── 03_global_scaling/
│   ├── global_scaling_params_reactions.tsv      (448 categories)
│   ├── global_scaling_params_ko.tsv             (684 categories)
│   ├── global_scaling_params_pathway.tsv        (198 categories)
│   └── qc_03_global_scaling_*.log
│
├── 04_env_scaling/
│   ├── env_scaling_params_reactions.tsv         (2,668 fits)
│   ├── env_scaling_params_ko.tsv                (4,617 fits)
│   ├── env_scaling_params_pathway.tsv           (1,420 fits)
│   ├── env_vs_global_Z_scores_reactions.tsv     (2,016 Z-scores)
│   ├── env_vs_global_Z_scores_ko.tsv            (3,568 Z-scores)
│   ├── env_vs_global_Z_scores_pathway.tsv       (1,372 Z-scores)
│   ├── category_Z_summary_reactions.tsv         (361 categories)
│   ├── category_Z_summary_ko.tsv                (574 categories)
│   ├── category_Z_summary_pathway.tsv           (190 categories)
│   └── qc_04_env_scaling_*.log
│
└── 05_kegg_labels/
    ├── kegg_term_labels_reactions.tsv           (In progress)
    ├── kegg_term_labels_ko.tsv                  (In progress)
    └── kegg_term_labels_pathway.tsv             (In progress)
```

---

## Example Analyses

### Find Most Environmentally Variable Reactions

```python
import pandas as pd

z_scores = pd.read_csv(
    'results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/category_Z_summary_reactions.tsv',
    sep='\t'
)

# Top 10 most variable
top_10 = z_scores.nlargest(10, 'Z_alpha_category')
print("Top 10 Most Environmentally Variable Reactions:")
print(top_10[['category', 'Z_alpha_category', 'Z_beta_category', 'n_envs_used']])
```

### Compare Scaling Exponents Across Data Types

```python
import pandas as pd

reactions = pd.read_csv('...global_scaling_params_reactions.tsv', sep='\t')
kos = pd.read_csv('...global_scaling_params_ko.tsv', sep='\t')
pathways = pd.read_csv('...global_scaling_params_pathway.tsv', sep='\t')

print(f"Reactions - Mean alpha: {reactions['alpha_global'].mean():.4f}")
print(f"KOs - Mean alpha: {kos['alpha_global'].mean():.4f}")
print(f"Pathways - Mean alpha: {pathways['alpha_global'].mean():.4f}")
```

### Analyze Specific Environment

```python
import pandas as pd

env_params = pd.read_csv(
    'results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/env_scaling_params_reactions.tsv',
    sep='\t'
)

# Get all fits for Aquatic environment
aquatic = env_params[env_params['environment'] == 'Aquatic']
print(f"Aquatic environment: {len(aquatic)} category fits")
print(f"Mean alpha: {aquatic['alpha_env'].mean():.4f}")
print(f"Mean R²: {aquatic['r_squared'].mean():.4f}")
```

---

## Documentation Files

- **README.md** (this file) - Overview and quick start
- **FINAL_STATUS_REPORT.md** - Comprehensive status and results
- **IMPLEMENTATION_PROGRESS.md** - Detailed implementation notes
- **REMAINING_WORK_SUMMARY.md** - What's left to do
- **PLOTTING_SCRIPTS_STATUS.md** - Plotting script information
- **04.5_plot_intermediate_figures_TEMPLATE.md** - Template for adapting plotting scripts

---

## Technical Notes

### Column Detection
The pipeline automatically detects KEGG data types:
- Reactions: Columns matching `^R\d{5}$`
- KOs: Columns matching `^K\d{5}$`
- Pathways: Columns matching `^(map|rn)\d{5}$`

### Prevalence Filtering
Two approaches supported:
1. **Pre-filtered files**: Uses `*_prev_95.txt` or `*_prev_99.txt` if available
2. **On-the-fly**: Filters during processing if pre-filtered files not found

### Wrapper Scripts
For scripts processing all three types, bash wrappers are provided:
- `04_fit_env_scaling_and_Z.sh` - Runs Script 04 for all types
- `05_map_kegg_labels.sh` - Runs Script 05 for all types
- `04.5_plot_minimal.sh` - Runs minimal plotting for all types

---

## Troubleshooting

### Issue: Script doesn't produce output
- Check that previous scripts have been run
- Verify input files exist in expected locations
- Check QC log files for errors

### Issue: Matplotlib not available
- Install: `pip install matplotlib pyparsing pillow`
- Or use custom plotting with available tools

### Issue: KEGG API timeout
- Script 05 uses caching to avoid repeated API calls
- Cached data stored in `results/3.5_KEGG_n_reaction_analyses/kegg_api_cache/`
- Can re-run if interrupted

---

## Citation & Methods

When using this pipeline, cite the original methods and note the adaptation to KEGG data:

> "We adapted the genome-wide scaling analysis pipeline to analyze KEGG reactions, ortholog groups (KOs), and metabolic pathways. For each KEGG category, we fitted power-law scaling relationships nc(g) = β × n(g)^α across all genomes, where nc(g) is the count of genes in category c and n(g) is the total gene count. We then computed environment-specific scaling parameters and Z-scores to identify categories with strong environmental variation."

---

## Support

For issues or questions:
1. Check the documentation files listed above
2. Review QC log files in output directories
3. Verify input data files are present and correctly formatted

---

**Pipeline Version**: 4.5 (KEGG Pathway Version)  
**Based On**: 4.0 (GO Term Version)  
**Last Updated**: December 7, 2025

