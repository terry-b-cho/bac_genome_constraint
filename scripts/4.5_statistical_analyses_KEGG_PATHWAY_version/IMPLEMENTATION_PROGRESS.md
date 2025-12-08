# Implementation Progress: KEGG Pathway Statistical Analyses Pipeline

## Status: PARTIALLY COMPLETE (Core Pipeline Functional)

Date: December 7, 2025

## Completed Components

### ✅ 1. Directory Structure
- Created: `scripts/4.5_statistical_analyses_KEGG_PATHWAY_version/`
- Created: `results/4.5_statistical_analyses_KEGG_PATHWAY_version/`
- All files copied from original GO version

### ✅ 2. Utility Functions (`prevalence_utils.py`)
**Status: COMPLETE**

New functions added:
- `detect_kegg_data_type(df)` - Auto-detects reactions/KOs/pathways from column names
- `get_kegg_columns(df, data_type)` - Extracts KEGG columns by type
- `filter_kegg_columns_by_prevalence()` - Filters KEGG terms by prevalence
- `get_kegg_input_files()` - Returns correct input file paths with prevalence support
- `get_data_type_suffix()` - Returns suffix for output files (_reactions, _ko, _pathway)

### ✅ 3. Script 01: Build Master Table
**Status: COMPLETE & TESTED**

- ✅ Processes all three data types (reactions, KOs, pathways)
- ✅ Generates outputs with appropriate suffixes
- ✅ Handles prevalence-filtered input files
- ✅ Tested in test mode and full mode
- ✅ All outputs verified

**Outputs Created:**
```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/01_master_table/
  - master_table_raw_reactions.tsv (3088 rows)
  - master_table_raw_ko.tsv (3088 rows)
  - master_table_raw_pathway.tsv (3088 rows)
  - master_table_high_quality_reactions.tsv (2209 rows)
  - master_table_high_quality_ko.tsv (2209 rows)
  - master_table_high_quality_pathway.tsv (2209 rows)
  - environment_counts_all_*.tsv (11 environments each)
  - qc_01_master_table_*.log
```

### ✅ 4. Script 02: Define Environment Cohorts
**Status: COMPLETE & TESTED**

- ✅ Processes all three data types
- ✅ Filters to environments with ≥20 genomes
- ✅ Generates outputs with suffixes
- ✅ Tested in test mode and full mode
- ✅ All outputs verified

**Outputs Created:**
```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts/
  - valid_environments_min20_reactions.tsv (8 environments)
  - valid_environments_min20_ko.tsv (8 environments)
  - valid_environments_min20_pathway.tsv (8 environments)
  - master_table_env_filtered_reactions.tsv/parquet (2164 genomes, 485 reactions)
  - master_table_env_filtered_ko.tsv/parquet (2164 genomes, 739 KOs)
  - master_table_env_filtered_pathway.tsv/parquet (2164 genomes, 202 pathways)
  - qc_02_env_cohorts_*.log
```

### ✅ 5. Script 03: Fit Global Scaling
**Status: COMPLETE & TESTED**

- ✅ Processes all three data types
- ✅ Fits power-law scaling: nc(g) = β × n(g)^α
- ✅ Generates outputs with suffixes
- ✅ Tested in test mode and full mode
- ✅ All outputs verified

**Outputs Created:**
```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/03_global_scaling/
  - global_scaling_params_reactions.tsv (448 categories)
  - global_scaling_params_ko.tsv (684 categories)
  - global_scaling_params_pathway.tsv (198 categories)
  - global_scaling_params_*.parquet
  - qc_03_global_scaling_*.log
```

**Key Results:**
- Reactions: 448/485 fitted (37 skipped due to low N)
- KOs: 684/739 fitted (55 skipped due to low N)
- Pathways: 198/202 fitted (4 skipped due to low N)

### ✅ 6. Script 04: Environment-Specific Fits & Z-Scores
**Status: COMPLETE & TESTED**

- ✅ Processes all three data types
- ✅ Fits per-environment scaling laws
- ✅ Computes Z-scores comparing env-specific to global
- ✅ Created wrapper script (04_fit_env_scaling_and_Z.sh) and single script
- ✅ Tested in test mode and full mode
- ✅ All outputs verified

**Outputs Created:**
```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/
  - env_scaling_params_reactions.tsv (2,668 env×category fits)
  - env_scaling_params_ko.tsv (4,617 env×category fits)
  - env_scaling_params_pathway.tsv (1,420 env×category fits)
  - env_vs_global_Z_scores_reactions.tsv (2,016 Z-scores)
  - env_vs_global_Z_scores_ko.tsv (3,568 Z-scores)
  - env_vs_global_Z_scores_pathway.tsv (1,372 Z-scores)
  - category_Z_summary_reactions.tsv (361 categories)
  - category_Z_summary_ko.tsv (574 categories)
  - category_Z_summary_pathway.tsv (190 categories)
  - qc_04_env_scaling_*.log
```

### 🔄 7. Script 05: Map KEGG Labels
**Status: IN PROGRESS (Running in background)**

- ✅ Created single script version (05_map_kegg_labels_single.py)
- ✅ Created wrapper script (05_map_kegg_labels.sh)
- ✅ Tested in test mode - works correctly
- 🔄 Currently running on full dataset (fetching from KEGG API with caching)
- ⏳ Expected completion time: ~10-30 minutes (depends on cache availability)

**Implementation:**
- Uses KEGG REST API to fetch term names and definitions
- Caches results in `results/3.5_KEGG_n_reaction_analyses/kegg_api_cache/`
- Processes: reactions (485 terms), KOs (739 terms), pathways (202 terms)

## Remaining Components

### ⏳ 8. Script 04.5: Plot Intermediate Figures
**Status: NOT STARTED**
**Complexity: HIGH (1616 lines)**

This script generates QC and intermediate diagnostic figures. Needs:
- Update all file paths to use KEGG data
- Add data type suffix to all figure filenames
- Update column identification logic
- Create wrapper to process all three types

**Estimated Time: 2-3 hours**

### ⏳ 9. Script 06: Make Scaling Figures
**Status: NOT STARTED**
**Complexity: HIGH (994 lines)**

This script generates publication-quality figures. Needs:
- Update all file paths to use KEGG data
- Add data type suffix to all figure filenames
- Update labels and titles for KEGG terminology
- Create wrapper to process all three types

**Estimated Time: 2-3 hours**

### ⏳ 10. Script 07: Scaling Extensions
**Status: NOT STARTED**
**Complexity: MEDIUM (431 lines)**

This script performs extended analyses (TF, mobile elements, nutrients). Needs:
- Update file paths
- Add data type suffixes
- Create wrapper

**Estimated Time: 1-2 hours**

### ⏳ 11. Documentation Updates
**Status: NOT STARTED**

Files to update:
- `CODE_WORKFLOW_AND_FIGURES_SUMMARY.md`
- `00_Figure_descriptions_n_methods.md`
- `descriptive_plan_and_required_sanity_check.md`
- `reproducible_data_description.md`

## How to Run Completed Scripts

### Full Pipeline (Scripts 01-04)
```bash
cd scripts/4.5_statistical_analyses_KEGG_PATHWAY_version

# Script 01: Build master tables
python 01_build_master_table_env.py

# Script 02: Define environment cohorts
python 02_define_env_cohorts.py

# Script 03: Fit global scaling
python 03_fit_global_scaling.py

# Script 04: Fit environment-specific scaling & Z-scores
bash 04_fit_env_scaling_and_Z.sh

# Script 05: Map KEGG labels (running in background)
bash 05_map_kegg_labels.sh
```

### With Prevalence Filtering
```bash
# Use --prevalence-threshold 95 or 99
python 01_build_master_table_env.py --prevalence-threshold 95
python 02_define_env_cohorts.py --prevalence-threshold 95
python 03_fit_global_scaling.py --prevalence-threshold 95
bash 04_fit_env_scaling_and_Z.sh --prevalence-threshold 95
bash 05_map_kegg_labels.sh --prevalence-threshold 95
```

### Test Mode
```bash
# Add --test-mode to any script
python 01_build_master_table_env.py --test-mode
```

## Key Implementation Details

### Data Type Detection
The pipeline automatically detects KEGG data types based on column patterns:
- **Reactions**: Columns starting with "R" + 5 digits (e.g., R00130)
- **KOs**: Columns starting with "K" + 5 digits (e.g., K00005)
- **Pathways**: Columns starting with "map" or "rn" (e.g., map00010, rn00010)

### File Naming Convention
All outputs include data type suffix:
- `*_reactions.tsv` - For reaction data
- `*_ko.tsv` - For KO data
- `*_pathway.tsv` - For pathway data

### Prevalence Filtering
The pipeline supports two approaches:
1. **Pre-filtered files**: Uses `*_prev_95.txt` or `*_prev_99.txt` if available
2. **On-the-fly filtering**: Filters during processing if pre-filtered files not found

## Next Steps

1. **Wait for Script 05 to complete** (~10-30 minutes)
2. **Adapt Script 04.5** (plotting intermediate figures)
3. **Adapt Script 06** (publication figures)
4. **Adapt Script 07** (extended analyses)
5. **Update documentation**

## Notes

- All core statistical analyses (Scripts 01-04) are functional and tested
- Scripts 06-07 (plotting) require significant adaptation due to size and complexity
- Consider whether all intermediate figures are needed or if a subset would suffice
- The plotting scripts may benefit from refactoring to reduce duplication

