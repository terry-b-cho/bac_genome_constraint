# Remaining Work Summary

## Current Status

### ✅ COMPLETED (Scripts 01-04)
The core statistical analysis pipeline is **fully functional** for all three KEGG data types:

1. **Script 01**: Build Master Table - ✅ COMPLETE
2. **Script 02**: Define Environment Cohorts - ✅ COMPLETE
3. **Script 03**: Fit Global Scaling - ✅ COMPLETE
4. **Script 04**: Environment-Specific Fits & Z-Scores - ✅ COMPLETE

### 🔄 IN PROGRESS
5. **Script 05**: Map KEGG Labels - Currently running in background (fetching from KEGG API)

### ⏳ REMAINING (Scripts 04.5, 06, 07)
These are large plotting scripts that require adaptation:

## Detailed Breakdown of Remaining Work

### Script 04.5: Plot Intermediate Figures (1616 lines)

**Purpose**: Generate QC and diagnostic plots for Scripts 01-04

**Required Changes:**
1. Update all input file paths to use KEGG data directory
2. Replace GO column detection with KEGG column detection
3. Add data type suffix to all figure filenames
4. Create wrapper to process all three types
5. Update plot titles and labels for KEGG terminology

**Key Sections to Modify:**
- Line 100-418: Script 01 QC plots (environment distributions, genome stats)
- Line 419-550: Script 02 QC plots (environment filtering)
- Line 551-952: Script 03 QC plots (global scaling fits, residuals)
- Line 953-1596: Script 04 QC plots (environment-specific fits, Z-scores)

**Approach:**
- Create single-type version: `04.5_plot_intermediate_figures_single.py`
- Create bash wrapper: `04.5_plot_intermediate_figures.sh`
- Each section needs minimal changes (mostly path updates and suffix additions)

**Estimated Time**: 2-3 hours

### Script 06: Make Scaling Figures (994 lines)

**Purpose**: Generate publication-quality figures

**Required Changes:**
1. Update all input file paths
2. Replace GO terminology with KEGG terminology in titles/labels
3. Add data type suffix to figure filenames
4. Create wrapper for all three types
5. Update legend labels and annotations

**Key Sections:**
- Multiple figure generation functions
- Complex multi-panel plots
- Statistical annotations

**Approach:**
- Create single-type version: `06_make_scaling_figures_single.py`
- Create bash wrapper: `06_make_scaling_figures.sh`
- Focus on core figures first, optional figures later

**Estimated Time**: 2-3 hours

### Script 07: Scaling Extensions (431 lines)

**Purpose**: Extended analyses for transcription factors, mobile elements, nutrients

**Required Changes:**
1. Update input file paths
2. Add data type suffixes
3. Create wrapper
4. May need to adapt TF/mobile element detection for KEGG context

**Approach:**
- Create single-type version: `07_scaling_extensions_single.py`
- Create bash wrapper: `07_scaling_extensions.sh`
- This script may need conceptual adaptation depending on analysis goals

**Estimated Time**: 1-2 hours

## Recommended Approach

### Option 1: Complete All Scripts (Comprehensive)
**Time Required**: 5-8 hours
- Adapt all plotting scripts
- Generate all figures for all data types
- Full documentation updates

### Option 2: Core Pipeline Only (Pragmatic)
**Time Required**: Current state + 1 hour
- Keep Scripts 01-04 as primary deliverable
- Create minimal plotting scripts for essential figures only
- Document that advanced plotting scripts can be adapted as needed

### Option 3: Incremental Completion
**Time Required**: Flexible
- Complete Script 05 (in progress)
- Adapt Script 04.5 for basic QC plots
- Leave Scripts 06-07 for future work with clear documentation

## Current Pipeline Capabilities

Even without the plotting scripts, the current pipeline provides:

### ✅ Statistical Outputs Available
- Master tables with QC filtering
- Environment cohort definitions
- Global scaling parameters (α, β, R², SE, CI)
- Environment-specific scaling parameters
- Z-scores for environment vs global comparisons
- Category-level Z-score summaries

### ✅ Data Types Processed
- Reactions: 485 terms → 448 with global fits → 361 with Z-scores
- KOs: 739 terms → 684 with global fits → 574 with Z-scores
- Pathways: 202 terms → 198 with global fits → 190 with Z-scores

### ✅ Environments Analyzed
8 environments with ≥20 genomes:
- Aquatic (598 genomes)
- Terrestrial (479 genomes)
- Mammals: Human (406 genomes)
- Plants (281 genomes)
- Mammals (198 genomes)
- Food production (103 genomes)
- Wastewater (63 genomes)
- Birds (36 genomes)

## How to Use Current Outputs

### Example: Analyzing Reaction Scaling

```python
import pandas as pd

# Load global scaling parameters
global_params = pd.read_csv(
    'results/4.5_statistical_analyses_KEGG_PATHWAY_version/03_global_scaling/global_scaling_params_reactions.tsv',
    sep='\t'
)

# Load Z-scores
z_scores = pd.read_csv(
    'results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/category_Z_summary_reactions.tsv',
    sep='\t'
)

# Find reactions with highest environmental variation
top_variable = z_scores.nlargest(10, 'Z_alpha_category')
print(top_variable[['category', 'Z_alpha_category', 'n_envs_used']])
```

### Custom Plotting

Users can create custom plots using the output data:
```python
import matplotlib.pyplot as plt
import numpy as np

# Plot global scaling for a specific reaction
reaction_id = 'R00130'
reaction_params = global_params[global_params['category'] == reaction_id].iloc[0]

# Load master table to plot raw data
master = pd.read_csv(
    'results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts/master_table_env_filtered_reactions.tsv',
    sep='\t'
)

# Create scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(np.log(master['genes_total']), np.log(master[reaction_id]), alpha=0.3)
plt.xlabel('log(Total Genes)')
plt.ylabel(f'log({reaction_id} Count)')
plt.title(f'Global Scaling: {reaction_id}')
plt.savefig(f'custom_plot_{reaction_id}.pdf')
```

## Conclusion

**The core statistical analysis pipeline (Scripts 01-04) is complete and functional.**

The remaining work involves:
1. Completing Script 05 (in progress, ~30 min)
2. Adapting plotting scripts (Scripts 04.5, 06, 07) - optional but recommended
3. Documentation updates

All essential statistical outputs are already generated and ready for analysis.



