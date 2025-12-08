# Plotting Scripts Status

## Issue: Matplotlib Dependencies

The plotting scripts (04.5, 06) require matplotlib with full dependencies, which is not currently available in the environment:

```
ModuleNotFoundError: No module named 'pyparsing'
```

## Solutions

### Option 1: Install Dependencies
```bash
pip install matplotlib pyparsing pillow
# or
conda install matplotlib pyparsing pillow
```

### Option 2: Use Existing Plotting Scripts
The original plotting scripts from `scripts/4_statistical_analyses_KEGG_PATHWAY_version/` can be manually adapted:

1. Copy the script
2. Update file paths to use KEGG data directory
3. Add data type suffix to figure names
4. Run for each data type separately

### Option 3: Custom Plotting in Jupyter/R
Use the generated TSV files to create custom plots:

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
global_params = pd.read_csv(
    'results/4.5_statistical_analyses_KEGG_PATHWAY_version/03_global_scaling/global_scaling_params_reactions.tsv',
    sep='\t'
)

# Create plots
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(global_params['alpha_global'], bins=50)
ax.set_xlabel('Alpha')
ax.set_ylabel('Count')
ax.set_title('Distribution of Scaling Exponents - Reactions')
plt.savefig('alpha_dist_reactions.pdf')
```

## Available Data for Plotting

All statistical outputs are available in TSV format:

### Global Scaling Parameters
```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/03_global_scaling/
  - global_scaling_params_reactions.tsv (448 categories)
  - global_scaling_params_ko.tsv (684 categories)
  - global_scaling_params_pathway.tsv (198 categories)
```

Columns: category, alpha_global, alpha_global_se, alpha_global_ci99_low, alpha_global_ci99_high, beta_global_log, beta_global_log_se, beta_global_ci99_low, beta_global_ci99_high, n_genomes_used, r_squared, p_value

### Environment-Specific Parameters
```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/
  - env_scaling_params_reactions.tsv (2,668 env×category fits)
  - env_scaling_params_ko.tsv (4,617 env×category fits)
  - env_scaling_params_pathway.tsv (1,420 env×category fits)
```

Columns: environment, category, alpha_env, alpha_env_se, alpha_env_ci99_low, alpha_env_ci99_high, beta_env_log, beta_env_log_se, n_genomes_used, r_squared, p_value

### Z-Scores
```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/
  - env_vs_global_Z_scores_reactions.tsv (2,016 Z-scores)
  - env_vs_global_Z_scores_ko.tsv (3,568 Z-scores)
  - env_vs_global_Z_scores_pathway.tsv (1,372 Z-scores)
  - category_Z_summary_reactions.tsv (361 categories)
  - category_Z_summary_ko.tsv (574 categories)
  - category_Z_summary_pathway.tsv (190 categories)
```

Z-score columns: environment, category, Z_alpha, Z_beta, n_genomes_used
Summary columns: category, Z_alpha_category, Z_beta_category, n_envs_used

### Master Tables
```
results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts/
  - master_table_env_filtered_reactions.tsv (2164 genomes × 485 reactions)
  - master_table_env_filtered_ko.tsv (2164 genomes × 739 KOs)
  - master_table_env_filtered_pathway.tsv (2164 genomes × 202 pathways)
```

## Recommended Figures to Generate

### Essential QC Figures
1. **Alpha distribution** - Histogram of scaling exponents
2. **Alpha vs R²** - Quality of fits
3. **Z-score distributions** - Environmental variation
4. **Top variable categories** - Bar plot of highest Z-scores
5. **Environment counts** - Sample sizes per environment

### Advanced Figures (Optional)
6. **Scaling plots per category** - Scatter plots with fitted lines
7. **Environment comparison** - Heatmaps of alpha values
8. **Residual plots** - Model diagnostics
9. **Confidence interval plots** - Uncertainty visualization

## Example Plotting Code

### Quick Visualization Script
```python
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def plot_kegg_results(data_type):
    """Generate essential plots for a KEGG data type."""
    BASE = Path("results/4.5_statistical_analyses_KEGG_PATHWAY_version")
    suffix = f"_{data_type}"
    
    # Load data
    global_params = pd.read_csv(BASE / f"03_global_scaling/global_scaling_params{suffix}.tsv", sep='\t')
    z_scores = pd.read_csv(BASE / f"04_env_scaling/category_Z_summary{suffix}.tsv", sep='\t')
    
    # Plot 1: Alpha distribution
    plt.figure(figsize=(10, 6))
    plt.hist(global_params['alpha_global'], bins=50, edgecolor='black')
    plt.xlabel('Alpha (Scaling Exponent)')
    plt.ylabel('Count')
    plt.title(f'Scaling Exponent Distribution - {data_type.upper()}')
    plt.savefig(f'alpha_dist_{data_type}.pdf')
    plt.close()
    
    # Plot 2: Top variable categories
    top_20 = z_scores.nlargest(20, 'Z_alpha_category')
    plt.figure(figsize=(12, 8))
    plt.barh(range(len(top_20)), top_20['Z_alpha_category'])
    plt.yticks(range(len(top_20)), top_20['category'])
    plt.xlabel('Z_alpha')
    plt.title(f'Top 20 Variable Categories - {data_type.upper()}')
    plt.tight_layout()
    plt.savefig(f'top_variable_{data_type}.pdf')
    plt.close()
    
    print(f"✓ Generated plots for {data_type}")

# Generate for all types
for dt in ['reactions', 'ko', 'pathway']:
    plot_kegg_results(dt)
```

## Status Summary

- **Core Pipeline (Scripts 01-04)**: ✅ COMPLETE
- **Label Mapping (Script 05)**: 🔄 IN PROGRESS
- **Plotting Scripts (04.5, 06, 07)**: ⏳ PENDING (requires matplotlib or manual adaptation)
- **Documentation**: ⏳ PENDING

The statistical analysis pipeline is fully functional. Plotting can be done separately using the generated TSV files.

