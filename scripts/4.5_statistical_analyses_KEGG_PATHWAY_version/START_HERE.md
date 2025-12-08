# START HERE: KEGG Pipeline Quick Reference

## ✅ What's Been Completed

The **core statistical analysis pipeline** for KEGG data is **fully functional**.

### Scripts Ready to Use

1. ✅ **Script 01**: Build Master Table - `python 01_build_master_table_env.py`
2. ✅ **Script 02**: Define Environment Cohorts - `python 02_define_env_cohorts.py`
3. ✅ **Script 03**: Fit Global Scaling - `python 03_fit_global_scaling.py`
4. ✅ **Script 04**: Environment-Specific Fits - `bash 04_fit_env_scaling_and_Z.sh`
5. 🔄 **Script 05**: Map KEGG Labels - `bash 05_map_kegg_labels.sh` (running)

### Results Generated

**81+ files** created in `results/4.5_statistical_analyses_KEGG_PATHWAY_version/`

- Master tables for 2,209 high-quality genomes
- Global scaling parameters for 1,330 KEGG categories
- Environment-specific parameters for 8 environments
- Z-scores identifying environmentally variable categories

---

## 🚀 Quick Start

### Run the Pipeline

```bash
cd scripts/4.5_statistical_analyses_KEGG_PATHWAY_version

# Run all core scripts (takes ~25 minutes)
python 01_build_master_table_env.py
python 02_define_env_cohorts.py
python 03_fit_global_scaling.py
bash 04_fit_env_scaling_and_Z.sh
```

### Test First (Recommended)

```bash
# Quick test run (~30 seconds)
python 01_build_master_table_env.py --test-mode
python 02_define_env_cohorts.py --test-mode
python 03_fit_global_scaling.py --test-mode
bash 04_fit_env_scaling_and_Z.sh --test-mode
```

---

## 📊 What You Get

### For Each Data Type (reactions, ko, pathway)

1. **Master Tables**
   - Raw data (3,088 genomes)
   - High-quality filtered (2,209 genomes)
   - Environment-filtered (2,164 genomes)

2. **Global Scaling Parameters**
   - Scaling exponent (α)
   - Intercept (β)
   - Standard errors and confidence intervals
   - R² and p-values

3. **Environment-Specific Parameters**
   - Per-environment scaling laws
   - Z-scores vs global
   - Category-level summaries

---

## 📖 Documentation

Read these files for details:

1. **README.md** - Overview and usage
2. **FINAL_STATUS_REPORT.md** - Complete status and results
3. **IMPLEMENTATION_COMPLETE_SUMMARY.md** - Technical details
4. **PLOTTING_SCRIPTS_STATUS.md** - Info about plotting scripts

---

## ⚠️ Important Notes

### Data Types Processed
- **Reactions**: 485 terms → 448 fitted → 361 with Z-scores
- **KOs**: 739 terms → 684 fitted → 574 with Z-scores
- **Pathways**: 202 terms → 198 fitted → 190 with Z-scores

### File Naming
All outputs include data type suffix:
- `*_reactions.tsv`
- `*_ko.tsv`
- `*_pathway.tsv`

### Plotting Scripts
Scripts 04.5, 06, 07 (plotting) are not yet adapted but:
- Templates provided for adaptation
- Minimal plotting script available
- All data available in TSV format for custom plotting

---

## 🔍 Quick Analysis Examples

### Find Top Variable Reactions

```bash
cd results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling
head -20 category_Z_summary_reactions.tsv | column -t
```

### View Global Scaling Parameters

```bash
cd results/4.5_statistical_analyses_KEGG_PATHWAY_version/03_global_scaling
head -20 global_scaling_params_reactions.tsv | column -t
```

### Check Environment Counts

```bash
cd results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts
cat valid_environments_min20_reactions.tsv
```

---

## 💡 Tips

1. **Always run test mode first** when trying new parameters
2. **Check QC log files** if something seems wrong
3. **Use prevalence filtering** to focus on common terms
4. **Scripts are independent per data type** - can run separately if needed

---

## 🎯 What to Do Next

### Option 1: Use Current Outputs
- Analyze the TSV files directly
- Create custom plots in Python/R
- No additional work needed

### Option 2: Complete Plotting Scripts
- Install matplotlib: `pip install matplotlib pyparsing pillow`
- Adapt Scripts 04.5, 06, 07 using templates
- Estimated time: 5-8 hours

### Option 3: Run with Different Parameters
- Try prevalence filtering: `--prevalence-threshold 95`
- Test on specific environments
- Customize QC thresholds

---

**🎉 Congratulations! The core KEGG statistical analysis pipeline is ready to use.**

See **README.md** for detailed usage instructions.
