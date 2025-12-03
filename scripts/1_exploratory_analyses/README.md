# Exploratory Analyses

This directory contains notebooks and scripts for preliminary data exploration and analysis.

## Environment Setup

Before running notebooks, set up the conda environment:

```bash
# Load required O2 modules
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate genome_constraint_envs_O2_py

# Verify installation
python -c "import pandas, numpy, matplotlib, seaborn; print('All packages available')"
```

## Scripts and Notebooks

### `01_preliminary_analysis.py` / `01_preliminary_analysis.ipynb`

Comprehensive preliminary analysis of bacterial genome size constraints, designed to support downstream modeling of how environment constrains genome size.

**Analysis Components:**

1. **Genome Size Distribution Analysis**
   - Histograms (linear and log scale)
   - Box plots and violin plots
   - Statistical summaries (mean, median, quartiles)

2. **Environment/Ecology Variation**
   - Genome size by ecosystem category (GOLD data)
   - Host-associated vs. free-living comparisons
   - Isolation source analysis
   - Statistical tests (Mann-Whitney U)

3. **Correlation Analyses**
   - Genome size vs. gene count (total and protein-coding)
   - Genome size vs. GC content
   - Scatter plots with correlation coefficients

4. **Environment Selection and Taxonomic Diversity**
   - Identification of environments with sufficient genomes (>=20)
   - Taxonomic diversity assessment per environment
   - Extreme genome size environments (smallest and largest)
   - Visualization of selected environments

5. **Genome Size Clustering Analysis**
   - Within-environment clustering patterns
   - Histogram analysis for large environments (>=50 genomes)
   - Identification of multiple survival strategies within environments

6. **Data Quality Verification**
   - KEGG analysis readiness (protein sequences, gene annotations)
   - Nutrient limitation analysis readiness (CDS sequences)
   - Data completeness assessment

7. **Downstream Analysis Preparation**
   - Summary tables for environment analysis
   - Environment selection guide with statistics
   - Ready-to-use datasets for regression/ML modeling

**Required packages:**
- pandas, numpy, matplotlib, seaborn, scipy
- All included in `genome_constraint_envs_O2_py` environment

**Data sources:**
- NCBI metadata: `data/ncbi/metadata/assembly_data_report_extracted.tsv`
- GOLD metadata: `data/gold/0_20251106_gold_metadata_*.csv`
- Dataset catalog: `data/ncbi/metadata/dataset_catalog_extracted.tsv`

**Output files:**
- `results/genome_size_distribution.png` - Distribution histograms
- `results/genome_size_boxplot.png` - Box and violin plots
- `results/genome_size_by_environment.png` - Environment comparisons
- `results/genome_size_vs_gene_count.png` - Correlation plots
- `results/genome_size_vs_gc_content.png` - GC content correlation
- `results/genome_size_extreme_environments.png` - Selected environments
- `results/genome_size_clustering_by_environment.png` - Clustering patterns
- `results/genome_environment_summary.tsv` - Summary table for downstream analysis
- `results/environment_selection_guide.tsv` - Environment statistics and recommendations

**Usage:**
```bash
# Activate environment (see above)
python 01_preliminary_analysis.py
```

**Key Findings:**
- 5,986 genomes analyzed (100% data completeness)
- 5,770 genomes successfully merged with GOLD environment data
- Strong correlations: Genome size vs genes (r=0.989-0.990)
- 16 environments identified with >=20 genomes suitable for detailed analysis
- Extreme environments: Insects (2.44 Mb) to Terrestrial (5.73 Mb)

## Running on O2

If running Jupyter on O2 compute nodes:

```bash
# Start Jupyter server
jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser

# Or use JupyterLab
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser
```

Then connect via SSH tunnel from your local machine:
```bash
ssh -L 8888:localhost:8888 your_username@o2.hms.harvard.edu
```

