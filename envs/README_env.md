# Environment Setup for HMS O2 Cluster

This directory contains conda environment specifications for the genome constraint analysis project on HMS O2.

## Overview

### Environment Files

1. **`genome_constraint_envs_O2_py.yaml`** - Original Python environment specification
2. **`genome_constraint_envs_O2_py_success.yaml`** - Python environment with successfully installed packages (use this for reference)
3. **`genome_constraint_envs_O2_R.yaml`** - Original R environment specification  
4. **`genome_constraint_envs_O2_R_success.yaml`** - R environment with successfully installed packages (use this for reference)
5. **`o2_avail_tools.txt`** - List of tools/modules available on O2 that should be loaded via `module load` rather than installed via conda

## Quick Start: Loading the Environment

### For Python Environment

```bash
# Step 1: Load required O2 modules
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0
module load cuda/12.8  # Optional, only needed for GPU/ML workloads

# Step 2: Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Step 3: Activate the Python environment
conda activate genome_constraint_envs_O2_py

# Step 4: Verify installation
python -c "import pandas, numpy, matplotlib, seaborn; print('Python environment ready!')"
```

### For R Environment

```bash
# Step 1: Load required O2 modules
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0
module load cuda/12.8  # Optional, only needed for GPU/ML workloads

# Step 2: Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Step 3: Activate the R environment
conda activate genome_constraint_envs_O2_R

# Step 4: Verify installation
Rscript -e "library(tidyverse); library(ggplot2); cat('R environment ready!\n')"
```

## Using Helper Scripts

For convenience, use the provided shell scripts:

```bash
# Activate Python environment
source envs/activate_py_env.sh

# Activate R environment
source envs/activate_R_env.sh
```

## O2 Module Tools (Load via `module load`)

These tools are available on O2 and should be loaded as modules rather than installed via conda:

```bash
# Bioinformatics tools
module load blast/2.16.0+
module load samtools/1.21
module load bedtools/2.31.0
module load bcftools/1.21
module load htslib/1.21
module load bwa/0.7.18
module load bowtie2/2.5.4
module load bowtie/1.3.1
module load fastqc/0.12.1
module load trimmomatic/0.39
module load sratoolkit/3.2.0
module load hisat2/2.2.1
module load star/2.7.11b
module load gatk/4.6.1.0
module load vcftools/0.1.16
module load snpEff/5.2f
module load homer/5.1
module load plink/1.90b7.7_20241022
module load ucsc-tools/475

# Core languages (if not using conda versions)
module load python/3.13.1
module load R/4.4.2
```

## Creating the Environment (First Time Setup)

If the environment doesn't exist yet:

### Python Environment

```bash
# Load O2 modules
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0
module load cuda/12.8

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Create environment from YAML file
conda env create -f envs/genome_constraint_envs_O2_py_success.yaml

# Activate the environment
conda activate genome_constraint_envs_O2_py
```

### R Environment

```bash
# Load O2 modules
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0
module load cuda/12.8

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Create environment from YAML file
conda env create -f envs/genome_constraint_envs_O2_R_success.yaml

# Activate the environment
conda activate genome_constraint_envs_O2_R
```

## Environment Details

### Python Environment (`genome_constraint_envs_O2_py`)

**Core packages:**
- Python 3.13.1
- Data science: numpy, pandas, scipy, matplotlib, seaborn
- Bioinformatics: biopython, pysam, gffutils, pyfaidx
- Machine learning: scikit-learn, xgboost, lightgbm
- Phylogeny: dendropy
- Workflow: snakemake, dask, joblib
- Jupyter: jupyter, ipython, notebook

**Known incompatibilities with Python 3.13:**
- `pybedtools` (requires Python <3.8)
- `pyvcf` (requires Python <3.13)
- `pybigwig` (requires Python <3.13)
- `ete3` (installed but has cgi module issues)
- `graphviz`/`pygraphviz` (installation failed)

### R Environment (`genome_constraint_envs_O2_R`)

**Core packages:**
- R 4.4.2
- Data manipulation: tidyverse, data.table
- Visualization: ggplot2, plotly, ggpubr
- Bioinformatics: Biostrings, GenomicRanges
- Phylogeny: ape, phytools
- Statistical modeling: lme4, mgcv, phylolm, caper
- Reporting: rmarkdown, knitr

**Known issues:**
- `bioconductor-rtracklayer` (installed but fails to load)
- `bioconductor-clusterprofiler` (R version compatibility)
- `r-caret` (installed but fails to load)

## For SLURM Job Scripts

When submitting jobs to the O2 cluster, include environment setup in your job script:

```bash
#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-04:00
#SBATCH --mem=4G
#SBATCH -c 1

# Load O2 modules
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0

# Initialize and activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate genome_constraint_envs_O2_py

# Your Python script here
python your_script.py
```

## Troubleshooting

### Conda command not found
- Make sure you've loaded the conda module: `module load conda/miniforge3/24.11.3-0`
- Initialize conda: `source $(conda info --base)/etc/profile.d/conda.sh`

### Environment not found
- Check if environment exists: `conda env list`
- If missing, create it using the instructions above

### Package import errors
- Verify package is in the `_success.yaml` file
- Some packages may have compatibility issues (see notes in YAML files)
- Try reinstalling: `conda install package_name`

### Module conflicts
- Don't load both conda Python and O2 Python module simultaneously
- Use conda Python from the environment, not `module load python/3.13.1`

