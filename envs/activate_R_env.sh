#!/bin/bash
# Script to activate R conda environment on HMS O2
# Usage: source envs/activate_R_env.sh

# Load required O2 modules
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0

# Optional: Load CUDA if needed for GPU/ML workloads
# module load cuda/12.8

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Activate the R environment
conda activate genome_constraint_envs_O2_R

# Verify activation
if [ $? -eq 0 ]; then
    echo "✓ R environment 'genome_constraint_envs_O2_R' activated"
    echo "  R version: $(R --version | head -n 1)"
    echo "  Conda environment: $CONDA_DEFAULT_ENV"
else
    echo "✗ Failed to activate R environment"
    echo "  Run: conda env create -f envs/genome_constraint_envs_O2_R_success.yaml"
    return 1
fi

