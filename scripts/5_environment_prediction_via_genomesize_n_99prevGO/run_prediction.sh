#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-04:00
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH -o logs/env_prediction_%j.out
#SBATCH -e logs/env_prediction_%j.err

# Load O2 modules
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0
module load cuda/12.8  # For GPU support (optional)

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Activate environment (adjust name if needed)
conda activate genome_constraint_envs_O2_py || {
    echo "ERROR: Conda environment not found. Please create it first."
    exit 1
}

# Change to project directory
cd /n/scratch/users/b/byc014/github/bac_genome_constraint

# Create logs directory
mkdir -p logs

# Run prediction script with all models
echo "Starting environment prediction analysis..."
python3 scripts/5_environment_prediction_via_genomesize_n_99prevGO/5_environment_prediction_via_genomesize_n_99prevGO.py \
    --prevalence-threshold 99 \
    --model all \
    --plot

echo "Analysis complete!"

