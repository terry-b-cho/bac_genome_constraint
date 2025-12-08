#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-04:00
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH -o logs/env_prediction_kegg_%j.out
#SBATCH -e logs/env_prediction_kegg_%j.err

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

# Parse command line arguments
PREVALENCE_THRESHOLD=${1:-99}
TEST_MODE=${2:-""}
MODEL=${3:-all}

echo "=================================================================================="
echo "KEGG Environment Prediction Pipeline"
echo "=================================================================================="
echo "Prevalence threshold: ${PREVALENCE_THRESHOLD}%"
echo "Test mode: ${TEST_MODE}"
echo "Model: ${MODEL}"
echo "=================================================================================="
echo ""

# Data types to process
DATA_TYPES=("reactions" "ko" "pathway")

# Process each data type
for data_type in "${DATA_TYPES[@]}"; do
    echo ""
    echo "=================================================================================="
    echo "Processing ${data_type} data type..."
    echo "=================================================================================="
    
    # Build command
    CMD="python3 scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version.py"
    CMD="${CMD} --data-type ${data_type}"
    CMD="${CMD} --prevalence-threshold ${PREVALENCE_THRESHOLD}"
    CMD="${CMD} --model ${MODEL}"
    CMD="${CMD} --plot"
    
    if [ "${TEST_MODE}" == "test" ]; then
        CMD="${CMD} --test-mode"
        echo "  Running in TEST MODE"
    fi
    
    echo "  Command: ${CMD}"
    echo ""
    
    # Run the command
    eval ${CMD}
    
    EXIT_CODE=$?
    if [ $EXIT_CODE -ne 0 ]; then
        echo ""
        echo "  ✗ ERROR: Failed to process ${data_type} data type (exit code: ${EXIT_CODE})"
        echo "  Continuing with next data type..."
    else
        echo ""
        echo "  ✓ Successfully completed ${data_type} data type"
    fi
    
    echo ""
done

echo "=================================================================================="
echo "All data types processed!"
echo "=================================================================================="
echo ""
echo "Results saved to: results/5.5_environment_prediction/"
echo "  Files have suffixes: _reactions, _ko, _pathway"
echo ""
