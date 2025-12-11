#!/bin/bash
# Wrapper to run Script 06 for all data types

TEST_MODE=""
PREVALENCE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --test-mode)
            TEST_MODE="--test-mode"
            shift
            ;;
        --prevalence-threshold)
            PREVALENCE="--prevalence-threshold $2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

echo ""
echo "================================================================================"
echo "Script 06: Make Scaling Figures - ALL DATA TYPES"
echo "================================================================================"
echo ""

for DATA_TYPE in reactions ko pathway; do
    echo ""
    echo "================================================================================"
    echo "Processing: ${DATA_TYPE} (normal)"
    echo "================================================================================"
    python 06_make_scaling_figures_single.py --data-type ${DATA_TYPE} ${TEST_MODE} ${PREVALENCE}
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to process ${DATA_TYPE} (normal)"
        exit 1
    fi
    
    echo ""
    echo "================================================================================"
    echo "Processing: ${DATA_TYPE} (r² filtered)"
    echo "================================================================================"
    python 06_make_scaling_figures_single.py --data-type ${DATA_TYPE} ${TEST_MODE} ${PREVALENCE} --rsq-filtered
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to process ${DATA_TYPE} (r² filtered)"
        exit 1
    fi
done

echo ""
echo "================================================================================"
echo "Script 06: All figures generated!"
echo "================================================================================"
