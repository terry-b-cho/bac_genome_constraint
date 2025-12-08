#!/bin/bash
# Wrapper to run minimal plotting for all data types

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
echo "Script 04.5: Minimal Plotting - ALL DATA TYPES"
echo "================================================================================"
echo ""

for DATA_TYPE in reactions ko pathway; do
    echo ""
    echo "Processing: ${DATA_TYPE}"
    python 04.5_plot_minimal_single.py --data-type ${DATA_TYPE} ${TEST_MODE} ${PREVALENCE}
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to process ${DATA_TYPE}"
        exit 1
    fi
done

echo ""
echo "================================================================================"
echo "Script 04.5: All plots generated!"
echo "================================================================================"
