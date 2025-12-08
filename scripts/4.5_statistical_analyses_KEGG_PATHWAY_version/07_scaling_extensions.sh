#!/bin/bash
# Wrapper to run Script 07 for all data types

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
echo "Script 07: Scaling Extensions - ALL DATA TYPES"
echo "================================================================================"
echo ""

for DATA_TYPE in reactions ko pathway; do
    echo ""
    echo "================================================================================"
    echo "Processing: ${DATA_TYPE}"
    echo "================================================================================"
    python 07_scaling_extensions_single.py --data-type ${DATA_TYPE} ${TEST_MODE} ${PREVALENCE}
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to process ${DATA_TYPE}"
        exit 1
    fi
done

echo ""
echo "================================================================================"
echo "Script 07: All extensions completed!"
echo "================================================================================"
