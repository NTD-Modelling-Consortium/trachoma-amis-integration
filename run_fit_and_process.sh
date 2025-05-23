#!/bin/bash

# Parse named arguments
while [ $# -gt 0 ]; do
    case "$1" in
    --id=*)
        ID="${1#*=}"
        ;;
    --failed_ids=*)
        FAILED_IDS="${1#*=}"
        ;;
    --amis_sigma=*)
        AMIS_SIGMA="${1#*=}"
        ;;
    --folder_id=*)
        FOLDER_ID="${1#*=}"
        ;;
    --help)
        echo "Usage: $0 [options]"
        echo ""
        echo "Required arguments:"
        echo "  --id=<id>              ID for fitting and processing"
        echo "  --folder_id=<folder>   Folder ID for realocation (e.g., 'source-data-20250220')"
        echo ""
        echo "Optional arguments:"
        echo "  --failed_ids=<ids>     Comma-separated list of failed IDs to skip"
        echo "  --amis_sigma=<number>       AMIS 'sigma' parameter, expects a floating point number"
        echo "  --help                 Show this help message"
        exit 0
        ;;
    *)
        echo "Error: Invalid argument $1"
        echo "Run '$0 --help' for usage information"
        exit 1
        ;;
    esac
    shift
done

# Check required arguments
if [ -z "$ID" ]; then
    echo "Error: --id is required"
    echo "Run '$0 --help' for usage information"
    exit 1
fi

if [ -z "$FOLDER_ID" ]; then
    echo "Error: --folder_id is required"
    echo "Run '$0 --help' for usage information"
    exit 1
fi

# Build command arguments for preprocess and realocate scripts
COMMON_ARGS=""
if [ ! -z "$FAILED_IDS" ]; then
    COMMON_ARGS="$COMMON_ARGS --failed_ids=$FAILED_IDS"
fi

export SLURM_ARRAY_TASK_ID=$ID

# Run the R scripts in sequence
echo "Running trachoma fitting..."
if [ ! -z "$AMIS_SIGMA" ]; then
    Rscript trachoma_fitting.R --id=$ID --amis_sigma=$AMIS_SIGMA || exit 1
else
    Rscript trachoma_fitting.R --id=$ID || exit 1
fi

echo "Running preprocessing for projections..."
Rscript preprocess_for_projections.R --id=$ID $COMMON_ARGS || exit 1

echo "Running realocation for projections..."
Rscript realocate_InputPars_MTP.R --id=$ID $COMMON_ARGS --folder_id=$FOLDER_ID || exit 1
