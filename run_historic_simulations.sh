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
    --folder_id=*)
        FOLDER_ID="${1#*=}"
        ;;
    --num_cores=*)
        NUM_CORES="${1#*=}"
        ;;
    --stop_importation)
        STOP_IMPORTATION=true
        ;;
    --help)
        echo "Usage: $0 [options]"
        echo ""
        echo "Required arguments:"
        echo "  --id=<id>              SLURM batch/task ID for fitting and historic simulations"
        echo "  --folder_id=<folder>   Folder for realocation (e.g., 'source-data-20250220')"
        echo ""
        echo "Optional arguments:"
        echo "  --failed_ids=<ids>     Comma-separated list of failed batch/task IDs to skip"
        echo "  --num_cores=<n>        Number of CPU cores to use for projections (default: 10)"
        echo "  --stop_importation     Stop importation of infections based on IU-specific year"
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

# Check if we're in the correct directory
EXPECTED_DIR="${TRACHOMA_AMIS_DIR}"
if [ -z "$EXPECTED_DIR" ]; then
    echo "Error: TRACHOMA_AMIS_DIR environment variable is not set"
    exit 1
fi

CURRENT_DIR="$(pwd)"
if [ "$CURRENT_DIR" != "$EXPECTED_DIR" ]; then
    echo "Error: This script must be run from '$EXPECTED_DIR', but it is running from '$CURRENT_DIR'."
    exit 1
fi

# Build common arguments for underlying scripts
COMMON_ARGS=""
if [ ! -z "$FAILED_IDS" ]; then
    COMMON_ARGS="$COMMON_ARGS --failed_ids=$FAILED_IDS"
fi

# Export task ID for sub-processes
export SLURM_ARRAY_TASK_ID=$ID

echo "Preparing histories and maps..."
Rscript prepare_histories_and_maps.R --id="$ID" || exit 1
Rscript prepare_histories_projections.R --id="$ID" || exit 1

# Run the fitting process
echo "Running trachoma fitting and preprocessing projections..."
./run_fit_and_process.sh --id="$ID" --folder_id="$FOLDER_ID" $COMMON_ARGS || exit 1

# Build arguments for historic simulation
PROJ_ARGS="--id=$ID --folder_id=$FOLDER_ID"
if [ ! -z "$NUM_CORES" ]; then
    PROJ_ARGS="$PROJ_ARGS --num_cores=$NUM_CORES"
fi
if [ ! -z "$STOP_IMPORTATION" ]; then
    PROJ_ARGS="$PROJ_ARGS --stop_importation"
fi

# Run the historic simulation
echo "Running historic simulations..."
./run_proj.sh $PROJ_ARGS || exit 1
