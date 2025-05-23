#!/bin/bash

# Parse named arguments
while [ $# -gt 0 ]; do
    case "$1" in
    --id=*)
        SLURM_TASK_ID="${1#*=}"
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
        echo "  --id=<id>                      Task ID for SLURM array"
        echo "  --folder_id=<folder>           Folder ID for outputs (e.g., 'sourcedata-20250220')"
        echo ""
        echo "Optional arguments:"
        echo "  --num_cores=<n>                Number of CPU cores to use (default: 10)"
        echo "  --stop_importation             Stop importation of infections based on IU-specific year"
        echo "  --help                         Show this help message"
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

# Check if required arguments are provided
if [ -z "$SLURM_TASK_ID" ]; then
    echo "Error: --slurm_task_id is required"
    echo "Run '$0 --help' for usage information"
    exit 1
fi

if [ -z "$FOLDER_ID" ]; then
    echo "Error: --folder_id is required"
    echo "Run '$0 --help' for usage information"
    exit 1
fi

# Build command with required arguments
export SLURM_ARRAY_TASK_ID=$SLURM_TASK_ID
CMD="python RunProjectionsTo2026.py --id=$SLURM_TASK_ID --folder_id=$FOLDER_ID"

# Add optional arguments if provided
if [ ! -z "$NUM_CORES" ]; then
    CMD="$CMD --num_cores=$NUM_CORES"
fi

if [ ! -z "$STOP_IMPORTATION" ]; then
    CMD="$CMD --stop_importation"
fi

# Execute the command
conda run --no-capture-output $CMD
