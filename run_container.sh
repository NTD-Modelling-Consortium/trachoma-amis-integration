#!/bin/bash
# run_container.sh - Simple Docker wrapper for Trachoma AMIS Pipeline

set -e

IMAGE="trachoma-amis-pipeline:latest"

# Base mounts for artefacts (always mount fitting, projections-prep, projections)
MOUNTS=(
    "-v" "$PWD/fitting/artefacts:/ntdmc/trachoma-amis-integration/fitting/artefacts"
    "-v" "$PWD/projections-prep/artefacts:/ntdmc/trachoma-amis-integration/projections-prep/artefacts" 
    "-v" "$PWD/projections/artefacts:/ntdmc/trachoma-amis-integration/projections/artefacts"
)

# Smart fitting-prep mounting: only mount if host directory is non-empty
# The container has built-in fitting-prep artefacts from the Docker build process.
# We only override them if the host directory contains files.
if [ "$(ls -A fitting-prep/artefacts 2>/dev/null)" ]; then
    echo "Using host fitting-prep artefacts (host directory is non-empty, will override container's built-in artefacts)"
    MOUNTS+=("-v" "$PWD/fitting-prep/artefacts:/ntdmc/trachoma-amis-integration/fitting-prep/artefacts")
else
    echo "Using container's built-in fitting-prep artefacts (host directory is empty, preserving container's pre-built artefacts)"
fi

# Generate unique container name
CONTAINER_NAME="trachoma-amis-pipeline-$(date +%s)"

# Run container and cleanup
docker run --name "$CONTAINER_NAME" "${MOUNTS[@]}" "$IMAGE" "$@"
docker rm "$CONTAINER_NAME" >/dev/null 2>&1