#!/bin/bash
# run-docker.sh

set -e # Exit on any error

IMAGE="trachoma-amis-pipeline:latest"
CONTAINER_NAME="trachoma-amis-pipeline-run-$(date +%s)"

# Define directory mappings
HOST_DIRS=(
    "./artefacts/fitting-prep"
    "./artefacts/fitting"
    "./artefacts/projections-prep"
    "./artefacts/projections"
)

CONTAINER_DIRS=(
    "/ntdmc/trachoma-amis-integration/fitting-prep/artefacts"
    "/ntdmc/trachoma-amis-integration/fitting/artefacts"
    "/ntdmc/trachoma-amis-integration/projections-prep/artefacts"
    "/ntdmc/trachoma-amis-integration/projections/artefacts"
)

# Cleanup function
cleanup() {
    echo "Cleaning up container..."
    docker rm "$CONTAINER_NAME" >/dev/null 2>&1 || true
}

# Set trap for cleanup on exit
trap cleanup EXIT

# Create host directories
for host_dir in "${HOST_DIRS[@]}"; do
    mkdir -p "$host_dir"
done

# Run container
echo "Running pipeline..."
if docker run --name "$CONTAINER_NAME" "$IMAGE" "$@"; then
    echo "Pipeline completed successfully"

    # Copy data to host only if successful
    echo "Copying results to host..."
    for i in "${!HOST_DIRS[@]}"; do
        host_dir="${HOST_DIRS[i]}"
        container_dir="${CONTAINER_DIRS[i]}"
        if docker cp "$CONTAINER_NAME:$container_dir/." "$host_dir/" 2>/dev/null; then
            echo "✓ Copied $container_dir -> $host_dir"
        else
            echo "⚠ Warning: Could not copy $container_dir"
        fi
    done

    echo "Results available in ./artefacts/"
else
    echo "❌ Pipeline failed!"
    echo "Container logs:"
    docker logs "$CONTAINER_NAME" | tail -50
    exit 1
fi
