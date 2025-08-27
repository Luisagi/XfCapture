#!/bin/bash

# Configuration
CORES=10
KRAKEN_MAX_JOBS=1
ALIGNMENT_MAX_JOBS=5
IQTREE_MAX_JOBS=3

# Check for auto flag
AUTO_MODE=false
if [[ "$1" == "--auto" ]]; then
    AUTO_MODE=true
    shift  # Remove --auto from arguments
fi


OUTPUT_DIR=$(grep "^  output_dir:" config.yaml | cut -d':' -f2 | tr -d ' "'"'"'')


echo "=== PHASE 1: Pipeline up to MLST + Checkpoints ==="

# Phase 1: Up to MLST (includes checkpoints)
snakemake \
    --cores $CORES \
    --use-conda \
    --resources kraken_jobs=$KRAKEN_MAX_JOBS \
    "$@"

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 1 failed"
    exit 1
fi

echo "Phase 1 completed. Checking successful samples..."

# Count successful samples
successful_count=$(grep -l "SUCCESS:" "${OUTPUT_DIR}/03.probes_reconstruction/checkpoint/"*_reconstruction_check.txt 2>/dev/null | wc -l)
total_count=$(ls "${OUTPUT_DIR}/03.probes_reconstruction/checkpoint/"*_reconstruction_check.txt 2>/dev/null | wc -l)

echo "Total samples: $total_count"
echo "Successful samples: $successful_count"

if [ $successful_count -eq 0 ]; then
    echo "No successful samples. Pipeline terminated."
    exit 0
fi

# Ask whether to continue (unless auto mode)
if [ "$AUTO_MODE" = true ]; then
    echo "Auto mode: Continuing with phylogenetic analysis..."
else
    read -p "Continue with phylogeny? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Pipeline terminated by user choice."
        exit 0
    fi
fi

echo "=== PHASE 2: Phylogenetic analysis ==="

# Phase 2: Phylogeny
snakemake phylogeny_phase \
    --cores $CORES \
    --use-conda \
    --resources iqtree_jobs=$IQTREE_MAX_JOBS alignment_jobs=$ALIGNMENT_MAX_JOBS \
    "$@"

if [ $? -eq 0 ]; then
    echo "Pipeline completed successfully!"
else
    echo "ERROR: Phase 2 failed"
    exit 1
fi