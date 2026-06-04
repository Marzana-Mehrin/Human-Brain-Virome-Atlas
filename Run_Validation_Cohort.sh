#!/bin/bash

# =========================
# VALIDATION COHORT PIPELINE
# =========================

set -euo pipefail

echo "Starting Validation Cohort Analysis"

# =========================
# SRR SAMPLE LIST
# =========================
SAMPLES=(
SRR21161944
SRR21161802
SRR21161882
SRR21161948
SRR21161942
SRR21161932
)

# =========================
# RUN PIPELINE
# =========================
for SAMPLE in "${SAMPLES[@]}"; do
    echo "=================================="
    echo "Processing Validation Sample: $SAMPLE"
    echo "=================================="

    bash blast_pipeline.sh "$SAMPLE"

    echo "Completed: $SAMPLE"
done

echo "ALL VALIDATION SAMPLES COMPLETED"
