#!/bin/bash

set -euo pipefail

# =========================
# DIAMOND FUNCTIONAL ANNOTATION PIPELINE
# =========================
# Software versions:
# DIAMOND v2.1.9
# RefSeq Viral Protein Database (downloaded 2026-05-19)
# =========================

# =========================
# INPUT / OUTPUT
# =========================

QUERY_DIR="../prodigal/proteins"
DB="/root/viral_refseq/viral_refseq.dmnd"

DIAMOND_OUT="../diamond"
RESULTS_OUT="../results"

mkdir -p "$DIAMOND_OUT" "$RESULTS_OUT"

# =========================
# DIAMOND ANNOTATION
# =========================

for f in "${QUERY_DIR}"/*.faa; do

    [ -e "$f" ] || continue

    sample=$(basename "$f" _proteins.faa)

    diamond blastp \
        -q "$f" \
        -d "$DB" \
        -o "${DIAMOND_OUT}/${sample}_diamond.m8" \
        -e 1e-10 \
        --sensitive \
        -k 1 \
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle

done

# =========================
# FILTER HIGH-CONFIDENCE HITS
# =========================

for f in "${DIAMOND_OUT}"/*_diamond.m8; do

    [ -e "$f" ] || continue

    sample=$(basename "$f" _diamond.m8)

    awk '$3 >= 80 && $4 >= 50' "$f" \
        > "${RESULTS_OUT}/${sample}_high_confidence.tsv"

done

echo "DIAMOND PIPELINE COMPLETE"
