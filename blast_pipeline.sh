#!/bin/bash

set -euo pipefail

# =========================
# BLAST PIPELINE (SAMPLE-LEVEL)
# Software:
# seqkit v2.3.0
# BLAST+ v2.12.0+
# =========================

# =========================
# INPUT CHECK
# =========================
if [ $# -lt 1 ]; then
    echo "Usage: bash blast_pipeline.sh <SAMPLE_ID>"
    exit 1
fi

SAMPLE=$1
CLEAN_FASTA="${SAMPLE}_clean.fasta"

echo "Processing: $SAMPLE"

if [ ! -f "$CLEAN_FASTA" ]; then
    echo "Missing file: $CLEAN_FASTA"
    exit 1
fi

# =========================
# STEP 1: SPLIT READS
# =========================
seqkit split2 -s 10000 "$CLEAN_FASTA" -O ${SAMPLE}_chunks

# =========================
# STEP 2: INITIAL VIRAL BLAST
# =========================
for f in ${SAMPLE}_chunks/*.fasta; do
    blastn -query "$f" \
    -db virus_db \
    -task megablast \
    -evalue 1e-10 \
    -outfmt "6 qseqid sseqid bitscore evalue" \
    -num_threads 1 \
    -out "${f%.fasta}_viral.tsv"
done

cat ${SAMPLE}_chunks/*_viral.tsv > ${SAMPLE}_viral_stage1.tsv

# =========================
# STEP 3: VIRAL CANDIDATES
# =========================
cut -f1 ${SAMPLE}_viral_stage1.tsv | sort | uniq > ${SAMPLE}_viral_ids.txt

seqkit grep -f ${SAMPLE}_viral_ids.txt "$CLEAN_FASTA" > ${SAMPLE}_viral_candidates.fasta

# =========================
# STEP 4: SECONDARY BLAST
# =========================
blastn -query ${SAMPLE}_viral_candidates.fasta \
-db virus_db \
-task megablast \
-evalue 1e-20 \
-outfmt "6 qseqid sseqid bitscore evalue" \
-num_threads 1 \
-out ${SAMPLE}_virus_stage2.tsv

blastn -query ${SAMPLE}_viral_candidates.fasta \
-db bacteria_db \
-task megablast \
-evalue 1e-20 \
-outfmt "6 qseqid sseqid bitscore evalue" \
-num_threads 1 \
-out ${SAMPLE}_bacteria_stage2.tsv

# =========================
# STEP 5: HUMAN FILTER
# =========================
seqkit split2 -s 5000 ${SAMPLE}_viral_candidates.fasta -O human_chunks_${SAMPLE}

for f in human_chunks_${SAMPLE}/*.fasta; do
    blastn -query "$f" \
    -db human_db \
    -task megablast \
    -evalue 1e-20 \
    -max_target_seqs 5 \
    -num_threads 2 \
    -out "${f%.fasta}_human.tsv"
done

# =========================
# STEP 6: LABELING
# =========================
awk '{print $0"\tvirus"}' ${SAMPLE}_virus_stage2.tsv > ${SAMPLE}_virus_labeled.tsv
awk '{print $0"\tbacteria"}' ${SAMPLE}_bacteria_stage2.tsv > ${SAMPLE}_bacteria_labeled.tsv
cat human_chunks_${SAMPLE}/*_human.tsv | awk '{print $0"\thuman"}' > ${SAMPLE}_human_labeled.tsv

cat ${SAMPLE}_virus_labeled.tsv \
    ${SAMPLE}_bacteria_labeled.tsv \
    ${SAMPLE}_human_labeled.tsv > ${SAMPLE}_stage2_labeled.tsv

# =========================
# STEP 7: BEST HIT SELECTION
# =========================
awk '{
if(!seen[$1] || $3>max[$1]){
    max[$1]=$3;
    line[$1]=$0;
    seen[$1]=1
}
}
END{for(i in line) print line[i]}' \
${SAMPLE}_stage2_labeled.tsv > ${SAMPLE}_stage2_best.tsv

# =========================
# STEP 8: FINAL FILTER
# =========================
awk '$5=="virus" && $3>=150' ${SAMPLE}_stage2_best.tsv > ${SAMPLE}_final_viral.tsv

cut -f2 ${SAMPLE}_final_viral.tsv | sort | uniq -c | sort -nr > ${SAMPLE}_viral_counts.txt
awk '{print $2"\t"$1}' ${SAMPLE}_viral_counts.txt > ${SAMPLE}_viral_counts_final.tsv

# =========================
# STEP 9: NCBI ANNOTATION
# =========================
cut -f1 ${SAMPLE}_viral_counts_final.tsv > ${SAMPLE}_acc.txt

if [ -s ${SAMPLE}_acc.txt ]; then
    esummary -db nucleotide -id $(paste -sd, ${SAMPLE}_acc.txt) > ${SAMPLE}_virus.xml

    xtract -input ${SAMPLE}_virus.xml \
    -pattern DocumentSummary \
    -element Caption,Title > ${SAMPLE}_accession_to_name.tsv
fi

# =========================
# STEP 10: CLEAN + MERGE
# =========================
awk 'BEGIN{OFS="\t"} {gsub(/\r/,""); gsub(/\.[0-9]+$/,"",$1); print $1,$2}' \
${SAMPLE}_viral_counts_final.tsv > counts_clean_${SAMPLE}.tsv

awk 'BEGIN{OFS="\t"} {gsub(/\r/,""); gsub(/\.[0-9]+$/,"",$1); name=$2; for(i=3;i<=NF;i++) name=name" "$i; print $1,name}' \
${SAMPLE}_accession_to_name.tsv > names_clean_${SAMPLE}.tsv

LC_ALL=C sort -k1,1 counts_clean_${SAMPLE}.tsv > a_${SAMPLE}.tsv
LC_ALL=C sort -k1,1 names_clean_${SAMPLE}.tsv > b_${SAMPLE}.tsv

join -t $'\t' a_${SAMPLE}.tsv b_${SAMPLE}.tsv > ${SAMPLE}_final_table.tsv

# =========================
# STEP 11: OUTPUT
# =========================
zip ${SAMPLE}_analysis.zip \
    ${SAMPLE}_final_table.tsv \
    ${SAMPLE}_accession_to_name.tsv 2>/dev/null

rm -rf ${SAMPLE}_chunks human_chunks_${SAMPLE}
echo "DONE: $SAMPLE"
