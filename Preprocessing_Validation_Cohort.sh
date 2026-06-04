#!/bin/bash

set -euo pipefail

SAMPLES=(
SRR21161944
SRR21161802
SRR21161882
SRR21161948
SRR21161942
SRR21161932
)

for SAMPLE in "${SAMPLES[@]}"; do

    base=$SAMPLE

    cutadapt \
        -a AGATCGGAAGAGC \
        -A AGATCGGAAGAGC \
        -q 20 -m 50 \
        -o ${base}_1_trimmed.fastq.gz \
        -p ${base}_2_trimmed.fastq.gz \
        ${base}_1.fastq.gz \
        ${base}_2.fastq.gz

    bowtie -p 6 \
        -x /mnt/e/wsl_work/human_index/hg38 \
        -1 ${base}_1_trimmed.fastq.gz \
        -2 ${base}_2_trimmed.fastq.gz \
        --un-conc ${base}_nonhuman.fastq \
        -S /dev/null

    if [ -f "${base}_nonhuman.1.fastq" ]; then

        seqtk seq -a ${base}_nonhuman.1.fastq > ${base}_1_clean.fasta
        seqtk seq -a ${base}_nonhuman.2.fastq > ${base}_2_clean.fasta

        cat ${base}_1_clean.fasta ${base}_2_clean.fasta > ${base}_clean.fasta

        rm -f ${base}_1_trimmed.fastq.gz ${base}_2_trimmed.fastq.gz
        rm -f ${base}_nonhuman.1.fastq ${base}_nonhuman.2.fastq
        rm -f ${base}_1.fastq.gz ${base}_2.fastq.gz

    fi

done
