#!/bin/bash

set -euo pipefail

# =========================
# DISCOVERY COHORT SRR IDS
# =========================
SAMPLES=(
SRR1424656 SRR1424657 SRR1424658 SRR1424659 SRR1424660
SRR1424661 SRR1424662 SRR1424663 SRR1424664 SRR1424665
SRR1424666 SRR1424667 SRR1424668 SRR1424669 SRR1424670
SRR1424671 SRR1424672 SRR1424673 SRR1424674 SRR1424675
SRR1424676 SRR1424677 SRR1424678 SRR1424679 SRR1424680
SRR1424681 SRR1424682 SRR1424683 SRR1424684 SRR1424685
SRR1424686 SRR1424687 SRR1424688 SRR1424689 SRR1424690
SRR1424691 SRR1424692 SRR1424693 SRR1424694 SRR1424695
SRR1424696 SRR1424697 SRR1424698 SRR1424699 SRR1424700
SRR1424701 SRR1424702 SRR1424703 SRR1424704 SRR1424705
SRR1424706 SRR1424707
)

# =========================
# DOWNLOAD RAW DATA
# =========================
for S in "${SAMPLES[@]}"; do
    echo "Downloading $S"

    wget -nc \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR142/${S: -3:1}0${S: -2}/${S}/${S}.fastq.gz
done

# =========================
# TRIMMING (Cutadapt)
# =========================
for f in *.fastq.gz; do
    cutadapt -a AGATCGGAAGAGC -q 20 -m 50 \
    -o "${f%.fastq.gz}_trimmed.fastq.gz" "$f"
done

# =========================
# ALIGNMENT (Bowtie)
# =========================
for f in *_trimmed.fastq.gz; do
    base=${f%_trimmed.fastq.gz}

    bowtie -p 4 -x human_index/hg38 \
    -U "$f" | samtools view -bS - > "${base}.bam"
done

# =========================
# REMOVE HUMAN READS
# =========================
for f in *.bam; do
    base=${f%.bam}
    samtools view -b -f 4 "$f" > "${base}_unmapped.bam"
done

# =========================
# CONVERT TO FASTA
# =========================
for f in *_unmapped.bam; do
    base=${f%_unmapped.bam}

    samtools fastq "$f" | seqtk seq -a - > "${base}_clean.fasta"
done

echo "PREPROCESSING COMPLETE"
