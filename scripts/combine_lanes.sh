#!/bin/bash
set -e

FASTQ_DIR="/Volumes/Data/RNASeq/allFastq"
OUTPUT_DIR="data/merged_fastq"

mkdir -p ${OUTPUT_DIR}

for sample in F1 F2 d1; do
    echo "Processing sample ${sample}..."
    
    cat ${FASTQ_DIR}/dropseq-${sample}_S*_L00{1,2,3,4}_R1_001.fastq.gz > ${OUTPUT_DIR}/${sample}_R1.fastq.gz
    cat ${FASTQ_DIR}/dropseq-${sample}_S*_L00{1,2,3,4}_R2_001.fastq.gz > ${OUTPUT_DIR}/${sample}_R2.fastq.gz
    
    echo "Finished merging lanes for sample ${sample}"
done
