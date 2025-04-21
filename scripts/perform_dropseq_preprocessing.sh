#!/bin/bash
set -e

INPUT_DIR="data/merged_fastq"
OUTPUT_DIR="results/dropseq_processed"
mkdir -p ${OUTPUT_DIR}

# Process each sample
for sample in F1 F2 d1; do
    echo "Processing sample ${sample}..."
    
    # 1. Convert FASTQ to uBAM
    picard FastqToSam \
        F1=${INPUT_DIR}/${sample}_R1.fastq.gz \
        F2=${INPUT_DIR}/${sample}_R2.fastq.gz \
        O=${OUTPUT_DIR}/${sample}_unaligned.bam \
        SAMPLE_NAME=${sample}

    # 2. Tag cell barcodes
    dropseq TagBamWithReadSequenceExtended \
        INPUT=${OUTPUT_DIR}/${sample}_unaligned.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_tagged_Cell.bam \
        SUMMARY=${OUTPUT_DIR}/${sample}_tagged_Cellular.bam_summary.txt \
        BASE_RANGE=1-12 \
        BASE_QUALITY=10 \
        BARCODED_READ=1 \
        DISCARD_READ=False \
        TAG_NAME=XC \
        NUM_BASES_BELOW_QUALITY=1

    rm ${OUTPUT_DIR}/${sample}_unaligned.bam  

    # 3. Tag molecular barcodes
    dropseq TagBamWithReadSequenceExtended \
        INPUT=${OUTPUT_DIR}/${sample}_tagged_Cell.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_tagged_CellMolecular.bam \
        SUMMARY=${OUTPUT_DIR}/${sample}_tagged_Molecular.bam_summary.txt \
        BASE_RANGE=13-20 \
        BASE_QUALITY=10 \
        BARCODED_READ=1 \
        DISCARD_READ=True \
        TAG_NAME=XM \
        NUM_BASES_BELOW_QUALITY=1

    rm ${OUTPUT_DIR}/${sample}_tagged_Cell.bam 

    # 4. Filter BAM
    dropseq FilterBam \
        TAG_REJECT=XQ \
        INPUT=${OUTPUT_DIR}/${sample}_tagged_CellMolecular.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_tagged_filtered.bam

    rm ${OUTPUT_DIR}/${sample}_tagged_CellMolecular.bam 

    # 5. Trim SMART adapter
    dropseq TrimStartingSequence \
        INPUT=${OUTPUT_DIR}/${sample}_tagged_filtered.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_trimmed_smart.bam \
        OUTPUT_SUMMARY=${OUTPUT_DIR}/${sample}_adapter_trimming_report.txt \
        SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
        MISMATCHES=0 \
        NUM_BASES=5

    rm ${OUTPUT_DIR}/${sample}_tagged_filtered.bam

    # 6. Trim polyA
    dropseq PolyATrimmer \
        INPUT=${OUTPUT_DIR}/${sample}_trimmed_smart.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_polyA_filtered.bam \
        OUTPUT_SUMMARY=${OUTPUT_DIR}/${sample}_polyA_trimming_report.txt \
        MISMATCHES=0 \
        NUM_BASES=6 \
        USE_NEW_TRIMMER=true

   rm ${OUTPUT_DIR}/${sample}_trimmed_smart.bam 

    # 7. Convert to FASTQ for STAR
    picard SamToFastq \
        INPUT=${OUTPUT_DIR}/${sample}_polyA_filtered.bam \
        FASTQ=${OUTPUT_DIR}/${sample}_for_star.fastq

    echo "Finished processing ${sample}"
done
