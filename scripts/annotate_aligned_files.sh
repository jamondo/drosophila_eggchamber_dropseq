#!/bin/bash
set -e

if [ $# -ne 3 ]; then
    echo "Usage: $0 <aligned_bam_dir> <tagged_bam_dir> <output_dir>"
    echo "Example: $0 results/star_aligned results/dropseq_processed results/final_dge"
    exit 1
fi

ALIGNED_DIR="$1"
TAGGED_DIR="$2"
OUTPUT_DIR="$3"
TMP_DIR="/Volumes/Data/RNASeq/tmp"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${TMP_DIR}  # make sure TMP_DIR exists

for aligned_bam in ${ALIGNED_DIR}/*_Aligned.sortedByCoord.out.bam; do
    sample=$(basename $aligned_bam _Aligned.sortedByCoord.out.bam)
    echo "Processing ${sample}..."

    # 1. Sort aligned BAM by queryname
    picard SortSam \
        I=${aligned_bam} \
        O=${OUTPUT_DIR}/${sample}_queryname_sorted.bam \
        SORT_ORDER=queryname \
        TMP_DIR=${TMP_DIR}

    # 2. Merge aligned BAM with cell/molecular barcodes
    picard MergeBamAlignment \
        REFERENCE_SEQUENCE=references/dmel-all-chromosome-r6.63.fasta \
        UNMAPPED_BAM=${TAGGED_DIR}/${sample}_polyA_filtered.bam \
        ALIGNED_BAM=${OUTPUT_DIR}/${sample}_queryname_sorted.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_merged.bam \
        INCLUDE_SECONDARY_ALIGNMENTS=false \
        PAIRED_RUN=false \
        TMP_DIR=${TMP_DIR}

    # 3. Tag with gene exons
    dropseq TagReadWithGeneExon \
        I=${OUTPUT_DIR}/${sample}_merged.bam \
        O=${OUTPUT_DIR}/${sample}_gene_tagged.bam \
        ANNOTATIONS_FILE=references/dmel-all-r6.63.gtf \
        TAG=GE \
        TMP_DIR=${TMP_DIR}

    # 4. Create Digital Gene Expression matrix
    dropseq DigitalExpression \
        I=${OUTPUT_DIR}/${sample}_gene_tagged.bam \
        O=${OUTPUT_DIR}/${sample}_dge.txt.gz \
        SUMMARY=${OUTPUT_DIR}/${sample}_dge_summary.txt \
        MIN_NUM_GENES_PER_CELL=100 \
        TMP_DIR=${TMP_DIR}

    # Clean up intermediate files
    rm ${OUTPUT_DIR}/${sample}_queryname_sorted.bam
    rm ${OUTPUT_DIR}/${sample}_merged.bam

    echo "Finished processing ${sample}"
done

