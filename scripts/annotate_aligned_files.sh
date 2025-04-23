#!/bin/bash
set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <aligned_bam_dir> <tagged_bam_dir>"
    echo "Example: $0 results/star_aligned results/dropseq_processed"
    exit 1
fi

ALIGNED_DIR="$1"
TAGGED_DIR="$2"
OUTPUT_DIR="results/final_dge"
mkdir -p ${OUTPUT_DIR}

for aligned_bam in ${ALIGNED_DIR}/*_Aligned.sortedByCoord.out.bam; do
    sample=$(basename $aligned_bam _Aligned.sortedByCoord.out.bam)
    echo "Processing ${sample}..."

    # 1. Sort aligned BAM by queryname
    picard SortSam \
        I=${aligned_bam} \
        O=${OUTPUT_DIR}/${sample}_queryname_sorted.bam \
        SORT_ORDER=queryname

    # 2. Merge aligned BAM with cell/molecular barcodes
    picard MergeBamAlignment \
        REFERENCE_SEQUENCE=references/genome.fasta \
        UNMAPPED_BAM=${TAGGED_DIR}/${sample}_polyA_filtered.bam \
        ALIGNED_BAM=${OUTPUT_DIR}/${sample}_queryname_sorted.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_merged.bam \
        INCLUDE_SECONDARY_ALIGNMENTS=false \
        PAIRED_RUN=false

    # 3. Tag with gene exons
    dropseq TagReadWithGeneExon \
        I=${OUTPUT_DIR}/${sample}_merged.bam \
        O=${OUTPUT_DIR}/${sample}_gene_tagged.bam \
        ANNOTATIONS_FILE=references/genes.gtf \
        TAG=GE

    # 4. Create Digital Gene Expression matrix
    dropseq DigitalExpression \
        I=${OUTPUT_DIR}/${sample}_gene_tagged.bam \
        O=${OUTPUT_DIR}/${sample}_dge.txt.gz \
        SUMMARY=${OUTPUT_DIR}/${sample}_dge_summary.txt \
        MIN_NUM_GENES_PER_CELL=100

    # Clean up intermediate files
    rm ${OUTPUT_DIR}/${sample}_queryname_sorted.bam
    rm ${OUTPUT_DIR}/${sample}_merged.bam

    echo "Finished processing ${sample}"
done

