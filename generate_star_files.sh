#!/bin/bash
set -e

mkdir -p references/star_index

STAR --runMode genomeGenerate \
    --genomeDir references/star_index \
    --genomeFastaFiles references/dmel-all-chromosome-r6.63.fasta \
    --sjdbGTFfile references/dmel-all-r6.63.gtf \
    --runThreadN 4 \
    --genomeSAindexNbases 12

