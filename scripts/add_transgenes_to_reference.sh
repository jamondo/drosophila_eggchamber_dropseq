#!/bin/bash
set -e

# Check if reference files exist
if [ ! -f references/dmel-all-chromosome-r6.63.fasta ] || [ ! -f references/dmel-all-r6.63.gtf ]; then
    echo "Error: Reference files not found. Run download_reference.sh first."
    exit 1
fi

# Create copies of reference files
cp references/dmel-all-chromosome-r6.63.fasta references/dmel-all-chromosome-r6.63.with_transgenes.fasta
cp references/dmel-all-r6.63.gtf references/dmel-all-r6.63.with_transgenes.gtf

# Add transgene sequences to FASTA
echo -e "\n>exogenous_genes" >> references/dmel-all-chromosome-r6.63.with_transgenes.fasta
grep -v "^>" reference_files/exogenous_gene_sequences.txt | tr -d '\n' >> references/dmel-all-chromosome-r6.63.with_transgenes.fasta
echo "" >> references/dmel-all-chromosome-r6.63.with_transgenes.fasta

# Calculate transgene lengths
dsred_len=$(grep -A1 "^>dsRedExpress" reference_files/exogenous_gene_sequences.txt | tail -n1 | tr -d '\n' | wc -c)
gal4_len=$(grep -A1 "^>Gal4" reference_files/exogenous_gene_sequences.txt | tail -n1 | tr -d '\n' | wc -c)
gfp_len=$(grep -A1 "^>LifeActGFP" reference_files/exogenous_gene_sequences.txt | tail -n1 | tr -d '\n' | wc -c)

# Convert symbols to names in GTF
sed -i '' -e 's/gene_symbol/gene_name/g' -e 's/transcript_symbol/transcript_name/g' references/dmel-all-r6.63.with_transgenes.gtf

# Add transgene annotations to GTF
printf "exogenous_genes\tCustom\tgene\t1\t%d\t.\t+\t.\tgene_id \"dsRedExpress\"; gene_name \"dsRedExpress\";\n" ${dsred_len} >> references/dmel-all-r6.63.with_transgenes.gtf
printf "exogenous_genes\tCustom\texon\t1\t%d\t.\t+\t.\tgene_id \"dsRedExpress\"; transcript_id \"dsRedExpress_t\"; gene_name \"dsRedExpress\";\n" ${dsred_len} >> references/dmel-all-r6.63.with_transgenes.gtf
printf "exogenous_genes\tCustom\tgene\t%d\t%d\t.\t+\t.\tgene_id \"Gal4\"; gene_name \"Gal4\";\n" $((dsred_len+1)) $((dsred_len+gal4_len)) >> references/dmel-all-r6.63.with_transgenes.gtf
printf "exogenous_genes\tCustom\texon\t%d\t%d\t.\t+\t.\tgene_id \"Gal4\"; transcript_id \"Gal4_t\"; gene_name \"Gal4\";\n" $((dsred_len+1)) $((dsred_len+gal4_len)) >> references/dmel-all-r6.63.with_transgenes.gtf
printf "exogenous_genes\tCustom\tgene\t%d\t%d\t.\t+\t.\tgene_id \"LifeActGFP\"; gene_name \"LifeActGFP\";\n" $((dsred_len+gal4_len+1)) $((dsred_len+gal4_len+gfp_len)) >> references/dmel-all-r6.63.with_transgenes.gtf
printf "exogenous_genes\tCustom\texon\t%d\t%d\t.\t+\t.\tgene_id \"LifeActGFP\"; transcript_id \"LifeActGFP_t\"; gene_name \"LifeActGFP\";\n" $((dsred_len+gal4_len+1)) $((dsred_len+gal4_len+gfp_len)) >> references/dmel-all-r6.63.with_transgenes.gtf

