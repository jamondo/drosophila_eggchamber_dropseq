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
cat <<EOF >> references/dmel-all-r6.63.with_transgenes.gtf
exogenous_genes	Custom	gene	1	${dsred_len}	.	+	.	gene_id "dsRedExpress"; gene_name "dsRedExpress";
exogenous_genes	Custom	exon	1	${dsred_len}	.	+	.	gene_id "dsRedExpress"; transcript_id "dsRedExpress_t"; gene_name "dsRedExpress";
exogenous_genes	Custom	gene	$((dsred_len+1))	$((dsred_len+gal4_len))	.	+	.	gene_id "Gal4"; gene_name "Gal4";
exogenous_genes	Custom	exon	$((dsred_len+1))	$((dsred_len+gal4_len))	.	+	.	gene_id "Gal4"; transcript_id "Gal4_t"; gene_name "Gal4";
exogenous_genes	Custom	gene	$((dsred_len+gal4_len+1))	$((dsred_len+gal4_len+gfp_len))	.	+	.	gene_id "LifeActGFP"; gene_name "LifeActGFP";
exogenous_genes	Custom	exon	$((dsred_len+gal4_len+1))	$((dsred_len+gal4_len+gfp_len))	.	+	.	gene_id "LifeActGFP"; transcript_id "LifeActGFP_t"; gene_name "LifeActGFP";
EOF
