#!/bin/bash
set -e

mkdir -p references

wget -P references https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.63_FB2025_02/fasta/dmel-all-chromosome-r6.63.fasta.gz
wget -P references https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.63_FB2025_02/gtf/dmel-all-r6.63.gtf.gz
gunzip references/*.gz

# Add transgenes
echo -e "\n>exogenous_genes" >> references/dmel-all-chromosome-r6.63.fasta
grep -v "^>" reference_files/exogenous_gene_sequences.txt | tr -d '\n' >> references/dmel-all-chromosome-r6.63.fasta
echo "" >> references/dmel-all-chromosome-r6.63.fasta

# Calculate gene positions and add to GTF
dsred_len=$(grep -A1 "^>dsRedExpress" reference_files/exogenous_gene_sequences.txt | tail -n1 | tr -d '\n' | wc -c)
gal4_len=$(grep -A1 "^>Gal4" reference_files/exogenous_gene_sequences.txt | tail -n1 | tr -d '\n' | wc -c)
gfp_len=$(grep -A1 "^>LifeActGFP" reference_files/exogenous_gene_sequences.txt | tail -n1 | tr -d '\n' | wc -c)

# Add GTF entries
cat << EOF >> references/dmel-all-r6.63.gtf

exogenous_genes	Custom	gene	1	${dsred_len}	.	+	.	gene_id "dsRedExpress"; gene_name "dsRedExpress"
exogenous_genes	Custom	exon	1	${dsred_len}	.	+	.	gene_id "dsRedExpress"; transcript_id "dsRedExpress_t"; gene_name "dsRedExpress"

exogenous_genes	Custom	gene	$((dsred_len+1))	$((dsred_len+gal4_len))	.	+	.	gene_id "Gal4"; gene_name "Gal4"
exogenous_genes	Custom	exon	$((dsred_len+1))	$((dsred_len+gal4_len))	.	+	.	gene_id "Gal4"; transcript_id "Gal4_t"; gene_name "Gal4"

exogenous_genes	Custom	gene	$((dsred_len+gal4_len+1))	$((dsred_len+gal4_len+gfp_len))	.	+	.	gene_id "LifeActGFP"; gene_name "LifeActGFP"
exogenous_genes	Custom	exon	$((dsred_len+gal4_len+1))	$((dsred_len+gal4_len+gfp_len))	.	+	.	gene_id "LifeActGFP"; transcript_id "LifeActGFP_t"; 
gene_name "LifeActGFP"
EOF
