Pipeline for processing Drop-seq single-cell RNA-seq data from Drosophila egg chambers.

## Setup

```bash
# Create conda environment
mamba env create -f environment.yml

# Activate environment
conda activate drosophila_dropseq

# Clone and build Drop-seq tools
bash scripts/setup_dropseq.sh

# Download Drosophila genome and add custom transgenes
bash scripts/download_reference.sh

# Generate STAR index
bash scripts/generate_star_files.sh

# Merge sequencing lanes for each sample
bash scripts/combine_lanes.sh

# Process samples through Drop-seq pipeline
bash scripts/perform_dropseq_preprocessing.sh
