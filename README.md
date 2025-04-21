Pipeline for processing Drop-seq single-cell RNA-seq data from Drosophila egg chambers.

## Setup

```bash
# Create conda environment
mamba env create -f environment.yml

# Activate environment
conda activate drosophila_dropseq

# Clone and build Drop-seq tools
bash scripts/setup_dropseq.sh
