##############
# Load required libraries
###############
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

###################
# Function to read DGE matrix and create Seurat object
##################
create_seurat_object <- function(file_path, sample_id) {
    counts <- read.table(file_path, header = TRUE, row.names = 1)
    seurat_obj <- CreateSeuratObject(
        counts = as.sparse(counts),
        project = "EggChamber",
        min.cells = 3,   
        min.features = 200 
    )
    
    seurat_obj$sample_id <- sample_id
    seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
    
    return(seurat_obj)
}

# Read in all samples
sample_paths <- c(
    "/Volumes/Data/RNASeq/drosophila_eggchamber_dropseq/results/final_dge/sample1_dge.txt.gz",
    "/Volumes/Data/RNASeq/drosophila_eggchamber_dropseq/results/final_dge/sample2_dge.txt.gz",
    "/Volumes/Data/RNASeq/drosophila_eggchamber_dropseq/results/final_dge/sample3_dge.txt.gz"
)

# Arbitrary labels coming out of the initial runs. They were all done on different days
sample_ids <- c("d1", "F1", "F2")  

# Create list of Seurat objects
seurat_list <- mapply(
    create_seurat_object,
    sample_paths,
    sample_ids,
    SIMPLIFY = FALSE
)

merged_seurat <- merge(seurat_list[[1]], 
                      y = seurat_list[2:length(seurat_list)], 
                      add.cell.ids = sample_ids)

#####################
# Calculate QC metrics
####################
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")


qc_plots <- function(seurat_obj) {
    
    p1 <- VlnPlot(seurat_obj, 
                  features = "nFeature_RNA", 
                  group.by = "sample_id", 
                  pt.size = 0) +
        ggtitle("Genes per Cell")

    p2 <- VlnPlot(seurat_obj, 
                  features = "nCount_RNA", 
                  group.by = "sample_id", 
                  pt.size = 0) +
        ggtitle("UMI Counts per Cell")

    p3 <- VlnPlot(seurat_obj, 
                  features = "percent.mt", 
                  group.by = "sample_id", 
                  pt.size = 0) +
        ggtitle("Percent Mitochondrial")

    # Scatter plots
    p4 <- FeatureScatter(seurat_obj, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA", 
                        group.by = "sample_id") +
        ggtitle("Genes vs UMIs")

    p5 <- FeatureScatter(seurat_obj, 
                        feature1 = "percent.mt", 
                        feature2 = "nFeature_RNA", 
                        group.by = "sample_id") +
        ggtitle("Genes vs Mt%")

    # Return combined plots
    return(wrap_plots(list(p1, p2, p3, p4, p5), ncol = 2))
}

print("QC plots before filtering:")
qc_plots(merged_seurat)

#####################################
# Print summary statistics per sample
####################################
print("Before filtering:")
print(table(merged_seurat$sample_id))
summary_stats <- merged_seurat@meta.data %>%
    group_by(sample_id) %>%
    summarise(
        Median_Genes = median(nFeature_RNA),
        Median_UMIs = median(nCount_RNA),
        Median_Mt_Percent = median(percent.mt),
        Cell_Count = n()
    )
print(summary_stats)

#####################
# Filter cells
#####################
merged_seurat <- subset(merged_seurat, 
    subset = nFeature_RNA > 500 &      # Minimum 500 genes per cell
             nFeature_RNA < 4000 &      # Maximum 4000 genes per cell
             percent.mt < 10 &          # Less than 10% mitochondrial
             nCount_RNA > 1000 &        # Minimum 1000 UMIs
             nCount_RNA < 15000         # Maximum 15000 UMIs
)

print("QC plots after filtering:")
qc_plots(merged_seurat)

print("After filtering:")
print(table(merged_seurat$sample_id))
summary_stats <- merged_seurat@meta.data %>%
    group_by(sample_id) %>%
    summarise(
        Median_Genes = median(nFeature_RNA),
        Median_UMIs = median(nCount_RNA),
        Median_Mt_Percent = median(percent.mt),
        Cell_Count = n()
    )
print(summary_stats)

#############################
# Normalize and scale data
#############################
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)

##########################
# Dimensionality reduction
#########################
merged_seurat <- RunPCA(merged_seurat)

# Show Elbow plot
ElbowPlot(merged_seurat, ndims = 50)

# Run UMAP
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)

# Display UMAP colored by sample
DimPlot(merged_seurat, reduction = "umap", group.by = "sample_id")

##########################
# Save Seurat
#########################