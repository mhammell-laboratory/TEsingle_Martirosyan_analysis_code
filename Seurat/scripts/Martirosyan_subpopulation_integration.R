#---------- Setting Up R Script -----------

#loading required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(presto) #for running findallmarker() imputations faster
library(sctransform)
library(glmGamPoi) #needed for sctransform to run sctransform to run faster
library(Matrix)
library(scales)
library(BPCells)
#setting Memory limit and bitmap type
options(future.globals.maxSize = 600 * 1024^3)

#get RDS object containing subsetted broad cell type (e.g. Neurons)
if(length(args) < 1){
    stop("Usage: Rscript Martirosyan_subpopulation_integration.R [celltype RDS object]", call.=FALSE)
}
infile <- args[1]
base <- tools:file_path_sans_ext(infile)

#setting output folders
rds_dir  <- "subpopulation/rds"
table_dir <- "subpopulation/tables"
plots_dir <- "subpopulation/plots"

#---------- Loading in Seurat Object -----------

#loading in subsetted broad cell type Seurat Object RDS
cat("Loading Seurat Object RDS...\n")
#change "BroadCellType" to the name of the broad cell type Seurat Object RDS file that you want to use for subpopulation integration
so <- readRDS(infile)

#---------- Subsetting Seurat Object based on nCount_RNA and nFeature_RNA -----------

#subsetting cells with at least 1000 UMIs per cell
cat("Subsetting Seurat Object based on nCount_RNA > 1000...\n")
so <- subset(so, subset = nCount_RNA > 1000)

#---------- Running SCTransform on SO -----------

#optional step:
#joining the RNA assay layers only if SO has multiple layers
cat("Joining RNA assay layers...\n")
so[["RNA"]] <- JoinLayers(so[["RNA"]])

#splitting the object by diagnosis
cat("Splitting Seurat Object by diagnosis...\n")
so[["RNA"]] = split(so[["RNA"]], f = so$diagnosis)

#performing scTransform on the splidt SO
cat("Running SCTransform on each subset of Seurat Object...\n")
so <- SCTransform(so, verbose = TRUE)

#saving SCTransformed RDS
cat("Saving Seurat Object sctransform RDS...\n")
outfile <- paste0(base,"_sct.rds")
saveRDS(so, file = (file.path(rds_dir, outfile)))

#---------- Performing PCA on SO SCT data -----------

#performing dimensionality reduction
cat("Running PCA on Seurat Object sct data...\n")
so <- RunPCA(so, verbose = FALSE, npcs = 30)

#---------- Integrating SO with CCA -----------

#creating all genes variable to use for integration
#genes with ":" in the name are TEs, so removing these from the integration features
cat("Creating all genes without TEs variable...\n")
all_genes <- rownames(so)[!grepl(":",rownames(so))]
cat("Number of genes for integration:", length(all_genes), "\n")

#integrating data with ccintegration
cat("Integrating layers with CCA integration...\n")
so <- IntegrateLayers(
  so, method = CCAIntegration, 
  orig.reduction = 'pca', normalization.method = 'SCT', 
  features = all_genes, new.reduction = 'integrated.cca', dims = 1:30)
joining layers after integration
cat("Joining layers after CCA integration...\n")
so[["RNA"]] <- JoinLayers(so[["RNA"]])

#saving RDS so can load in later
cat("Saving Seurat object with CCA integration...\n")
outfile <- paste0(base,"_sct_cca.rds")
saveRDS(so, file = (file.path(rds_dir, outfile)))

#---------- Finding Neighbors and Clusters with Leiden -----------

#finding neighbors with integrated.cca reduction, dimensions 1:30
cat("Finding neighbors...\n")
so <- FindNeighbors(so, reduction = "integrated.cca", dims = 1:30, verbose = FALSE)

#finding clusters with Leiden
cat("Finding clusters with Leiden...\n")
so <- FindClusters(so, resolution = 0.5, algorithm = 4, cluster.name = "leiden_subpopulation", verbose = TRUE)

#dimension reduction with UMAP
cat("Running UMAP...\n")
so <- RunUMAP(so, reduction = "integrated.cca", reduction.name = "umap.cca", dims = 1:30)
DefaultAssay(so) = "RNA"

#setting identity to leiden_subpopulation
cat("Setting active identity to leiden_subpopulation and creating UMAP...\n")
Idents(so) <- "leiden_subpopulation"
#creating a plot of all the labeled clusters from leiden
so_leiden_UMAP <- DimPlot(so, reduction = "umap.cca", shuffle = TRUE, label = TRUE)
#saving UMAP
outplot <- paste0(base,"_Leiden_UMAP.png")
ggsave(filename = file.path(plots_dir, outplot),
       plot = so_leiden_UMAP,
       width = 10,
       height = 6,
       dpi = 300)

#create UMAPs for different metadata categories
metadata_categories <- c("diagnosis", "gender", "patientID")
for (category in metadata_categories) {
  if (category %in% colnames(so@meta.data)) {
    cat(paste("Setting active identity to", category, "and creating UMAP...\n"))
    # Create UMAP plot
    plot <- DimPlot(so, reduction = "umap.cca", group.by = category, label = TRUE, repel = TRUE, pt.size = 0.5) +
      ggtitle(paste("TITLE_OF_UMAP:", category)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Save plot
    outplot <- paste0(base,"_", category,"_UMAP.png")
    ggsave(filename = file.path(plots_dir, outplot),
           plot = plot,
           width = 10,
           height = 6,
           dpi = 300)
  }
  else {
    cat(paste("Warning:", category, "not found in metadata\n"))
  }
}

#saving RDS with leiden subpopulation clusters
cat("Saving Seurat object with Leiden subpopulation clusters...\n")
Idents(so) <- "leiden_subpopulation"
outfile <- paste0(base,"_integrated.rds")
saveRDS(so, file = (file.path(rds_dir, outfile)))

#---------- Seurat FindMarkers() to Define Clusters -----------

#features to use to find cluster markers
cat("Using only genes to define cluster markers...\n")
all_genes <- rownames(so)[grep(":",rownames(so),invert=T)]
Idents(so) <- "leiden_subpopulation"

#Finding all markers DE in all clusters with LR
cat("Finding cluster markers using LR method...\n")
so_leidenmarkers <- FindAllMarkers(so, test.use= "LR", features = all_genes, min.pct= 0.1, logfc.threshold = 0.5)

#adding column to indicate if gene is sig up, down, or not DE then arranging the df so that it orders genes by cluster, then by the smallest p_val_adj, and finally by the absoulte value of (avg_log2FC, large to small)
cat("Arranging markers and adding regulation column...\n")
so_leidenmarkers <- so_leidenmarkers %>%
  arrange(cluster, p_val_adj, desc(abs(avg_log2FC))) %>%
  mutate(regulation = case_when(
    avg_log2FC > 1 & p_val_adj < 0.05 ~ "up",
    avg_log2FC < -1 & p_val_adj < 0.05 ~ "down",
    TRUE ~ "none"))
head(so_leidenmarkers)

#saving table of the top genes up and down regulated in each cluster
cat("Saving Leiden Cluster Markers...\n")
outfile <- paste0(base,"_LeidenMarkers.csv")
write.csv(so_leidenmarkers, file.path(table_dir, outfile), row.names = FALSE)
