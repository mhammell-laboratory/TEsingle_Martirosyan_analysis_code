#---------- Setting Up R Script -----------

#loading required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(scales)
library(tibble)
library(forcats)
library(tidyr)
#setting Memory limit and bitmap type
options(future.globals.maxSize = 600 * 1024^3)

#setting input & output folders
rds_dir  <- "subpopulation/rds"
table_dir <- "subpopulation/tables"
plots_dir <- "subpopulation/plots"
de_dir <- "subpopulations/diffExp"
dataset <- "Martirosyan_Neuron"

#---------- Loading in Seurat Object and Leiden Cluster Markers -----------

#loading in Integrated subpopulation Seurat Object RDS
cat("Loading Seurat Object RDS...\n")
#change "BroadCellType" to the name of the integrated broad cell type Seurat Object RDS file 
#that you want to use for subpopulation identification
infile = paste0(dataset, "_integrated.rds")
so <- readRDS(file.path(working_dir, infile))
#setting identity of SO to leiden_subpopulation
Idents(so) = "leiden_subpopulation"

#---------- Renaming Leiden Clusters as Predicted Subpopulations -----------

#defining new metadata column with renamed leiden clusters based on subpopulation identification
#example code: renaming leiden clusters to subpopulation cell types for neuronal broad cell type
so$cluster_subpopulation <- plyr::mapvalues(
  Idents(so),
  from = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
  to = c("DopaN", "Other_Neurons", "Doublet (ODC)", "DopaN", "Excitatory", "Inhibitory", "Doublet (ODC)", "Doublet (ODC)", "Inhibitory", "Doublet (Astro)", "DopaN", "Inhibitory", "Excitatory"))

#saving Seurat object with predicted subpopulation clusters
outfile <- paste0(dataset,"_subpopulations.rds")
saveRDS(so, file = file.path(rds_dir, outfile))

#---------- Defining Categories For DE Analysis -----------

#defining category for DE analysis
cat("Setting categories...\n")
#"cluster_subpopulation" is the metadata column name of the subpopulation cell type
category <- paste(as.vector(so$diagnosis), as.vector(so$cluster_subpopulation), sep="_")
so[["category"]] <- category
Idents(so) <- "category"

#---------- Running PrepSCTFindMarkers Function -----------

#need to run PrepSCTFindMarkers before running FindMarkers if used SCTransform normalization method
cat("Running PrepSCTFindMarkers Function...\n")
so <- PrepSCTFindMarkers(so)

#---------- Running FindMarkers Function -----------

# For each cell type of interest, run FindMarkers() to compare disease to control cells
cat("Running FindMarkers based on cell type and diagnosis...\n")
#subpopulation cell types for each broad cell type:
#neuron cell types: DopaN, Inhibitory, Excitatory, Other_Neurons
#astrocyte cell types: Homeostatic, Intermediate, Reactive, Disease_associated, Inflammatory, Metallothionein
#microglia cell types: Homeostatic, Cytokine_response, Interferon_response, Antigen_presenting, Disease_associated
#oligodendrocyte cell types: Resting, Myelinating, Immature, Stress_response, Ox_stress

# extract list of cell subpopulations
subpopulation_list = levels(factor(so$cluster_subpopulation))
subpopulation_list = subpopulation_list[grep("Doublet",subpopulation_list, invert=TRUE)]
for(cell in subpopulation_list) {
  disease <- paste0("Parkinsons_",cell)
  control <- paste0("Control_",cell)
  #run analysis using likelihood ratio test with a logistic regression model
  diff <- FindMarkers(object = so, ident.1 = disease, ident.2 = control, test.use = "LR")
  #in outfile names, change "broadcelltype" to the name of the broad cell type you are analyzing
  outfile <- paste0(dataset, "_", cell,"_PDvCtrl_LR_results.txt")
  write.table(diff, file.path(de_dir, outfile)), sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
  sig <- diff[!is.na(diff$p_val_adj) & diff$p_val_adj < 0.05,]
  sig <- sig %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  mutate(regulation = case_when(
    avg_log2FC > 0.01 & p_val_adj < 0.05 ~ "up",
    avg_log2FC < -0.01 & p_val_adj < 0.05 ~ "down",
    TRUE ~ "none"
    ))
  sig$cell_type <- cell
  outfile <- paste0(dataset,"_", cell, "_PDvCtrl_LR_sighits.txt")
  write.table(sig, file.path(de_dir, outfile)), sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
}

#saving RDS with DE results
cat("Saving Seurat Object...\n")
outfile = paste0(dataset,"_DE.rds")
saveRDS(so, file = file.path(rds_dir, outfile))

==========================

#loading in csv file of FindAllMarker() leiden cluster markers csv
infile = paste0(dataset,"_LeidenMarkers.csv")
leiden_clustermarkers <- read.csv(file = file.path(table_dir, infile), header = TRUE)

#creating table with only the upregulated cluster markers
leiden_clustermarkers_up <- leiden_clustermarkers %>%
  filter(regulation == "up") %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))

#---------- Creating General and Subpopulation Marker Heatmaps -----------
#this is the first pass to annotating the leiden subpopulation clusters
#the general marker heatmap helps with identifying potential doublet clusters
#the subpopulation marker heatmap helps further define specific subpopulation of the broad cell types found in the dataset

#creating a heatmap of general cell type markers
#loading in general cell type markers file
general_markers <- read.table(file.path(table_dir, "general_markers_file.txt"), sep="\t", row.names = NULL, header = F)
#creating a heatmap of general cell type markers based on features found in the SO
if (!is.null(general_markers)) {
  #creating a variable of available general markers in the so and printing how many are found
  available_general_markers <- rownames(general_markers)[rownames(general_markers) %in% rownames(so)]
  cat("Using", length(available_general_markers), "out of", nrow(general_markers), "general markers found in dataset\n")
  #creating heatmap of general markers based on available markers found in so
  if (length(available_general_markers) > 0) {
    general_marker_heatmap <- DoHeatmap(object = subset(so, downsample = 100), features = available_general_markers, size = 3) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
    theme(axis.text.y = element_text(size = 8))
  }
}
#saving general markers heatmap
ggsave(filename = file.path(plots_dir, "Martirosyan_BroadCellType_general_heatmap_YYYYMMDD.png"),
         plot = general_marker_heatmap,
         width = 10,
         height = 8,
         dpi = 300)

#creating a heatmap of subpopulation cell type markers based on previously published datasets
#loading in subpopulation cell type markers file
subpopulation_markers <- read.table(file.path(table_dir, "broad_celltype_subpopulation_markers_file.txt"), sep="\t", row.names = NULL, header = F)
#creating a heatmap of subpopulation cell type markers based on features found in the SO
if (!is.null(subpopulation_markers)) {
  #creating a variable of available subpopulation markers in the so and printing how many are found
  available_subpopulation_markers <- rownames(subpopulation_markers)[rownames(subpopulation_markers) %in% rownames(so)]
  cat("Using", length(available_subpopulation_markers), "out of", nrow(subpopulation_markers), "subpopulation markers found in dataset\n")
#creating heatmap of subpopulation markers based on available markers found in so
  if (length(available_subpopulation_markers) > 0) {
    subpopulation_marker_heatmap <- DoHeatmap(object = subset(so, downsample = 100), features = available_subpopulation_markers, size = 3) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
    theme(axis.text.y = element_text(size = 8))
  }
}
#saving subpopulation markers heatmap
ggsave(filename = file.path(plots_dir, "Martirosyan_BroadCellType_subpopulation_heatmap_YYYYMMDD.png"),
         plot = subpopulation_marker_heatmap,
         width = 10,
         height = 8,
         dpi = 300)

#---------- Using Gene Lists to Identify Leiden Subpopulation Clusters -----------

#this is the second pass to annotating the leiden subpopulation clusters
#using published gene lists, I looked to see which clusters had markers that overlapped with these gene lists

#example code: using dopaminergic neuronal genes to identify which leiden clusters are dopaminergic neurons
leiden_clustermarkers_up[leiden_clustermarkers_up$gene %in% c("TH", "SLC6A3", "KCNJ6", "NR4A2", "LMX1B", "LMX1A", "FOLR1", "FOXA2", "DBH", "SLC6A2", "PITX3", "SMAD3", "NEUROD6", "SLC18A2", "DDC", "SLC18A3"),]
#listing leiden cluster numbers that have these markers enriched
#Clusters: 1,2,4,11

#example code: using GABAergic neuronal genes to identify which leiden clusters are GABAergic neurons
leiden_clustermarkers_up[leiden_clustermarkers_up$gene %in% c("GABBR1", "GABBR2", "GAD1", "GAD2", "SLC6A1", "SLC32A1", "DLX1"),]
#listing leiden cluster numbers that have these markers enriched
#Clusters: 6,9,12

#loading in specific subpopulation gene lists and checking for overlap with leiden cluster markers
#example code: using a Ramcao reactive astrocytes genes list to identify which leiden clusters are reactive ramcao astrocytes
astro_ramcao <- read.table(file.path(table_dir_dir, "Astro_ramcao.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#identifying which genes overlap with upregulated leiden cluster markers
intersect(astro_ramcao, leiden_clustermarkers_up$gene)
#listing genes that are found to be upregulated in leiden cluster markers table
#"EMP1", "B3GNT5", "HSPB1", "CD14", "CD44", "HMOX1", "OSMR", "GAP43", "SBNO2", "GCNT2", "CP", "GADD45B", "SERPINA3", "AKAP12", "GFAP"
#using these genes to identify which leiden clusters are ramcao reactive astrocytes
leiden_clustermarkers_up[leiden_clustermarkers_up$gene %in% c("EMP1", "B3GNT5", "HSPB1", "CD14", "CD44", "HMOX1", "OSMR", "GAP43", "SBNO2", "GCNT2", "CP", "GADD45B", "SERPINA3", "AKAP12", "GFAP"),]
#listing leiden cluster numbers that have these markers enriched
#Clusters: 8, 9, 10, 13, 15

#---------- Using FeaturePlots, Violin Plots, and Dot Plots to Identify Leiden Subpopulation Clusters -----------

#this is the third pass to annotating the leiden subpopulation clusters
#using Seurat's FeaturePlots and Violin plots to visualize expression of specific genes and seeing if they are
#enriched in specific leiden clusters to help identify subpopulation of cell types

#FeaturePlots
#example code: using TH marker to identify which leiden clusters are dopaminergic neurons
TH_feature <- FeaturePlot(so, features = 'TH', reduction = "umap.cca", order = TRUE, label = TRUE)
#example code: using TH and SLC6A3 (DAT1) markers to identify which leiden clusters are dopaminergic neurons
#parameter blend = TRUE allows for visualization of co-expression of two genes
TH_SLC6A3 <- FeaturePlot(so, features = c("TH", "SLC6A3"), reduction = "umap.cca", label = TRUE, blend = TRUE, 
cols = c("lightgrey", "green", "red"), ncol=4)

#Violin Plots
#example code: using violin plots to visualize expression of dopaminergic neuronal markers in leiden clusters
DopaN_VlnPlot <- wrap_plots(lapply(c("TH", "SLC6A3", "LMX1B", "LMX1A"), function(g) VlnPlot(so, features = g) + NoLegend()), ncol = 2) +
  plot_annotation(title = "DopaN Markers", theme = theme(plot.title = element_text(hjust = 0.5)))
#example code: using violin plots to visualize expression of GABAergic neuronal markers in leiden clusters
GABA_VlnPlot <- wrap_plots(lapply(c("GAD1", "GAD2", "SLC6A1", "GABBR1"), function(g) VlnPlot(so, features = g) + NoLegend()), ncol = 2) +
  plot_annotation(title = "GABAergic Neuronal Markers", theme = theme(plot.title = element_text(hjust = 0.5)))

#DotPlots
#example code: creating dotplots to visualize average expression and percent expression
#of dopaminergic neuronal markers in leiden clusters
DopaN_dotplot <- DotPlot(so, cols = c("cornflowerblue", "red"), features = c("TH", "SLC6A3", "LMX1B", "LMX1A")) +
RotatedAxis()
#example code: dotplot visualization of average expression and percent expression of GABAergic neuronal markers in leiden clusters
GABA_dotplot <- DotPlot(so, cols = c("cornflowerblue", "red"), features = c("GAD1", "GAD2", "SLC6A1", "GABBR1")) +
RotatedAxis()

