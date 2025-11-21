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
de_dir <- "subpopulations/diffExp"
dataset <- "Martirosyan_Micro"

#---------- Loading in Seurat Object and Leiden Cluster Markers -----------

#loading in Integrated subpopulation Seurat Object RDS
cat("Loading Seurat Object RDS...\n")
#loading in the subpopulation integrated dataset
infile <- paste0(dataset, "_integrated.rds")
so <- readRDS(file.path(rds_dir, infile))
Idents(so) <- "leiden_subpopulation"

#---------- Renaming Leiden Clusters as Predicted Subpopulations -----------

#defining new metadata column with renamed leiden clusters based on subpopulation identification
so$cluster_subpopulation <- plyr::mapvalues(
  Idents(so),
  from = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
  to = c("Homeostatic", "Doublet (Neurons)", "Cytokine_response", "Disease_associated", "Interferon_response", "Homeostatic", "Cytokine_response", "Doublet (ODC)", "Antigen_presenting", "Interferon_response", "Antigen_presenting", "Doublet (Macrophage)", "Homeostatic", "Doublet (Astro)"))

#saving Seurat object with predicted subpopulation clusters
Idents(so) <- "cluster_subpopulation"
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
# extract list of cell subpopulations
subpopulation_list <- levels(factor(so$cluster_subpopulation))
subpopulation_list <- subpopulation_list[grep("Doublet",subpopulation_list, invert=TRUE)]
subpopulation_sighits_list <- list()

for(i in 1:length(subpopulation_list)) {
  cell <- subpopulation_list[i]
  disease <- paste0("Parkinsons_",cell)
  control <- paste0("Control_",cell)
  cat("Running Differential Analysis on", cell, "subpopulation")
  #run analysis using likelihood ratio test with a logistic regression model
  diff <- FindMarkers(object = so, ident.1 = disease, ident.2 = control, test.use = "LR")
  outfile <- paste0(dataset, "_", cell,"_PDvCtrl_LR_results.txt")
  write.table(diff, file.path(de_dir, outfile), sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
  sig <- diff[!is.na(diff$p_val_adj) & diff$p_val_adj < 0.05,]
  sig <- sig %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    mutate(regulation = case_when(
      avg_log2FC > 0.01 & p_val_adj < 0.05 ~ "up",
      avg_log2FC < -0.01 & p_val_adj < 0.05 ~ "down",
      TRUE ~ "none"))
  sig$cell_type <- cell
  outfile <- paste0(dataset,"_", cell, "_PDvCtrl_LR_sighits.txt")
  write.table(sig, file.path(de_dir, outfile), sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
  sig$gene <- rownames(sig)
  rownames(sig) <- NULL
  sig <- sig[,c(ncol(sig),2:ncol(sig)-1)]
  subpopulation_sighits_list[[i]] <- sig
}

subpopulation_sighits <- dplyr::bind_rows(subpopulation_sighits_list)
outfile <- paste0(dataset,"_subpopulation_PDvCtrl_LR_sighits.txt")
write.table(subpopulation_sighits, file.path(de_dir, outfile), sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)

#saving RDS with DE results
cat("Saving Seurat Object...\n")
outfile <- paste0(dataset,"_DE.rds")
saveRDS(so, file = file.path(rds_dir, outfile))
