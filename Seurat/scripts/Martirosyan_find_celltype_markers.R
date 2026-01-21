#!/bin/env Rscript

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

#---------- Loading in Seurat Object and Grouping Cell Types -----------

#loading in so
so <- readRDS("Martirosyan_2024_integrated.rds")

#changing idents to celltype so then can group 
Idents(so) <-  so$celltype

#if needed grouping cell types into condensed groups
so$celltype_condensed <- plyr::mapvalues(Idents(so),
  from = c("Oligo", "Astro", "Exc", "DopaN", "Micro", "OPC", "Inh", "OtherNeuron", "Unk", "Pericyte_Endothelial_VLMC", "Tcell" ),
  to = c("Oligo", "Astro", "Neurons", "Neurons", "Micro", "OPCs", "Neurons", "Neurons", "T_Cells", "Vascular", "T_Cells"))

#"celltype_condensed" is the metadata column with the condensed cell type groups
Idents(so) <-  so$celltype_condensed

#reordering celltypes
Idents(so) <- factor(so$celltype_condensed, levels = c("Neurons",
                                                       "Astro",
                                                       "Micro",
                                                       "Oligo",
                                                       "OPCs",
                                                       "Vascular",
                                                       "T_Cells"))

#---------- Running FindMarkers Function -----------

# For each cell type of interest, run FindMarkers() to compare against other cell types
#creating a variable for all broad cell types within seurat object
broad_celltypes <- unique(so$celltype_condensed)
sighits_list <- list()
results_list <- list()

cat("Running FindMarkers based on cell type...\n")
for(cell_type in broad_celltypes) {
  cat("Processing cell type:", cell_type, "\n")
  other_types <- setdiff(broad_celltypes, cell_type)
  
  #Run analysis using likelihood ratio test with a logistic regression model
  diff = FindMarkers(object = so, ident.1 = cell_type, ident.2 = other_types, test.use = "LR")
  diff$broad_celltype <- cell_type
  diff$gene <- rownames(diff)
  diff <- diff %>%
    mutate(regulation = case_when(
      avg_log2FC > 0 & p_val_adj < 0.05 ~ "up",
      avg_log2FC < 0 & p_val_adj < 0.05 ~ "down",
      TRUE ~ "none"))
  diff <- diff[, c("gene", setdiff(names(diff), "gene"))]
  results_list[[cell_type]] <- diff
  outfile = paste0("differential_analysis/",cell_type, "_cellvcell_LR_results.txt")
  write.table(diff, outfile, sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
  sig = diff[!is.na(diff$p_val_adj) & diff$p_val_adj < 0.05,]
  outfile = paste0("differential_analysis/",cell_type, "_cellvcell_LR_sighits.txt")
  write.table(sig, outfile, sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
  sig$gene <- rownames(sig)
  rownames(diff) <- NULL
  rownames(sig) <- NULL
  results_list[[cell_type]] <- diff
  sighits_list[[cell_type]] <- sig
}

sighits <- dplyr::bind_rows(sighits_list)
results <- dplyr::bind_rows(results_list)
sighits_outfile <- paste0("differential_analysis/Martirosyan_cellvcell_LR_sighits.txt")
results_outfile <- paste0("differential_analysis/Martirosyan_cellvcell_LR_results.txt")
write.table(sighits, file.path(sighits_outfile), sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(results, file.path(results_outfile), sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
