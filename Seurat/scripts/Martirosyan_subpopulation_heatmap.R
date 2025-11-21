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

if(length(args) < 1){
    stop("Usage: Rscript Martirosyan_subpopulation_integration.R [Neuron/Astro/Micro/Oligo]", call.=FALSE)
}
celltype <- args[1]


#setting input & output folders
rds_dir  <- "subpopulation/rds"
plots_dir <- "subpopulation/plots"
marker_dir <- "subpopulations/markers"
dataset <- paste0("Martirosyan_",celltype)

#---------- Loading in Seurat Object and Leiden Cluster Markers -----------

#loading in Integrated subpopulation Seurat Object RDS
cat("Loading Seurat Object RDS...\n")
#loading in the subpopulation integrated dataset
infile <- paste0(dataset, "_integrated.rds")
so <- readRDS(file.path(rds_dir, infile))
Idents(so) <- "leiden_subpopulation"

#---------- Creating General and Subpopulation Marker Heatmaps -----------

#the general marker heatmap helps with identifying potential doublet clusters
#creating a heatmap of general cell type markers
#loading in general cell type markers file
general_markers <- read.table(file.path(marker_dir,"official_markers.txt"), sep="\t", row.names = 1, header = F)
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
outplot <- paste0(dataset,"_general_markers_heatmap.png")
ggsave(filename = file.path(plots_dir, outplot),
       plot = general_marker_heatmap,
       width = 10,
       height = 8,
       dpi = 300)

#the subpopulation marker heatmap helps further define specific subpopulations of the broad cell types found in the dataset
#creating a heatmap of subpopulation cell type markers based on previously published datasets
#loading in subpopulation cell type markers file
subpopulation_markers <- read.table(file.path(marker_dir, "astrocyte_markers.txt"), sep="\t", row.names = 1, header = F)
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
outplot <- paste0(dataset,"_subpopulation_markers_heatmap.png")
ggsave(filename = file.path(plots_dir, outplot),
       plot = subpopulation_marker_heatmap,
       width = 10,
       height = 8,
       dpi = 300)
