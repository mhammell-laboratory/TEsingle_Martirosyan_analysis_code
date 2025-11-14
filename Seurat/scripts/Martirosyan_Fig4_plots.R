#!/bin/env Rscript

#loading in necessary packages to run Seurat, analyses, and plot graphs
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
library(ggvenn)
library(ComplexUpset)
#setting Memory limit and bitmap type
options(future.globals.maxSize = 600 * 1024^3)

rds_dir  <- "subpopulation/rds"
de_dir <- "subpopulation/diffExp"
plots_dir <- "figures"

## Figure 4: Martirosyan Dataset Overview

### Figure 4B: Martirosyan UMAP

#loading in so
so <- readRDS("Martirosyan_2024_integrated.rds"))

#if needed grouping cell types into condensed groups
so$celltype_condensed <- plyr::mapvalues(
  Idents(so), from = c("Oligo", "Astro", "Exc", "DopaN", "Micro", "OPC", "Inh", "OtherNeuron", "Unk", "Pericyte_Endothelial_VLMC", "Tcell"), 
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

#setting cell type colors
celltype_cols <- c("Neurons" = "#DF939A",
                  "Astro" = "#75B2C7",
                  "Micro" = "#726FB9",
                  "Oligo" = "#55A35D",
                  "OPCs" = "#D59F66",
                  "Vascular" = "#E3C775",
                  "T_Cells"= "#73736F")

#parameters for displaying UMAP
so_UMAP <- DimPlot(so, cols = celltype_cols, reduction = "umap.cca", shuffle = TRUE, label = TRUE, label.size = 5) +
  labs(color = "Cell Types") +
  theme(legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14, face = "plain"))
#so_UMAP

#saving UMAP plot
outplot <- "Fig4B.png"
ggsave(filename = file.path(plots_dir, outplot),
         plot = so_UMAP,
         width = 10,
         height = 8,
         dpi = 300)


### Figure 4C: Martirosyan Cell Type Marker Dot Plot

#setting expression colors
expression_cols <- c("snow2","#39738D")

#setting celltype specific marker genes used in Martirosyan paper
celltypes_genes <- c("SNAP25", "SYT1", "GFAP", "AQP4", "P2RY12", "CD74", "MBP", "MOBP", "PDGFRA", "VCAN", "DCN", "FLT1", "CD2", "THEMIS")

#rearranging order of celltypes to be displayed
so$celltype_condensed <- factor(so$celltype_condensed,
                                levels = c("Neurons",
                                           "Astro",
                                           "Micro",
                                           "Oligo",
                                           "OPCs",
                                           "Vascular",
                                           "T_Cells"))

#parameters used to make dotplot
celltype_dotplot <- DotPlot(object = so, features = celltypes_genes, group.by = "celltype_condensed", scale.by = "size", dot.scale = 10, cols = expression_cols) +
  labs(x = "Genes", y = "Cell Types") +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(face = "plain", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "plain", size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14, face = "plain"),
        panel.spacing = unit(1, "lines"),
        aspect.ratio = 0.4)
#celltype_dotplot

outplot <- "Fig4C.png"
ggsave((file.path(plots_dir, outplot)),
       plot = celltype_dotplot,
       width = 10,
       height = 4,
       dpi = 300)

### Figure 4D: Martirosyan Cell Proportions Horizontal Bar Graph

#creating data needed to make proportion bar graph
#creating a variable that holds all control patientIDs
control_ids <- unique(so$patientID[so$diagnosis == "Control"])
#creating a variable that holds all PD patientIDs
parkinsons_ids <- unique(so$patientID[so$diagnosis == "Parkinsons"])

#creating a logical index for control cells
control_cells <- so$patientID %in% control_ids
#creating a count table that only contain control cells for each cell type
control_cell_counts_table <- table(so$patientID[control_cells], so$celltype_condensed[control_cells])
#finding the mean number of cells for each control cell type
control_celltype_mean <- colMeans(control_cell_counts_table)

#creating a logical index for PD cells
pd_cells <- so$patientID %in% parkinsons_ids
#creating a count table that only contain PD cells for each cell type
pd_cell_counts_table <- table(so$patientID[pd_cells], so$celltype_condensed[pd_cells])
#finding the mean number of cells for each PD cell type
pd_celltype_mean <- colMeans(pd_cell_counts_table)

#creating long dataframe
mean_prop_df <- data.frame(celltype = names(control_celltype_mean), 
                           Control = control_celltype_mean, Parkinsons = pd_celltype_mean) %>%
  pivot_longer(cols = c(Control, Parkinsons), names_to = "diagnosis", values_to = "mean_count")

#calculating total mean count across diagnoses for each cell type
mean_prop_df <- mean_prop_df %>%
  group_by(celltype) %>%
  mutate(celltype_total = sum(mean_count)) %>%
  ungroup()
#calculating percentage contribution of each diagnosis to each cell type's total mean count
mean_prop_df <- mean_prop_df %>%
  group_by(diagnosis) %>%
  mutate(percent = mean_count / celltype_total * 100) %>%
  ungroup()

#reordering celltype
mean_prop_df$celltype <- factor(mean_prop_df$celltype, levels = rev(c("Neurons",
                                                                        "Astro",
                                                                        "Micro",
                                                                        "Oligo",
                                                                        "OPCs",
                                                                        "Vascular",
                                                                        "T_Cells")))

#parameters uesed to create the horizontal bar plot showing mean count for each cell type and diagnosis contribution to each cell type
celltype_mean_composition_plot <- ggplot(mean_prop_df, aes(x = percent, y = celltype, fill = diagnosis)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "#73736F") +
  geom_text(aes(label = round(mean_count, 0), 
                x = ifelse(diagnosis == "Parkinsons", 0, 100), 
                hjust = ifelse(diagnosis == "Parkinsons", -0.2, 1.2)),
            color = "white", size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Control" = "#75B2C7", "Parkinsons" = "#CF3B3D")) +
  labs(title = "Cell Type Composition by Diagnosis (Mean per Sample)",
       x = "\nOrigin %",
       y = "Cell Type",
       fill = "") +
  scale_x_continuous(breaks = c(0, 50, 100), labels = c("0", "50", "100")) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "bold"))
#celltype_mean_composition_plot

outplot <- "Fig4D.png"
ggsave((file.path(plots_dir, outplot)),
  plot = celltype_mean_composition_plot,
  width = 10,
  height = 8,
  dpi = 300)

### Figure 4E: Martirosyan Cell Type Pie Chart

#creating df for pie chart using total cell type counts
celltype_counts_df <- as.data.frame(table(so$celltype_condensed))
#renaming columns in the df
colnames(celltype_counts_df) <- c("celltype", "total_counts")

#setting cell type colors
celltype_cols <- c("Neurons" = "#DF939A",
                  "Astro" = "#75B2C7",
                  "Micro" = "#726FB9",
                  "Oligo" = "#55A35D",
                  "OPCs" = "#D59F66",
                  "Vascular" = "#E3C775",
                  "T_Cells"= "#73736F")

#reordering cell types
celltype_counts_df$celltype <- factor(celltype_counts_df$celltype, 
                                          levels = c("Neurons",
                                                     "Astro",
                                                     "Micro",
                                                     "Oligo",
                                                     "OPCs",
                                                     "Vascular",
                                                     "T_Cells"))

#adding total cell count to df from the dataset
celltype_counts_df$totalcells <- ncol(so)

#calculating the percent of each cell type within the dataset
celltype_counts_df$percent_celltype <- round((celltype_counts_df$total_counts / celltype_counts_df$totalcells) *100, 1)
#adding labels to be used of the legend
celltype_counts_df <- celltype_counts_df %>%
  mutate(legend_label = paste0(celltype, " (", percent_celltype, "%)"))
#reordering df so that it matches factor order
celltype_counts_df <- celltype_counts_df %>%
  arrange(celltype)

#parameters used to make cell type composition pie chart
celltype_composition_pie <- ggplot(celltype_counts_df, 
                                   aes(x = "", y = percent_celltype, fill = celltype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = celltype_cols, 
                   labels = celltype_counts_df$legend_label,
                   name = "Cell Types", drop = FALSE) +
  labs(fill = "Cell Types") +
  theme_void() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 14, face = "plain"),
        legend.title = element_text(size = 18, face = "bold"))
#celltype_composition_pie

outplot <- "Fig4E.png"
ggsave((file.path(plots_dir, outplot)),
  plot = celltype_composition_pie,
  width = 10,
  height = 8,
  dpi = 300)

### Figure 4F: Martirosyan DE Genes Comparison Venn Diagram

#loading in martirosyan DEGs table
martirosyan_degs <- read.csv("path/to/martirosyan/DEG/table")
#creating a variable for all DEGs found in martirosyan et al
martirosyan_degs_genes <- unique(martirosyan_degs$gene)

#loading in TEsingle DEGs table
TEsingle_degs <- read.csv("path/to/TEsingle/DEG/table")
#creating a variable for all DEGs found in martirosyan et al
TEsingle_genes <- unique(TEsingle_degs$gene)

#parameters used to make venn diagram
all_degs <- ggvenn(list(Martirosyan = martirosyan_degs_genes, TEsingle = TEsingle_genes),
  fill_color = c("#55A35D", "#D59F66"),
  stroke_size = 0.5, set_name_size = 8, auto_scale = FALSE, show_percentage = FALSE) +
  ggtitle("Shared and Unique DEGs in all Cell Types\n") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA))
#all_degs

outplot = "Fig4F.png"
ggsave(filename = file.path(plots_dir, outplot),
         plot = all_degs,
         width = 10,
         height = 8,
         dpi = 300)


