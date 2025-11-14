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


## Supplementary Figure: ODCs Subpopulation

### Supp Figure A: ODCs Subpopulation UMAP


#loading in microglia subpopulation so
odc_so <- readRDS(file.path(rds_dir, "Martirosyan_Oligo_subpopulations.rds"))

#setting subpopulation colors
celltype_cols <- c("Resting" = "#346b50",
                   "Myelinating" = "#408759",
                   "Immature" = "#55A35D",
                   "Stress_associated" = "#80b285",
                   "Ox_stress" = "#A9C0AC")

#parameters used to create UMAP
odc_UMAP <- DimPlot(odc_so, cols = celltype_cols, reduction = "umap.cca", shuffle = TRUE, label = TRUE, label.size = 5) +
  labs(color = "Cell Types") +
  theme(legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14, face = "plain"))
#odc_UMAP

outplot <- "SuppFigXA.png"
ggsave(filename = file.path(plots_dir, outplot),
         plot = odc_UMAP,
         width = 10,
         height = 8,
         dpi = 300)

### Supp Figure B: ODCs Subpopulation Marker Violin Plots

#parameters used to create Violin Plots
odc_VlnPlot <- wrap_plots(lapply(c("OPALIN", "OLIG1", "DNAJC3", "HMOX2", "MOG"), function(g) 
  VlnPlot(odc_so, features = g, pt.size = 0, cols = celltype_cols, y.max = 5) + 
  NoLegend() + 
  ggtitle("") +
  ylab(g) +
  xlab("") +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, 
                                    vjust = 0.5, face = "bold", size = 18))), ncol = 3)
#odc_VlnPlot

outplot <- "SuppFigXB.png"
ggsave(filename = file.path(plots_dir, outplot),
         plot = odc_VlnPlot,
         width = 8,
         height = 8,
         dpi = 300)


### Supp Figure C: ODCs Subpopulation Horizontal DE TE Loci Bar Plot


#loading in subpopulation DE Genes and TE table
master_DE_odc_df  <- read.table(file.path(de_dir,"Martirosyan_Oligo_subpopulation_PDvCtrl_LR_sighits.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#adding TE_element information column
#this can be extracted from the TE name
master_DE_odc_df <- master_DE_odc_df %>%
  mutate(TE_element = case_when(
    grepl(":", gene) ~ sapply(strsplit(master_DE_odc_df$gene, ":"), function(x) x[4]),
    TRUE ~ NA_character_))

#adding TE_family information column
master_DE_odc_df <- master_DE_odc_df %>%
  mutate(TE_family = case_when(
    grepl(":", gene) ~ sapply(strsplit(master_DE_odc_df$gene, ":"), function(x) x[3]),
    TRUE ~ NA_character_))

#adding TE_subfamily information column
master_DE_odc_df <- master_DE_odc_df %>%
  mutate(TE_subfamily = case_when(
    grepl(":", gene) ~ sapply(strsplit(master_DE_odc_df$gene, ":"), function(x) x[2]),
    TRUE ~ NA_character_))

#grouping TE_element into broad TE groups
master_DE_odc_df <- master_DE_odc_df %>%
  mutate(
    TE_group = case_when(
      TE_element == "SINE" ~ "SINEs",
      TE_element == "Retroposon" & grepl("^SVA", TE_family) ~ "SINEs",
      TE_element == "LINE"  ~ "LINEs",
      TE_element == "LTR" & grepl("^ERV", TE_family) ~ "ERVs",
      TE_element == "LTR" & grepl("Gypsy", TE_family) ~ "ERVs",
      TE_element == "DNA"  ~ "Other",
      TE_element == "Satellite"  ~ "Other",
      TE_element == "Unknown"  ~ "Other",
      is.na(TE_element) ~ NA_character_, 
      TRUE ~ NA_character_))

#setting regulation direction colors
regulation_cols <- c("down" = "#346b50", "up" = "#A9C0AC")
#setting regulation lables for legend
regulation_label <- c("up" = "Up", "down" = "Down")

#creating df for horizontal barplot
#creating df that only contain DE TEs
master_DE_TE_odc_df <- master_DE_odc_df %>%
  filter(grepl(":", master_DE_odc_df$gene))

#creating df that counts how many TE loci are DE within each subpopulation
TE_family_counts <- master_DE_TE_odc_df %>%
  group_by(cell_type, regulation, TE_group) %>%
  summarise(count = n(), .groups = "drop")
#counting number of DE TE loci and flipping downregulated loci for plotting
TE_family_counts <- TE_family_counts %>%
  group_by(cell_type, regulation) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(count = ifelse(regulation == "down", -count, count))

#reordering regulation
TE_family_counts$regulation <- factor(TE_family_counts$regulation,
                                     levels = c("up", "down"))
#reordering subpopulation
TE_family_counts$cell_type <- factor(TE_family_counts$cell_type,
                                     levels = rev(c("Resting", 
                                                "Myelinating", 
                                                "Immature",
                                                "Stress_associated",
                                                "Ox_stress")))

#parameters used to make the horizontal bar plot
total_TE_DE_Counts <- ggplot(TE_family_counts, aes(x = count, y = cell_type, fill = regulation)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = abs(count)),
      hjust = ifelse(TE_family_counts$count < 0, 1.1, -0.1),
      size = 5, fontface = "bold") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#73736F") +
    scale_fill_manual(values = regulation_cols, labels = regulation_label) +
    labs(x = "\nTotal TE Loci Counts",
      y = "Cell Type Subpopulations\n",
      fill = "Direction") +
    theme_classic() +
    theme(axis.text.x = element_text(face = "plain", size = 12), color = "black",
    axis.text.y = element_text(face = "plain", size = 14, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14, face = "plain"),
    plot.title = element_text(hjust = 0.5))
#total_TE_DE_Counts

outplot <- "SuppFigXC.png"
ggsave(filename = file.path(plots_dir, outplot),
       plot = total_TE_DE_Counts,
       width = 12,
       height = 6,
       dpi = 300)

### Supp Figure D: ODCs Subpopulation Upset Plot

#loading package used to make the upset plot
#not compatible with ggplot2 version 4.0.0

#creating variables to list TE loci specific to the cell type subpopulation
resting_TE_loci <- unique(master_DE_TE_odc_df$gene[master_DE_TE_odc_df$cell_type == "Resting"])
myelinating_TE_loci <- unique(master_DE_TE_odc_df$gene[master_DE_TE_odc_df$cell_type == "Myelinating"])
immature_TE_loci <- unique(master_DE_TE_odc_df$gene[master_DE_TE_odc_df$cell_type == "Immature"])
stress_TE_loci <- unique(master_DE_TE_odc_df$gene[master_DE_TE_odc_df$cell_type == "Stress_associated"])
ox_TE_loci <- unique(master_DE_TE_odc_df$gene[master_DE_TE_odc_df$cell_type == "Ox_stress"])

#creating named list to be used for upset plot
te_lists <- list(Resting = resting_TE_loci,
                 Myelinating = myelinating_TE_loci,
                 Immature = immature_TE_loci,
                 Stress_associated = stress_TE_loci,
                 Ox_stress = ox_TE_loci)
 
#converting list to the format ComplexUpset needs
all_TEs <- unique(unlist(te_lists))
upset_data <- sapply(te_lists, function(x) as.numeric(all_TEs %in% x))
rownames(upset_data) <- all_TEs
upset_data <- as.data.frame(upset_data)

#parameters used to make upset plot
upset_plot <- upset(upset_data, colnames(upset_data), name = 'TE Loci', width_ratio = 0.3,
                    stripes='white', 
                    base_annotations = list('Intersection size' = 
                                              intersection_size(bar_number_threshold = 1)),
                    set_sizes = FALSE, 
                    themes=upset_default_themes(text=element_text(color='black', 
                                                                  face = "bold", 
                                                                  size = 18)))
#upset_plot

outplot <- "SuppXD.png"
ggsave(filename = file.path(plots_dir, outplot),
       plot = upset_plot,
       width = 10,
       height = 6,
       dpi = 300)

### Supp Figure XE: ODCs Subpopulation DE TE Loci Volcano Plot

#creating variable for all the cell type subpopulation 
celltype_subpopulation <- unique(master_DE_odc_df$cell_type)

#creating an empty plot list to hold all subpopulation volcano plots
volcano_plots <- list()
#creating a loop to make each subpopulation volcano plot
for (cell in celltype_subpopulation) {
  #creating df for one subpopulation
  volcanoplot_df <- master_DE_odc_df %>%
    filter(cell_type == cell) %>%
    mutate(regulation = case_when(avg_log2FC > 0 ~ "up",
                                  avg_log2FC < 0 ~ "down",
                                  TRUE ~ "neutral"),
           #creating a variable to highlight TE loci, which contain ":" in the name
           is_highlighted <- grepl(":", gene),
           #creating a color group to ....
           color_group <- case_when(
             avg_log2FC > 1 & regulation == "up" ~ "up",
             avg_log2FC < -1 & regulation == "down" ~ "down",
             TRUE ~ "other")) %>%
    mutate(color_group = factor(color_group, levels = c("up", "down", "other"))) %>%
    arrange(regulation)
  
  #creating separate df for highlighted and non-highlighted points
  non_highlighted <- volcanoplot_df %>%
    filter(!is_highlighted)
  
  #creating separate dfs for each TE element group to color them differently
  sines <- volcanoplot_df %>% filter(is_highlighted & TE_element == "SINEs")
  lines <- volcanoplot_df %>% filter(is_highlighted & TE_element == "LINEs")
  ervs <- volcanoplot_df %>% filter(is_highlighted & TE_element == "ERVs")
  other_tes <- volcanoplot_df %>% filter(is_highlighted & TE_element == "Other")
  
  #parameters used to make the volcano plots
  p <- ggplot(volcanoplot_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_vline(xintercept = c(-1, 1), col = "dimgrey", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "dimgrey", linetype = 'dashed') +
    #plotting non-highlighted points first as the background layer with a smaller dot size
    geom_point(data = non_highlighted, color = "gray", size = 1.5) +
    #plotting highlighted points by TE element group with a larger dot size
    geom_point(data = sines, aes(color = "SINEs"), size = 2) +
    geom_point(data = lines, aes(color = "LINEs"), size = 2) +
    geom_point(data = ervs, aes(color = "ERVs"), size = 2) +
    geom_point(data = other_tes, aes(color = "Other"), size = 2) +
    scale_color_manual(values = c("SINEs" = "#55A35D",
                                  "LINEs" = "#E3C775",
                                  "ERVs" = "#39738D",
                                  "Other" = "#726FB9"),
                       breaks = c("LINEs", "SINEs", "ERVs", "Other"),
                       name = "TE Element") +
    labs(x = expression("log"[2]*"FC"),
         y = expression("-log"[10]*"adj p-value"),
         title = paste('DE TEs in', cell, master_DE_odc_df$broad_celltype, "\n")) +
    scale_x_continuous(breaks = seq(-10, 10, 2)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          legend.title = element_text(face = "bold", size = 14),
          legend.key = element_rect(fill = "transparent", colour = NA),
          legend.text = element_text(size = 12),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.major = element_line(color = "lightgray", linewidth = 0.3),
          panel.grid.minor = element_line(color = "lightgray", linewidth = 0.1, linetype = "dotted"))
  #storing plots into the plot list
  volcano_plots[[cell]] <- p
}

#calling plots to visualize each subpopulation
#volcano_plots[["Resting"]]
#volcano_plots[["Myelinating"]]
#volcano_plots[["Immature"]]
#volcano_plots[["Stress_associated"]]
#volcano_plots[["Ox_stress"]]

for(subpop in c("Resting","Myelinating","Immature","Stress_associated","Ox_stress")){
    outplot = paste0("SuppFigXE_",subpop,".png")
    ggsave(filename = file.path(plots_dir, outplot),
       plot = volcano_plots[[subpop]],
       width = 10,
       height = 6,
       dpi = 300)
}


### Supp Figure XF: ODCs Subpopulation DE PD v Ctrl TE Loci Dot Plot

#list of specific TE loci to plot
TE_loci <- c("AluY-dup35121:AluY:Alu:SINE", "AluYa5-dup1888:AluYa5:Alu:SINE", "AluYb8-dup2669:AluYb8:Alu:SINE", "AluYc-dup2455:AluYc:Alu:SINE", "AluYc3-dup184:AluYc3:Alu:SINE", "AluYe5-dup4334:AluYe5:Alu:SINE", "AluYg6-dup1733:AluYg6:Alu:SINE", "AluYk2-dup8577:AluYk2:Alu:SINE", "AluYk3-dup6184:AluYk3:Alu:SINE", "AluYm1-dup2128:AluYm1:Alu:SINE", "SVA-B-dup645:SVA-B:SVA:Retroposon", "SVA-B-dup681:SVA-B:SVA:Retroposon", "SVA-D-dup874:SVA-D:SVA:Retroposon", "SVA-E-dup1146:SVA-E:SVA:Retroposon", "SVA-E-dup253:SVA-E:SVA:Retroposon", "L1PA3-dup1904:L1PA3:L1:LINE", "L1PA3-dup8421:L1PA3:L1:LINE", "HERV9N-int-dup3:HERV9N-int:ERV1:LTR", "HERVK9-int-dup536:HERVK9-int:ERVK:LTR", "LTR13-dup271:LTR13:ERVK:LTR", "MER94-dup3820:MER94:hAT-Blackjack:DNA")
#creating a short-handed name for the TE loci
short_name <- sub(":.*", "", TE_loci)

#parameters used to make the TE loci dotplot
odc_dotplot <- DotPlot(odc_so, features = TE_loci, group.by = "cluster_subpopulation", dot.scale = 5, split.by = "diagnosis", cols = rev(diagnosis_cols)) +
  scale_x_discrete(labels = short_name) +
  scale_y_discrete(limits = c("Resting_Control", "Resting_Parkinsons", "Myelinating_Control", "Myelinating_Parkinsons", "Immature_Control", "Immature_Parkinsons", "Stress_associated_Control", "Stress_associated_Parkinsons", "Ox_stress_Control", "Ox_stress_Parkinsons")) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(size = 14, hjust = 1, angle = 90),
        axis.text.y = element_text(size = 12, angle = 45),
        panel.spacing = unit(0.2, "lines"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        aspect.ratio = 0.2,
        legend.text = element_text(size = 14, face = "plain"),
        legend.title = element_text(size = 14, face = "bold"))
#odc_dotplot

outplot <- "SuppFigXF.png"
ggsave(filename = file.path(plots_dir, outplot),
       plot = odc_dotplot,
       width = 10,
       height = 4,
       dpi = 300)

### Supp Figure XG: Stress Associated ODC Subpopulation Specific TE Loci FeaturePlot

#creating custom UMAP to display subpopulation enriched TE loci similar to Seurat's FeaturePlot
#obtaining the UMAP coordinates to be used for graphing
umap_coords <- Embeddings(odc_so, reduction = "umap.cca")

#obtaining TE loci expression data
#specifying which TE loci to ge expression data for
TE_loci <- "SVA-E-dup253:SVA-E:SVA:Retroposon"
#creating a short-hand name for the TE loci
short_name <- sub(":.*", "", TE_loci)
#getting the TE loci expression data
TE_expression <- GetAssayData(odc_so, assay = "RNA", slot = "data")[TE_loci,]

#creating a variable for cluster identities
#if needed, change idents to the subpopulation labels
#Idents(odc_so) <- "cluster_subpopulation"
clusters <- Idents(odc_so)

#creating the df needed for making the plot
umap_df <- data.frame(UMAP_1 = umap_coords[,1], #x-axis umap coordinates
                        UMAP_2 = umap_coords[,2], #y-axis umap coordinates
                        Expression = as.numeric(TE_expression),
                        Detected = TE_expression > 0,
                        Cluster = as.factor(clusters),
                        row.names = names(TE_expression))

#removing rows with NA values
umap_df <- umap_df[complete.cases(umap_df), ]

#checks to make sure that umap df information were obtained correctly
cat("Plot data dimensions:", dim(umap_df), "\n")
cat("UMAP_1 range:", range(umap_df$UMAP_1), "\n")
cat("UMAP_2 range:", range(umap_df$UMAP_2), "\n")

#creating separate layers to define order for plotting points
#plotting first layer of dots which are cells that do not express the TE loci
noexpression_cells <- umap_df[umap_df$Expression == 0, ]
#plotting second layer of dots which are cells that do express the TE loci
expressing_cells <- umap_df[umap_df$Expression > 0, ]
#this is ordered by expression level:
#lowest expression plotted first then highest expression plotted last
expressing_cells <- expressing_cells[order(expressing_cells$Expression), ]

#parameters used to create the UMAP
subpopulation_specific_TE_FP <- ggplot() +
  #first plotting all non-expressing cells with a smaller size
  geom_point(data = noexpression_cells,
             aes(x = UMAP_1, y = UMAP_2),
             color = "#E2E2E2", size = 0.25, alpha = 1) +
  #nex plotting expressing cells based on expression gradient
  geom_point(data = expressing_cells,
             aes(x = UMAP_1, y = UMAP_2, 
                 color = Expression),
             size = 0.75, alpha = 1) +
  #defining expression gradient cols
  scale_color_gradient(low = "snow2", high = "#005D8F", 
                      name = "Expression") +
  guides(size = "none") +
  theme_void() +
  theme(axis.line = element_line(color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text = element_text(color = "black", size = 12),
        axis.title.y = element_text(angle = 90),
        axis.title = element_text(color = "black", size = 12, margin = margin(10, 10, 10, 10)),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),  # Facet panel titles
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(title = paste(short_name, "\n"), x = "umapcca_1", y = "umapcca_2")
#subpopulation_specific_TE_FP

outplot <- "SuppFigXG.png"
ggsave(filename = file.path(plots_dir, outplot),
         plot = subpopulation_specific_TE_FP,
         width = 10,
         height = 8,
         dpi = 300)
