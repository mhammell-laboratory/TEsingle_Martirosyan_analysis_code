library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

scdata = readRDS("Martirosyan_2024_integrated.rds")

markers = read.table("official_markers.txt", sep="\t", row.names = 1, header = F)

leiden = as.vector(scdata$Leiden)
celltype = as.vector(scdata$Leiden)

celltype[leiden == 1 | leiden == 9 | leiden == 15] = "DopaN"
celltype[leiden == 10] = "Pericyte_Endothelial_VLMC"
celltype[leiden == 11] = "Exc"
celltype[leiden == 12] = "OtherNeuron"
celltype[leiden == 20] = "Inh"
celltype[leiden == 6 | leiden == 7 | leiden == 16 | leiden == 19] = "Astro"
celltype[leiden == 2 | leiden == 4 | leiden == 5 | leiden == 13] = "Oligo"
celltype[leiden == 14 | leiden == 23] = "Unknown"
celltype[leiden == 17] = "Tcell"
celltype[leiden == 3 | leiden == 18 | leiden == 21] = "Micro"
celltype[leiden == 8 | leiden == 22] = "OPC"

scdata[["celltype"]] = celltype

pdf("Martirosyan_labelled_clusters.pdf")
Idents(scdata) = "Leiden"
DimPlot(scdata, label = TRUE, repel = TRUE)
DoHeatmap(subset(scdata, downsample = 100), features = rownames(markers)) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
Idents(scdata) = "celltype"
DimPlot(scdata, label = TRUE, repel = TRUE)
DoHeatmap(subset(scdata, downsample = 100), features = rownames(markers)) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
Idents(scdata) = "diagnosis"
DimPlot(scdata, label = TRUE, repel = TRUE)
Idents(scdata) = "gender"
DimPlot(scdata, label = TRUE, repel = TRUE)
dev.off()

saveRDS(scdata,"Martirosyan_2024_integrated.rds")

for(cell in c("Neuron","Astro","Micro","OPC","Oligo")){
    if(cell == "Neuron"){
        cell = c("DopaN","Exc","Inh","OtherNeuron")
    }
    subdata = subset(x = scdata, idents = cell)
    output = paste0("cellclusters/Martirosyan_",cell,".rds")
    saveRDS(subdata,output)
    metadata = subdata@meta.data
    output = paste0("cellclusters/Martirosyan_",cell,"_metadata.txt")
    write.table(metadata,output,row.names=T,col.names=T,sep="\t",quote=F)
}

category = paste(as.vector(scdata$diagnosis), as.vector(scdata$celltype), sep="_")
scdata[["category"]] = category
Idents(scdata) = "category"

for(cell in c("DopaN","Astro","Micro", "Exc","Inh", "OtherNeuron", "Oligo","OPC")){
    disease = paste0("Parkinsons_",cell)
    control = paste0("Control_",cell)
    diff = FindMarkers(object = scdata, ident.1 = disease, ident.2 = control, test.use = "LR")
    outfile = paste0("differential_analysis/",cell, "_PDvCtrl_LR_results.txt")
    write.table(diff, outfile, sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
    sig = diff[!is.na(diff$p_val_adj) & diff$p_val_adj < 0.05,]
    outfile = paste0("differential_analysis/",cell, "_PDvCtrl_LR_sighits.txt")
    write.table(sig, outfile, sep="\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
}
