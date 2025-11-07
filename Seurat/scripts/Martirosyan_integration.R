library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# 1Tb max vector size
options(future.globals.maxSize = 1000 * 1024 ^ 3)

## Creating Seurat object with unnormalized counts
data = Read10X(data.dir = "input", gene.column = 1)
scdata = CreateSeuratObject(counts = data, project = "Martirosyan_2024", min.cells = 3, min.features = 200)
rm(data)
 
## Get MT percentage - Not useful for T2T
#scdata[["percent.mt"]] = PercentageFeatureSet(scdata,pattern="^MT-")
 
## Add metadata
libdata = read.csv("annotations/Martirosyan_barcode_annotations.csv", header=F,row.names=1)
sampleID = libdata[match(rownames(scdata@meta.data), rownames(libdata)),1]
patientID = libdata[match(rownames(scdata@meta.data), rownames(libdata)),3]
diagnosis = libdata[match(rownames(scdata@meta.data), rownames(libdata)),4]
age = libdata[match(rownames(scdata@meta.data), rownames(libdata)),5]
gender = libdata[match(rownames(scdata@meta.data), rownames(libdata)),6]
pmi = libdata[match(rownames(scdata@meta.data), rownames(libdata)),7]
rin = libdata[match(rownames(scdata@meta.data), rownames(libdata)),8]
midLewy = libdata[match(rownames(scdata@meta.data), rownames(libdata)),9]
amygLewy = libdata[match(rownames(scdata@meta.data), rownames(libdata)),10]
fcxLewy = libdata[match(rownames(scdata@meta.data), rownames(libdata)),11]
cerad = libdata[match(rownames(scdata@meta.data), rownames(libdata)),12]
braak = libdata[match(rownames(scdata@meta.data), rownames(libdata)),13]
scdata[["sampleID"]] = sampleID
scdata[["patientID"]] = patientID
scdata[["diagnosis"]] = diagnosis
scdata[["age"]] = age
scdata[["gender"]] = gender
scdata[["pmi"]] = pmi
scdata[["rin"]] = rin
scdata[["midLewy"]] = midLewy
scdata[["amygLewy"]] = amygLewy
scdata[["fcxLewy"]] = fcxLewy
scdata[["cerad"]] = cerad
scdata[["braak"]] = braak

## Add EmptyDrops and DoubletFinder results
preprocess = read.csv("annotations/Martirosyan_ED_DF_output.csv", header=T,row.names = 1)
ED = preprocess[match(rownames(scdata@meta.data), rownames(preprocess)),4]
DF = preprocess[match(rownames(scdata@meta.data), rownames(preprocess)),6]
scdata[["EmptyDropsFDR"]] = ED
scdata[["DoubletFinder"]] = DF

write.table(scdata@meta.data, "annotations/Martirosyan_full_metadata.csv",sep=",",row.names = T, col.names = T, quote = F)

## Filter out empty drops and doublets
emptyDropsThreshold = 0.001
print(paste("Original matrix size:",ncol(scdata)))
scdata = subset(x = scdata, subset = EmptyDropsFDR <= emptyDropsThreshold & DoubletFinder == "Singlet")
print(paste("Filtered matrix size:",ncol(scdata)))

## Split and preprocess data
scdata[["RNA"]] = split(scdata[["RNA"]], f = scdata$diagnosis)
scdata = NormalizeData(scdata)
scdata = FindVariableFeatures(scdata)
all.genes = rownames(scdata)[grep(":",rownames(scdata),invert=T)]
scdata = ScaleData(scdata, features = all.genes)
scdata = RunPCA(scdata, npcs = 30)

## Integrate data
scdata = IntegrateLayers(object = scdata, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", features = all.genes)

scdata[["RNA"]] = JoinLayers(scdata[["RNA"]])
saveRDS(scdata,"Martirosyan_2024_integrated.rds")

scdata = FindNeighbors(scdata, reduction = "integrated.cca", dims = 1:30)
scdata = FindClusters(scdata, resolution = 0.5, algorithm = 4, method = "igraph", cluster.name = "Leiden")
scdata = RunUMAP(scdata, reduction = "integrated.cca", reduction.name = "umap.cca", dims = 1:30)

DefaultAssay(scdata) = "RNA"
saveRDS(scdata,"Martirosyan_2024_integrated.rds")

markers = read.table("celltype_markers.txt", sep="\t", row.names = 1, header = F)

pdf("Martirosyan_Leiden_clusters.pdf")
Idents(scdata) = "Leiden"
DimPlot(scdata, label = TRUE, repel = TRUE)
DoHeatmap(subset(scdata, downsample = 100), features = rownames(markers)) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
Idents(scdata) = "diagnosis"
DimPlot(scdata, label = TRUE, repel = TRUE)
Idents(scdata) = "gender"
DimPlot(scdata, label = TRUE, repel = TRUE)
dev.off()
