#!/bin/env Rscript

if(! require("BiocManager")){
    install.packages("BiocManager")
}

required_packages <- c("Seurat","presto","patchwork","glmGamPoi","ggrepel","ggplot2","dpylr,"BPCells")
for(package in required_packages){
    if(! require(package)){
        BiocManager::install(package)
    }
}
