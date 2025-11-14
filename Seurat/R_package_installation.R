#!/bin/env Rscript

if(! require("BiocManager")){
    install.packages("BiocManager")
}

required_packages <- c("Seurat","sctransform","ggrepel","ggvenn","ComplexUpset","BPCells","glmGamPoi","tidyr","Matrix","dplyr","forcats","ggplot2","patchwork","presto","scales","tibble")
for(package in required_packages){
    if(! require(package)){
        BiocManager::install(package)
    }
}
