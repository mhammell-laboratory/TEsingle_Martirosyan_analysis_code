#!/bin/env Rscript

if(! require("BiocManager")){
    install.packages("BiocManager")
}

required_packages <- c("Seurat","DropletUtils","DoubletFinder")
for(package in required_packages){
    if(! require(package)){
        BiocManager::install(package)
    }
}
