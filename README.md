# TEsingle analysis of Martirosyan dataset
Code used to analyze Martirosyan et al. dataset with TEsingle

## Overview
This repository contains code used to analyze the [Martirosyan et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38245794/) dataset with TEsingle.

The pipeline is divided into four portions:
1. Setting up STAR index and running TEsingle.
2. Running emptyDrops and DoubletFinder on TEsingle output.
3. Aggregating TEsingle runs into a single matrix.
4. Running downstream analyses with Seurat.

The raw data for this single nuclei dataset is available on GEO ([GSE243639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243639)).
