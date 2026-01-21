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

## Installation
### Dependencies
- [STAR](https://github.com/alexdobin/STAR) v2.7.11b : [installation instructions](https://github.com/alexdobin/STAR/blob/master/README.md)
- [TEsingle](https://github.com/mhammell-laboratory/TEsingle) v1.0 : [installation_instructions](https://github.com/mhammell-laboratory/TEsingle/blob/main/README.rst)
- Perl >= v5 (tested on v5.28.0) Unix/Linux and MacOSX have perl installed. Please see [here](https://www.perl.org/get.html) for Windows.
- R >= v4 (tested on v4.3.3) : [installation instructions](https://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-R-be-installed_003f)
- Seurat >= v5 (tested on v5.3.0) : [installation instructions](https://satijalab.org/seurat/articles/install.html)
- DropletUtils (tested on v1.20.0) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/QC/R_package_installation.R)
- DoubletFinder (tested on v2.0.6) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/QC/R_package_installation.R)
- sctranform (tested on v0.4.2) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- BPCells (tested on v0.3.0) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- presto (tested on v1.0.0) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- patchwork (tested on v1.3.0) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- glmGamPoi (tested on v1.12.1) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- ggrepel (tested on v0.9.6) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- ggvenn (tested on v0.1.19) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- ComplexUpset (tested on v1.3.3) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- dpylr (tested on v1.1.4) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- tidyr (tested on v1.3.1) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)

### Obtaining analysis pipeline code
```
$ git clone https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code.git
```
## Folder structure
- `LICENSE`: BSD 3-clause license
- `README.md`: README file
- `QC`: Contains scripts to perform emptyDrops and DoubletFinder on TEsingle output, as well as a package installation script and `sessionInfo()` output.
    - `src`: additional script for QC.
- `Seurat`: Organized to enable data analysis with Seurat, with input and output folders pre-made. Also contains a package installation script and `sessionInfo()` output.
    - `annotations`: Annotation files and QC output for the dataset should be placed here.
    - `cellclusters`: R objects and associated metadata for each subsetted broad cell types (e.g. Neurons) are output here.
    - `differential_analysis`: Output from Parkinson's disease vs control differential analysis for broad cell types are saved here.
    - `figures`: Contains the figures generated from the Martirosyan analysis.
    - `input`: The `barcodes.tsv.gz`, `features.tsv.gz` and `matrix.mtx.gz` from the aggregated TEsingle runs should be stored here.
    - `scripts`: Contains all the R scripts used in the analysis.
    - `subpopulation`: Contains the subpopulation analysis output.
        - `diffExp`: Output from Parkinson's disease vs control differential analysis for cell subpopulations are saved here.
        - `markers`: Subpopulation marker lists are located here.
        - `plots`: Plots from subpopulation analysis are saved here.
        - `rds`: R objects from subpopulation analysis are saved here.
        - `tables`: Leiden markers from subpopulation analysis are saved here.
- `TEsingle`: Contains scripts to set up the STAR index and to run `STAR` and `TEsingle`.
- `aggregate_output`: Contains scripts to aggregate TEsingle output into a single matrix.
    - `src`: additional scripts/files for aggregating TEsingle outputs.

## How to use the pipeline
### Quantification with TEsingle
#### Initial setup
##### System requirments
- CPU: 10
- Memory: 3G per CPU (30G total)
- Allowed time: up to 12 hours

You will need to download the FASTQ files from GEO ([GSE243639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243639)) using software such as the [NCBI SRA Toolkit](https://hpc.nih.gov/apps/sratoolkit.html). 

You will then need to generate the STAR index for alignment using the provided `setup_STAR_index.sh` script:
```
# Run locally
$ sh /path/to/setup_STAR_index.sh T2T_STAR_index
# Submit to SLURM
$ sbatch /path/to/setup_STAR_index/sh T2T_STAR_index
```
This will create a folder called `T2T_STAR_index` that should contain the T2T genome FASTA, both gene and TE GTF required for TEsingle, a barcode whitelist (from 10x Genomics v3) and a STAR index for alignment.

#### Running TEsingle
##### System requirments
- CPU: 10
- Memory: 30G per CPU (300G total)
- Allowed time: up to 5 days

You will then need to run TEsingle on each sample:
```
# Run locally
$ sh /path/to/run_STAR_TEsingle.sh /path/to/T2T_STAR_index [Read 1 FASTQ] [Read 2 FASTQ]
# Submit to SLURM
$ sbatch /path/to/run_STAR_TEsingle.sh /path/to/T2T_STAR_index [Read 1 FASTQ] [Read 2 FASTQ]
```
Each run will generate three files: 
- `XXXX_TEsingle.annots`: Contains the features used in the quantification
- `XXXX_TEsingle.cbcs`: Contains the barcodes detected in the library
- `XXXX_TEsingle.mtx`: Matrix containing counts

### Library QC post quantification
##### System requirments
- CPU: 1
- Memory: 100G
- Allowed time: up to 5 days

After running TEsingle, you will need to run emptyDrops and DoubletFinder to identify empty droplets and doublets from the datasets. You will need to provide the path to the `XXXX_TEsingle.mtx` file, with the assumption that the other two files (`.annots` and `.cbcs`) are in the same folder.
```
# Run locally
$ sh /path/to/emptyDrops_DoubletFinder_QC.sh /path/to/XXXX_TEsingle.mtx
# Submit to SLURM
$ sbatch /path/to/emptyDrops_DoubletFinder_QC.sh /path/to/XXXX_TEsingle.mtx
```
The code will create a folder (`XXXX_TEsingle`, depends on the file name) containing the expected output (`barcodes.tsv`, `features.tsv` and `matrix.mtx`) commonly found from Cell Ranger output. This will allow import into Seurat/emptyDrops/DoubletFinder for the QC step, resulting in a file (`XXXX_TEsingle_QC_output.txt`) in the `preprocess` folder.

### Aggregating output
##### System requirments
- CPU: 10
- Memory: 1G per CPU (10G total)
- Allowed time: up to 12 hours

After the QC step, the TEsingle runs and the QC output need to be aggregated into a single dataset for use downstream by Seurat.

To aggregate the TEsingle runs, run the following script:
```
# Run locally
$ sh /path/to/aggregate_TEsingle_runs.sh /path/to/Seurat/input /path/to/*_TEsingle.mtx
# Submit to SLURM
$ sbatch /path/to/aggregate_TEsingle_runs.sh /path/to/Seurat/input /path/to/*_TEsingle.mtx
```
This code will combine the barcodes and counts for each TEsingle run (i.e. each sample) and put them into an aggregated matrix. It is highly recommended that you set the output folder into the `input` subfolder in the `Seurat` folder, as the downstream Seurat analyses work within that folder structure.

To generate annotations for the downstream analysis, run the following script:
```
# Run locally
$ sh /path/to/make_aggregated_annotations.sh /path/to/Seurat/input /path/to/preprocess/*_TEsingle_QC_output.txt
# Submit to SLURM
$ sbatch /path/to/make_aggregated_annotations.sh /path/to/Seurat/input /path/to/preprocess/*_TEsingle_QC_output.txt
```
This code will take the Martirosyan sample annotation and assign it to each detected barcode (`Martirosyan_barcode_annotations.csv`). It will also combine the QC results into a single file for subsequent loading into Seurat (`Martirosyan_ED_DF_output.csv`). It is recommended that these two files are copied/moved into the `annotations` subfolder in the `Seurat` folder.

### Seurat analysis
The code for the analysis using Seurat is designed to work in the folder structure created in the [`Seurat`](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat) subfolder of the GitHub repository. The R scripts should be run in the base folder, which will take inputs from relevant folders (`input` and `annotations`) and generate output either in the base folder, or various subfolders (e.g. `differential_analysis`)

```
$ echo $PWD
/path/to/Seurat
$ Rscript scripts/Martirosyan_integration.R
$ Rscript scripts/Martirosyan_celltype_DE.R
$ Rscript scripts/Martirosyan_subpopulation_integration.R cellclusters/Martirosyan_Neuron.rds
$ Rscript scripts/Martirosyan_neuronal_subpopulation.R
$ Rscript scripts/Martirosyan_subpopulation_integration.R cellclusters/Martirosyan_Astro.rds
$ Rscript scripts/Martirosyan_astrocytic_subpopulation.R
$ Rscript scripts/Martirosyan_subpopulation_integration.R cellclusters/Martirosyan_Micro.rds
$ Rscript scripts/Martirosyan_microglial_subpopulation.R
$ Rscript scripts/Martirosyan_subpopulation_integration.R cellclusters/Martirosyan_Oligo.rds
$ Rscript scripts/Martirosyan_oligodendrocytic_subpopulation.R
```
#### Data integration
##### System requirments
- CPU: 1
- Memory: 1Tb
- Allowed time: up to 5 days

This step takes the combined matrix from the TEsingle runs, and adds the annotation and QC results for the identified barcodes to the metadata. It will then filter out barcodes that were determined as empty droplets (FDR >0.001) or suspected doublets, perform integration based on diagnosis (PD or control), perform dimenional reduction (PCA and UMAP) and find clusters of cell with similar profile (using Leiden algorithm). 
```
$ Rscript scripts/Martirosyan_integration.R
```
The code outputs a R object of the integrated dataset and a UMAP plot PDF showing the Leiden clusters, diagnosis and gender (for QC), and a heatmap of broad cell type markers for the Leiden clusters (to enable cluster annotation).

#### Cell type annotation, subsetting and differential analysis
##### System requirments
- CPU: 1
- Memory: 600Gb
- Allowed time: up to 5 days

This step takes the integrated Seurat R object, annotates the Leiden clusters to broad cell types, subset the broad cell types into separate R objects (for subpopulation analysis) and runs differential analysis between Parkinson's Disease and controls for many of the broad cell types.
```
$ Rscript scripts/Martirosyan_celltype_DE.R
```
The code outputs a PDF showing the annotated Leiden clusters and their marker expression, R objects containing only barcodes corresponding to various broad cell types in the `cellclusters` subfolder, and differential analysis output (logistic regression) in the `differential_analysis` subfolder (all features tested, `[celltype]_PDvCtrl_LR_results.txt`, and significant hits (FDR < 0.05), `[celltype]_PDvCtrl_LR_sighits.txt`).

#### Find celltype markers based on called population
##### System requirements
- CPU: 1
- Memory: 600Gb
- Allowed time: up to 5 days

This step takes the integrated Seurat R object with the broad cell type annotations, and compare each cell type to the others to obtain a list of markers specific to each annotation.
```
$ Rscript scripts/Martirosyan_find_celltype_markers.R
```
The code outputs the cell markers (determined by differential expression using logistric regression in the `differential_analysis` subfolder (all features tested, `[celltype]_cellvcell_LR_results.txt`, and significant hits (FDR < 0.05), `[celltype]_cellvcell_LR_sighits.txt`). Combined files for all features tested, `Martirosyan_cellvcell_LR_results.txt`, and all significant hits (FDR < 0.05), `Martirosyan_cellvcell_LR_sighits.txt` are also generated in the `differential_analysis` subfolder.

#### Subpopulation integration
##### System requirments
- CPU: 1
- Memory: 600Gb
- Allowed time: up to 5 days

This step takes the subsetted R objects, re-normalize the counts (using `sctransform`), re-integrate based on diagnosis, perform dimensional reduction and find clusters of similar expression profile (using Leiden algorithm).
```
$ Rscript scripts/Martirosyan_subpopulation_integration.R cellclusters/Martirosyan_Neuron.rds
$ Rscript scripts/Martirosyan_subpopulation_integration.R cellclusters/Martirosyan_Astro.rds
$ Rscript scripts/Martirosyan_subpopulation_integration.R cellclusters/Martirosyan_Micro.rds
$ Rscript scripts/Martirosyan_subpopulation_integration.R cellclusters/Martirosyan_Oligo.rds
```
The code outputs RDS objects at various stages of the pipeline into the `subpopulation/rds` subfolder, UMAP plots showing Leiden clusters and various metadata categories (for QC) into the `subpopulation/plots` subfolder, and a CSV containing differential "markers" for each of the identified Leiden cluster in the `subpopuation/tables` subfolder.

#### Marker heatmaps of cell subpopulations
##### System requirments
- CPU: 1
- Memory: 600Gb
- Allowed time: up to 5 days

This step takes the integrated subpopulation RDS object in the `subpopulation/rds` subfolder and generate expression heatmaps of cell type/subpopulation markers.
```
$ Rscript scripts/Martirosyan_subpopulation_heatmap.R Neuron
$ Rscript scripts/Martirosyan_subpopulation_heatmap.R Astro
$ Rscript scripts/Martirosyan_subpopulation_heatmap.R Micro
$ Rscript scripts/Martirosyan_subpopulation_heatmap.R Oligo
```
The code outputs heatmap plots into the `subpopulation/plots` subfolder. 

#### Subpopulation differential analysis
##### System requirments
- CPU: 1
- Memory: 600Gb
- Allowed time: up to 5 days

This step takes the integrated cell type subset R objects, labels the Leiden clusters as a particular cell subpopulation, and perform differential analyses (logistic regression) in Parkinson's Disease vs control for each subpopulation.
```
$ Rscript scripts/Martirosyan_neuronal_subpopulation.R
$ Rscript scripts/Martirosyan_astrocytic_subpopulation.R
$ Rscript scripts/Martirosyan_microglial_subpopulation.R
$ Rscript scripts/Martirosyan_oligodendrocytic_subpopulation.R
```
The code outputs labelled RDS objects at various stages of the pipeline into the `subpopulation/rds` subfolder, and differential analysis results into the `subpopuation/diffExp` subfolder.

### Figure generation
##### System requirments
- CPU: 1
- Memory: 600Gb
- Allowed time: up to 12 hours

The following scripts are used to generate figures for the paper, specifically figures 4B to 4F, 5A to 5G, 6A to 6G, 7A to 7G and supplementary figures XA to XG. The scripts should be run in the base folder used for Seurat analysis, and will output PNG of the figures into the `figures` subfolder.
```
$ Rscript scripts/Martirosyan_Fig4_plots.R
$ Rscript scripts/Martirosyan_Fig5_plots.R
$ Rscript scripts/Martirosyan_Fig6_plots.R
$ Rscript scripts/Martirosyan_Fig7_plots.R
$ Rscript scripts/Martirosyan_SuppFig_plots.R
```

## Limitations
This pipeline has been designed for testing a specific version of STAR, TEsingle and Seurat (see versions used in the dependency section). Newer versions of the software may have changed parameters and output, and could lead to different results. 

Cell type annotation was done manually by comparing Seurat Leiden clusters wth common neuronal & glial markers, as well as some previous published subpopulation markers (see our [paper]() for details). An automated approach with more comprehensive marker lists may give slightly different annotations, especially for subpopulations not previously well characterized.

### Using QC and aggregate output scripts on SLURM
The QC and aggregate output scripts depend on several supporting code/files in the `src` subfolder. When submitting to SLURM, the `src` subfolder might be unlinked from the folder containing the accuracy scripts, leading to the following errors:
```
Fatal error: cannot open file '.../emptyDrops_DoubletFinder_QC.R': No such file or directory
Can't open perl script ".../TEsingle_aggregate.pl": No such file or directory
Can't open perl script ".../annotate_barcode.pl": No such file or directory
```
In order to fix this, you may need to change the following lines:
- [`emptyDrops_DoubletFinder_QC.sh`](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/blob/main/QC/emptyDrops_DoubletFinder_QC.sh#L15)
- [`aggregate_TEsingle_runs.sh`](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/blob/main/aggregate_output/aggregate_TEsingle_runs.sh#L18)
- [`make_aggregated_annotations.sh`](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/blob/main/aggregate_output/make_aggregated_annotations.sh#L19)

from
```
SCRIPTDIR=$(dirname $0)
```
to
```
SCRIPTDIR=/path/to/QC                 # For QC script
SCRIPTDIR=/path/to/aggregate_output   # For aggregate output scripts
```

## Citation
To be provided

## License
The code in this repository is distributed under the BSD 3-clause license per ASAP Open Access (OA) policy, which facilitates the rapid and free exchange of scientific ideas and ensures that ASAP-funded research fund can be leveraged for future discoveries.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

A copy of BSD 3-clause licence is included along with the software, and can be accessed [here](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/blob/main/LICENSE).

## Acknowledgments

- Contributors: Esther Cheng, Talitha Forcier, Oliver Tam & Molly Gale Hammell

This research was funded in whole or in part by Aligning Science Across Parkinson’s (ASAP-000520) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF).
