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
#### Data integration
##### System requirments
- CPU: 1
- Memory: 1Tb
- Allowed time: up to 5 days

## Limitations
This pipeline has been designed for testing a specific version of STAR, TEsingle and Seurat (see versions used in the dependency section). Newer versions of the software may have changed parameters and output, and could lead to different results. 

Cell type annotation was done manually by comparing Seurat Leiden clusters wth common neuronal & glial markers, as well as some previous published subpopulation markers (see our [paper]() for details). An automated approach with more comprehensive marker lists may give slightly different annotations, especially for subpopulations not previously well characterized.

### Using QC and aggregate output scripts on SLURM
The QC and aggregate output scripts depend on several supporting code/files in the src subfolder. When submitting to SLURM, the src subfolder might be unlinked from the folder containing the accuracy scripts, leading to the following errors:
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

This research was funded in whole by Aligning Science Across Parkinson’s (ASAP-000520) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF).
