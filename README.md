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
- presto (tested on v1.0.0) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- patchwork (tested on v1.3.0) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- glmGamPoi (tested on v1.12.1) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- ggrepel (tested on v0.9.6) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- dpylr (tested on v1.1.4) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)
- BPCells (tested on v0.3.0) : [installation code](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/tree/main/Seurat/R_package_installation.R)

### Obtaining analysis pipeline code
```
$ git clone https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code.git
```

## Folder structure
- `LICENSE`: BSD 3-clause license
- `README.md`: README file
- `QC`: Contains scripts to perform emptyDrops and DoubletFinder on TEsingle output, as well as a package installation script and `sessionInfo()` output.
- `Seurat`: Organized to enable data analysis with Seurat, with input and output folders pre-made.
    - `annotations`: Annotation files and QC output for the dataset should be placed here.
    - `cellclusters`: R objects and associated metadata for each subsetted broad cell types (e.g. Neurons) are output here.
    - `differential_analysis`: Output from Parkinson's disease vs control differential analysis for broad cell types are saved here.
    - `input`: The `barcodes.tsv.gz`, `features.tsv.gz` and `matrix.mtx.gz` from the aggregated TEsingle runs should be stored here.
    - `scripts`: Contains all the R scripts used in the analysis.
    - `subpopulation`: Contains the subpopulation analysis output.
        - `diffExp`: Output from Parkinson's disease vs control differential analysis for cell subpopulations are saved here.
        - `plots`: Plots from subpopulation analysis are saved here.
        - `rds`: R objects from subpopulation analysis are saved here.
        - `tables`: Leiden markers from subpopulation analysis are saved here.
- `TEsingle`: Contains scripts to set up the STAR index and to run `STAR` and `TEsingle`.
- `aggregate_output`: Contains scripts to aggregate TEsingle output into a single matrix.
    - `src`: additional scripts/files for 

## How to use the pipeline

## Limitations
This pipeline has been designed for testing a specific version of STAR, TEsingle and Seurat (see versions used in the dependency section). Newer versions of the software may have changed parameters and output, and could lead to different results.

### Using aggregate output scripts on SLURM
The aggreate output scripts depend on several supporting code/files in the src subfolder. When submitting to SLURM, the src subfolder might be unlinked from the folder containing the accuracy scripts, leading to the following errors:
```
Can't open perl script ".../TEsingle_aggregate.pl": No such file or directory
Can't open perl script ".../annotate_barcode.pl": No such file or directory
```
In order to fix this, you may need to change the following lines:
- [`aggregate_TEsingle_runs.sh`](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/blob/main/aggregate_output/aggregate_TEsingle_runs.sh#L18)
- [`make_aggregated_annotations.sh`](https://github.com/mhammell-laboratory/TEsingle_Martirosyan_analysis_code/blob/main/aggregate_output/make_aggregated_annotations.sh#L19)

from
```
SCRIPTDIR=$(dirname $0)
```
to
```
SCRIPTDIR=/path/to/aggregate_output
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
