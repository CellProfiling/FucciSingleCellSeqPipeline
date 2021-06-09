# Snakemake analysis pipeline for FUCCI single-cell RNA-Seq data

## Single-cell proteogenomic analysis

This repository contains the _snakemake_ pipeline for analyzing the RNA sequencing data for ~1k single cells. The results of this single-cell RNA-Seq analysis provide a transcriptomic context to a proteomic analysis based on immunofluorescence staining of ~200k individual cells. For the code used to perform that single-cell proteogenomic analysis of the human cell cycle, please see the [CellProfiling/SingleCellProteogenomics](https://github.com/CellProfiling/SingleCellProteogenomics) repository.

## Single-cell sequencing files

The single-cell RNA-Seq data is available at GEO SRA under project number [GSE146773](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146773).

This data is downloaded automatically in this pipeline.

## Updating the Ensembl version

The genome and Ensembl versions are located at the top of the file `Snakefile`.
These can be updated, and the references will be downloaded automatically.

## Usage

1) Clone repository and initialize submodules: `git clone --recurse-submodules https://github.com/CellProfiling/FucciSingleCellSeqPipeline.git && cd FucciSingleCellSeqPipeline/workflow`
1) Install conda: https://docs.conda.io/en/latest/miniconda.html
2) Install snakemake using conda: `conda install -c conda-forge snakemake-minimal`
4) Run the workflow: `snakemake --use-conda --cores 24 --resources mem_mb=100000`, where you can subsitute the max number of cores and max memory allocation. At least 54 GB of free memory should be available.

## Usage on cluster

1) Clone repository and initialize submodules: `git clone --recurse-submodules https://github.com/CellProfiling/FucciSingleCellSeqPipeline.git && cd FucciSingleCellSeqPipeline/workflow`
2) If 

## Citation

Mahdessian, D.\*; Cesnik, A. J.\*; Gnann, C.; Danielsson, F.; Stenström, L.; Arif, M.; Zhang, C.; Le, T.; Johansson, F.; Shutten, R.; Bäckström, A.; Axelsson, U.; Thul, P.; Cho, N. H.; Carja, O.; Uhlén, M.; Mardinoglu, A.; Stadler, C.; Lindskog, C.; Ayoglu, B.; Leonetti, M. D.; Pontén, F.; Sullivan, D. P.; Lundberg, E. “Spatiotemporal dissection of the cell cycle with single cell proteogenomics.” Nature, 2021, 590, 649–654. \*Contributed equally. https://www.nature.com/articles/s41586-021-03232-9
