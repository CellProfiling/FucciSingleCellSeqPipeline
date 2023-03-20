# Snakemake analysis pipeline for FUCCI single-cell RNA-Seq data

## Single-cell proteogenomic analysis

This repository contains the _snakemake_ pipeline for analyzing the RNA sequencing data for ~1k single cells. The results of this single-cell RNA-Seq analysis provide a transcriptomic context to a proteomic analysis based on immunofluorescence staining of ~200k individual cells. For the code used to perform that single-cell proteogenomic analysis of the human cell cycle, please see the [CellProfiling/SingleCellProteogenomics](https://github.com/CellProfiling/SingleCellProteogenomics) repository.

## Single-cell sequencing files

The single-cell RNA-Seq data is available at GEO SRA under project number [GSE146773](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146773).

This data is downloaded automatically in this pipeline.

## Updating the Ensembl version

The genome and Ensembl versions are located at the top of the file `workflow/config/FucciSingleCell.yaml`.
These can be updated, and the references will be downloaded automatically.

## Usage

1) Clone repository and initialize submodules: `git clone --recurse-submodules https://github.com/CellProfiling/FucciSingleCellSeqPipeline.git && cd FucciSingleCellSeqPipeline/workflow`
1) Install conda: https://docs.conda.io/en/latest/miniconda.html
2) Create and activate setup environment: `conda env create -n fuccisetup -f envs/setup.yaml && conda activate fuccisetup`
4) Run the workflow: `snakemake --use-conda --conda-frontend mamba --cores 24 --resources mem_mb=100000`, where you can subsitute the max number of cores and max memory allocation. At least 54 GB of free memory should be available.

## Usage on cluster

In place of installing conda, you may need to activate it as a module, such as by `module load conda` and then follow the instructions to initialize it.

Adapt `config/cluster_config.yaml` for your needs.

In place of the last step above, you can use the scheduler like this:
`snakemake -j 500 --cores 16 --cluster-config config/cluster_config.yaml --latency-wait 60 --keep-going --use-conda --conda-frontend mamba --cluster "sbatch -A {cluster.account} -t {cluster.time} -N {cluster.nodes} --cpus-per-task {threads} -p {cluster.partition}"`

## Usage on protected access cluster

1) Clone repository and initialize submodules on your local machine: `git clone --recurse-submodules https://github.com/CellProfiling/FucciSingleCellSeqPipeline.git && cd FucciSingleCellSeqPipeline/workflow`
2) Install conda: https://docs.conda.io/en/latest/miniconda.html
3) Create and activate setup environment: `conda env create -n fuccisetup -f envs/setup.yaml && conda activate fuccisetup`
4) If running the pipeline on protected access computer, predownload files by running `snakemake -j 16 ../results/setup.txt` on a machine with internet access.
5) Make a tarball of the project with `cd ../.. && tar -cxvf FucciSingleCellSeqPipeline.zip FucciSingleCellSeqPipeline` and transfer it to the protected access cluster.
6) Load conda as a module on the protected access cluster, such as with `module load conda`, and follow the instructions to activate it.
7) Create and activate setup environment: `conda env create -n fuccisetup -f envs/setup.yaml && conda activate fuccisetup`
8) Adapt `config/cluster_config.yaml` for your needs.
9) Use the scheduler from snakemake like this:
`snakemake -j 500 --cores 16 --cluster-config config/cluster_config.yaml --latency-wait 60 --keep-going --use-conda --conda-frontend mamba --cluster "sbatch -A {cluster.account} -t {cluster.time} -N {cluster.nodes} --cpus-per-task {threads} -p {cluster.partition}"`

## Citation

Mahdessian, D.\*; Cesnik, A. J.\*; Gnann, C.; Danielsson, F.; Stenström, L.; Arif, M.; Zhang, C.; Le, T.; Johansson, F.; Shutten, R.; Bäckström, A.; Axelsson, U.; Thul, P.; Cho, N. H.; Carja, O.; Uhlén, M.; Mardinoglu, A.; Stadler, C.; Lindskog, C.; Ayoglu, B.; Leonetti, M. D.; Pontén, F.; Sullivan, D. P.; Lundberg, E. “Spatiotemporal dissection of the cell cycle with single cell proteogenomics.” Nature, 2021, 590, 649–654. \*Contributed equally. https://www.nature.com/articles/s41586-021-03232-9
