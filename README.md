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
2) Install snakemake using conda: `conda install -c conda-forge -c bioconda snakemake-minimal`
4) Run the workflow: `snakemake --use-conda --cores 24 --resources mem_mb=100000`, where you can subsitute the max number of cores and max memory allocation. At least 54 GB of free memory should be available.

## Usage on cluster

In place of the last step above, you can use the scheduler like this:
`snakemake -j 500 --cores 16 --cluster-config cluster_config.yaml --latency-wait 60 --keep-going --use-conda --cluster "sbatch -A {cluster.account} -t {cluster.time} -N {cluster.nodes} --cpus-per-task {threads} -p {cluster.partition}"`

Replace 99 with the number of cores specified above in workflow/rules/align.smk and workflow/rules/quant.smk.

Where `cluster_config.yaml` may look like this:
```
__default__:
    account: sens2020535
    partition: core
    time: 2-0 # time limit for each job
    nodes: 1
    ntasks-per-node: 16 #Request n cores be allocated per node.
    output: ../results/slurmout/spritz-%j.out
    error: ../results/slurmerr/spritz-%j.err
```

## Usage on protected access cluster

1) Clone repository and initialize submodules: `git clone --recurse-submodules https://github.com/CellProfiling/FucciSingleCellSeqPipeline.git && cd FucciSingleCellSeqPipeline/workflow`
1) Install conda: https://docs.conda.io/en/latest/miniconda.html
2) Install snakemake using conda: `conda install -c conda-forge -c bioconda snakemake-minimal`
2) If running the pipeline on protected access computer, predownload files by running `snakemake -j 16 ../results/setup.txt` on a machine with internet access.
4) Make a tarball of the project with `cd ../.. && tar -cxvf FucciSingleCellSeqPipeline.zip FucciSingleCellSeqPipeline`

## Citation

Mahdessian, D.\*; Cesnik, A. J.\*; Gnann, C.; Danielsson, F.; Stenström, L.; Arif, M.; Zhang, C.; Le, T.; Johansson, F.; Shutten, R.; Bäckström, A.; Axelsson, U.; Thul, P.; Cho, N. H.; Carja, O.; Uhlén, M.; Mardinoglu, A.; Stadler, C.; Lindskog, C.; Ayoglu, B.; Leonetti, M. D.; Pontén, F.; Sullivan, D. P.; Lundberg, E. “Spatiotemporal dissection of the cell cycle with single cell proteogenomics.” Nature, 2021, 590, 649–654. \*Contributed equally. https://www.nature.com/articles/s41586-021-03232-9
