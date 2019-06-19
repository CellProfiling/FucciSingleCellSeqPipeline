# Snakemake analysis pipeline for scRNA-Seq data

## Usage

1) install conda: https://docs.conda.io/en/latest/miniconda.html
2) create the conda environment: `conda env create --file environment.yaml --name cellquant`
3) activate the conda environment: `conda activate cellquant`
4) run the workflow: recommended command is `snakemake --cores 24 --resources mem_mb=100000`, where you can subsitute the max number of cores and max memory allocation. The memory allocation should be at least 50000 MB if possible. It might work with 32000 MB, but no guarantees.