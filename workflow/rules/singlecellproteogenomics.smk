SCP_VERSION="1.2"

rule download_scp:
    input:
        quant=[
            "../results/quant/IsoformToGene.csv.gz",
            "../results/quant/Counts.csv",
            "../results/quant/IdsToNames.csv",
            "../results/quant/Tpms.csv",
            "../results/quant/Counts_Isoforms.csv",
            "../results/quant/IdsToNames_Isoforms.csv",
            "../results/quant/Tpms_Isoforms.csv"
        ],
        velocity=[
            "../results/velocity/a.loom",
            "../results/velocity/a.obs_names.csv",
        ],
    output: directory("../results/SingleCellProteogenomics/newinputs/RNAData/")
    conda: "../envs/downloads.yaml"
    log: "../results/downoad_scp.log"
    benchmark: "../results/downoad_scp.benchmark"
    params:
        folder=lambda w, output: os.path.dirname(output[0]),
        url=f"https://github.com/CellProfiling/SingleCellProteogenomics/archive/refs/tags/{SCP_VERSION}.tar.gz"
    shell:
        #"(wget -c {params.url} -O - | tar --strip-components=1 -xzC {output} && "
        "(mkdir -p {output} && cp {input.quant} {input.velocity} {output}) &> {log}"

rule SingleCellProteogenomics_run:
    input: "../results/SingleCellProteogenomics/newinputs/RNAData/"
    output:
        protein="../results/final/ProteinPseudotimePlotting.csv.gz",
        rna="../results/final/RNAPseudotimePlotting.csv.gz",
    conda: "../../results/SingleCellProteogenomics/workflow/envs/enviro.yaml"
    log: "../results/SingleCellProteogenomics.log"
    benchmark: "../results/SingleCellProteogenomics.log"
    params: scp_results=lambda w, input: os.path.dirname(input[0]),
    threads: 1
    shell:
        "(cd ../results/SingleCellProteogenomics/workflow && snakemake -j {threads} && "
        "cp {params.scp_results}/ProteinPseudotimePlotting.csv.gz {output.protein} && "
        "cp {params.scp_results}/RNAPseudotimePlotting.csv.gz {output.rna}) &> {log}"
