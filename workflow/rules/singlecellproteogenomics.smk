import os

INPUTID="149ICTtieYjuKWZoLwRLzimwff0n6eWqw"

rule SingleCellProteogenomics_copyResults:
    input:
        quant=[
            "../results/quant/Counts.csv",
            "../results/quant/Counts.csv.ercc.csv",
            "../results/quant/Tpms.csv",
            "../results/quant/Tpms.csv.ercc.csv",
            "../results/quant/Counts_Isoforms.csv",
            "../results/quant/Counts_Isoforms.csv.ercc.csv",
            "../results/quant/Tpms_Isoforms.csv",
            "../results/quant/Tpms_Isoforms.csv.ercc.csv",
        ],
        ids=[
            "../results/quant/IsoformToGene.csv.gz",
            "../results/quant/IdsToNames.csv",
            "../results/quant/IdsToNames_Isoforms.csv",
        ],
        velocity=[
            "../results/velocity/a.loom",
            "../results/velocity/a.obs_names.csv",
        ],
    output: directory("../results/newinputs/RNAData/")
    conda: "../envs/downloads.yaml"
    log: "../results/SingleCellProteogenomics_copyResults.log"
    benchmark: "../results/SingleCellProteogenomics_copyResults.benchmark"
    shell:
        "(mkdir -p {output} && cp {input.quant} {input.ids} {input.velocity} {output}) &> {log}"

rule downloadInputs:
    output: directory("../results/input")
    conda: "../../SingleCellProteogenomics/workflow/envs/download.yaml"
    log: "../results/SingleCellProteogenomics_downloadInputs.log"
    params: outdir=lambda w, output: os.path.dirname(output[0]),
    shell:
        "(gdown -O {output}.zip"
        f" \"https://drive.google.com/uc?export=download&id={INPUTID}\""
        " && unzip -d {params.outdir} {output}.zip) 2> {log}"

rule ProteinCellCycleClusters:
    input:
        "../results/input/",
        "../results/newinputs/RNAData/"
    output: "../results/output/pickles/mockbulk_phases.npy"
    conda: "../../SingleCellProteogenomics/workflow/envs/enviro.yaml"
    log: "../results/output/1_ProteinCellCycleClusters.log"
    threads: workflow.cores # use a whole node
    shell: "cd ../results && python ../SingleCellProteogenomics/1_ProteinCellCycleClusters.py &> {log}"

rule ProteinFucciPsuedotime:
    input: "../results/output/pickles/mockbulk_phases.npy"
    output: "../results/output/ProteinPseudotimePlotting.csv.gz"
    conda: "../../SingleCellProteogenomics/workflow/envs/enviro.yaml"
    log: "../results/output/2_ProteinFucciPsuedotime.log"
    threads: workflow.cores # use a whole node
    shell: "cd ../results && python ../SingleCellProteogenomics/2_ProteinFucciPsuedotime.py &> {log}"

rule RNAFucciPseudotime:
    input: "../results/output/ProteinPseudotimePlotting.csv.gz"
    output: "../results/output/RNAPseudotimePlotting.csv.gz"
    conda: "../../SingleCellProteogenomics/workflow/envs/enviro.yaml"
    log: "../results/output/3_RNAFucciPseudotime.log"
    threads: workflow.cores # use a whole node
    shell: "cd ../results && python ../SingleCellProteogenomics/3_RNAFucciPseudotime.py &> {log}"

rule TemporalDelay:
    input: "../results/output/RNAPseudotimePlotting.csv.gz"
    output: "../results/output/diff_max_pol.csv"
    conda: "../../SingleCellProteogenomics/workflow/envs/enviro.yaml"
    log: "../results/output/4_TemporalDelay.log"
    threads: workflow.cores # use a whole node
    shell: "cd ../results && python ../SingleCellProteogenomics/4_TemporalDelay.py &> {log}"

rule ProteinProperties:
    input: "../results/output/diff_max_pol.csv"
    output: "../results/output/upstreamKinaseResults.csv"
    conda: "../../SingleCellProteogenomics/workflow/envs/enviro.yaml"
    log: "../results/output/5_ProteinProperties.log"
    threads: workflow.cores # use a whole node
    shell: "cd ../results && python ../SingleCellProteogenomics/5_ProteinProperties.py &> {log}"

rule SingleCellProteogenomics_final:
    input:
        protein="../results/output/ProteinPseudotimePlotting.csv.gz",
        rna="../results/output/RNAPseudotimePlotting.csv.gz",
        kinase="../results/output/upstreamKinaseResults.csv" # prompt full run
    output:
        protein="../results/final/ProteinPseudotimePlotting.csv.gz",
        rna="../results/final/RNAPseudotimePlotting.csv.gz",
    conda: "../../SingleCellProteogenomics/workflow/envs/enviro.yaml"
    log: "../results/SingleCellProteogenomics_run.log"
    benchmark: "../results/SingleCellProteogenomics_run.benchmark"
    threads: 1
    shell:
        "(cp {input.protein} {output.protein} && "
        "cp {input.rna} {output.rna}) &> {log}"
