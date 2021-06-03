rule velocity_analysis:
    '''Run the velocity analysis'''
    input:
        unload_genome="../results/align/unloaded_2pass",
        gtf=f"{GTF}.fix.gtf",
        bams=expand("../results/align/{sra}Aligned.sortedByCoord.out.bam", sra=config['sra'])
    output: "../results/velocity/a.loom"
    conda: "../envs/velo.yaml"
    log: "../results/velocity.log"
    benchmark: "../results/velocity.benchmark"
    params:
        outfolder=lambda w, output: os.path.dirname(output[0]),
        sampleid="a" # used this historically, so just keeping it consistent
    shell:
        "velocyto run-smartseq2 --outputfolder {params.outfolder}"
        " --sampleid {params.sampleid} {input.bams} {input.gtf} &> {log}"

rule srr_lookup:
    '''Lookup the SRX accessions for the SRRs'''
    output: "../results/srr_lookup.txt"
    conda: "../envs/velo.yaml"
    log: "../results/srr_lookup.log"
    benchmark: "../results/srr_lookup.benchmark"
    params:
        commandpart1="esearch -db sra -query",
        commandpart2="| esummary | xtract -pattern DocumentSummary -element Run@acc Experiment@acc"
    shell: "(" + "\n".join(["{params.commandpart1} " + sra + " {params.commandpart2} >> {output}" for sra in config['sra']]) + ") 2> {log}"

rule make_aobs:
    '''Make the observation name file expected'''
    input:
        series_matrix="../resources/GSE146773_series_matrix.txt",
        srr_lookup="../results/srr_lookup.txt",
        loom="../results/velocity/a.loom"
    output: "../results/velocity/a.obs_names.csv"
    log: "../results/velocity/a.obs_names.log"
    benchmark: "../results/velocity/a.obs_names.benchmark"
    conda: "../envs/velo.yaml"
    shell: "python scripts/make_a_obs.py {input.loom} {input.srr_lookup} {input.series_matrix} {output} 2> {log}"
