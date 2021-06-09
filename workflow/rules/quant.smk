rule filter_ercc:
    '''Filter the empty ERCC entries'''
    input: "../resources/ERCC.fa"
    output: "../resources/ERCC.filtered.fa"
    log: "../resources/ERCC.filtered.log"
    benchmark: "../resources/ERCC.filtered.benchmark"
    conda: "../envs/quant.yaml"
    shell: "grep -v ERCC_ID {input} > {output} 2> {log}"

rule rsem_reference:
    '''Create an RSEM reference'''
    input:
        efa="../resources/ERCC.filtered.fa",
        gfa=FA,
        gtf=f"{GTF}.fix.gtf"
    output:
        grp=f"{RSEM_REF_FOLDER}RsemReference.grp",
        fa=f"{RSEM_REF_FOLDER}RsemReference.transcripts.fa",
        suffix=f"{RSEM_REF_FOLDER}SA"
    params: prefix=lambda w, output: output.grp[:-4]
    threads: 99 # can do fewer without STAR references
    log: "../resources/ensembl/prepare-reference.log"
    benchmark: "../resources/ensembl/prepare-reference.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
        " --star" # TODO: figure out how to use transcript STAR output
        " \"{input.efa}\",\"{input.gfa}\" {params.prefix}) &> {log}"

rule rsem_calculations:
    '''Align to transcripts with STAR and quantify with RSEM'''
    input:
        grp=f"{RSEM_REF_FOLDER}RsemReference.grp",
        suffix=f"{RSEM_REF_FOLDER}SA",
        #bam="../results/align/{sra}Aligned.toTranscriptome.out.bam", #single end only
        fq="../results/fastq/{sra}.trim.fastq.gz"
    output:
        "../results/quant/{sra}.isoforms.results",
        "../results/quant/{sra}.genes.results",
        directory("../results/quant/{sra}.stat"),
    params:
        prefix=lambda w, input: input.grp[:-4],
        results=lambda w, output: os.path.dirname(os.path.dirname(output[0])),
        quant=lambda w, output: os.path.dirname(output[0]),
        options="--no-bam-output --single-cell-prior" # --bam"
    resources: mem_mb=30000
    threads: 8
    log: "../results/quant/{sra}calculate-expression.log"
    benchmark: "../results/quant/{sra}calculate-expression.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(rsem-calculate-expression --num-threads {threads}"
        " --star"
        " {params.options} <(zcat {input.fq})" #{input.bam}"
        " {params.prefix}"
        " {params.quant}/{wildcards.sra}) &> {log}"

rule make_gene_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("../results/quant/{sra}.genes.results", sra=config['sra']),
        gff=f"{GFF}.fix.gff3",
        series_matrix="../resources/GSE146773_series_matrix.txt",
        srr_lookup="../results/srr_lookup.txt",
    output:
        counts="../results/quant/Counts.csv",
        counts_ercc="../results/quant/Counts.csv.ercc.csv",
        names="../results/quant/IdsToNames.csv",
        tpms="../results/quant/Tpms.csv",
        tpms_ercc="../results/quant/Tpms.csv.ercc.csv",
        ids="../results/quant/IsoformToGene.csv.gz",
    conda: "../envs/quant.yaml"
    log: "../results/quant/Counts.log"
    benchmark: "../results/quant/Counts.benchmark"
    shell:
        "python scripts/make_rsem_dataframe.py genes {input.gff} {input.srr_lookup} {input.series_matrix}"
        " {output.counts} {output.tpms} {output.names} {output.ids} &> {log}"

rule make_isoform_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("../results/quant/{sra}.isoforms.results", sra=config['sra']),
        gff=f"{GFF}.fix.gff3",
        series_matrix="../resources/GSE146773_series_matrix.txt",
        srr_lookup="../results/srr_lookup.txt",
    output:
        counts="../results/quant/Counts_Isoforms.csv",
        counts_ercc="../results/quant/Counts_Isoforms.csv.ercc.csv",
        names="../results/quant/IdsToNames_Isoforms.csv",
        tpms="../results/quant/Tpms_Isoforms.csv",
        tpms_ercc="../results/quant/Tpms_Isoforms.csv.ercc.csv",
        ids="../results/quant/IsoformToGene_Isoforms.csv.gz",
    conda: "../envs/quant.yaml"
    log: "../results/quant/Counts_Isoforms.log"
    benchmark: "../results/quant/Counts_Isoforms.benchmark"
    shell:
        "python scripts/make_rsem_dataframe.py isoforms {input.gff} {input.srr_lookup} {input.series_matrix}"
        " {output.counts} {output.tpms} {output.names} {output.ids} &> {log}"
