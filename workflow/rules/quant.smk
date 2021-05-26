rule filter_ercc:
    '''Filter the empty ERCC entries'''
    input: "../resources/ERCC.fa"
    output: "../resources/ERCC.filtered.fa"
    log: "../resources/ERCC.filtered.log"
    benchmark: "../resources/ERCC.filtered.benchmark"
    conda: "../envs/quant.yaml"
    shell: "grep -v ERCC_ID {input} > {output} 2> {log}"

rule rsem_star_genome:
    '''Create an RSEM reference with STAR indices'''
    input:
        efa="../resources/ERCC.filtered.fa",
        gfa=GENOME_FA,
        gff=f"{ENSEMBL_GFF}.fix.gff3"
    output:
        gtf=f"{REFSTAR_FOLDER}RsemStarReference.gtf",
        suffix = f"{REFSTAR_FOLDER}SA"
    params:
        prefix=lambda w, output: output.gtf[:-4]
    threads: 99
    resources: mem_mb=60000
    log: "../resources/ensembl/prepare-reference.log"
    benchmark: "../resources/ensembl/prepare-reference.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star --gff3 {input.gff}"
        " \"{input.efa}\",\"{input.gfa}\" {params.prefix}) &> {log}"

rule rsem_star_align:
    '''Align to transcripts with STAR and quantify with RSEM'''
    input:
        suffix=f"{REFSTAR_FOLDER}SA",
        gtf=f"{REFSTAR_FOLDER}RsemStarReference.gtf",
        fastq=expand("../results/fastq/{sra}.trim.fastq.gz", sra=config['sra']) #single end only
    output:
        "../results/quant/{sra}.isoforms.results",
        "../results/quant/{sra}.genes.results",
        "../results/quant/{sra}.time",
        directory("../results/quant/{sra}.stat"),
    params:
        prefix=lambda w, input: input.gtf[:-4],
        results=lambda w, output: os.path.dirname(os.path.dirname(output[0])),
        quant=lambda w, output: os.path.dirname(output[0]),
    resources: mem_mb=30000
    threads: 8
    log: "../results/quant/{sra}calculate-expression.log"
    benchmark: "../results/quant/{sra}calculate-expression.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(rsem-calculate-expression --time --no-bam-output --star" # --calc-ci" probably not using confidence intervals here
        " --num-threads {threads} <(zcat {params.results}/{wildcards.sra}.trim.fastq.gz)"
        " {params.prefix}"
        " {params.quant}/{wildcards.sra}) &> {log}"

rule make_gene_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("../results/quant/{sra}.genes.results", sra=config['sra']),
        gff=f"{ENSEMBL_GFF}.fix.gff3"
    output:
        counts="../results/quant/Counts.csv",
        names="../results/quant/IdsToNames.csv",
        tpms="../results/quant/Tpms.csv"
    conda: "../envs/quant.yaml"
    log: "../results/quant/Counts.log"
    benchmark: "../results/quant/Counts.benchmark"
    shell:
        "python scripts/make_rsem_dataframe.py genes {input.gff} "
        "{output.counts} {output.tpms} {output.names} &> {log}"

rule make_isoform_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("../results/quant/{sra}.isoforms.results", sra=config['sra']),
        gff=f"{ENSEMBL_GFF}.fix.gff3"
    output:
        counts="../results/quant/Counts_Isoforms.csv",
        names="../results/quant/IdsToNames_Isoforms.csv",
        tpms="../results/quant/Tpms_Isoforms.csv"
    conda: "../envs/quant.yaml"
    log: "../results/quant/Counts_Isoforms.log"
    benchmark: "../results/quant/Counts_Isoforms.benchmark"
    shell:
        "python scripts/make_rsem_dataframe.py isoforms {input.gff} "
        "{output.counts} {output.tpms} {output.names} &> {log}"
