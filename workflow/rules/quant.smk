rule rsem_star_genome:
    '''Create an RSEM reference with STAR indices'''
    input:
        efa="../resources/ERCC.fa",
        gfa=GENOME_FA,
        gff=f"{ENSEMBL_GFF}.fix.gff3"
    output:
        f"{REFSTAR_PREFIX}.gtf",
        suffix = f"{REFSTAR_FOLDER}SA"
    params: prefix=lambda w, output: output.suffix.strip("SA")
    threads: 99
    resources: mem_mb=60000
    log: "../resources/ensembl/prepare-reference.log"
    conda: "../envs/quant.yaml"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star --gff3 {input.gff}"
        " \"{input.efa}\",\"{input.gfa}\" {params.prefix}) 2> {log}"

rule rsem_star_align:
    '''Align to transcripts with STAR and quantify with RSEM'''
    input:
        suffix=f"{REFSTAR_FOLDER}SA",
        gtf=f"{REFSTAR_PREFIX}.gtf",
        fastq=expand("../results/fastq/{sra}.trim.fastq.gz", sra=config['sra']) #single end only
    output:
        "../results/quant/{sra}.isoforms.results",
        "../results/quant/{sra}.genes.results",
        "../results/quant/{sra}.time",
        directory("../results/{sra}.stat"),
    params: prefix=lambda w, input: input.suffix.strip("SA")
    resources: mem_mb=50000
    threads: 16
    log: "../results/{sra}calculate-expression.log"
    benchmark: "../results/{sra}calculate-expression.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" probably not using confidence intervals here
        " --num-threads {threads} <(zcat ../../results/{wildcards.sra}.trim.fastq.gz)"
        " {params.prefix}"
        " ../../results/quant/{wildcards.sra}) &> {log}"

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
        "python scripts/make_rsem_dataframe.py genes {input.gff} {output.counts} {output.tpms} {output.names}"

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
        "python scripts/make_rsem_dataframe.py isoforms {input.gff} {output.counts} {output.tpms} {output.names}"
