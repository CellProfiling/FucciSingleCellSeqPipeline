rule rsem_star_genome:
    '''Create an RSEM reference with STAR indices'''
    input:
        efa="../../resources/ERCC.fa",
        gfa=ENSEMBL_FA,
        gff=ENSEMBL_GFF + ".fix.gff3"
    output:
        REFSTAR_PREFIX + ".gtf",
        suffix = REFSTAR_FOLDER + "SA"
    threads: 99
    resources: mem_mb=60000
    log: "../../resources/ensembl/prepare-reference.log"
    conda: "../envs/environment.yaml"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star --gff3 {input.gff}"
        " \"{input.efa}\",\"{input.gfa}\" " + REFSTAR_PREFIX +
        ") 2> {log}"

rule rsem_star_align:
    '''Align to transcripts with STAR and quantify with RSEM'''
    input:
        suffix=REFSTAR_FOLDER + "SA",
        gtf=REFSTAR_PREFIX + ".gtf",
        fastq=expand("../../results/fastq/{sra}.trim.fastq.gz", sra=config['sra']) #single end only
    output:
        "../../results/quant/{sra}.isoforms.results",
        "../../results/quant/{sra}.genes.results",
        "../../results/quant/{sra}.time",
        directory("../../results/{sra}.stat"),
    resources: mem_mb=50000
    threads: 12
    log: "../../results/{sra}calculate-expression.log"
    shell:
        "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" probably not using confidence intervals here
        " --num-threads {threads} <(zcat ../../results/{wildcards.sra}.trim.fastq.gz) " + REFSTAR_PREFIX + " ../../results/quant/{wildcards.sra}) &> {log}"

rule make_gene_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("../../results/quant/{sra}.genes.results", sra=config['sra']),
        gff=ENSEMBL_GFF + ".fix.gff3"
    output:
        counts="../../results/quant/Counts.csv",
        names="../../results/quant/IdsToNames.csv",
        tpms="../../results/quant/Tpms.csv"
    shell:
        "python scripts/make_rsem_dataframe.py genes {input.gff} {output.counts} {output.tpms} {output.names}"

rule make_isoform_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("../../results/quant/{sra}.isoforms.results", sra=config['sra']),
        gff=ENSEMBL_GFF + ".fix.gff3"
    output:
        counts="../../results/quant/Counts_Isoforms.csv",
        names="../../results/quant/IdsToNames_Isoforms.csv",
        tpms="../../results/quant/Tpms_Isoforms.csv"
    shell:
        "python scripts/make_rsem_dataframe.py isoforms {input.gff} {output.counts} {output.tpms} {output.names}"
