rule rsem_star_genome:
    '''Create an RSEM reference with STAR indices'''
    input:
        efa="data/ERCC.fa",
        gfa=FA,
        gff=GFF3 + ".fix.gff3"
    output:
        REFSTAR_PREFIX + ".gtf",
        suffix = REFSTAR_FOLDER + "SA"
    threads: 99
    resources: mem_mb=60000
    log: "ensembl/prepare-reference.log"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star --gff3 {input.gff}"
        " \"{input.efa}\",\"{input.gfa}\" " + REFSTAR_PREFIX +
        ") 2> {log}"

rule rsem_star_align:
    '''Align to transcripts with STAR and quantify with RSEM'''
    input:
        suffix=REFSTAR_FOLDER + "SA",
        gtf=REFSTAR_PREFIX + ".gtf",
        fastq=expand(FQ_FOLDER + "{fq}.fastq.gz", fq=FQ_PREFIXES) #single end only
    output:
        "output/{fq}.isoforms.results",
        "output/{fq}.genes.results",
        "output/{fq}.time",
        directory("output/{fq}.stat"),
    resources: mem_mb=50000
    threads: 12
    log: "output/{fq}calculate-expression.log"
    shell:
        "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" probably not using confidence intervals here
        " --num-threads {threads} <(zcat " + FQ_FOLDER + "{wildcards.fq}.fastq.gz) " + REFSTAR_PREFIX + " output/{wildcards.fq}) &> {log}"

rule make_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("output/{fq}.genes.results", fq=FQ_PREFIXES),
        gff=GFF3 = ".fix.gff3"
    output:
        counts="output/Counts.csv",
        names="output/IdsToNames.csv",
        tpms="output/Tpms.csv"
    shell:
        "python scripts/make_rsem_dataframe.py {input.gff} {output.counts} {output.tpms} {output.names}"

# rule rsem_genome:
#     '''Create an RSEM reference. For use with STAR aligned BAMs.'''
#     input:
#         efa="data/ERCC.fa",
#         gfa=FA,
#         gff=GFF3 + ".fix.gff3"
#     output:
#         REF_PREFIX + ".gtf",
#     threads: 12
#     resources: mem_mb=60000
#     log: "ensembl/prepare-reference.log"
#     shell:
#         "(rsem-prepare-reference --num-threads {threads} --gff3 {input.gff}"
#         " \"{input.efa}\",\"{input.gfa}\" " + REF_PREFIX +
#         ") 2> {log}"
#
# rule rsem_bam:
#     input:
#         REF_PREFIX + ".gtf",
#         temp("temp/unloaded_2pass"), # trigger unload star reference
#         "output/{fq}" + unique_tag() + "Aligned.sortedByCoord.out.bam"
#     output:
#         "output/{fq}" + unique_tag() + ".isoforms.results",
#         "output/{fq}" + unique_tag() + ".genes.results",
#         "output/{fq}" + unique_tag() + ".time",
#         directory("output/{fq}.stat"),
#     resources: mem_mb=50000
#     threads: 12
#     log: "output/{fq}calculate-expression.log"
#     shell:
#         "(rsem-calculate-expression --no-bam-output --time" # --calc-ci" probably not using confidence intervals here
#         " --num-threads {threads} --alignments output/{wildcards.fq}" + unique_tag() + "Aligned.sortedByCoord.out.bam " +
#         REF_PREFIX + " output/{wildcards.fq}) &> {log}"

# rule make_rsem_dataframe:
#     input: expand("output/{fq}" + unique_tag() + ".genes.results", fq=FQ_PREFIXES)
#     output:
#         counts="output/counts" + unique_tag() + ".tsv",
#         tpm="output/tpm" + unique_tag() + ".tsv"
#     shell:
#         "echo \"python scripts/make_rsem_dataframe.py 4 {output.counts} {input}\" && "
#         "python scripts/make_rsem_dataframe.py 4 {output.counts} {input} && "
#         "python scripts/make_rsem_dataframe.py 5 {output.tpm} {input}"
