RAM_MB_REQ = 50000 #mb
RAM_B_REQ = RAM_MB_REQ * 1000 # BYTES

rule star_genome_generate:
    input:
        efa="data/ERCC.fa",
        gfa=FA,
        gff=GFF3 + ".fix.gff3"
    output: STAR_REF_FOLDER + "/SA"
    params: genomedir=STAR_REF_FOLDER
    threads: 99
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genomedir} "
        "--genomeFastaFiles {input.efa} {input.gfa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"

rule load_star_genome_firstpass:
    input: STAR_REF_FOLDER + "/SA"
    output: temp("output/loaded_firstpass")
    params: genomedir=STAR_REF_FOLDER
    resources: mem_mb = RAM_MB_REQ
    shell: "STAR --genomeLoad LoadAndExit --genomeDir {params} && touch {output}"

rule star_firstpass:
    input:
        temp("output/loaded_firstpass"),
        suffix=STAR_REF_FOLDER + "/SA",
        fastq=FQ_FOLDER + "{fq}.fastq.gz" #single end only
    output: "output/{fq}SJ.out.tab"
    threads: 6
    params:
        junctions="--outFilterIntronMotifs RemoveNoncanonical", # adds XS tag to all alignments that contain a splice junction
        bam="--outSAMtype None",
        genomedir=STAR_REF_FOLDER
    shell:
        "STAR --genomeLoad LoadAndKeep --genomeDir {params.genomedir}"
        " --runThreadN {threads} {params.junctions} {params.bam}"
        " --outFileNamePrefix output/{wildcards.fq}"
        " --readFilesIn <(zcat " + FQ_FOLDER + "{wildcards.fq}.fastq.gz)"

rule unload_firstpass_genome:
    input:
        STAR_REF_FOLDER + "/SA",
        jj=expand("output/{fq}SJ.out.tab", fq=FQ_PREFIXES)
    output: temp("temp/unloaded_firstpass")
    params: genomedir=STAR_REF_FOLDER
    shell:
        "STAR --genomeLoad Remove --genomeDir {params.genomedir} && "
        "touch {output}"

rule process_first_pass:
    input:
        "temp/unloaded_firstpass",
        jj=expand("output/{fq}SJ.out.tab", fq=FQ_PREFIXES)
    output: "output/combinedSJ.out.tab"
    shell:
        "awk 'BEGIN {{OFS=\"\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";}} "
        "{{if($5>0){{print $1,$2,$3,strChar[$4]}}}}' {input.jj} | "
        "grep -v 'MT' >> {output}"

rule star_genome_generate_secondpass:
    input:
        efa="data/ERCC.fa",
        gfa=FA,
        gff=GFF3 + ".fix.gff3",
        jj="output/combinedSJ.out.tab"
    output: STAR_REF_FOLDER + unique_tag() + "/SA"
    params: genomedir=STAR_REF_FOLDER + unique_tag()
    threads: 99
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genomedir} "
        "--genomeFastaFiles {input.efa} {input.gfa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"
        "--limitSjdbInsertNsj 1200000 --sjdbFileChrStartEnd {input.jj}"

rule load_star_genome_2pass:
    input: STAR_REF_FOLDER + unique_tag() + "/SA"
    output: temp("output/loaded_2pass")
    params: genomedir=STAR_REF_FOLDER + unique_tag()
    resources: mem_mb = RAM_MB_REQ
    shell: "STAR --genomeLoad LoadAndExit --genomeDir {params} && touch {output}"

rule star_2pass:
    input:
        temp("output/loaded_2pass"),
        genomedir=STAR_REF_FOLDER + unique_tag() + "/SA",
        fastq=FQ_FOLDER + "{fq}.fastq.gz" #single end only
    output: "output/{fq}" + unique_tag() + "Aligned.sortedByCoord.out.bam"
    threads: 6
    params:
        junctions="--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical", # adds XS tag to all alignments that contain a splice junction
        bam="--outSAMtype BAM SortedByCoordinate --outBAMcompression 10 --limitBAMsortRAM " + str(RAM_B_REQ),
        gatk="--outSAMattrRGline ID:1 PU:platform  PL:illumina SM:sample LB:library --outSAMmapqUnique 60",
        genomedir=STAR_REF_FOLDER + unique_tag()
    shell:
        "STAR --runMode alignReads --genomeLoad LoadAndKeep --genomeDir {params.genomedir}"
        " --runThreadN {threads} {params.junctions} {params.bam} {params.gatk}"
        " --outFileNamePrefix output/{wildcards.fq}" + unique_tag() +
        " --readFilesIn <(zcat " + FQ_FOLDER + "{wildcards.fq}.fastq.gz)"

rule unload_2pass_genome:
    input:
        STAR_REF_FOLDER + unique_tag() + "/SA",
        expand("output/{fq}" + unique_tag() + "Aligned.sortedByCoord.out.bam", fq=FQ_PREFIXES)
    output: temp("temp/unloaded_2pass")
    params: genomedir=STAR_REF_FOLDER + unique_tag()
    shell:
        "STAR --genomeLoad Remove --genomeDir {params.genomedir} && "
        "touch {output}"
