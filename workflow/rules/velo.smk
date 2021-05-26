RAM_MB_REQ = 50000 #mb
RAM_B_REQ = RAM_MB_REQ * 1000 # BYTES

rule star_genome_generate:
    input:
        efa="../resources/ERCC.filtered.fa",
        gfa=GENOME_FA,
        gff=f"{ENSEMBL_GFF}.fix.gff3"
    output: f"{STAR_REF_FOLDER}/SA"
    params: genomedir=lambda w, output: os.path.dirname(output[0])
    threads: 99
    log: f"{STAR_REF_FOLDER}.genome_generate.log"
    benchmark: f"{STAR_REF_FOLDER}.genome_generate.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genomedir}"
        " --genomeFastaFiles {input.efa} {input.gfa} --sjdbGTFfile {input.gff}"
        " --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 &> {log}"

rule load_star_genome_firstpass:
    input: f"{STAR_REF_FOLDER}/SA"
    output: temp("../results/loaded_firstpass")
    params: genomedir=lambda w, input: os.path.dirname(input[0])
    resources: mem_mb = RAM_MB_REQ
    log: f"{STAR_REF_FOLDER}.genome_load.log"
    benchmark: f"{STAR_REF_FOLDER}.genome_load.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --genomeLoad LoadAndExit --genomeDir {params.genomedir} && "
        "touch {output}) &> {log}"

rule star_firstpass:
    input:
        loadedgenome="../results/loaded_firstpass",
        suffix=f"{STAR_REF_FOLDER}/SA",
        fastq="../results/fastq/{sra}.trim.fastq.gz", #single end only
    output: "../results/align/SJ1st/{sra}SJ.out.tab"
    threads: 4
    params:
        outfolder=lambda w, output: os.path.dirname(output[0]),
        junctions="--outFilterIntronMotifs RemoveNoncanonical", # adds XS tag to all alignments that contain a splice junction
        bam="--outSAMtype None",
        sjfilter=" --outSJfilterReads Unique", # for JUM
        genomedir=lambda w, input: os.path.dirname(input.suffix)
    log: "../results/align/SJ1st/{sra}SJ.firstpass.log"
    benchmark: "../results/align/SJ1st/{sra}SJ.firstpass.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --genomeLoad LoadAndKeep --genomeDir {params.genomedir}"
        " --runThreadN {threads} {params.junctions} {params.bam} {params.sjfilter}"
        " --outFileNamePrefix {params.outfolder}/{wildcards.sra}"
        " --readFilesIn <(zcat {input.fastq})) &> {log}"

rule unload_firstpass_genome:
    input:
        f"{STAR_REF_FOLDER}/SA",
        jj=expand("../results/align/SJ1st/{sra}SJ.out.tab", sra=config['sra'])
    output: temp("temp/unloaded_firstpass")
    params: genomedir=lambda w, input: os.path.dirname(input[0])
    log: f"{STAR_REF_FOLDER}.genome_unload.log"
    benchmark: f"{STAR_REF_FOLDER}.genome_unload.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --genomeLoad Remove --genomeDir {params.genomedir} && "
        "touch {output}) &> {log}"

rule star_genome_generate_secondpass:
    input:
        efa="../resources/ERCC.filtered.fa",
        gfa=GENOME_FA,
        gff=f"{ENSEMBL_GFF}.fix.gff3",
        jj=expand("../results/align/SJ1st/{sra}SJ.out.tab", sra=config['sra'])
    output: f"{STAR_REF_FOLDER}SecondPass/SA"
    params: genomedir=lambda w, output: os.path.dirname(output[0])
    log: f"{STAR_REF_FOLDER}SecondPass.genome_generate.log"
    benchmark: f"{STAR_REF_FOLDER}SecondPass.genome_generate.benchmark"
    conda: "../envs/quant.yaml"
    threads: 99
    shell:
        "(STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genomedir}"
        " --genomeFastaFiles {input.efa} {input.gfa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"
        " --limitSjdbInsertNsj 1200000 --sjdbFileChrStartEnd {input.jj}) &> {log}"

rule load_star_genome_2pass:
    input: f"{STAR_REF_FOLDER}SecondPass/SA"
    output: temp("../results/loaded_2pass")
    params: genomedir=lambda w, input: os.path.dirname(input[0])
    resources: mem_mb=RAM_MB_REQ
    log: f"{STAR_REF_FOLDER}SecondPass.genome_load.log"
    benchmark: f"{STAR_REF_FOLDER}SecondPass.genome_load.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --genomeLoad LoadAndExit --genomeDir {params.genomedir} && "
        "touch {output}) &> {log}"

rule star_2pass:
    input:
        "../results/loaded_2pass",
        genome_sa=f"{STAR_REF_FOLDER}SecondPass/SA",
        fastq="../results/fastq/{sra}.trim.fastq.gz" #single end only
    output:
        sj="../results/align/{sra}SJ.out.tab",
        bam="../results/align/{sra}Aligned.sortedByCoord.out.bam",
        log="../results/align/{sra}Log.out",
        final="../results/align/{sra}Log.final.out",
        progress="../results/align/{sra}Log.progress.out",
    threads: 6
    params:
        junctions="--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical", # adds XS tag to all alignments that contain a splice junction
        bam=f"--outSAMtype BAM SortedByCoordinate --outBAMcompression 10 --limitBAMsortRAM {str(RAM_B_REQ)}",
        gatk="--outSAMattrRGline ID:1 PU:platform  PL:illumina SM:sample LB:library --outSAMmapqUnique 60",
        sjfilter="--outSJfilterReads Unique", # for JUM
        genomedir=lambda w, input: os.path.dirname(input.genome_sa),
        outfolder=lambda w, output: os.path.dirname(output.bam)
    log: "../results/align/{sra}.secondpass.log"
    benchmark: "../results/align/{sra}.secondpass.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --runMode alignReads --genomeLoad LoadAndKeep --genomeDir {params.genomedir} {params.sjfilter}"
        " --runThreadN {threads} {params.junctions} {params.bam} {params.gatk}"
        " --outFileNamePrefix {params.outfolder}/{wildcards.sra}"
        " --readFilesIn <(zcat {input.fastq})) &> {log}"

rule unload_2pass_genome:
    input:
        f"{STAR_REF_FOLDER}SecondPass/SA",
        bam=expand("../results/align/{sra}Aligned.sortedByCoord.out.bam", sra=config['sra']),
    output: temp("temp/unloaded_2pass")
    params: genomedir=lambda w, input: os.path.dirname(input[0])
    log: f"{STAR_REF_FOLDER}SecondPass.genome_unload.log"
    benchmark: f"{STAR_REF_FOLDER}SecondPass.genome_unload.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --genomeLoad Remove --genomeDir {params.genomedir} && "
        "touch {output}) &> {log}"

rule velocity_analysis:
    '''Run the velocity analysis'''
    input:
        "temp/unloaded_2pass",
        gtf=ENSEMBL_GTF,
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
        " --sampleid {params.sampleid} {input.bams} {input.gtf}"
