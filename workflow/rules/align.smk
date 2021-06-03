RAM_MB_REQ = 50000 #mb
RAM_B_REQ = RAM_MB_REQ * 1000 # BYTES

rule star_genome_generate:
    input:
        efa="../resources/ERCC.filtered.fa",
        gfa=FA,
        gtf=f"{GTF}.fix.gtf"
    output: f"{STAR_REF_FOLDER}/SA"
    params:
        genomedir=lambda w, output: os.path.dirname(output[0]),
        size="--genomeSAindexNbases 12" if FA == TEST_GENOME_FA else "",
        sjoptions="--sjdbOverhang 42" # 43bp reads for this experiment
    threads: 99
    log: f"{STAR_REF_FOLDER}.genome_generate.log"
    benchmark: f"{STAR_REF_FOLDER}.genome_generate.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genomedir}"
        " --genomeFastaFiles {input.efa} {input.gfa} --sjdbGTFfile {input.gtf}"
        " {params.sjoptions} {params.size}) &> {log}"

rule load_star_genome_firstpass:
    input: suffix=f"{STAR_REF_FOLDER}/SA"
    output: temp("../results/align/loaded_firstpass")
    params: genomedir=lambda w, input: os.path.dirname(input.suffix)
    resources: mem_mb = RAM_MB_REQ
    log: "../results/align/loaded_firstpass.log"
    benchmark: "../results/align/loaded_firstpass.benchmark"
    conda: "../envs/quant.yaml"
    shell: "(STAR --genomeLoad LoadAndExit --genomeDir {params} && touch {output}) &> {log}"

rule star_firstpass:
    input:
        "../results/align/loaded_firstpass",
        suffix=f"{STAR_REF_FOLDER}/SA",
        fastq="../results/fastq/{sra}.trim.fastq.gz",
    output: sj="../results/align/SJ1st/{sra}SJ.out.tab"
    threads: 8
    params:
        bam="--outSAMtype None",
        genomedir=lambda w, input: os.path.dirname(input.suffix),
        outfolder=lambda w, output: os.path.dirname(output.sj)
    log: "../results/align/SJ1st/{sra}SJ.log"
    benchmark: "../results/align/SJ1st/{sra}SJ.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --genomeLoad LoadAndKeep --genomeDir {params.genomedir}"
        " --runThreadN {threads} {params.bam} "
        " --outFileNamePrefix {params.outfolder}/{wildcards.sra}"
        " --readFilesIn <(zcat {input.fastq}) ) &> {log}"

rule unload_firstpass_genome:
    input:
        suffix=f"{STAR_REF_FOLDER}/SA",
        jj=expand("../results/align/SJ1st/{sra}SJ.out.tab", sra=config['sra'])
    output: temp("../results/align/unloaded_firstpass")
    params: genomedir=lambda w, input: os.path.dirname(input.suffix)
    log: "../results/align/unloaded_firstpass.log"
    benchmark: "../results/align/unloaded_firstpass.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --genomeLoad Remove --genomeDir {params.genomedir} && "
        "touch {output}) &> {log}"

rule star_genome_generate_secondpass:
    input:
        unload="../results/align/unloaded_firstpass",
        efa="../resources/ERCC.fa",
        gfa=FA,
        gtf=f"{GTF}.fix.gtf",
        jj=expand("../results/align/SJ1st/{sra}SJ.out.tab", sra=config['sra'])
    output: suffix=f"{STAR_REF_FOLDER}SecondPass/SA"
    params:
        genomedir=lambda w, output: os.path.dirname(output.suffix),
        size="--genomeSAindexNbases 12" if FA == TEST_GENOME_FA else "",
        sjoptions="--sjdbOverhang 42 --limitSjdbInsertNsj 1200000" # 43bp reads for this experiment
    threads: 99
    log: f"{STAR_REF_FOLDER}SecondPass.generate.log"
    benchmark: f"{STAR_REF_FOLDER}SecondPass.generate.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genomedir}"
        " --genomeFastaFiles {input.efa} {input.gfa} --sjdbGTFfile {input.gtf} {params.sjoptions}"
        " --sjdbFileChrStartEnd {input.jj} {params.size}) &> {log}"

rule load_star_genome_2pass:
    input: suffix=f"{STAR_REF_FOLDER}SecondPass/SA"
    output: temp("../results/align/loaded_2pass")
    params: genomedir=lambda w, input: os.path.dirname(input.suffix),
    resources: mem_mb = RAM_MB_REQ
    log: "../results/align/loaded_2pass.log"
    benchmark: "../results/align/loaded_2pass.benchmark"
    conda: "../envs/quant.yaml"
    shell: "(STAR --genomeLoad LoadAndExit --genomeDir {params} && touch {output}) &> {log}"

rule star_2pass:
    input:
        "../results/align/loaded_2pass",
        suffix=f"{STAR_REF_FOLDER}SecondPass/SA",
        fastq="../results/fastq/{sra}.trim.fastq.gz"
    output:
        sj="../results/align/{sra,[A-Z,0-9]+}SJ.out.tab",
        bam="../results/align/{sra}Aligned.sortedByCoord.out.bam",
        transcript="../results/align/{sra}Aligned.toTranscriptome.out.bam",
        log="../results/align/{sra}Log.out",
        final="../results/align/{sra}Log.final.out",
        progress="../results/align/{sra}Log.progress.out",
    wildcard_constraints: sra="[A-Z0-9]+"
    threads: 8
    params:
        mode="--runMode alignReads --genomeLoad LoadAndKeep",
        bam=f"--outSAMtype BAM SortedByCoordinate --outBAMcompression 10 --limitBAMsortRAM {str(RAM_B_REQ)}",
        transcripts=f"--quantMode TranscriptomeSAM",
        outfolder=lambda w, output: os.path.dirname(output.sj),
        genomedir=lambda w, input: os.path.dirname(input.suffix),
    log: "../results/align/{sra}Align.log"
    benchmark: "../results/align/{sra}Align.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR {params.mode} --genomeDir {params.genomedir}"
        " --runThreadN {threads} {params.bam} {params.transcripts}"
        " --outFileNamePrefix {params.outfolder}/{wildcards.sra}"
        " --readFilesIn <(zcat {input.fastq}) ) &> {log}"

rule unload_2pass_genome:
    input:
        suffix=f"{STAR_REF_FOLDER}SecondPass/SA",
        bam=expand("../results/align/{sra}Aligned.sortedByCoord.out.bam", sra=config['sra']),
    output: temp("../results/align/unloaded_2pass")
    params: genomedir=lambda w, input: os.path.dirname(input.suffix)
    log: "../results/align/unloaded_2pass.log"
    benchmark: "../results/align/unloaded_2pass.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "(STAR --genomeLoad Remove --genomeDir {params.genomedir} && "
        "touch {output}) &> {log}"
