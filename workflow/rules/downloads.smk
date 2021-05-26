PROTOCOL = "http"

rule download_ensembl_references:
    output:
        gfa=GENOME_FA,
        gff3=ENSEMBL_GFF,
        pfa=f"../resources/ensembl/{REF}.pep.all.fa",
    params:
        primary=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/dna/{REF}.dna.primary_assembly.fa.gz",
        toplevel=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/dna/{REF}.dna.toplevel.fa.gz",
        gff=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/gff3/{SPECIES_LOWER}/{REF}.{ENSEMBL_VERSION}.gff3.gz",
        pep=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/pep/{REF}.pep.all.fa.gz",
    benchmark: "../resources/ensembl/downloads.benchmark"
    log: "../resources/ensembl/downloads.log"
    conda: "../envs/downloads.yaml"
    shell:
        "((wget -O - {params.primary} || wget -O - {params.toplevel}) | gunzip -c - > {output.gfa} && "
        "wget -O - {params.gff} | gunzip -c - > {output.gff3} && "
        "wget -O - {params.pep} | gunzip -c - > {output.pfa}) 2> {log}"

rule fix_gff3_for_rsem:
    '''This script changes descriptive notes in column 4 to "gene" if a gene row, and it also adds ERCCs to the gene model'''
    input: ENSEMBL_GFF
    output: f"{ENSEMBL_GFF}.fix.gff3"
    log: f"{ENSEMBL_GFF}.fix.gff3.log"
    benchmark: f"{ENSEMBL_GFF}.fix.gff3.benchmark"
    conda: "../envs/downloads.yaml"
    shell: "python scripts/fix_gff3_for_rsem.py {input} {output} &> {log}"

rule prefetch_sras_se:
    '''Prefetch SRA from GEO SRA'''
    output: temp("../results/sra/{sra}/{sra}.sra")
    benchmark: "../results/sra/{sra}/{sra}.benchmark"
    log: "../results/sra/{sra}/{sra}.log"
    params: outdir=lambda w, output: os.path.dirname(os.path.dirname(output[0]))
    conda: "../envs/downloads.yaml"
    shell:
        "prefetch {wildcards.sra}"
        " --output-directory {params.outdir} &> {log}"

rule split_sras_se:
    input: "../results/sra/{sra}/{sra}.sra"
    output: "../results/fastq/{sra}.fastq" # independent of pe/se
    benchmark: "../results/fastq/{sra}.benchmark"
    log: "../results/fastq/{sra}.log"
    params: outdir=lambda w, output: os.path.dirname(output[0])
    conda: "../envs/downloads.yaml"
    shell:
        "fastq-dump -I --outdir {params.outdir} --split-files {input} && "
        "mv {params.outdir}/{wildcards.sra}_1.fastq {output} &> {log}"

rule fastp_sra_se:
    '''Trim adapters, read quality filtering, make QC outputs'''
    input: "../results/fastq/{sra}.fastq",
    output:
        fq="../results/fastq/{sra}.trim.fastq.gz",
        html="../results/fastq/{sra}.trim.html",
        json="../results/fastq/{sra}.trim.json",
    threads: 2
    log: "../results/fastq/{sra}.trim.log"
    conda: "../envs/downloads.yaml"
    params:
        quality=20,
        title="{sra}"
    shell:
        "fastp -q {params.quality}"
        " -i {input} -o {output.fq}"
        " -h {output.html} -j {output.json}"
        " -w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"
