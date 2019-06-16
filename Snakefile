GENOME_VERSION = "GRCh38"
ENSEMBL_VERSION = "GRCh38.96"
GENOME_FA = f"ensembl/Homo_sapiens.{GENOME_VERSION}.dna.primary_assembly.fa"
ENSEMBL_GFF = f"ensembl/Homo_sapiens.{ENSEMBL_VERSION}.gff3"
REF_PREFIX = f"ensembl/Homo_sapiens.{ENSEMBL_VERSION}/RsemStarReference."

configfile: "config.yaml"

rule all:
    input:
        expand("output/{fq}.genes.results", fq=config["fq"])

rule download_ensembl_references:
    output:
        gfa=GENOME_FA,
        gff=ENSEMBL_GFF,
    log: "ensembl/downloads.log"
    shell:
        "(wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens.{GENOME_VERSION}.dna.primary_assembly.fa.gz | "
        "gunzip -c > {output.gfa} && "
        "wget -O - ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.{ENSEMBL_VERSION}.gff3.gz | "
        "gunzip -c > {output.gff}) 2> {log}"

rule rsem_star_genome:
    input:
        efa="data/ERCC.fa",
        gfa=GENOME_FA,
        gff=ENSEMBL_GFF
    output:
        suffix = REF_PREFIX + "SA"
    threads: 12
    resources: mem_mb=60000
    shell:
        "rsem-prepare-reference --num-threads {threads} --star"
        " --gff3 {input.gff} {input.efa} {input.gfa} " + REF_PREFIX

rule rsem_star_align:
    input:
        suffix=REF_PREFIX + "SA",
        fq=expand("{fq}", fq=config["fq"]) #single end only
    output:
        "output/{fq}.isoforms.results",
        "output/{fq}.genes.results",
        "output/{fq}.transcript.bam",
        "output/{fq}.time",
        "output/{fq}.stat",
    log:
    resources: mem_mb=60000
    threads: 12
    shell:
        "rsem-calculate-expression --no-bam-output --time --calc-ci --star"
        " --num-threads {threads} {input.fq} " + REF_PREFIX
