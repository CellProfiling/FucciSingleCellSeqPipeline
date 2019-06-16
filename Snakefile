GENOME_VERSION = "GRCh38"
ENSEMBL_VERSION = "96"
GENEMODEL_VERSION = GENOME_VERSION + "." + ENSEMBL_VERSION
GENOME_FA = f"ensembl/Homo_sapiens.{GENOME_VERSION}.dna.primary_assembly.fa"
ENSEMBL_GFF = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}.gff3"
REF_PREFIX = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}/RsemStarReference"
REF_FOLDER = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}/"

FQ_FOLDER = "/mnt/e/ProjectsActive/tonsil/data_spritzsnake/trimmed/"
# (FQ_PREFIXES,) = glob_wildcards(FQ_FOLDER + "{fq}.fastq.gz")
(FQ_PREFIXES,) = glob_wildcards(FQ_FOLDER + "{fq}.fastq.gz")

rule all:
    input: "output/counts.tsv"

# rule clean:
#     shell: "rm -rf ensembl output"

rule download_ensembl_genome:
    output:
        gfa=GENOME_FA,
    log: "ensembl/downloads_fa.log"
    shell:
        "(wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.fa.gz | "
        "gunzip -c > {output.gfa}) 2> {log}"

rule download_ensembl_genemodel:
    output:
        gff=ENSEMBL_GFF,
    log: "ensembl/download_gff.log"
    shell:
        "(wget -O - ftp://ftp.ensembl.org/pub/release-" + ENSEMBL_VERSION + "/gff3/homo_sapiens/Homo_sapiens." + GENEMODEL_VERSION + ".gff3.gz | "
        "gunzip -c > {output.gff}) 2> {log}"

rule fix_gff3_for_rsem:
    input: ENSEMBL_GFF
    output: ENSEMBL_GFF + ".fix.gff3"
    shell: "python scripts/fix_gff3_for_rsem.py {input} {output}"

rule rsem_star_genome:
    input:
        efa="data/ERCC.fa",
        gfa=GENOME_FA,
        gff=ENSEMBL_GFF + ".fix.gff3"
    output:
        REF_PREFIX + ".gtf",
        suffix = REF_FOLDER + "SA"
    threads: 99
    resources: mem_mb=60000
    log: "ensembl/prepare-reference.log"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star"
        " --gff3 {input.gff} \"{input.efa}\",\"{input.gfa}\" " + REF_PREFIX +
        ") 2> {log}"

rule rsem_star_align:
    input:
        suffix=REF_FOLDER + "SA",
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
        " --num-threads {threads} --star-gzipped-read-file " + FQ_FOLDER + "{wildcards.fq}.fastq.gz " + REF_PREFIX + " output/{wildcards.fq}) &> {log}"

rule make_rsem_dataframe:
    input: expand("output/{fq}.genes.results", fq=FQ_PREFIXES)
    output:
        counts="output/counts.tsv",
        tpm="output/tpm.tsv"
    shell:
        "echo \"python scripts/make_rsem_dataframe.py 4 {output.counts} {input}\" && "
        "python scripts/make_rsem_dataframe.py 4 {output.counts} {input} && "
        "python scripts/make_rsem_dataframe.py 5 {output.tpm} {input}"

# 190616: This isn't working. Nothing is getting aligned even in cases where I know there
# should be good alignments. Should I align with STAR and then import the alignments for quant
# within RSEM?
