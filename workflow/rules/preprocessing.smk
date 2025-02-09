rule trim3:
    input:
        "resources/fastq/{experiment}/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/{experiment}/{accession}_trim3.fastq"
    threads:
        32
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/{experiment}/{accession}_trim3.log"
    message:
        "3'-end adapter trimming for {wildcards.experiment} {wildcards.accession} via Trimmomatic."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "ILLUMINACLIP:resources/adapters.fa:3:20:10 SLIDINGWINDOW:4:20 MINLEN:18"


rule read_collapse:
    input:
        icSHAPE_dir = rules.download_icSHAPE_git.output,
        fq = rules.trim3.output
    output:
        "resources/preprocessed_fastq/{experiment}/{accession}_trim3_nodup.fastq"
    conda:
        "../envs/git.yaml"
    message:
        "Remove PCR duplicates from {wildcards.experiment} {wildcards.accession}."
    shell:
        "perl {input.icSHAPE_dir}/scripts/readCollapse.pl -U {input.fq} -o {output}"


rule trim5:
    input:
        rules.read_collapse.output
    output:
        "resources/preprocessed_fastq/{experiment}/{accession}.fastq"
    params:
        "resources/preprocessed_fastq/{experiment}/{accession}_trim3"
    threads:
        32
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/{experiment}/{accession}_trim5.log"
    message:
        "5'-end adapter trimming for {wildcards.experiment} {wildcards.accession} via Trimmomatic."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "HEADCROP:17 MINLEN:20 && "
        "rm {params}*"


###############################################################
# FASTQC
###############################################################

rule fastqc_raw:
    input:
        "resources/fastq/{experiment}/{accession}.fastq.gz"
    output:
        "results/fastqc/raw/{experiment}/{accession}_fastqc.html"
    params:
        "results/fastqc/raw/{experiment}"
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/fastqc/raw/{experiment}/{accession}.log"
    message:
        "Quality control check on the raw sequencing data of {wildcards.experiment} {wildcards.accession}."
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "{input} "
        "&> {log}"


rule fastqc_preprocessed:
    input:
        rules.trim5.output
    output:
        "results/fastqc/preprocessed/{experiment}/{accession}_fastqc.html"
    params:
        "results/fastqc/preprocessed/{experiment}"
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/fastqc/preprocessed/{experiment}/{accession}.log"
    message:
        "Quality control check on the preprocessed sequencing data of {wildcards.experiment} {wildcards.accession}."
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "{input} "
        "&> {log}"