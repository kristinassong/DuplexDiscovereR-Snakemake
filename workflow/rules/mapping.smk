###############################################################
# CoCo -- Corrected annotation
###############################################################

rule coco_ca:
    input:
        gtf = config["genome_gtf"],
        coco_dir = rules.download_CoCo_git.output
    output:
        gtf_corrected = "resources/correct_annotation.gtf"
    conda:
        "../envs/coco.yaml"
    message:
        "Generate corrected annotation using CoCo."
    shell:
        "export PATH=$PWD/{input.coco_dir}/bin:$PATH && "
        "python3 {input.coco_dir}/bin/coco.py ca {input.gtf} -o {output.gtf_corrected}"


###############################################################
# Map reads to the genome
###############################################################

rule STAR_index:
    input:
        fasta = config["genome_fasta"],
        gtf = rules.coco_ca.output.gtf_corrected
    output:
        directory("results/STAR/index")
    conda:
        "../envs/mapping.yaml"
    threads:
        16
    log:
        "results/logs/STAR/index.log"
    message:
        "Generate genome index files using STAR."
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "--genomeSAindexNbases 14 "
        "--genomeChrBinNbits 18 "
        "--sjdbGTFfile {input.gtf} "
        "&> {log}"


rule STAR_align:
    input:
        fq = rules.trim5.output,
        idx = rules.STAR_index.output
    output:
        "results/STAR/{experiment}/{accession}/{accession}_Aligned.sortedByCoord.out.bam"
    params:
        out_prefix = "results/STAR/{experiment}/{accession}/{accession}_"
    conda:
        "../envs/mapping.yaml"
    threads:
        16
    log:
        "results/logs/STAR/{experiment}/{accession}_align.log"
    message:
        "Align {wildcards.experiment} {wildcards.accession} reads to the reference genome using STAR."
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {input.idx} "
        "--readFilesIn {input.fq} "
        "--outFileNamePrefix {params.out_prefix} "
        "--runThreadN {threads} --genomeLoad NoSharedMemory --outReadsUnmapped Fastx  "
        "--outFilterMultimapNmax 10 --outFilterScoreMinOverLread 0 "
        "--outSAMattributes All "
        "--outSAMtype BAM Unsorted SortedByCoordinate --alignIntronMin 1 --scoreGap 0 "
        "--scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 "
        "--scoreGenomicLengthLog2scale -1 --chimOutType Junctions WithinBAM HardClip "
        "--chimSegmentMin 5 --chimJunctionOverhangMin 5 --chimScoreJunctionNonGTAG 0 "
        "-- chimScoreDropMax 80 --chimNonchimScoreDropMin 20 --chimMultimapNmax 5 "
        "--limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 "
        "--quantMode GeneCounts "
        "&> {log}"


###############################################################
# DuplexDiscovereR
###############################################################

rule dna_to_rna_fasta:
    input:
        dna = config["genome_fasta"]
    output:
        rna = "resources/genome_rna.fa"
    conda:
        "../envs/r.yaml"
    message:
        "Convert DNA fasta to RNA fasta for DuplexDiscoverer."
    script:
        "../scripts/dna_to_rna_fasta.R"


rule DuplexDiscoverer:
    input:
        star_bam = rules.STAR_align.output,
        gtf = rules.coco_ca.output.gtf_corrected,
        fa = rules.dna_to_rna_fasta.output.rna
    output:
        dg = "results/duplexdiscoverer/{experiment}/{accession}_dg.tsv"
    params:
        sample = "{experiment}_{accession}",
        table = "STAR", # STAR or bedpe
        lib = "SE", # SE or PE
        sj = 'results/STAR/{experiment}/{accession}/{accession}_SJ.out.tab',
        chim = 'results/STAR/{experiment}/{accession}/{accession}_Chimeric.out.junction',
        gc = 'results/STAR/{experiment}/{accession}/{accession}_ReadsPerGene.out.tab',
    conda:
        "../envs/r.yaml"
    message:
        "Run DuplexDiscovereR to detect RNA-RNA interactions."
    script:
        "../scripts/DuplexDiscoverer.R"


rule merge_dgs:
    input:
        paris = expand(rules.DuplexDiscoverer.output.dg,accession=config['PARIS'],experiment=['PARIS']),
        paris2 = expand(rules.DuplexDiscoverer.output.dg,accession=config['PARIS2'],experiment=['PARIS2']),
        splash = expand(rules.DuplexDiscoverer.output.dg,accession=config['SPLASH'],experiment=['SPLASH']),
        ligr_seq = expand(rules.DuplexDiscoverer.output.dg,accession=config['LIGR_seq'],experiment=['LIGR_seq'])
    params:
        samples = "resources/samples.txt",
        outdir = "results/duplexdiscoverer",
        gtf = rules.coco_ca.output.gtf_corrected
    output:
        dg = "results/duplexdiscoverer/all_dgs.tsv"
    conda:
        "../envs/coco.yaml"
    message:
        "Combine all DGs into one file."
    script:
        "../scripts/merge_dgs.py"