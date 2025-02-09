rule wget_fq:
    input:
        "resources/sra_ids.txt"
    output:
        "resources/fastq/download_status.txt"
    message:
        "Download FASTQ files from GEO via wget."
    shell:
        "echo \'Downloading fastq files from GEO...\' > {output} && "
        "while IFS=\'\t\' read exp url; do "
        "date >> {output}; "
        "echo $exp $url >> {output}; "
        "wget -P resources/fastq/$exp $url; "
        "done < {input} && "
        "echo \'Download complete\' >> {output}"


rule download_CoCo_git:
    output:
        directory('resources/CoCo')
    conda:
        '../envs/git.yaml'
    message:
        "Download CoCo git repository."
    shell:
        'git clone git@github.com:scottgroup/coco.git {output}'


rule download_icSHAPE_git:
    output:
        directory('resources/icSHAPE')
    conda:
        '../envs/git.yaml'
    message:
        "Download icSHAPE git repository."
    shell:
        "git clone git@github.com:qczhang/icSHAPE.git {output}"