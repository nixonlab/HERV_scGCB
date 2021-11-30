######################################################
## Analyze droplet single-cell RNA sequencing data
## with STARsolo
######################################################
## to do:
## input as function for cDNA and barcodes
## a scheme to name output directories
## see if I can integrate PEP here

rule starsolo_alignment:
    """
    Align sequencing reads from a 10x V3 single-cell RNA-seq experiment using STARsolo
    """
    conda:
        "../envs/star.yaml"
    input:
        cDNA = "runs/{s}/{s}_2.fastq.gz",
        barcodes = "runs/{s}/{s}_1.fastq.gz",
        genome = config['indexes']['star'],
        whitelist = "refs/downloads/whitelist.10x.v2.txt"
    output:
        "results/starsolo_algn/{s}/{s}_GDC38.Aligned.sortedByCoord.out.bam"
    params:
        out_prefix="results/starsolo_algn/{s}/{s}_GDC38.",
        cb_start=config["cellbarcode_start"],
        cb_length=config["cellbarcode_length"],
        umi_start=config["umi_start"],
        umi_length=config["umi_length"],
        max_multimap=config["max_multimap"]
    benchmark: "benchmarks/star_alignment/{s}_star_alignment.tsv"
    threads: config['star_alignment_threads']
    shell:
        '''
        #--- STARsolo (turned on by --soloType CB_UMI_Simple)
        STAR\
            --runThreadN {threads}\
            --genomeDir {input.genome}\
            --readFilesIn {input.cDNA} {input.barcodes}\
            --readFilesCommand gunzip -c\
            --soloType CB_UMI_Simple\
            --soloCBwhitelist {input.whitelist}\
            --outFilterMultimapNmax {params.max_multimap}\
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM\
            --outSAMtype BAM SortedByCoordinate\
            --outFileNamePrefix {params.out_prefix}
        '''

rule samtools_collate:
    conda: "../envs/utils.yaml"
    output: "results/starsolo_algn/{s}/{s}_GDC38.collated.out.bam"
    input: "results/starsolo_algn/{s}/{s}_GDC38.Aligned.sortedByCoord.out.bam"
    benchmark: "benchmarks/samtools_collate/{s}_samtools_collate.tsv"
    params:
        tmpdir = config['fasterq_dump_tmp']
    threads: config['samtools_collate_threads']
    shell:
        '''
        samtools collate\
        {input}\
        -o {output}\
        -@ {threads}\
        {params.tmpdir}
        '''
