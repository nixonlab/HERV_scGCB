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
        cDNA = "samples/{samid}_2.fastq.gz",
        barcodes = "samples/{samid}_1.fastq.gz",
        genome = config['indexes']['star'],
        whitelist = lambda wc: SAMPLE_WHITELIST[wc.samid][0]
    output:
        "results/starsolo_algn/{samid}/{samid}_GDC38.Aligned.sortedByCoord.out.bam",
        "results/starsolo_algn/{samid}/{samid}_GDC38.Solo.out/Gene/filtered/barcodes.tsv"
    params:
        out_prefix="results/starsolo_algn/{samid}/{samid}_GDC38.",
        cb_start=config["cellbarcode_start"],
        cb_length=config["cellbarcode_length"],
        umi_start=config["umi_start"],
        umi_length=config["umi_length"],
        max_multimap=config["max_multimap"]
    benchmark: "benchmarks/star_alignment/{samid}_star_alignment.tsv"
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
            --soloBarcodeReadLength 0\
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM\
            --outSAMtype BAM SortedByCoordinate\
            --outFileNamePrefix {params.out_prefix}
        '''

rule stellarscope_cellsort:
    conda:
        "../envs/telescope.yaml"
    output:
        sorted_bam = "results/stellarscope_cellsort/{samid}/{samid}.Aligned.sortedByCB.bam"
    input:
        alignment_bam=rules.starsolo_alignment.output[0],
        filtered_barcodes=rules.starsolo_alignment.output[1]
    benchmark:
        "benchmarks/stellarscope_cellsort/{samid}_stellarscope_cellsort.tsv"
    log:
        "results/stellarscope_cellsort/{samid}/stellarscope_cellsort.log"
    threads: 4
    shell:
        '''
stellarscope cellsort --ncpu {threads} --outfile {output.sorted_bam} {input.alignment_bam} {input.filtered_barcodes}
        '''
