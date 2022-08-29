#! /usr/bin/env python
# -*- coding utf-8 -*-

rule stellarscope:
    conda:
        "../envs/telescope.yaml"
    output:
        "results/telescope/{samid}/{samid}-TE_counts.mtx"
    input:
        bam = "results/starsolo_algn/{samid}/{samid}_GDC38.collated.out.bam",
        annotation = rules.telescope_annotation.output,
        barcode = "results/starsolo_algn/{samid}/{samid}_GDC38.Solo.out/Gene/filtered/barcodes.tsv"
    benchmark: "benchmarks/telescope/{samid}_telescope.tsv"
    log:
        "results/telescope/{samid}/telescope.log"
    threads: config['telescope_threads']
    params:
        tmpdir = config['local_tmp'],
        out = "results/telescope/{samid}",
        exp_tag = "{samid}"
    shell:
        """

        getrss() {{
            local cmd=$1
            local param=$2
            ps -C $cmd -o rss= -o args= | grep "$param" | awk '{{$1=int(100 * $1/1024/1024)/100"GB";}}{{ print $1;}}' | while read v; do echo "Memory usage (RSS) for $cmd (param: $param): $v"; done
        }}

        while true; do getrss stellarscope {wildcards.samid}; sleep 5; done &

        stellarscope assign\
        --updated_sam\
        --exp_tag {params.exp_tag}\
        --outdir {params.out}\
         {input[0]}\
         {input[1]}\
         --barcodefile {input[2]}\
         2>&1 | tee {log[0]}

        """

rule cell_sort:
    conda:
        "../envs/samtools.yaml"
    input:
        bam = "results/starsolo_algn/{samid}/{samid}_GDC38.collated.out.bam",
        barcode = "results/starsolo_algn/{samid}/{samid}_GDC38.Solo.out/Gene/filtered/barcodes.tsv"
    output: "results/telescope/{samid}/{samid}-Aligned.sortedByCB.out.bam"
    threads: config['telescope_threads']
    params:
        tmpdir = config['cell_sort_tmp'],
        outfile = "results/telescope/{samid}-Aligned.sortedByCB.out.bam"
    log: "results/telescope/{samid}/cellsort.log"
    shell:
        """
        samtools view -u -F 4 -D CB:{input.barcode} {input.bam} | samtools sort -@ {threads} -n -t CB > {output}
        """

rule stellarscope_individual:
    conda: "../envs/telescope.yaml"
    input:
        bam = "results/telescope/{samid}/{samid}-Aligned.sortedByCB.out.bam",
        annotation = rules.telescope_annotation.output,
        barcode = "results/starsolo_algn/{samid}/{samid}_GDC38.Solo.out/Gene/filtered/barcodes.tsv"
    output: "results/telescope_individual/{samid}/{samid}_individual-TE_counts.mtx"
    threads: config['telescope_threads']
    params:
        tmpdir = config['local_tmp'],
        out = "results/telescope_individual/{samid}",
        exp_tag = "{samid}_individual"
    log: "results/telescope_individual/{samid}/telescope.log"
    shell:
        """
        stellarscope assign\
        --pooling_mode individual\
        --updated_sam\
        --exp_tag {params.exp_tag}\
        --outdir {params.out}\
         {input[0]}\
         {input[1]}\
         --barcodefile {input[2]}\
         2>&1 | tee {log[0]}
        """

rule sample_complete:
    input:
        rules.stellarscope.output,
        rules.cell_sort.output,
        rules.stellarscope_individual.output
    output:
        touch("results/completed/{samid}_completed.txt")
