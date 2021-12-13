#! /usr/bin/env python
# -*- coding utf-8 -*-

rule stellarscope:
    conda:
        "../envs/telescope.yaml"
    output:
        "results/telescope/{s}/{s}-TE_counts.tsv",
        "results/telescope/{s}/{s}-updated.bam",
        "results/telescope/{s}/{s}-other.bam"
    input:
        bam = "results/starsolo_algn/{s}/{s}_GDC38.collated.out.bam",
        annotation = rules.telescope_annotation.output,
        barcode = "results/starsolo_algn/{s}/{s}_GDC38.Solo.out/Gene/filtered/barcodes.tsv"
    benchmark: "benchmarks/telescope/{s}_telescope.tsv"
    log:
        "results/telescope/{s}/telescope.log"
    threads: config['telescope_threads']
    params:
        tmpdir = config['local_tmp'],
        out = "results/telescope/{s}",
        exp_tag = "{s}"
    shell:
        """

        getrss() {{
            local cmd=$1
            local param=$2
            ps -C $cmd -o rss= -o args= | grep "$param" | awk '{{$1=int(100 * $1/1024/1024)/100"GB";}}{{ print $1;}}' | while read v; do echo "Memory usage (RSS) for $cmd (param: $param): $v"; done
        }}

        while true; do getrss stellarscope {wildcards.s}; sleep 5; done &

        stellarscope assign\
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
        rules.stellarscope.output
    output:
        touch("results/completed/{s}_completed.txt")
