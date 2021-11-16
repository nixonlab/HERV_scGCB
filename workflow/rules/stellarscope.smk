#! /usr/bin/env python
# -*- coding utf-8 -*-

rule stellarscope:
    conda: "../envs/telescope.yaml"
    output:
        "results/telescope/{s}/{s}_telescope.report.tsv",
        "results/telescope/{s}/{s}_telescope.updated.bam",
        "results/telescope/{s}/{s}_telescope.other.bam"
    input:
        bam = "results/starsolo_algn/{s}/{s}_GDC38.collated.out.bam",
        annotation = rules.telescope_annotation.output
    benchmark: "benchmarks/telescope/{s}_telescope.tsv"
    log:
        "results/telescope/{s}/telescope.log"
    threads: config['telescope_threads']
    params:
        tmpdir = config['local_tmp']
    shell:
        """
        telescope sc assign\
         {input[0]}\
         {input[1]}\
         2>&1 | tee {log[0]}
        """

rule sample_complete:
    input:
        rules.stellarscope.output
    output:
        touch("results/completed/{s}_completed.txt")
