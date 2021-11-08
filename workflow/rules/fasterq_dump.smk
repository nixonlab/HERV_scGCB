#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download FASTQs
"""

rule fasterq_dump:
    conda:
        "../envs/utils.yaml"
    output:
        "runs/{s}/{s}_1.fastq",
        "runs/{s}/{s}_2.fastq"
    params:
        tmpdir = config['tmpdir'],
        outdir = "runs/{s}"
    threads: min(20, snakemake.utils.available_cpu_count())
    log: "runs/{s}/fasterq_sra_to_fastq.log"
    shell:
        """
        fasterq-dump -e {threads} --temp {params.tmpdir} -O {params.outdir} {wildcards.s} &> {log[0]}
        """

rule conversion_complete:
    input:
        expand("runs/{s}/{s}_1.fastq", s=RUNS)
