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
        tmpdir = config['fasterq_dump_tmp'],
        outdir = "runs/{s}"
    threads: 8
    resources:
        mem_mb = 10000, disk_mb = 60000
    log: "runs/{s}/fasterq_sra_to_fastq.log"
    shell:
        """
        fasterq-dump -e {threads} --temp {params.tmpdir} -O {params.outdir} {wildcards.s} &> {log[0]}
        """

rule gzip:
    output:
        R1_in = "runs/{s}/{s}_1.fastq.gz",
        R2_in = "runs/{s}/{s}_2.fastq.gz"
    input:
        R1_in = "runs/{s}/{s}_1.fastq",
        R2_in = "runs/{s}/{s}_2.fastq"
    shell:
        """
        gzip {input.R1_in}
        gzip {input.R2_in}
        """

rule conversion_complete:
    input:
        expand("runs/{s}/{s}_1.fastq.gz", s=RUNS)
