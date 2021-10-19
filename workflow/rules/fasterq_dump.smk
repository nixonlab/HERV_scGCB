#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download FASTQs
"""

rule fasterq_dump:
    conda:
        "../envs/utils.yaml"
    output:
        "runs/{run_acc}/{run_acc}_1.fastq",
        "runs/{run_acc}/{run_acc}_2.fastq"
    params:
        tmpdir = config['tmpdir'],
        outdir = "runs/{run_acc}"
    threads: min(20, snakemake.utils.available_cpu_count())
    log: "runs/{run_acc}/fasterq_sra_to_fastq.log"
    shell:
        """
        fasterq-dump -e {threads} --temp {params.tmpdir} -O {params.outdir} {wildcards.run_acc} &> {log[0]}
        """

rule conversion_complete:
    input:
        expand("runs/{run_acc}/{run_acc}_1.fastq", run_acc=RUNS)
