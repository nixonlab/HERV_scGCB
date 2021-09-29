#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert SRA files to fastq
"""

localrules: sra_to_fastq

rule sra_to_fastq:
    conda:
        "../envs/utils.yaml"
    output:
        "runs/{run_acc}/{run_acc}_1.fastq",
        "runs/{run_acc}/{run_acc}_2.fastq"
    input: config['indir']
    params:
        tmpdir = config['tmpdir'],
        outdir = "runs/{run_acc}"
    threads: min(20, snakemake.utils.available_cpu_count())
    shell:
        """
        fasterq-dump -e {threads} --temp {params.tmpdir} --outdir {params.rundir} --split-3 {input} &> {log[0]}
        """
    log:
        "runs/{run_acc}/fasterq_sra_to_fastq.log"

localrules: conversion_complete
rule conversion_complete:
    input:
        expand("runs/{run_acc}/{run_acc}_1.fastq", run_acc=RUNS)
