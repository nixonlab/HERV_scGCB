#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download data from SRA
"""
localrules: sra_download

print('here')

rule sra_download:
    output:
        "runs/{run_acc}/{run_acc}_1.fastq",
        "runs/{run_acc}/{run_acc}_2.fastq"
    log:
        "runs/{run_acc}/fasterq-dump.log"
    params:
        keyfile = config['dbgap_key'] if 'dbgap_key' in config else 'nonexistent',
        tmpdir = config['tmpdir'],
        rundir = "runs/{run_acc}"
    conda:
        "../envs/utils.yaml"
    threads: min(20, snakemake.utils.available_cpu_count())
    shell:
        '''
        mkdir -p {params.rundir}

        #--- Download using fasterq-dump
        if [[ ! -e "{params.keyfile}" ]]; then
            fasterq-dump\
              -e {threads}\
              -t {params.tmpdir}\
              -O {params.rundir}\
              {wildcards.run_acc}\
            &> {log[0]}
        else
            CDIR=$PWD
            NCBIDIR=$(vdb-config --import {params.keyfile} | tail -n+2 | cut -d"'" -f2)
            cd $NCBIDIR
            fasterq-dump\
              -e {threads}\
              -t {params.tmpdir}\
              {wildcards.run_acc}\
            &> $CDIR/{log[0]}
            cd $CDIR
            mv $NCBIDIR/{wildcards.run_acc}* {params.rundir}/
        fi
        '''

localrules: sample_complete
rule sample_complete:
    input:
        expand("runs/{run_acc}/{run_acc}_1.fastq", run_acc=RUNS)
