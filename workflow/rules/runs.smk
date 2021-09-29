#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download data from SRA
"""
localrules: sra_download

rule sra_download:
    output:
        "runs/{run_acc}/unmapped.bam"
    log:
        "runs/{run_acc}/fasterq-dump.log"
    params:
        keyfile = config['dbgap_key'] if 'dbgap_key' in config else 'nonexistent',
        tmpdir = config['tmpdir'],
        rundir = "runs/{run_acc}"
    conda:
        "envs/utils.yaml"
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
              {run_acc}\
            &> {log[0]}
        else
            CDIR=$PWD
            NCBIDIR=$(vdb-config --import {params.keyfile} | tail -n+2 | cut -d"'" -f2)
            cd $NCBIDIR
            fasterq-dump\
              -e {threads}\
              -t {params.tmpdir}\
              {run_acc}\
            &> $CDIR/{log[0]}
            cd $CDIR
            mv $NCBIDIR/{run_acc}* {params.rundir}/
        fi
        '''

localrules: sample_complete
rule sample_complete:
    input:
        rules.sra_download.output,
    output:
        lambda wc: touch("runs/{run_acc}/completed.txt", run_acc=SAMPLE_RUN[wc.samp_acc])
