#! /usr/bin/env python
# -*- coding: utf-8 -*-

wildcard_constraints:
    sra_run = "SRR[0-9]{6,9}"

"""
Download SRA file 
"""
rule prefetch_sra:
    conda: 
        "../envs/sra_data.yaml"
    output:
        temp("results/{dataset}/sra/{samp}/{sra_run}.sra")
    log:
        "results/{dataset}/sra/{samp}/prefetch_{sra_run}.log"
    shell:
        '''
	prefetch\
            --max-size 50G\
	    {wildcards.sra_run}\
	    --output-file {output[0]}\
	    &> {log}
	'''

#fastq-dump adds SRA ID to each read in the file, to avoid it, use –-origfmt
"""
dump SRA to FASTQ
"""
rule sra_to_fastq:
    conda:
        "../envs/sra_data.yaml"
    output:
        "results/{dataset}/fastq/{samp}/{sra_run}_1.fastq.gz",
        "results/{dataset}/fastq/{samp}/{sra_run}_2.fastq.gz"
    input:
        rules.prefetch_sra.output
    params:
        layout = lambda wc: runs.loc[wc.sra_run]['layout']
    threads: 
        config['sra_to_fastq_threads']
    log:
        "results/{dataset}/fastq/{samp}/sra_to_fastq_{sra_run}.log"
    shell:
        '''
        tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sra_run}.XXXXXX)
	parallel-fastq-dump\
            --threads {threads}\
            --tmpdir $tdir\
	    --outdir $tdir\
	    --split-files\
            --origfmt\
            --skip-technical\
            --sra-id {input}\
	    &> {log}

        mkdir -p $(dirname {output[0]})
        pigz -p {threads} -c $tdir/{wildcards.sra_run}_1.fastq > {output[0]}
        if [[ "{params.layout}" == "PAIRED" ]]; then
            pigz -p {threads} -c $tdir/{wildcards.sra_run}_2.fastq > {output[1]}
        else
            touch {output[1]}
        fi

        rm -rf $tdir

        #chmod 0440 {output}
        #chmod 0550 $(dirname {output[0]})
	'''

