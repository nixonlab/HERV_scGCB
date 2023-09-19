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
        temp("results/sra/{sra_run}/{sra_run}.sra")
    log:
        "results/sra/{sra_run}/prefetch_sra_{sra_run}.log"
    shell:
        '''
	prefetch\
	    {wildcards.sra_run}\
	    --output-file {output[0]}\
	    &> {log}
	'''

# function to figure out how many fastq output files there are
#
# f-strings are used to embed expressions inside string literals
# The curly braces {} inside the string indicate where the value of wildcards.sra_run should be inserted into the string
#def get_fastq_output(wildcards):
#    layout = runs.loc[wildcards.sra_run]['layout']
#    if layout == "SINGLE":
#        return f"results/fastq/{wildcards.sra_run}/{wildcards.sra_run}_1.fastq.gz"
#    elif layout == "PAIRED":
#        return (
#            f"results/fastq/{wildcards.sra_run}/{wildcards.sra_run}_1.fastq.gz",
#            f"results/fastq/{wildcards.sra_run}/{wildcards.sra_run}_2.fastq.gz"
#        )

#fastq-dump adds SRA ID to each read in the file, to avoid it, use –-origfmt
"""
dump SRA to FASTQ
"""
rule sra_to_fastq:
    conda:
        "../envs/sra_data.yaml"
    output:
        temp("results/fastq/{sra_run}/{sra_run}_1.fastq.gz"),
        temp("results/fastq/{sra_run}/{sra_run}_2.fastq.gz")
    input:
        rules.prefetch_sra.output
    params:
        layout = lambda wc: runs.loc[wc.sra_run]['layout']
    threads: 
        config['sra_to_fastq_threads']
    log:
        "results/fastq/{sra_run}/sra_to_fastq_{sra_run}.log"
    shell:
        '''
        tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.sra_run}.XXXXXX)
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


