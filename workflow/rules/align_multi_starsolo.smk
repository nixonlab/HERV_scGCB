#! /usr/bin/env python
# -*- coding: utf-8 -*-

def get_input_files(wildcards):
    _sruns = samples.loc[wildcards.samp]['runid']
    _layouts = runs.loc[_sruns]['layout']

    if not all(_ == 'PAIRED' for _ in _layouts):
        msg = f'All runs for sample "{wildcards.samp}" should be PAIRED for this rule'
        raise WorkflowError(msg)

    return expand("results/{{dataset}}/fastq/{{samp}}/{sra_run}_{r}.fastq.gz", sra_run=_sruns, r=[1,2])

# this function creates the actual strings with the fastq files in a way STAR takes them
#
def star_readfiles_in(wildcards, input, output):
    _sruns = samples.loc[wildcards.samp]['runid']

    r1_files = [f"results/{wildcards.dataset}/fastq/{wildcards.samp}/{sra_run}_1.fastq.gz" for sra_run in _sruns]
    r2_files = [f"results/{wildcards.dataset}/fastq/{wildcards.samp}/{sra_run}_2.fastq.gz" for sra_run in _sruns]

    r1 = ','.join(r1_files)
    r2 = ','.join(r2_files)

    bcmate = config['barcode_mate'][meta_table.loc[wildcards.samp].chemistry]
    if bcmate == 'R1':
        readBC = r1
        readTX = r2
    elif bcmate == 'R2':
        readBC = r2
        readTX = r1
    return f'{readTX} {readBC}'

rule align_multi_starsolo:
    output:
        "results/{dataset}/align_multi_starsolo/{samp}/Aligned.sortedByCoord.out.bam",
        "results/{dataset}/align_multi_starsolo/{samp}/Solo.out/Gene/filtered/matrix.mtx",
        "results/{dataset}/align_multi_starsolo/{samp}/Solo.out/Gene/filtered/barcodes.tsv",
        "results/{dataset}/align_multi_starsolo/{samp}/Solo.out/Gene/filtered/features.tsv"
    input:
        get_input_files
    params:
        readfiles_in = star_readfiles_in,
        common_args = lambda wc: ' '.join(f'--{k} {v}' for k,v in config['align_multi_starsolo'][wc.dataset].items()),
        whitelist = lambda wc: config['whitelist'][wc.dataset]
    threads: 18
    conda: "../envs/star.yaml"
    shell:
        '''
        mkdir -p $(dirname {output[0]})
        STAR\
            --runThreadN {threads}\
            --genomeDir {config[star_index]}\
            --readFilesIn {params.readfiles_in}\
            --soloCBwhitelist {params.whitelist}\
            {params.common_args}\
            --outFileNamePrefix $(dirname {output[0]})/
        '''

