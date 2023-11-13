#! /usr/bin/env python
# -*- coding: utf-8 -*-

# take the runid associated with the sample and check that all the runs have the same layout
#
def get_fq_r1(wildcards):
    _sruns = samples.loc[wildcards.samp]['runid']
    _layouts = runs.loc[_sruns]['layout']
    if not all(_layouts.iloc[0] == _ for _ in _layouts):
        msg = f'All runs for sample "{wildcards.samp}" do not have the same layout (PAIRED or SINGLE)'
        raise WorkflowError(msg)
    if _layouts.iloc[0] == 'PAIRED':
        return expand("results/{{dataset}}/fastq/{{samp}}/{sra_run}_1.fastq.gz", sra_run=_sruns)
    else:
        raise WorkflowError(f'Invalid layout: "{_layouts[0]}"')

def get_fq_r2(wildcards):
    _sruns = samples.loc[wildcards.samp]['runid']
    return expand("results/{{dataset}}/fastq/{{samp}}/{sra_run}_2.fastq.gz", sra_run=_sruns)

# this function creates the actual strings with the fastq files in a way STAR takes them
#
def star_readfiles_in(wildcards, input, output):
    bcmate = config['barcode_mate'][meta_table.loc[wildcards.samp].chemistry]
    if bcmate == 'R1':
        readBC = ','.join(input['R1'])
        readTX = ','.join(input['R2'])
    elif bcmate == 'R2':
        readBC = ','.join(input['R2'])
        readTX = ','.join(input['R1'])
    return f'{readTX}  {readBC}'

rule align_multi_starsolo:
    output:
        "results/{dataset}/align_multi_starsolo/{samp}/Aligned.sortedByCoord.out.bam",
        "results/{dataset}/align_multi_starsolo/{samp}/Solo.out/Gene/filtered/matrix.mtx",
        "results/{dataset}/align_multi_starsolo/{samp}/Solo.out/Gene/filtered/barcodes.tsv",
        "results/{dataset}/align_multi_starsolo/{samp}/Solo.out/Gene/filtered/features.tsv"
    input:
        R1 = get_fq_r1, 
        R2 = get_fq_r2
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

