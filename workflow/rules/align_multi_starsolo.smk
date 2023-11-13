#! /usr/bin/env python
# -*- coding: utf-8 -*-

# take the runid associated with the sample and check that all the runs have the same layout
#
def fastq_input(wc):
    _sruns = samples.loc[wc.samp]['runid']
    _layouts = runs.loc[_sruns]['layout']
    if not all(_layouts.iloc[0] == _ for _ in _layouts):
        msg = f'All runs for sample "{wc.samp}" do not have the same layout (PAIRED or SINGLE)'
        raise WorkflowError(msg)
    if _layouts.iloc[0] == 'SINGLE':
        return expand("results/{dataset}/fastq/{samp}/{sra_run}_1.fastq.gz", sra_run=_sruns)
    elif _layouts.iloc[0] == 'PAIRED':
        return expand("results/{dataset}/fastq/{samp}/{sra_run}_{rn}.fastq.gz", sra_run=_sruns, rn=[1,2])
    else:
        raise WorkflowError(f'Invalid layout: "{_layouts[0]}"')

def get_fq_r1(wc):
    _sruns = samples.loc[wc.samp]['runid']
    ret1 = ','.join(expand("results/fastq/{runid}/{runid}_1.fastq.gz", runid=_sruns))
    
    return expand("results/{dataset}/fastq/{samp}/{sra_run}_1.fastq.gz", sra_run=_sruns)

def get_fq_r2(wc):
    _sruns = samples.loc[wc.samp]['runid']
    return expand("results/{dataset}/fastq/{samp}/{sra_run}_2.fastq.gz", sra_run=_sruns)

# this function creates the actual strings with the fastq files in a way STAR takes them
#
def star_readfiles_in(wildcards, input, output):
    _sruns = samples.loc[wc.samp]['runid']
    r1 = ','.join(expand("results/{dataset}/fastq/{samp}/{sra_run}_1.fastq.gz", sra_run=_sruns))
    r2 = ','.join(expand("results/{dataset}/fastq/{samp}/{sra_run}_2.fastq.gz", sra_run=_sruns))

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
        fastq_input
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

