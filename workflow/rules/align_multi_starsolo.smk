#! /usr/bin/env python
# -*- coding: utf-8 -*-

# this function makes sure that the fastq files exist (which is why we request it as input)
#
# take the runid associated with the sampleid
# and check that all the runs for the sample have the same layout
#
# if the layout is single, the function returns a list of file paths to the single end fastqs
# if the layout is paired, it returns a list of file paths to the paired end fastqs
def fastq_input(wildcards):
    _sruns = samples.loc[wildcards.samp]['runid']
    _layouts = runs.loc[_sruns]['layout']
    if not all(_layouts[0] == _ for _ in _layouts):
        msg = f'All runs for sample "{wildcards.samp}" do not have the same layout (PAIRED or SINGLE)'
        raise WorkflowError(msg)
    
    if _layouts[0] == 'SINGLE':
        return expand("results/{dataset}/fastq/{samp}/{sra_run}_1.fastq.gz", runid=_sruns)
    elif _layouts[0] == 'PAIRED':
        return expand("results/{dataset}/fastq/{samp}/{sra_run}_1.fastq.gz", runid=_sruns, rn=[1,2])
    else:
        raise WorkflowError(f'Invalid layout: "{_layouts[0]}"')

# WE SHOULD CHECK HERE THAT SAMPLES ARE assay = rna_seq



# this function creates the actual strings with the fastq files in a way STAR takes them
#
def star_readfiles_in(wildcards, input, output):
    bcmate = config['barcode_mate'][meta_table.loc[wildcards.samp].chem_meta]
        

def readfiles_param(wildcards):
    _sruns = samples.loc[wildcards.sample_id]['runid']
    _layouts = runs.loc[_sruns]['layout']
    if not all(_layouts[0] == _ for _ in _layouts):
        msg = f'All runs for sample "{wildcards.sample_id}" do not have the same layout (PAIRED or SINGLE)'
        raise WorkflowError(msg)
    
    ret1 = ','.join(expand("results/fastq/{runid}/{runid}_1.fastq.gz", runid=_sruns))
    if _layouts[0] == 'SINGLE':
        return ret1
    elif _layouts[0] == 'PAIRED':
        ret2 = ','.join(expand("results/fastq/{runid}/{runid}_2.fastq.gz", runid=_sruns))
        return f'{ret1} {ret2}'
    else:
        raise WorkflowError(f'Invalid layout: "{_layouts[0]}"')

# this function gets the values for all the STAR parameters that are defined in the input
#
def star_common_args(wildcards):
    return ' '.join(f'--{k} {v}' for k,v in config['star_alignment'].items())

rule align_multi_starsolo:
    output:
        "results/RNAseq/star_alignment/{sample_id}/Aligned.sortedByCoord.out.bam",
        "results/RNAseq/star_alignment/{sample_id}/ReadsPerGene.out.tab",
        "results/RNAseq/star_alignment/{sample_id}/Log.final.out",
        "results/RNAseq/star_alignment/{sample_id}/SJ.out.tab"
    input:
        fastq_input,
        genomeDir = config['star_alignment']['genomeDir']
    params:
        readfiles = readfiles_param,
        common_args = star_common_args
    threads: 18
    conda: "../envs/star.yaml"
    shell:
        '''
tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.sample_id}.XXXXXX)

STAR\
  --readFilesIn {params.readfiles}\
  --outFileNamePrefix $(dirname {output[0]})/\
  --runThreadN {threads}\
  {params.common_args}
        '''

#The reason we check that all the runids from one sample have the same layout is because STAR
#expects that all the reads from a paired-end sample are in separate files 
#(one for each mate)
#while all the reads from a single-end sample are in a single file
#If the layout is mixed, with some runs being paired-end and others single-end,
#the STAR command would be tricky,
#as the reads from the paired-end runs will not be correctly paired.
