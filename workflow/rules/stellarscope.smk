#! /usr/bin/env python
# -*- coding utf-8 -*-

reassign_modes = [
    "best_exclude",
    "best_conf",
    "best_random",
    "best_average",
    "initial_unique",
    "initial_random",
    "total_hits"
]

wildcard_constraints:
    pmode = "pseudobulk|individual|celltype"

rule stellarscope_cellsort:
    output:
        protected("results/{dataset}/stellarscope/{samp}/Aligned.sortedByCB.bam")
    input:
        aln_bam = rules.align_multi_starsolo.output[0],
        passBC_tsv = rules.scrublet.output[2]
    conda:
        "../envs/stellarscope_dev.yaml"
    threads: 12
    params:
        tempdir = config['tmpdir']
    shell:
        '''
        mkdir -p $(dirname {output[0]})

        stellarscope cellsort\
            --nproc {threads}\
            --tempdir {params.tempdir}\
            --outfile {output[0]}\
            {input.aln_bam}\
            {input.passBC_tsv}
        '''

def get_strand_arg(wildcards):
    stranded = lambda wc: samples.loc[wildcards.samp]['stranded']
    if stranded == "U": # handle the unstranded case
        return ""
    return f'--stranded_mode {stranded}'

rule stellarscope_load:
    output:
        "results/{dataset}/stellarscope/{samp}/stload-checkpoint.load_alignment.pickle",
        "results/{dataset}/stellarscope/{samp}/stload-checkpoint.dedup_umi.pickle",
        "results/{dataset}/stellarscope/{samp}/stload-stats.final.tsv"
    input:
        bamfile = rules.stellarscope_cellsort.output,
        passBC_tsv = rules.scrublet.output[2],
        gtffile = config['stellarscope']['retro_gtf']
    conda:
        "../envs/stellarscope_dev.yaml"
    benchmark:
        "benchmarks/{dataset}/{samp}.stload.tsv"
    log:
        "results/{dataset}/stellarscope/{samp}/stload.log"
    threads: 1
    params:
        strandarg = get_strand_arg 
    shell:
        '''
        stellarscope assign\
            --exp_tag stload\
            --outdir $(dirname {output[0]})\
            --skip_em\
            --updated_sam\
            {params.strandarg}\
            --whitelist {input.passBC_tsv}\
            {input.bamfile}\
            {input.gtffile}\
            2>&1 | tee {log[0]}
        '''
 
rule stellarscope_resume:
    output:
        "results/{dataset}/stellarscope/{samp}/{pmode}-TE_counts.mtx",
        "results/{dataset}/stellarscope/{samp}/{pmode}-barcodes.tsv",
        "results/{dataset}/stellarscope/{samp}/{pmode}-features.tsv",
        "results/{dataset}/stellarscope/{samp}/{pmode}-checkpoint.final.pickle",
        "results/{dataset}/stellarscope/{samp}/{pmode}-stats.final.tsv",
        expand("results/{{dataset}}/stellarscope/{{samp}}/{{pmode}}-TE_counts.{rmode}.mtx",
            rmode = reassign_modes[1:]
        )
    input:
        rules.stellarscope_load.output[1],
        rules.annotate_celltypes_azimuth.output["tonsilL2_tsv"]
    conda:
        "../envs/stellarscope_dev.yaml"
    benchmark:
        "benchmarks/{dataset}/{samp}.{pmode}.stellarscope_resume.tsv"
    log:
        "results/{dataset}/stellarscope/{samp}/{pmode}.log"
    threads: 12
    params:
        strandarg = get_strand_arg
    shell:
        '''
        stellarscope resume\
            --exp_tag {wildcards.pmode}\
            --outdir $(dirname {output[0]})\
            --pooling_mode {wildcards.pmode}\
            --celltype_tsv {input[1]}\
            --updated_sam\
            --use_every_reassign_mode\
            --nproc {threads}\
            --max_iter 500\
            {input[0]}\
            2>&1 | tee {log[0]}
        '''

# Since we are hardcoding it really didnt make sense to leave the L3 for tonsil here,
# but maybe in the future we will work with some reference that has up to L3
def get_celltype_tsv(wildcards):
    if wildcards.lvl == 'l1':
        return rules.annotate_celltypes_azimuth.output["tonsilL1_tsv"]
    elif wildcards.lvl == 'l3':
        return rules.annotate_celltypes_azimuth.output["tonsilL3_tsv"]    

rule stellarscope_resume_celltype:
    """
    Run stellarscope resume for other celltype annotation levels (l1 and l3)
    """
    output:
         "results/{dataset}/stellarscope/{samp}/celltype.{lvl}-TE_counts.mtx",
         "results/{dataset}/stellarscope/{samp}/celltype.{lvl}-barcodes.tsv",
         "results/{dataset}/stellarscope/{samp}/celltype.{lvl}-features.tsv",
         "results/{dataset}/stellarscope/{samp}/celltype.{lvl}-checkpoint.final.pickle",
         "results/{dataset}/stellarscope/{samp}/celltype.{lvl}-stats.final.tsv",
         expand("results/{{dataset}}/stellarscope/{{samp}}/celltype.{{lvl}}-TE_counts.{rmode}.mtx",
                 rmode = reassign_modes[1:]
        )
    input:
        rules.stellarscope_load.output[1],
        get_celltype_tsv
    wildcard_constraints:
        lvl = "(l1|l3)"
    conda:
        "../envs/stellarscope_dev.yaml"
    benchmark:
        "benchmarks/{dataset}/{samp}.celltype.{lvl}.stellarscope_resume.tsv"
    log:
        "results/{dataset}/stellarscope/{samp}/celltype.{lvl}.log"
    threads: 12
    params:
        strandarg = get_strand_arg
    shell:
        '''
        stellarscope resume\
            --exp_tag celltype.{wildcards.lvl}\
            --outdir $(dirname {output[0]})\
            --pooling_mode celltype\
            --celltype_tsv {input[1]}\
            --updated_sam\
            --use_every_reassign_mode\
            --nproc {threads}\
            --max_iter 500\
            {input[0]}\
            2>&1 | tee {log[0]}
        '''

localrules: stellarscope_poolmodes
rule stellarscope_poolmodes:
    output:
        "results/{dataset}/stellarscope/{samp}/poolmodes.txt"
    input:
        "results/{dataset}/stellarscope/{samp}/individual-TE_counts.mtx",
        "results/{dataset}/stellarscope/{samp}/pseudobulk-TE_counts.mtx",
        "results/{dataset}/stellarscope/{samp}/celltype-TE_counts.mtx",
        "results/{dataset}/stellarscope/{samp}/celltype.l1-TE_counts.mtx"
    shell:
        '''
        for f in {input}; do
        echo $(basename $f) | sed 's/-TE_counts.mtx//'
        done > {output[0]}
        '''

