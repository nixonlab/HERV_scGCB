#! /usr/bin/env python
# -*- coding utf-8 -*-

def get_strand_cmd(wc):
    # input functions are passed a single parameter, wildcards
    # wildcards is a namespace so refer to the variable with a "."
    if wc.smode == "U": # handle the unstranded case
        return ""
    # all other options (R, F, RF, and FR) will just be passed as-is
    return f'--stranded_mode {wc.smode}' # this is a python f-string

rule stellarscope_pseudobulk:
    conda: "../envs/telescope.yaml"
    output:
        "results/telescope_pseudobulk/{samid}_{smode}/{samid}_{smode}_pseudobulk-TE_counts.mtx"
    input:
        bam = rules.stellarscope_cellsort.output,
        annotation = rules.telescope_annotation.output,
        barcodes = rules.starsolo_alignment.output[1]
    benchmark:
        "benchmarks/telescope_pseudobulk/{samid}_{smode}_telescope_pseudobulk.tsv"
    log:
        "results/telescope_pseudobulk/{samid}_{smode}/telescope.log"
    threads:
        config['telescope_threads']
    wildcard_constraints:
        smode="U|R|F|RF|FR"
    params:
        tmpdir = config['local_tmp'],
        out = "results/telescope_pseudobulk/{samid}_{smode}",
        exp_tag = "{samid}_{smode}_pseudobulk",
        stranded_mode = get_strand_cmd # the function name. this is a "callable" object that takes 1 argument, wildcards.
    shell:
        '''
stellarscope\
 assign\
 --exp_tag {params.exp_tag}\
 --updated_sam\
 --outdir {params.out}\
 --use_every_reassign_mode\
 {params.stranded_mode}\
 --whitelist {input.barcodes}\
 --pooling_mode pseudobulk\
 {input.bam}\
 {input.annotation}\
 2>&1 | tee {log[0]}
	'''

# when you request targets then you can specify the reps or stranded mode in the filename
# Below will launch 5 runs of stellarscope. Three runs will use no argument for stranded
# mode and two will use --stranded_mode F. The outputs will all be to different
# directories so we can compare.

rule telescope_complete:
    output:
        touch("results/completed/{samid}_completed.txt")
    input:
        "results/telescope_pseudobulk/{samid}_U/{samid}_U_pseudobulk-TE_counts.mtx",
        "results/telescope_pseudobulk/{samid}_F/{samid}_F_pseudobulk-TE_counts.mtx"
