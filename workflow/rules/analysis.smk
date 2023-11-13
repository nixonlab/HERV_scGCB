#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Downstream analysis """

wildcard_constraints:
    exptag = "pseudobulk|individual|celltype(.l\d)*"

rule create_seurat_obj:
    output:
        'results/{dataset}/analysis/{samp}/{exptag}.orig.h5seurat'
    input:
        stellar_counts_mtx = "results/{dataset}/stellarscope/{samp}/{exptag}-TE_counts.mtx",
        stellar_barcodes_tsv = "results/{dataset}/stellarscope/{samp}/{exptag}-barcodes.tsv",
        stellar_features_tsv = "results/{dataset}/stellarscope/{samp}/{exptag}-features.tsv",
        star_counts_mtx = rules.align_multi_starsolo.output[1],
        star_barcodes_tsv = rules.align_multi_starsolo.output[2],
        star_features_tsv = rules.align_multi_starsolo.output[3]
    params:
        gtype_rds = config['stellarscope']['gtype_rds'],
        retro_rds = config['stellarscope']['retro_rds'],
        remove_nofeat = True,
        min_cells = 0,
        min_features = 0
    threads: 1
    container:
        'docker://hreypar/scope-r:latest'
    script:
        '../analysis/create_seurat_obj.R'

rule map_azimuth:
    output:
        'results/{dataset}/analysis/{samp}/{exptag}.map.h5seurat'
    input:
        'results/{dataset}/analysis/{samp}/{exptag}.orig.h5seurat'
    params:
        refname = lambda wc: samples.loc[wc.samp]['azimuth_ref']
    threads: 1
    container:
        'docker://hreypar/scope-r:latest'
    script:
        '../analysis/map_azimuth.R'

rule aggregate_features:
    output:
        'results/{dataset}/analysis/{samp}/{exptag}.agg.h5seurat'
    input:
        'results/{dataset}/analysis/{samp}/{exptag}.map.h5seurat'
    params:
        refname = lambda wc: samples.loc[wc.samp]['azimuth_ref']
    threads: 1
    container:
        'docker://hreypar/scope-r:latest'
    script:
        '../analysis/aggregate_features.R'

rule sctransform:
    output:
        'results/{dataset}/analysis/{samp}/{exptag}.sct.h5seurat'
    input:
        'results/{dataset}/analysis/{samp}/{exptag}.agg.h5seurat'
    params:
        refname = lambda wc: samples.loc[wc.samp]['azimuth_ref']
    threads: 1
    container:
        'docker://hreypar/scope-r:latest'
    script:
        '../analysis/sctransform.R'

