#!/usr/bin/env python
# -*- coding: utf-8 -*-

rule cell_qc:
    """ 
    Quality control using canonical genes (STAR output)
    """
    output:
        passBC_tsv = 'results/{dataset}/canonicalQC/{samp}/barcodes.pass.tsv',
        seurat_rds = 'results/{dataset}/canonicalQC/{samp}/tx_seurat.cell_qc.rds'
    input:
        counts_mtx = rules.align_multi_starsolo.output[1],
        barcodes_tsv = rules.align_multi_starsolo.output[2],
        features_tsv = rules.align_multi_starsolo.output[3]
#        inst_deps = rules.install_deps_Renv.output
    container:
        'docker://hreypar/scope-r:latest'
    script:
        '../scripts/cell_qc.R'

rule annotate_celltypes_azimuth:
    output:
        tonsilL1_tsv = 'results/{dataset}/canonicalQC/{samp}/tonsilref.celltype.l1.tsv',
        tonsilL2_tsv = 'results/{dataset}/canonicalQC/{samp}/tonsilref.celltype.l2.tsv',
        mapped_seurat_rds = 'results/{dataset}/canonicalQC/{samp}/tx_seurat.mapped.rds' 
    input:
        seurat_rds = rules.cell_qc.output['seurat_rds']
#        inst_deps = rules.install_deps_Renv.output
    params:
        thresh = config['annotate_celltypes_azimuth']['filter_threshold'],
        refname = lambda wc: samples.loc[wc.samp]['azimuth_ref']
    container:
        'docker://hreypar/scope-r:latest'
    script:
        '../scripts/annotate_azimuth.R'

localrules: seurat_to_mtx_forscr
rule seurat_to_mtx_forscr:
    input:
        seurat_rds = rules.annotate_celltypes_azimuth.output['mapped_seurat_rds']
    output:
        counts_mtx = 'results/{dataset}/canonicalQC/{samp}/tx_seurat.mapped.counts.mtx',
        features_tsv = 'results/{dataset}/canonicalQC/{samp}/tx_seurat.mapped.features.tsv',
        barcodes_tsv = 'results/{dataset}/canonicalQC/{samp}/tx_seurat.mapped.barcodes.tsv'
    container:
        'docker://hreypar/scope-r:latest'
    script:
        '../scripts/seurat_to_mtx.R'

rule scrublet:
    input:
        counts_mtx = rules.seurat_to_mtx_forscr.output['counts_mtx'],
        features_tsv = rules.seurat_to_mtx_forscr.output['features_tsv'],
        barcodes_tsv = rules.seurat_to_mtx_forscr.output['barcodes_tsv']
    output:
        'results/{dataset}/canonicalQC/{samp}/scrublet/SCR.predictions_all.tsv',
        'results/{dataset}/canonicalQC/{samp}/scrublet/SCR.predictions_single.tsv',
        'results/{dataset}/canonicalQC/{samp}/scrublet/SCR.barcodes.pass.tsv',
        'results/{dataset}/canonicalQC/{samp}/scrublet/SCR.summary.txt'
    params:
        doublet_rate = 0.08
    conda:
        '../envs/scrublet.yaml'
    script:
        '../scripts/run_scrublet.py'

rule doubletfinder:
    input:
        seurat_rds = rules.annotate_celltypes_azimuth.output['mapped_seurat_rds']
    output:
         'results/{dataset}/canonicalQC/{samp}/doubletfinder/DF.predictions_all.tsv',
         'results/{dataset}/canonicalQC/{samp}/doubletfinder/DF.predictions_single.tsv',
         'results/{dataset}/canonicalQC/{samp}/doubletfinder/DF.barcodes.pass.tsv',
         'results/{dataset}/canonicalQC/{samp}/doubletfinder/DF.summary.tsv'        
    params:
        doublet_rate = 0.08,
        cell_idents = "predicted.celltype.l1"
    threads: 8
    container:
        'docker://hreypar/scope-r:latest'
    script:
        '../scripts/run_doubletfinder.R'

rule gather_doublets:
    input:
        rules.scrublet.output,
        rules.doubletfinder.output
    output:
        touch('results/{dataset}/canonicalQC/{samp}_doublets.completed')

