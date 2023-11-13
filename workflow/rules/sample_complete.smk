#!/usr/bin/env python
# -*- coding: utf-8 -*-

s_rna_gex = dsets.loc['rna_gex_PRJNA587486']['sample']
s_cite_gex = dsets.loc['cite_gex_PRJNA587486']['sample']

rule complete_sample:
    input:
        expand('results/{dataset}/analysis/{samp}/{exptag}.sct.h5seurat',
            dataset = 'rna_gex_PRJNA587486',
            samp = s_rna_gex,
            exptag=['pseudobulk','celltype','celltype.l1']
        ),
        expand('results/{dataset}/analysis/{samp}/{exptag}.sct.h5seurat',
            dataset = 'cite_gex_PRJNA587486',
            samp = s_cite_gex,
            exptag=['pseudobulk','celltype','celltype.l1']
        )
    output:
        touch('results/{dataset}/completed.txt')

