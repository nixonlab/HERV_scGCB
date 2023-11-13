#!/usr/bin/env python
# -*- coding: utf-8 -*-

rule complete_sample:
    input:
        expand('results/{{dataset}}/analysis/{{samp}}/{exptag}.sct.h5seurat',
            exptag=['pseudobulk','celltype','celltype.l1']
        )
    output:
        touch('results/{dataset}/completed/{samp}.txt')

