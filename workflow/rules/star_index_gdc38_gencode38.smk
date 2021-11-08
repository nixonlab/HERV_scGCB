#! /usr/bin/env python
# -*- coding utf-8 -*-

rule star_index_gdc38_gencode38:
    conda:
        "../envs/star.yaml"
    input:
        genome = config['sequences']['genome'],
        annotation_gtf_gz = config['annotations']['gencode']
    output:
        directory(config['indexes']['star'])
    params:
        sjdbOverhang = config['splice_junction_overhang']
    threads: config['star_index_threads']
    resources:
        mem_mb=config['star_index_mem_mb']
    shell:
        """
        tdir=$(mktemp -d {config[local_tmp]}/{rule}.XXXXXX)

        pigz -dc {input.annotation_gtf_gz} > $tdir/gencode.v38.annotation.gtf
        pigz -dc {input.genome} > $tdir/genome.fa

        STAR\
            --runThreadN {threads}\
            --runMode genomeGenerate\
            --genomeDir {output}\
            --outFileNamePrefix {output}\
            --genomeFastaFiles $tdir/genome.fa\
            --sjdbGTFfile $tdir/gencode.v38.annotation.gtf\
            --sjdbOverhang {params.sjdbOverhang}
        """
