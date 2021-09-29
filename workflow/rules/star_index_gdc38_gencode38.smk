rule star_index_gdc38_gencode38:
    input:
        fasta_gz = 'databases/remotefiles/GRCh38.d1.vd1.fa.tar.gz',
        gtf_gz = 'databases/remotefiles/gencode.v38.annotation.gtf.gz'
    output:
        directory("databases/star_index_GDCHG38_gencode38")
    params:
        sjdbOverhang = config['splice_junction_overhang']
    conda:
        "../envs/star.yaml"
    threads: workflow.cores
    shell:
        '''
        tdir=$(mktemp -d {config[local_tmp]}/{rule}.XXXXXX)

        tar -xvzf {input.fasta_gz} -C $tdir

	pigz -dc {input.gtf_gz} > $tdir/gencode.v38.annotation.gtf

        STAR\
            --runThreadN {threads}\
            --runMode genomeGenerate\
            --genomeDir {output}\
            --outFileNamePrefix {output}\
            --genomeFastaFiles $tdir/GRCh38.d1.vd1.fa\
            --sjdbGTFfile $tdir/gencode.v38.annotation.gtf\
            --sjdbOverhang {params.sjdbOverhang}
            '''
