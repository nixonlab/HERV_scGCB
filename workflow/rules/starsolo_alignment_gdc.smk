######################################################
## Analyze droplet single-cell RNA sequencing data
## with STARsolo
######################################################
## to do:
## input as function for cDNA and barcodes
## a scheme to name output directories
## see if I can integrate PEP here

rule decompress_whitelist:

rule starsolo_alignment:
	"""
	Align sequencing reads from a 10x V3 single-cell RNA-seq experiment using STARsolo
	"""
	input:
		cDNA = "runs/{run_acc}/{run_acc}_2.fastq",
		barcodes = "runs/{run_acc}/{run_acc}_1.fastq",
		genome = "databases/star_index_GDCHG38_gencode38",
		whitelist_gz = "databases/remotefiles/whitelist.10x.v3.gz"
	output:
		"results/{run_acc}/{run_acc}_GDC38.Aligned.sortedByCoord.out.bam"
	params:
		out_prefix="results/{run_acc}/{run_acc}_GDC38.",
		cb_start=config["cellbarcode_start"],
		cb_length=config["cellbarcode_length"],
		umi_start=config["umi_start"],
		umi_length=config["umi_length"],
		max_multimap=config["max_multimap"]
	conda:
		"../envs/star.yaml"
	threads: worflow.cores
	shell:
		'''
		pigz {input.whitelist_gz}
	
		#--- STARsolo (turned on by --soloType CB_UMI_Simple)
		STAR\
			--runThreadN {threads}\
			--genomeDir {input.genome}\
			--readFilesIn {input.cDNA} {input.barcodes}\
			--readFilesCommand gunzip -c\
			--soloType CB_UMI_Simple\
			--soloCBwhitelist databases/remotefiles/whitelist.10x.v3\
			--soloCBstart {params.cb_start}\
			--soloCBlen {params.cb_length}\
			--soloUMIstart {params.umi_start}\
			--soloUMIlen {params.umi_length}\
			--outFilterMultimapNmax {params.max_multimap}\
			--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM\
			--outSAMtype BAM SortedByCoordinate\
			--outFileNamePrefix {params.out_prefix}
		'''
