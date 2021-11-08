######################################################
## Analyze droplet single-cell RNA sequencing data
## with STARsolo
######################################################
## to do:
## input as function for cDNA and barcodes
## a scheme to name output directories
## see if I can integrate PEP here

#rule decompress_whitelist:
#	input: "databases/remotefiles/whitelist.10x.v2.txt.gz"
#	output: "databases/remotefiles/whitelist.10x.v2.txt"
#	shell:
#		"gunzip {input}"

rule starsolo_alignment:
	"""
	Align sequencing reads from a 10x V3 single-cell RNA-seq experiment using STARsolo
	"""
	conda:
		"../envs/star.yaml"
	input:
		cDNA = "runs/{run_acc}/{run_acc}_2.fastq",
		barcodes = "runs/{run_acc}/{run_acc}_1.fastq",
		genome = config['indexes']['star'],
		whitelist = "refs/downloads/whitelist.10x.v2.txt"
	output:
		"results/starsolo_algn/{run_acc}/{run_acc}_GDC38.Aligned.sortedByCoord.out.bam"
	params:
		out_prefix="results/starsolo_algn/{run_acc}/{run_acc}_GDC38.",
		cb_start=config["cellbarcode_start"],
		cb_length=config["cellbarcode_length"],
		umi_start=config["umi_start"],
		umi_length=config["umi_length"],
		max_multimap=config["max_multimap"]
    benchmark: "benchmarks/star_alignment/{run_acc}_star_alignment.tsv"
	threads: config['star_alignment_threads']
	shell:
		'''
		#--- STARsolo (turned on by --soloType CB_UMI_Simple)
		STAR\
			--runThreadN {threads}\
			--genomeDir {input.genome}\
			--readFilesIn {input.cDNA} {input.barcodes}\
			--soloType CB_UMI_Simple\
			--soloCBwhitelist {input.whitelist}\
			--soloCBstart {params.cb_start}\
			--soloCBlen {params.cb_length}\
			--soloUMIstart {params.umi_start}\
			--soloUMIlen {params.umi_length}\
			--outFilterMultimapNmax {params.max_multimap}\
			--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM\
			--outSAMtype BAM SortedByCoordinate\
			--outFileNamePrefix {params.out_prefix}
		'''
