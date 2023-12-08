#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(scater)
})

################################################################################
#   Utility Functions
################################################################################
cell_qc <- function(the.seurat) {
    require(Seurat)
    require(scater)
    
    the.seurat[['percent.mt']] <- Seurat::PercentageFeatureSet(
        the.seurat,
        features=rownames(the.seurat)[grepl('^MT-', rownames(the.seurat))]
    )
    
    qc.ncount_rna <- scater::isOutlier(
        the.seurat$nCount_RNA,
        log = TRUE,
        type = "both"
    )
    qc.nfeature_rna <- scater::isOutlier(
        the.seurat$nFeature_RNA,
        log = TRUE,
        type = "both"
    )
    qc.percent_mt <- scater::isOutlier(
        the.seurat$percent.mt,
        type="higher"
    )
    
    thresh <- data.frame(
        ncount = attr(qc.ncount_rna, "thresholds"),
        nfeature = attr(qc.nfeature_rna, "thresholds"),
        mt = attr(qc.percent_mt, "thresholds")
    )
    
    subset(
        the.seurat,
        subset = nCount_RNA > thresh["lower", "ncount"] &
            nCount_RNA < thresh["higher", "ncount"] &
            nFeature_RNA >  thresh["lower", "nfeature"] &
            nFeature_RNA < thresh["higher", "nfeature"] &
            percent.mt < thresh["higher", "mt"]
    )
}

################################################################################
#   Main Function
################################################################################
cell_qc.run <- function(
    counts.mtx, barcodes.tsv, features.tsv, outdir
) {
    # Load gene data as Seurat object
    cat("Loading count matrix...\n")
    print(sobj <- Seurat::CreateSeuratObject(
        Seurat::ReadMtx(
            mtx = counts.mtx,
            cells = barcodes.tsv,
            features = features.tsv
        )
    ))
    
    # Subset cells by RNA count, feature count, and MT content
    print("Filtering count matrix...")
    print(sobj <- cell_qc(sobj))
    
    # Save Seurat object with mapping information
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(sobj, file.path(outdir, 'tx_seurat.cell_qc.rds'))
    
    # Write passing barcodes to file
    outfile <- file.path(outdir, 'barcodes.pass.tsv')
    rownames(sobj@meta.data) %>%
        write.table(file=outfile, quote=F, row.names=F, col.names=F)
}


################################################################################
#   Parse arguments
################################################################################
if (sys.nframe() == 0L) {
    if(exists("snakemake")) {
        counts.mtx <- snakemake@input[["counts_mtx"]]
        barcodes.tsv <- snakemake@input[["barcodes_tsv"]]        
        features.tsv <- snakemake@input[["features_tsv"]]
        outdir <- dirname(snakemake@output[[1]])
    } else {
        args = commandArgs(trailingOnly=TRUE)
        usage <- 'USAGE:\n  cell_qc.R [-h] counts.dir outdir'
        if(args[1] == '-h' | args[1] == '--help') {
            cat(usage, '\n')
            q("no")
        }
        if(length(args) != 2)
            stop(usage, call.=FALSE)
        
        counts.mtx <- file.path(args[1], 'matrix.mtx')
        barcodes.tsv <- file.path(args[1], 'barcodes.tsv')
        features.tsv <- file.path(args[1], 'features.tsv')
        stopifnot(file.exists(counts.mtx, barcodes.tsv, features.tsv))
        outdir <- args[2]
    }
    cell_qc.run(counts.mtx, barcodes.tsv, features.tsv, outdir)
    cat("#--- cell_qc.R completed. ---#\n")
}
