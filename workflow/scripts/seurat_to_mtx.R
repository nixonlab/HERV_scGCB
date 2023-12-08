#! /usr/bin/env Rscript --vanilla

library(Seurat)

################################################################################
#   Parse arguments
################################################################################
if(exists("snakemake")) {
    seurat.rds <- snakemake@input[['seurat_rds']]
    counts.mtx <- snakemake@output[['counts_mtx']]
    features.tsv <- snakemake@output[['features_tsv']]
    barcodes.tsv <- snakemake@output[['barcodes_tsv']]
    
    assay.str <- ifelse('assay' %in% names(snakemake@params),
        snakemake@params[['assay']], 'RNA'
    )
    slot.str <- ifelse('slot' %in% names(snakemake@params),
        snakemake@params[['slot']], 'counts'
    )
} else {
    args = commandArgs(trailingOnly=TRUE)
    usage <- 'USAGE:\n  seurat_to_mtx.R [-h] seurat.rds [output.h5seurat]'
    if(length(args)==0 | length(args) > 2)
        stop(usage, call.=FALSE)
    if(args[1] == '-h' | args[1] == '--help') {
        cat(usage, '\n')
        q("no")
    }
}

################################################################################
#   Main
################################################################################
sobj <- readRDS(seurat.rds)
countmat <- slot(sobj[[assay.str]], slot.str)

outdir <- dirname(counts.mtx)
dir.create(outdir, recursive=TRUE, showWarnings = FALSE)

Matrix::writeMM(countmat, file=counts.mtx)
write(rownames(countmat), file=features.tsv)
write(colnames(countmat), file=barcodes.tsv)