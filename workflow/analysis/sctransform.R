#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(SeuratDisk)
})

################################################################################
#   Utility Functions
################################################################################
rv_preserve_ltr <- function(
    the.seurat, 
    assay='RNA', 
    target.ltr=120, 
    default.rv.th = 0.9
) {
    # Run SCTransform with default threshold to get residual variance
    the.seurat <- SCTransform(the.seurat,
                              assay=assay,
                              new.assay.name=gsub('^RNA','SCT', assay),
                              vst.flavor="v2",
                              variable.features.n = NULL,
                              return.only.var.genes = FALSE
    )
    # get the SCT model
    sct.model <- the.seurat[[gsub('^RNA','SCT', assay)]]@SCTModel.list[[1]]
    
    # Get residual variance of LTRs
    mf <- the.seurat[[assay]]@meta.features
    tmp <- sct.model@feature.attributes %>%
        rownames_to_column('feature_name') %>%
        mutate(feature.class = mf[feature_name, 'feature.class']) %>%
        filter(feature.class == 'LTR')
    
    # Return if not enough LTR
    if(nrow(tmp)<target.ltr) {
        cat('Not enough LTR features to meet threshold\n')
        return(default.rv.th)
    }
    tmp %<>%
        slice_max(residual_variance, n=target.ltr) %>%
        arrange(residual_variance)
    
    top.rv <- tmp$residual_variance
    cat("Using residual variance threshold: ", top.rv[1], '\n')
    return(top.rv[1])
}

################################################################################
#   Parse arguments
################################################################################
if(exists("snakemake")) {
    infile <- snakemake@input[[1]]
    outfile <- snakemake@output[[1]]
    workers <- snakemake@threads
} else {
    args = commandArgs(trailingOnly=TRUE)
    usage <- 'USAGE:\n  sctransform.R [-h] input.h5seurat [output.h5seurat]'
    if(length(args)==0 | length(args) > 2)
        stop(usage, call.=FALSE)
    if(args[1] == '-h' | args[1] == '--help') {
        cat(usage, '\n')
        q("no")
    }
    
    infile <- args[1]
    if(length(args) == 2) {
        outfile <- args[2]
    } else {
        outfile <- ifelse(
            grepl('.h5seurat$', infile, ignore.case=T),
            gsub('.h5seurat$', '.sct.h5seurat', infile, ignore.case=T),
            gsub('$', '.sct.h5seurat', infile, ignore.case=T)
        )
    }
}

if(workers > 1) {
    cat('   using', workers, 'workers\n')
    library(future)
    plan("multicore", workers = workers)
}
################################################################################
#   Main
################################################################################
hfile <- SeuratDisk::Connect(infile)
hfile$index()

rna_assays <- hfile[['assays']]$names[grepl('^RNA', hfile[['assays']]$names)]
sobj <- SeuratDisk::LoadH5Seurat(hfile,
                                 assays=rna_assays,
                                 reductions=FALSE,
                                 graphs=FALSE,
                                 neighbors=FALSE,
                                 images=FALSE
)

# Get a threshold for residual variance that preserves LTR features
threshold.var <- rv_preserve_ltr(sobj, target.ltr=120)

for(the.assay in rna_assays) {
    cat('SCTransform for assay:', the.assay, '...\n')
    is_small <- nrow(sobj[[the.assay]]) < 3000
    if(is_small) {
        sobj <- SCTransform(
            sobj,
            assay=the.assay, 
            new.assay.name=gsub('^RNA','SCT', the.assay),
            vst.flavor="v2",
            variable.features.n = 3000,
            return.only.var.genes = FALSE
        )
    } else {
        sobj <- SCTransform(
            sobj,
            assay=the.assay, 
            new.assay.name=gsub('^RNA','SCT', the.assay),
            vst.flavor="v2",
            variable.features.n = NULL,
            variable.features.rv.th = threshold.var,
            return.only.var.genes = FALSE
        )
    }
}

sobj <- AppendData(infile, sobj)
DefaultAssay(sobj) <- 'SCT'
SaveH5Seurat(sobj, filename=outfile, overwrite=T)

cat("#--- sctransform.R completed. ---#\n")
