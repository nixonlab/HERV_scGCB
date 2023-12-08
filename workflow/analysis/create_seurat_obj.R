#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(SeuratDisk)
})

source('workflow/analysis/01-helpers.R')

################################################################################
#   Parse input
################################################################################
if(exists("snakemake")) {
    dataset <- snakemake@wildcards[['dataset']]
    samp <- snakemake@wildcards[['samp']]
    exptag <- snakemake@wildcards[['exptag']]

    gtype.rds <- snakemake@params[["gtype_rds"]]
    retro.rds <- snakemake@params[["retro_rds"]]
    remove.nofeat <- snakemake@params[['remove_nofeat']]
    min.cells <- snakemake@params[['min_cells']]
    min.features <- snakemake@params[['min_features']]

    outdir.star <- dirname(snakemake@input[["star_counts_mtx"]])
    outdir.stellar <- dirname(snakemake@input[["stellar_counts_mtx"]])

    outh5 <- snakemake@output[[1]]
    
    workers <- snakemake@threads
} else {
    cat("Not using snakemake\n")
    dataset <- 'greenleaf'
    samp <- 'Tonsil_1a_gex'
    exptag <- 'pseudobulk'

    gtype.rds <- 'resources/gencode.v38/metadata.gid_gtype.rds'
    retro.rds <- 'resources/retro.hg38.v1.rds'
    remove.nofeat <- TRUE
    min.cells <- 0
    min.features <- 0
    
    outdir.star <- file.path('results', dataset, 'align_multi_starsolo', samp, 'Solo.out/Gene/filtered')
    outdir.stellar <- file.path('results', dataset, 'stellarscope', samp)

    outh5 <- file.path('results', dataset, 'downstream', samp, 'orig.h5seurat')
    workers <- 1
}

if(workers > 1) {
    cat('   using', workers, 'workers\n')
    library(future)
    plan("multicore", workers = workers)
}

################################################################################
#   Main
################################################################################
#--- Load count data as sparse matrix
cat("Loading count data...\n")
gmat <- Seurat::ReadMtx(
    mtx = file.path(outdir.star, 'matrix.mtx'),
    cells = file.path(outdir.star, 'barcodes.tsv'),
    features = file.path(outdir.star, 'features.tsv')
)

tmat <- Seurat::ReadMtx(
    mtx = file.path(outdir.stellar, paste0(exptag, '-TE_counts.mtx')),
    cells = file.path(outdir.stellar, paste0(exptag, '-barcodes.tsv')),
    features = file.path(outdir.stellar, paste0(exptag, '-features.tsv')),
    feature.column=1
)

#--- Remove no feature
if(remove.nofeat) {
    rem <- which(rownames(tmat) == '__no_feature')
    if(length(rem)>0) tmat <- tmat[-rem, ]
}

#--- Change underscore to hyphen
rownames(gmat) <- gsub('_', '-', rownames(gmat))
rownames(tmat) <- gsub('_', '-', rownames(tmat))

#--- Harmonize columns (barcodes)
bcint <- intersect(colnames(gmat), colnames(tmat))
gmat <- gmat[ , bcint]
tmat <- tmat[ , bcint]

stopifnot(all(colnames(gmat) == colnames(tmat)))

#--- Feature table
cat("Loading feature metadata...\n")
gid_gtype <- readRDS(gtype.rds)
gfeat <- read.table(
        file.path(outdir.star, 'features.tsv'),
        sep='\t',
        col.names = c('id', 'name', 'assay')
    ) %>%
    dplyr::mutate(
        feature.class = 'CG',
        feature.gtype = gid_gtype[id, 1],
        name.unique = gsub('_', '-', rownames(gmat))
    ) %>%
    dplyr::select(name.unique, id, name, feature.class, feature.gtype, assay) %>%
    tibble::column_to_rownames("name.unique")

tfeat <- readRDS(retro.rds) %>%
    dplyr::mutate(
        name.unique = gsub('_', '-', locus),
        id = locus,
        name = locus,
        feature.class = te_class,
        feature.gtype = family,
        assay = 'Gene Expression'        
    ) %>%
    dplyr::select(name.unique, id, name, feature.class, feature.gtype, assay) %>%
    tibble::column_to_rownames("name.unique")

tfeat <- tfeat[rownames(tmat), ]

features <- rbind(gfeat, tfeat) %>% dplyr::mutate(dplyr::across(3:5, as.factor))
rm(gfeat, tfeat)

#--- Combine matrices
cat("Creating Seurat object...\n")
counts <- rbind(gmat, tmat)
stopifnot(all(rownames(counts) == rownames(features)))

sobj <- Seurat::CreateSeuratObject(
    counts,
    project=paste0(samp, '.', exptag),
    min.cells=min.cells,
    min.features=min.features
)
sobj[['RNA']] <- Seurat::AddMetaData(sobj[['RNA']], features)
sobj <- cell_summary_data(sobj)

stopifnot(all(rownames(sobj) == rownames(counts)))
stopifnot(all(colnames(sobj) == colnames(counts)))

#--- Create output
dir.create(dirname(outh5), recursive=TRUE, showWarnings = FALSE)
SeuratDisk::SaveH5Seurat(sobj, filename=outh5, overwrite=T)

cat("#--- create_seurat_obj.R completed. ---#\n")
