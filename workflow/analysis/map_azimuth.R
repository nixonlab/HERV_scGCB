#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(SeuratData)
    library(SeuratDisk)
})

################################################################################
#   Parse arguments
################################################################################
if(exists("snakemake")) {
    dataset <- snakemake@wildcards[['dataset']]
    samp <- snakemake@wildcards[['samp']]
    exptag <- snakemake@wildcards[['exptag']]

    infile <- snakemake@input[[1]]
    outh5 <- snakemake@output[[1]]
    
    refname <- snakemake@params[["refname"]]
    workers <- snakemake@threads
} else {
    cat("Not using snakemake\n")
    dataset <- 'greenleaf'
    samp <- 'Tonsil_1a'
    exptag <- 'pseudobulk'
    
    infile <- 'results/pbmc3p/analysis/pbmc3p.500/pseudobulk.orig.h5seurat'
    outh5 <- 'results/pbmc3p/analysis/pbmc3p.500/pseudobulk.map.h5seurat'
    
    refname <- 'tonsilref'
    workers <- 1
}

if(workers > 1) {
    cat('   using', workers, 'workers\n')
    library(future)
    plan("multicore", workers = workers)
}

outdir <- dirname(outh5)
dir.create(outdir, recursive=TRUE, showWarnings = FALSE)

################################################################################
#   Main
################################################################################
#--- Load reference
cat("Loading reference data...\n")
SeuratData::InstallData(refname)
reference <- SeuratData::LoadData(refname, type = "azimuth")$map
annotation.levels <- names(slot(object = reference, name = "meta.data"))
annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
do.adt <- "ADT" %in% names(reference@assays)

if(do.adt) cat("    ...including ADT")

#--- Load Seurat object
cat("Loading Seurat object...\n")
hfile <- SeuratDisk::Connect(infile)
hfile$index()

sobj <- SeuratDisk::LoadH5Seurat(hfile)

# Run Azimuth
cat("Running Azimuth...\n")

### This was needed because of expired https cert
# options(url.method='libcurl')
# httr::set_config(httr::config(ssl_verifypeer = 0L))

print(sobj.mapped <- Azimuth::RunAzimuth(
    sobj,
    reference = refname,
    do.adt = do.adt
))

# Adding Azimuth results using Seurat::AddAzimuthResults requires
# an Azimuth mapping scores file, which is an RDS file containing a list
# with:
#   results$umap : Dimensionality reduction (UMAP) as a `DimReduc` object
#   results$pred.df : `data.frame` with predicted celltypes and scores
#   results$impADT :  `Assay` object with imputed values for antibody derived
#                      tags (optional)
results <- list()

if('impADT' %in% Seurat::Assays(sobj.mapped)) {
    results$impADT <- sobj.mapped[['impADT']]
}

if('ref.umap' %in% Seurat::Reductions(sobj.mapped)) {
    results$umap <- sobj.mapped[['ref.umap']]
}

results$pred.df <- sobj.mapped@meta.data %>%
    tibble::rownames_to_column('cell') %>%
    dplyr::select(
        cell,
        dplyr::starts_with('predicted.'),
        mapping.score
    ) %>%
    as.data.frame

az.outfile <- file.path(outdir, paste0(exptag, '.tmpAzimuthResults.rds'))
saveRDS(results, file = az.outfile)

# Add Azimuth results
sobj <- Seurat::AddAzimuthResults(sobj, filename=az.outfile)
file.remove(az.outfile)

# Predicted celltype, score, and mapping.score columns are added to the
# metadata but are all NA (bug?). Add it manually here
stopifnot(all(rownames(sobj@meta.data) == results$pred.df$cell))
for (n in annotation.levels) {
    pcol <- paste0('predicted.', n)
    scol <- paste0('predicted.', n, '.score')
    sobj@meta.data[[pcol]] <- results$pred.df[[pcol]]
    sobj@meta.data[[scol]] <- results$pred.df[[scol]]
}
sobj@meta.data[['mapping.score']] <- results$pred.df[['mapping.score']]

# Update assays
for (n in annotation.levels) {
    aname <- paste0('prediction.score.', n)
    sobj[[aname]] <- sobj.mapped@assays[[aname]]
}

# Save Seurat object with mapping information
SeuratDisk::SaveH5Seurat(sobj, filename=outh5, overwrite=T)

cat("#--- azimuth_mapping.R completed. ---#\n")
