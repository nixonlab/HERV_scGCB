#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(SeuratData)
    library(Azimuth)
})


################################################################################
#   Utility Functions
################################################################################
azimuth_ref_info <- function(refname) {
    # Load reference
    SeuratData::InstallData(refname)
    reference <- SeuratData::LoadData(refname, type = "azimuth")$map
    annotation.levels <- names(slot(object = reference, name = "meta.data"))
    annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
    has.adt <- "ADT" %in% names(reference@assays)
    return(list(annotation.levels, has.adt))
}

my.AddAzimuthResults <- function(the.seurat, the.seurat.mapped, lvls, tempdir) {
    # Adding Azimuth results using Seurat::AddAzimuthResults requires
    # an Azimuth mapping scores file, which is an RDS file containing a list
    # with:
    #   results$umap : Dimensionality reduction (UMAP) as a `DimReduc` object
    #   results$pred.df : `data.frame` with predicted celltypes and scores
    #   results$impADT :  `Assay` object with imputed values for antibody derived
    #                      tags (optional)
    results <- list()

    if('impADT' %in% Seurat::Assays(the.seurat.mapped)) {
        results$impADT <- the.seurat.mapped[['impADT']]
    }

    if('ref.umap' %in% Seurat::Reductions(the.seurat.mapped)) {
        results$umap <- the.seurat.mapped[['ref.umap']]
    }

    results$pred.df <- the.seurat.mapped@meta.data %>%
        tibble::rownames_to_column('cell') %>%
        dplyr::select(
            cell,
            dplyr::starts_with('predicted.'),
            mapping.score
        ) %>%
        as.data.frame

    az.outfile <- file.path(tempdir, paste0('tmpAzimuthResults.rds'))
    saveRDS(results, file = az.outfile)

    # Add Azimuth results
    the.seurat <- Seurat::AddAzimuthResults(the.seurat, filename=az.outfile)
    file.remove(az.outfile)

    # # Another way to get annotation levels
    # lvls <- the.seurat@meta.data %>%
    #     dplyr::select(dplyr::starts_with('predicted.')) %>%
    #     dplyr::select(-dplyr::ends_with('.score')) %>%
    #     names %>%
    #     str_remove('predicted.')

    # Predicted celltype, score, and mapping.score columns are added to the
    # metadata but are all NA (bug?). Add it manually here
    stopifnot(all(rownames(the.seurat@meta.data) == results$pred.df$cell))
    for (n in lvls) {
        pcol <- paste0('predicted.', n)
        scol <- paste0('predicted.', n, '.score')
        the.seurat@meta.data[[pcol]] <- results$pred.df[[pcol]]
        the.seurat@meta.data[[scol]] <- results$pred.df[[scol]]
    }
    the.seurat@meta.data[['mapping.score']] <- results$pred.df[['mapping.score']]

    # Update assays
    for (n in lvls) {
        aname <- paste0('prediction.score.', n)
        the.seurat[[aname]] <- the.seurat.mapped@assays[[aname]]
    }

    # Return Seurat object with mapping information
    return(the.seurat)
}

################################################################################
#   Main Function
################################################################################

annotate_azimuth.run <- function(seurat.rds, refname, score.th, outdir) {
    # Load seurat object
    sobj <- readRDS(seurat.rds)
    
    # Get reference info
    refinfo <- azimuth_ref_info(refname)
    annotation.levels <- refinfo[[1]]
    do.adt <- refinfo[[2]]
    
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
    
    # Add Azimuth results to Seurat object
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    sobj <- my.AddAzimuthResults(sobj, sobj.mapped, annotation.levels, outdir)
    saveRDS(sobj, file.path(outdir, 'tx_seurat.mapped.rds'))
    
    # Output celltype annotation files
    for (n in annotation.levels) {
        outfile <- file.path(
            outdir,
            paste0(refname, '.', n, '.tsv')
        )
        sobj@meta.data %>% tibble::rownames_to_column('cell') %>%
            dplyr::select(cell, paste0('predicted.', n)) %>%
            write_tsv(file=outfile, quote='none', col_names = F)
    }
    
    # Output filtered celltype annotation files
    for (n in annotation.levels) {
        outfile <- file.path(
            outdir,
            paste0(refname, '.', n, '.filt.tsv')
        )
        sobj@meta.data %>%
            dplyr::filter(.data[[paste0('predicted.', n, ".score")]] > score.th) %>%
            tibble::rownames_to_column('cell') %>%
            dplyr::select(cell, paste0('predicted.', n)) %>%
            write_tsv(file=outfile, quote='none', col_names = F)
    }
}
    
################################################################################
#   Parse arguments
################################################################################
if (sys.nframe() == 0L) {
    if(exists("snakemake")) {
        seurat.rds <- snakemake@input[["seurat_rds"]]
        refname <- snakemake@params[['refname']]
        score.th <- snakemake@params[['thresh']]
        outdir <- dirname(snakemake@output[[1]])
    } else {
        args = commandArgs(trailingOnly=TRUE)
        usage <- 'USAGE:\n  annotate_azimuth.R [-h] seurat.rds refname outdir'
        if(args[1] == '-h' | args[1] == '--help') {
            cat(usage, '\n')
            q("no")
        }
        if(length(args) != 3)
            stop(usage, call.=FALSE)
        
        seurat.rds <- args[1]
        refname <- args[2]
        score.th <- 0.8
        outdir <- args[3]
        stopifnot(file.exists(seurat.rds))
    }
    annotate_azimuth.run(seurat.rds, refname, score.th, outdir)
    cat("#--- annotate_azimuth.R completed. ---#\n")
}

# 
#     
# if(exists("snakemake")) {
#     outdir.star <- dirname(snakemake@input[["counts_mtx"]])
#     refname <- snakemake@params[["refname"]]
#     outdir <- dirname(snakemake@output[[1]])
#     thresh <- snakemake@params[["thresh"]]
# } else {
#     cat("Not using snakemake\n")
#     outdir.star <- 'results/pbmc3p/align_multi_starsolo/pbmc3p.500/Solo.out/Gene/filtered'
#     refname <- 'pbmcref'
#     outdir <- 'results/pbmc3p/celltype_annotation/pbmc3p.500'
#     thresh <- 0.8
# }
# 
#     
#     
#     
#     
# }
# # Load reference
# SeuratData::InstallData(refname)
# reference <- SeuratData::LoadData(refname, type = "azimuth")$map
# annotation.levels <- names(slot(object = reference, name = "meta.data"))
# annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
# annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
# annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
# do.adt <- "ADT" %in% names(reference@assays)
# 
# if(do.adt) print("Including ADT")

dir.create(outdir, recursive=TRUE, showWarnings = FALSE)

# # Load gene data as Seurat object
# print("Loading count matrix...")
# print(sobj <- Seurat::CreateSeuratObject(
#     Seurat::ReadMtx(
#         mtx = file.path(outdir.star, 'matrix.mtx'),
#         cells = file.path(outdir.star, 'barcodes.tsv'),
#         features = file.path(outdir.star, 'features.tsv')
#     )
# ))
# 
# # Subset cells by RNA count, feature count, and MT content
# print("Filtering count matrix...")
# print(sobj <- cell_qc(sobj))
# 
# # Write passing barcodes to file
# outfile <- file.path(outdir, 'barcodes.pass.tsv')
# rownames(sobj@meta.data) %>%
#     write.table(file=outfile, quote=F, row.names=F, col.names=F)


# Run Azimuth
cat("Running Azimuth...\n")

### This was needed because of expired https cert
# options(url.method='libcurl')
# httr::set_config(httr::config(ssl_verifypeer = 0L))

# print(sobj.mapped <- Azimuth::RunAzimuth(
#     sobj,
#     reference = refname,
#     do.adt = do.adt
# ))


# Adding Azimuth results using Seurat::AddAzimuthResults requires
# an Azimuth mapping scores file, which is an RDS file containing a list
# with:
#   results$umap : Dimensionality reduction (UMAP) as a `DimReduc` object
#   results$pred.df : `data.frame` with predicted celltypes and scores
#   results$impADT :  `Assay` object with imputed values for antibody derived
#                      tags (optional)
# results <- list()
# 
# if('impADT' %in% Seurat::Assays(sobj.mapped)) {
#     results$impADT <- sobj.mapped[['impADT']]
# }
# 
# if('ref.umap' %in% Seurat::Reductions(sobj.mapped)) {
#     results$umap <- sobj.mapped[['ref.umap']]
# }
# 
# results$pred.df <- sobj.mapped@meta.data %>%
#     tibble::rownames_to_column('cell') %>%
#     dplyr::select(
#         cell,
#         dplyr::starts_with('predicted.'),
#         mapping.score
#     ) %>%
#     as.data.frame
# 
# az.outfile <- file.path(outdir, paste0('tmpAzimuthResults.rds'))
# saveRDS(results, file = az.outfile)
# 
# # Add Azimuth results
# sobj <- Seurat::AddAzimuthResults(sobj, filename=az.outfile)
# 
# # # Another way to get annotation levels
# # annotation.levels <- sobj@meta.data %>%
# #     dplyr::select(dplyr::starts_with('predicted.')) %>%
# #     dplyr::select(-dplyr::ends_with('.score')) %>%
# #     names %>%
# #     str_remove('predicted.')
# 
# # Predicted celltype, score, and mapping.score columns are added to the
# # metadata but are all NA (bug?). Add it manually here
# stopifnot(all(rownames(sobj@meta.data) == results$pred.df$cell))
# for (n in annotation.levels) {
#     pcol <- paste0('predicted.', n)
#     scol <- paste0('predicted.', n, '.score')
#     sobj@meta.data[[pcol]] <- results$pred.df[[pcol]]
#     sobj@meta.data[[scol]] <- results$pred.df[[scol]]
# }
# sobj@meta.data[['mapping.score']] <- results$pred.df[['mapping.score']]
# 
# # Update assays
# for (n in annotation.levels) {
#     aname <- paste0('prediction.score.', n)
#     sobj[[aname]] <- sobj.mapped@assays[[aname]]
# }
# 
# # Save Seurat object with mapping information
# saveRDS(sobj, file.path(outdir, 'tx_seurat.mapped.rds'))

# # Output celltype annotation files
# for (n in annotation.levels) {
#     outfile <- file.path(
#         outdir,
#         paste0(refname, '.', n, '.tsv')
#     )
#     sobj@meta.data %>% tibble::rownames_to_column('cell') %>%
#         dplyr::select(cell, paste0('predicted.', n)) %>%
#         write_tsv(file=outfile, quote='none', col_names = F)
# }
# 
# # Output filtered celltype annotation files
# for (n in annotation.levels) {
#     outfile <- file.path(
#         outdir,
#         paste0(refname, '.', n, '.filt.tsv')
#     )
#     sobj@meta.data %>%
#         dplyr::filter(.data[[paste0('predicted.', n, ".score")]] > thresh) %>%
#         tibble::rownames_to_column('cell') %>%
#         dplyr::select(cell, paste0('predicted.', n)) %>%
#         write_tsv(file=outfile, quote='none', col_names = F)
# }

cat("#--- annotate_azimuth.R completed. ---#\n")
