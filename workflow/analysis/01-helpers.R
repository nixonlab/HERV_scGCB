#! /usr/bin/env Rscript

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

cell_summary_data <- function(sobj) {
    # Feature metadata
    fmeta <- sobj[['RNA']]@meta.features

    # Select mitochondrial features and calculate mt percentage
    mt_feats <- grepl('^MT-', fmeta$name)
    sobj[['percent.mt']] <-
        Seurat::PercentageFeatureSet(sobj, features = rownames(fmeta)[mt_feats])

    # Select HERV features and calculate HERV percentage
    herv_feats <- !is.na(fmeta$feature.class) & fmeta$feature.class == 'LTR'
    sobj[['percent.HERV']] <-
        Seurat::PercentageFeatureSet(sobj, features = rownames(fmeta)[herv_feats])

    # Select L1 features and calculate L1 percentage
    l1_feats <- !is.na(fmeta$feature.class) & fmeta$feature.class == 'LINE'
    sobj[['percent.L1']] <-
        Seurat::PercentageFeatureSet(sobj, features = rownames(fmeta)[l1_feats])

    # Select TE features and calculate TE percentage
    te_feats <- herv_feats | l1_feats
    sobj[['percent.TE']] <-
        Seurat::PercentageFeatureSet(sobj, features = rownames(fmeta)[te_feats])
    sobj
}


#' Bind a list of vectors into a sparse column-oriented matrix
#'
#' @param lst
#'
#' @return
#' @export
#'
#' @examples
sv.cbind <- function (lst) {
    require(Matrix)
    input <- lapply( lst, as, "dsparseVector" )
    thelength <- unique(sapply(input,length))
    stopifnot( length(thelength)==1 )

    Matrix::sparseMatrix(
        x=unlist(lapply(input,slot,"x")),
        i=unlist(lapply(input,slot,"i")),
        p=c(0,cumsum(sapply(input,function(x){length(x@x)}))),
        dims=c(thelength,length(input))
    )
}


#' Add counts from multiple features to an aggregate feature
#'
#' @param sobj
#'
#' @return
#' @export
#'
#' @examples
aggregate_features <- function(sobj, group_map, assay=sobj@active.assay) {
    require(Matrix)
    # Get the counts
    countmat <- sobj[[assay]]@counts

    # Remove features in group_map that do not appear in countmat
    group_map <- group_map[group_map[,1] %in% rownames(countmat),]

    # Separate counts that are to be aggregated from those to leave alone
    countmat.ret <- countmat[!(rownames(countmat) %in% group_map[,1]), ]
    countmat.sel <- countmat[(rownames(countmat) %in% group_map[,1]), ]

    # Reorder group_map to match countmat.sel
    rownames(group_map) <- group_map[,1]
    group_map <- group_map[rownames(countmat.sel), ]
    stopifnot(all(rownames(countmat.sel) == group_map[,1]))

    groups <- unique(group_map[,2])
    agg <- lapply(groups,
                  function(grp) {
                      Matrix::colSums(countmat.sel[group_map[,2]==grp, ], sparseResult=T)
                  }) %>%
        sv.cbind() %>%
        Matrix::t()

    agg@Dimnames[[1]] <- groups
    agg@Dimnames[[2]] <- colnames(countmat.sel)

    # ?? problem if countmat.ret has no rows?
    rbind2(countmat.ret, agg)
}

#' Determine mapping from HERV locus to family
#'
#' HERV names are extracted from the Seurat object and the HERV family is
#' assumed to be the part of the locus name preceding `sep`.
#'
#' @param sobj Seurat object with HERVs
#' @param sep Character separating family name from locus identifier. This is
#'    an underscore in the Telescope annotation and is replaced with a hyphen
#'    by Seurat. Default is '-'.
#' @param exclude_prefix Gene name prefixes to exclude from the mapping.
#'
#' @return A dataframe with HERV loci in column 1 and the corresponding HERV
#'    family in column 2.
#' @export
#'
#' @examples
get_herv_family_mapping <- function(sobj, sep = '-', exclude_prefix = c('ENSG', 'L1FLnI', 'L1FLI', 'L1ORF2')) {
    exclude_re <- paste0('^(', paste(exclude_prefix, collapse = '|'), ')')
    ret <- data.frame(locus = rownames(sobj)) %>%
        dplyr::filter(!grepl(exclude_re, locus)) %>%
        tidyr::separate(locus, c('family'), sep = sep, remove = F, extra = 'drop')
    ret
}


#' Assign celltypes using celltype prediction score thresholds
#'
#' @param scores Matrix of celltype prediction scores, rows are celltype labels
#'    and columns are cells. Found in the data slot of Azimuth prediction assay,
#'    i.e. `seuratobj[['prediction.score.celltype.l1']]@data`
#' @param conf.thresh Threshold for prediction score, which ranges from 0 to 1.
#' @param missing Label to use for cells that do not meet threshold.
#'    Default is "_low_
#'
#' @return A vector with the celltype assignment for each cell if threshold is
#'     met.
#' @export
#'
#' @examples
#' l1.conf <- assign_celltypes_conf(sobj[['prediction.score.celltype.l1']]@data)
#'
assign_celltypes_conf <- function(scores, conf.thresh=0.95, missing='_low') {
    tmp <- do.call(rbind, lapply(1:ncol(scores), function(col) {
        idx <- which(scores[,col] > conf.thresh)
        if(length(idx)==1) {
            return(c(rownames(scores)[idx], scores[idx, col]))
        } else {
            return(c('_low', 0.))
        }
    })) %>% data.frame
    names(tmp) <- c('celltype','score')
    tmp$celltype
}


#--- Plotting
umap_guides <- function(trunc_upper=unit(1, "cm")) {
    axis <- ggh4x::guide_axis_truncated(
        trunc_lower = unit(0, "npc"),
        trunc_upper = trunc_upper
    )
    guides(x=axis, y=axis)
}

umap_guide_theme <- function(...) {
    theme(
        axis.line = element_line(arrow = arrow(length = unit(0.25, "cm"))),
        axis.title = element_text(hjust = 0, size=6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
    )
}




