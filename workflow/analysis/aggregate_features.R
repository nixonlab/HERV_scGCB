#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(SeuratDisk)
})

################################################################################
#   Utility Functions
################################################################################
source('workflow/analysis/01-helpers.R')

################################################################################
#   Parse arguments
################################################################################
if(exists("snakemake")) {
    infile <- snakemake@input[[1]]
    outfile <- snakemake@output[[1]]
} else {
    args = commandArgs(trailingOnly=TRUE)
    usage <- 'USAGE:\n  aggregate_features.R [-h] input.h5seurat [output.h5seurat]'
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
            gsub('.h5seurat$', '.agg.h5seurat', infile, ignore.case=T),
            gsub('$', '.agg.h5seurat', infile, ignore.case=T)
        )
    }
}

################################################################################
#   Main
################################################################################
hfile <- SeuratDisk::Connect(infile)
hfile$index()
sobj <- SeuratDisk::LoadH5Seurat(hfile,
                                 assays=list(RNA=c('counts')),
                                 reductions=FALSE,
                                 graphs=FALSE,
                                 neighbors=FALSE,
                                 images=FALSE
                                 )

#--- Add pcts assay
sobj[['pcts']] <- Seurat::CreateAssayObject(data = t(as.matrix(
    sobj@meta.data[ , c('percent.mt','percent.HERV','percent.L1','percent.TE')]
)))

#--- Add aggregated assays
mf <- sobj[['RNA']]@meta.features

#--- RNAclass
# Aggreated by Feature class: 
# c('CG', 'LTR', 'LINE') # CG is canonical gene
map.class <- data.frame(
    feature=rownames(mf),
    group=as.character(mf$feature.class)
)

table(map.class$group)

sobj[['RNAclass']] <- Seurat::CreateAssayObject(
    counts = aggregate_features(sobj, map.class)
)

stopifnot(sum(sobj[['RNAclass']]@counts) == sum(sobj[['RNA']]@counts))

mfrn <- rownames(sobj[['RNAclass']]@meta.features)
sobj[['RNAclass']]@meta.features$agg <- TRUE
sobj[['RNAclass']]@meta.features$nfeat <- table(map.class$group)[mfrn]
sobj[['RNAclass']]@meta.features$desc <- c('CG'='canonical', 'LTR'='HERV', 'LINE'='L1')[mfrn]

#--- RNAgtype
# Aggregated by gene biotype:
# c('protein-coding', 'lncRNA', 'sncRNA', 'pseudogene', 
#   'IG', 'TR', 'IG-pseudo', 'TR-pseudo', 'LINE', 'LTR')
map.gtype <- data.frame(
    feature=rownames(mf),
    group=ifelse(
        mf$feature.class=='CG',
        as.character(mf$feature.gtype), 
        as.character(mf$feature.class)
    )
)

# Recode biotypes so there are fewer
recode.gtype <- data.frame(
    row.names=sort(unique(map.gtype$group)),
    recode=sort(unique(map.gtype$group))
)
lncRNA.btypes <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding",  
                   "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA", 
                   "bidirectional_promoter_lncrna", "lncRNA")
sncRNA.btypes <- c("snRNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "miRNA", "ribozyme", "sRNA", 
                   "scaRNA", "vaultRNA", "vault_RNA", "scRNA")
recode.gtype['protein_coding', 'recode'] <- 'protein-coding'
recode.gtype[rownames(recode.gtype) %in% lncRNA.btypes, 'recode'] <- 'lncRNA'
recode.gtype[rownames(recode.gtype) %in% sncRNA.btypes, 'recode'] <- 'sncRNA'

recode.gtype[grepl('_pseudogene$', rownames(recode.gtype)), 'recode'] <- 'pseudogene'

recode.gtype[grepl('^IG.+_gene$', rownames(recode.gtype)), 'recode'] <- 'IG'
recode.gtype[grepl('^TR.+_gene$', rownames(recode.gtype)), 'recode'] <- 'TR'
recode.gtype[grepl('^IG.+_pseudogene$', rownames(recode.gtype)), 'recode'] <- 'IG-pseudo'
recode.gtype[grepl('^TR.+_pseudogene$', rownames(recode.gtype)), 'recode'] <- 'TR-pseudo'

map.gtype$recode <- recode.gtype[map.gtype$group, 'recode']
table(map.gtype$recode)

sobj[['RNAgtype']] <- Seurat::CreateAssayObject(
    counts=aggregate_features(sobj, map.gtype[,c('feature','recode')])
)

stopifnot(sum(sobj[['RNAgtype']]@counts) == sum(sobj[['RNA']]@counts))

mfrn <- rownames(sobj[['RNAgtype']]@meta.features)
nfeat.table <- table(map.gtype$recode)
sobj[['RNAgtype']]@meta.features$agg <- TRUE
sobj[['RNAgtype']]@meta.features$nfeat <- nfeat.table[mfrn]
sobj[['RNAgtype']]@meta.features$desc <- ""


#--- RNAfam
# HERV features are aggregated by family, others are left as-is
map.hervfam <- data.frame(
    feature=rownames(mf),
    group=ifelse(
        mf$feature.class %in% c('LTR','LINE'),
        as.character(mf$feature.gtype),
        rownames(mf)
    )
)

# ?? further aggregate LINE (now L1FLI, L1FLnI, L1)

# remove CG from map
map.hervfam <- map.hervfam[mf$feature.class %in% c('LTR','LINE'), ]
table(map.hervfam[ , 'group'])

sobj[['RNAfam']] <- Seurat::CreateAssayObject(
    counts=aggregate_features(sobj, map.hervfam)
)

stopifnot(sum(sobj[['RNAfam']]@counts) == sum(sobj[['RNA']]@counts))

mfrn <- rownames(sobj[['RNAfam']]@meta.features)
nfeat.table <- table(map.hervfam[ , 'group'])
for(n in names(sobj[['RNA']]@meta.features)) {
    sobj[['RNAfam']]@meta.features[,n] <- ifelse(
        mfrn %in% rownames(sobj[['RNA']]@meta.features), 
        as.character(sobj[['RNA']]@meta.features[mfrn, n]), 
        NA
    )
}
sobj[['RNAfam']]@meta.features$agg <- ifelse(mfrn %in% names(nfeat.table), 'agg', NA)
sobj[['RNAfam']]@meta.features$nfeat <- ifelse(mfrn %in% names(nfeat.table), nfeat.table[mfrn], NA)

sobj <- SeuratDisk::AppendData(infile, sobj)
DefaultAssay(sobj) <- 'RNA'
SeuratDisk::SaveH5Seurat(sobj, filename=outfile, overwrite=T)

cat("#--- aggregate_features.R completed. ---#\n")
