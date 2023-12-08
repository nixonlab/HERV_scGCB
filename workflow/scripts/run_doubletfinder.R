#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyverse)
library(Seurat) 

# Check if DoubletFinder package is installed, if not, install it
if (!requireNamespace("DoubletFinder", quietly = TRUE))
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', quiet=TRUE)

library(DoubletFinder)

# Check if the script is run using Snakemake workflow, set input/output directories and parameters accordingly
if(exists("snakemake")) {
    seurat.rds <- snakemake@input[["seurat_rds"]]
    outdir <- dirname(snakemake@output[[1]])
    doublet.rate <- snakemake@params[["doublet_rate"]]
    cell.idents <- snakemake@params[['cell_idents']]
} else {
    cat("Not using snakemake\n")
    seurat.rds <- 'results/canonicalQC/Tonsil_1a/tx_seurat.mapped.rds'
    outdir <- 'results/canonicalQC/Tonsil_1a/doubletfinder'
    doublet.rate <- 0.08
    cell.idents <- 'predicted.celltype.l1'
}

# Read Seurat object from RDS file
sobj <- readRDS(seurat.rds)

# Preprocess data if SCT not applied
if(!('SCT' %in% Assays(sobj))) {
    sobj <- Seurat::SCTransform(
        sobj,
        vst.flavor="v2",
        variable.features.n = 5000,
        return.only.var.genes = FALSE
    )
}
sobj <- Seurat::RunPCA(sobj)

# Determine optimal pK values
#
# pK is the PC neighborhood size used to compute pANN, expressed as a 
# proportion of the merged real-artificial data. No default is set, as pK
# should be adjusted for each scRNA-seq dataset. Optimal pK values can be
# determined using mean-variance-normalized bimodality coefficient (BCmvn)
sweep.list <- DoubletFinder::paramSweep_v3(sobj, PCs = 1:50, sct=TRUE)
sweep.stats <- DoubletFinder::summarizeSweep(sweep.list)

pdf(file.path(outdir, 'DF.BCmvn_pK.pdf'), paper='USr')
bcmvn <- DoubletFinder::find.pK(sweep.stats)
dev.off()

emp.pK <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric), 'pK']
emp.pK <- as.numeric(as.character(emp.pK))

cat(sprintf(
    'Maximum BCmvn = %.3f found for pK = %.3f',
    max(bcmvn$BCmetric),
    emp.pK
))


# Determine optimal nExp values
#
# The total number of doublet predictions produced. This value can best be 
# estimated from cell loading densities into the 10X/Drop-Seq device, and 
# adjusted according to the estimated proportion of homotypic doublets.
stopifnot(cell.idents %in% names(sobj@meta.data))
homotypic.prop <- DoubletFinder::modelHomotypic(sobj@meta.data[[cell.idents]])
cat(sprintf(
    'Expected homotypic doublet proportion (identities="%s"): %.3f%%\n',
    cell.idents,
    homotypic.prop*100
))

nExp_poi <- round(doublet.rate * nrow(sobj@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

cat(sprintf(
    'Expecting %d doublets for multiplet rate of %.1f%%. \n',
    nExp_poi.adj,
    doublet.rate*100
))

# Run DoubletFinder with empirical pK and nExp calculated using the provided
# multiplet rate, adjusted for homotypic doublets
sobj.dd <- DoubletFinder::doubletFinder_v3(
    sobj,
    PCs = 1:50,
    pN = 0.25,
    pK = emp.pK,
    nExp = nExp_poi.adj,
    sct = TRUE
)
sobj.md <- sobj.dd@meta.data
score.col <- names(sobj.md)[grepl('^pANN', names(sobj.md))][1]
class.col <- names(sobj.md)[grepl('^DF.classifications', names(sobj.md))][1]

# Sweep nExp ± 2 percent
sweep.doublet.rate <- seq(doublet.rate-0.02, doublet.rate+0.02, by=0.005)
sweep.nExp_poi <- round(sweep.doublet.rate * nrow(sobj.md))
sweep.nExp_poi.adj <- round(sweep.nExp_poi * (1 - homotypic.prop))

for(cur.nExp in sweep.nExp_poi.adj) {
    sobj.dd <- DoubletFinder::doubletFinder_v3(
        sobj.dd,
        PCs = 1:50,
        pN = 0.25,
        pK = emp.pK,
        nExp = cur.nExp,
        sct = TRUE,
        reuse.pANN = score.col
    )
}

# Extract DoubletFinder predictions and scores
sobj.md <- sobj.dd@meta.data
sweep.cols <- names(sobj.md)[grepl('^DF.classifications', names(sobj.md))]
# remove the default
sweep.cols <- sweep.cols[sweep.cols != class.col]

DFpred <- sobj.md %>% 
    select(DF.classifications = class.col[1], score.col[1], all_of(sweep.cols))

cat(sprintf(
    'Renamed default column "%s" to "DF.classifications"\n',
    class.col
))

# Generate density plots for DoubletFinder classifications
plts <- lapply(c('DF.classifications', sweep.cols), function(cc) {
    ggplot(DFpred, aes(x=.data[[score.col]], color=.data[[cc]])) + 
        geom_density() + 
        theme_minimal() +
        ggtitle(cc) +
        NoLegend()
})

# Save the density plots as PDF
pdf(file.path(outdir, 'DF.score_density.pdf'), paper='USr')
plts[[1]]
(plts[[2]] + plts[[3]]) /  (plts[[4]] + plts[[5]])
(plts[[6]] + plts[[7]]) /  (plts[[8]] + plts[[9]])
dev.off()

# Generate summary table for DoubletFinder predictions
summary.df <- lapply(c('DF.classifications', sweep.cols),
                     function(col) table(DFpred[, col])
                     ) %>%
    bind_rows() %>%
    mutate(cname = c(class.col, sweep.cols)) %>%
    mutate(default = cname == class.col) %>%
    tidyr::separate(cname,
                    c('x', 'pN', 'pK', 'nExp'),
                    sep = '_',
                    remove = FALSE) %>%
    dplyr::mutate(across(pN:nExp, as.numeric)) %>%
    dplyr::select(-x) %>%
    dplyr::mutate(
        Doublet = as.numeric(Doublet),
        Singlet = as.numeric(Singlet)
    ) %>%
    dplyr::arrange(nExp)

summary.df[summary.df$default, 'cname'] <- 'DF.classifications'

print(summary.df)

# Write DoubletFinder summary to file
outfile <- file.path(outdir, 'DF.summary.tsv')
summary.df %>% readr::write_tsv(outfile, quote = 'none')

# Write DoubletFinder predictions to file (all)
outfile <- file.path(outdir, 'DF.predictions_all.tsv')
DFpred %>% 
    tibble::rownames_to_column('barcode') %>%
    readr::write_tsv(file=outfile, quote='none')

# Write DoubletFinder predictions to file (singlet only)
outfile <- file.path(outdir, 'DF.predictions_single.tsv')
DFpred %>%
    tibble::rownames_to_column('barcode') %>%    
    dplyr::filter(DF.classifications == 'Singlet') %>%
    readr::write_tsv(file=outfile, quote='none')

# Write passing barcodes to file
outfile <- file.path(outdir, 'DF.barcodes.pass.tsv')
DFpred %>%
    dplyr::filter(DF.classifications == 'Singlet') %>%
    rownames() %>%
    write.table(file=outfile, quote=F, row.names=F, col.names=F)

