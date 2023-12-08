#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Import necessary libraries and modules
from typing import Optional # Import the Optional type hint for optional function arguments
from numbers import Number
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
from os import path
import pandas as pd
import sys


def scrublet_wrapper(
        counts_fn: str,
        features_fn: str,
        barcodes_fn: str,
        outdir: str,
        fcol: int = 0,
        bcol: int = 0,
        expected_doublet_rate: Optional[float] = None, 
        sim_doublet_ratio: Number = 2.0,
    ):
    """
    Wrapper function for Scrublet doublet detection.
    Parameters:
        counts_fn (str): File path to the input count matrix.
        features_fn (str): File path to the input features (genes) in tab-separated format.
        barcodes_fn (str): File path to the input barcodes in tab-separated format.
        outdir (str): Output directory where the results and plots will be saved.
        fcol (int, optional): Column index for features in the input file (default is 0).
        bcol (int, optional): Column index for barcodes in the input file (default is 0).
        expected_doublet_rate (float, optional): Expected doublet rate for the dataset (default is None).
        sim_doublet_ratio (Number, optional): Ratio for simulating doublets (default is 2.0).
    """
    
    # Read the counts data and transpose it
    counts_csc = scipy.io.mmread(counts_fn).T.tocsc()
    
    # Read features and barcodes data into Pandas Series
    with open(features_fn, 'r') as fh:
        features = pd.Series(l.strip('\n').split('\t')[fcol] for l in fh)

    with open(barcodes_fn, 'r') as fh:
        barcodes = pd.Series(l.strip('\n').split('\t')[bcol] for l in fh)

    # Check for incongruent input data dimensions
    assert counts_csc.shape == (len(barcodes), len(features)), 'Incongruent input'    
    
    # Calculate expected doublet rate if not provided
    if expected_doublet_rate is None:
        ## This is the calculation for 10x doublet rate but will be different for other platforms
        expected_doublet_rate = counts_csc.shape[0]/1000 * 0.008 

    print(f'Using expected doublet rate = {expected_doublet_rate}', file=sys.stderr)
    
    # Initialize Scrublet object with input data and parameters
    scrub = scr.Scrublet(
        counts_csc,
        expected_doublet_rate = expected_doublet_rate,
        sim_doublet_ratio = sim_doublet_ratio,
    )

    # Scrublet analysis to identify doublets
    _scores, _preds = scrub.scrub_doublets(
        min_counts = 3, # Minimum total UMI counts for a cell to be considered
        min_cells = 3, # Minimum number of cells a gene must be detected in to be considered
        min_gene_variability_pctl = 85,  # Minimum gene variability percentile for analysis
        n_prin_comps = 30 # Number of principal components for PCA analysis
    )

    ''' Plotting and saving '''
    # Create output directory if it doesn't exist
    os.makedirs(outdir, exist_ok = True)

    # Plot and save doublet score histogram
    scrub.plot_histogram()
    plt.savefig(path.join(outdir,'SCR.doublet_score_histogram.png'))

    print('Running UMAP...', file=sys.stderr)

    # Perform UMAP embedding and plot UMAP visualization
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.', file=sys.stderr)

    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig(os.path.join(outdir,'SCR.UMAP.png'))

    # Create results dataframe and save the predictions
    results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
    scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
    df = pd.DataFrame({
        'barcode': barcodes,
        'scrublet_score': pd.Series(scrub.doublet_scores_obs_),
        'pbool': pd.Series(scrub.predicted_doublets_)
    }).set_index('barcode')
    
    # Label scrublet predictions as 'Singlet' or 'Doublet' and save the results in different files
    df['scrublet_prediction'] = 'Singlet'
    df.loc[df['pbool'], 'scrublet_prediction'] = 'Doublet'
    
    # Save the results in separate output files
    with open(path.join(outdir, 'SCR.predictions_all.tsv'), 'w') as outh:
        df[['scrublet_prediction','scrublet_score']].to_csv(outh, sep='\t')
    
    with open(path.join(outdir, 'SCR.predictions_single.tsv'), 'w') as outh:
        df.loc[~df['pbool']][['scrublet_prediction','scrublet_score']].to_csv(outh, sep='\t')
    
    with open(path.join(outdir, 'SCR.barcodes.pass.tsv'), 'w') as outh:
        pd.Series(df.loc[~df['pbool']].index).to_csv(outh, index=False, header=False)
    
    ### Make summary of scrublet predictions and save it
    summary = pd.DataFrame(df.scrublet_prediction.value_counts())
    with open(path.join(outdir, 'SCR.summary.txt'), 'w') as outh:
        summary.to_csv(outh, sep = "\t")

# Check if the script is being run from the command line or as part of a Snakemake workflow
try:
    _ = snakemake.input
    _rate = None
    if 'doublet_rate' in snakemake.params.keys():
        _rate = float(snakemake.params['doublet_rate'])
    
    # Call the scrublet_wrapper function with input from Snakemake workflow
    scrublet_wrapper(
        counts_fn = snakemake.input['counts_mtx'],
        features_fn = snakemake.input['features_tsv'],
        barcodes_fn = snakemake.input['barcodes_tsv'],
        outdir = path.dirname(snakemake.output[0]),
        expected_doublet_rate = _rate
    )

# If Snakemake variables are not defined, run the function with default input data and parameters
except NameError:
    if __name__ == '__main__':
        input_dir = 'results/pbmc3p/celltype_annotation/pbmc3p.20k/cell_qc'
        scrublet_wrapper(
            counts_fn = path.join(input_dir, 'counts.mtx'),
            features_fn = path.join(input_dir, 'features.tsv'),
            barcodes_fn = path.join(input_dir, 'barcodes.tsv'),
            outdir = 'results/pbmc3p/celltype_annotation/pbmc3p.20k/scrublet',
            expected_doublet_rate = 0.08
        )
