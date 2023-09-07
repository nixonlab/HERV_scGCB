#! /usr/bin/env python

# Definition of dataframes that will be used in other parts of the workflow
#
#import pandas as pd

# read in the SRA run table CSV file and store in dataframe
# use rename method on the dataframe to change the column names
sra_runtable = pd.read_csv(
    config['sra_runtable']['filename']
    ).rename(columns={
        config['sra_runtable']['sampid']: 'sampid',
	config['sra_runtable']['runid']: 'runid',
	config['sra_runtable']['layout']: 'layout',
        config['sra_runtable']['assay']: 'assay'
    })

# group the rows of the sra_runtable dataframe by sampid and apply the list function to runid column
# (because there could be more than one run by sample)
# the list function takes an iterable (here a column of the df) and returns a list containing all the elements of the iterable
# 
# this returns a new dataframe where for each sample, the runid column contains a list of all the SRA run IDs associated with that sample
samples = pd.DataFrame(sra_runtable.groupby('sampid')['runid'].apply(list))

# set runid as the index of sra_runtable and store the resulting dataframe in runs 
# (runs is a new dataframe that has the same rows as sra_runtable but with runid column as the index)
# this means that we could use the .loc method to access rows/cols in the dataframe by index
runs = sra_runtable.set_index('runid')

# get samples grouped by assay
assays = pd.DataFrame(sra_runtable.groupby('assay')['sampid'].apply(list))

