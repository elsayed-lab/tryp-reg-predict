#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
combine_kmer_counts.py
keith hughitt (khughitt@umd.edu)

Combines output from k-mer counts across different features into a single
matrix.
"""
import os
import glob
import pandas as pd

"""Main script body"""
# dataframe to store result in
result = pd.DataFrame()

# get a list of all genes
genes = set([os.path.splitext(os.path.basename(x))[0] for x in snakemake.input])

# iterate over features, one gene at a time
for i, gene_id in enumerate(genes):
    # data.frame to store results for a single gene
    dat = pd.DataFrame()

    # get all kmer input files for specified gene
    infiles = [x for x in snakemake.input if gene_id + ".txt" in x]

    for infile in infiles:
        # skip empty files; can happen with sequences containing all N's
        if os.stat(infile).st_size == 0:
            continue

        # feature type
        FEATURE_IDX = 2
        feature = infile.split('/')[FEATURE_IDX]

        # load k-mers as a single column dataframe
        kmers = pd.read_csv(infile, sep=' ', index_col=0, header=None)

        # update row and column names
        kmers.columns = [gene_id]
        kmers.index = feature + '_' + kmers.index

        # add rows to gene-specific data frame
        dat = dat.append(kmers)

    # add column to result table
    result = result.join(dat, how='outer')

# transpose matrix so that rows=genes and cols=features
result = result.T

# convert gene id index to a column
feature_names = sorted(result.columns)
result['gene'] = result.index
result = result.reindex_axis(['gene'] +  feature_names, axis=1)

# save result to a file
result.fillna(0).to_csv(snakemake.output[0], index=False)

