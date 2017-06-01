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

def main():
    """Main script body"""
    # dataframe to store result in
    result = pd.DataFrame()

    # TESTING
    dat = pd.DataFrame()

    # iterate over features
    for filepath in snakemake.input:
        # input_dir = os.path.dirname(filepath)

        # # feature name and k-mer size
        # base_dir, kmer_size = os.path.split(input_dir)
        # feature = os.path.split(base_dir)[-1]
        
        # # iterate over genes
        # input_files = glob.glob(os.path.join(input_dir, '*.txt'))

        # intermediate dataframe
        # dat = pd.DataFrame()

        # print("Parsing kmer counts from %s" % input_dir)

        # for infile in input_files:
            # skip empty files

        # testing
        infile = filepath

        if os.stat(infile).st_size == 0:
            continue

        # gene id
        gene_id = os.path.splitext(os.path.basename(infile))[0]

        # load k-mers as a single column dataframe
        kmers = pd.read_csv(infile, sep=' ', index_col=0, header=None)

        # update row and column names
        kmers.columns = [gene_id]
        kmers.index = feature + '_' + kmers.index

        # add column to result table
        dat = dat.join(kmers, how='outer')

        # append to result dataframe
        # result = result.append(dat)

    # TESTING
    result = dat

    # transpose matrix so that rows=genes and cols=features
    result = result.T

    # convert gene id index to a column
    feature_names = sorted(result.columns)
    result['gene'] = result.index
    result = result.reindex_axis(['gene'] +  feature_names, axis=1)

    # save result to a file
    result.fillna(0).to_csv(snakemake.output[0], index=False)

if __name__ == '__main__':
    main()
