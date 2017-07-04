#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get_gene_feature_sequences.py
keith hughitt (khughitt@umd.edu)

Generates gene-level FASTA files of CDS, UTR, and intergenic region sequences
"""
import os
import pandas as pd

# load gene feature statistics and sequences
feature = snakemake.params['feature']

features_filepath = snakemake.input[feature]
features = pd.read_csv(features_filepath)

# For UTR's, if predicted UTR length is very short, use the "static" assumed
# UTR length values instead (e.g. for T. cruzi, 75nt 5' UTR and 125nt 3' UTR)
if snakemake.params['feature'] in ['utr5', 'utr3']:
    short_utrs = features.static_seq.apply(len) > features.seq.apply(len)

    features.loc[short_utrs, 'seq'] = features[short_utrs].static_seq

    # Not used here...
    # features.loc[short_utrs, 'gc']  = features[short_utrs].static_gc
    # features.loc[short_utrs, 'ct']  = features[short_utrs].static_ct

# remove sequences with mostly N's
ratio_n = features.seq.str.count('N') / features.seq.apply(len)
features = features[ratio_n < 0.5]


# save individual gene sequences to FASTA files
# with open(snakemake.output, 'w') as cluster_outfile:
for entry in features.itertuples():
    # output filepath
    gene_filepath = snakemake.output[0].replace('__snakemake_dynamic__', entry.gene)

    # write gene sequence to fasta file
    with open(gene_filepath, 'w') as gene_outfile:
        gene_outfile.write('>%s\n' % entry.gene)
        gene_outfile.write('%s' % entry.seq)

