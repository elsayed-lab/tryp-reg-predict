#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get_gene_feature_sequences.py
keith hughitt (khughitt@umd.edu)

Generates gene-level FASTA files of CDS, UTR, and intergenic region sequences
"""
import os
import random
import pandas as pd

# load gene feature statistics and sequences
feature = snakemake.params['feature']

if feature in ['5utr', '3utr']:
    # 5' and 3'UTR stats
    features_filepath = snakemake.config[feature + '_stats']
else:
    # CDS, intergenic regions
    features_filepath = snakemake.input[feature]

features = pd.read_csv(features_filepath)

# For UTR's, if predicted UTR length is very short, use the "static" assumed
# UTR length values instead (e.g. for T. cruzi, 75nt 5' UTR and 125nt 3' UTR)
if snakemake.params['feature'] in ['5utr', '3utr']:
    short_utrs = features.static_seq.apply(len) > features.seq.apply(len)

    features.loc[short_utrs, 'seq'] = features[short_utrs].static_seq
    features.loc[short_utrs, 'gc']  = features[short_utrs].static_gc
    features.loc[short_utrs, 'ct']  = features[short_utrs].static_ct

# remove sequences with mostly N's
ratio_n = features.seq.str.count('N') / features.seq.apply(len)
features = features[ratio_n < 0.5]

# dataframe indices
GENE_IDX = 1
SEQ_IDX = 2

# save individual gene sequences to FASTA files
with open(snakemake.output.positive, 'w') as cluster_outfile:
    for entry in features.itertuples():
        gene_filepath = snakemake.output.gene.replace('__snakemake_dynamic__', entry[GENE_IDX])

        # write gene sequence to fasta file
        with open(gene_filepath, 'w') as gene_outfile:
            gene_outfile.write('>%s\n' % entry[GENE_IDX])
            gene_outfile.write('%s\n' % entry[SEQ_IDX])

