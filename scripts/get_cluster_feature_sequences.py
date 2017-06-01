#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get_cluster_feature_sequences.py
keith hughitt (khughitt@umd.edu)

Generates FASTA files containing the collection of CDS, UTR, and intergenic
region sequences, as well as a "negative" set of randomly selected sequences
from genes _not_ in the cluster.
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

# load cluster memship csv
cluster = pd.read_csv(snakemake.input.clusters)
cluster_genes = list(cluster.gene)

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

# feature sequences for genes in module (positive)
seqs = features[features.gene.isin(cluster_genes)][['gene', 'seq']]

# save cluster sequences to FASTA files
with open(snakemake.output.positive, 'w') as cluster_outfile:
    for entry in seqs.itertuples():
        # append to cluster multi fasta
        cluster_outfile.write('>%s\n' % entry[GENE_IDX])
        cluster_outfile.write('%s\n' % entry[SEQ_IDX])

# feature sequences for genes not in module (negative)

# grab N random sequences NOT in the module of interest
noncluster_genes = [x for x in list(features.gene) if x not in cluster_genes]

n = snakemake.config['extreme_num_neg_seqs']
negative_genes = random.sample(noncluster_genes, n)

seqs = features[features.gene.isin(negative_genes)][['gene', 'seq']]

# save positive sequences to FASTA file
with open(snakemake.output.negative, 'w') as fp:
    for entry in seqs.itertuples():
        fp.write('>%s\n' % entry[GENE_IDX])
        fp.write('%s\n' % entry[SEQ_IDX])

