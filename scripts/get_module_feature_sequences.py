#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get_module_feature_sequences.py
keith hughitt (khughitt@umd.edu)

Generates FASTA files containing the collection of CDS, UTR, and intergenic
region sequences, as well as a "negative" set of randomly selected sequences
from genes _not_ in the module.
"""
import os
import pandas as pd

# load gene/module mapping
module_assignments = pd.read_table(snakemake.config['module_assignments'])

# Exclude ungrouped 'grey' module genes
module_assignments = module_assignments[module_assignments.color != 'grey']

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

    features[short_utrs].seq = features[short_utrs].static_seq
    features[short_utrs].gc = features[short_utrs].static_gc
    features[short_utrs].ct = features[short_utrs].static_ct

# remove sequences with mostly N's
ratio_n = features.seq.str.count('N') / features.seq.apply(len)
features = features[ratio_n < 0.5]

# dataframe indices
GENE_IDX = 1
SEQ_IDX = 2

# iterate over co-expression modules and generate FASTA files containing
# sequences and negative (not in module) sequences for each gene in the module
for module in module_assignments.color.unique():
    # feature sequences for genes in module (positive)
    base_dir = os.path.join(snakemake.config['output_dir'],
                            snakemake.params['build_dir'])

    outfile = os.path.join(base_dir, module + ".fa")

    genes = module_assignments[module_assignments.color == module].gene_id
    seqs = features[features.gene.isin(genes)][['gene', 'seq']]

    # save positive sequences to FASTA file
    with open(outfile, 'w') as fp:
        for entry in seqs.itertuples():
            fp.write('>%s\n' % entry[GENE_IDX])
            fp.write('%s\n' % entry[SEQ_IDX])

    # feature sequences for genes not in module (negative)
    base_dir = os.path.join(base_dir, 'negative')
    outfile = os.path.join(base_dir, module + ".fa")

    # grab N random sequences NOT in the module of interest
    n = snakemake.config['extreme_num_neg_seqs']
    genes = module_assignments[module_assignments.color != module].sample(n).gene_id
    seqs = features[features.gene.isin(genes)][['gene', 'seq']]

    # save positive sequences to FASTA file
    with open(outfile, 'w') as fp:
        for entry in seqs.itertuples():
            fp.write('>%s\n' % entry[GENE_IDX])
            fp.write('%s\n' % entry[SEQ_IDX])

