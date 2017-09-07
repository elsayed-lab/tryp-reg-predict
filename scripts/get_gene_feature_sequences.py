#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get_gene_feature_sequences.py
keith hughitt (khughitt@umd.edu)

Generates gene-level FASTA files of CDS, UTR, and intergenic region sequences
"""
import os
import pandas as pd

for i, features_filepath in enumerate(snakemake.input):
    # Determine feature type
    infile = os.path.basename(features_filepath)

    if infile == '5utr_stats.csv':
        ftype = '5utr'
        out_prefix = 'build/sequences/utr5/genes'
    elif infile == '3utr_stats.csv':
        ftype = '3utr'
        out_prefix = 'build/sequences/utr3/genes'
    elif infile == 'gene_stats_cds.csv':
        ftype = 'cds'
        out_prefix = 'build/sequences/cds/genes'
    elif infile == 'gene_stats_downstream_intergenic_region.csv':
        ftype = 'downstream_intergenic'
        out_prefix = 'build/sequences/downstream_intergenic_region/genes'
    elif infile == 'gene_stats_upstream_intergenic_region.csv':
        ftype = 'upstream_intergenic'
        out_prefix = 'build/sequences/upstream_intergenic_region/genes'

    # load gene feature statistics and sequences
    features = pd.read_csv(features_filepath)

    # For UTR's, if predicted UTR length is very short, use the "static" assumed
    # UTR length values instead (e.g. for T. cruzi, 75nt 5' UTR and 125nt 3' UTR)
    if ftype in ['5utr', '3utr']:
        short_utrs = features.static_seq.apply(len) > features.seq.apply(len)
        features.loc[short_utrs, 'seq'] = features[short_utrs].static_seq

    # exclude genes with incomplete information
    features = features[features.gene.isin(snakemake.params['genes'])]

    # save individual gene sequences to FASTA files
    for entry in list(features.itertuples()):
        # skip genes without complete UTR information
        if entry.gene not in snakemake.params['genes']:
            continue
            
        # output filepath
        outfile = os.path.join(out_prefix, '%s.fa' % entry.gene)

        # write gene sequence to fasta file
        with open(outfile, 'w') as gene_outfile:
            gene_outfile.write('>%s\n' % entry.gene)
            gene_outfile.write('%s' % entry.seq)

