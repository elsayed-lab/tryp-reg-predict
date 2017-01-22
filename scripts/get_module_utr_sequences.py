#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generates FASTA files containing the collection of UTR sequences for which the
UTR lengths are known, as well as a "negative" set of randomly selected
sequences from genes _not_ in the module.
"""
import os
import pandas as pd

# load gene/module mapping
module_assignments = pd.read_table(snakemake.input['module_assignments'])

# load UTR statistics and sequences
utr_stats_filepath = snakemake.config[snakemake.params['utr'] + '_stats']
utr_stats = pd.read_csv(utr_stats_filepath)

# If predicted UTR length is very short, use the "static" assumed UTR length
# values instead (e.g. for T. cruzi, 75nt 5' UTR and 125nt 3' UTR)
short_utrs = utr_stats.static_seq.apply(len) > utr_stats.seq.apply(len)

utr_stats[short_utrs].seq = utr_stats[short_utrs].static_seq
utr_stats[short_utrs].gc = utr_stats[short_utrs].static_gc
utr_stats[short_utrs].ct = utr_stats[short_utrs].static_ct

# dataframe indices
GENE_IDX = 1
SEQ_IDX = 2

# iterate over modules and generate FASTA files containing module UTR
# sequences and negative (not in module) UTR sequences
for module in module_assignments.color.unique():
    # utr sequences for genes in module (positive)
    base_dir = os.path.join(snakemake.config['output_dir'],
                            snakemake.params['build_dir'])

    outfile = os.path.join(base_dir, module + ".fa")

    genes = module_assignments[module_assignments.color == module].gene_id
    seqs = utr_stats[utr_stats.gene.isin(genes)][['gene', 'seq']]

    # save positive sequences to FASTA file
    with open(outfile, 'w') as fp:
        for entry in seqs.itertuples():
            fp.write('>%s\n' % entry[GENE_IDX])
            fp.write('%s\n' % entry[SEQ_IDX])

    # utr sequences for genes not in module (negative)
    base_dir = os.path.join(base_dir, 'negative')
    outfile = os.path.join(base_dir, module + ".fa")

    # grab N random sequences NOT in the module of interest
    n = snakemake.config['extreme_num_neg_seqs']
    genes = module_assignments[module_assignments.color != module].sample(n).gene_id
    seqs = utr_stats[utr_stats.gene.isin(genes)][['gene', 'seq']]

    # save positive sequences to FASTA file
    with open(outfile, 'w') as fp:
        for entry in seqs.itertuples():
            fp.write('>%s\n' % entry[GENE_IDX])
            fp.write('%s\n' % entry[SEQ_IDX])

