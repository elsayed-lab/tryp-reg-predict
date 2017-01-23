#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
compute_cai.py
keith hughitt (khughitt@umd.edu)

Snakemake script for computing the codon adaptation index (CAI) for a set of
sequences.
"""
from Bio import SeqIO

# load input CDS FASTA sequenc
infile = snakemake.config['input_cds_fasta']

# open file to store results in
fp = open(snakemake.output[0], 'w')
fp.write('gene,length,gc,ct,seq')

# parse CDS FASTA entries
for seq in SeqIO.parse(infile, format="fasta"):
    # compute CDS statistics
    seq_id = seq.id.split(':')[0]

    # Get GC- and CT-richness
    gc_richness = round((seq.seq.count('G') + seq.seq.count('C')) / len(seq), 3)
    ct_richness = round((seq.seq.count('C') + seq.seq.count('T')) / len(seq), 3)

    row = '\n{},{},{},{},{}'.format(seq_id, len(seq), gc_richness, ct_richness, str(seq.seq))
    fp.write(row)

fp.close()
