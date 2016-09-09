#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Snakemake script for computing the codon adaptation index (CAI) for a set of
sequences.
"""
from Bio import SeqIO
from Bio.SeqUtils import CodonUsage
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex

# load input CDS FASTA sequenc
infile = snakemake.config['input_cds_fasta']

# create CAI index
cai = CodonAdaptationIndex()
cai.generate_index(infile)

# open file to store results in
fp = open('build/cai.csv', 'w')
fp.write('gene_id,cai')

# compute CAI for each CDS
for seq in SeqIO.parse(infile, format="fasta"):
   gene_cai = cai.cai_for_gene(str(seq.seq))
   fp.write('\n%s,%f' % (seq.id, gene_cai))

fp.close()
