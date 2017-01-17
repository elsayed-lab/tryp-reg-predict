#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Converts CMFinder Stockholm motif alignment output to a PWM

CMFinder outputs motifs as stockholm-formatted MSA, including all of the input
sequences. Since it is often the case that some of the input sequences did not
contain the motif, and thus have a low weight assignment in the output, these
will be excluded for the purposes of computing the PWM.

Usage:

    cmfinder_motif_to_pwm.py /path/to/input.meme
"""
import os
import sys
import pandas as pd

# check for valid input filepath
infile = sys.argv[1]

if not os.path.isfile(infile):
    sys.exit("Invalid input filepath specified.")

outfile = infile.replace('h1.1', 'meme')

# weight threshold for determining which sequences to use
THRESHOLD = 0.5

# stockholm format columns
STOCKHOLM_GENE_IND = 1

# open input motif alignment and iterate over lines
fp = open(infile, 'rt')

# keep track of which input sequences contained the motif
motif_genes = []
motif_seqs = []

for line in fp:
    # parse sequence weights
    if ' WT ' in line:
        weight = float(line.split().pop())

        if weight >= THRESHOLD:
            gene_id = line.split()[STOCKHOLM_GENE_IND]
            motif_genes.append(gene_id)
    elif  line != '\n' and line != '//\n' and not line.startswith('#'):
        # parse relevant sequences
        gene_id, seq = line.split()

        if gene_id in motif_genes:
            motif_seqs.append(seq)

# convert to a pandas dataframe
df = pd.DataFrame([list(x) for x in motif_seqs])
pfm = df.apply(pd.Series.value_counts)

# drop gap row
pfm = pfm.iloc[1:,]

# assign uniform probabilities to gap regions
pfm.loc[:,pd.isnull(pfm).all()] = 0.25

# convert remaining nans to zeros
pfm = pfm.fillna(0)

# convert to PWM
pwm = pfm / pfm.shape[1]

# generate MEME version of output; since gaps aren't supported, we will
# simply replace them with uniform probability entries
with open(outfile, 'wt') as out:
    out.write("MEME version 4.9.0\n\n")
    out.write("ALPHABET= ACGT\n\n")
    out.write("strands: + -\n\n")
    out.write("Background letter frequencies (from uniform background):\n")
    out.write("A 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n")
    out.write("MOTIF M1 O1\n\n")
    out.write("letter-probability matrix: alength= 4 w= %d" % (pwm.shape[1]))

    # write pwm entries
    for row in pwm.T.iterrows():
        values = tuple(row[1].values.tolist())
        entry = "%0.3f %0.3f %0.3f %0.3f\n" % (values)
        out.write(entry)

    out.write('\n')

