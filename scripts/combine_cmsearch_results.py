#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
combine_cmsearch_results.py
keith hughitt (khughitt@umd.edu)

Combines output from multiple cmsearch runs results into a single gene x motif
count matrix.
"""
import gzip
import pandas as pd

def main():
    """Main script execution"""
    # iterate over cmsearch results and generate a corresponding count matrix
    result = None

    for filepath in snakemake.input:
        # parse cmsearch results for motifs from a single module
        df = parse_cmsearch_results(filepath)

        # update results matrix
        if result is None:
            result = df
        else:
            # merge counts into result matrix
            result = pd.concat([result, df], axis=1)

    # save results
    result.fillna(0).to_csv(snakemake.output[0])

def parse_cmsearch_results(filepath):
    """Parses file containing one or more cmsearch run results"""
    # open cmsearch results
    fp = gzip.open(filepath, 'rt')

    # column containing gene indentifiers in cmsearch table output
    GENE_ID_INDEX = 5

    # use a dictionary of dictionaries to keep track of motif/gene counts
    motifs = {}

    # iterate over lines and parse relevant entries
    for line in fp:
        # motif
        if line.startswith('# query'):
            motif = "%s:%s" % (filepath.replace('build/motifs/', ''),
                               line.split().pop())

            if motif not in motifs:
                motifs[motif] = {}
        # gene
        elif line.startswith('  ('):
            gene = line.split()[GENE_ID_INDEX]

            # update counts
            if gene in motifs[motif]:
                motifs[motif][gene] += 1
            else:
                motifs[motif][gene] = 1

    return pd.DataFrame.from_dict(motifs).fillna(0)

if __name__ == "__main__":
    main()

