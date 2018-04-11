Prediction of Trypanosomatid Regulatory Elements
================================================

Overview
--------

### Features

- Sequence motifs
    - [X] 5' UTR
    - [X] 3' UTR
    - [X] Upstream gene's 3' UTR
    - [X] Downstream genes 5' UTR
    - [X] Upstream intergenic region
    - [X] Downstream intergenic region
    - [X] CDS
- Sequence composition
    - [X] 5' UTR GC/CT composition
    - [X] 3' UTR GC/CT composition
    - [X] CDS GC/CT composition
    - [X] Polypyrimidine tract GC/CT composition
    - [X] Kmer counts
- Sequence lengths
    - [X] 5' UTR length
    - [X] 3' UTR length
    - [X] Polypyrimidine tract length
    - [X] Interenic region / inter-CDS length
- Other
    - [X] CDS codon adaptation index (CAI)

Installation
------------

The Trypanosomatid Regulatory Elements prediction pipeline makes use a number
of different R and Python packages, as well as several standalone tools.

Below is a list of all of the requirements needed to run this pipeline.

### Requirements

#### Software requirements

- [CMFinder](http://bio.cs.washington.edu/CMfinderWeb/CMfinderDownload.pl) (0.2)
- [Infernal](http://eddylab.org/infernal/) (1.1.2+)
- [EXTREME](https://github.com/uci-cbcl/EXTREME) (2.0)
- [Jellyfish]() ()

#### Python requirements

- [Python 3](https://www.python.org/downloads/)
- [Snakemake](https://snakemake.readthedocs.io/en/latest/) (3.10.0+)
- [Biopython](http://biopython.org/wiki/Biopython) (1.6.8+)
- [Numpy]() ()
- [Pandas](http://pandas.pydata.org/) (0.19.2+)
- [PyYAML](http://pyyaml.org/)(3.12+)

#### R requirements

- [Bioconductor](https://bioconductor.org) (3.3+)
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (2.40.0+)
- [caret](http://topepo.github.io/caret/index.html) (6.0-73+)
- [tibble](https://cran.r-project.org/web/packages/tibble/index.html) (1.2+)
- [seqLogo](https://www.bioconductor.org/packages/release/bioc/html/seqLogo.html) (1.38+)

Note: In order to avoid running out of memory during execution, the
hierarchical clustering portion of the EXTREME script
`run_consensus_clusering_using_wm.pl` may need to be edited to increase the
value `Xmx`, e.g.: `-Xmx10000m`.

Usage
-----

TODO: describe software for predicting UTR boundaries, etc.

```sh
snakemake --configfile settings/config.yml combine_motifs
```


