Prediction of Trypanosomatid Regulatory Elements
================================================

Installation
------------

The Trypanosomatid Regulatory Elements prediction pipeline makes use a number
of different R and Python packages, as well as several standalone tools.

Below is a list of all of the requirements needed to run this pipeline.

### Requirements

#### Software requirements

- [CMFinder](http://bio.cs.washington.edu/CMfinderWeb/CMfinderDownload.pl) (0.2)
- [Infernal](http://eddylab.org/infernal/) (1.1.2+)

#### Python requirements

- [Python 3](https://www.python.org/downloads/)
- [Snakemake](https://snakemake.readthedocs.io/en/latest/) (3.10.0+)
- [Biopython](http://biopython.org/wiki/Biopython) (1.6.8+)
- [Pandas](http://pandas.pydata.org/) (0.19.2+)

#### R requirements

- [Bioconductor](https://bioconductor.org) (3.3+)
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (2.40.0+)
- [caret](http://topepo.github.io/caret/index.html) (6.0-73+)
- [seqLogo](https://www.bioconductor.org/packages/release/bioc/html/seqLogo.html) (1.38+)

Usage
-----

```sh
snakemake --configfile settings/config.yml combine_motifs
```


