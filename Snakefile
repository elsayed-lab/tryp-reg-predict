"""
Pipeline for the prediction of trypanosomatid regulatory elements
Keith Hughitt
September 02, 2016
"""
import os
import yaml
import pandas as pd

# specify build/output directory
workdir: config['output_dir']

# get co-expression module names
df = pd.read_table(config['module_assignments'])

# get a list of the co-expression modules
MODULES = list(df.color.unique())

# TEMP: limit to first 3 modules during testing...
MODULES = MODULES[:3]

#
# snakemake directives
#
rule create_training_set:
    input:
        utr5="build/motifs/extreme/5utr-filtered.csv",
        utr3="build/motifs/extreme/3utr-filtered.csv"
    output:
        "build/model/training_set.csv"
    shell:
        "touch {output}"
    # script:
    #     "scripts/create_training_set.py"

rule filter_and_combine_motifs:
    input:
        utr5=expand("build/motifs/extreme/5utr/{module}.csv", module=MODULES),
        utr3=expand("build/motifs/extreme/3utr/{module}.csv", module=MODULES)
    output:
        utr5="build/motifs/extreme/5utr-filtered.csv",
        utr3="build/motifs/extreme/3utr-filtered.csv"

"""
Compute Codon Adaptation Index for each CDS in the genome.
"""
rule compute_cai:
    output:
        "build/cai.csv"
    script:
        "scripts/compute_cai.py"

rule get_module_5utr_sequences:
    input:
        module_assignments=config['module_assignments'],
        utr_stats=config['5utr_stats']
    output:
        expand("build/sequences/5utr/{module}.fa", module=MODULES),
        expand("build/sequences/5utr/negative/{module}.fa", module=MODULES)
    params:
        build_dir='build/sequences/5utr'
    script:
        "scripts/get_module_utr_sequences.py"

rule get_module_3utr_sequences:
    input:
        module_assignments=config['module_assignments'],
        utr_stats=config['3utr_stats']
    output:
        expand("build/sequences/3utr/{module}.fa", module=MODULES),
        expand("build/sequences/3utr/negative/{module}.fa", module=MODULES)
    params:
        build_dir='build/sequences/3utr'
    script:
        "scripts/get_module_utr_sequences.py"

rule detect_5utr_motifs_extreme:
    input:
        utr5_seqs=expand("build/sequences/5utr/{module}.fa", module=MODULES),
        utr5_neg_seqs=expand("build/sequences/5utr/negative/{module}.fa", module=MODULES)
    output:
        expand("build/motifs/extreme/5utr/{module}.csv", module=MODULES)
    params:
        build_dir='build/motifs/extreme/5utr'
    shell:
        """
        python2 {config.extreme_dir}/src/GappedKmerSearch.py -l 8 -ming 0 -maxg 10 -minsites 5 GM12878_NRSF_ChIP.fasta GM12878_NRSF_ChIP_shuffled.fasta GM12878_NRSF_ChIP.words
        perl {config.extreme_dir}/src/run_consensus_clusering_using_wm.pl GM12878_NRSF_ChIP.words 0.3
        python2 {config.extreme_dir}/src/Consensus2PWM.py GM12878_NRSF_ChIP.words.cluster.aln GM12878_NRSF_ChIP.wm
        """

rule detect_3utr_motifs_extreme:
    input:
        utr3_seqs=expand("build/sequences/3utr/{module}.fa", module=MODULES),
        utr3_neg_seqs=expand("build/sequences/3utr/negative/{module}.fa", module=MODULES)
    output:
        expand("build/motifs/extreme/3utr/{module}.csv", module=MODULES)
    params:
        build_dir='build/motifs/extreme/3utr'
    shell:
        "touch {output}"

# vim: ft=python
