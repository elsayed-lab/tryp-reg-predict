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

"""
Combine motifs, filtering out redundant or low information ones.
"""
rule combine_motifs:
    input:
        utr5=expand("build/motifs/extreme/5utr/{module}.words", module=MODULES)
        # utr3=expand("build/motifs/extreme/3utr/{module}.words", module=MODULES)
    output:
        utr5="build/motifs/extreme/5utr-filtered.csv"
        # utr3="build/motifs/extreme/3utr-filtered.csv"
    shell:
        "touch {output}"

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
        "build/sequences/5utr/{module}.fa",
        "build/sequences/5utr/negative/{module}.fa"
    params:
        build_dir='build/sequences/5utr'
    script:
        "scripts/get_module_utr_sequences.py"

rule get_module_3utr_sequences:
    input:
        module_assignments=config['module_assignments'],
        utr_stats=config['3utr_stats']
    output:
        "build/sequences/3utr/{module}.fa",
        "build/sequences/3utr/negative/{module}.fa"
    params:
        build_dir='build/sequences/3utr'
    script:
        "scripts/get_module_utr_sequences.py"

rule detect_5utr_motifs_extreme:
    input:
        utr5_seqs="build/sequences/5utr/{module}.fa",
        utr5_neg_seqs="build/sequences/5utr/negative/{module}.fa"
    output:
        "build/motifs/extreme/5utr/{module}.words"
    params:
        build_dir="build/motifs/extreme/5utr",
        word_file="build/motifs/extreme/5utr/{module}.words",
        weight_matrix="build/motifs/extreme/5utr/{module}.wm"
    shell:
        """
        python2 {config[extreme_dir]}/src/GappedKmerSearch.py \
            -l {config[extreme_half_length]} \
            -ming {config[extreme_min_gap_size]} \
            -maxg {config[extreme_max_gap_size]} \
            -minsites {config[extreme_min_sites]} \
            {input.utr5_seqs} \
            {input.utr5_neg_seqs} \
            {params.word_file}

        perl {config[extreme_dir]}/src/run_consensus_clusering_using_wm.pl \
            {params.word_file} {config[extreme_clustering_threshold]}

        python2 {config[extreme_dir]}/src/Consensus2PWM.py \
            {params.word_file}.cluster.aln \
            {params.weight_matrix}

        for i in `seq 1 {config[extreme_max_motif_seeds]}`; do
            if [[ -n $(grep ">cluster$i\s" {params.weight_matrix}) ]]; then
                python2 {config[extreme_dir]}/src/EXTREME.py \
                    {input.utr5_seqs} \
                    {input.utr5_neg_seqs} \
                    --saveseqs \
                    {params.weight_matrix} \
                    $i
                fi
        done
        """

# vim: ft=python
