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
MODULES = MODULES[:50]

# EXTREME runs
EXTREME_RUNS = config['extreme_settings'].keys()

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
        expand("build/motifs/extreme/{run}/{utr}/{module}/{module}.finished",
               utr=['5utr', '3utr'], run=EXTREME_RUNS, module=MODULES)
    output:
        "build/motifs/extreme/5utr-filtered.csv",
        "build/motifs/extreme/3utr-filtered.csv"
    shell:
        "touch {output}"
    # script:
    #     "scripts/combine_motifs.R"

"""
Compute Codon Adaptation Index for each CDS in the genome.
"""
rule compute_cai:
    output:
        "build/cai.csv"
    script:
        "scripts/compute_cai.py"

"""
Generates FASTA files for all 5' and 3'UTRs
"""
rule get_module_utr_sequences:
    input:
        module_assignments=config['module_assignments']
    output:
        "build/sequences/{utr}/{module}.fa",
        "build/sequences/{utr}/negative/{module}.fa"
    params:
        utr='{utr}',
        build_dir='build/sequences/{utr}'
    script:
        "scripts/get_module_utr_sequences.py"

"""
Performs RNA motif detection using EXTREME
http://www.ncbi.nlm.nih.gov/pubmed/24532725
"""
rule detect_utr_motifs_extreme:
    input:
        utr_seqs="build/sequences/{utr}/{module}.fa",
        utr_neg_seqs="build/sequences/{utr}/negative/{module}.fa"
    output:
        "build/motifs/extreme/{run}/{utr}/{module}/{module}.finished"
    params:
        half_length=lambda wildcards: config['extreme_settings'][wildcards.run]['half_length'],
        ming=lambda wildcards: config['extreme_settings'][wildcards.run]['min_gap_size'],
        maxg=lambda wildcards: config['extreme_settings'][wildcards.run]['max_gap_size'],
        min_sites=lambda wildcards: config['extreme_settings'][wildcards.run]['min_sites'],
        clustering_threshold=lambda wildcards: config['extreme_settings'][wildcards.run]['clustering_threshold'],
        max_motif_seeds=lambda wildcards: config['extreme_settings'][wildcards.run]['max_motif_seeds']
    shell:
        """
        cd $(dirname {output})

        echo {params.ming}

        python2 {config[extreme_dir]}/src/GappedKmerSearch.py \
            -l {params.half_length} \
            -ming {params.ming} \
            -maxg {params.maxg} \
            -minsites {params.min_sites} \
            {config[output_dir]}/{input.utr_seqs} \
            {config[output_dir]}/{input.utr_neg_seqs} \
            {wildcards.module}.words

        perl {config[extreme_dir]}/src/run_consensus_clusering_using_wm.pl \
            {wildcards.module}.words \
            {params.clustering_threshold}

        python2 {config[extreme_dir]}/src/Consensus2PWM.py \
            {wildcards.module}.words.cluster.aln \
            {wildcards.module}.wm
        
        for i in `seq 1 {params.max_motif_seeds}`; do
            if [[ -n $(grep ">cluster$i\s" {wildcards.module}.wm) ]]; then
                python2 {config[extreme_dir]}/src/EXTREME.py \
                    {config[output_dir]}/{input.utr_seqs} \
                    {config[output_dir]}/{input.utr_neg_seqs} \
                    --saveseqs \
                    {wildcards.module}.wm \
                    $i
            fi
        done

        touch {wildcards.module}.finished
        """

# vim: ft=python
