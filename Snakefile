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

# exclude unassigned (grey) genes
df = df[df.color != 'grey']

# get a list of the co-expression modules
MODULES = list(df.color.unique())

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
Scans co-expression clusters for presence of motifs detected using EXTREME
"""
rule count_motifs_extreme:
    input:
        expand("build/motifs/extreme/{run}/{utr}/{module}/{module}.finished",
               utr=['5utr', '3utr'], run=EXTREME_RUNS, module=MODULES)
    output:
        "build/motifs/extreme-motif-counts.csv"
    script:
        "scripts/count_motifs_extreme.R"

"""
Combines motif search output from cmsearch into a single table
"""
rule combine_cmsearch_results:
    input: 
        expand("build/motifs/cmfinder/{utr}/{module}/{module}.cmsearch.gz",
               utr=['5utr', '3utr'], module=MODULES)
    output:
        "build/motifs/cmfinder-motif-counts.csv",
    script:
        "scripts/combine_cmsearch_results.py"

"""
Scans co-expression clusters for presence of motifs detected using CMFinder
"""
rule count_motifs_cmfinder:
    input:
        "build/motifs/cmfinder/{utr}/{module}/{module}.finished"
    output:
        "build/motifs/cmfinder/{utr}/{module}/{module}.cmsearch.gz"
    shell:
        """
        # full path to output file
        outfile={config[output_dir]}/{output}

        cd $(dirname {input})

        # for each detected motif, iterate over all UTR sequences and count the
        # occurences of that motif
        for motif in *.cm; do
            for utr_seqs in {config[output_dir]}/build/sequences/*/*.fa; do
                cmsearch --noali \
                         -E {config[cmfinder_settings][evalue_cutoff]} \
                         ${{motif}} \
                         ${{utr_seqs}} >> ${{outfile/.gz//}}
            done
        done

        # compress results
        gzip ${{outfile}}
        """

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
Performs RNA motif detection using CMFinder
"""
rule detect_utr_motifs_cmfinder:
    input:
        utr_seqs="build/sequences/{utr}/{module}.fa"
    output:
        "build/motifs/cmfinder/{utr}/{module}/{module}.finished"
    shell:
        """
        # create output directory and copy cluster utr sequences
        cd $(dirname {output})
        cp {config[output_dir]}/{input.utr_seqs} .

        # run cmfinder
        cmfinder.pl -f {config[cmfinder_settings][motif_fraction]} \
                    -m {config[cmfinder_settings][min_length]} \
                    -b $(basename {input.utr_seqs})

        # rebuilt and calibrate the motif covariance models using infernal
        for motif in *.motif.*; do
            cmbuild ${{motif}}.cm ${{motif}}
            cmcalibrate ${{motif}}.cm
        done

        rm latest.cm

        touch {wildcards.module}.finished
        """

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
