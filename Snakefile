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

# TESTING
MODULES = MODULES[:5]

# EXTREME runs
EXTREME_RUNS = config['extreme_settings'].keys()

"""
Scans co-expression clusters for presence of motifs detected using EXTREME
"""
rule count_motifs_extreme:
    input:
        expand("build/motifs/extreme/{run}/{feature}/{module}/{module}.finished",
               feature=['5utr', '3utr', 'cds', 'downstream_intergenic_region',
                        'upstream_intergenic_region'], 
               run=EXTREME_RUNS, module=MODULES)
    output:
        "build/features/motif_counts_extreme.csv"
    params:
        build_dir='build/motifs'
    script:
        "scripts/count_motifs_extreme.R"

"""
Combines motif search output from cmsearch into a single table
"""
rule combine_cmsearch_results:
    input: 
        expand("build/motifs/cmfinder/{feature}/{module}/{module}.cmsearch.gz",
               feature=['5utr', '3utr', 'cds', 'downstream_intergenic_region',
                        'upstream_intergenic_region'], 
               module=MODULES)
    output:
        "build/features/motif_counts_cmfinder.csv",
    script:
        "scripts/combine_cmsearch_results.py"

rule copy_input_features:
    input:
        config['polypyrimidine_stats'],
        config['5utr_stats'],
        config['3utr_stats']
    output:
        'build/features/polypyrimidine_tracts.csv',
        'build/features/5utr_stats.csv',
        'build/features/3utr_stats.csv'
    shell:
        """
        for x in {input}; do
            outfile=$(basename $x)
            echo "cp $x build/features/$outfile"
            cp $x build/features/$outfile
        done
        """

"""
Scans co-expression clusters for presence of motifs detected using CMFinder
"""
rule count_motifs_cmfinder:
    input:
        "build/motifs/cmfinder/{feature}/{module}/{module}.finished"
    output:
        "build/motifs/cmfinder/{feature}/{module}/{module}.cmsearch.gz"
    shell:
        """
        # full path to output file
        outfile={config[output_dir]}/{output}

        cd $(dirname {input})

        # for each detected motif, iterate over all feature sequences and count
        # the occurences of that motif
        if ls *.cm 1>/dev/null 2>&1; then
            for motif in *.cm; do
                for feature_seqs in {config[output_dir]}/build/sequences/*/*.fa; do
                    cmsearch --noali \
                            -E {config[cmfinder_settings][evalue_cutoff]} \
                            ${{motif}} \
                            ${{feature_seqs}} >> ${{outfile/.gz/}}
                done
            done

            # compress results
            gzip ${{outfile/.gz/}}
        else
            # if no motifs were detected, create a dummy output file to satisfy
            # snakemake
            touch ${{outfile}}
        fi

        """

"""
Compute Codon Adaptation Index for each CDS in the genome.
"""
rule compute_cai:
    output:
        "build/features/gene_features_cai.csv"
    script:
        "scripts/compute_cai.py"

"""
Generates table containing sequences and basic stats for CDS's
"""
rule generate_cds_stats:
    output:
        "build/features/gene_stats_cds.csv"
    script:
        "scripts/generate_cds_stats.py"

"""
Compute upstream and downstream intergenic region statistics for each gene
"""
rule generate_intergenic_stats:
    output:
        downstream="build/features/gene_stats_downstream_intergenic_region.csv",
        upstream="build/features/gene_stats_upstream_intergenic_region.csv"
    run:
        # load intergenic stats
        df = pd.read_csv(config['intergenic_stats'])

        # upstream intergenic regions
        df['gene'] = df['left_gene'].where(df['strand'] == 1, df['right_gene'])
        upstream = df[['gene', 'inter_cds_length', 'intergenic_length', 'gc',
                       'ct', 'seq']]
        upstream.to_csv(output.upstream, index=False)

        # downstream intergenic regions
        df['gene'] = df['right_gene'].where(df['strand'] == 1, df['left_gene'])
        downstream = df[['gene', 'inter_cds_length', 'intergenic_length',
                         'gc', 'ct', 'seq']]
        downstream.to_csv(output.downstream, index=False)

"""
Determine upstream and downstream neighbors for each gene.

NOTE: no longer needed: R used directly instead.
"""
rule determine_neighboring_genes:
    output:
        'build/features/gene_neighbors.csv'
    run:
        # load intergenic stats
        df = pd.read_csv(config['intergenic_stats'])

        # upstream genes for each gene
        upstream_df = pd.concat({
            'gene':          df['right_gene'].where(df['strand'] == 1, df['left_gene']),
            'upstream_gene': df['left_gene'].where(df['strand'] == 1, df['right_gene'])
        }, axis=1)

        # downstream genes for each gene
        downstream_df = pd.concat({
            'gene':          df['left_gene'].where(df['strand'] == 1, df['right_gene']),
            'downstream_gene': df['right_gene'].where(df['strand'] == 1, df['left_gene'])
        }, axis=1)

        result = pd.merge(upstream_df, downstream_df, on='gene', how='outer')
        result.to_csv(output[0])

"""
Generates FASTA files for all 5' and 3' UTRs, intergenic regions, and CDS's
"""
rule get_coexpression_module_feature_sequences:
    input:
        cds=rules.generate_cds_stats.output[0],
        downstream_intergenic_region=rules.generate_intergenic_stats.output.downstream,
        upstream_intergenic_region=rules.generate_intergenic_stats.output.upstream
    output:
        "build/sequences/{feature}/{module}.fa",
        "build/sequences/{feature}/negative/{module}.fa"
    params:
        feature='{feature}',
        build_dir='build/sequences/{feature}'
    script:
        "scripts/get_module_feature_sequences.py"

"""
Splits up module-wide multifasta sequence files into separate FASTA files for
each gene. This is useful for counting k-mers at the gene level using jellyfish
"""
rule split_multifasta_files:
    input:
        # rules.get_coexpression_module_feature_sequences.output[0]
        expand("build/sequences/{feature}/{module}.fa",
               feature=['5utr', '3utr', 'cds', 
                       'downstream_intergenic_region',
                       'upstream_intergenic_region'], 
               module=MODULES)
    output:
        dynamic("build/sequences/genes/{gene}.fa")
    shell:
        """
        for input_file in {input}; do
            for x in `grep '^>' ${{input_file}} | cut -c2-`; do
                outfile="build/sequences/genes/${{x}}.fa"
                awk -v seq="${{x}}" -v RS='>' '$1 == seq {{print RS $0}}' ${{input_file}} > ${{outfile}}
            done
        done
        """

rule count_kmers:
    input:
        # "build/sequences/genes/{gene}.fa"
        rules.split_multifasta_files.output
    output:
        "build/kmers/{gene}.txt"
    shell:
        """
        jellyfish count -m4 -s256 -t4 -o {output}.jf {input}
        jellyfish dump -c {output}.jf > {output}
        """

rule test_kmer_counts:
    input: 
        dynamic("build/kmers/{gene}.txt")

"""
Performs RNA motif detection using CMFinder
"""
rule detect_feature_motifs_cmfinder:
    input:
        feature_seqs="build/sequences/{feature}/{module}.fa"
    output:
        "build/motifs/cmfinder/{feature}/{module}/{module}.finished"
    shell:
        """
        # create output directory and copy cluster feature sequences
        cd $(dirname {output})
        cp {config[output_dir]}/{input.feature_seqs} .

        # run cmfinder
        cmfinder.pl -f {config[cmfinder_settings][motif_fraction]} \
                    -m {config[cmfinder_settings][min_length]} \
                    -b $(basename {input.feature_seqs})

        # rebuilt and calibrate the motif covariance models using infernal
        if ls *.motif.* 1>/dev/null 2>&1; then
            for motif in *.motif.*; do
                cmbuild ${{motif}}.cm ${{motif}}
                cmcalibrate ${{motif}}.cm
            done

            rm -f latest.cm
        fi

        touch {wildcards.module}.finished
        """

"""
Performs RNA motif detection using EXTREME
http://www.ncbi.nlm.nih.gov/pubmed/24532725
"""
rule detect_feature_motifs_extreme:
    input:
        feature_seqs="build/sequences/{feature}/{module}.fa",
        feature_neg_seqs="build/sequences/{feature}/negative/{module}.fa"
    output:
        "build/motifs/extreme/{run}/{feature}/{module}/{module}.finished"
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
            {config[output_dir]}/{input.feature_seqs} \
            {config[output_dir]}/{input.feature_neg_seqs} \
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
                    {config[output_dir]}/{input.feature_seqs} \
                    {config[output_dir]}/{input.feature_neg_seqs} \
                    --saveseqs \
                    {wildcards.module}.wm \
                    $i
            fi
        done

        touch {wildcards.module}.finished
        """

#
# snakemake directives
#
rule create_training_set:
    input:
        extreme=rules.count_motifs_extreme.output[0],
        cmfinder=rules.combine_cmsearch_results.output[0],
        cai=rules.compute_cai.output[0],
        cds=rules.generate_cds_stats.output[0],
        downstream_intergenic_region=rules.generate_intergenic_stats.output.downstream,
        upstream_intergenic_region=rules.generate_intergenic_stats.output.upstream
    output:
        "build/model/training_set.csv"
    script:
        "scripts/create_training_set.R"

# vim: ft=python
