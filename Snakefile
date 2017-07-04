"""
Pipeline for the prediction of trypanosomatid regulatory elements
Keith Hughitt
September 02, 2016
"""
import os
import yaml
import pandas as pd

# specify build/output directory
workdir: os.path.join(config['output_dir'], config['version'])

# work-around for ambiguity
ruleorder: get_gene_feature_sequences > get_coexpression_cluster_feature_sequences 

# wildcard fills
EXTREME_RUNS = config['extreme_settings'].keys()
FEATURES = ['utr5', 'utr3', 'cds', 'downstream_intergenic_region', 'upstream_intergenic_region']

###############################
# Gene structure & composition
###############################
"""
Copies UTR and polypyrimidine tract statistics generated in separate pipelines
into the features dir for use here.
"""
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
Compute Codon Adaptation Index for each CDS in the genome.
"""
rule compute_cai:
    output:
        "build/features/gene_features_cai.csv"
    script:
        "scripts/compute_cai.py"

#############
# Clustering
#############
"""
Clusters genes based on their co-expression profiles
"""
rule generate_coexpression_clusters:
    output:
        dynamic("build/clusters/{cluster}.csv")
    script:
        "scripts/generate_coexpression_clusters.R"

############
# Sequences
############
"""
Generates cluster-level FASTA files for all 5' and 3' UTRs, intergenic regions, and CDS's
"""
rule get_coexpression_cluster_feature_sequences:
    input:
        utr5='build/features/5utr_stats.csv',
        utr3='build/features/3utr_stats.csv',
        cds="build/features/gene_stats_cds.csv",
        downstream_intergenic_region="build/features/gene_stats_downstream_intergenic_region.csv",
        upstream_intergenic_region="build/features/gene_stats_upstream_intergenic_region.csv",
        clusters="build/clusters/{cluster}.csv"
    output:
        positive="build/sequences/{feature}/{cluster}.fa",
        negative="build/sequences/{feature}/negative/{cluster}.fa"
    params:
        feature='{feature}',
        cluster='{cluster}'
    script:
        "scripts/get_cluster_feature_sequences.py"

"""
Generates gene-level FASTA files for all 5' and 3' UTRs, intergenic regions, and CDS's
"""
rule get_gene_feature_sequences:
    input:
        utr5='build/features/5utr_stats.csv',
        utr3='build/features/3utr_stats.csv',
        cds="build/features/gene_stats_cds.csv",
        downstream_intergenic_region="build/features/gene_stats_downstream_intergenic_region.csv",
        upstream_intergenic_region="build/features/gene_stats_upstream_intergenic_region.csv"
    output:
        dynamic("build/sequences/{feature}/genes/{gene}.fa")
    params:
        feature='{feature}'
    script:
        "scripts/get_gene_feature_sequences.py"

rule foo:
    input:
        dynamic(expand("build/sequences/{feature}/genes/{{gene}}.fa",
                       feature=FEATURES))

# """
# Splits up module-wide multifasta sequence files into separate FASTA files for
# each gene. This is useful for counting k-mers at the gene level using jellyfish
# """
# rule split_multifasta_files:
#     input:
#         "build/sequences/{feature}/{cluster}.fa"
#     output:
#         dynamic("build/sequences/{{feature}}/genes/{gene}.fa")
#     shell:
#         """
#         for input_file in {input}; do
#             for x in `grep '^>' ${{input_file}} | cut -c2-`; do
#                 echo $x
#                 awk -v seq="${{x}}" -v RS='>' '$1 == seq {{print RS $0}}' ${{input_file}} > {output}
#             done
#         done
#         """

##############
# Kmer counts
##############
"""
Uses Jellyfish to count kmers in different gene features (CDS, UTR's, etc.)
"""
rule count_kmers:
    input:
        "build/sequences/{feature}/genes/{gene}.fa"
    output:
        "build/kmers/{feature}/{kmer_size}/{gene}.txt"
    shell:
        """
        # indir=`dirname {input}`
        # outdir=`dirname {output}`

        # for gene_fasta in ${{indir}}/*.fa; do
        # intermediate and output filepaths
        # file_prefix=`basename ${{gene_fasta/.fa/}}`
        # jf_counts=${{outdir}}/${{file_prefix}}.jf
        jf_counts={output}.jf

        # count kmers and save output to plaintext file
        hash_size=`echo $((4 * {wildcards.kmer_size}))`
        jellyfish count -m {wildcards.kmer_size} -s ${{hash_size}} -t3 -o ${{jf_counts}} ${{gene_fasta}}
        jellyfish dump -c ${{jf_counts}} > {output}

        # delete intermediate file
        rm ${{jf_counts}}
        # done

        # touch {output}
        """

"""
Combines kmer counts for various features and size of k into a single CSV file.
"""
rule combine_kmer_counts:
    input:
        dynamic(expand("build/kmers/{feature}/{kmer_size}/{{gene}}.txt",
                       feature=FEATURES, 
                       kmer_size=['3', '4', '5']))
    output:
        "build/features/kmer_counts.csv"
    script:
        "scripts/combine_kmer_counts.py"

###########################
# Motif detection: EXTREME
###########################
"""
Scans co-expression clusters for presence of motifs detected using EXTREME
"""
rule count_motifs_extreme:
    input:
        dynamic(expand("build/motifs/extreme/{run}/{feature}/{{cluster}}/{{cluster}}.finished",
                        feature=FEATURES, run=EXTREME_RUNS))
    output:
        "build/features/motif_counts_extreme.csv"
    params:
        build_dir='build/motifs'
    script:
        "scripts/count_motifs_extreme.R"

"""
Performs RNA motif detection using EXTREME
http://www.ncbi.nlm.nih.gov/pubmed/24532725
"""
rule detect_feature_motifs_extreme:
    input:
        feature_seqs="build/sequences/{feature}/{cluster}.fa",
        feature_neg_seqs="build/sequences/{feature}/negative/{cluster}.fa"
    output:
        "build/motifs/extreme/{run}/{feature}/{cluster}/{cluster}.finished"
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
            {wildcards.cluster}.words

        perl {config[extreme_dir]}/src/run_consensus_clusering_using_wm.pl \
            {wildcards.cluster}.words \
            {params.clustering_threshold}

        python2 {config[extreme_dir]}/src/Consensus2PWM.py \
            {wildcards.cluster}.words.cluster.aln \
            {wildcards.cluster}.wm
        
        for i in `seq 1 {params.max_motif_seeds}`; do
            if [[ -n $(grep ">cluster$i\s" {wildcards.cluster}.wm) ]]; then
                python2 {config[extreme_dir]}/src/EXTREME.py \
                    {config[output_dir]}/{input.feature_seqs} \
                    {config[output_dir]}/{input.feature_neg_seqs} \
                    --saveseqs \
                    {wildcards.cluster}.wm \
                    $i
            fi
        done

        touch {wildcards.cluster}.finished
        """

###########################
# Motif detection: CMFinder
###########################

"""
Performs RNA motif detection using CMFinder
"""
rule detect_feature_motifs_cmfinder:
    input:
        feature_seqs="build/sequences/{feature}/{cluster}.fa"
    output:
        "build/motifs/cmfinder/{feature}/{cluster}/{cluster}.finished"
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

        touch {wildcards.cluster}.finished
        """

"""
Scans co-expression clusters for presence of motifs detected using CMFinder
"""
rule count_motifs_cmfinder:
    input:
        "build/motifs/cmfinder/{feature}/{cluster}/{cluster}.finished"
    output:
        "build/motifs/cmfinder/{feature}/{cluster}/{cluster}.cmsearch.gz"
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
Combines motif search output from cmsearch into a single table
"""
rule combine_cmsearch_results:
    input: 
        dynamic(expand("build/motifs/cmfinder/{feature}/{{cluster}}/{{cluster}}.cmsearch.gz",
                feature=FEATURES))
    output:
        "build/features/motif_counts_cmfinder.csv",
    script:
        "scripts/combine_cmsearch_results.py"

############################
# Training set construction
############################
"""
Combined multiple feature types and gene cluster assignments into a single
training set.
"""
rule create_training_set:
    input:
        extreme=rules.count_motifs_extreme.output[0],
        cmfinder=rules.combine_cmsearch_results.output[0],
        cai=rules.compute_cai.output[0],
        cds=rules.generate_cds_stats.output[0],
        downstream_intergenic_region=rules.generate_intergenic_stats.output.downstream,
        upstream_intergenic_region=rules.generate_intergenic_stats.output.upstream,
        kmers=rules.combine_kmer_counts.output[0]
    output:
        "build/model/training_set.csv"
    script:
        "scripts/create_training_set.R"

# vim: ft=python
