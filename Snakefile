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

##########################################################################
# Determine list of genes with required information for each feature type
##########################################################################

# function to filter out entries with mostly N's
def filter_n(features):
    ratio_n = features.seq.str.count('N') / features.seq.apply(len)
    return features[ratio_n < 0.5]

ig = pd.read_csv(config['intergenic_stats'])
utr5 = pd.read_csv(config['5utr_stats'])
utr3 = pd.read_csv(config['3utr_stats'])

# upstream intergenic region
ig['gene'] = ig['left_gene'].where(ig['strand'] == 1, ig['right_gene'])
upstream_ids = filter_n(ig)['gene']

# downstream intergenic regions
ig['gene'] = ig['right_gene'].where(ig['strand'] == 1, ig['left_gene'])
downstream_ids = filter_n(ig)['gene']

# For UTRs, if predicted UTR length is very short, use the "static" assumed
# UTR length values instead (e.g. for T. cruzi, 75nt 5' UTR and 125nt 3' UTR)
short_5utrs = utr5.static_seq.apply(len) > utr5.seq.apply(len)
utr5.loc[short_5utrs, 'seq'] = utr5[short_5utrs].static_seq

short_3utrs = utr3.static_seq.apply(len) > utr3.seq.apply(len)
utr3.loc[short_3utrs, 'seq'] = utr3[short_3utrs].static_seq

utr5_ids = list(utr5.gene)
utr3_ids = list(utr3.gene)

with open(config['input_cds_fasta']) as fp:
    # parse gene ids from FASTA header entries, e.g.:
    # >TcCLB.410961.10:mRNA | organism=Trypanosoma_cruzi...
    fasta_gene_ids = []

    for line in fp.readlines():
        if line.startswith('>'):
            fasta_gene_ids.append(line.split(':')[0][1:])

# list of all genes
GENES = list(set(utr5_ids).intersection(utr3_ids)
                          .intersection(fasta_gene_ids)
                          .intersection(upstream_ids)
                          .intersection(downstream_ids))

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
If a pre-existing gene clustering is provided in the config file, it will be
used, otherwise hierarhical clustering will be performed.
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
        positive="build/sequences/{feature}/clusters/{cluster}.fa",
        negative="build/sequences/{feature}/clusters-negative/{cluster}.fa"
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
        expand("build/sequences/{feature}/genes/{gene}.fa", gene=GENES, feature=FEATURES)
    params:
        genes=GENES
    script:
        "scripts/get_gene_feature_sequences.py"

# rule foo:
#     input:
#         dynamic(expand("build/sequences/{feature}/genes/{{gene}}.fa",
#                        feature=FEATURES))

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
        jf_counts={output}.jf

        # count kmers and save output to plaintext file
        hash_size=`echo $((4 * {wildcards.kmer_size}))`
        jellyfish count -m {wildcards.kmer_size} -s ${{hash_size}} -t2 -o ${{jf_counts}} {input}
        jellyfish dump -c ${{jf_counts}} > {output}

        # delete intermediate file
        rm ${{jf_counts}}
        """

"""
Combines kmer counts for various features and size of k into a single CSV file.
"""
rule combine_kmer_counts:
    input:
        expand("build/kmers/{feature}/{kmer_size}/{gene}.txt",
                feature=FEATURES,
                kmer_size=['3', '4'],
                gene=GENES)
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
        feature_seqs="build/sequences/{feature}/clusters/{cluster}.fa",
        feature_neg_seqs="build/sequences/{feature}/clusters-negative/{cluster}.fa"
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

        outdir={config[output_dir]}/{config[version]}

        python2 {config[extreme_dir]}/src/GappedKmerSearch.py \
            -l {params.half_length} \
            -ming {params.ming} \
            -maxg {params.maxg} \
            -minsites {params.min_sites} \
            $outdir/{input.feature_seqs} \
            $outdir/{input.feature_neg_seqs} \
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
                    $outdir/{input.feature_seqs} \
                    $outdir/{input.feature_neg_seqs} \
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
        feature_seqs="build/sequences/{feature}/clusters/{cluster}.fa"
    output:
        "build/motifs/cmfinder/{feature}/{cluster}/{cluster}.finished"
    shell:
        """
        # create output directory and copy cluster feature sequences
        cd $(dirname {output})
        cp {config[output_dir]}/{config[version]}/{input.feature_seqs} .
        
        # echo "----DEBUGGING----"
        # pwd
        # cmd="cmfinder.pl -f {config[cmfinder_settings][motif_fraction]} \
        #             -m {config[cmfinder_settings][min_length]} \
        #             -b $(basename {input.feature_seqs})"
        # echo $cmd
        # echo "-----------------"

        # run cmfinder
        cmfinder.pl -f {config[cmfinder_settings][motif_fraction]} \
                    -m {config[cmfinder_settings][min_length]} \
                    -b $(basename {input.feature_seqs}) 2>/dev/null

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
        outdir={config[output_dir]}/{config[version]}
        outfile=$outdir/{output}

        cd $(dirname {input})

        # for each detected motif, iterate over all feature sequences and count
        # the occurences of that motif
        if ls *.cm 1>/dev/null 2>&1; then
            for motif in *.cm; do
                for feature_seqs in $outdir/build/sequences/*/clusters/*.fa; do
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

# DEBUGGING
# rule foo:
#     input:

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
        # cai=rules.compute_cai.output[0],
        cds=rules.generate_cds_stats.output[0],
        downstream_intergenic_region=rules.generate_intergenic_stats.output.downstream,
        upstream_intergenic_region=rules.generate_intergenic_stats.output.upstream,
        kmers=rules.combine_kmer_counts.output[0]
    output:
        "build/model/training_set.csv"
    script:
        "scripts/create_training_set.R"

"""
Trains Random Forest models for random subsets of clusters and features, and
determines variable importance for each model. Features which never have a
sufficiently high variable importance score are filtered out to generate a
smaller, more informative training set.
"""
rule filter_training_set:
    input:
        "build/model/training_set.csv"
    output:
        "build/model/training_set_filtered.csv"
    script:
        "scripts/create_filtered_training_set.R"

# vim: ft=python
