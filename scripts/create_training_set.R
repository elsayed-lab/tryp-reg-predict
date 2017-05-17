#!/usr/bin/env Rscript
###############################################################################
# 
# create_training_set.R
# keith hughitt (khughitt@umd.edu)
#
# Combines motif counts and other gene structure and composition variables into
# a single dataset and filters out highly correlated or low-information
# features.
#
###############################################################################
library('caret')
library('dplyr')
library('rtracklayer')

feature_dir <- '/cbcb/nelsayed-scratch/keith/reg-predict/tcruzi/build/features'

inputs <- list(
    '5utr_stats'=read.csv(file.path(feature_dir, '5utr_stats.csv')) %>%
        select(gene, utr5_len=length, utr5_gc=gc, utr5_ct=ct),
    '3utr_stats'=read.csv(file.path(feature_dir, '3utr_stats.csv')) %>%
        select(gene, utr3_len=length, utr3_gc=gc, utr3_ct=ct),
    'polypyr_tracts'=read.csv(file.path(feature_dir, 'polypyrimidine_tracts.csv')) %>%
        select(gene, sl_dist, polya_dist, tract_length, 
               polypyr_gc=gc, polypyr_ct=ct)
)

# Debugging 
#snakemake_input <- list(
#    extreme=file.path(feature_dir, 'motif_counts_extreme.csv'),
#    cmfinder=file.path(feature_dir, 'motif_counts_cmfinder.csv'),
#    cai=file.path(feature_dir, 'gene_features_cai.csv'),
#    cds=file.path(feature_dir, 'gene_stats_cds.csv'),
#    downstream_intergenic_region=file.path(feature_dir, 'gene_stats_upstream_intergenic_region.csv'),
#    upstream_intergenic_region=file.path(feature_dir, 'gene_stats_downstream_intergenic_region.csv')
#)

# retrieve snakemake input files using named keys; the snakemake@input list
# includes each values twice; once with a numeric key and once using a 
# string key, if specified.
input_keys <- c('extreme', 'cmfinder', 'cai', 'cds',
                'downstream_intergenic_region', 'upstream_intergenic_region')

snakemake_input <- list()
for (x in input_keys) {
    snakemake_input[[x]] <- snakemake@input[[x]]
}

# create a combined list of input feature files
inputs <- c(lapply(snakemake_input, read.csv), inputs)

# remove unndeeded columns
inputs[['cds']] <- inputs[['cds']] %>%
    select(gene, cds_length=length, cds_gc=gc, cds_ct=ct)
inputs[['upstream_intergenic_region']] <- inputs[['upstream_intergenic_region']] %>%
    select(gene, upstream_inter_cds_length=inter_cds_length,
           upstream_intergenic_length=intergenic_length, upstream_gc=gc,
           upstream_ct=ct)
inputs[['downstream_intergenic_region']] <- inputs[['downstream_intergenic_region']] %>%
    select(gene, downstream_inter_cds_length=inter_cds_length,
           downstream_intergenic_length=intergenic_length, downstream_gc=gc,
           downstream_ct=ct)

# load gene annotations
gff <- import.gff3('/cbcb/lab/nelsayed/ref_data/tcruzi_clbrener_esmeraldo-like/annotation/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like.gff')
genes <- gff[gff$type == 'gene']

# Combine input training variable sources, including genes with missing data
dat <- Reduce(function(...) merge(..., by='gene', all=TRUE), inputs)
rownames(dat) <- dat$gene
dat <- dat %>% select(-gene)

# Remove genes with a significant amount of missing feature data; often this
# is due to large number of N's on either side of the CDS
num_missing <- apply(dat, 1, function(x) { sum(is.na(x))})
gene_mask <- num_missing < (0.1 * ncol(dat))
print(sprintf("Removing %d/%d genes with many missing features", sum(!gene_mask), nrow(dat)))
dat <- dat[gene_mask,]

# QUESTION: do we only want to consider genes for which both UTR boundaries
# could be determined?..

# remove any features with many NA's; because UTR's could only be detected for
# some genes, there are many NA's for UTR length, etc., however, we still want
# to keep these fields.
#feature_mask <- apply(dat, 2, function (x) { sum(is.na(x)) < 100 })
#print(sprintf("Removing %d/%d features with many missing values", sum(!feature_mask), ncol(dat)))

# DEBUGGING
#print("Features being removed:")
#x <- apply(dat, 2, function(x) { sum(is.na(x)) })
#print(names(x[x > 0]))

#dat <- dat[,feature_mask]

# Load gene neighbor information and extend training set to include features
# of upstream and downstream genes
upstream <- follow(genes)
downstream <- precede(genes)

# add placeholder NA's to matrix representing upstream and downstream gene
# features
num_features <- ncol(dat)
dat <- cbind(dat, matrix(NA, nrow=nrow(dat), ncol=(2 * num_features)))

# determine column indices where upstream and downstream neighbor entries will go
upstream_cols <- (num_features + 1):(num_features * 2)
downstream_cols <- ((num_features * 2) + 1):(num_features * 3)

colnames(dat)[upstream_cols] <- paste0(colnames(dat)[1:num_features], '_upstream_neighbor')
colnames(dat)[downstream_cols] <- paste0(colnames(dat)[1:num_features], '_downstream_neighbor')

# iterate over genes and for each one that exists in the training set, copy
# the feature columns from its up- and down-stream neighbors, if they exist.
for (i in 1:length(genes)) {
    gid <- genes$ID[i]

    # skip genes that are not in the training set inputs
    if (!gid %in% rownames(dat)) {
        next
    }

    # add features from upstream genes, if they exist
    if (!is.na(upstream[i])) {
        upstream_gid <- genes$ID[upstream[i]]

        if (upstream_gid %in% rownames(dat)) {
            upstream_features <- dat[rownames(dat) == upstream_gid, 1:num_features]
            dat[rownames(dat) == gid, upstream_cols] <- upstream_features 
        }

    }

    # add features from downstream genes, if they exist
    if (!is.na(downstream[i])) {
        downstream_gid <- genes$ID[downstream[i]]

        if (downstream_gid %in% rownames(dat)) {
            downstream_features <- dat[rownames(dat) == downstream_gid, 1:num_features]
            dat[rownames(dat) == gid, downstream_cols] <- downstream_features 
        }

    }
}

# generate a feature correlation matrix
print("Generating feature correlation matrix")
cor_mat <- cor(dat, method='spearman', use='pairwise.complete.obs')

# Replace the few remaining NA's (much less than 1% of total values in testing)
# with 0's to allow the caret findCorrelation function to be used
cor_mat[is.na(cor_mat)] <- 0

# Use the caret findCorrelation function to find and remove highly correlated
# variables
redundant_features <- findCorrelation(cor_mat, cutoff=0.5, verbose=FALSE)
print(sprintf("Removing %d/%d features which have high redundancy",
              length(redundant_features), ncol(cor_mat)))

dat <- dat[,-redundant_features]
print(sprintf("Remaining features: %d", ncol(dat)))

write.csv(dat, file=snakemake@output[[1]])

#cor_mat <- cor(dat, method='spearman')
#cor_mat[is.na(cor_mat)] <- 0
