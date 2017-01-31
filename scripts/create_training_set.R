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

##/TEMP DEV####################################################################
feature_dir <- '/cbcb/nelsayed-scratch/keith/reg-predict/tcruzi/build/features'

inputs <- list(
    read.csv(file.path(feature_dir, 'motif_counts_extreme.csv')),
    read.csv(file.path(feature_dir, 'motif_counts_cmfinder.csv')),
    #read.csv(file.path(feature_dir, 'gene_features_cai.csv')),
    read.csv(file.path(feature_dir, '5utr_stats.csv')) %>%
        select(gene, utr5_len=length, utr5_gc=gc, utr5_ct=ct),
    read.csv(file.path(feature_dir, '3utr_stats.csv')) %>%
        select(gene, utr3_len=length, utr3_gc=gc, utr3_ct=ct),
    read.csv(file.path(feature_dir, 'polypyrimidine_tracts.csv')) %>%
        select(gene, sl_dist, polya_dist, tract_length, 
               polypyr_gc=gc, polypyr_ct=ct),
    read.csv(file.path(feature_dir, 'gene_stats_cds.csv')) %>% 
        select(gene, cds_length=length, cds_gc=gc, cds_ct=ct), 
    read.csv(file.path(feature_dir, 'gene_stats_upstream_intergenic_region.csv')) %>%
        select(gene, upstream_inter_cds_length=inter_cds_length, 
               upstream_intergenic_length=intergenic_length, upstream_gc=gc, upstream_ct=ct),
    read.csv(file.path(feature_dir, 'gene_stats_downstream_intergenic_region.csv')) %>%
        select(gene, downstream_inter_cds_length=inter_cds_length, 
               downstream_intergenic_length=intergenic_length, downstream_gc=gc, downstream_ct=ct)
)

###TEMP DEV\###################################################################

# including gene with missing data
#dat <- Reduce(function(...) merge(..., by='gene', all=TRUE), inputs)

dat <- Reduce(function(...) merge(..., by='gene'), inputs)
rownames(dat) <- dat$gene
dat <- dat %>% select(-gene)

# Use the caret findCorrelation function to find and remove highly correlated
# variables

# 5' UTR
#cor_mat <- cor(dat, method='spearman', use='pairwise.complete.obs')
cor_mat <- cor(dat, method='spearman')

# remove any entries with many NA's
mask <- apply(cor_mat, 1, function (x) { sum(is.na(x)) < 100 })
cor_mat <- cor_mat[mask,mask]


motifs_to_ignore <- findCorrelation(cor_mat, cutoff=0.5, verbose=FALSE)
dat <- dat[,-motifs_to_ignore]

print(sprintf("Ignoring %d/%d highly-correlated 5' UTR motifs...",
              length(motifs_to_ignore), nrow(cor_mat)))
print(sprintf("Remaining motifs: %d", ncol(dat)))

cor_mat <- cor(dat, method='spearman')
cor_mat[is.na(cor_mat)] <- 0

