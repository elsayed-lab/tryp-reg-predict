#!/usr/bin/env Rscript
###############################################################################
# 
# create_filtered_training_set.R 
# keith hughitt (khughitt@umd.edu)
#
# Trains Random Forest models for random subsets of clusters and features, and
# determines variable importance for each model. Features which never have a
# sufficiently high variable importance score are filtered out to generate a
# smaller, more informative training set.
#
###############################################################################
suppressMessages(library("randomForest"))
suppressMessages(library("dplyr"))
suppressMessages(library("tibble"))

set.seed(1)

train <- read.csv(snakemake@input[['unfiltered']], row.names=1)

# output files for feature importance and model error rate stats
feature_importance_file <- sub('training_set_filtered.csv',
                               'feature_importance.csv', snakemake@output[[1]])

model_error_rate_file <- sub('training_set_filtered.csv',
                             'model_class_error_rates.csv', snakemake@output[[1]])

# number of itereations
NUM_ITERATIONS <- 50000

# number of clusters to use when training a random forest
NUM_CLUSTERS <- 2

# ratio of features to genes to randomly select for each model, prior to
# filtering out mostly zero-valued columns; after zero- and NA-filtering, 
# the ratio will be much smaller (e.g. from 35x ~> 2-6x)
FEATURE_TO_GENE_RATIO <- 35

# Number of high-performing features to preserve
NUM_FEATURES_TO_KEEP <- 500

#
# generate_train_subset
#
# Creates a subset including only genes from k randomly selected clusters
#
generate_train_subset <- function(train, clusters, k=2, balanced=TRUE,
                                  min_nonzero_values=5, max_ratio_na=0.25) {
    # randomly select k clusters
    ind <- clusters$cluster %in% sample(unique(clusters$cluster), k)
    clusters_subset <- clusters[ind,]

    # balance training set if requested
    if (balanced) {
        # balance based on smallest cluster size
        n = min(table(clusters_subset$cluster))

        clusters_subset <- clusters_subset %>% 
            group_by(cluster) %>% 
            sample_n(n) %>%
            ungroup
        
        gene_subset <- as.character(clusters_subset$gene )
        response <- as.factor(clusters_subset$cluster)
    } else {
        gene_subset <- as.character(clusters$gene[ind])
        response <- as.factor(clusters$cluster[ind])
    }

    train_subset <- train[row.names(train) %in% gene_subset,]

    # discard features with only a few non-zero entries
    num_nonzeros <- apply(train_subset, 2, function(x) { sum(x != 0, na.rm=TRUE) })
    feature_mask <- num_nonzeros >= min_nonzero_values
    print(sprintf("Removing %d/%d features with little information relating to selected clusters.",
                sum(!feature_mask), length(feature_mask)))
    train_subset <- train_subset[,feature_mask]

    # remove features that are missing for more than 25% of genes; note that if
    # this cutoff is set too low, some useful features such as 5' and 3' UTR length
    # will be excluded.
    feature_cutoff <- max_ratio_na * nrow(train_subset)
    feature_mask <- apply(train_subset, 2, function(x) { sum (is.na(x)) } ) <= feature_cutoff
    print(sprintf("Removing %d/%d features with many NA values.",
                sum(!feature_mask), length(feature_mask)))
    train_subset <- train_subset[,feature_mask]

    dat <- cbind(train_subset, response)

    return(dat)
}

###############################################################################
# Feature selection
#
# Approach: Randomly select a subset of k clusters of genes and n features, and
# construct a random forest model to predict cluster membership. Save the top x
# features, as measured by variable importance and repeat the process some
# number of times. After sufficient iterations have been performed, filter out
# all features which never showed a high variable importance.
###############################################################################

# create lists to keep track of maximum and total variable importance for each
# feature, as well as the number of times that feature was randomly selected
total_var_importance <- setNames(vector("list", ncol(train)), colnames(train))
total_var_importance[1:length(total_var_importance)] <- 0

times_included <- total_var_importance
max_var_importance <- total_var_importance

# Also keep track of clusters tested and average class error for each model
clusters_tested <- c()
oob_error_rates <- c()

# TEMP
#clusters <- read.csv('/cbcb/nelsayed-scratch/keith/reg-predict/tcruzi/3.0/build/clusters/clusters.csv')
#clusters <- read.csv(snakemake@input[['clusters']])
clusters <- read.csv(file.path('build', 'features', 'coex_clusters.csv'))
message('WD:')
message(getwd())
clusters <- clusters[clusters$gene %in% rownames(train),]

# Exclude clusters containing only a single element (can happen due to
# post-clustering filtering during training set construction)
clusters_to_remove <- names(table(clusters$cluster)[table(clusters$cluster) < 2])
clusters <- clusters[!clusters$cluster %in% clusters_to_remove,]

# TESTING 2017/09/20
# For now, let's also exclude very small clusters to avoid biasing our feature
# selection
clusters_to_remove <- names(table(clusters$cluster)[table(clusters$cluster) <= 5])
clusters <- clusters[!clusters$cluster %in% clusters_to_remove,]

for (i in 1:NUM_ITERATIONS) {
    message(sprintf("Feature selection round %d/%d", i, NUM_ITERATIONS))

    # randomly select k clusters
    ind <- clusters$cluster %in% sample(unique(clusters$cluster), NUM_CLUSTERS)
    gene_subset <- as.character(clusters$gene[ind])
    response <- as.factor(clusters$cluster[ind])

    train_subset <- train[row.names(train) %in% gene_subset,]

    # randomly select features
    num_features <- FEATURE_TO_GENE_RATIO * nrow(train_subset)
    train_subset <- train_subset[,sample(ncol(train_subset), min(num_features, ncol(train_subset)))]

    # discard features with only a few non-zero entries
    num_nonzeros <- apply(train_subset, 2, function(x) { sum(x != 0, na.rm=TRUE) })
    feature_mask <- num_nonzeros >= 5
    train_subset <- train_subset[,feature_mask]

    # remove features that are missing for more than 15% of genes
    feature_cutoff <- 0.15 * nrow(train_subset)
    feature_mask <- apply(train_subset, 2, function(x) { sum (is.na(x)) } ) <= feature_cutoff
    train_subset <- train_subset[,feature_mask]

    dat <- cbind(train_subset, response)

    # If NA's still remain, impute values before training
    if (sum(is.na(dat)) > 0) {
        dat <- rfImpute(response ~ ., dat)
    }

    # train a random forest model
    rf <- randomForest(response ~ ., data=dat)

    #top_features <- names(head(rf$importance[order(rf$importance, decreasing=TRUE),], 3))

    clusters_tested <- c(clusters_tested, paste0(unique(response), collapse='_'))
    oob_error_rates <- c(oob_error_rates, mean(rf$confusion[,'class.error']))
    
    # update variable importance counters
    for (x  in rownames(rf$importance)) {
        total_var_importance[[x]] <- total_var_importance[[x]] + rf$importance[x,1]
        times_included[[x]] <- times_included[[x]] + 1

        if(rf$importance[x,1] > max_var_importance[[x]]) {
            max_var_importance[[x]] <- rf$importance[x,1]
        }
    }
}

# save feature selection results
features <- data.frame(total_variable_importance=unlist(total_var_importance),
                       max_variable_importance=unlist(max_var_importance),
                       times_included=unlist(times_included)) %>%
    rownames_to_column('feature') %>%
    mutate(average_variable_importance=total_variable_importance / times_included) %>%
    arrange(desc(average_variable_importance))
write.csv(features, file=feature_importance_file, quote=FALSE, row.names=FALSE)

# save model performance results
error_rates <- data.frame(clusters=clusters_tested, error_rate=oob_error_rates)
write.csv(error_rates, file=model_error_rate_file, quote=FALSE, row.names=FALSE)

# keep only the top N features with the highest average variable importance
# during feature selection
top_features <- features$feature[1:NUM_FEATURES_TO_KEEP]
train <- train[,colnames(train) %in% c('response', top_features)]

# Save results
write.csv(train, file=snakemake@output[[1]], quote=FALSE)

