#!/usr/bin/env Rscript
###############################################################################
# 
# generate_coexpression_clusters.R
# keith hughitt (khughitt@umd.edu)
#
# Clusters genes using cor-dist (weighted combination of pearson correlation
# and euclidean distance), hierarchical clustering, and dynamicTreeCut.
#
###############################################################################
if (snakemake@config[['verbose']]) {
    message("Output:")
    message(snakemake@output)
}

###############################################################################
#
# Approach 1: Use pre-existing clustering
#
###############################################################################
if (snakemake@config[['gene_clusters']] != '') {
    message(sprintf("Using existing gene clustering located at: %s",
                    snakemake@config[['gene_clusters']]))
    # CSV
    if (endsWith(snakemake@config[['gene_clusters']], 'csv')) {
        clusters <- read.csv(snakemake@config[['gene_clusters']])
    } else {
        # Tab-delimited (consensus network)
        clusters <- read.delim(snakemake@config[['gene_clusters']], sep='\t')
        colnames(clusters) <- c('gene', 'color', 'cluster', 'description')

        # drop description (not needed and contains commas)
        clusters <- clusters[,c('gene', 'color', 'cluster')]
    }
} else {
###############################################################################
#
# Approach 2: Generate clusters using hierarchical clustering
#
###############################################################################
    # load libraries
    suppressMessages(library("dplyr"))
    suppressMessages(library("dynamicTreeCut"))
    suppressMessages(library("reshape2"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("tibble"))

    # make reproducible
    set.seed(1)

    # load count table
    expr <- as.matrix(read.csv(snakemake@config[['count_table']], row.names=1))

    # filter out any genes with zero variance (will result in NA's in cor matrix)
    zero_var_mask <- apply(expr, 1, var) != 0
    expr <- expr[zero_var_mask,]

    # size factor normalization (cpm)
    expr <- prop.table(expr, 2) * 1000000

    # generate a similarity matrix using cor-dist (weighted combination of pearson correltion and
    # euclidean distance)
    cor_matrix  <- cor(t(expr))
    dist_matrix <- log1p(as.matrix(dist(expr, diag=TRUE, upper=TRUE)))
    dist_matrix <- 1 - (dist_matrix / max(dist_matrix))

    similarity_matrix <- sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)

    # shift to range (0,1)
    similarity_matrix <- ((1 + similarity_matrix) / 2)

    # construct a hierarchical clustering denodrogram and split into separate
    # clusters
    gene_tree <- hclust(as.dist(1 - similarity_matrix), method="average")
    cluster_labels <- cutreeDynamicTree(dendro=gene_tree,
                                        minModuleSize=snakemake@config[['min_cluster_size']],
                                        deepSplit=snakemake@config[['deep_split']])
    clusters <- data.frame(gene=row.names(expr), cluster=cluster_labels)

    # generate diagnostic plot of expression clusters
    expr_flat <- melt(expr)
    colnames(expr_flat) <- c('gene', 'sample', 'value')
    expr_flat <- merge(expr_flat, clusters, by='gene')

    plt_dat <- expr_flat %>% 
        group_by(cluster, sample) %>%
        summarize(expr=mean(value))

    # save plot to disk
    plot_filepath <- sub('__snakemake_dynamic__.csv', 'clusters.png', 
                        snakemake@output[[1]])
    png(plot_filepath)
    ggplot(plt_dat, aes(x=sample, y=log(expr), group=cluster, color=factor(cluster))) +
        geom_line()
    dev.off()

    # save clustering results as a single file
    #outfile <- sub('__snakemake_dynamic__', 'clusters', snakemake@output[[1]])
    outfile <- 'features/coex_clusters.csv'
    write.csv(clusters, file=outfile, quote=FALSE, row.names=FALSE)
}

# save clustering results to separate files
for (cluster in unique(clusters$cluster)) {
    # skip unassigned genes
    if (cluster == 0) {
        next
    }

    # save cluster membership list to csv
    outfile <- sub('__snakemake_dynamic__', cluster, snakemake@output[[1]])

    write.csv(clusters[clusters$cluster == cluster,], file=outfile, 
              quote=FALSE, row.names=FALSE)
}

