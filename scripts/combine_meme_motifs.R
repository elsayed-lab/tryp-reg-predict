#!/usr/bin/env Rscript
###############################################################################
# Loads detected motifs and creates a utr-motif matrix quantifying the
# presence of motifs for each 5' and 5'UTR.
#
# Highly correlated motifs and those with low variance across genes are
# filtered out to obtain a reduced set of informative motifs.
#
###############################################################################
library('Biostrings')
library('seqLogo')
library('caret')
options(stringsAsFactors=FALSE)

###############################################################################
#
# load_meme_pwm
# 
# Loads a MEME 4.90 formatted PWM.
#
# Creates a matrix in a format usable by the Biostrings matchPWM and countPWM
# functions.
#
# For parsing output from earlier versions of MEME, see:
# - http://cran.r-project.org/web/packages/MEET/index.html
#
###############################################################################
load_meme_pwm <- function(filepath) {
    # Read file contents
    lines <- readLines(filepath)

    # Matrix starts on line 14 and goes until next to last line
    MATRIX_START <- 14
    MATRIX_END   <- length(lines) - 1

    # Get matrix lines
    lines <- lines[MATRIX_START:MATRIX_END]

    # Create matrix
    pwm <- matrix(nrow=length(lines), ncol=4)

    i <- 1
    for (line in lines) {
        parts <- unlist(strsplit(line, '\\s'))
        pwm[i,] <- as.numeric(parts[parts != ""])
        i <- i + 1
    }
    colnames(pwm) <- c("A", "C", "G", "T")
    return(t(pwm))
}

# list of all motifs detected across all modules and all runs
build_dir <- dirname(snakemake@output$utr5)

infiles_5utr <- Sys.glob(file.path(build_dir, '*', '5utr', '*', '*', '*.meme'))
infiles_3utr <- Sys.glob(file.path(build_dir, '*', '3utr', '*', '*', '*.meme'))

# load motif PWMs
motifs_5utr <- list()

for (infile in infiles_5utr) {
    # load motif
    motifs_5utr[[infile]] <- load_meme_pwm(infile)

    # generate sequence logo
    png(sub('.meme', '.png', infile))
    seqLogo(motifs_5utr[[infile]])
    dev.off()
}

motifs_3utr <- list()

for (infile in infiles_3utr) {
    motifs_3utr[[infile]] <- load_meme_pwm(infile)

    # generate sequence logo
    png(sub('.meme', '.png', infile))
    seqLogo(motifs_3utr[[infile]])
    dev.off()
}

#
# Quantify gene-level motif presence
#

# load 5'UTR sequences
utr5_seqs <- read.csv(snakemake@config[['5utr_stats']])

# create an empty data frame to store counts
utr5_counts <- NULL

# iterate over 5'UTR sequences
for (i in 1:nrow(utr5_seqs)) {
    # get gene UTR sequence
    utr <- utr5_seqs[i,]$seq

    # vector to store motif counts for gene
    entries <- c()

    message(sprintf("Counting 5'UTR motif occurrences for gene %d/%d",
                    i, nrow(utr5_seqs)))

    # Iterate over motifs and count occurrences of each motif in the UTR
    for (motif in motifs_5utr) {
        entries <- append(entries, countPWM(motif, utr))
    }

    # Add row to data frame
    utr5_counts <- rbind(utr5_counts, entries)
}

# Convert to a data frame
utr5_counts <- data.frame(utr5_counts)
rownames(utr5_counts) <- utr5_seqs$gene
colnames(utr5_counts) <- names(motifs_5utr)

# load 3'UTR sequences
utr3_seqs <- read.csv(snakemake@config[['3utr_stats']])

# create an empty data frame to store counts
utr3_counts <- NULL

# iterate over 3'UTR sequences
for (i in 1:nrow(utr3_seqs)) {
    # get gene UTR sequence
    utr <- utr3_seqs[i,]$seq

    # vector to store motif counts for gene
    entries <- c()

    message(sprintf("Counting 3'UTR motif occurrences for gene %d/%d",
                    i, nrow(utr3_seqs)))

    # Iterate over motifs and count occurrences of each motif in the UTR
    for (motif in motifs_3utr) {
        entries <- append(entries, countPWM(motif, utr))
    }

    # Add row to data frame
    utr3_counts <- rbind(utr3_counts, entries)
}

# Convert to a data frame
utr3_counts <- data.frame(utr3_counts)
rownames(utr3_counts) <- utr3_seqs$gene
colnames(utr3_counts) <- names(motifs_3utr)

# Use the caret findCorrelation function to find and remove highly correlated
# variables

# 5' UTR
cor_mat_5utr <- cor(utr5_counts, method='spearman')
cor_mat_5utr[is.na(cor_mat_5utr)] <- 0

motifs_to_ignore <- findCorrelation(cor_mat_5utr, cutoff=0.5, verbose=FALSE)
utr5_counts <- utr5_counts[,-motifs_to_ignore]

print(sprintf("Ignoring %d/%d highly-correlated 5' UTR motifs...",
              length(motifs_to_ignore), nrow(cor_mat_5utr)))
print(sprintf("Remaining motifs: %d", ncol(utr5_counts)))

cor_mat_5utr <- cor(utr5_counts, method='spearman')
cor_mat_5utr[is.na(cor_mat_5utr)] <- 0

# 3' UTR
cor_mat_3utr <- cor(utr3_counts, method='spearman')
cor_mat_3utr[is.na(cor_mat_3utr)] <- 0

motifs_to_ignore <- findCorrelation(cor_mat_3utr, cutoff=0.5, verbose=FALSE)
utr3_counts <- utr3_counts[,-motifs_to_ignore]

print(sprintf("Ignoring %d/%d highly-correlated 3' UTR motifs...",
              length(motifs_to_ignore), nrow(cor_mat_3utr)))
print(sprintf("Remaining motifs: %d", ncol(utr3_counts)))

# save output
write.csv(utr5_counts, file=snakemake@output$utr5)
write.csv(utr3_counts, file=snakemake@output$utr3)

