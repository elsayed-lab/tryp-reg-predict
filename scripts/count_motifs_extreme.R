#!/usr/bin/env Rscript
###############################################################################
#
# count_motifs_extreme.R
# keith hughitt (khughitt@umd.edu)
#
# Loads detected motifs and creates a utr-motif matrix quantifying the
# presence of motifs for each gene feature: 5' UTR, 3' UTR, CDS, and integenic
# region.
#
###############################################################################
library('Biostrings')
library('seqLogo')
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

    message(sprintf('running load_meme_pwm: %s', filepath))

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
build_dir <- dirname(snakemake@output[[1]])

infiles_5utr <- Sys.glob(file.path(build_dir, 'extreme', '*', '5utr', '*', '*', '*.meme'))
infiles_3utr <- Sys.glob(file.path(build_dir, 'extreme', '*', '3utr', '*', '*', '*.meme'))

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
colnames(utr5_counts) <- sub('MEMEoutput.meme', '', 
                             sub('build/motifs/', '', names(motifs_5utr)))

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
colnames(utr3_counts) <- sub('MEMEoutput.meme', '', 
                             sub('build/motifs/', '', names(motifs_3utr)))

# normalize gene entries across utr5/utr3 data frames
utr5_missing <- row.names(utr3_counts)[!row.names(utr3_counts) %in% row.names(utr5_counts)]
utr3_missing <- row.names(utr5_counts)[!row.names(utr5_counts) %in% row.names(utr3_counts)]

if (length(utr5_missing) > 0) {
    utr5_missing_entries <- matrix(0, nrow=length(utr5_missing), ncol=ncol(utr5_counts))
    rownames(utr5_missing_entries) <- utr5_missing
    colnames(utr5_missing_entries) <- colnames(utr5_counts)
    utr5_counts <- rbind(utr5_counts, utr5_missing_entries)
}

if (length(utr3_missing) > 0) {
    utr3_missing_entries <- matrix(0, nrow=length(utr3_missing), ncol=ncol(utr3_counts))
    rownames(utr3_missing_entries) <- utr3_missing
    colnames(utr3_missing_entries) <- colnames(utr3_counts)
    utr3_counts <- rbind(utr3_counts, utr3_missing_entries)
}


# combine results and save output
results <- merge(utr5_counts, utr3_counts, by='row.names')
rownames(results) <- results$Row.names
results <- results[,!colnames(results) == 'Row.names']

write.csv(results, file=snakemake@output[[1]])

