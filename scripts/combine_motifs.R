#!/usr/bin/env Rscript
###############################################################################
# Loads detected motifs and creates a utr-motif matrix quantifying the
# presence of motifs for each 5' and 5'UTR.
#
# Highly correlated motifs and those with low variance across genes are
# filtered out to obtain a reduced set of informative motifs.
#
###############################################################################

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
load_meme_pwm = function(filepath) {
    # Read file contents
    lines = readLines(filepath)

    # Matrix starts on line 14 and goes until next to last line
    MATRIX_START = 14
    MATRIX_END   = length(lines) - 1

    # Get matrix lines
    lines = lines[MATRIX_START:MATRIX_END]

    # Create matrix
    pwm = matrix(nrow=length(lines), ncol=4)

    i = 1
    for (line in lines) {
        parts = unlist(strsplit(line, '\\s'))
        pwm[i,] = as.numeric(parts[parts != ""])
        i = i + 1
    }
    colnames(pwm) = c("A", "C", "G", "T")
    return(t(pwm))
}
