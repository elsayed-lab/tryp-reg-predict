#!/usr/bin/env sh
#
# EXTREME wrapper for tryp regulatory element prediction pipeline
#
python2 ../src/fasta-dinucleotide-shuffle.py -f GM12878_NRSF_ChIP.fasta > GM12878_NRSF_ChIP_shuffled.fasta
python2 ../src/GappedKmerSearch.py -l 8 -ming 0 -maxg 10 -minsites 5 GM12878_NRSF_ChIP.fasta GM12878_NRSF_ChIP_shuffled.fasta GM12878_NRSF_ChIP.words
perl ../src/run_consensus_clusering_using_wm.pl GM12878_NRSF_ChIP.words 0.3
python2 ../src/Consensus2PWM.py GM12878_NRSF_ChIP.words.cluster.aln GM12878_NRSF_ChIP.wm
