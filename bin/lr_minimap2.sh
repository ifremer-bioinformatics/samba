#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Map long reads (Nanopore/PacBio) to SILVA              ##
##                                                                           ##
###############################################################################

args=("$@")
fastq=${args[0]}
sam=${args[1]}
manufacturer=${args[2]}
database=${args[3]}
cpus=${args[4]}
logcmd=${args[5]}

#Mapping of Nanopore/PacBio reads against database and filter with Samtools
#Note: -K option is used to limit ram usage
#Note: database indexed with minimap2, kmer 15
cmd="minimap2 -t ${cpus} -K 25M -ax ${manufacturer} -L ${database} ${fastq} | samtools view -F0xe80 > ${sam}"
echo $cmd > $logcmd
eval $cmd
