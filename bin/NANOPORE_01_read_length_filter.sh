#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Length filtering of raw nanopore reads                 ##
##                                                                           ##
###############################################################################

args=("$@")
CPUS=${args[0]}
L_MIN=${args[1]}
L_MAX=${args[2]}
FASTQ=${args[3]}
LOGCMD=${args[4]}

SAMPLE=${FASTQ%%.*}

# Run seqkit
CMD="seqkit seq -m ${L_MIN} -M ${L_MAX} ${FASTQ} | pigz -p ${CPUS} > ${SAMPLE}.filtered.fastq.gz"
echo ${CMD} > ${LOGCMD}
eval ${CMD}
