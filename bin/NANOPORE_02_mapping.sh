#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Mapping Nanopore reads against taxonomic database      ##
##                    using minimap2                                         ##
##                                                                           ##
###############################################################################

args=("$@")
CPUS=${args[0]}
BATCH_SIZE=${args[1]}
PRESET=${args[2]}
DB=${args[3]}
FASTQ=${args[4]}
BAM_OUTPUT=${args[5]}
BAM_INDEX_LOG=${args[6]}
LOGCMD=${args[7]}

# Run minimap2
CMD="minimap2 -t ${CPUS} -K ${BATCH_SIZE} -ax ${PRESET} -L ${DB} ${FASTQ} | samtools view -h -F0xe00 | samtools sort -o ${BAM_OUTPUT} -O bam -"
echo ${CMD} > ${LOGCMD}
eval ${CMD}
