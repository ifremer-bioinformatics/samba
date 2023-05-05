#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Identify optimal optimizing microbiome rRNA gene       ##
##                    trimming parameters for DADA2                          ##
##                                                                           ##
###############################################################################

args=("$@")

RAW_DATA_DIR=${args[0]}
AMPLICON_LENGTH=${args[1]}
F_PRIMER=${args[2]}
R_PRIMER=${args[3]}
FIGARO_OUTPUT=${args[4]}
LOGCMD=${args[5]}

# Retrieve primer length
LENGTH_PRIMER_F=$(echo ${F_PRIMER} | awk '{print length}')
LENGTH_PRIMER_R=$(echo ${R_PRIMER} | awk '{print length}')

# Run FIGARO
CMD="figaro -i ${RAW_DATA_DIR} -o ${FIGARO_OUTPUT} -a ${AMPLICON_LENGTH} -f ${LENGTH_PRIMER_F} -r ${LENGTH_PRIMER_R}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}
