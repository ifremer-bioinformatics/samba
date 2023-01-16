#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Clustering of AVSs using swarm                         ##
##                                                                           ##
###############################################################################

args=("$@") 

CPUS=${args[0]}
ASV_SEQS_FASTA=${args[1]}
SWARM_ASV_SEQS_FA=${args[2]}
TMP_SWARM_CLUSTER_LIST=${args[3]}
FINAL_SWARM_ASV_SEQS_FA=${args[4]}
LOGCMD=${args[5]}

# Run swarm
CMD="swarm -t ${CPUS} -f -l swarm.log -w ${SWARM_ASV_SEQS_FA} ${ASV_SEQS_FASTA} > ${TMP_SWARM_CLUSTER_LIST} "
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Reformat swarm sequence header FASTA file
CMD="cut -d '_' -f1 ${SWARM_ASV_SEQS_FA} > ${FINAL_SWARM_ASV_SEQS_FA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
