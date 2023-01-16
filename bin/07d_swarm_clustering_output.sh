#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Reformat swarm output results for future analyses      ##
##                                                                           ##
###############################################################################

args=("$@")

SWARM_TABLE_TSV=${args[0]}
SWARM_TABLE_BIOM=${args[1]}
SWARM_TABLE_QZA=${args[2]}
SWARM_ASV_SEQS_FA=${args[3]}
SWARM_ASV_SEQS_QZA=${args[4]}
SWARM_TABLE_QZV=${args[5]}
METADATA=${args[6]}
SWARM_CLUSTERING_OUTPUT=${args[7]}
SWARM_ASV_SEQS_QZV=${args[8]}
LOGCMD=${args[9]}

# Convert TSV swarm table to BIOM format
CMD='biom convert -i ${SWARM_TABLE_TSV} -o ${SWARM_TABLE_BIOM} --to-hdf5 --table-type="OTU table"'
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert swarm table to QIIME 2 format
CMD="qiime tools import --input-path ${SWARM_TABLE_BIOM} --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ${SWARM_TABLE_QZA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert FASTA file of the swarm reference ASVs to QIIME 2 format
CMD="qiime tools import --input-path ${SWARM_ASV_SEQS_FA} --output-path ${SWARM_ASV_SEQS_QZA} --type 'FeatureData[Sequence]'"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export all data to an QIIME 2 html report
CMD="qiime feature-table summarize --i-table ${SWARM_TABLE_QZA} --o-visualization ${SWARM_TABLE_QZV} --m-sample-metadata-file ${METADATA} ;
qiime tools export --input-path ${SWARM_TABLE_QZA} --output-path ${SWARM_CLUSTERING_OUTPUT} ;
qiime tools export --input-path ${SWARM_TABLE_QZV} --output-path ${SWARM_CLUSTERING_OUTPUT} ;
qiime feature-table tabulate-seqs --i-data ${SWARM_ASV_SEQS_QZA} --o-visualization ${SWARM_ASV_SEQS_QZV} ;
qiime tools export --input-path ${SWARM_ASV_SEQS_QZV} --output-path ${SWARM_CLUSTERING_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
