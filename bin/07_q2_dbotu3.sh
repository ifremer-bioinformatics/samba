#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Clustering of AVSs to reduce PCR errors                ##
##                                                                           ##
###############################################################################

args=("$@") 

TABLE=${args[0]}
REP_SEQS=${args[1]}
GENETIC_CRITERION=${args[2]}
ABUNDANCE_CRITERION=${args[3]}
PVALUE_CRITERION=${args[4]}
dbOTU3_SEQS_QZA=${args[5]}
dbOTU3_TABLE_QZA=${args[6]}
dbOTU3_DETAILS=${args[7]}
dbOTU3_TABLE_QZV=${args[8]}
METADATA=${args[9]}
dbOTU3_SEQS_QZV=${args[10]}
dbOTU3_OUTPUT=${args[11]}
LOGCMD=${args[12]}

# Run dbOTU3
CMD="qiime dbotu-q2 call-otus --verbose --i-table ${TABLE} --i-sequences ${REP_SEQS} --p-gen-crit ${GENETIC_CRITERION} --p-abund-crit ${ABUNDANCE_CRITERION} --p-pval-crit ${PVALUE_CRITERION} --o-representative-sequences ${dbOTU3_SEQS_QZA} --o-dbotu-table ${dbOTU3_TABLE_QZA} > ${dbOTU3_DETAILS}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Generate visual and tabular summaries of the DADA2 ASV table
CMD="qiime feature-table summarize --verbose --i-table ${dbOTU3_TABLE_QZA} --o-visualization ${dbOTU3_TABLE_QZV} --m-sample-metadata-file ${METADATA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Generate tabular view of ASV sequences
CMD="qiime feature-table tabulate-seqs --verbose --i-data ${dbOTU3_SEQS_QZA} --o-visualization ${dbOTU3_SEQS_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export all data to an QIIME 2 html report
CMD="qiime tools export --input-path ${dbOTU3_SEQS_QZV} --output-path ${dbOTU3_OUTPUT} ;
qiime tools export --input-path ${dbOTU3_TABLE_QZA} --output-path ${dbOTU3_OUTPUT} ;
qiime tools export --input-path ${dbOTU3_TABLE_QZV} --output-path ${dbOTU3_OUTPUT} ;
cp ${dbOTU3_DETAILS} ${dbOTU3_OUTPUT}/"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
