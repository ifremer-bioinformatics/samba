#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: All Dada2 processes                                    ##
##                                                                           ##
###############################################################################

args=("$@") 

SINGLE_END=${args[0]}
DADA2_INPUT_DATA=${args[1]}
F_TRIM_LEFT=${args[2]}
R_TRIM_LEFT=${args[3]}
F_TRUNC_LENGTH=${args[4]}
R_TRUNC_LENGTH=${args[5]}
TRUNC_Q=${args[6]}
F_MAX_EE=${args[7]}
R_MAX_EE=${args[8]}
POOLING_METHOD=${args[9]}
CHIMERA_METHOD=${args[10]}
CPUS=${args[11]}
REP_SEQS_QZA=${args[12]}
TABLE_QZA=${args[13]}
DADA2_STATS_QZA=${args[14]}
DADA2_STATS_QZV=${args[15]}
TABLE_QZV=${args[16]}
METADATA=${args[17]}
REP_SEQS_QZV=${args[18]}
DADA2_OUTPUT=${args[19]}
LOGCMD=${args[20]}

if ${SINGLE_END}
then
    options="qiime dada2 denoise-single --p-trim-left ${F_TRIM_LEFT} --p-trunc-len ${F_TRUNC_LENGTH} --p-max-ee ${F_MAX_EE}"
else
    options="qiime dada2 denoise-paired --p-trim-left-f ${F_TRIM_LEFT} --p-trim-left-r ${R_TRIM_LEFT} --p-trunc-len-f ${F_TRUNC_LENGTH} --p-trunc-len-r ${R_TRUNC_LENGTH} --p-max-ee-f ${F_MAX_EE} --p-max-ee-r ${R_MAX_EE}"
fi

# Run DADA2
CMD="${options} --verbose --i-demultiplexed-seqs ${DADA2_INPUT_DATA} --p-trunc-q ${TRUNC_Q} --p-pooling-method ${POOLING_METHOD} --p-chimera-method ${CHIMERA_METHOD} --p-n-threads ${CPUS} --o-representative-sequences ${REP_SEQS_QZA} --o-table ${TABLE_QZA} --o-denoising-stats ${DADA2_STATS_QZA}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Generate a tabular view of DADA2 process stats
CMD="qiime metadata tabulate --verbose --m-input-file ${DADA2_STATS_QZA} --o-visualization ${DADA2_STATS_QZV}"
echo ${CMD} >> ${LOGCMD} 
eval ${CMD}

# Generate visual and tabular summaries of the DADA2 ASV table
CMD="qiime feature-table summarize --verbose --i-table ${TABLE_QZA} --o-visualization ${TABLE_QZV} --m-sample-metadata-file ${METADATA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Generate tabular view of ASV sequences
CMD="qiime feature-table tabulate-seqs --verbose --i-data ${REP_SEQS_QZA} --o-visualization ${REP_SEQS_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export all data to an QIIME 2 html report
CMD="qiime tools export --input-path ${REP_SEQS_QZV} --output-path ${DADA2_OUTPUT} ; 
qiime tools export --input-path ${TABLE_QZA} --output-path ${DADA2_OUTPUT} ;
qiime tools export --input-path ${TABLE_QZV} --output-path ${DADA2_OUTPUT} ;
qiime tools export --input-path ${DADA2_STATS_QZV} --output-path ${DADA2_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
