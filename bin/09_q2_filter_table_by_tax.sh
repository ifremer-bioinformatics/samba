#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Filtering ASV table & sequences based on taxonomy      ##
##                                                                           ##
###############################################################################

args=("$@")

ASV_TABLE=${args[0]}
TAX_QZA=${args[1]}
TAX_TO_EXCLUDE=${args[2]}
FILTERED_TAX_TABLE_QZA=${args[3]}
ASV_SEQS=${args[4]}
FILTERED_TAX_SEQS_QZA=${args[5]}
FILTERED_TAX_TABLE_QZV=${args[6]}
METADATA=${args[7]}
FILTERED_TAX_SEQS_QZV=${args[8]}
FILTERED_TAX_OUTPUT=${args[9]}
TAX_TSV=${args[10]}
FILTERED_TAX_TABLE_BIOM=${args[11]}
FILTERED_TAX_TABLE_TSV=${args[12]}
LOGCMD=${args[13]}

# Filtering ASV table based on taxonomy
CMD="qiime taxa filter-table --i-table ${ASV_TABLE} --i-taxonomy ${TAX_QZA} --p-exclude ${TAX_TO_EXCLUDE} --o-filtered-table ${FILTERED_TAX_TABLE_QZA}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}
    
# Filtering ASV sequences based on taxonomy
CMD="qiime taxa filter-seqs --i-sequences ${ASV_SEQS} --i-taxonomy ${TAX_QZA} --p-exclude ${TAX_TO_EXCLUDE} --o-filtered-sequences ${FILTERED_TAX_SEQS_QZA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export all data to an QIIME 2 html report
CMD="qiime feature-table summarize --i-table ${FILTERED_TAX_TABLE_QZA} --o-visualization ${FILTERED_TAX_TABLE_QZV} --m-sample-metadata-file ${METADATA} ;
qiime feature-table tabulate-seqs --i-data ${FILTERED_TAX_SEQS_QZA} --o-visualization ${FILTERED_TAX_SEQS_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

CMD="qiime tools export --input-path ${FILTERED_TAX_TABLE_QZV} --output-path ${FILTERED_TAX_OUTPUT} ;
qiime tools export --input-path ${FILTERED_TAX_SEQS_QZV} --output-path ${FILTERED_TAX_OUTPUT} ;
qiime tools export --input-path ${FILTERED_TAX_TABLE_QZA} --output-path ${FILTERED_TAX_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Add taxonomy to count table (biom format)
CMD="biom add-metadata -i ${FILTERED_TAX_OUTPUT}/feature-table.biom --observation-metadata-fp ${TAX_TSV} -o ${FILTERED_TAX_TABLE_BIOM} --sc-separated taxonomy"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert biom table to tabular
CMD="biom convert -i ${FILTERED_TAX_TABLE_BIOM} -o ${FILTERED_TAX_TABLE_TSV} --to-tsv --header-key taxonomy ;
cp ${FILTERED_TAX_TABLE_TSV} ${FILTERED_TAX_OUTPUT}/"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
