#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Filtering ASV table & sequences based on taxonomy      ##
##                                                                           ##
###############################################################################

args=("$@")
FILTERING_TYPE=${args[0]}
ASV_TABLE=${args[1]}
TAX_QZA=${args[2]}
TAX_TO_FILTER=${args[3]}
FILTERED_TAX_TABLE_QZA=${args[4]}
ASV_SEQS=${args[5]}
FILTERED_TAX_SEQS_QZA=${args[6]}
FILTERED_TAX_TABLE_QZV=${args[7]}
METADATA=${args[8]}
FILTERED_TAX_SEQS_QZV=${args[9]}
FILTERED_TAX_OUTPUT=${args[10]}
TAX_TSV=${args[11]}
FILTERED_TAX_TABLE_BIOM=${args[12]}
FILTERED_TAX_TABLE_TSV=${args[13]}
LOGCMD=${args[14]}

if [ "${FILTERING_TYPE}" = "exclude" ]
then
    # Exclude ASV from the ASV table based on taxonomy
    CMD="qiime taxa filter-table --i-table ${ASV_TABLE} --i-taxonomy ${TAX_QZA} --p-exclude ${TAX_TO_FILTER} --o-filtered-table ${FILTERED_TAX_TABLE_QZA}"
    echo ${CMD} > ${LOGCMD}
    eval ${CMD}
    
    # Exclude ASV sequences based on taxonomy
    CMD="qiime taxa filter-seqs --i-sequences ${ASV_SEQS} --i-taxonomy ${TAX_QZA} --p-exclude ${TAX_TO_FILTER} --o-filtered-sequences ${FILTERED_TAX_SEQS_QZA}"
    echo ${CMD} >> ${LOGCMD}
    eval ${CMD}
elif [ "${FILTERING_TYPE}" = "include" ]
then
    # Only include ASV from the ASV table based on taxonomy
    CMD="qiime taxa filter-table --i-table ${ASV_TABLE} --i-taxonomy ${TAX_QZA} --p-include ${TAX_TO_FILTER} --o-filtered-table ${FILTERED_TAX_TABLE_QZA}"
    echo ${CMD} > ${LOGCMD}
    eval ${CMD}

    # Only include ASV sequences based on taxonomy
    CMD="qiime taxa filter-seqs --i-sequences ${ASV_SEQS} --i-taxonomy ${TAX_QZA} --p-include ${TAX_TO_FILTER} --o-filtered-sequences ${FILTERED_TAX_SEQS_QZA}"
    echo ${CMD} >> ${LOGCMD}
    eval ${CMD}
fi

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
