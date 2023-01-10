#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Filtering ASV table & sequences based on taxonomy      ##
##                                                                           ##
###############################################################################

args=("$@")

FILTER_BY_ID=${args[0]}
SAMPLES_TO_REMOVE=${args[1]}
ASV_TABLE=${args[2]}
ID_FILTERED_TABLE_QZA=${args[3]}
ASV_SEQS=${args[4]}
ID_FILTERED_SEQS_QZA=${args[5]}
FILTER_BY_FREQUENCY=${args[6]}
MIN_FREQUENCY_SAMPLE=${args[7]}
MIN_FREQUENCY_ASV=${args[8]}
CONTINGENCY_ASV=${args[9]}
FREQ_FILTERED_TABLE_QZA=${args[10]}
FREQ_FILTERED_SEQS_QZA=${args[11]}
FINAL_FILTERED_TABLE_QZV=${args[12]}
METADATA=${args[13]}
FINAL_FILTERED_SEQS_QZV=${args[14]}
FILTERED_DATA_OUTPUT=${args[15]}
TAX_TSV=${args[16]}
FILTERED_DATA_TABLE_BIOM=${args[17]}
FILTERED_DATA_TABLE_TSV=${args[18]}
LOGCMD=${args[19]}

# Filtering data based on sample ID
if "${FILTER_BY_ID}"
then
    echo SampleID > list_samples_to_remove.tsv
    echo "${SAMPLES_TO_REMOVE}" >> list_samples_to_remove.tsv
    if [[ "${SAMPLES_TO_REMOVE}" =~ "," ]]
    then 
        sed -i 's/,/\n/g' list_samples_to_remove.tsv
    fi
    CMD="qiime feature-table filter-samples --i-table ${ASV_TABLE} --m-metadata-file list_samples_to_remove.tsv --o-filtered-table ${ID_FILTERED_TABLE_QZA} --p-filter-empty-features --p-exclude-ids ;
    qiime feature-table filter-seqs --i-data ${ASV_SEQS} --i-table ${ID_FILTERED_TABLE_QZA} --o-filtered-data ${ID_FILTERED_SEQS_QZA}"
    echo ${CMD} > ${LOGCMD}
    eval ${CMD}
fi

# Filtering data based total-frequency and contingency of samples and ASVs
if "${FILTER_BY_FREQUENCY}"
then
    if "${FILTER_BY_ID}"
    then
        ASV_TABLE=${ID_FILTERED_TABLE_QZA}
    fi
    CMD="qiime feature-table filter-samples --i-table ${ASV_TABLE} --p-min-frequency ${MIN_FREQUENCY_SAMPLE}  --o-filtered-table tmp_frequency_filtered_table.qza --p-filter-empty-features ;
    qiime feature-table filter-features --i-table tmp_frequency_filtered_table.qza --p-min-frequency ${MIN_FREQUENCY_ASV} --p-min-samples ${CONTINGENCY_ASV} --o-filtered-table ${FREQ_FILTERED_TABLE_QZA} ;
    qiime feature-table filter-seqs --i-data ${ASV_SEQS} --i-table ${FREQ_FILTERED_TABLE_QZA} --o-filtered-data ${FREQ_FILTERED_SEQS_QZA}"
    echo ${CMD} >> ${LOGCMD}
    eval ${CMD}
fi

# Export all data to an QIIME 2 html report
if "${FILTER_BY_FREQUENCY}"
then
    cp ${FREQ_FILTERED_TABLE_QZA} final_asv_table_filtered.qza
    cp ${FREQ_FILTERED_SEQS_QZA} final_asv_seqs_filtered.qza
else
    cp ${ID_FILTERED_TABLE_QZA} final_asv_table_filtered.qza
    cp ${ID_FILTERED_SEQS_QZA} final_asv_seqs_filtered.qza
fi

CMD="qiime feature-table summarize --i-table final_asv_table_filtered.qza --o-visualization ${FINAL_FILTERED_TABLE_QZV} --m-sample-metadata-file ${METADATA} ;
qiime feature-table tabulate-seqs --i-data final_asv_seqs_filtered.qza --o-visualization ${FINAL_FILTERED_SEQS_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

CMD="qiime tools export --input-path ${FINAL_FILTERED_TABLE_QZV} --output-path ${FILTERED_DATA_OUTPUT} ;
qiime tools export --input-path ${FINAL_FILTERED_SEQS_QZV} --output-path ${FILTERED_DATA_OUTPUT} ;
qiime tools export --input-path final_asv_table_filtered.qza --output-path ${FILTERED_DATA_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Add taxonomy to count table (biom format)
CMD="biom add-metadata -i ${FILTERED_DATA_OUTPUT}/feature-table.biom --observation-metadata-fp ${TAX_TSV} -o ${FILTERED_DATA_TABLE_BIOM} --sc-separated taxonomy"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert biom table to tabular
CMD="biom convert -i ${FILTERED_DATA_TABLE_BIOM} -o ${FILTERED_DATA_TABLE_TSV} --to-tsv --header-key taxonomy ;
cp ${FILTERED_DATA_TABLE_TSV} ${FILTERED_DATA_OUTPUT}/"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
