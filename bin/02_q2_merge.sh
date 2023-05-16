#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Merge multiple metabarcoding sequencing runs           ##
##                                                                           ##
###############################################################################

args=("$@")
TABLE_DIR=${args[0]}
MERGED_TABLE_QZA=${args[1]}
REPSEQ_DIR=${args[2]}
MERGED_REPSEQ_QZA=${args[3]}
MERGED_TABLE_QZV=${args[4]}
METADATA=${args[5]}
MERGED_REPSEQ_QZV=${args[6]}
MERGE_OUTPUT=${args[7]}
LOGCMD=${args[8]}

# Find all tables of the different runs in the merge directory
options=""
for table in ${TABLE_DIR}/*
do
   options="$options --i-tables $table "
done

# Merge tables with QIIME 2
CMD="qiime feature-table merge $options --o-merged-table ${MERGED_TABLE_QZA}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Find all ASV representative sequences of the diffrent runs in the merge directory
options=""
for repseq in ${REPSEQ_DIR}/*
do
   options="$options --i-data $repseq "
done

# Merge ASV representative sequences with QIIME 2
CMD="qiime feature-table merge-seqs $options --o-merged-data ${MERGED_REPSEQ_QZA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Generate visual and tabular summaries of the MERGED ASV table
CMD="qiime feature-table summarize --verbose --i-table ${MERGED_TABLE_QZA} --o-visualization ${MERGED_TABLE_QZV} --m-sample-metadata-file ${METADATA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Generate tabular view of ASV sequences
CMD="qiime feature-table tabulate-seqs --verbose --i-data ${MERGED_REPSEQ_QZA} --o-visualization ${MERGED_REPSEQ_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export all data to an QIIME 2 html report
CMD="qiime tools export --input-path ${MERGED_REPSEQ_QZV} --output-path ${MERGE_OUTPUT};
qiime tools export --input-path ${MERGED_TABLE_QZV} --output-path ${MERGE_OUTPUT} ;
qiime tools export --input-path ${MERGED_TABLE_QZA} --output-path ${MERGE_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert BIOM table to TSV
CMD="biom convert -i ${MERGE_OUTPUT}/feature-table.biom -o ${MERGE_OUTPUT}/merged_asv_table.tsv --to-tsv ;
sed -i '1d' ${MERGE_OUTPUT}/merged_asv_table.tsv ;
sed -i 's/#OTU ID/ASV_ID/g' ${MERGE_OUTPUT}/merged_asv_table.tsv"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
