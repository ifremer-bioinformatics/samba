#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Reformat microDecon results for future analyses        ##
##                                                                           ##
###############################################################################

args=("$@")

DECONTAM_TABLE_TSV=${args[0]}
DECONTAM_TABLE_BIOM=${args[1]}
DECONTAM_TABLE_QZA=${args[2]}
ASV_SEQS=${args[3]}
DECONTAM_ASV_SEQS_FA=${args[4]}
DECONTAM_ASV_SEQS_QZA=${args[5]}
DECONTAM_TABLE_QZV=${args[6]}
METADATA=${args[7]}
FILTER_CONTAMINANTS_OUTPUT=${args[8]}
DECONTAM_ASV_SEQS_QZV=${args[9]}
LOGCMD=${args[10]}

# Convert TSV decontaminated table to BIOM format
CMD='biom convert -i ${DECONTAM_TABLE_TSV} -o ${DECONTAM_TABLE_BIOM} --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy'
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert decontaminated table to QIIME 2 format
CMD="qiime tools import --input-path ${DECONTAM_TABLE_BIOM} --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ${DECONTAM_TABLE_QZA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Retrieve non-contaminant ASV IDs
CMD="cut -d \$'\t' -f1 ${DECONTAM_TABLE_TSV} | sed '1d' > decontaminated_ASV_ID.txt"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Build FASTA file of the non-contaminant ASVs
CMD="seqtk subseq ${ASV_SEQS} decontaminated_ASV_ID.txt > ${DECONTAM_ASV_SEQS_FA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert FASTA file of the non-contaminant ASVs to QIIME 2 format
CMD="qiime tools import --input-path ${DECONTAM_ASV_SEQS_FA} --output-path ${DECONTAM_ASV_SEQS_QZA} --type 'FeatureData[Sequence]'"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export all data to an QIIME 2 html report
CMD="qiime feature-table summarize --i-table ${DECONTAM_TABLE_QZA} --o-visualization ${DECONTAM_TABLE_QZV} --m-sample-metadata-file ${METADATA} ;
qiime tools export --input-path ${DECONTAM_TABLE_QZA} --output-path ${FILTER_CONTAMINANTS_OUTPUT} ;
qiime tools export --input-path ${DECONTAM_TABLE_QZV} --output-path ${FILTER_CONTAMINANTS_OUTPUT} ;
qiime feature-table tabulate-seqs --i-data ${DECONTAM_ASV_SEQS_QZA} --o-visualization ${DECONTAM_ASV_SEQS_QZV} ;
qiime tools export --input-path ${DECONTAM_ASV_SEQS_QZV} --output-path ${FILTER_CONTAMINANTS_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
