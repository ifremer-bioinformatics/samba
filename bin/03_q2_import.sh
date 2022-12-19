#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Check integrity of raw data                            ##
##                                                                           ##
###############################################################################

args=("$@")
SINGLE_END=${args[0]}
MANIFEST=${args[1]}
DATA_IMPORTED_QZA=${args[2]}
DATA_IMPORTED_QZV=${args[3]}
OUTPUT_EXPORT=${args[4]}
LOGCMD=${args[5]}

if ${SINGLE_END}; then
    options="--type 'SampleData[SequencesWithQuality]' --input-format SingleEndFastqManifestPhred33V2"
else
    options="--type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2" 
fi

# Import all samples FASTQ files listed in manifest to QIIME2 format

CMD="qiime tools import --input-path ${MANIFEST} --output-path ${DATA_IMPORTED_QZA} $options"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Summarize counts per sample for all samples

CMD="qiime demux summarize --verbose --i-data ${DATA_IMPORTED_QZA} --o-visualization ${DATA_IMPORTED_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export html report
CMD="qiime tools export --input-path ${DATA_IMPORTED_QZV} --output-path ${OUTPUT_EXPORT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
