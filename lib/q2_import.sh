#!/usr/bin/env bash
## Command run by nextflow process :
### q2_import.sh ${manifest} data.qza data.qzv import_output > q2_import.log 2>&1

#inputs
manifest=$1
dataqza=$2
dataqzv=$3
import_output=$4

#Import all samples paired-end files listed in manifest to qiime data structure
qiime tools import --input-path $manifest --output-path $dataqza --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2

#Summarize counts per sample for all samples
qiime demux summarize --verbose --i-data $dataqza --o-visualization $dataqzv

#Export html report
qiime tools export --input-path $dataqzv --output-path $import_output
