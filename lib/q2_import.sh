#!/usr/bin/env bash
## Command run by nextflow process :
### q2_import.sh ${manifest} data.qza data.qzv import_output complete_import_cmd > q2_import.log 2>&1
# store arguments in a special array 
args=("$@")

manifest=${args[0]}
dataqza=${args[1]}
dataqzv=${args[2]}
import_output=${args[3]}
logcmd=${args[4]}

#Import all samples paired-end files listed in manifest to qiime data structure
cmd="qiime tools import --input-path $manifest --output-path $dataqza --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2"
echo $cmd > $logcmd
eval $cmd

#Summarize counts per sample for all samples
cmd="qiime demux summarize --verbose --i-data $dataqza --o-visualization $dataqzv"
echo $cmd >> $logcmd
eval $cmd

#Export html report
cmd="qiime tools export --input-path $dataqzv --output-path $import_output"
echo $cmd >> $logcmd
eval $cmd
