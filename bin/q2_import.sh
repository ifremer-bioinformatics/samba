#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Check integrity of raw data                            ##
##                                                                           ##
###############################################################################

args=("$@")
singleEnd=${args[0]}
manifest=${args[1]}
dataqza=${args[2]}
dataqzv=${args[3]}
import_output=${args[4]}
logcmd=${args[5]}

#Import all samples files listed in manifest to qiime data structure
if ${singleEnd}; then
    cmdoptions="--type 'SampleData[SequencesWithQuality]' --input-format SingleEndFastqManifestPhred33V2"
else
    cmdoptions="--type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2" 
fi
cmd="qiime tools import --input-path $manifest --output-path $dataqza $cmdoptions"
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
