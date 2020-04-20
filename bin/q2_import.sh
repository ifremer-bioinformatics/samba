#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Script name: q2_import.sh                                               ####
##                                                                           ##
## Purpose of script: Check integrity of raw data                            ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2020-01-31                                                 ####
##                                                                           ##
## Email: samba-sebimer@ifremer.fr                                         ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## This program is free software: you can redistribute it and/or modify it   ##
## under the terms of the GNU Affero General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or         ##
## (at your option) any later version.                                       ## 
##                                                                           ##
## License at https://www.gnu.org/licenses/agpl-3.0.txt                      ##
##                                                                           ##
## This program is distributed in the hope that it will be useful, but       ##
## WITHOUT ANY WARRANTY; without even the implied warranty of                ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                      ##
## See the GNU Affero General Public License for more details.               ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Command run by nextflow :
### q2_import.sh ${manifest} data.qza data.qzv import_output complete_import_cmd > q2_import.log 2>&1

# Arguments 
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
