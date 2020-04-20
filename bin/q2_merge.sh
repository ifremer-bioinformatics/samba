#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Script name: q2_merge.sh                                                ####
##                                                                           ##
## Purpose of script:  Merge ASV tables and repseq sequences                 ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2020-03-23                                               ####
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

# Arguments 
args=("$@")
table_dir=${args[0]}
repseq_dir=${args[1]}
merged_table=${args[2]}
merged_seq=${args[3]}
output_merge=${args[4]}
logcmd=${args[5]}

#find all tables in table_dir directory
cmdoptions=""
for table in $table_dir/*
do
   cmdoptions="$cmdoptions --i-tables $table "
done

#merge tables with qiime2
cmd="qiime feature-table merge $cmdoptions --o-merged-table $merged_table"
echo $cmd > $logcmd
eval $cmd

#find all repseq in repseq_dir directory
cmdoptions=""
for repseq in $repseq_dir/*
do
   cmdoptions="$cmdoptions --i-data $repseq "
done

#merge repseqs with qiime2
cmd="qiime feature-table merge-seqs $cmdoptions --o-merged-data $merged_seq"
echo $cmd >> $logcmd
eval $cmd

#export results
cmd="qiime tools export --input-path $merged_table --output-path $output_merge"
echo $cmd >> $logcmd
eval $cmd


