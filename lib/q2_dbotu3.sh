#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Script name: q2_dbotu3.sh                                               ####
##                                                                           ##
## Purpose of script: Clustering of AVSs to reduce PCR errors                ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2020-02-25                                               ####
## Modified on: 2020-02-25                                                 ####
##                                                                           ##
## Email: samba-sebimer@ifremer.fr                                         ####
##                                                                           ##
## Copyright (c) SeBiMER, february-2020                                    ####
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
## q2_dbotu3.sh ${table} ${seqs} ${metadata_dbotu3} dbotu3_details.txt dbotu3_seqs.qza dbotu3_seqs.qzv dbotu3_table.qza dbotu3_table.qzv dbotu3_output ${params.dbotu3.gen-crit} ${params.dbotu3.abund-crit} ${params.dbotu3.pval-crit} completecmd > q2_dbotu3.log 2>&1

# Arguments 
args=("$@") 

table=${args[0]}
seqs=${args[1]}
metadata=${args[2]}
dbotu3_details=${args[3]}
dbotu3_seqsqza=${args[4]}
dbotu3_seqsqzv=${args[5]}
dbotu3_tableqza=${args[6]}
dbotu3_tableqzv=${args[7]}
dbotu3_output=${args[8]}
gen-crit=${args[9]}
abund-crit=${args[10]}
pval-crit=${args[11]}
logcmd=${args[12]}

# Run dbotu2
cmd="qiime dbotu-q2 call-otus --verbose --i-table $table --i-sequences $seqs --p-gen-crit $gen-crit --p-abund-crit $abund-crit --p-pval-crit $pval-crit --o-representative-sequences $dbotu3_seqsqza --o-dbotu-table $dbotu3_tableqza > $dbotu3_details"
echo $cmd > $logcmd
eval $cmd

# Generate visual and tabular summaries of a feature table
cmd="qiime feature-table summarize --verbose --i-table $dbotu3_tableqza --o-visualization $dbotu3_tableqzv --m-sample-metadata-file $metadata"
echo $cmd >> $logcmd
eval $cmd

# Generate tabular view of feature identifier to sequence mapping
cmd="qiime feature-table tabulate-seqs --verbose --i-data $dbotu3_seqsqza --o-visualization $dbotu3_seqsqzv"
echo $cmd >> $logcmd
eval $cmd

# Export data to html
cmd="qiime tools export --input-path $dbotu3_seqsqzv --output-path $dbotu3_output"
echo $cmd >> $logcmd
eval $cmd

cmd="qiime tools export --input-path $dbotu3_tableqzv --output-path $dbotu3_output"
echo $cmd >> $logcmd
eval $cmd

cmd="qiime tools export --input-path $dbotu3_tableqza --output-path $dbotu3_output"
echo $cmd >> $logcmd
eval $cmd

# Export dbotu3 details to outdir
cmd="cp $dbotu3_details $dbotu3_output/"
echo $cmd >> $logcmd
eval $cmd
