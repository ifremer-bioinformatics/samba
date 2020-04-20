#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Script name: q2_ANCOM.sh                                                ####
##                                                                           ##
## Purpose of script: Differential abundance testing with ANCOM              ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2020-04-20                                               ####
## Modified on: 2020-04-20                                                 ####
##                                                                           ##
## Email: samba-sebimer@ifremer.fr                                         ####
##                                                                           ##
## Copyright (c) SeBiMER, april-2020                                       ####
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
## Command run by snakemake :
### q2_ANCOM.sh ${table} compo_table.qza ${metadata} ${params.ancom_var} ancom_${params.ancom_var}.qzv export_ancom_${params.ancom_var} ${taxonomy} collapsed_table_family.qza compo_table_family.qza ancom_${params.ancom_var}_family.qzv export_ancom_${params.ancom_var}_family collapsed_table_genus.qza compo_table_genus.qza ancom_${params.ancom_var}_genus.qzv export_ancom_${params.ancom_var}_genus completecmd > q2_ANCOM.log 2>&1

# Arguments
args=("$@")
table=${args[0]}
compo_table=${args[1]}
metadata=${args[2]}
criteria=${args[3]}
output_global=${args[4]}
export_global=${args[5]}
taxonomy=${args[6]}
collaped_table_family=${args[7]}
compo_table_family=${args[8]}
output_family=${args[9]}
export_family=${args[10]}
collaped_table_genus=${args[11]}
compo_table_genus=${args[12]}
output_genus=${args[13]}
export_genus=${args[14]}
logcmd=${args[15]}

# Build a composition table artefact file
cmd="qiime composition add-pseudocount --i-table ${table} --o-composition-table ${compo_table}"
echo $cmd > $logcmd
eval $cmd

# ANCOM analysis based on the user-selected variable
cmd="qiime composition ancom --i-table ${compo_table} --m-metadata-file ${metadata} --m-metadata-column ${criteria} --p-transform-function log --p-difference-function f_statistic --o-visualization ${output_global} ;
qiime tools export --input-path ${output_global} --output-path ${export_global}"
echo $cmd >> $logcmd
eval $cmd

# ANCOM analysis based on the user-selected variable at the family level
cmd="qiime taxa collapse --i-table ${table} --i-taxonomy ${taxonomy} --p-level 5 --o-collapsed-table ${collaped_table_family} ;
qiime composition add-pseudocount --i-table ${collaped_table_family} --o-composition-table ${compo_table_family} ;
qiime composition ancom --i-table ${compo_table_family} --m-metadata-file ${metadata} --m-metadata-column ${criteria} --p-transform-function log --p-difference-function f_statistic --o-visualization ${output_family} ;
qiime tools export --input-path ${output_family} --output-path ${export_family}"
echo $cmd >> $logcmd
eval $cmd

# ANCOM analysis based on the user-selected variable at the genus level
cmd="qiime taxa collapse --i-table ${table} --i-taxonomy ${taxonomy} --p-level 6 --o-collapsed-table ${collaped_table_genus} ;
qiime composition add-pseudocount --i-table ${collaped_table_genus} --o-composition-table ${compo_table_genus} ;
qiime composition ancom --i-table ${compo_table_genus} --m-metadata-file ${metadata} --m-metadata-column ${criteria} --p-transform-function log --p-difference-function f_statistic --o-visualization ${output_genus} ;
qiime tools export --input-path ${output_genus} --output-path ${export_genus}"
echo $cmd >> $logcmd
eval $cmd
