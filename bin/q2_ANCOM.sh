#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Differential abundance testing with ANCOM              ##
##                                                                           ##
###############################################################################

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
