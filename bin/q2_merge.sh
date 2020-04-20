#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script:  Merge ASV tables and repseq sequences                 ##
##                                                                           ##
###############################################################################

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


