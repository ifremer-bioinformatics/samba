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
metadata=${args[4]}
merged_tableqzv=${args[5]}
output_merge=${args[6]}
merged_seqqzv=${args[7]}
logcmd=${args[8]}

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

#Generate visual and tabular summaries of a feature table
cmd="qiime feature-table summarize \
    --verbose \
    --i-table $merged_table \
    --o-visualization $merged_tableqzv \
    --m-sample-metadata-file $metadata"
echo $cmd >> $logcmd
eval $cmd

#Generate tabular view of feature identifier to sequence mapping
cmd="qiime feature-table tabulate-seqs \
    --verbose \
    --i-data $merged_seq \
    --o-visualization $merged_seqqzv"
echo $cmd >> $logcmd
eval $cmd

#Export data to html
cmd="qiime tools export \
    --input-path $merged_seqqzv \
    --output-path $output_merge"
echo $cmd >> $logcmd
eval $cmd

cmd="qiime tools export \
    --input-path $merged_tableqzv \
    --output-path $output_merge"
echo $cmd >> $logcmd
echo $cmd 
eval $cmd

cmd="qiime tools export \
    --input-path $merged_table \
    --output-path $output_merge"
echo $cmd >> $logcmd
eval $cmd
