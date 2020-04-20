#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: All Dada2 processes                                    ##
##                                                                           ##
###############################################################################

args=("$@") 

singleEnd=${args[0]}
trimmed_data=${args[1]}
metadata=${args[2]}
rep_seqsqza=${args[3]}
rep_seqsqzv=${args[4]}
tableqza=${args[5]}
tableqzv=${args[6]}
statsqza=${args[7]}
statsqzv=${args[8]}
dada2_output=${args[9]}
trimLeft=${args[10]}
trimRigth=${args[11]}
FtruncLen=${args[12]}
RtruncLen=${args[13]}
FmaxEE=${args[14]}
RmaxEE=${args[15]}
minQ=${args[16]}
chimeras=${args[17]}
cpus=${args[18]}
logcmd=${args[19]}

#Run dada2 : denoises paired-end sequences, dereplicates them and filters chimeras
if ${singleEnd}; then
    cmdoptions="qiime dada2 denoise-single --p-trim-left $trimLeft --p-trunc-len $FtruncLen --p-max-ee $FmaxEE"
else
    cmdoptions="qiime dada2 denoise-paired --p-trim-left-f $trimLeft --p-trim-left-r $trimRigth --p-trunc-len-f $FtruncLen --p-trunc-len-r $RtruncLen --p-max-ee-f $FmaxEE --p-max-ee-r $RmaxEE" 
fi
cmd="$cmdoptions \
    --verbose \
    --i-demultiplexed-seqs $trimmed_data \
    --p-trunc-q $minQ \
    --p-chimera-method $chimeras \
    --p-n-threads $cpus \
    --o-representative-sequences $rep_seqsqza \
    --o-table $tableqza  \
    --o-denoising-stats $statsqza"
echo $cmd > $logcmd
eval $cmd

#Generate a tabular view of taxonomy metadata
cmd="qiime metadata tabulate \
    --verbose \
    --m-input-file $statsqza \
    --o-visualization $statsqzv"
echo $cmd >> $logcmd 
eval $cmd

#Generate visual and tabular summaries of a feature table
cmd="qiime feature-table summarize \
    --verbose \
    --i-table $tableqza \
    --o-visualization $tableqzv \
    --m-sample-metadata-file $metadata"
echo $cmd >> $logcmd
echo $cmd
eval $cmd

#Generate tabular view of feature identifier to sequence mapping
cmd="qiime feature-table tabulate-seqs \
    --verbose \
    --i-data $rep_seqsqza \
    --o-visualization $rep_seqsqzv"
echo $cmd >> $logcmd
eval $cmd

#Export data to html
cmd="qiime tools export \
    --input-path $rep_seqsqzv \
    --output-path $dada2_output"
echo $cmd >> $logcmd
eval $cmd

cmd="qiime tools export \
    --input-path $tableqzv \
    --output-path $dada2_output"
echo $cmd >> $logcmd
echo $cmd 
eval $cmd

cmd="qiime tools export \
    --input-path $statsqzv \
    --output-path $dada2_output"
echo $cmd >> $logcmd
eval $cmd

cmd="qiime tools export \
    --input-path $tableqza \
    --output-path $dada2_output"
echo $cmd >> $logcmd
eval $cmd

