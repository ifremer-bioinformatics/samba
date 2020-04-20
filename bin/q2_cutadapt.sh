#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Remove primers in raw sequences                        ##
##                                                                           ##
###############################################################################


args=("$@")

singleEnd=${args[0]}
cpus=${args[1]}
imported_data=${args[2]}
primerF=${args[3]}
primerR=${args[4]}
errorRate=${args[5]}
overlap=${args[6]}
data_trimmedqza=${args[7]}
data_trimmedqzv=${args[8]}
trimmed_output=${args[9]}
logcmd=${args[10]}

#Run cutadapt
if ${singleEnd}; then
    cmdoptions="qiime cutadapt trim-single --p-front $primerF"
else
    cmdoptions="qiime cutadapt trim-paired --p-front-f $primerF --p-front-r $primerR" 
fi
cmd="$cmdoptions \
    --verbose \
    --p-cores $cpus \
    --i-demultiplexed-sequences $imported_data \
    --p-error-rate $errorRate \
    --p-discard-untrimmed --p-match-read-wildcards \
    --p-overlap $overlap \
    --o-trimmed-sequences $data_trimmedqza"

echo $cmd > $logcmd
eval $cmd

#Summarize counts per sample for all samples
cmd="qiime demux summarize \
    --verbose \
    --i-data $data_trimmedqza \
    --o-visualization $data_trimmedqzv"
echo $cmd >> $logcmd
eval $cmd

#Export html report
cmd="qiime tools export \
    --input-path $data_trimmedqzv \
    --output-path $trimmed_output"
echo $cmd >> $logcmd
eval $cmd

