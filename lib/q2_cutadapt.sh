#!/usr/bin/env bash
## Command run by nextflow process :
### ${baseDir}/lib/q2_cutadapt.sh ${task.cpus} ${imported_data} ${params.cutadapt.primerF} ${params.cutadapt.primerR} \ 
###                               ${params.cutadapt.errorRate} ${params.cutadapt.overlap} data_trimmed.qza data_trimmed.qzv \
###                               trimmed_output completecmd > q2_cutadapt.log 2>&1
# store arguments in a special array 
args=("$@")

cpus=${args[0]}
imported_data=${args[1]}
primerF=${args[2]}
primerR=${args[3]}
errorRate=${args[4]}
overlap=${args[5]}
data_trimmedqza=${args[6]}
data_trimmedqzv=${args[7]}
trimmed_output=${args[8]}
logcmd=${args[9]}

#Run cutadapt
cmd="qiime cutadapt trim-paired \
    --verbose \
    --p-cores $cpus \
    --i-demultiplexed-sequences $imported_data \
    --p-front-f $primerF \
    --p-front-r $primerR \
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

