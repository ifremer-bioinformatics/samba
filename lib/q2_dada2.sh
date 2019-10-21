#!/usr/bin/env bash
## Command run by nextflow process :
### ${baseDir}/lib/q2_dada2_2.sh ${trimmed_data} ${metadata} rep_seqs.qza rep_seqs.qzv table.qza table.qzv stats.qza stats.qzv dada2_output ${params.dada2.trim3F} ${params.dada2.trim3R} ${params.dada2.trunclenF} ${params.dada2.trunclenR} ${params.dada2.maxee_f} ${params.dada2.maxee_r} ${params.dada2.minqual} ${params.dada2.chimeras} ${task.cpus} completecmd > q2_dada2.log 2>&1
# store arguments in a special array 
args=("$@") 

trimmed_data=${args[0]}
metadata=${args[1]}
rep_seqsqza=${args[2]}
rep_seqsqzv=${args[3]}
tableqza=${args[4]}
tableqzv=${args[5]}
statsqza=${args[6]}
statsqzv=${args[7]}
dada2_output=${args[8]}
trim3F=${args[9]}
trim3R=${args[10]}
trunclenF=${args[11]}
trunclenR=${args[12]}
maxee_f=${args[13]}
maxee_r=${args[14]}
minqual=${args[15]}
chimeras=${args[16]}
cpus=${args[17]}
logcmd=${args[18]}

#Run dada2 : denoises paired-end sequences, dereplicates them and filters chimeras
cmd="qiime dada2 denoise-paired \
    --verbose \
    --i-demultiplexed-seqs $trimmed_data \
    --p-trim-left-f $trim3F \
    --p-trim-left-r $trim3R \
    --p-trunc-len-f $trunclenF \
    --p-trunc-len-r $trunclenR \
    --p-max-ee-f $maxee_f \
    --p-max-ee-r $maxee_r \
    --p-trunc-q $minqual \
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

