#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Clustering of AVSs to reduce PCR errors                ##
##                                                                           ##
###############################################################################

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
gen_crit=${args[9]}
abund_crit=${args[10]}
pval_crit=${args[11]}
logcmd=${args[12]}

# Run dbotu2
cmd="qiime dbotu-q2 call-otus --verbose --i-table $table --i-sequences $seqs --p-gen-crit $gen_crit --p-abund-crit $abund_crit --p-pval-crit $pval_crit --o-representative-sequences $dbotu3_seqsqza --o-dbotu-table $dbotu3_tableqza > $dbotu3_details"
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
