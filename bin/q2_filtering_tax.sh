#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Filtering ASV table & seq based on taxonomy            ##
##                                                                           ##
###############################################################################
args=("$@")
asv_table=${args[0]}
asv_tax=${args[1]}
tax_to_exclude=${args[2]}
filtered_table=${args[3]}
asv_seq=${args[4]}
filtered_seq=${args[5]}
filtered_table_qzv=${args[6]}
metadata=${args[7]}
filtered_tax_output=${args[8]}
filtered_seq_qzv=${args[9]}
asv_taxo_tsv=${args[10]}
final_filtered_table_biom=${args[11]}
final_filtered_table_tsv=${args[12]}
logcmd=${args[13]}

#Filtering ASV table based on taxonomy
cmd="qiime taxa filter-table \
  --i-table $asv_table \
  --i-taxonomy $asv_tax \
  --p-exclude $tax_to_exclude \
  --o-filtered-table $filtered_table"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi

#Filtering ASV sequences based on taxonomy
cmd="qiime taxa filter-seqs \
  --i-sequences $asv_seq \
  --i-taxonomy $asv_tax \
  --p-exclude $tax_to_exclude \
  --o-filtered-sequences $filtered_seq"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi

#Export data
cmd="qiime feature-table summarize \
    --i-table $filtered_table \
    --o-visualization $filtered_table_qzv \
    --m-sample-metadata-file $metadata"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi

cmd="qiime tools export \
    --input-path $filtered_table_qzv \
    --output-path $filtered_tax_output"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi

cmd="qiime feature-table tabulate-seqs \
    --i-data $filtered_seq \
    --o-visualization $filtered_seq_qzv"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi

cmd="qiime tools export \
    --input-path $filtered_seq_qzv \
    --output-path $filtered_tax_output"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi

cmd="qiime tools export \
    --input-path $filtered_table \
    --output-path $filtered_tax_output"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi

#Add taxonomy to count table (biom format)
cmd="biom add-metadata \
    -i $filtered_tax_output/feature-table.biom \
    --observation-metadata-fp $asv_taxo_tsv \
    -o $final_filtered_table_biom \
    --sc-separated taxonomy"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi

#Convert biom table to tabular
cmd="biom convert \
    -i $final_filtered_table_biom \
    -o $final_filtered_table_tsv \
    --to-tsv \
    --header-key taxonomy"
echo $cmd >> $logcmd
eval $cmd
if [ ! $? -eq 0 ]; then exit 1; fi
