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
tax_to_include=${args[3]}
filtered_table=${args[4]}
asv_seq=${args[5]}
filtered_seq=${args[6]}
filtered_table_qzv=${args[7]}
metadata=${args[8]}
filtered_tax_output=${args[9]}
filtered_seq_qzv=${args[10]}
asv_taxo_tsv=${args[11]}
final_filtered_table_biom=${args[12]}
final_filtered_table_tsv=${args[13]}
logcmd=${args[14]}

if [$tax_to_exclude != "none"]
then
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
else
    #Filtering ASV table based on taxonomy
    cmd="qiime taxa filter-table \
      --i-table $asv_table \
      --i-taxonomy $asv_tax \
      --p-include $tax_to_include \
      --o-filtered-table $filtered_table"
    echo $cmd >> $logcmd
    eval $cmd
    if [ ! $? -eq 0 ]; then exit 1; fi

    #Filtering ASV sequences based on taxonomy
    cmd="qiime taxa filter-seqs \
      --i-sequences $asv_seq \
      --i-taxonomy $asv_tax \
      --p-include $tax_to_include \
      --o-filtered-sequences $filtered_seq"
    echo $cmd >> $logcmd
    eval $cmd
    if [ ! $? -eq 0 ]; then exit 1; fi
fi

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
