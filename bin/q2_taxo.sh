#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Taxonomic assignation of ASVs                          ##
##                                                                           ##
###############################################################################
args=("$@")

cpus=${args[0]}
extract_db=${args[1]}
Fprimer=${args[2]}
Rprimer=${args[3]}
confidence=${args[4]}
data_repseqs=${args[5]}
taxoqza=${args[6]}
taxoqzv=${args[7]}
taxo_output=${args[8]}
asv_taxo_tsv=${args[9]}
dbotu3_summary=${args[10]}
final_asv_taxo_biom=${args[11]}
final_asv_taxo_tsv=${args[12]}
database=${args[13]}
seqs_db_filtered=${args[14]}
logcmd=${args[15]}

if [ "$extract_db" = true ]; then
    seqs_db=${args[16]}
    taxa_db=${args[17]}
    #Train the classifier
    cmd="qiime feature-classifier extract-reads \
        --i-sequences $seqs_db \
        --p-f-primer $Fprimer \
        --p-r-primer $Rprimer \
        --o-reads $seqs_db_filtered"
    echo $cmd > $logcmd
    eval $cmd
    
    cmd="qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads $seqs_db_filtered \
        --i-reference-taxonomy $taxa_db \
        --o-classifier $database"
    echo $cmd >> $logcmd
    eval $cmd

else 
    # Keep complete database
    database=${args[16]}
fi

#Run RDP Classifier for taxonomy assignment
cmd="qiime feature-classifier classify-sklearn \
    --p-n-jobs $cpus \
    --p-confidence $confidence \
    --i-classifier $database \
    --i-reads $data_repseqs \
    --o-classification $taxoqza"
echo $cmd >> $logcmd
eval $cmd

#Generate a tabular view of taxonomy metadata
cmd="qiime metadata tabulate \
    --m-input-file $taxoqza 
    --o-visualization $taxoqzv"
echo $cmd >> $logcmd
eval $cmd

#Export data
cmd="qiime tools export \
    --input-path $taxoqzv \
    --output-path $taxo_output"
echo $cmd >> $logcmd
eval $cmd

#Rename tabular taxonomy file and modify header
cmd="mv $taxo_output/metadata.tsv $asv_taxo_tsv"
echo $cmd >> $logcmd
eval $cmd

cmd="sed -i '1,2d' $asv_taxo_tsv"
echo $cmd >> $logcmd
eval $cmd

cmd="sed -i '1 i\\#OTUID\ttaxonomy\tconfidence' $asv_taxo_tsv"
echo $cmd >> $logcmd
eval $cmd

#Add taxonomy to count table (biom format)
cmd="biom add-metadata \
    -i $dbotu3_summary/feature-table.biom \
    --observation-metadata-fp $asv_taxo_tsv \
    -o $final_asv_taxo_biom \
    --sc-separated taxonomy"
echo $cmd >> $logcmd
eval $cmd

#Convert biom table to tabular
cmd="biom convert \
    -i $final_asv_taxo_biom \
    -o $final_asv_taxo_tsv \
    --to-tsv \
    --header-key taxonomy"
echo $cmd >> $logcmd
eval $cmd
