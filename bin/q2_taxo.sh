#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Script name: q2_taxo.sh                                                 ####
##                                                                           ##
## Purpose of script: Taxonomic assignation of ASVs                          ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2020-01-31                                                 ####
##                                                                           ##
## Email: samba-sebimer@ifremer.fr                                         ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## This program is free software: you can redistribute it and/or modify it   ##
## under the terms of the GNU Affero General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or         ##
## (at your option) any later version.                                       ## 
##                                                                           ##
## License at https://www.gnu.org/licenses/agpl-3.0.txt                      ##
##                                                                           ##
## This program is distributed in the hope that it will be useful, but       ##
## WITHOUT ANY WARRANTY; without even the implied warranty of                ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                      ##
## See the GNU Affero General Public License for more details.               ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Command run by nextflow :
##  q2_taxo.sh ${task.cpus} ${params.taxo.seqs_db} ${params.taxo.taxa_db} ${params.cutadapt.primerF} ${params.cutadapt.primerR} ${params.taxo.confidence} ${data_repseqs} taxonomy.qza taxonomy.qzv taxo_output ASV_taxonomy.tsv ${summary} Final_ASV_table_with_taxonomy.biom Final_ASV_table_with_taxonomy.tsv taxonomic_database.qza seqs_db_amplicons.qza completecmd > q2_taxo.log 2>&1

# Arguments 
args=("$@")

cpus=${args[0]}
seqs_db=${args[1]}
taxa_db=${args[2]}
database_no_extract=${args[3]}
extract_db=${args[4]}
Fprimer=${args[5]}
Rprimer=${args[6]}
confidence=${args[7]}
data_repseqs=${args[8]}
taxoqza=${args[9]}
taxoqzv=${args[10]}
taxo_output=${args[11]}
asv_taxo_tsv=${args[12]}
dbotu3_summary=${args[13]}
final_asv_taxo_biom=${args[14]}
final_asv_taxo_tsv=${args[15]}
database=${args[16]}
seqs_db_filtered=${args[17]}
logcmd=${args[18]}

if [ "$extract_db" = true ]; then
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
    database=$database_no_extract
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
