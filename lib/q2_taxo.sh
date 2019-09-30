#!/usr/bin/env bash
## Command run by nextflow process :
###  ${baseDir}/lib/q2_taxo.sh ${task.cpus} ${params.taxo.db_seqs} ${params.taxo.db_tax} ${params.cutadapt.primerF} ${params.cutadapt.primerR} ${params.taxo.confidence} ${data_repseqs} taxonomy.qza taxonomy.qzv taxo_output ASV_taxonomy.tsv ${dada2_summary} Final_ASV_table_with_taxonomy.biom Final_ASV_table_with_taxonomy.tsv taxonomic_database.qza db_seqs_amplicons.qza completecmd > q2_taxo.log 2>&1
# store arguments in a special array 
args=("$@")

cpus=${args[0]}
db_seqs=${args[1]}
db_tax=${args[2]}
Fprimer=${args[3]}
Rprimer=${args[4]}
confidence=${args[5]}
data_repseqs=${args[6]}
taxoqza=${args[7]}
taxoqzv=${args[8]}
taxo_output=${args[9]}
asv_taxo_tsv=${args[10]}
dada2_summary=${args[11]}
final_asv_taxo_biom=${args[12]}
final_asv_taxo_tsv=${args[13]}
database=${args[14]}
db_seqs_filtered=${args[15]}
logcmd=${args[16]}

#Train the classifier
cmd="qiime feature-classifier extract-reads \
    --i-sequences $db_seqs \
    --p-f-primer $Fprimer \
    --p-r-primer $Rprimer \
    --o-reads $db_seqs_filtered"
echo $cmd > $logcmd
eval $cmd

cmd="qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads $db_seqs_filtered \
    --i-reference-taxonomy $db_tax \
    --o-classifier $database"
echo $cmd >> $logcmd
eval $cmd

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
    -i $dada2_summary/feature-table.biom \
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
