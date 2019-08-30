#!/usr/bin/env bash
## Command run by nextflow process :
### ${baseDir}/lib/q2_taxo.sh ${task.cpus} ${params.taxo.confidence} ${params.taxo.database} ${data_repseqs} taxonomy.qza taxonomy.qzv taxo_output ASV_taxonomy.tsv ${dada2_summary} Final_ASV_table_with_taxonomy.biom Final_ASV_table_with_taxonomy.tsv completecmd > q2_taxo.log 2>&1
# store arguments in a special array 
args=("$@")

cpus=${args[0]}
confidence=${args[1]}
database=${args[2]}
data_repseqs=${args[3]}
taxoqza=${args[4]}
taxoqzv=${args[5]}
taxo_output=${args[6]}
asv_taxo_tsv=${args[7]}
dada2_summary=${args[8]}
final_asv_taxo_biom=${args[9]}
final_asv_taxo_tsv=${args[10]}
logcmd=${args[11]}

#Run RDP Classifier for taxonomy assignment
cmd="qiime feature-classifier classify-sklearn \
    --p-n-jobs $cpus \
    --p-confidence $confidence \
    --i-classifier $database \
    --i-reads $data_repseqs \
    --o-classification $taxoqza"
echo $cmd > $logcmd
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
