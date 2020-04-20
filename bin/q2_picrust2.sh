#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Make functional predictions from metabarcoding data    ##
##                                                                           ##
###############################################################################

args=("$@")

table=${args[0]}
seqs=${args[1]}
output=${args[2]}
ncpus=${args[3]}
method=${args[4]}
nsti=${args[5]}
logcmd=${args[6]}

#Import all samples paired-end files listed in manifest to qiime data structure
cmd="qiime picrust2 full-pipeline --i-table $table --i-seq $seqs --output-dir $output --p-threads $ncpus --p-hsp-method $method --p-max-nsti $nsti --p-highly-verbose --verbose"
echo $cmd > $logcmd
eval $cmd

#Get summary information regarding the prediction files
cmd="qiime feature-table summarize --i-table $output/pathway_abundance.qza --o-visualization $output/pathway_abundance.qzv ;
qiime tools export --input-path $output/pathway_abundance.qzv --output-path $output/pathway_abundance_visu ;
qiime feature-table summarize --i-table $output/ec_metagenome.qza --o-visualization $output/ec_metagenome.qzv ;
qiime tools export --input-path $output/ec_metagenome.qzv --output-path $output/ec_metagenome_visu ;
qiime feature-table summarize --i-table $output/ko_metagenome.qza --o-visualization $output/ko_metagenome.qzv ;
qiime tools export --input-path $output/ko_metagenome.qzv --output-path $output/ko_metagenome_visu ;"
echo $cmd >> $logcmd
eval $cmd

#Get summary of predictions for stats
cmd="qiime tools export --input-path $output/pathway_abundance.qza --output-path $output/pathway_abundance_exported ;
biom convert -i $output/pathway_abundance_exported/feature-table.biom -o $output/pathway_abundance_exported/pathway_abundance_predictions.tsv --to-tsv ;
qiime tools export --input-path $output/ec_metagenome.qza --output-path $output/ec_metagenome_exported ;
biom convert -i $output/ec_metagenome_exported/feature-table.biom -o $output/ec_metagenome_exported/ec_metagenome_predictions.tsv --to-tsv ;
qiime tools export --input-path $output/ko_metagenome.qza --output-path $output/ko_metagenome_exported ;
biom convert -i $output/ko_metagenome_exported/feature-table.biom -o $output/ko_metagenome_exported/ko_metagenome_predictions.tsv --to-tsv"
echo $cmd >> $logcmd
eval $cmd

#Add descriptions EC
cmd="sed -i '1d' $output/ec_metagenome_exported/ec_metagenome_predictions.tsv ;
sed -i 's/#OTU ID/EC_ID/g' $output/ec_metagenome_exported/ec_metagenome_predictions.tsv ;
add_descriptions.py -i $output/ec_metagenome_exported/ec_metagenome_predictions.tsv -m EC -o $output/ec_metagenome_exported/ec_metagenome_predictions_with-descriptions.tsv "
echo $cmd >> $logcmd
eval $cmd

#Add descriptions KO
cmd="sed -i '1d' $output/ko_metagenome_exported/ko_metagenome_predictions.tsv ;
sed -i 's/#OTU ID/KO_ID/g' $output/ko_metagenome_exported/ko_metagenome_predictions.tsv ;
add_descriptions.py -i $output/ko_metagenome_exported/ko_metagenome_predictions.tsv -m KO -o $output/ko_metagenome_exported/ko_metagenome_predictions_with-descriptions.tsv "
echo $cmd >> $logcmd
eval $cmd

#Add descriptions METACYC
cmd="sed -i '1d' $output/pathway_abundance_exported/pathway_abundance_predictions.tsv ;
sed -i 's/#OTU ID/METACYC_ID/g' $output/pathway_abundance_exported/pathway_abundance_predictions.tsv ;
add_descriptions.py -i $output/pathway_abundance_exported/pathway_abundance_predictions.tsv -m METACYC -o $output/pathway_abundance_exported/pathway_abundance_predictions_with-descriptions.tsv"
echo $cmd >> $logcmd
eval $cmd
