#!/usr/bin/env bash
## Command run by nextflow process :
### ${baseDir}/lib/prepare_data_for_stats.sh ${metadata} ${biom_tsv} ASV_table_with_taxo_for_stats q2_metadata complete_cmd_prepare_stats > prepare_data_for_stats.log 2&>1
# store arguments in a special array 
args=("$@")

metadata=${args[0]}
biom_tsv=${args[1]}
asv_table=${args[2]}
metadata_stats=${args[3]}
logcmd=${args[4]}

cmd="cp $metadata $metadata_stats; 
sed '0,/#/s/#//' $metadata_stats;
sed -i '1s/-/_/g' $metadata_stats;
cp $biom_tsv $asv_table; 
sed -i '1d' $asv_table; 
sed -i 's/#OTU ID/ASV_ID/g' $asv_table; 
sed -i 's/D_.__//g' $asv_table"

#sed -i 's/#SampleID/SampleID/g' $metadata_stats; 
echo $cmd > $logcmd
eval $cmd
