#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Reformat ASV table to use it in statistical analysis   ##
##                                                                           ##
###############################################################################

args=("$@")

metadata=${args[0]}
biom_tsv=${args[1]}
asv_table=${args[2]}
metadata_stats=${args[3]}
microDecon=${args[4]}

cmd="cp $metadata $metadata_stats;
     sed -i '1s/#//' $metadata_stats;
     sed -i '1s/-/_/g' $metadata_stats;
     cp $biom_tsv $asv_table;"
eval $cmd

if [ $microDecon = "true" ]
then
   cmd="sed -i 's/D_0__//g' $asv_table;
   sed -i 's/ D_.__//g' $asv_table"
   eval $cmd
else
   cmd="sed -i '1d' $asv_table; 
   sed -i 's/#OTU ID/ASV_ID/g' $asv_table; 
   sed -i 's/D_0__//g' $asv_table;
   sed -i 's/ D_.__//g' $asv_table"
   eval $cmd
fi
