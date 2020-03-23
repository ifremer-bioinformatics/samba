#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Script name: prepare_data_for_stats.sh                                  ####
##                                                                           ##
## Purpose of script: Reformat ASV table to use it in statistical analysis   ##
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
## prepare_data_for_stats.sh ${metadata} ${biom_tsv} ASV_table_with_taxo_for_stats q2_metadata ${params.microDecon_enable} > prepare_data_for_stats.log 2&>1

# Arguments 
args=("$@")

metadata=${args[0]}
biom_tsv=${args[1]}
asv_table=${args[2]}
metadata_stats=${args[3]}
microDecon=${args[4]}

if [ $microDecon = "false" ]
then
   cmd="cp $metadata $metadata_stats; 
   sed -i '1s/#//' $metadata_stats;
   sed -i '1s/-/_/g' $metadata_stats;
   cp $biom_tsv $asv_table; 
   sed -i '1d' $asv_table; 
   sed -i 's/#OTU ID/ASV_ID/g' $asv_table; 
   sed -i 's/D_0__//g' $asv_table;
   sed -i 's/ D_.__//g' $asv_table"
   eval $cmd
else
   cmd="cp $metadata $metadata_stats; 
   sed -i '1s/#//' $metadata_stats;
   sed -i '1s/-/_/g' $metadata_stats;
   cp $biom_tsv $asv_table; 
   sed -i 's/D_0__//g' $asv_table;
   sed -i 's/ D_.__//g' $asv_table"
   eval $cmd
fi
