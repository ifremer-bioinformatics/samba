#!/usr/bin/env Rscript

###############################################################################
##                                                                           ##
## Script name: microDecon.R                                               ####
##                                                                           ##
## Purpose of script: Filter contaminants from ASV table                     ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2020-02-28                                               ####
## Modified on: 2020-03-05                                                 ####
##                                                                           ##
## Email: samba-sebimer@ifremer.fr                                         ####
##                                                                           ##
## Copyright (c) SeBiMER, february-2020                                    ####
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
##                                                                           ##
## Notes: This part of the script performs decontamination of your data      ##
##        based on control samples using microDecon R package                ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Load up the needed packages ####
requiredPackages = c("microDecon","dplyr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

# @@@@@@@@@@@@@@@@@@@@@@@@ #
#                          #
# Decontamination analysis #
#                          #
# @@@@@@@@@@@@@@@@@@@@@@@@ #
args = commandArgs(trailingOnly=TRUE)

ASV_table = args[1]
control-list =  unlist(strsplit(args[[2]],","))
blanks = noquote(args[3])
samples = noquote(args[4])
output_table = args[5]
output_abund_removed = args[6]
output_ASVs_removed = args[7]

# Load ASV table
TABLE = read.table(ASV_table,h=T,sep="\t",check.names=FALSE)
TABLE_control = TABLE %>% select("ASV_ID",c(control-list), everything())  
head(TABLE_control)

# microDecon
result = decon(TABLE_control, numb.blanks=as.numeric(blanks), numb.ind=as.numeric(samples), taxa=TRUE, runs=2, thresh=1)

# Output reformat
decontaminated_table = result$decon.table[-2]
final_decontaminated_table = decontaminated_table[rowSums(decontaminated_table[,2:(length(decontaminated_table)-1)])!=0, ]

# Save output
write.table(final_decontaminated_table, file=output_table, row.names=F,col.names=T,quote=F,sep="\t")
write.table(result$reads.removed, file=output_abund_removed,row.names=F,col.names=T,quote=F,sep="\t")
write.table(result$OTUs.removed, file=output_ASVs_removed,row.names=F,col.names=T,quote=F,sep="\t")
