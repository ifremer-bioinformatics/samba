#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Filter contaminants from ASV table                     ##
##                                                                           ##
###############################################################################

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
control_list =  unlist(strsplit(args[[2]],","))
blanks = noquote(args[3])
samples = noquote(args[4])
output_table = args[5]
output_abund_removed = args[6]
output_ASVs_removed = args[7]

# Load ASV table
TABLE = read.table(ASV_table,h=T,sep="\t",check.names=FALSE)
TABLE = TABLE[,colSums(TABLE != 0) >0]
control_list = control_list[control_list %in% colnames(TABLE)]
TABLE_control = TABLE %>% select("ASV_ID",c(control_list), everything())  
blanks = length(control_list)
samples = (length(TABLE_control)-2)-blanks

# microDecon
result = decon(TABLE_control, numb.blanks=as.numeric(blanks), numb.ind=as.numeric(samples), taxa=TRUE, runs=2, thresh=1)

# Output reformat
decontaminated_table = result$decon.table[-2]
final_decontaminated_table = decontaminated_table[rowSums(decontaminated_table[,2:(length(decontaminated_table)-1)])!=0, ]

# Save output
write.table(final_decontaminated_table, file=output_table, row.names=F,col.names=T,quote=F,sep="\t")
write.table(result$reads.removed, file=output_abund_removed,row.names=F,col.names=T,quote=F,sep="\t")
write.table(result$OTUs.removed, file=output_ASVs_removed,row.names=F,col.names=T,quote=F,sep="\t")
