#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Summarise ANCOM-BC results                             ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("stringr","dplyr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

args = commandArgs(trailingOnly=TRUE)

var_name = args[1]
lfc_file = args[2]
qvalue_file = args[3]
ancombc_summary_output= args[4]

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#                                #
# Merge ANCOM-BC multiple result #
# into a unique file             #
#                                #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

## Import ANCOM-BC results ####
lfc = read.csv(lfc_file, check.names=F)
colnames(lfc) = str_remove_all(colnames(lfc), var_name) 
lfc = lfc[, c(2,4:length(lfc))]
colnames(lfc) = c("taxon", paste(colnames(lfc)[2:length(lfc)],"lfc",sep="_"))

padjvalue = read.csv(qvalue_file, check.names=F)
colnames(padjvalue) = str_remove_all(colnames(padjvalue), var_name)
padjvalue = padjvalue[, c(2,4:length(padjvalue))]
padjvalue_signif = padjvalue
colnames(padjvalue) = c("taxon", paste(colnames(padjvalue)[2:length(padjvalue)],"padjvalue",sep="_"))
rownames(padjvalue_signif) = padjvalue_signif$taxon
padjvalue_signif = padjvalue_signif %>% select(-taxon)
padjvalue_signif[padjvalue_signif <= 0.001] <- "***"
padjvalue_signif[(padjvalue_signif < 0.01) & (padjvalue_signif > 0.001)] <- "**"
padjvalue_signif[(padjvalue_signif < 0.05) & (padjvalue_signif > 0.01)] <- "*"
padjvalue_signif[padjvalue_signif > 0.05] <- "ns"
colnames(padjvalue_signif) = paste(colnames(padjvalue_signif), "signif", sep="_")
padjvalue = merge(padjvalue, padjvalue_signif, by.x="taxon", by.y="row.names")

## Merge files ####
merged_ancombc_results = merge(lfc, padjvalue, by="taxon")
merged_ancombc_results = merged_ancombc_results %>% select(taxon, order(colnames(merged_ancombc_results)))

## Save output ####
write.table(merged_ancombc_results, ancombc_summary_output, col.names=T, row.names=F, sep="\t", quote=F)
