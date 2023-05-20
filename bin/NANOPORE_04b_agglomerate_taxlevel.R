#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: agglomerate the phyloseq object at different taxonomic ##
##                    levels                                                 ##
##                                                                           ##
###############################################################################

## Load up the packages needed ####
requiredPackages = c("phyloseq", "speedyseq", "dplyr", "stringr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

agglomerate_phyloseq <- function(phyloseq_rds_raw_all, phyloseq_rds_raw_assigned, tax_level) {
  
  # Input data
  RAW_PHYLOSEQ_ALL = readRDS(phyloseq_rds_raw_all)
  RAW_PHYLOSEQ_ASSIGNED = readRDS(phyloseq_rds_raw_assigned)
  level = tax_level
  LEVEL = toupper(level)

  # Agglomerate data at the defined taxonomic level for ALL ASSIGNATION
  assign(paste("PHYLOSEQ_", LEVEL, "_ALL_ASSIGNATION", sep=""), tax_glom(RAW_PHYLOSEQ_ALL, level))
  TABLE_ALL = data.frame(otu_table(eval(parse(text = paste("PHYLOSEQ_", LEVEL, "_ALL_ASSIGNATION", sep="")))), check.names=FALSE)
  TABLE_ALL$taxonomy = data.frame(tax_table(eval(parse(text=paste("PHYLOSEQ_", LEVEL, "_ALL_ASSIGNATION", sep="")))))[,level]
  TABLE_ALL = data.frame(TABLE_ALL %>% group_by(taxonomy) %>% summarise(across(where(is.numeric), list(sum=sum), na.rm=TRUE)))
  rownames(TABLE_ALL) = TABLE_ALL$taxonomy
  TABLE_ALL = TABLE_ALL %>% select(-taxonomy)
  colnames(TABLE_ALL) = str_replace_all(colnames(TABLE_ALL), "_sum", "")
  write.table(TABLE_ALL, paste("asv_table_all_assignation_", level, ".tsv", sep=""), col.names=T, row.names=T, sep="\t", quote=F)
  saveRDS(eval(parse(text=paste("PHYLOSEQ_", LEVEL, "_ALL_ASSIGNATION", sep=""))), file=paste("phyloseq_all_assignation_", level, ".rds", sep=""))

  # Agglomerate data at the defined taxonomic level for ASSIGNED
  assign(paste("PHYLOSEQ_", LEVEL, "_ONLY_ASSIGNED", sep=""), tax_glom(RAW_PHYLOSEQ_ASSIGNED, level))
  TABLE_ASSIGNED = data.frame(otu_table(eval(parse(text=paste("PHYLOSEQ_", LEVEL, "_ONLY_ASSIGNED", sep="")))), check.names=FALSE)
  TABLE_ASSIGNED$taxonomy = data.frame(tax_table(eval(parse(text=paste("PHYLOSEQ_", LEVEL, "_ONLY_ASSIGNED", sep="")))))[,level]
  TABLE_ASSIGNED = data.frame(TABLE_ASSIGNED %>% group_by(taxonomy) %>% summarise(across(where(is.numeric), list(sum=sum), na.rm=TRUE)))
  rownames(TABLE_ASSIGNED) = TABLE_ASSIGNED$taxonomy
  TABLE_ASSIGNED = TABLE_ASSIGNED %>% select(-taxonomy)
  colnames(TABLE_ASSIGNED) = str_replace_all(colnames(TABLE_ASSIGNED), "_sum", "")
  write.table(TABLE_ASSIGNED, paste("asv_table_only_assigned_", level, ".tsv", sep=""), col.names=T, row.names=T, sep="\t", quote=F)
  saveRDS(eval(parse(text=paste("PHYLOSEQ_", LEVEL, "_ONLY_ASSIGNED", sep=""))), file=paste("phyloseq_only_assigned_", level, ".rds", sep=""))
  
}

main <- function() {
  # Get arguments from RScript command line
  args = commandArgs(trailingOnly=TRUE)
  phyloseq_rds_all = args[1]
  phyloseq_rds_assigned = args[2]
  tax_level = args[3]
  agglomerate_phyloseq(phyloseq_rds_all, phyloseq_rds_assigned, tax_level)
}

if (!interactive()) {
  main()
}
