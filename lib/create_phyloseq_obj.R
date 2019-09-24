#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Script name: R script for the qiime2_snakemake workflow                 ####
##                                                                           ##
## Purpose of script: Automated Statistical Analyzes of Metabarcoding Data   ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Author: Dr. Cyril NOEL                                                  ####
##         Bioinformatics engineer                                           ##
##         SeBiMER, Ifremer                                                  ##
##                                                                           ##
## Date Created: 2019-08-29                                                ####
##                                                                           ##
## Copyright (c) Cyril NOEL, aout-2019                                     ####
## Email: cyril.noel@ifremer.fr                                            ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs only the creation of a phyloseq   ##
##        object                                                             ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


#### _____________________________________________________________________ ####

## Load up the packages needed ####
library("dplyr")
library("stringr")
library("phyloseq")
#library("ggplot2")
#library("RColorBrewer")
#library("svglite")
#library("tidyr")
#library("gridExtra")
#library("egg")

create_phyloseq_obj <- function(phyloseq_rds, biom_tsv, metadata) {
    #Input data
    rawASVtable = read.table(biom_tsv, h=T, sep="\t", dec=".", check.names=FALSE)
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    
    #Reformatting of the input data
    abund = rawASVtable[,1:length(rawASVtable)-1]
    row.names(abund) = abund$ASV_ID
    abund = abund %>% select (-ASV_ID)
    tax = rawASVtable[,c(1,length(names(rawASVtable)))]
    tax = data.frame(tax$ASV_ID,do.call(rbind,list(str_split_fixed(tax$taxonomy, ";",7))))
    colnames(tax) = c("ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    row.names(tax) = tax$ASV_ID
    tax = as.matrix(tax %>% select (-ASV_ID))
    
    ## Construction of the phyloseq object ####
    ABUND = otu_table(abund,taxa_are_rows=TRUE)
    TAX = tax_table(tax)
    METADATA = sample_data(metadata)
    PHYLOSEQ = phyloseq(ABUND,TAX,METADATA)
    saveRDS(PHYLOSEQ, file=phyloseq_rds)
}

main <- function() {
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    phyloseq_rds = args[1]
    biom_tsv = args[2]
    metadata = args[3]
    create_phyloseq_obj(phyloseq_rds, biom_tsv, metadata)
}

if (!interactive()) {
        main()
}
