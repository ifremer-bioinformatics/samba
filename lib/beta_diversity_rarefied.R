#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Script name: R script for the SAMBA-nextflow workflow                   ####
##                                                                           ##
## Purpose of script: Automated Statistical Analyses of Metabarcoding Data   ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Cyril NOEL and Laure QUINTRIC                                  ####
##         Bioinformatics engineers                                          ##
##         SeBiMER, Ifremer                                                  ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2019-10-23                                                 ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## Emails: cyril.noel@ifremer.fr and laure.quintric@ifremer.fr             ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs the beta diversity based on the   ##
##        rarefied abundance table                                           ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Load up the packages needed ####
library("dplyr")
library("stringr")
library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("svglite")
library("tidyr")
library("gridExtra")
library("egg")
library("vegan")
library("stringr")

betadiversity_rarefied <- function (PHYLOSEQ, final_rarefied_ASV_table_with_taxonomy, distance, criteria, samples_ordination_plot_rarefied, metadata) {

    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script: analysis of the beta diversity **           ####
    
    #~~~~~~~~~~~~~~~#
    # Rarefied data #
    #~~~~~~~~~~~~~~~#
    PHYLOSEQ_rarefied = rarefy_even_depth(PHYLOSEQ, sample.size=min(sample_sums(PHYLOSEQ)),rngseed=1000,replace=TRUE,trimOTUs=TRUE)
    rarefied_table = cbind(as.data.frame(otu_table(PHYLOSEQ_rarefied)),as.data.frame(tax_table(PHYLOSEQ_rarefied)))
    write.table(rarefied_table,final_rarefied_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=".",quote=F) 

    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    ord_rarefied = ordinate(PHYLOSEQ_rarefied, "NMDS", distance, trymax = 1000)

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]
    
    group_rarefied = get_variable(PHYLOSEQ_rarefied, criteria)
    anosim_result_rarefied = anosim(distance(PHYLOSEQ_rarefied,"bray"),group_rarefied, permutations = 999)

    ## /3\ Sample analysis ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, sample ordination plot name, width of graph, heigth of graph, graph title
    sample_nmds(PHYLOSEQ_rarefied, ord_rarefied, criteria, color_samples, anosim_result_rarefied, samples_ordination_plot_rarefied, 12, 10, "NMDS on rarefied data")
}

main <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1]) 
    final_rarefied_ASV_table_with_taxonomy = args[2]
    distance = args[3]
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[4], "-", "_")
    samples_ordination_plot_rarefied = args[5]
    metadata = args[6]
    workflow_dir = args[7]
    if (!exists("sample_nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity_rarefied(PHYLOSEQ, final_rarefied_ASV_table_with_taxonomy, distance, criteria, samples_ordination_plot_rarefied, metadata)
}

if (!interactive()) {
        main()
}
