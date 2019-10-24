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
library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("svglite")
library("tidyr")
library("gridExtra")
library("egg")
library("metagenomeSeq")
library("vegan")
library("stringr")

betadiversity_css <- function (PHYLOSEQ, final_css_ASV_table_with_taxonomy, distance, criteria, samples_ordination_plot_css, metadata) {

    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script: analysis of the beta diversity **           ####
    
    #~~~~~~~~~~~~~~~~~~~#
    # CSS normalization #
    #~~~~~~~~~~~~~~~~~~~#
    css = phyloseq_to_metagenomeSeq(PHYLOSEQ)
    p = cumNormStatFast(css)
    css = cumNorm(css, p=p)
    css_norm_factor = normFactors(css)          
    CSS_TABLE = MRcounts(css, norm = T)
    PHYLOSEQ_css = PHYLOSEQ
    otu_table(PHYLOSEQ_css) = otu_table(CSS_TABLE,taxa_are_rows=TRUE)
    CSS_normalized_table = cbind(as.data.frame(otu_table(PHYLOSEQ_css)),as.data.frame(tax_table(PHYLOSEQ_css)))
    write.table(CSS_normalized_table,final_css_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=",") 
    
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    ord_css = ordinate(PHYLOSEQ_css,"NMDS", distance, trymax = 1000)

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]
    
    group_css = get_variable(PHYLOSEQ_css, criteria)
    anosim_result_css = anosim(distance(PHYLOSEQ_css,"bray"),group_css, permutations = 999)
    
    ### /3\ Sample analysis ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, sample ordination plot name, width of graph, heigth of graph, graph title
    sample_nmds(PHYLOSEQ_css, ord_css, criteria, color_samples, anosim_result_css, samples_ordination_plot_css, 12, 10, "NMDS on CSS normalized data")
}

main <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    final_css_ASV_table_with_taxonomy = args[2]
    distance = args[3]
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[4], "-", "_")
    samples_ordination_plot_css = args[5]
    metadata = args[6]
    workflow_dir = args[7]
    if (!exists("sample_nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity_css(PHYLOSEQ, final_css_ASV_table_with_taxonomy, distance, criteria, samples_ordination_plot_css, metadata)
}

if (!interactive()) {
        main()
}

