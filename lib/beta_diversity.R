#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Script name: R script for the qiime2_nextflow workflow                 ####
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
## Copyright (c) Cyril NOEL, august-2019                                     ####
## Email: cyril.noel@ifremer.fr                                            ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs the beta diversity based on the   ##
##        non-normalized abundance table                                     ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

#### _____________________________________________________________________ ####

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

betadiversity <- function(PHYLOSEQ, distance, criteria, samples_ordination_plot, metadata) {

    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script: analysis of the beta diversity **           ####
    
    #~~~~~~~~~~~~~~~~~~~~~#
    # Non-normalized data #
    #~~~~~~~~~~~~~~~~~~~~~#
    
    ## /1\ Ordination process ####
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    ord = ordinate(PHYLOSEQ,"NMDS",distance,trymax = 1000)

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]

    group = get_variable(PHYLOSEQ, criteria)
    anosim_result = anosim(distance(PHYLOSEQ,"bray"),group, permutations = 999)
    
    ## /3\ Sample analysis ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, sample ordination plot name, width of graph, heigth of graph
    sample_nmds(PHYLOSEQ, ord, criteria, color_samples, anosim_result, samples_ordination_plot, 12, 10, "NMDS on non-normalized data")
}

main  <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    distance = args[2]
    criteria = args[3]
    samples_ordination_plot = args[4]
    metadata = args[5]
    workflow_dir = args[6]
    if (!exists("sample_nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity(PHYLOSEQ, distance, criteria, samples_ordination_plot, metadata)
}

if (!interactive()) {
        main()
}

