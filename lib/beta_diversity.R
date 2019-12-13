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
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2019-12-11                                                 ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## Emails: cyril.noel@ifremer.fr and laure.quintric@ifremer.fr             ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs the beta diversity based on the   ##
##        non-normalized abundance table                                     ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

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
library("dendextend")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#										#
# Function to standardize ordination analysis and plot for each distance matrix #
#										#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

 
betadiversity <- function(PHYLOSEQ, distance, metadata, criteria, nmds, pcoa, method_hc, plot_hc){
    
    #~~~~~~~~~~~~~~~~~~~~~#
    # Non-normalized data #
    #~~~~~~~~~~~~~~~~~~~~~#
    
    ## Ordination process ####
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    
    if (distance == "jaccard" | distance == "bray") {
    ord_nmds = ordinate(PHYLOSEQ,"NMDS",distance,trymax = 1000)
    ord_pcoa = ordinate(PHYLOSEQ,"PCoA",distance,trymax = 1000)
    }
    else {
    ord_nmds = ordinate(PHYLOSEQ,"NMDS",distance)
    ord_pcoa = ordinate(PHYLOSEQ,"PCoA",distance)
    }

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]

    group = get_variable(PHYLOSEQ, criteria)
    anosim_result = anosim(distance(PHYLOSEQ,distance),group, permutations = 999)

    ## Sample analysis ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, ordination plot name, width of graph, heigth of graph
    plot.nmds(PHYLOSEQ, ord_nmds, criteria, color_samples, anosim_result, nmds, distance, 12, 10, paste("NMDS on non-normalized data","based on",distance,"distance",sep=" "))
    plot.pcoa(PHYLOSEQ, ord_pcoa, criteria, color_samples, anosim_result, pcoa, distance, 12, 10, paste("MDS-PCoA on non-normalized data","based on",distance,"distance",sep=" "))
    
    ## Hierarchical clsutering ####    
    dist = distance(PHYLOSEQ, distance, type="samples")
    hc = hclust(dist, method = method_hc)
    dendro = as.dendrogram(hc)
    group = data.frame(PHYLOSEQ@sam_data[,criteria])[,1]
    n_group = length(unique(group))
    cols = color_vector[1:n_group]
    col_group = cols[group]
    plot.hc(dendro, group, cols, col_group, method_hc, plot_hc, distance, 12, 10)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#                                             #
# Ordination process for each distance matrix #
#                                             #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

main_jaccard  <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    distance = "jaccard"
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[2], "-", "_")
    metadata = args[3]
    workflow_dir = args[4]
    nmds = args[5]
    pcoa = args[6]
    method_hc = args[7]
    plot_hc = args[8]
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity(PHYLOSEQ, distance, metadata, criteria, nmds, pcoa, method_hc, plot_hc)
}

if (!interactive()) {
        main_jaccard()
}


main_bray  <- function(){
    # Get arguments from RScript command lin
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    distance = "bray"
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[2], "-", "_")
    metadata = args[3]
    workflow_dir = args[4]
    nmds = args[5]
    pcoa = args[6]
    method_hc = args[7]
    plot_hc = args[8]
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity(PHYLOSEQ, distance, metadata, criteria, nmds, pcoa, method_hc, plot_hc)
}

if (!interactive()) {
        main_bray()
}

main_unifrac  <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    distance = "unifrac"
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[2], "-", "_")
    metadata = args[3]
    workflow_dir = args[4]
    nmds = args[5]
    pcoa = args[6]
    method_hc = args[7]
    plot_hc = args[8]
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity(PHYLOSEQ, distance, metadata, criteria, nmds, pcoa, method_hc, plot_hc)
}

if (!interactive()) {
        main_unifrac()
}

main_wunifrac  <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    distance = "wunifrac"
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[2], "-", "_")
    metadata = args[3]
    workflow_dir = args[4]
    nmds = args[5]
    pcoa = args[6]
    method_hc = args[7]
    plot_hc = args[8]
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity(PHYLOSEQ, distance, metadata, criteria, nmds, pcoa, method_hc, plot_hc)
}

if (!interactive()) {
        main_wunifrac()
}
