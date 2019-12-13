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
## Modified on: 2019-12-13                                                 ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## Emails: cyril.noel@ifremer.fr and laure.quintric@ifremer.fr             ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs the beta diversity (NMDS, PCoa &  ##
##        Hierachical Clustering) on the non-normalized ASV table based on   ##
##        four distance matrices (Jaccard, Bray-Curtis, UniFrac & Weighted   ##
##        UniFrac) 							     ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

## Install (if necessary) and load up the needed packages ####
requiredPackages_CRAN = c("dplyr","stringr","ggplot2","RColorBrewer","svglite","tidyr","gridExtra","egg","vegan","dendextend","BiocManager")
for(package in requiredPackages_CRAN){
  if(!require(package,character.only = TRUE)) install.packages(package)
  library(package,character.only = TRUE)
}

requiredPackages_BIOCONDUCTOR = c("phyloseq")
for(package in requiredPackages_BIOCONDUCTOR){
  if(!require(package,character.only = TRUE)) BiocManager::install(package)
  library(package,character.only = TRUE)
}

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

    ## Ordination plots ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, ordination plot name, distance, width of graph, heigth of graph
    plot.nmds(PHYLOSEQ, ord_nmds, criteria, color_samples, anosim_result, nmds, distance, 12, 10, paste("NMDS on non-normalized data","based on",distance,"distance",sep=" "))
    plot.pcoa(PHYLOSEQ, ord_pcoa, criteria, color_samples, anosim_result, pcoa, distance, 12, 10, paste("MDS-PCoA on non-normalized data","based on",distance,"distance",sep=" "))
    
    ## Hierarchical clustering ####    
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
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[2], "-", "_")
    metadata = args[3]
    workflow_dir = args[4]
    nmds = args[5]
    pcoa = args[6]
    method_hc = args[7]
    plot_hc = args[8]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
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
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[2], "-", "_")
    metadata = args[3]
    workflow_dir = args[4]
    nmds = args[5]
    pcoa = args[6]
    method_hc = args[7]
    plot_hc = args[8]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
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
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[2], "-", "_")
    metadata = args[3]
    workflow_dir = args[4]
    nmds = args[5]
    pcoa = args[6]
    method_hc = args[7]
    plot_hc = args[8]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
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
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[2], "-", "_")
    metadata = args[3]
    workflow_dir = args[4]
    nmds = args[5]
    pcoa = args[6]
    method_hc = args[7]
    plot_hc = args[8]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity(PHYLOSEQ, distance, metadata, criteria, nmds, pcoa, method_hc, plot_hc)
}

if (!interactive()) {
        main_wunifrac()
}
