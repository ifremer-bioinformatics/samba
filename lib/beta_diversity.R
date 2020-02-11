#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Script name: R script for the SAMBA-nextflow workflow                   ####
##                                                                           ##
## Purpose of script: Automated Statistical Analyses of Metabarcoding Data   ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2020-01-31                                                 ####
##                                                                           ##
## Emails: samba-sebimer@ifremer.fr                                        ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## This program is free software: you can redistribute it and/or modify it   ##
## under the terms of the GNU Affero General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or         ##
## (at your option) any later version.                                       ## 
##									                                         ##
## License at https://www.gnu.org/licenses/agpl-3.0.txt 		             ##
## 									                                         ##
## This program is distributed in the hope that it will be useful, but       ##
## WITHOUT ANY WARRANTY; without even the implied warranty of                ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                      ##
## See the GNU Affero General Public License for more details.               ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs the beta diversity (NMDS, PCoa &  ##
##        Hierachical Clustering) on the non-normalized ASV table based on   ##
##        four distance matrices (Jaccard, Bray-Curtis, UniFrac & Weighted   ##
##        UniFrac)                    					                     ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Load up the needed packages ####
requiredPackages = c("dplyr","stringr","ggplot2","RColorBrewer","svglite","tidyr","gridExtra","egg","vegan","dendextend","BiocManager", "phyloseq")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#										                                        #
# Function to standardize ordination analysis and plot for each distance matrix #
#										                                        #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

 
betadiversity <- function(PHYLOSEQ, distance, metadata, variance_significance_tests, criteria, nmds, pcoa, method_hc, plot_hc, plot_pie){
    
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

    ## Statistic about variance explained ####
    ### using user-specified variable #### 
    group = get_variable(PHYLOSEQ, criteria)
    adonis_result = adonis(distance(PHYLOSEQ,distance) ~ group, permutations = 9999)

    ### using all available variables ####
    all_var = c()
    for (var in sample_variables(PHYLOSEQ) ) {
      l = length(levels(as.factor(get_variable(PHYLOSEQ,var))))
      if(l > 1 && l < nsamples(PHYLOSEQ)){
        all_var <- cbind(all_var,var)
      }
    }

    variables = paste(collapse =" + ", all_var )
    
    sink(file = paste(variance_significance_tests,distance,".txt",sep="") , type = "output")
      f  = paste("distance(PHYLOSEQ,distance)"," ~ ", variables)
      cat(sep = "", "###############################################################\n",
                "#Perform Adonis test on multiple variables: ",variables," using the ",distance," distance matrix")
      adonis_all=adonis(as.formula(f), data=metadata, perm = 9999)
      print(adonis_all)
      cat("\n\n")
    sink()
    
    ### Explained variance graphs ####
    ExpVar_perc = adonis_all$aov.tab$R2[-length(adonis_all$aov.tab$R2)]*100
    ExpVar_name = rownames(adonis_all$aov.tab)[-length(rownames(adonis_all$aov.tab))]
    ExpVar_piedata = data.frame(ExpVar_name,ExpVar_perc)
    ExpVar_piedata = ExpVar_piedata[order(ExpVar_piedata$ExpVar_perc),]
    ExpVar_pielabels = sprintf("%s = %3.1f%s", ExpVar_piedata$ExpVar_name,ExpVar_piedata$ExpVar_perc, "%")

    plot.pie(ExpVar_piedata$ExpVar_perc, ExpVar_pielabels, distance, plot_pie, 12, 10)

    ## Ordination plots ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, adonis result, ordination plot name, distance, width of graph, heigth of graph
    plot.nmds(PHYLOSEQ, ord_nmds, criteria, color_samples, adonis_result, nmds, distance, 12, 10, paste("NMDS on non-normalized data","based on",distance,"distance",sep=" "))
    plot.pcoa(PHYLOSEQ, ord_pcoa, criteria, color_samples, adonis_result, pcoa, distance, 12, 10, paste("MDS-PCoA on non-normalized data","based on",distance,"distance",sep=" "))
    
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
    variance_significance_tests=args[9]
    plot_pie = args[10]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity(PHYLOSEQ, distance, metadata, variance_significance_tests, criteria, nmds, pcoa, method_hc, plot_hc, plot_pie)
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
    variance_significance_tests=args[9]
    plot_pie = args[10]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity(PHYLOSEQ, distance, metadata, variance_significance_tests, criteria, nmds, pcoa, method_hc, plot_hc, plot_pie)
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
    variance_significance_tests=args[9]
    plot_pie = args[10]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity(PHYLOSEQ, distance, metadata, variance_significance_tests, criteria, nmds, pcoa, method_hc, plot_hc, plot_pie)
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
    variance_significance_tests=args[9]
    plot_pie = args[10]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity(PHYLOSEQ, distance, metadata, variance_significance_tests, criteria, nmds, pcoa, method_hc, plot_hc, plot_pie)
}

if (!interactive()) {
        main_wunifrac()
}
