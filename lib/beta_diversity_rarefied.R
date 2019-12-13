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
##        Hierachical Clustering) on the rarefied ASV table based on four    ##
##        distance matrices (Jaccard, Bray-Curtis, UniFrac & Weighted        ##
##        UniFrac)                                                           ##   
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
library("dendextend")

# @@@@@@@@@@@@@@@@@@@ #
#                     #
# Rarefaction process #
#                     #
# @@@@@@@@@@@@@@@@@@@ #

args = commandArgs(trailingOnly=TRUE)
PHYLOSEQ = readRDS(args[1])
final_rarefied_ASV_table_with_taxonomy = args[2]

PHYLOSEQ_rarefied = rarefy_even_depth(PHYLOSEQ, sample.size=min(sample_sums(PHYLOSEQ)),rngseed=1000,replace=TRUE,trimOTUs=TRUE)
rarefied_table = cbind(as.data.frame(otu_table(PHYLOSEQ_rarefied)),as.data.frame(tax_table(PHYLOSEQ_rarefied)))
write.table(rarefied_table,final_rarefied_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=".",quote=F)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#										#
# Function to standardize ordination analysis and plot for each distance matrix #
#										#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

betadiversity_rarefied <- function (PHYLOSEQ_rarefied, distance, metadata, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc) {
   
    #~~~~~~~~~~~~~~~#
    # Rarefied data #
    #~~~~~~~~~~~~~~~#

    ## Ordination process ####
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    
    if (distance == "jaccard" | distance == "bray") {
      ord_rarefied_nmds = ordinate(PHYLOSEQ_rarefied, "NMDS", distance, trymax = 1000)
      ord_rarefied_pcoa = ordinate(PHYLOSEQ_rarefied, "PCoA", distance, trymax = 1000)
    }
    else {
      ord_rarefied_nmds = ordinate(PHYLOSEQ_rarefied, "NMDS", distance)
      ord_rarefied_pcoa = ordinate(PHYLOSEQ_rarefied, "PCoA", distance)
    }

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]
    
    group_rarefied = get_variable(PHYLOSEQ_rarefied, criteria)
    anosim_result_rarefied = anosim(distance(PHYLOSEQ_rarefied,distance),group_rarefied, permutations = 999)

    ## Sample analysis ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, ordination plot name, width of graph, heigth of graph, graph title
    plot.nmds(PHYLOSEQ_rarefied, ord_rarefied_nmds, criteria, color_samples, anosim_result_rarefied, nmds_rarefied, distance, 12, 10, paste("NMDS on rarefied data","based on",distance,"distance",sep=" "))
    plot.pcoa(PHYLOSEQ_rarefied, ord_rarefied_pcoa, criteria, color_samples, anosim_result_rarefied, pcoa_rarefied, distance, 12, 10, paste("MDS-PCoA on rarefied data","based on",distance,"distance",sep=" "))

    ## Hierarchical clsutering ####
    dist = distance(PHYLOSEQ_rarefied, distance, type="samples")
    hc = hclust(dist, method = method_hc)
    dendro = as.dendrogram(hc)
    group = data.frame(PHYLOSEQ_rarefied@sam_data[,criteria])[,1]
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

main_jaccard <- function(){
    PHYLOSEQ_rarefied = PHYLOSEQ_rarefied 
    distance = "jaccard"
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_rarefied = args[6]
    pcoa_rarefied = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity_rarefied(PHYLOSEQ_rarefied, distance, metadata, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc)
}

if (!interactive()) {
        main_jaccard()
}

main_bray <- function(){
    PHYLOSEQ_rarefied = PHYLOSEQ_rarefied
    distance = "bray"
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_rarefied = args[6]
    pcoa_rarefied = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity_rarefied(PHYLOSEQ_rarefied, distance, metadata, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc)
}

if (!interactive()) {
        main_bray()
}

main_unifrac <- function(){
    PHYLOSEQ_rarefied = PHYLOSEQ_rarefied
    distance = "unifrac"
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_rarefied = args[6]
    pcoa_rarefied = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity_rarefied(PHYLOSEQ_rarefied, distance, metadata, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc)
}

if (!interactive()) {
        main_unifrac()
}

main_wunifrac <- function(){
    PHYLOSEQ_rarefied = PHYLOSEQ_rarefied
    distance = "wunifrac"
    # get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_rarefied = args[6]
    pcoa_rarefied = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity_rarefied(PHYLOSEQ_rarefied, distance, metadata, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc)
}

if (!interactive()) {
        main_wunifrac()
}
