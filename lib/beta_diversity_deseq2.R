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
##        Hierachical Clustering) on DESeq2 normalized ASV table based on    ##
##        four distance matrices (Jaccard, Bray-Curtis, UniFrac & Weighted   ##
##        UniFrac) 	 						     ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Install (if necessary) and load up the needed packages ####
requiredPackages_CRAN = c("dplyr","stringr","ggplot2","RColorBrewer","svglite","tidyr","gridExtra","egg","vegan","dendextend","BiocManager")
for(package in requiredPackages_CRAN){
  if(!require(package,character.only = TRUE)) install.packages(package)
  library(package,character.only = TRUE)
}

requiredPackages_BIOCONDUCTOR = c("phyloseq"; "DESeq2")
for(package in requiredPackages_BIOCONDUCTOR){
  if(!require(package,character.only = TRUE)) BiocManager::install(package)
  library(package,character.only = TRUE)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#                              #
# DESeq2 normalization process #
#                              #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

# Get arguments from RScript command line
args = commandArgs(trailingOnly=TRUE)

PHYLOSEQ = readRDS(args[1])
final_deseq2_ASV_table_with_taxonomy = args[2]

deseq2 = phyloseq_to_deseq2(PHYLOSEQ , ~ 1)
gm_mean = function(x, na.rm=TRUE) {
                   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
                   }
geoMeans = apply(counts(deseq2), 1, gm_mean)
deseq2 = estimateSizeFactors(deseq2,geoMeans=geoMeans)
deseq2 = estimateDispersions(deseq2,fitType="parametric",maxit=1000)
deseq2_vst = getVarianceStabilizedData(deseq2)
PHYLOSEQ_deseq2 = PHYLOSEQ
otu_table(PHYLOSEQ_deseq2) = otu_table(deseq2_vst,taxa_are_rows=TRUE)
otu_table(PHYLOSEQ_deseq2)[otu_table(PHYLOSEQ_deseq2) <0 ] <-0
DESeq2_normalized_table = cbind(as.data.frame(otu_table(PHYLOSEQ_deseq2)),as.data.frame(tax_table(PHYLOSEQ_deseq2)))
write.table(DESeq2_normalized_table,final_deseq2_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=",")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#										#
# Function to standardize ordination analysis and plot for each distance matrix #
#										#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

betadiversity_deseq2 <- function (PHYLOSEQ_deseq2, distance, metadata, criteria, nmds_deseq2, pcoa_deseq2, method_hc, plot_hc) { 
  
    #~~~~~~~~~~~~~~~~~~~~~~#
    # DESeq2 normalization #
    #~~~~~~~~~~~~~~~~~~~~~~#
    
    ## Ordination process ####     
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    
    if (distance == "jaccard" | distance == "bray") {
      ord_deseq2_nmds = ordinate(PHYLOSEQ_deseq2, "NMDS", distance,trymax = 1000)
      ord_deseq2_pcoa = ordinate(PHYLOSEQ_deseq2, "PCoA", distance,trymax = 1000)
    }
    else {
      ord_deseq2_nmds = ordinate(PHYLOSEQ_deseq2, "NMDS", distance)
      ord_deseq2_pcoa = ordinate(PHYLOSEQ_deseq2, "PCoA", distance)
    }

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]
    
    group_deseq2 = get_variable(PHYLOSEQ_deseq2, criteria)
    PHYLOSEQ_deseq2_dist = phyloseq::distance(PHYLOSEQ_deseq2, distance)
    anosim_result_deseq2 = anosim(PHYLOSEQ_deseq2_dist,group_deseq2, permutations = 999)
    
    ## Sample analysis ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, ordination plot name, distance, width of graph, heigth of graph, graph title
    plot.nmds(PHYLOSEQ_deseq2, ord_deseq2_nmds, criteria, color_samples, anosim_result_deseq2, nmds_deseq2, distance, 12, 10, paste("NMDS on deseq2 normalized data","based on",distance,"distance",sep=" "))
    plot.pcoa(PHYLOSEQ_deseq2, ord_deseq2_pcoa, criteria, color_samples, anosim_result_deseq2, pcoa_deseq2, distance, 12, 10, paste("MDS-PCoA on deseq2 normalized data","based on",distance,"distance",sep=" "))

    ## Hierarchical clsutering ####    
    hc = hclust(PHYLOSEQ_deseq2_dist, method = method_hc)
    dendro = as.dendrogram(hc)
    group = data.frame(PHYLOSEQ_deseq2@sam_data[,criteria])[,1]
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
    PHYLOSEQ_deseq2 = PHYLOSEQ_deseq2
    distance = "jaccard"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_deseq2 = args[6]
    pcoa_deseq2 = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_deseq2(PHYLOSEQ_deseq2, distance, metadata, criteria, nmds_deseq2, pcoa_deseq2, method_hc, plot_hc)
}

if (!interactive()) {
        main_jaccard()
}


main_bray <- function(){
    PHYLOSEQ_deseq2 = PHYLOSEQ_deseq2
    distance = "bray"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_deseq2 = args[6]
    pcoa_deseq2 = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_deseq2(PHYLOSEQ_deseq2, distance, metadata, criteria, nmds_deseq2, pcoa_deseq2, method_hc, plot_hc)
}

if (!interactive()) {
        main_bray()
}

main_unifrac <- function(){
    PHYLOSEQ_deseq2 = PHYLOSEQ_deseq2
    distance = "unifrac"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_deseq2 = args[6]
    pcoa_deseq2 = args[7] 
    method_hc = args[8]
    plot_hc = args[9]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_deseq2(PHYLOSEQ_deseq2, distance, metadata, criteria, nmds_deseq2,  pcoa_deseq2, method_hc, plot_hc)
}

if (!interactive()) {
        main_unifrac()
}

main_wunifrac <- function(){
    PHYLOSEQ_deseq2 = PHYLOSEQ_deseq2
    distance = "wunifrac"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_deseq2 = args[6]
    pcoa_deseq2 = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_deseq2(PHYLOSEQ_deseq2, distance, metadata, criteria, nmds_deseq2, pcoa_deseq2, method_hc, plot_hc)
}

if (!interactive()) {
        main_wunifrac()
}
