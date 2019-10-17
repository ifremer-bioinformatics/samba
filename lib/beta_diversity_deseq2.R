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
##        normalized abundance table using DESeq2                            ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

#### _____________________________________________________________________ ####

## Load up the packages needed ####
library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("svglite")
library("tidyr")
library("gridExtra")
library("egg")
library("DESeq2")
library("vegan")

betadiversity_deseq2 <- function (PHYLOSEQ, final_deseq2_ASV_table_with_taxonomy, distance, criteria, samples_ordination_plot_deseq2, metadata) { 

    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script: analysis of the beta diversity **           ####
    
    #~~~~~~~~~~~~~~~~~~~~~~#
    # DESeq2 normalization #
    #~~~~~~~~~~~~~~~~~~~~~~#
    deseq2 = phyloseq_to_deseq2(PHYLOSEQ , ~ 1)
    gm_mean = function(x, na.rm=TRUE){
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
    
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    ord_deseq2 = ordinate(PHYLOSEQ_deseq2, "NMDS", distance,trymax = 1000)

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]
    
    group_deseq2 = get_variable(PHYLOSEQ_deseq2, criteria)
    PHYLOSEQ_deseq2_dist = phyloseq::distance(PHYLOSEQ_deseq2, "bray")
    anosim_result_deseq2 = anosim(PHYLOSEQ_deseq2_dist,group_deseq2, permutations = 999)
    
    ### /3\ Sample analysis ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, sample ordination plot name, width of graph, heigth of graph, graph title
    sample_nmds(PHYLOSEQ_deseq2, ord_deseq2, criteria, color_samples, anosim_result_deseq2, samples_ordination_plot_deseq2, 12, 10, "NMDS on deseq2 normalized data")
}

main <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    final_deseq2_ASV_table_with_taxonomy = args[2]
    distance = args[3]
    criteria = args[4] 
    samples_ordination_plot_deseq2 = args[5]
    metadata = args[6]
    workflow_dir = args[7]
    if (!exists("sample_nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    betadiversity_deseq2(PHYLOSEQ, final_deseq2_ASV_table_with_taxonomy, distance, criteria, samples_ordination_plot_deseq2, metadata)
}

if (!interactive()) {
        main()
}
