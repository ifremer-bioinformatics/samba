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

betadiversity <- function(PHYLOSEQ, ASV_ordination_plot, distance, replicats, ASV_ordination_plot_wrapped, samples_ordination_plot, split_graph_ordination_plot, metadata) {

    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script: analysis of the beta diversity **           ####
    
    #~~~~~~~~~~~~~~~~~~~~~#
    # Non-normalized data #
    #~~~~~~~~~~~~~~~~~~~~~#
    
    ## /1\ Ordination process ####
    ord = ordinate(PHYLOSEQ,"NMDS",distance,trymax = 1000)
    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)

    ## /2\ ASV analysis ####
    
    ## ______ NMDS ####
#    plot_ordination(PHYLOSEQ,ord,type="taxa",color="Class") + 
#      theme_classic() +
#      geom_point(size=4) +
#      scale_color_manual(values=color_ord_class) +
#      theme(legend.text=element_text(size=13)) +
#      theme(legend.title=element_text(size=14,face="bold")) +
#      theme(axis.text=element_text(size=12,color="black")) +
#      annotate(geom="text",x=min(ord$points[,1]),y=max(ord$points[,1]),label=paste("Stress:",round(ord$stress,4),sep=" "))
#    ggsave(filename=ASV_ordination_plot,width=20,height=9)
#    
#    ## ______ wrapped NMDS ####
#    plot_ordination(PHYLOSEQ,ord,type="taxa",color="Class") + 
#      theme_classic() +
#      geom_point(size=3) +
#      theme(legend.position="none") +
#      theme(axis.text=element_text(size=12,color="black")) +
#      theme(strip.text.x=element_text(size=13,face="bold",color="blue")) +
#      scale_color_manual(values=color_ord_class) +
#      facet_wrap(~Class,4)
#    ggsave(filename=ASV_ordination_plot_wrapped,width=13,height=10)
#    
    ## /3\ Sample analysis ####
    criteria = "ExtractMethod"    
    color_ord_class = color_vector[1:length(unique(PHYLOSEQ@tax_table@.Data[,4]))]
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]
    group = get_variable(PHYLOSEQ, criteria)
    anosim_result = anosim(distance(PHYLOSEQ,"bray"),group, permutations = 999)
    
    #plot_ordination(PHYLOSEQ,ord,type="samples",color=replicats) +
    #plot_ordination(PHYLOSEQ,ord,type="samples",color=criteria, shape="Group") +
    plot_ordination(PHYLOSEQ,ord,type="samples",color=criteria) +
      theme_classic() +
      geom_point(size=3) +
      geom_text(aes(label=rownames(sample_data(PHYLOSEQ))),col="black",size=2.5,vjust=2,hjust=1) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      #labs(color="Polymer", shape="Group") +
      labs(color=criteria) +
      theme(axis.text=element_text(size=12,color="black")) +
      scale_fill_manual(values=alpha(color_samples,0.4)) +
      scale_color_manual(values=color_samples) +
      stat_ellipse(geom="polygon",alpha=0.1,type="t",aes(fill=ExtractMethod)) +
      labs(caption = paste("Stress:",round(ord$stress,4),
          "\nANOSIM statistic R:",round(anosim_result$statistic,4),
          paste("\nAnosim based on ", criteria,": p-value"),anosim_result$signif,sep=" ")) 
    ggsave(filename=samples_ordination_plot,width=12,height=10)
    
    ## /4\ Split graphic ####
#    plot_ordination(PHYLOSEQ,ord,type="split",color="Class") +
#      theme_classic() +
#      geom_text(size=8,aes(label=Group),vjust=-0.8) +
#      theme(legend.text=element_text(size=20)) +
#      theme(legend.title=element_text(size=22)) +
#      theme(strip.text.x=element_text(size=22,face="bold",color="blue")) +
#      theme(axis.text=element_text(size=18,color="black")) +
#      theme(axis.title=element_text(size=19)) + 
#      geom_point(size=5) +
#      scale_color_manual(values=c("black",color_ord_class))
#    ggsave(filename=split_graph_ordination_plot,width=20,height=13)
}

main  <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    ASV_ordination_plot = args[2]
    distance = args[3]
    replicats = args[4]
    ASV_ordination_plot_wrapped = args[5]
    samples_ordination_plot = args[6]
    split_graph_ordination_plot = args[7]
    metadata = args[8]
    betadiversity(PHYLOSEQ, ASV_ordination_plot, distance, replicats, ASV_ordination_plot_wrapped, samples_ordination_plot, split_graph_ordination_plot, metadata)
}

if (!interactive()) {
        main()
}

