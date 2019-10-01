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
## Copyright (c) Cyril NOEL, aout-2019                                     ####
## Email: cyril.noel@ifremer.fr                                            ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs the beta diversity based on the   ##
##        rarefied abundance table                                           ##   
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

betadiversity_rarefied <- function (PHYLOSEQ, final_rarefied_ASV_table_with_taxonomy, ASV_ordination_plot_rarefied, distance, replicats, ASV_ordination_plot_wrapped_rarefied, samples_ordination_plot_rarefied, split_graph_ordination_plot_rarefied, metadata) {

    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script: analysis of the beta diversity **           ####
    
    #~~~~~~~~~~~~~~~#
    # Rarefied data #
    #~~~~~~~~~~~~~~~#
    PHYLOSEQ_rarefied = rarefy_even_depth(PHYLOSEQ, sample.size=min(sample_sums(PHYLOSEQ)),rngseed=1000,replace=TRUE,trimOTUs=TRUE)
    rarefied_table = cbind(as.data.frame(otu_table(PHYLOSEQ_rarefied)),as.data.frame(tax_table(PHYLOSEQ_rarefied)))
    write.table(rarefied_table,final_rarefied_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=".",quote=F) 
    
    ## /1\ Ordination process ####
    ord_rarefied = ordinate(PHYLOSEQ_rarefied, "NMDS", distance, trymax = 1000)
    
    ## /2\ ASV analysis ####
    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_ord_class = color_vector[1:length(unique(PHYLOSEQ_rarefied@tax_table@.Data[,3]))]
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    
    ## ______ NMDS ####
    plot_ordination(PHYLOSEQ_rarefied,ord_rarefied,type="taxa",color="Class") + 
      theme_classic() +
      geom_point(size=4) +
      scale_color_manual(values=color_ord_class) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      theme(axis.text=element_text(size=12,color="black")) +
      annotate(geom="text",x=min(ord_rarefied$points[,1]),y=max(ord_rarefied$points[,1]),label=paste("Stress:",round(ord_rarefied$stress,4),sep=" "))
    ggsave(filename=ASV_ordination_plot_rarefied,width=13,height=9)
    
    ## ______ wrapped NMDS ####
    plot_ordination(PHYLOSEQ_rarefied,ord_rarefied,type="taxa",color="Class") + 
      theme_classic() +
      geom_point(size=3) +
      theme(legend.position="none") +
      theme(strip.text.x=element_text(size=13,face="bold",color="blue")) +
      scale_color_manual(values=color_ord_class) +
      facet_wrap(~Class,4) +
      theme(axis.text=element_text(size=12,color="black")) 
    ggsave(filename=ASV_ordination_plot_wrapped_rarefied,width=13,height=10)
    
    ## /3\ Sample analysis ####
    color_samples = color_vector[1:length(levels(metadata[,replicats]))]
    group_rarefied = get_variable(PHYLOSEQ_rarefied,replicats)
    anosim_result_rarefied = anosim(distance(PHYLOSEQ_rarefied,"bray"),group_rarefied, permutations = 999)
    plot_ordination(PHYLOSEQ_rarefied,ord_rarefied,type="samples",color=replicats) +
      theme_classic() +
      geom_point(size=3) +
      geom_text(aes(label=rownames(sample_data(PHYLOSEQ_rarefied))),col="black",size=2.5,vjust=2,hjust=1) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      labs(color="Species") +
      theme(axis.text=element_text(size=12,color="black")) +
      scale_fill_manual(values=alpha(color_samples,0.4)) +
      scale_color_manual(values=color_samples) +
      annotate(geom="text",x=min(ord_rarefied$points[,1]),y=max(ord_rarefied$points[,1]),label=paste("Stress:",round(ord_rarefied$stress,4),sep=" ")) +
      stat_ellipse(geom="polygon",alpha=0.1,type="t",aes(fill=Species)) +
      annotate(geom="text",x=min(ord_rarefied$points[,1]),y=max(ord_rarefied$points[,1])-0.3,label=paste("Anosim (based on species) : p-value",anosim_result_rarefied$signif,sep=" "))
    ggsave(filename=samples_ordination_plot_rarefied,width=12,height=10)

   print("plot ordination")
    
    ## /4\ Split graphic ####
    
    plot_ordination(PHYLOSEQ_rarefied,ord_rarefied,type="split",color="Class") +
      theme_classic() +
      geom_text(size=8,aes(label=Labels_NMDS),vjust=-0.8) +
      theme(legend.text=element_text(size=20)) +
      theme(legend.title=element_text(size=22)) +
      theme(strip.text.x=element_text(size=22,face="bold",color="blue")) +
      theme(axis.text=element_text(size=18,color="black")) +
      theme(axis.title=element_text(size=19)) +
      geom_point(size=5) +
      scale_color_manual(values=c("black",color_ord_class))
    ggsave(filename=split_graph_ordination_plot_rarefied,width=20,height=13)
    
    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** End of the script **                                                 ####
}

main <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1]) 
    final_rarefied_ASV_table_with_taxonomy = args[2]
    ASV_ordination_plot_rarefied = args[3]
    distance = args[4]
    replicats = args[5]
    ASV_ordination_plot_wrapped_rarefied = args[6]
    samples_ordination_plot_rarefied = args[7]
    split_graph_ordination_plot_rarefied = args[8]
    metadata = args[9]
    betadiversity_rarefied(PHYLOSEQ, final_rarefied_ASV_table_with_taxonomy, ASV_ordination_plot_rarefied, distance, replicats, ASV_ordination_plot_wrapped_rarefied, samples_ordination_plot_rarefied, split_graph_ordination_plot_rarefied, metadata)
}

if (!interactive()) {
        main()
}
