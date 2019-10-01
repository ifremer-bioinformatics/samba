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

betadiversity_deseq2 <- function (PHYLOSEQ, final_deseq2_ASV_table_with_taxonomy, ASV_ordination_plot_deseq2, ASV_ordination_plot_wrapped_deseq2, samples_ordination_plot_deseq2, split_graph_ordination_plot_deseq2, distance, replicats, metadata) {

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
    
    ### /1\ Ordination process ####
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    ord_deseq2 = ordinate(PHYLOSEQ_deseq2, "NMDS", distance,trymax = 1000)
    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,replicats]))]
    color_ord_class = color_vector[1:length(unique(PHYLOSEQ_deseq2@tax_table@.Data[,3]))]
    
    ### /2\ ASV analysis ####
    
    ## ______ NMDS ####
    plot_ordination(PHYLOSEQ_deseq2,ord_deseq2,type="taxa",color="Class") +
      theme_classic() +
      geom_point(size=4) +
      scale_color_manual(values=color_ord_class) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      theme(axis.text=element_text(size=12,color="black")) +
      annotate(geom="text",x=min(ord_deseq2$points[,1]),y=max(ord_deseq2$points[,1]),label=paste("Stress:",round(ord_deseq2$stress,4),sep=" "))
    ggsave(filename=ASV_ordination_plot_deseq2,width=13,height=9)
    
    ### ______ wrapped NMDS ####
    plot_ordination(PHYLOSEQ_deseq2,ord_deseq2,type="taxa",color="Class") +
      theme_classic() +
      geom_point(size=3) +
      theme(legend.position="none") +
      theme(strip.text.x=element_text(size=13,face="bold",color="blue")) +
      scale_color_manual(values=color_ord_class) +
      facet_wrap(~Class,4) +
      theme(axis.text=element_text(size=12,color="black"))
    ggsave(filename=ASV_ordination_plot_wrapped_deseq2,width=13,height=10)
    
    ### /3\ Sample analysis ####
    group_deseq2 = get_variable(PHYLOSEQ_deseq2, replicats)
    PHYLOSEQ_deseq2_dist = phyloseq::distance(PHYLOSEQ_deseq2, "bray")
    anosim_result_deseq2 = anosim(PHYLOSEQ_deseq2_dist,group_deseq2, permutations = 999)
    
    plot_ordination(PHYLOSEQ_deseq2,ord_deseq2,type="samples",color=replicats) +
      theme_classic() +
      geom_point(size=3) +
      geom_text(aes(label=rownames(sample_data(PHYLOSEQ_deseq2))),col="black",size=2.5,vjust=2,hjust=1) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      labs(color="Species") +
      theme(axis.text=element_text(size=12,color="black")) +
      scale_fill_manual(values=alpha(color_samples,0.4)) +
      scale_color_manual(values=color_samples) +
      stat_ellipse(geom="polygon",alpha=0.1,type="t",aes(fill=Species)) +
      annotate(geom="text",x=min(ord_deseq2$points[,1])+0.01,y=max(ord_deseq2$points[,1]),label=paste("Stress:",round(ord_deseq2$stress,4),
                                                                                             "\nANOSIM statistic R:",round(anosim_result_deseq2$statistic,4),
                                                                                             "\nAnosim (based on Sample_Location) : p-value",anosim_result_deseq2$signif,sep=" "))

    ggsave(filename=samples_ordination_plot_deseq2,width=12,height=10)
    
    ### /4\ Split graphic ####
    
    plot_ordination(PHYLOSEQ_deseq2,ord_deseq2,type="split",color="Class") +
      theme_classic() +
      geom_text(size=8,aes(label=Labels_NMDS),vjust=-0.8) +
      theme(legend.text=element_text(size=20)) +
      theme(legend.title=element_text(size=22)) +
      theme(strip.text.x=element_text(size=22,face="bold",color="blue")) +
      theme(axis.text=element_text(size=18,color="black")) +
      theme(axis.title=element_text(size=19)) +
      geom_point(size=5) +
      scale_color_manual(values=c("black",color_ord_class))
    ggsave(filename=split_graph_ordination_plot_deseq2,width=20,height=13)
    
    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** End of the script **                                                 ####
}

main <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    final_deseq2_ASV_table_with_taxonomy = args[2]
    ASV_ordination_plot_deseq2 = args[3]
    ASV_ordination_plot_wrapped_deseq2 = args[4]
    samples_ordination_plot_deseq2 = args[5]
    split_graph_ordination_plot_deseq2 = args[6]
    distance = args[7]
    replicats = args[8] 
    metadata = args[9]
    betadiversity_deseq2(PHYLOSEQ, final_deseq2_ASV_table_with_taxonomy, ASV_ordination_plot_deseq2, ASV_ordination_plot_wrapped_deseq2, samples_ordination_plot_deseq2, split_graph_ordination_plot_deseq2, distance, replicats, metadata)
}

if (!interactive()) {
        main()
}
