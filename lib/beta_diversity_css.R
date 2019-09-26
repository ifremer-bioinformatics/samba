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
##        normalized abundance table using the CSS algorythm performed in    ##
##        the metagenomeSeq R package                                        ##   
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
library("metagenomeSeq")
library("vegan")

betadiversity_css <- function (PHYLOSEQ, final_css_ASV_table_with_taxonomy, ASV_ordination_plot_css, ASV_ordination_plot_wrapped_css, samples_ordination_plot_css, split_graph_ordination_plot_css, distance, replicats, metadata) {

    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script: analysis of the beta diversity **           ####
    
    #~~~~~~~~~~~~~~~~~~~#
    # CSS normalization #
    #~~~~~~~~~~~~~~~~~~~#
    css = phyloseq_to_metagenomeSeq(PHYLOSEQ)
    p = cumNormStatFast(css)
    css = cumNorm(css, p=p)
    css_norm_factor = normFactors(css)          
    
    CSS_TABLE = MRcounts(css, norm = T)
    
    PHYLOSEQ_css = PHYLOSEQ
    otu_table(PHYLOSEQ_css) = otu_table(CSS_TABLE,taxa_are_rows=TRUE)
    CSS_normalized_table = cbind(as.data.frame(otu_table(PHYLOSEQ_css)),as.data.frame(tax_table(PHYLOSEQ_css)))
    
    write.table(CSS_normalized_table,final_css_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=",") 
    
    ### /1\ Ordination process ####
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    ord_css = ordinate(PHYLOSEQ_css,"NMDS", distance, trymax = 1000)
    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = sample(color_vector,length(levels(metadata[,replicats])))
    color_ord_class = sample(color_vector,length(unique(PHYLOSEQ_css@tax_table@.Data[,3])))
    
    ### /2\ ASV analysis ####
    
    ## ______ NMDS ####
    plot_ordination(PHYLOSEQ_css,ord_css,type="taxa",color="Class") +
      theme_classic() +
      geom_point(size=3) +
      scale_color_manual(values=color_ord_class) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      theme(axis.text=element_text(size=12,color="black")) +
      annotate(geom="text",x=min(ord_css$points[,1]),y=max(ord_css$points[,1]),label=paste("Stress:",round(ord_css$stress,4),sep=" "))
    ggsave(filename=ASV_ordination_plot_css,width=13,height=9)
    
    ## ______ wrapped NMDS ####
    plot_ordination(PHYLOSEQ_css,ord_css,type="taxa",color="Class") +
      theme_classic() +
      geom_point(size=3) +
      theme(legend.position="none") +
      theme(strip.text.x=element_text(size=13,face="bold",color="blue")) +
      scale_color_manual(values=color_ord_class) +
      facet_wrap(~Class,4) +
      theme(axis.text=element_text(size=12,color="black"))
    ggsave(filename=ASV_ordination_plot_wrapped_css,width=13,height=10)
    
    ### /3\ Sample analysis ####
    group_css = get_variable(PHYLOSEQ_css, replicats)
    anosim_result_css = anosim(distance(PHYLOSEQ_css,"bray"),group_css, permutations = 999)
    
    plot_ordination(PHYLOSEQ_css,ord_css,type="samples",color=replicats) +
      theme_classic() +
      geom_point(size=4) +
      geom_text(aes(label=rownames(sample_data(PHYLOSEQ_css))),col="black",size=2.5,vjust=2,hjust=1) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      labs(color="Species") +
      theme(axis.text=element_text(size=12,color="black")) +
      scale_fill_manual(values=alpha(color_samples,0.4)) +
      scale_color_manual(values=color_samples) +
      stat_ellipse(geom="polygon",alpha=0.1,type="t",aes(fill=Species)) +
      annotate(geom="text",x=min(ord_css$points[,1])+0.01,y=max(ord_css$points[,1]),label=paste("Stress:",round(ord_css$stress,4),
                                                                                             "\nANOSIM statistic R:",round(anosim_result_css$statistic,4),
                                                                                             "\nAnosim (based on Sample_Location) : p-value",anosim_result_css$signif,sep=" "))

    ggsave(filename=samples_ordination_plot_css,width=12,height=10)
    
    ### /4\ Split graphic ####
    
    plot_ordination(PHYLOSEQ_css,ord_css,type="split",color="Class") +
      theme_classic() +
      geom_text(size=8,aes(label=Labels_NMDS),vjust=-0.8) +
      theme(legend.text=element_text(size=20)) +
      theme(legend.title=element_text(size=22)) +
      theme(strip.text.x=element_text(size=22,face="bold",color="blue")) +
      theme(axis.text=element_text(size=18,color="black")) +
      theme(axis.title=element_text(size=19)) +
      geom_point(size=5) +
      scale_color_manual(values=c("black",color_ord_class))
    ggsave(filename=split_graph_ordination_plot_css,width=20,height=13)
    
    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** End of the script **                                                 ####
}

main <- function(){
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    final_css_ASV_table_with_taxonomy = args[2]
    ASV_ordination_plot_css = args[3]
    ASV_ordination_plot_wrapped_css = args[4]
    samples_ordination_plot_css = args[5]
    split_graph_ordination_plot_css = args[6]
    distance = args[7]
    replicats = args[8]
    metadata = args[9]
    betadiversity_css(PHYLOSEQ, final_css_ASV_table_with_taxonomy, ASV_ordination_plot_css, ASV_ordination_plot_wrapped_css, samples_ordination_plot_css, split_graph_ordination_plot_css, distance, replicats, metadata)
}

if (!interactive()) {
        main()
}

