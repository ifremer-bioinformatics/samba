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
## Modified on: 2019-10-23                                                 ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## Emails: cyril.noel@ifremer.fr & laure.quintric@ifremer.fr               ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs only the alpha diversity and      ##
##        taxonomic diversity based on your data                             ##   
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

alphadiversity <- function(PHYLOSEQ, alpha_div_plots, barplot_relabund_phylum, barplot_relabund_family, barplot_relabund_genus, heatmap_class, heatmap_family, heatmap_genus, threshold, distance, group){
    
    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    alpha_rich = estimate_richness(PHYLOSEQ,measures=c("Observed","Chao1","Shannon","InvSimpson"))
    library("microbiome")
    evenness = evenness(PHYLOSEQ,"pielou")
    detach("package:microbiome", unload=TRUE)
    alpha_rich$Pielou = evenness$pielou
    df = data.frame(alpha_rich,sample_data(PHYLOSEQ))
    df2 = gather(df,key="Measure",value="Value",Observed,Chao1,Shannon,InvSimpson,Pielou)
    df2$Measure = factor(df2$Measure,levels=c("Observed","Chao1","InvSimpson","Shannon","Pielou"))
    
    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script : Analysis of the alpha diversity **         ####
    
    #### /1\ Alpha diversity ####
    
    ## ___ Process of the analysis ####
    plot_alpha_global=ggplot(df2, aes_string(x=group,y="Value")) +
      facet_wrap(~Measure, scale="free") +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text=element_text(size=12,color="black")) +
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)) +
      theme(strip.text.x=element_text(size = 14,face="bold",color="blue")) +
      theme(strip.background = element_rect(fill="grey")) +
      theme(axis.title.x=element_text(size=14,face="bold",vjust=-1)) +
      xlab("Samples") +
      theme(axis.title.y=element_text(size=14,face="bold"))
    plot1_alpha = plot_alpha_global %+% subset(df2 , Measure == "Observed" | Measure == "Chao1" | Measure == "InvSimpson") + labs(x=NULL)
    plot2_alpha = plot_alpha_global %+% subset(df2 , Measure == "Shannon" | Measure == "Pielou")
    
    final_alpha_plot = arrangeGrob(grobs=lapply(list(plot1_alpha,plot2_alpha),set_panel_size,width=unit(10,"cm"),height=unit(10,"cm")))
    ggsave(filename=alpha_div_plots,final_alpha_plot, width=14, height=14)

    ## ___ Statistical significance of the index ####
    anova_data = cbind(sample_data(PHYLOSEQ), alpha_rich) 

    #Anova on Observed richness
    anova.Observed=aov(Observed~group,anova_data)
    anova.Observed.res=summary(anova.Observed)

    #Anova on Chao1 index
    anova.Chao1=aov(Chao1~group,anova_data)
    anova.Chao1.res=summary(anova.Chao1)

    #Anova on Shannon index
    anova.Shannon=aov(Shannon~group,anova_data)
    anova.Shannon.res=summary(anova.Shannon)

    #Anova on Inv. Simpson index
    anova.InvSimpson=aov(InvSimpson~group,anova_data)
    anova.InvSimpson.res=summary(anova.InvSimpson)

    #Anova on Pielou index
    anova.Pielou=aov(Pielou~group,anova_data)
    anova.Pielou.res=summary(anova.Pielou)

    #Output of significance test
    index.list=list(Anova.Observed=anova.Observed.res,Anova.Chao1=anova.Chao1.res,Anova.Shannon=anova.Shannon.res,Anova.InvSimpson=anova.InvSimpson.res,Anova.Pielou=anova.Pielou.res)
    capture.output(print(index.list),file=paste(FIGURES.alpha,"index_significance_tests.txt",sep=""))

    #### /2\ Taxonomic diversity ####
    
    ## ___ Barplot representation ####
    ## ______ at the phylum level ####
    Relabund_phylum = PHYLOSEQ %>%
      tax_glom(taxrank="Phylum") %>%
      transform_sample_counts(function(x){x/sum(x)*100}) %>%
      psmelt() %>%
      arrange(Phylum)
    
    color_phylum = color_vector[1:length(levels(Relabund_phylum$Phylum))]
    
    ggplot(Relabund_phylum,aes(x=Sample,y=Abundance,fill=Phylum)) +
      geom_bar(stat = "identity",position="fill") +
      facet_wrap(group, nrow=1, scale="free") +
      theme_classic() +
      scale_y_continuous(expand=c(0,0),labels=c("0","25","50","75","100")) +
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black",size=10)) +
      theme(axis.title.x=element_text(vjust=-1,color="black",size=11)) +
      theme(axis.text.y=element_text(hjust=0.8,color="black",size=10)) +
      theme(axis.title.y=element_text(size=11)) +
      theme(legend.text=element_text(size=12)) +
      theme(legend.title=element_text(size=13,face="bold")) +
      xlab("Samples") +
      scale_fill_manual(values=color_phylum)
    ggsave(filename=barplot_relabund_phylum,width=20,height=12)


    ## ______ at the family level ####
    Relabund_family = PHYLOSEQ %>%
      tax_glom(taxrank="Family") %>%
      transform_sample_counts(function(x){x/sum(x)*100}) %>%
      psmelt() %>%
      filter(Abundance > threshold) %>%
      arrange(Family)
    
    color_family = color_vector[1:length(levels(Relabund_family$Family))]
    
    ggplot(Relabund_family,aes(x=Sample,y=Abundance,fill=Family)) +
      geom_bar(stat = "identity",position="fill") +
      facet_wrap(group, nrow=1, scale="free") +
      theme_classic() +
      scale_y_continuous(expand=c(0,0),labels=c("0","25","50","75","100")) +
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black",size=10)) +
      theme(axis.title.x=element_text(vjust=-1,color="black")) +
      theme(axis.text.y=element_text(hjust=0.8,color="black",size=10)) +
      theme(axis.title.y=element_text(size=11)) +
      theme(legend.text=element_text(size=12)) +
      theme(legend.title=element_text(size=13,face="bold")) +
      xlab("Samples") +
      scale_fill_manual(values=color_family)
    ggsave(filename=barplot_relabund_family,width=20,height=12)

    ## ______ at the genus level ####
    Relabund_genus = PHYLOSEQ %>%
      tax_glom(taxrank="Genus") %>%
      transform_sample_counts(function(x){x/sum(x)*100}) %>%
      psmelt() %>%
      filter(Abundance > threshold) %>%
      arrange(Genus)
    
    color_genus = color_vector[1:length(levels(Relabund_genus$Genus))]
    
    ggplot(Relabund_genus, aes(x = Sample, y = Abundance, fill = Genus)) +
      geom_bar(stat = "identity",position="fill") +
      facet_wrap(group, nrow=1, scale="free") +
      theme_classic() +
      scale_y_continuous(labels=c("0","25","50","75","100"),expand=c(0,0)) +
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black",size=10)) +
      theme(axis.title.x=element_text(vjust=-1,color="black")) +
      theme(axis.text.y=element_text(hjust=0.8,color="black",size=10)) +
      theme(axis.title.y=element_text(size=11)) +
      theme(legend.text=element_text(size=12)) +
      theme(legend.title=element_text(size=13,face="bold")) +
      xlab("Samples") +
      scale_fill_manual(values=color_genus)
    ggsave(filename=barplot_relabund_genus,width=20,height=12)
}

main <- function() {
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    threshold = args[2]
    distance = args[3]
    alpha_div_plots = args[4]
    barplot_relabund_phylum = args[5]
    barplot_relabund_family = args[6]
    barplot_relabund_genus = args[7]
    heatmap_class = args[8]
    heatmap_family = args[9]
    heatmap_genus = args[10]
    #get group variable an replace "-" by "_"
    group = str_replace(args[11], "-", "_")
    #Run alpha diversity calculations
    alphadiversity(PHYLOSEQ, alpha_div_plots, barplot_relabund_phylum, barplot_relabund_family, barplot_relabund_genus, heatmap_class, heatmap_family, heatmap_genus, threshold, distance, group)
}
if (!interactive()) {
        main()
}
