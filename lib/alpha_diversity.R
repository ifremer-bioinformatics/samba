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
## Notes: This part of the script performs only the alpha diversity and      ##
##        taxonomic diversity based on your data                             ##   
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

alphadiversity <- function(biom_tsv, metadata, alpha_div_plots, barplot_relabund_phylum, barplot_relabund_family, barplot_relabund_genus, heatmap_class, heatmap_family, heatmap_genus, threshold, distance){
    #Input data
    rawASVtable = read.table(biom_tsv, h=T, sep="\t", dec=".", check.names=FALSE)
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    
    #Reformatting of the input data
    abund = rawASVtable[,1:length(rawASVtable)-1]
    row.names(abund) = abund$ASV_ID
    abund = abund %>% select (-ASV_ID)
    tax = rawASVtable[,c(1,length(names(rawASVtable)))]
    tax = data.frame(tax$ASV_ID,do.call(rbind,list(str_split_fixed(tax$taxonomy, ";",7))))
    colnames(tax) = c("ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    row.names(tax) = tax$ASV_ID
    tax = as.matrix(tax %>% select (-ASV_ID))
    
    ## Construction of the phyloseq object ####
    ABUND = otu_table(abund,taxa_are_rows=TRUE)
    TAX = tax_table(tax)
    METADATA = sample_data(metadata)
    PHYLOSEQ = phyloseq(ABUND,TAX,METADATA)
    
    #### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
    ## ** Beginning of the script : Analysis of the alpha diversity **         ####
    
    #### /1\ Alpha diversity ####
    
    ## ___ Process of the analysis ####
    alpha_rich = estimate_richness(PHYLOSEQ,measures=c("Observed","Chao1","Shannon","InvSimpson"))
    library("microbiome")
    evenness = evenness(PHYLOSEQ,"pielou")
    detach("package:microbiome", unload=TRUE)
    alpha_rich$Pielou = evenness$pielou
    df = data.frame(alpha_rich,sample_data(PHYLOSEQ))
    df2 = gather(df,key="Measure",value="Value",Observed,Chao1,Shannon,InvSimpson,Pielou)
    df2$Measure = factor(df2$Measure,levels=c("Observed","Chao1","InvSimpson","Shannon","Pielou"))
    plot_alpha_global=ggplot(df2, aes(x=Replicats,y=Value)) +
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
    ggsave(filename=alpha_div_plots,final_alpha_plot, width=14, height=11)

    #### /2\ Taxonomic diversity ####
    
    ## ___ Barplot representation ####
    
    ## ______ at the phylum level ####
    Relabund_phylum = PHYLOSEQ %>%
      tax_glom(taxrank="Phylum") %>%
      transform_sample_counts(function(x){x/sum(x)*100}) %>%
      psmelt() %>%
      arrange(Phylum)
    
    color_phylum = sample(color_vector,length(levels(Relabund_phylum$Phylum)))
    
    ggplot(Relabund_phylum,aes(x=Sample,y=Abundance,fill=Phylum)) +
      geom_bar(stat = "identity",position="fill") +
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
    ggsave(filename=barplot_relabund_phylum,width=12,height=10)


    ## ______ at the family level ####
    Relabund_family = PHYLOSEQ %>%
      tax_glom(taxrank="Family") %>%
      transform_sample_counts(function(x){x/sum(x)*100}) %>%
      psmelt() %>%
      filter(Abundance > threshold) %>%
      arrange(Family)
    
    color_family = sample(color_vector,length(unique(Relabund_family$Family)))
    
    ggplot(Relabund_family,aes(x=Sample,y=Abundance,fill=Family)) +
      geom_bar(stat = "identity",position="fill") +
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
    ggsave(filename=barplot_relabund_family,width=12,height=10)

    ## ______ at the genus level ####
    Relabund_genus = PHYLOSEQ %>%
      tax_glom(taxrank="Genus") %>%
      transform_sample_counts(function(x){x/sum(x)*100}) %>%
      psmelt() %>%
      filter(Abundance > threshold) %>%
      arrange(Genus)
    
    color_genus = sample(color_vector,length(unique(Relabund_genus$Genus)))
    
    ggplot(Relabund_genus, aes(x = Sample, y = Abundance, fill = Genus)) +
      geom_bar(stat = "identity",position="fill") +
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
    ggsave(filename=barplot_relabund_genus,width=12,height=10)

    ## ___ Heatmap representation ####
    
    ## ______ at the class level ####
    heatmap_data_class = tax_glom(PHYLOSEQ, taxrank="Class")
    plot_heatmap(heatmap_data_class,"NMDS",distance,sample.label="Labels_NMDS",taxa.label="Class",taxa.order="Class",trans=NULL, low="#000033",high="#CCFF66",na.value="#000033") +
     theme_classic() +
     theme(axis.text.x=element_text(angle=90,size=14,vjust=0.3,color="black")) +
     theme(axis.text.y=element_text(size=14,color="black")) +
     theme(axis.title.y=element_text(size=16,face="bold")) +
     theme(legend.title=element_text(size=14)) +
     theme(legend.text=element_text(size=12))
    ggsave(filename=heatmap_class,width=15,height=4)
    
    ## ______ at the family level ####
    heatmap_data_family = tax_glom(PHYLOSEQ, taxrank="Family")
    plot_heatmap(heatmap_data_family, "NMDS",distance,sample.label="Labels_NMDS",taxa.label="Family",taxa.order="Family",trans=NULL,low="#000033", high="#CCFF66",na.value="#000033") +
     theme_classic() +
     theme(axis.text.x=element_text(angle=90,size=14,vjust=0.3,color="black")) +
     theme(axis.text.y=element_text(size=14,color="black")) +
     theme(axis.title.y=element_text(size=16,face="bold")) +
     theme(legend.title=element_text(size=14)) +
     theme(legend.text=element_text(size=12))
    ggsave(filename=heatmap_family,width=15,height=13)
    
    ## ______ at the genus level ####
    heatmap_data_genus = tax_glom(PHYLOSEQ, taxrank="Genus")
    plot_heatmap(heatmap_data_genus,"NMDS",distance,sample.label="Labels_NMDS",taxa.label="Genus",taxa.order="Genus",trans=NULL,low="#000033", high="#CCFF66",na.value="#000033") +
     theme_classic() +
     theme(axis.text.x=element_text(angle=90,size=14,vjust=0.3,color="black")) +
     theme(axis.text.y=element_text(size=14,color="black")) +
     theme(axis.title.y=element_text(size=16,face="bold")) +
     theme(legend.title=element_text(size=14)) +
     theme(legend.text=element_text(size=12))
    ggsave(filename=heatmap_genus,width=15,height=20)


}


main <- function() {
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    project_name = args[1]
    biom_tsv = args[2]
    metadata = args[3]
    threshold = args[4]
    distance = args[5]
    alpha_div_plots = args[6]
    barplot_relabund_phylum = args[7]
    barplot_relabund_family = args[8]
    barplot_relabund_genus = args[9]
    heatmap_class = args[10]
    heatmap_family = args[11]
    heatmap_genus = args[12]
    #Run alpha diversity calculations
    alphadiversity(biom_tsv, metadata, alpha_div_plots, barplot_relabund_phylum, barplot_relabund_family, barplot_relabund_genus, heatmap_class, heatmap_family, heatmap_genus, threshold, distance)

}
if (!interactive()) {
        main()
}

#### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
## ** End of the script **                                                 ####

#save.image(file=paste(SCRIPT,RData_name,sep=""))
#file.copy(from="../../../libs/Rscript_qiime2_snakemake.R",to=paste(SCRIPT,Rscript,sep=""))
