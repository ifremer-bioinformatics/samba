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
##         Bioinformatics engineers                                          ##
##         SeBiMER, Ifremer                                                  ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2020-01-10                                                 ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## Emails: cyril.noel@ifremer.fr and laure.quintric@ifremer.fr             ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script contains all functions used for make       ##
##        beta diversity plots                                               ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

plot.nmds <- function(PHYLOSEQ, ord_nmds, criteria, color_samples, adonis_result, nmds, distance, width, height, graph_title) {
    ## Sample ordination - NMDS ####
    plot_ordination(PHYLOSEQ,ord_nmds,type="samples",color=criteria, title=graph_title) +
      theme_classic() +
      geom_point(size=3) +
      geom_text(aes(label=rownames(sample_data(PHYLOSEQ))),col="black",size=2.5,vjust=2,hjust=1) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      labs(color=criteria) +
      theme(axis.text=element_text(size=12,color="black")) +
      scale_fill_manual(values=alpha(color_samples,0.4)) +
      scale_color_manual(values=color_samples) +
      stat_ellipse(geom="polygon",alpha=0.1,type="t",aes_string(fill=criteria)) +
      labs(caption = paste("Stress:",round(ord_nmds$stress,4),
                           "\nAdonis statistic R:",round(adonis_result$aov.tab$R2[1]*100,2),
                     paste("\nAdonis based on ", criteria,": p-value"),adonis_result$aov.tab$`Pr(>F)`[1],sep=" "))
    ggsave(filename=paste(nmds,"_",distance,".svg",sep=""), device="svg", width = width, height = height)
    ggsave(filename=paste(nmds,"_",distance,".png",sep=""), device="png", width = width, height = height)
}

plot.pcoa <- function(PHYLOSEQ, ord_pcoa, criteria, color_samples, adonis_result, pcoa, distance, width, height, graph_title) {
    ## Sample ordination - MDS-PCoA ####
    plot_ordination(PHYLOSEQ,ord_pcoa,type="samples",color=criteria, title=graph_title) +
      theme_classic() +
      geom_point(size=3) +
      geom_text(aes(label=rownames(sample_data(PHYLOSEQ))),col="black",size=2.5,vjust=2,hjust=1) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      labs(color=criteria) +
      theme(axis.text=element_text(size=12,color="black")) +
      scale_fill_manual(values=alpha(color_samples,0.4)) +
      scale_color_manual(values=color_samples) +
      stat_ellipse(geom="polygon",alpha=0.1,type="t",aes_string(fill=criteria)) +
      labs(caption = paste("\nAdonis statistic R:",round(adonis_result$aov.tab$R2[1]*100,2),
                     paste("\nAdonis based on ", criteria,": p-value"),adonis_result$aov.tab$`Pr(>F)`[1],sep=" "))
    ggsave(filename=paste(pcoa,"_",distance,".svg",sep=""), device="svg", width = width, height = height)
    ggsave(filename=paste(pcoa,"_",distance,".png",sep=""), device="png", width = width, height = height)
}

plot.hc <- function(dendro, group, cols, col_group, method_hc, plot_hc, distance, width, height) {
    ## Sort GROUP color palette according to dend ####
    color = col_group[order.dendrogram(dendro)]
    ## Plot dendrogram ####
    svglite(paste(plot_hc,"_",distance,".svg",sep=""), width = width, height = height)
    plot = dendro %>% set("labels_colors", color) %>% plot(main = paste("Hierarchical clustering with the", method_hc, "method", "based on", distance, "distance", sep=" "))
    legend("topright", legend = levels(group), fill = cols, cex = 0.8, horiz=FALSE, border="white",box.lty=0)
    dev.off()
    png(filename=paste(plot_hc,"_",distance,".png",sep=""), res=150, width = 2000, height = 1200)
    plot = dendro %>% set("labels_colors", color) %>% plot(main = paste("Hierarchical clustering with the", method_hc, "method", "based on", distance, "distance", sep=" "))
    legend("topright", legend = levels(group), fill = cols, cex = 0.8, horiz=FALSE, border="white", box.lty=0)
    dev.off()

}

plot.pie <- function(ExpVar_perc, labels, distance, plot_pie, width, height) {
    svglite(paste(plot_pie,distance,".svg",sep=""), width = width, height = height)
    pie(ExpVar_perc,
      labels=labels,
      clockwise=TRUE,
      radius=1,
      col=brewer.pal(length(ExpVar_perc),"Set1"),
      border="white",
      cex=1,
      main=paste("Percentage of variance explained by each variable for",distance,"matrix",sep=" "))
    dev.off()
    png(filename=paste(plot_pie,distance,".png",sep=""), res=150, width = 2000, height = 1200)
    pie(ExpVar_perc,
      labels=labels,
      clockwise=TRUE,
      radius=1,
      col=brewer.pal(length(ExpVar_perc),"Set1"),
      border="white",
      cex=1,
      main=paste("Percentage of variance explained by each variable for",distance,"matrix",sep=" "))
    dev.off()
}
