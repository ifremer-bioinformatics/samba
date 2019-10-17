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
## Notes: general functions for beta diversity graphics                      ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

#### _____________________________________________________________________ ####

## Load up the packages needed ####

sample_nmds <- function(PHYLOSEQ, ord, criteria, color_samples, anosim_result, samples_ordination_plot, width, height, graph_title) {
    ## /3\ Sample analysis ####
    plot_ordination(PHYLOSEQ,ord,type="samples",color=criteria, title=graph_title) +
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
      labs(caption = paste("Stress:",round(ord$stress,4),
          "\nANOSIM statistic R:",round(anosim_result$statistic,4),
          paste("\nAnosim based on ", criteria,": p-value"),anosim_result$signif,sep=" "))
    ggsave(filename=samples_ordination_plot,width=width,height=height)
}
