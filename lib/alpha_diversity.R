#!/usr/bin/env Rscript

###############################################################################
##                                                                           ##
## Script name: R script for the SAMBA-nextflow workflow                   ####
##                                                                           ##
## Purpose of script: Automated Statistical Analyses of Metabarcoding Data   ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2020-01-24                                                 ####
##                                                                           ##
## Email: samba-sebimer@ifremer.fr                                         ####
## 									     ##
## Copyright (c) SeBiMER, august-2019                                      ####
## This program is free software: you can redistribute it and/or modify it   ##
## under the terms of the GNU Affero General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or         ##
## (at your option) any later version.                                       ## 
##									     ##
## License at https://www.gnu.org/licenses/agpl-3.0.txt 		     ##
## 									     ##
## This program is distributed in the hope that it will be useful, but       ##
## WITHOUT ANY WARRANTY; without even the implied warranty of                ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                      ##
## See the GNU Affero General Public License for more details.               ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs only the alpha diversity and      ##
##        taxonomic diversity based on your data                             ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Load up the needed packages ####
requiredPackages = c("dplyr","stringr","ggplot2","svglite","vegan","RColorBrewer","tidyr","gridExtra","egg","reshape2","BiocManager","phyloseq", "microbiome")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#						   #
# Function to standardize alpha diversity analysis #
#						   #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

alphadiversity <- function(PHYLOSEQ, alpha_div_plots, index_significance_tests, barplot_phylum, barplot_class, barplot_order, barplot_family, barplot_genus, kingdom, taxa-nb, distance, group, plot_rarefaction) {
    
    # ~~~~~~~~~~~~~~~ #
    # Alpha diversity #
    # ~~~~~~~~~~~~~~~ #

    ## Calcul of diversity indexes (Observed, Chao1, Shannon, InvSimpson, Pielou) ####

    alpha_rich = estimate_richness(PHYLOSEQ,measures=c("Observed","Chao1","Shannon","InvSimpson"))
    
    evenness = evenness(PHYLOSEQ,"pielou")
    detach("package:microbiome", unload=TRUE)
    
    alpha_rich$Pielou = evenness$pielou
    df = data.frame(alpha_rich,sample_data(PHYLOSEQ))
    df2 = gather(df,key="Measure",value="Value",Observed,Chao1,Shannon,InvSimpson,Pielou)
    df2$Measure = factor(df2$Measure,levels=c("Observed","Chao1","InvSimpson","Shannon","Pielou"))
    
    ## Alpha diversity plots ####
    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))

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
    ggsave(filename=paste(alpha_div_plots,".svg",sep=""),final_alpha_plot, device="svg", width=14, height=14)
    ggsave(filename=paste(alpha_div_plots,".png",sep=""),final_alpha_plot, device="png", width=14, height=14)

    ## Rarefaction curve ####
    rarec <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                   label = TRUE, cols = c(rep('red', nrow(x) / 2), rep('blue', nrow(x) / 2)), ...) {
      tot <- rowSums(x)
      S <- specnumber(x)
      nr <- nrow(x)
      out <- lapply(seq_len(nr), function(i) {
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i]) 
          n <- c(n, tot[i])
        drop(rarefy(x[i, ], n))
      })
      Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
      Smax <- sapply(out, max)
      plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
           type = "n", ...)
      rare <- min(rowSums(x))
      abline(v = rare, lwd = 2, lty=2, col="red")
      for (ln in seq_len(length(out))) {
        N <- attr(out[[ln]], "Subsample")
        lines(N, out[[ln]], col = "black", ...)
      }
      if (label) {
        ordilabel(cbind(tot, S), labels = rownames(x), ...)
      }
      invisible(out)
    }

    svglite(paste(plot_rarefaction,".svg",sep=""), width = 14, height = 10)
    rarec(t(otu_table(PHYLOSEQ)),step=20, cex = 0.6)
    dev.off()

    png(filename=paste(plot_rarefaction,".png",sep=""), res=150, width = 2000, height = 1200)
    rarec(t(otu_table(PHYLOSEQ)),step=20, cex = 0.6)
    dev.off()

    ## Statistical significance of indexes  ####
    anova_data = cbind(sample_data(PHYLOSEQ), alpha_rich)
    index = c("Observed","Chao1","Shannon","InvSimpson","Pielou")
    anova_data$Depth = sample_sums(PHYLOSEQ)
    variables = c()

    for (var in sample_variables(PHYLOSEQ) ) {
      l = length(levels(as.factor(get_variable(PHYLOSEQ,var))))
      if(l > 1 && l < nsamples(PHYLOSEQ)){
        variables <- cbind(variables,var)
      }
    }

    variables = paste(sep=" + ", "Depth", paste(collapse =" + ", variables ))

    sink(file = index_significance_tests , type = "output")
    for (m in index){
      f <- paste(m," ~ ", variables)
      cat(sep = "", "###############################################################\n
                     #Perform ANOVA on ",m,", which effects are significant\n
                     anova.",m," <-aov( ",f,", anova_data)\n
                     summary(anova.",m,")\n")
      anova_res <- aov( as.formula(f), anova_data)
      res <- summary(anova_res)
      print(res)
      cat("\n\n")
    }
    sink()

    # ~~~~~~~~~~~~~~~~~~~ #
    # Taxonomic diversity #
    # ~~~~~~~~~~~~~~~~~~~ #

    ## Variable definition ####
    taxaSet1 = unlist(strsplit(kingdom, " "))
    color_bar = color_vector[1:taxa-nb]
    
    ## Barplot representation ####
    #### at the phylum level ####
    composition(PHYLOSEQ, "Kingdom", taxaSet1, "Phylum", taxa-nb, fill="Phylum", group, color_bar, barplot_phylum) 

    #### at the class level ####
    composition(PHYLOSEQ, "Kingdom", taxaSet1, "Class", taxa-nb, fill="Class", group, color_bar, barplot_class)

    #### at the order level ####
    composition(PHYLOSEQ, "Kingdom", taxaSet1, "Order", taxa-nb, fill="Order", group, color_bar, barplot_order)

    #### at the family level ####
    composition(PHYLOSEQ, "Kingdom", taxaSet1, "Family", taxa-nb, fill="Family", group, color_bar, barplot_family)

    #### at the genus level ####
    composition(PHYLOSEQ, "Kingdom", taxaSet1, "Genus", taxa-nb, fill="Genus", group, color_bar, barplot_genus)
}

# @@@@@@@@@@@@@@@@@@@@@@@@ #
#                          #
# Alpha diversity analyses #
#                          #
# @@@@@@@@@@@@@@@@@@@@@@@@ #

main <- function() {
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    PHYLOSEQ = readRDS(args[1])
    distance = args[2]
    alpha_div_plots = args[3]
    kingdom = args[4]
    taxa-nb = args[5]
    barplot_phylum = args[6]
    barplot_class = args[7]
    barplot_order = args[8]
    barplot_family = args[9]
    barplot_genus = args[10]
    # Get group variable an replace "-" by "_"
    group = str_replace(args[11], "-", "_")
    index_significance_tests = args[12]
    workflow_dir = args[13]
    plot_rarefaction = args[14]
    # Get plot_composition function
    if (!exists("composition", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/barplot_graph_functions.R")))
    # Run alpha diversity analyses
    alphadiversity(PHYLOSEQ, alpha_div_plots, index_significance_tests, barplot_phylum, barplot_class, barplot_order, barplot_family, barplot_genus, kingdom, taxa-nb, distance, group, plot_rarefaction)
}
if (!interactive()) {
        main()
}
