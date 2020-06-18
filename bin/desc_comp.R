#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script:  Descriptive comparisons with UpsetR                   ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("phyloseq","UpSetR","svglite","dplyr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

args = commandArgs(trailingOnly=TRUE)

PHYLOSEQ = readRDS(args[1])
PHYLOSEQ_merged = merge_samples(PHYLOSEQ, args[2])

sets = data.frame(t(otu_table(PHYLOSEQ_merged)))
sets[sets > 0] <- 1

svg(paste(args[3],".svg",sep=""), width = 12, height = 10)
upset(sets, nsets=ncol(sets), number.angles = 0, text.scale = c(1.6,1.3,1.6,1.3,1.3,1.3),point.size = 2, line.size = 1, mainbar.y.label = "Number of ASVs", sets.x.label = "Total ASV")
dev.off()

png(paste(args[3],".png",sep=""), res=150, width=1600, height = 1400)
upset(sets, nsets=ncol(sets),number.angles = 0, text.scale = c(1.6,1.3,1.6,1.3,1.3,1.3),point.size = 2, line.size = 1, mainbar.y.label = "Number of ASVs", sets.x.label = "Total ASV")
dev.off()

