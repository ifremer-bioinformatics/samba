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
## Creation Date: 2020-04-02                                               ####
## Modified on: 2020-04-02                                                 ####
##                                                                           ##
## Emails: samba-sebimer@ifremer.fr                                        ####
##                                                                           ##
## Copyright (c) SeBiMER, april-2020                                       ####
## This program is free software: you can redistribute it and/or modify it   ##
## under the terms of the GNU Affero General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or         ##
## (at your option) any later version.                                       ## 
##                                                                           ##
## License at https://www.gnu.org/licenses/agpl-3.0.txt                      ##
##                                                                           ##
## This program is distributed in the hope that it will be useful, but       ##
## WITHOUT ANY WARRANTY; without even the implied warranty of                ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                      ##
## See the GNU Affero General Public License for more details.               ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs identification og all possible    ##
##        logical relationships between several sets of present in your      ##
##        ASV table based on your specified variable                         ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Load up the needed packages ####
requiredPackages = c("phyloseq","UpSetR","svglite","dplyr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

# @@@@@@@@@@@@ #
#              #
# Set analysis #
#              #
# @@@@@@@@@@@@ #
args = commandArgs(trailingOnly=TRUE)

PHYLOSEQ = readRDS(args[1])
PHYLOSEQ_merged = merge_samples(PHYLOSEQ, args[2])

sets = data.frame(t(otu_table(PHYLOSEQ_merged)))
sets[sets > 0] <- 1

PHYLOSEQ_relabund = transform_sample_counts(PHYLOSEQ, function(x) x / sum(x) * 100)
PHYLOSEQ_melt = psmelt(PHYLOSEQ_relabund)
PHYLOSEQ_melt = PHYLOSEQ_melt[,c(1,3)]
PHYLOSEQ_melt[PHYLOSEQ_melt==0] <- NA
PHYLOSEQ_melt_value = PHYLOSEQ_melt[complete.cases(PHYLOSEQ_melt),]
PHYLOSEQ_melt_mean = data.frame(aggregate(Abundance ~ OTU, PHYLOSEQ_melt_value , mean))

sets_abund = merge(PHYLOSEQ_melt_mean,sets,by.x="OTU",by.y="row.names")
row.names(sets_abund)=sets_abund[,1]
final_sets = sets_abund %>% select (-OTU)

svglite(paste(args[3],".svg",sep=""), width = 12, height = 10)
upset(final_sets, nsets=ncol(final_sets), boxplot.summary = "Abundance",number.angles = 0, text.scale = c(1.6,1.3,1.6,1.3,1.3,1.3),point.size = 2, line.size = 1, mainbar.y.label = "Number of ASVs", sets.x.label = "Total ASV")
dev.off()

png(paste(args[3],".png",sep=""), res=150, width=1600, height = 1400)
upset(final_sets, nsets=ncol(final_sets), boxplot.summary = "Abundance",number.angles = 0, text.scale = c(1.6,1.3,1.6,1.3,1.3,1.3),point.size = 2, line.size = 1, mainbar.y.label = "Number of ASVs", sets.x.label = "Total ASV")
dev.off()