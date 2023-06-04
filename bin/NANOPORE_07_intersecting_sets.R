#!/usr/bin/env Rscript
################################################################################
##                                                                            ##
## Purpose of script: Visualize set intersections between multiple conditions ##
##                    of a tested variable                                    ##
##                                                                            ##
################################################################################

## Get arguments from RScript command line ####
args = commandArgs(trailingOnly=TRUE)

## Load up the needed packages ####
requiredPackages = c("phyloseq","UpSetR")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

## Phyloseq object import ####
if (args[2] == "silva") {
    PHYLOSEQ_Phylum = readRDS(paste0("phyloseq_", args[1], "_Phylum.rds"))
    PHYLOSEQ_var_Phylum = merge_samples(PHYLOSEQ_Phylum, args[3])
}

PHYLOSEQ_Class = readRDS(paste0("phyloseq_", args[1], "_Class.rds"))
PHYLOSEQ_var_Class = merge_samples(PHYLOSEQ_Class, args[3])

PHYLOSEQ_Order = readRDS(paste0("phyloseq_", args[1], "_Order.rds"))
PHYLOSEQ_var_Order = merge_samples(PHYLOSEQ_Order, args[3])

PHYLOSEQ_Family = readRDS(paste0("phyloseq_", args[1], "_Family.rds"))
PHYLOSEQ_var_Family = merge_samples(PHYLOSEQ_Family, args[3])

PHYLOSEQ_Genus = readRDS(paste0("phyloseq_", args[1], "_Genus.rds"))
PHYLOSEQ_var_Genus = merge_samples(PHYLOSEQ_Genus, args[3])

PHYLOSEQ_Species = readRDS(paste0("phyloseq_", args[1], "_Species.rds"))
PHYLOSEQ_var_Species = merge_samples(PHYLOSEQ_Species, args[3])

if (args[2] == "pr2-4" || args[2] == "pr2-5") {

    PHYLOSEQ_Supergroup = readRDS(paste0("phyloseq_", args[1], "_Supergroup.rds"))
    PHYLOSEQ_var_Supergroup = merge_samples(PHYLOSEQ_Supergroup, args[3])

    PHYLOSEQ_Division = readRDS(paste0("phyloseq_", args[1], "_Division.rds"))
    PHYLOSEQ_var_Division = merge_samples(PHYLOSEQ_Division, args[3])

}

if (args[2] == "pr2-5") {

    PHYLOSEQ_Subdivision = readRDS(paste0("phyloseq_", args[1], "_Subdivision.rds"))
    PHYLOSEQ_var_Subdivision = merge_samples(PHYLOSEQ_Subdivision, args[3])

}

## Data processing for UpSetR ####

if (args[2] == "silva") {
    UpSetR_data_Phylum = data.frame(t(otu_table(PHYLOSEQ_var_Phylum)))
    UpSetR_data_Phylum[UpSetR_data_Phylum > 0] <- 1
    png(paste0(args[3],"_UpSetR_Phylum.png"), res=150, width=1600, height = 1400)
    upset(UpSetR_data_Phylum, nsets=ncol(UpSetR_data_Phylum) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label="Number of Phylum", sets.x.label = "Total Phylum by condition")
    dev.off()
}

common_tax_level = c("Class","Order","Family","Genus", "Species")
for (rank in common_tax_level) {
    UpSetR_data_rank = data.frame(t(otu_table(eval(parse(text=paste0("PHYLOSEQ_var_", rank))))))
    UpSetR_data_rank[UpSetR_data_rank > 0] <- 1
    png(paste0(args[3],"_UpSetR_", rank, ".png"), res=150, width=1600, height = 1400)
    upset(UpSetR_data_rank, nsets=ncol(UpSetR_data_rank) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label=paste("Number of", rank, sep=" "), sets.x.label = paste("Total", rank, "by condition", sep=" "))
    dev.off()
}

if (args[2] == "pr2-4" || args[3] == "pr2-5") {
    pr2_tax_level = c("Supergroup","Division")
    for (rank in pr2_tax_level) {
        UpSetR_data_rank = data.frame(t(otu_table(eval(parse(text=paste0("PHYLOSEQ_var_", rank))))))
        UpSetR_data_rank[UpSetR_data_rank > 0] <- 1
        png(paste0(args[3],"_UpSetR_", rank, ".png"), res=150, width=1600, height = 1400)
        upset(UpSetR_data_rank, nsets=ncol(UpSetR_data_rank) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label=paste("Number of", rank, sep=" "), sets.x.label = paste("Total", rank, "by condition", sep=" "))
        dev.off()
    }
}

if (args[2] == "pr2-5") {
    UpSetR_data_Subdivision = data.frame(t(otu_table(PHYLOSEQ_var_Subdivision)))
    UpSetR_data_Subdivision[UpSetR_data_Subdivision > 0] <- 1
    png(paste0(args[3],"_UpSetR_Subdivision..png"), res=150, width=1600, height = 1400)
    upset(UpSetR_data_Subdivision, nsets=ncol(UpSetR_data_Subdivision) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label="Number of Subdivision", sets.x.label = "Total Subdivision by condition")
    dev.off()
}
