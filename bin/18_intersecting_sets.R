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
PHYLOSEQ = readRDS(args[1])

## Merge samples according to the user-selected variable ####
PHYLOSEQ_var = merge_samples(PHYLOSEQ, args[2])

## Group ASV at the different taxonomic level ####
if (args[3] == "silva") {
    PHYLOSEQ_var_Phylum = tax_glom(PHYLOSEQ_var, "Phylum")
}
PHYLOSEQ_var_Class = tax_glom(PHYLOSEQ_var, "Class")
PHYLOSEQ_var_Order = tax_glom(PHYLOSEQ_var, "Order")
PHYLOSEQ_var_Family = tax_glom(PHYLOSEQ_var, "Family")
PHYLOSEQ_var_Genus = tax_glom(PHYLOSEQ_var, "Genus")
PHYLOSEQ_var_Species = tax_glom(PHYLOSEQ_var, "Species")
if (args[3] == "pr2-4" || args[3] == "pr2-5") {
    PHYLOSEQ_var_Supergroup = tax_glom(PHYLOSEQ_var, "Supergroup")
    PHYLOSEQ_var_Division = tax_glom(PHYLOSEQ_var, "Division")
}
if (args[3] == "pr2-5") {
    PHYLOSEQ_var_Subdivision = tax_glom(PHYLOSEQ_var, "Subdivision")
}

## Data processing for UpSetR ####
UpSetR_data_asv = data.frame(t(otu_table(PHYLOSEQ_var)))
UpSetR_data_asv[UpSetR_data_asv > 0] <- 1
png(paste0(args[2],"_UpSetR_ASV.png"), res=150, width=1600, height = 1400)
upset(UpSetR_data_asv, nsets=ncol(UpSetR_data_asv) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label="Number of ASVs", sets.x.label = "Total ASV by condition")
dev.off()

if (args[3] == "silva") {
    UpSetR_data_Phylum = data.frame(t(otu_table(PHYLOSEQ_var_Phylum)))
    UpSetR_data_Phylum[UpSetR_data_Phylum > 0] <- 1
    png(paste0(args[2],"_UpSetR_Phylum.png"), res=150, width=1600, height = 1400)
    upset(UpSetR_data_Phylum, nsets=ncol(UpSetR_data_Phylum) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label="Number of Phylum", sets.x.label = "Total Phylum by condition")
    dev.off()
}

common_tax_level = c("Class","Order","Family","Genus", "Species")
for (rank in common_tax_level) {
    UpSetR_data_rank = data.frame(t(otu_table(eval(parse(text=paste0("PHYLOSEQ_var_", rank))))))
    UpSetR_data_rank[UpSetR_data_rank > 0] <- 1
    png(paste0(args[2],"_UpSetR_", rank, ".png"), res=150, width=1600, height = 1400)
    upset(UpSetR_data_rank, nsets=ncol(UpSetR_data_rank) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label=paste("Number of", rank, sep=" "), sets.x.label = paste("Total", rank, "by condition", sep=" "))
    dev.off()
}

if (args[3] == "pr2-4" || args[3] == "pr2-5") {
    pr2_tax_level = c("Supergroup","Division")
    for (rank in pr2_tax_level) {
        UpSetR_data_rank = data.frame(t(otu_table(eval(parse(text=paste0("PHYLOSEQ_var_", rank))))))
        UpSetR_data_rank[UpSetR_data_rank > 0] <- 1
        png(paste0(args[2],"_UpSetR_", rank, ".png"), res=150, width=1600, height = 1400)
        upset(UpSetR_data_rank, nsets=ncol(UpSetR_data_rank) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label=paste("Number of", rank, sep=" "), sets.x.label = paste("Total", rank, "by condition", sep=" "))
        dev.off()
    }
}

if (args[3] == "pr2-5") {
    UpSetR_data_Subdivision = data.frame(t(otu_table(PHYLOSEQ_var_Subdivision)))
    UpSetR_data_Subdivision[UpSetR_data_Subdivision > 0] <- 1
    png(paste0(args[2],"_UpSetR_Subdivision..png"), res=150, width=1600, height = 1400)
    upset(UpSetR_data_Subdivision, nsets=ncol(UpSetR_data_Subdivision) ,keep.order=T, matrix.color="#00609B", shade.color="#FAE500", main.bar.color="#00609B", sets.bar.color="#00609B", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8), point.size=3, line.size=1.2, mainbar.y.label="Number of Subdivision", sets.x.label = "Total Subdivision by condition")
    dev.off()
}
