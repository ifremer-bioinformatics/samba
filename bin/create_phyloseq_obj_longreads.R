#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Automated Statistical Analyzes of Metabarcoding Data   ##
##                                                                           ##
###############################################################################

## Load up the packages needed ####
requiredPackages = c("phyloseq","dplyr","stringr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

create_phyloseq_obj <- function(phyloseq_rds, nanopore_count_table, metadata, rank, kingdom, final_table) {
    #Input data
    rawtable = read.table(nanopore_count_table, h=T, sep="\t", dec=".", check.names=FALSE, quote="")
    rawtable = rawtable %>% select (c(-Assignation))
    rawtable$Taxonomy = str_split_fixed(rawtable$Taxonomy, ";",7)[,as.numeric(rank)]
    rawtable_rank = rawtable %>% select (-Read_id) %>% group_by(Taxonomy) %>% summarise_each(funs(sum))
    rawtable_rank[1,1] = "Unclassified"
    rawtable_rank = data.frame(rawtable_rank)
    write.table(rawtable_rank,final_table, sep="\t", dec=",", col.names=T,row.names=T,quote=F)
    nanopore_metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)

    #Reformatting of the input data
    nanopore_abund = rawtable_rank
    row.names(nanopore_abund) = nanopore_abund$Taxonomy
    nanopore_abund = nanopore_abund %>% select (-Taxonomy)
    nanopore_tax = data.frame(Kingdom=rep(kingdom,length(rawtable_rank$Taxonomy)),rawtable_rank$Taxonomy)
    if (as.numeric(rank) == 1) {
        rank_name = "Kingdom"
    } else if (as.numeric(rank) == 2) {
        rank_name = "Phylum"
    } else if (as.numeric(rank) == 3) {
        rank_name = "Class"
    } else if (as.numeric(rank) == 4) {
        rank_name = "Order"
    } else if (as.numeric(rank) == 5) {
        rank_name = "Family"
    } else if (as.numeric(rank) == 6) {
        rank_name = "Genus"
    } else {
        rank_name = "Species"
    }
    colnames(nanopore_tax) = c("Kingdom",rank_name)
    row.names(nanopore_tax) = nanopore_tax[[rank_name]]
    nanopore_tax = as.matrix(nanopore_tax)

    ## Construction of the phyloseq object ####
    nanopore_abund = nanopore_abund[,colSums(nanopore_abund) > 0]
    NANOPORE_ABUND = otu_table(nanopore_abund,taxa_are_rows=TRUE)
    NANOPORE_TAX = tax_table(nanopore_tax)
    nanopore_metadata = nanopore_metadata[rownames(nanopore_metadata) %in% colnames(nanopore_abund), ]
    write.table(nanopore_metadata, "metadata_stats.tsv", col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)
    NANOPORE_METADATA = sample_data(nanopore_metadata)
    NANOPORE_PHYLOSEQ = phyloseq(NANOPORE_ABUND,NANOPORE_TAX,NANOPORE_METADATA)
    saveRDS(NANOPORE_PHYLOSEQ, file=phyloseq_rds)
}

main <- function() {
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    phyloseq_rds = args[1]
    nanopore_count_table = args[2]
    metadata = args[3]
    rank = noquote(args[4])
    kingdom = args[5]
    final_table = args[6]
    create_phyloseq_obj(phyloseq_rds, nanopore_count_table, metadata, rank, kingdom, final_table)
}

if (!interactive()) {
        main()
}
