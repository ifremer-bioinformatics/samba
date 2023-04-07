#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Automated Statistical Analyzes of Metabarcoding Data   ##
##                                                                           ##
###############################################################################

## Load up the packages needed ####
requiredPackages = c("dplyr","stringr","phyloseq")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

create_phyloseq_obj <- function(raw_asv_table, asv_table_phyloseq, raw_metadata, metadata_phyloseq, tree, phyloseq_rds) {
    #Input data
    raw_asv_table = read.table(raw_asv_table, h=T, sep="\t", dec=".", check.names=FALSE, quote="")
    write.table(raw_asv_table, asv_table_phyloseq, sep="\t", dec=",", col.names=T, row.names=F, quote=F)
    raw_metadata = read.table(raw_metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    final_metadata = data.frame(raw_metadata[rownames(raw_metadata) %in% colnames(raw_asv_table[,-length(raw_asv_table)]), ])
    write.table(final_metadata, metadata_phyloseq, col.names=NA, row.names=TRUE, sep="\t",quote=FALSE)

    #Reformatting of the input data
    abund = raw_asv_table[,1:length(raw_asv_table)-1]
    row.names(abund) = abund$ASV_ID
    abund = abund %>% select (-ASV_ID)
    tax = raw_asv_table[,c(1,length(names(raw_asv_table)))]
    tax = data.frame(tax$ASV_ID,do.call(rbind,list(str_split_fixed(tax$taxonomy, ";",7))))
    colnames(tax) = c("ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    row.names(tax) = tax$ASV_ID
    tax = as.matrix(tax %>% select (-ASV_ID))

    ## Construction of the phyloseq object ####
    abund = abund[,colSums(abund) > 0]
    ABUND = otu_table(abund,taxa_are_rows=TRUE)
    TAX = tax_table(tax)
    METADATA_phyloseq = sample_data(final_metadata)
    TREE = read_tree(tree)
    PHYLOSEQ = phyloseq(ABUND,TAX,METADATA_phyloseq, TREE)
    PHYLOSEQ = subset_taxa(PHYLOSEQ, Kingdom != "Unassigned")

    ## Reformat taxonomy
    not_formated_tax_table = data.frame(tax_table(PHYLOSEQ))
    not_formated_tax_table = data.frame(apply(not_formated_tax_table, 2, function(x) gsub(" ", "", x)))
    not_formated_tax_table = data.frame(apply(not_formated_tax_table, 2, function(x) gsub("uncultured$", "", x)))
    not_formated_tax_table = data.frame(apply(not_formated_tax_table, 2, function(x) gsub("Unknown_Family", "", x)))
    
    # Genus level
    if(length(not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Family!="",]$Genus) != 0) {
      not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Family!="",]$Genus <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Family!="",]$Family, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Order!="",]$Genus) != 0) {
      not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Order!="",]$Genus <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Order!="",]$Order, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Class!="",]$Genus) != 0) {
      not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Class!="",]$Genus <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Class!="",]$Class, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Phylum!="",]$Genus) != 0) {
      not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Phylum!="",]$Genus <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Genus=="" & not_formated_tax_table$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Kingdom=="Bacteria",]$Kingdom) > length(not_formated_tax_table[not_formated_tax_table$Kingdom=="Eukaryota",]$Kingdom)) {
      unknown = "Unknown Bacteria"
    } else {
      unknown = "Unknown Eukaryota"
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Genus=="",]$Genus) > 0) {
      not_formated_tax_table[not_formated_tax_table$Genus=="",]$Genus <- unknown
    }
    
    # Family level
    if(length(not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Order!="",]$Family) != 0) {
      not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Order!="",]$Family <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Order!="",]$Order, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Class!="",]$Family) != 0) {
      not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Class!="",]$Family <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Class!="",]$Class, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Phylum!="",]$Family) != 0) {
      not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Phylum!="",]$Family <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Family=="" & not_formated_tax_table$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Family=="",]$Family) > 0) {
      not_formated_tax_table[not_formated_tax_table$Family=="",]$Family <- unknown
    }
    
    # Order level
    if(length(not_formated_tax_table[not_formated_tax_table$Order=="" & not_formated_tax_table$Class!="",]$Order) != 0) {
      not_formated_tax_table[not_formated_tax_table$Order=="" & not_formated_tax_table$Class!="",]$Order <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Order=="" & not_formated_tax_table$Class!="",]$Class, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Order=="" & not_formated_tax_table$Phylum!="",]$Order) != 0) {
      not_formated_tax_table[not_formated_tax_table$Order=="" & not_formated_tax_table$Phylum!="",]$Order <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Order=="" & not_formated_tax_table$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Order=="",]$Order) > 0) {
      not_formated_tax_table[not_formated_tax_table$Order=="",]$Order <- unknown
    }
    
    # Class level
    if(length(not_formated_tax_table[not_formated_tax_table$Class=="" & not_formated_tax_table$Phylum!="",]$Class) != 0) {
      not_formated_tax_table[not_formated_tax_table$Class=="" & not_formated_tax_table$Phylum!="",]$Class <- paste("Unknown", not_formated_tax_table[not_formated_tax_table$Class=="" & not_formated_tax_table$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(not_formated_tax_table[not_formated_tax_table$Class=="",]$Class) > 0) {
      not_formated_tax_table[not_formated_tax_table$Class=="",]$Class <- unknown
    }
    
    # Phylum level
    if(length(not_formated_tax_table[not_formated_tax_table$Phylum=="",]$Phylum) > 0) {
      not_formated_tax_table[not_formated_tax_table$Phylum=="",]$Phylum <- unknown
    }
    
    # Import new tax table in the phyloseq object
    final_formated_tax_table = as.matrix(not_formated_tax_table)
    tax_table(PHYLOSEQ) = tax_table(final_formated_tax_table)

    ## Save final phyloseq object
    saveRDS(PHYLOSEQ, file=phyloseq_rds)

}

main <- function() {
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    raw_asv_table = args[1]
    asv_table_phyloseq = args[2]
    raw_metadata = args[3]
    metadata_phyloseq = args[4]
    tree = args[5]
    phyloseq_rds = args[6]
    create_phyloseq_obj(raw_asv_table, asv_table_phyloseq, raw_metadata, metadata_phyloseq, tree, phyloseq_rds)
}

if (!interactive()) {
        main()
}
