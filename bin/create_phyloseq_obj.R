#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Automated Statistical Analyzes of Metabarcoding Data   ##
##                                                                           ##
###############################################################################

## Load up the packages needed ####
library("dplyr")
library("stringr")
library("phyloseq")
library("phangorn")

create_phyloseq_obj <- function(phyloseq_rds, biom_tsv, metadata, microDecon, control, tree, remove_sample, sample_to_remove) {
    #Input data
    rawASVtable = read.table(biom_tsv, h=T, sep="\t", dec=".", check.names=FALSE, quote="")
    METADATA = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)

    if (microDecon == "true") {
        control_list =  unlist(strsplit(control,","))
        METADATA = METADATA[!rownames(METADATA) %in% control_list, ]
        write.table(METADATA, metadata, col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)
    }

    if (remove_sample == "true") {
        sample_to_remove_list =  unlist(strsplit(sample_to_remove,","))
        rawASVtable = rawASVtable[, !colnames(rawASVtable) %in% sample_to_remove_list ]
        METADATA = METADATA[!rownames(METADATA) %in% sample_to_remove_list, ]
        write.table(rawASVtable, biom_tsv, col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)
        write.table(METADATA, metadata, col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)
    }

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
    abund = abund[,colSums(abund) > 0]
    ABUND = otu_table(abund,taxa_are_rows=TRUE)
    TAX = tax_table(tax)
    var_data = data.frame(METADATA[rownames(METADATA) %in% colnames(abund), ])
    colnames(var_data) = colnames(METADATA)
    rownames(var_data) = rownames(METADATA)[rownames(METADATA) %in% colnames(abund)]
    METADATA = var_data
    write.table(METADATA, metadata, col.names=TRUE, row.names=TRUE, sep="\t",quote=FALSE)
    METADATA_phyloseq = sample_data(METADATA)
    TREE = read_tree(tree)
    TREE_ROOTED = phangorn::midpoint(TREE)
    PHYLOSEQ = phyloseq(ABUND,TAX,METADATA_phyloseq, TREE_ROOTED)
    saveRDS(PHYLOSEQ, file=phyloseq_rds)
}

main <- function() {
    # Get arguments from RScript command line
    args = commandArgs(trailingOnly=TRUE)
    phyloseq_rds = args[1]
    biom_tsv = args[2]
    metadata = args[3]
    microDecon = args[4]
    control = args[[5]]
    tree = args[6]
    remove_sample = args[7]
    sample_to_remove = args[[8]]
    create_phyloseq_obj(phyloseq_rds, biom_tsv, metadata, microDecon, control, tree, remove_sample, sample_to_remove)
}

if (!interactive()) {
        main()
}
