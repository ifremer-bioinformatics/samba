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

create_phyloseq <- function(phyloseq_rds, nanopore_count_table, metadata, final_table_all_assignation, final_table_only_assigned, db, tax2filter) {
  
  # Input data
  nanopore_metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
  rawtable = read.table(nanopore_count_table, h=T, sep="\t", dec=".", check.names=FALSE, quote="")
  rawtable = rawtable %>% select (c(-Identity, -Coverage))
  tax_to_filter = unlist(strsplit(tax2filter,","))
  if(length(tax_to_filter) > 0) {
      tax_to_filter_collapsed = paste(tax_to_filter, collapse='|')
      rawtable = rawtable[!grepl(tax_to_filter_collapsed, rawtable$Taxonomy),]
  }
  if(db == "silva") {
    rawtable = data.frame(rawtable[,1:(length(rawtable)-2)], rawtable[,length(rawtable)], do.call(rbind, list(str_split_fixed(rawtable$Taxonomy, ";",7))), check.names=FALSE)
    colnames(rawtable)[(length(rawtable)-7):length(rawtable)] = c("Assignation", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Reformat taxonomy
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub(" ", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^uncultured$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^unidentified_marine$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^uncultured_marine$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^unidentified_eubacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^unidentified$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^Unknown_Family$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^alpha_proteobacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^Alphaproteobacteria_bacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^bacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^bacterium_endosymbiont$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^bacterium_enrichment$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^beta_proteobacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^bioreactor_metagenome$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^candidate_division$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^delta_proteobacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^denitrifying_bacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^endosymbiont_of$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^gamma_proteobacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^groundwater_metagenome$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^gut_metagenome$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^human_gut$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^iron-reducing_bacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^iron-reducing_enrichment$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^metagenome$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^microbial_mat$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^mine_drainage$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^museum_specimen$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^nitrogen-fixing_bacterium$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^primary_endosymbiont$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^rock_porewater$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^saltmarsh_clone$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^synthetic_construct$", "", x)), check.names=FALSE)
  rawtable = data.frame(apply(rawtable, 2, function(x) gsub("^Incertae_Sedis$", "", x)), check.names=FALSE)
  rawtable$Species = str_replace_all(rawtable$Species, "^wastewater_metagenome$", "")
  if(length(rawtable[rawtable$Kingdom=="Bacteria",]$Kingdom) > length(rawtable[rawtable$Kingdom=="Eukaryota",]$Kingdom)) {
    unknown = "Unknown Bacteria"
  } else {
    unknown = "Unknown Eukaryota"
  }
  
    ## Species level
    if(length(rawtable[rawtable$Species=="" & rawtable$Genus!="",]$Species) != 0) {
      rawtable[rawtable$Species=="" & rawtable$Genus!="",]$Species <- paste("Unknown", rawtable[rawtable$Species=="" & rawtable$Genus!="",]$Genus, sep=" ")
    }
    if(length(rawtable[rawtable$Species=="" & rawtable$Family!="",]$Species) != 0) {
      rawtable[rawtable$Species=="" & rawtable$Family!="",]$Species <- paste("Unknown", rawtable[rawtable$Species=="" & rawtable$Family!="",]$Family, sep=" ")
    }
    if(length(rawtable[rawtable$Species=="" & rawtable$Order!="",]$Species) != 0) {
      rawtable[rawtable$Species=="" & rawtable$Order!="",]$Species <- paste("Unknown", rawtable[rawtable$Species=="" & rawtable$Order!="",]$Order, sep=" ")
    }
    if(length(rawtable[rawtable$Species=="" & rawtable$Class!="",]$Species) != 0) {
      rawtable[rawtable$Species=="" & rawtable$Class!="",]$Species <- paste("Unknown", rawtable[rawtable$Species=="" & rawtable$Class!="",]$Class, sep=" ")
    }
    if(length(rawtable[rawtable$Species=="" & rawtable$Phylum!="",]$Species) != 0) {
      rawtable[rawtable$Species=="" & rawtable$Phylum!="",]$Species <- paste("Unknown", rawtable[rawtable$Species=="" & rawtable$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(rawtable[rawtable$Species=="",]$Species) > 0) {
      rawtable[rawtable$Species=="",]$Species <- unknown
    }
  
    ## Genus level
    if(length(rawtable[rawtable$Genus=="" & rawtable$Family!="",]$Genus) != 0) {
      rawtable[rawtable$Genus=="" & rawtable$Family!="",]$Genus <- paste("Unknown", rawtable[rawtable$Genus=="" & rawtable$Family!="",]$Family, sep=" ")
    }
    if(length(rawtable[rawtable$Genus=="" & rawtable$Order!="",]$Genus) != 0) {
      rawtable[rawtable$Genus=="" & rawtable$Order!="",]$Genus <- paste("Unknown", rawtable[rawtable$Genus=="" & rawtable$Order!="",]$Order, sep=" ")
    }
    if(length(rawtable[rawtable$Genus=="" & rawtable$Class!="",]$Genus) != 0) {
      rawtable[rawtable$Genus=="" & rawtable$Class!="",]$Genus <- paste("Unknown", rawtable[rawtable$Genus=="" & rawtable$Class!="",]$Class, sep=" ")
    }
    if(length(rawtable[rawtable$Genus=="" & rawtable$Phylum!="",]$Genus) != 0) {
      rawtable[rawtable$Genus=="" & rawtable$Phylum!="",]$Genus <- paste("Unknown", rawtable[rawtable$Genus=="" & rawtable$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(rawtable[rawtable$Genus=="",]$Genus) > 0) {
      rawtable[rawtable$Genus=="",]$Genus <- unknown
    }
  
    ## Family level
    if(length(rawtable[rawtable$Family=="" & rawtable$Order!="",]$Family) != 0) {
      rawtable[rawtable$Family=="" & rawtable$Order!="",]$Family <- paste("Unknown", rawtable[rawtable$Family=="" & rawtable$Order!="",]$Order, sep=" ")
    }
    if(length(rawtable[rawtable$Family=="" & rawtable$Class!="",]$Family) != 0) {
      rawtable[rawtable$Family=="" & rawtable$Class!="",]$Family <- paste("Unknown", rawtable[rawtable$Family=="" & rawtable$Class!="",]$Class, sep=" ")
    }
    if(length(rawtable[rawtable$Family=="" & rawtable$Phylum!="",]$Family) != 0) {
      rawtable[rawtable$Family=="" & rawtable$Phylum!="",]$Family <- paste("Unknown", rawtable[rawtable$Family=="" & rawtable$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(rawtable[rawtable$Family=="",]$Family) > 0) {
      rawtable[rawtable$Family=="",]$Family <- unknown
    }
  
    ## Order level
    if(length(rawtable[rawtable$Order=="" & rawtable$Class!="",]$Order) != 0) {
      rawtable[rawtable$Order=="" & rawtable$Class!="",]$Order <- paste("Unknown", rawtable[rawtable$Order=="" & rawtable$Class!="",]$Class, sep=" ")
    }
    if(length(rawtable[rawtable$Order=="" & rawtable$Phylum!="",]$Order) != 0) {
      rawtable[rawtable$Order=="" & rawtable$Phylum!="",]$Order <- paste("Unknown", rawtable[rawtable$Order=="" & rawtable$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(rawtable[rawtable$Order=="",]$Order) > 0) {
      rawtable[rawtable$Order=="",]$Order <- unknown
    }
  
    ## Class level
    if(length(rawtable[rawtable$Class=="" & rawtable$Phylum!="",]$Class) != 0) {
      rawtable[rawtable$Class=="" & rawtable$Phylum!="",]$Class <- paste("Unknown", rawtable[rawtable$Class=="" & rawtable$Phylum!="",]$Phylum, sep=" ")
    }
    if(length(rawtable[rawtable$Class=="",]$Class) > 0) {
      rawtable[rawtable$Class=="",]$Class <- unknown
    }
  
    ## Phylum level
    if(length(rawtable[rawtable$Phylum=="",]$Phylum) > 0) {
      rawtable[rawtable$Phylum=="",]$Phylum <- unknown
    }

    ## Kingdom level
    if(length(rawtable[rawtable$Kingdom=="",]$Kingdom) > 0) {
      if(unknown == "Unknown Bacteria") {
        rawtable[rawtable$Kingdom=="",]$Kingdom <- "Bacteria"
      } else {
        rawtable[rawtable$Kingdom=="",]$Kingdom <- "Eukaryota"
      }
    }

    ## Manage uncultured species
    if(length(rawtable[grepl("uncultured", rawtable$Species),]$Species) > 0) {
      rawtable[grepl("uncultured", rawtable$Species),]$Species <- paste("Uncultured", rawtable[grepl("uncultured", rawtable$Species),]$Genus, sep=" ")
    }
    rawtable = data.frame(apply(rawtable, 2, function(x) gsub("Uncultured Unknown", "Uncultured", x)), check.names=FALSE)
  } else if(db == "pr2-4") {
    rawtable = data.frame(rawtable[,1:(length(rawtable)-2)], rawtable[,length(rawtable)], do.call(rbind, list(str_split_fixed(rawtable$Taxonomy, ";",8))), check.names=FALSE)
    colnames(rawtable)[(length(rawtable)-8):length(rawtable)] = c("Assignation", "Domain", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")
  } else if(db == "pr2-5") {
    rawtable = data.frame(rawtable[,1:(length(rawtable)-2)], rawtable[,length(rawtable)], do.call(rbind, list(str_split_fixed(rawtable$Taxonomy, ";",9))), check.names=FALSE)
    colnames(rawtable)[(length(rawtable)-9):length(rawtable)] = c("Assignation", "Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
  }
  
  # Subset data according to assignation category
  rawtable_all_assignation = rawtable %>% select(-Assignation)
  write.table(rawtable_all_assignation, final_table_all_assignation, sep="\t", dec=",", col.names=T, row.names=T, quote=F)
  
  rawtable_only_assigned = subset(rawtable, Assignation == "Assigned")
  rawtable_only_assigned = rawtable_only_assigned %>% select(-Assignation)
  write.table(rawtable_only_assigned, final_table_only_assigned, sep="\t", dec=",", col.names=T, row.names=T, quote=F)
  
  # Reformat data for phyloseq
  nanopore_abund_all_assignation = rawtable_all_assignation
  row.names(nanopore_abund_all_assignation) = nanopore_abund_all_assignation$Read_id
  if(db == "silva") {
    nanopore_abund_all_assignation = nanopore_abund_all_assignation %>% select(-c(Read_id,Kingdom,Phylum,Class,Order,Family,Genus,Species))
  } else if(db == "pr2-4") {
  nanopore_abund_all_assignation = nanopore_abund_all_assignation %>% select(-c(Read_id,Domain,Supergroup,Division,Class,Order,Family,Genus,Species))
  } else if(db == "pr2-5") {
    nanopore_abund_all_assignation = nanopore_abund_all_assignation %>% select(-c(Read_id,Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species))
  }
  nanopore_abund_all_assignation = data.frame(row.names=rownames(nanopore_abund_all_assignation), lapply(nanopore_abund_all_assignation, as.numeric), check.names=FALSE)
  if(db == "silva") {
    nanopore_tax_all_assignation = rawtable_all_assignation %>% select(c(Read_id,Kingdom,Phylum,Class,Order,Family,Genus,Species))
  } else if(db == "pr2-4") {
    nanopore_tax_all_assignation = rawtable_all_assignation %>% select(c(Read_id,Domain,Supergroup,Division,Class,Order,Family,Genus,Species))
  } else if(db == "pr2-5") {
    nanopore_tax_all_assignation = rawtable_all_assignation %>% select(c(Read_id,Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species))
  }
  row.names(nanopore_tax_all_assignation) = nanopore_tax_all_assignation$Read_id
  nanopore_tax_all_assignation = nanopore_tax_all_assignation %>% select(-Read_id)
  nanopore_tax_all_assignation = as.matrix(nanopore_tax_all_assignation)
  
  nanopore_abund_only_assigned = rawtable_only_assigned
  row.names(nanopore_abund_only_assigned) = nanopore_abund_only_assigned$Read_id
  if(db == "silva") {
    nanopore_abund_only_assigned = nanopore_abund_only_assigned %>% select(-c(Read_id,Kingdom,Phylum,Class,Order,Family,Genus,Species))
  } else if(db == "pr2-4") {
    nanopore_abund_only_assigned = nanopore_abund_only_assigned %>% select(-c(Read_id,Domain,Supergroup,Division,Class,Order,Family,Genus,Species))
  } else if(db == "pr2-5") {
    nanopore_abund_only_assigned = nanopore_abund_only_assigned %>% select(-c(Read_id,Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species))
  }
  nanopore_abund_only_assigned = data.frame(row.names=rownames(nanopore_abund_only_assigned), lapply(nanopore_abund_only_assigned, as.numeric), check.names=FALSE)
  if(db == "silva") {
    nanopore_tax_only_assigned = rawtable_only_assigned %>% select(c(Read_id,Kingdom,Phylum,Class,Order,Family,Genus,Species))
  } else if(db == "pr2-4") {
    nanopore_tax_only_assigned = rawtable_only_assigned %>% select(c(Read_id,Domain,Supergroup,Division,Class,Order,Family,Genus,Species))
  } else if(db == "pr2-5") {
    nanopore_tax_only_assigned = rawtable_only_assigned %>% select(c(Read_id,Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species))
  } 
  row.names(nanopore_tax_only_assigned) = nanopore_tax_only_assigned$Read_id
  nanopore_tax_only_assigned = nanopore_tax_only_assigned %>% select(-Read_id)
  nanopore_tax_only_assigned = as.matrix(nanopore_tax_only_assigned)
  
  # Construction of the phyloseq object
  NANOPORE_ABUND_ALL_ASSIGNATION = otu_table(nanopore_abund_all_assignation, taxa_are_rows=TRUE)
  NANOPORE_TAX_ALL_ASSIGNATION = tax_table(nanopore_tax_all_assignation)
  
  NANOPORE_ABUND_ONLY_ASSIGNED = otu_table(nanopore_abund_only_assigned, taxa_are_rows=TRUE)
  NANOPORE_TAX_ONLY_ASSIGNED = tax_table(nanopore_tax_only_assigned)
  
  NANOPORE_METADATA = sample_data(nanopore_metadata)
  
  NANOPORE_PHYLOSEQ_ALL_ASSIGNATION = phyloseq(NANOPORE_ABUND_ALL_ASSIGNATION, NANOPORE_TAX_ALL_ASSIGNATION, NANOPORE_METADATA)
  saveRDS(NANOPORE_PHYLOSEQ_ALL_ASSIGNATION, file=paste(phyloseq_rds,"all_assignation.rds", sep=""))
  
  NANOPORE_PHYLOSEQ_ONLY_ASSIGNED = phyloseq(NANOPORE_ABUND_ONLY_ASSIGNED, NANOPORE_TAX_ONLY_ASSIGNED, NANOPORE_METADATA)
  saveRDS(NANOPORE_PHYLOSEQ_ONLY_ASSIGNED, file=paste(phyloseq_rds,"only_assigned.rds", sep=""))

}

main <- function() {
  # Get arguments from RScript command line
  args = commandArgs(trailingOnly=TRUE)
  phyloseq_rds = args[1]
  nanopore_count_table = args[2]
  metadata = args[3]
  final_table_all_assignation = args[4]
  final_table_only_assigned = args[5]
  db = args[6]
  tax2filter = args[7]
  create_phyloseq(phyloseq_rds, nanopore_count_table, metadata, final_table_all_assignation, final_table_only_assigned, db, tax2filter)
}

if (!interactive()) {
  main()
}
