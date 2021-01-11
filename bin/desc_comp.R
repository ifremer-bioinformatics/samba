#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script:  Descriptive comparisons with UpsetR                   ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("phyloseq","UpSetR","svglite","dplyr","stringr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

args = commandArgs(trailingOnly=TRUE)

PHYLOSEQ = readRDS(args[1])
PHYLOSEQ_merged = merge_samples(PHYLOSEQ, args[2])

sets = data.frame(t(otu_table(PHYLOSEQ_merged)),check.names=F)
sets[sets > 0] <- 1
sets$Taxonomy = data.frame(tax_table(PHYLOSEQ_merged))[,args[4]]
sets$Taxonomy = as.character(sets$Taxonomy)
sets$Taxonomy[sets$Taxonomy==""] <- "Unknown"
sets$Taxonomy[sets$Taxonomy==" Ambiguous_taxa"] <- "Unknown"
sets$Taxonomy = str_replace(sets$Taxonomy, "Candidatus", "C.")
sets[] = lapply(sets, function(x) replace(x, grep("uncultured", x), "Unknown"))
cols = names(sets[,1:length(sets)-1])
sets[cols] <- sapply(sets[cols],as.numeric)
sets$Taxonomy = as.factor(sets$Taxonomy)

count_taxa = data.frame(summary(sets$Taxonomy))
colnames(count_taxa) = "count"
count_taxa$name=row.names(count_taxa)
count_taxa = count_taxa %>% arrange(desc(count))
row.names(count_taxa) = count_taxa$name

top_taxa = row.names(count_taxa)[1:10]

sets$Taxonomy[!(sets$Taxonomy %in% top_taxa)] <- NA
sets$Taxonomy = as.character(sets$Taxonomy)
sets[is.na(sets)] <- "Others"

color_vector = c("antiquewhite2","darkblue","chartreuse3","coral2","darkkhaki","dodgerblue","deeppink","red","seagreen1","orange")

i<-0
taxa_list<-list()
if ("Unknown" %in% top_taxa) {
    vectorUniqueValue <- c("Unknown","Others",unique(sets$Taxonomy[!(sets$Taxonomy %in% c("Unknown","Others"))]))
    colors = c("black","grey45",color_vector)
} else {
    vectorUniqueValue <- c("Others",unique(sets$Taxonomy[!(sets$Taxonomy %in% "Others")]))
    colors = c("grey45",color_vector)
}
while ( length(vectorUniqueValue)>0 ){
  i<-i+1
  taxa_list[[i]]<-list(query = elements, params = list("Taxonomy",as.character(vectorUniqueValue)), color = colors[i], active = T, query.names=as.character(vectorUniqueValue)[1])
  vectorUniqueValue<-vectorUniqueValue[-1]
}

svg(paste(args[3],".svg",sep=""), width = 12, height = 10)
upset(sets, queries = taxa_list, query.legend="top", nsets=ncol(sets), matrix.color = "blue", shade.color = "white",  main.bar.color = "black", sets.bar.color = "blue", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8),point.size = 3, line.size = 1.2, mainbar.y.label = "Number of ASVs by set", sets.x.label = "Total ASV by condition")
dev.off()

png(paste(args[3],".png",sep=""), res=150, width=1600, height = 1400)
upset(sets, queries = taxa_list, query.legend="top", nsets=ncol(sets), matrix.color = "blue", shade.color = "white",  main.bar.color = "black", sets.bar.color = "blue", number.angles=0, text.scale = c(1.6,1.5,1.6,1.5,1.5,1.8),point.size = 3, line.size = 1.2, mainbar.y.label = "Number of ASVs by set", sets.x.label = "Total ASV by condition")
dev.off()
