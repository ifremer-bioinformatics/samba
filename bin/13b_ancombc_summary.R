#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Summarise ANCOM-BC results                             ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("stringr","dplyr", "reshape2" ,"ggplot2", "svglite")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

args = commandArgs(trailingOnly=TRUE)

var_name = args[1]
lfc_file = args[2]
qvalue_file = args[3]
ancombc_summary_output= args[4]
metadata_file = args[5]
ancombc_heatmap = args[6]

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#                                #
# Merge ANCOM-BC multiple result #
# into a unique file             #
#                                #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

## Import ANCOM-BC results ####
lfc = read.csv(lfc_file, check.names=F)
colnames(lfc) = str_remove_all(colnames(lfc), var_name)
values_tested = colnames(lfc %>% select(-c(taxon,feature_id,`(Intercept)`)))
lfc = lfc[, c(2,4:length(lfc))]
colnames(lfc) = c("taxon", paste(colnames(lfc)[2:length(lfc)],"lfc",sep="_"))

padjvalue = read.csv(qvalue_file, check.names=F)
colnames(padjvalue) = str_remove_all(colnames(padjvalue), var_name)
padjvalue = padjvalue[, c(2,4:length(padjvalue))]
padjvalue_signif = padjvalue
colnames(padjvalue) = c("taxon", paste(colnames(padjvalue)[2:length(padjvalue)],"padjvalue",sep="_"))
rownames(padjvalue_signif) = padjvalue_signif$taxon
padjvalue_signif = padjvalue_signif %>% select(-taxon)
padjvalue_signif[padjvalue_signif <= 0.001] <- "***"
padjvalue_signif[(padjvalue_signif < 0.01) & (padjvalue_signif > 0.001)] <- "**"
padjvalue_signif[(padjvalue_signif < 0.05) & (padjvalue_signif > 0.01)] <- "*"
padjvalue_signif[padjvalue_signif > 0.05] <- "ns"
colnames(padjvalue_signif) = paste(colnames(padjvalue_signif), "signif", sep="_")
padjvalue = merge(padjvalue, padjvalue_signif, by.x="taxon", by.y="row.names")

## Merge files ####
merged_ancombc_results = merge(lfc, padjvalue, by="taxon")
merged_ancombc_results = merged_ancombc_results %>% select(taxon, order(colnames(merged_ancombc_results)))

## Save output ####
write.table(merged_ancombc_results, ancombc_summary_output, col.names=T, row.names=F, sep="\t", quote=F)

## Create a heatmap output ####
metadata = read.table(metadata_file, h=T, sep="\t")
values_var = unique(metadata[,var_name])
ref_value = setdiff(values_var,values_tested)
melted_merged_ancombc_results = melt(merged_ancombc_results)
lfc_melted_merged_ancombc_results = melted_merged_ancombc_results[grep("lfc", melted_merged_ancombc_results$variable),]
filtered_melted_merged_ancombc_results = subset(lfc_melted_merged_ancombc_results, lfc_melted_merged_ancombc_results[,2:(length(lfc_melted_merged_ancombc_results)-2)] != "ns")
filtered_melted_merged_ancombc_results = na.omit(filtered_melted_merged_ancombc_results)

for(var in values_tested) {
  assign(paste("table_",var, sep=""), filtered_melted_merged_ancombc_results[grep(var, filtered_melted_merged_ancombc_results$variable),])

}
for(df_id in ls(pattern="^table_")) {
  var_column_signif = paste(str_remove(df_id, "table_"), "_signif", sep="")
  df = get(noquote(df_id))
  assign(paste("signif_value_",df_id, sep=""), df %>% mutate(value_signif = case_when(df[var_column_signif] != 'ns' ~ df$value)))
}
for(signif_id in ls(pattern="^signif_value_table")) {
  var_column_signif = paste(str_remove(signif_id, "signif_value_table_"), "_signif", sep="")
  df = get(noquote(signif_id))
  assign(paste("final_signif_value_",signif_id, sep=""), df %>% mutate(sb_signif = case_when(df$value_signif != 'NA' ~ df[,var_column_signif])))
}
list_signif_table = lapply(ls(pattern="final_signif_value_"), get)
only_signif_filtered_melted_merged_ancombc_results = bind_rows(list_signif_table)

heatmap_y_label = str_replace_all(str_remove(unique(only_signif_filtered_melted_merged_ancombc_results$variable), "_lfc"), "_" , " ")
min_lfc = min(only_signif_filtered_melted_merged_ancombc_results$value)-1
max_lfc = max(only_signif_filtered_melted_merged_ancombc_results$value)+1
ggplot(only_signif_filtered_melted_merged_ancombc_results, aes(x=taxon, y=variable, fill=value_signif)) +
  coord_fixed() +
  geom_tile(color="white", lwd=1.5, linetype=1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  theme(plot.background = element_rect(fill="white", colour=NA)) +
  theme(panel.grid.major=element_blank(), panel.border=element_blank()) +
  theme(axis.text = element_text(color="black", size=12)) +
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle=70, hjust=1)) +
  labs(fill="LFC", title=paste("ANCOM-BC:", str_replace_all(var_name, "_", " "), "vs", str_replace_all(ref_value, "_", " "), sep=" ")) +
  scale_fill_gradient2(limits=c(min_lfc,max_lfc), high="navyblue", low="red", mid="white", midpoint=0, na.value="white") +
  scale_y_discrete(labels=heatmap_y_label) +
  geom_text(aes(label=sb_signif))
ggsave(paste(ancombc_heatmap, ".svg", sep=""), device="svg")
ggsave(paste(ancombc_heatmap, ".png", sep="" ), device="png")
