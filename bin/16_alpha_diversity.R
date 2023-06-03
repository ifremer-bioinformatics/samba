#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Alpha diversity analysis from ILLUMINA data            ##
##                                                                           ##
###############################################################################

## Get arguments from RScript command line ####
args = commandArgs(trailingOnly=TRUE)

## Load up the needed packages ####
requiredPackages = c("RColorBrewer", "phyloseq", "reshape2", "dplyr", "tidyr", "ggplot2","ggpubr", "stringr", "rstatix","grid", "vegan")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}
color_set = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))

## Functions ####

format_data_pie <- function(PHYLOSEQ, env_var, taxa, unknown, taxa_nb) {
  PHYLOSEQ_GROUP <- merge_samples(PHYLOSEQ, env_var)
  asvtab <- t(otu_table(PHYLOSEQ_GROUP))
  asvtab <- as(asvtab, "matrix")
  asvtab <- apply(asvtab, 2, function(x) x / sum(x))
  ## Reformat ASV table (using melt function of the reshape2 package)
  mdf <- melt(asvtab, varnames=c("ASV", env_var))
  colnames(mdf)[3] <- "Abundance"
  ## Add taxonomic information and replace NA and unclassified Unknown
  taxtab <- tax_table(PHYLOSEQ_GROUP)
  taxtab <- as(taxtab, "matrix")
  taxtab[, taxa][is.na(taxtab[, taxa])] <- unknown
  taxtab[, taxa][taxtab[, taxa] %in% c("", "unclassified", "Unclassified", "uncultured", "Uncultured")] <- unknown
  taxtab <- data.frame(ASV = rownames(taxtab), taxtab)
  mdf <- merge(mdf, taxtab, by="ASV")
  ## Aggregate at Genus level
  formula_mdf_agg <- paste(paste("Abundance", env_var, sep="~"), taxa, sep="+")
  mdf <- aggregate(as.formula(formula_mdf_agg), data=mdf, FUN=sum)
  formula_total_abund <- paste("Abundance", taxa, sep="~")
  total_abundance <- aggregate(as.formula(formula_total_abund), data=mdf, FUN=sum)
  ## Keep only taxa_nb top taxa and aggregate the rest as "Other"
  ordered_abund <- total_abundance[ order(total_abundance[, "Abundance"], decreasing=TRUE),]
  min_top_abund <- round(ordered_abund[10, "Abundance"] / sum(ordered_abund$Abundance) * 100, 2)
  list_ordered_abund <- as.character(ordered_abund[, taxa])
  list_ordered_abund <- list_ordered_abund[list_ordered_abund != unknown]
  top <- list_ordered_abund[1:min(length(list_ordered_abund), taxa_nb)]
  mdf[, taxa] <- as.character(mdf[ , taxa])
  ii <- (mdf[, taxa] %in% c(top, unknown))
  Others_lab <- paste("Others (", taxa, " <", min_top_abund, "%)", sep="")
  mdf[!ii , taxa] <- Others_lab
  mdf <- aggregate(as.formula(formula_mdf_agg), data=mdf, FUN=sum)
  mdf[, taxa] <- factor(mdf[, taxa], levels=c(sort(top), Others_lab, unknown))
  ## Sort the entries by abundance to produce nice stacked bars in ggplot2
  mdf <- mdf[ order(mdf[, taxa], mdf$Abundance, decreasing=TRUE), ]
  return(mdf)
}

ggplot_pie <- function(PHYLOSEQ, env_var, taxa, unknown, taxa_nb, color_pie, pie_data_filename) {
  formated_data <- format_data_pie(PHYLOSEQ, env_var, taxa, unknown, taxa_nb)
  formated_data$Abundance <- formated_data$Abundance*100
  Others_lab <- as.character(unique(formated_data[grepl("Others", formated_data[, taxa]), taxa]))
  ## Management of pie colors
  if (!is.null(taxa) && any(c(Others_lab, unknown) %in% unique(formated_data[, taxa]))) {
    taxa4color <- as.character(unique(formated_data[, taxa]))
    taxa4color <- taxa4color[ ! taxa4color %in% c(Others_lab, unknown)]
    colvals <- c(color_pie[1:taxa_nb], "grey45", "black")
    names(colvals) <- c(taxa4color, Others_lab, unknown)
  } else {
    taxa4color <- as.character(unique(formated_data[, taxa]))
    colvals <- color_pie[1:taxa_nb]
    names(colvals) <- taxa4color
  }
  group_formula <- c(env_var, taxa, "Abundance")
  formated_data <- formated_data %>% group_by_at(group_formula) %>% summarise()
  dcast_formula <- paste(taxa, env_var, sep="~")
  pie_data <- dcast(formated_data, as.formula(dcast_formula))
  write.table(pie_data, paste0("pie_data_", env_var, "_", taxa, ".tsv"), col.names=T, row.names=F, sep="\t", quote=F)
  
  ggplot(formated_data, aes_string(x=1, y="Abundance", fill=taxa)) +
    geom_bar(position="stack", stat="identity", width=2) +
    coord_polar("y", start=0) +
    theme_void() +
    facet_wrap(as.formula(paste("~", env_var))) +
    scale_fill_manual(values=colvals) +
    theme(strip.text.x = element_text(size=12, face="bold")) +
    theme(strip.text.y = element_text(size=12, face="bold")) +
    theme(legend.box.spacing = unit(1, "cm")) +
    theme(legend.position = "right") +
    theme(legend.title = element_text(size=14, face="bold")) +
    theme(legend.text = element_text(size=12)) +
    theme(plot.background = element_rect(fill="white", color="white")) +
    theme(panel.background = element_rect(fill="white", color="white"))
  ggsave(filename=paste("pie_",env_var,"_",taxa,".svg",sep=""), device="svg", width=14, height=10)
  ggsave(filename=paste("pie_",env_var,"_",taxa,".png",sep=""), device="png", width=14, height=10)
  
}

format_data_barplot <- function(PHYLOSEQ, taxa, unknown, taxa_nb) {
  asvtab <- otu_table(PHYLOSEQ)
  asvtab <- as(asvtab, "matrix")
  asvtab <- apply(asvtab, 2, function(x) x / sum(x))
  ## Reformat ASV table (using melt function of the reshape2 package)
  mdf <- melt(asvtab, varnames=c("ASV", "Samples"))
  colnames(mdf)[3] <- "Abundance"
  ## Add taxonomic information and replace NA and unclassified Unknown
  taxtab <- as(tax_table(PHYLOSEQ), "matrix")
  taxtab[, taxa][is.na(taxtab[, taxa])] <- unknown
  taxtab[, taxa][taxtab[, taxa] %in% c("", "unclassified", "Unclassified", "uncultured", "Uncultured")] <- unknown
  taxtab <- data.frame(ASV = rownames(taxtab), taxtab)
  mdf <- merge(mdf, taxtab, by="ASV")
  ## Aggregate at Genus level
  formula_mdf_agg <- paste("Abundance ~ Samples", taxa, sep="+")
  mdf <- aggregate(as.formula(formula_mdf_agg), data=mdf, FUN=sum)
  formula_total_abund <- paste("Abundance", taxa, sep="~")
  total_abundance <- aggregate(as.formula(formula_total_abund), data=mdf, FUN=sum)
  ## Keep only taxa_nb top class and aggregate the rest as "Other"
  ordered_abund <- total_abundance[ order(total_abundance[, "Abundance"], decreasing=TRUE),]
  min_top_abund <- round(ordered_abund[10, "Abundance"] / sum(ordered_abund$Abundance) * 100, 2)
  list_order_abund <- as.character(ordered_abund[, taxa])
  list_order_abund <- list_order_abund[list_order_abund != unknown]
  top <- list_order_abund[1:min(length(list_order_abund), taxa_nb)]
  mdf[, taxa] <- as.character(mdf[ , taxa])
  ii <- (mdf[, taxa] %in% c(top, unknown))
  Others_lab <- paste("Others (", taxa, " <", min_top_abund, "%)", sep="") 
  mdf[!ii , taxa] <- Others_lab
  mdf <- aggregate(as.formula(formula_mdf_agg), data=mdf, FUN=sum)
  mdf[, taxa] <- factor(mdf[, taxa], levels=c(sort(top), Others_lab, unknown))
  ## Add metadata
  metadata <- as(sample_data(PHYLOSEQ), "data.frame")
  metadata$Samples <- sample_names(PHYLOSEQ)
  mdf <- merge(mdf, metadata, by.x = "Samples")
  ## Sort the entries by abundance to produce nice stacked bars in ggplot
  mdf <- mdf[ order(mdf[ , taxa], mdf$Abundance, decreasing = TRUE), ]
  return(mdf)
}

ggplot_barplot <- function(PHYLOSEQ, env_var, taxa, unknown, taxa_nb, color_bar, abundance_data_filename) {
  formated_data <- format_data_barplot(PHYLOSEQ, taxa, unknown, taxa_nb)
  formated_data$Abundance <- formated_data$Abundance*100
  Others_lab <- as.character(unique(formated_data[grepl("Others", formated_data[, taxa]), taxa]))
  ## Management of pie colors
  if (!is.null(taxa) && any(c(Others_lab, unknown) %in% unique(formated_data[, taxa]))) {
    taxa4color <- as.character(unique(formated_data[, taxa]))
    taxa4color <- taxa4color[ ! taxa4color %in% c(Others_lab, unknown)]
    colvals <- c(color_bar[1:taxa_nb], "grey45", "black")
    names(colvals) <- c(taxa4color, Others_lab, unknown)
  } else {
    taxa4color <- as.character(unique(formated_data[, taxa]))
    colvals <- color_bar[1:taxa_nb]
    names(colvals) <- taxa4color
  }
  sample_abund <- formated_data[, c("Samples", taxa, "Abundance")]
  dcast_formula <- paste(taxa, "Samples", sep="~")
  abundance_data <- dcast(sample_abund, as.formula(dcast_formula))
  write.table(abundance_data, paste0("barplot_data_", env_var, "_", taxa, ".tsv"), col.names=T, row.names=F, sep="\t", quote=F)
  
  ggplot(formated_data, aes_string(x="Samples", y="Abundance", fill=taxa)) +
    geom_bar(position="stack", stat="identity", width=0.6) +
    theme_minimal() +
    facet_wrap(as.formula(paste("~", env_var)), scales = "free_x", nrow = 1) +
    scale_fill_manual(values=colvals) +
    theme(strip.text=element_text(size=18, color="#FAE500")) +
    theme(strip.background = element_rect(colour="#00609B", fill="#00609B", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size=12, face="bold")) +
    theme(strip.text.y = element_text(size=12, face="bold")) +
    theme(legend.box.spacing = unit(1, "cm")) +
    theme(legend.position = "right") +
    theme(legend.title = element_text(size=14, face="bold")) +
    theme(legend.text = element_text(size=12)) +
    theme(plot.background = element_rect(fill="white", color="white")) +
    theme(panel.background = element_rect(fill="white", color="white")) +
    theme(panel.grid = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(color="black", size=20)) +
    labs(x="", y="Relative abundance (%)") +
    theme(axis.text.x = element_text(color="black", size=14, angle=60, hjust=1, vjust=1.2)) +
    theme(axis.text.y = element_text(color="black", size=14))
  ggsave(filename=paste("barplot_",env_var,"_",taxa,".svg",sep=""), device="svg", width=14, height=10)
  ggsave(filename=paste("barplot_",env_var,"_",taxa,".png",sep=""), device="png", width=14, height=10)
  
}

# ~~~~~~~~~~~~~~~ #
# Alpha diversity #
# ~~~~~~~~~~~~~~~ #

## Phyloseq objects import ####
PHYLOSEQ = readRDS(args[1])

## Calcul of diversity indexes (Chao1, Shannon, InvSimpson) ####
alpha_rich = estimate_richness(PHYLOSEQ, measures=c("Observed", "Chao1", "ACE", "Shannon", "InvSimpson"))
evenness = microbiome::evenness(PHYLOSEQ, "pielou")
df_alpha_div = cbind(alpha_rich, evenness)

write.table(df_alpha_div, args[2], col.names=NA, row.names=TRUE, sep="\t", dec=",", quote=FALSE)

df_alpha_div_sdf = data.frame(df_alpha_div, sample_data(PHYLOSEQ))
df_alpha_div_sdf = gather(df_alpha_div_sdf, key="Measure", value="Value", Observed, Chao1, ACE, Shannon, InvSimpson, pielou)
df_alpha_div_sdf$Measure = factor(df_alpha_div_sdf$Measure, levels=c("Observed", "Chao1", "ACE", "Shannon", "InvSimpson", "pielou"))

## Alpha diversity plots according to user-selected variable ####
test_formula = paste("Value", args[3], sep="~")
var_value = length(unique(df_alpha_div_sdf[, args[3]]))
stats_var = df_alpha_div_sdf %>% group_by(Measure) %>% t_test(as.formula(test_formula)) %>% adjust_pvalue(method="bonferroni") %>% add_significance() %>% add_y_position(scales="free")
var_color = adjustcolor(color_set[2:(var_value+1)], alpha.f=0.8)

alpha_bxp_var = ggboxplot(df_alpha_div_sdf, x=args[3], y="Value", fill=args[3], facet.by="Measure", add="jitter", scales="free_y", add.params=list(size=1)) +
  stat_pvalue_manual(stats_var, hide.ns=TRUE, label="{signif(p.adj,1)}{p.adj.signif}", label.size=4, tip.length = 0.02) +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.1))) +
  scale_fill_manual(values=var_color) +
  labs(y="Index values", fill=args[3]) +
  theme_bw() +
  theme(axis.text = element_text(size=14, color="black")) +
  theme(axis.text.x = element_text(hjust=1, vjust=1, angle=60)) +
  theme(strip.text.x = element_text(size=16, face="bold", color="black")) +
  theme(strip.background = element_rect(fill="grey85")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=16, face="bold")) +
  theme(legend.title = element_text(size=16)) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.position = "bottom")

ggsave(filename=paste(args[4],"_", args[3],".svg", sep=""), alpha_bxp_var, device="svg", width=16, height=10)
ggsave(filename=paste(args[4],"_", args[3],".png", sep=""), alpha_bxp_var, device="png", width=16, height=10)

## Rarefaction curve ####
asv_matrix = as.matrix(t(data.frame(otu_table(PHYLOSEQ), check.names=F)))
metadata_df = data.frame(sample_data(PHYLOSEQ))
metadata_df$sample = as.factor(rownames(metadata_df))
rarefaction_curve = rarecurve(asv_matrix, step=20, tidy=TRUE)
rarefaction_curve_metadata = left_join(rarefaction_curve, metadata_df, by=c("Site" = "sample"))

rarefaction_value = min(rowSums(asv_matrix))

FINAL_rarefaction_curve = ggplot(rarefaction_curve_metadata) +
  geom_line(aes_string(x="Sample", y="Species", group="Site", colour=args[3])) +
  theme_bw() +
  facet_wrap(as.formula(paste("~", args[3])), nrow=3) +
  theme(strip.text.x = element_text(size=16, face="bold", color="black")) +
  theme(strip.background = element_rect(fill="grey85")) +
  scale_color_manual(values=color_set) +
  theme(axis.text = element_text(size=14, color="black")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=16, face="bold")) +
  geom_vline(xintercept=rarefaction_value,lty=2) +
  labs(y="Taxa count") +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0,0.1,0.1,0, "cm"))

ggsave(filename=paste(args[5],"_",args[3],".svg", sep=""), FINAL_rarefaction_curve, device="svg", width=14, height=12)
ggsave(filename=paste(args[5],"_",args[3],".png", sep=""), FINAL_rarefaction_curve, device="png", width=14, height=12)

# ~~~~~~~~~~~~~~~~~~~ #
# Taxonomy inspection #
# ~~~~~~~~~~~~~~~~~~~ #

## Pie graph with sample grouped by env_var 
tmp_df <- data.frame(tax_table(PHYLOSEQ))
if(args[7] == "silva") {
  if(length(tmp_df[tmp_df$Kingdom=="Bacteria",]$Kingdom) > length(tmp_df[tmp_df$Kingdom=="Eukaryota",]$Kingdom)) {
    unknown = "Unknown Bacteria"
  } else {
    unknown = "Unknown Eukaryota"
  }
}
if(args[7] == "pr2-4" || args[7] == "pr2-5") {
    unknown = "Unknown Eukaryota"
}

if(args[7] == "silva") {
  ggplot_pie(PHYLOSEQ, args[3], "Phylum", unknown, args[6], color_set)
}
ggplot_pie(PHYLOSEQ, args[3], "Class", unknown, args[6], color_set)
ggplot_pie(PHYLOSEQ, args[3], "Order", unknown, args[6], color_set)
ggplot_pie(PHYLOSEQ, args[3], "Family", unknown, args[6], color_set)
ggplot_pie(PHYLOSEQ, args[3], "Genus", unknown, args[6], color_set)
ggplot_pie(PHYLOSEQ, args[3], "Species", unknown, args[6], color_set)
if(args[7] == "pr2-4" || args[7] == "pr2-5") {
  ggplot_pie(PHYLOSEQ, args[3], "Supergroup", unknown, args[6], color_set)
  ggplot_pie(PHYLOSEQ, args[3], "Division", unknown, args[6], color_set)
}
if(args[7] == "pr2-5") {
  ggplot_pie(PHYLOSEQ, args[3], "Subdivision", unknown, args[6], color_set)
}

## Taxonomic barplot for all samples facetted by env_var
if(args[7] == "silva") {
  ggplot_barplot(PHYLOSEQ, args[3], "Phylum", unknown, args[6], color_set)
}
ggplot_barplot(PHYLOSEQ, args[3], "Class", unknown, args[6], color_set)
ggplot_barplot(PHYLOSEQ, args[3], "Order", unknown, args[6], color_set)
ggplot_barplot(PHYLOSEQ, args[3], "Family", unknown, args[6], color_set)
ggplot_barplot(PHYLOSEQ, args[3], "Genus", unknown, args[6], color_set)
ggplot_barplot(PHYLOSEQ, args[3], "Species", unknown, args[6], color_set)
if(args[7] == "pr2-4" || args[7] == "pr2-5") {
  ggplot_barplot(PHYLOSEQ, args[3], "Supergroup", unknown, args[6], color_set)
  ggplot_barplot(PHYLOSEQ, args[3], "Division", unknown, args[6], color_set)
}
if(args[7] == "pr2-5") {
  ggplot_barplot(PHYLOSEQ, args[3], "Subdivision", unknown, args[6], color_set)
}
