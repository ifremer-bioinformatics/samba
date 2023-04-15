#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Alpha diversity analysis from NANOPORE data            ##
##                                                                           ##
###############################################################################

## Get arguments from RScript command line ####
args = commandArgs(trailingOnly=TRUE)

## Load up the needed packages ####
requiredPackages = c("RColorBrewer", "phyloseq", "dplyr", "tidyr", "ggplot2","ggpubr", "stringr", "rstatix","grid", "vegan")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}
color_set = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))

# ~~~~~~~~~~~~~~~ #
# Alpha diversity #
# ~~~~~~~~~~~~~~~ #

## Phyloseq object import ####
PHYLOSEQ = readRDS(args[1])

## Aggregate sequence at each taxa level ####
PHYLOSEQ_PHYLUM = tax_glom(PHYLOSEQ, "Phylum")
PHYLOSEQ_CLASS = tax_glom(PHYLOSEQ, "Class")
PHYLOSEQ_ORDER = tax_glom(PHYLOSEQ, "Order")
PHYLOSEQ_FAMILY = tax_glom(PHYLOSEQ, "Family")
PHYLOSEQ_GENUS = tax_glom(PHYLOSEQ, "Genus")
PHYLOSEQ_SPECIES = tax_glom(PHYLOSEQ, "Species")

## Calcul of diversity indexes (Chao1, Shannon, InvSimpson) ####
alpha_rich = estimate_richness(PHYLOSEQ_SPECIES, measures=c("Observed", "Chao1", "ACE", "Shannon", "InvSimpson"))
evenness = microbiome::evenness(PHYLOSEQ_SPECIES, "pielou")
df_alpha_div = cbind(alpha_rich, evenness)

write.table(df_alpha_div, args[2], col.names=NA, row.names=TRUE, sep="\t", dec=",", quote=FALSE)

df_alpha_div_sdf = data.frame(df_alpha_div, sample_data(PHYLOSEQ_SPECIES))
df_alpha_div_sdf = gather(df_alpha_div_sdf, key="Measure", value="Value", Observed, Chao1, ACE, Shannon, InvSimpson, pielou)
df_alpha_div_sdf$Measure = factor(df_alpha_div_sdf$Measure, levels=c("Observed", "Chao1", "ACE", "Shannon", "InvSimpson", "pielou"))

## Alpha diversity plots according to user-selected variable ####
test_formula = paste("Value", args[3], sep="~")
var_value = length(unique(df_alpha_div_sdf[, args[3]]))
stats_var = df_alpha_div_sdf %>% group_by(Measure) %>% t_test(as.formula(test_formula)) %>% adjust_pvalue(method="bonferroni") %>% add_significance() %>% add_xy_position(x = "supp")
var_color = adjustcolor(color_set[2:var_value+1)], alpha.f=0.8)

alpha_bxp_var = ggboxplot(df_alpha_div_sdf, x=args[3], y="Value", fill=args[3], facet.by="Measure", add="jitter", scales="free_y", add.params=list(size=1)) +
  stat_pvalue_manual(stats_var, hide.ns=TRUE, label="{signif(p.adj,1)}{p.adj.signif}", label.size=4) +
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
asv_matrix = as.matrix(t(data.frame(otu_table(PHYLOSEQ_SPECIES), check.names=F)))
metadata_df = data.frame(sample_data(PHYLOSEQ_SPECIES))
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
