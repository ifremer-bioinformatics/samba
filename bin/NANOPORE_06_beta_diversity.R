#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Build ordination beta diversity analyses               ## 
##                                                                           ##
###############################################################################

## Get arguments from RScript command line ####
args = commandArgs(trailingOnly=TRUE)

## Load up the needed packages ####
requiredPackages = c("phyloseq", "metagMisc", "DESeq2", "vegan","ggplot2","stringr","metagenomeSeq","gridExtra","dplyr","grid", "pairwiseAdonis")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~ #
# Phyloseq input loading #
# ~~~~~~~~~~~~~~~~~~~~~~ #
PHYLOSEQ <- readRDS(args[1])

# ~~~~~~~~~~~~~~~~~~~~~~ #
# Normalisation          #
# ~~~~~~~~~~~~~~~~~~~~~~ #
if(args[2] == "rarefaction") {
    PHYLOSEQ_norm <- rarefy_even_depth(PHYLOSEQ, sample.size=min(sample_sums(PHYLOSEQ)), rngseed=1000, replace=FALSE, trimOTUs=TRUE)
    rarefied_table <- cbind(as.data.frame(otu_table(PHYLOSEQ_norm)), as.data.frame(tax_table(PHYLOSEQ_norm)))
    write.table(rarefied_table, "asv_table_for_stats_rarefied.tsv", sep="\t", col.names=NA, row.names=T, dec=".", quote=F)
} else if(args[2] == "css") {
    PHYLOSEQ_norm <- phyloseq_transform_css(PHYLOSEQ, norm=TRUE, log=TRUE)
    css_table <- cbind(as.data.frame(otu_table(PHYLOSEQ_norm)), as.data.frame(tax_table(PHYLOSEQ_norm)))
    write.table(css_table, "asv_table_for_stats_css_norm.tsv", sep="\t", col.names=NA, row.names=T, dec=".", quote=F)
} else if(args[2] == "deseq2") {
    deseq2 <- phyloseq_to_deseq2(PHYLOSEQ , ~ 1)
    gm_mean <- function(x, na.rm=TRUE) {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
              }
    geoMeans <- apply(counts(deseq2), 1, gm_mean)
    deseq2 <- estimateSizeFactors(deseq2, geoMeans=geoMeans)
    deseq2 <- estimateDispersions(deseq2, fitType="parametric", maxit=1000)
    deseq2_vst <- getVarianceStabilizedData(deseq2)
    deseq2_vst <- deseq2_vst + abs(min(deseq2_vst))
    PHYLOSEQ_norm <- PHYLOSEQ
    otu_table(PHYLOSEQ_norm) = otu_table(deseq2_vst, taxa_are_rows=TRUE)
    deseq2_table = cbind(as.data.frame(otu_table(PHYLOSEQ_norm)), as.data.frame(tax_table(PHYLOSEQ_norm)))
    write.table(deseq2_table, "asv_table_for_stats_deseq2_norm.tsv", sep="\t", col.names=NA, row.names=T, dec=".", quote=F)
}

# ~~~~~~~~~~~~~~~~~~~~~~ #
# Ordination calculation #
# ~~~~~~~~~~~~~~~~~~~~~~ #

index = c("jaccard", "bray")
ord_type = c("NMDS", "PCoA")

for(idx in index) { 
    for(type in ord_type) {
        assign(paste0("ordination_", type , "_", idx), ordinate(PHYLOSEQ_norm, type, idx, k=2, trymax=100))
        if (type == "NMDS" && eval(parse(text=paste0("ordination_", type, "_", idx)))$stress > 0.15) {
          assign(paste0("ordination_", type, "-", idx), ordinate(PHYLOSEQ_norm, type, idx, k=3, trymax=100))    
        }
        if (type == "NMDS" && eval(parse(text=paste0("ordination_", type, "_", idx)))$stress > 0.15) {
          assign(paste0("ordination_", type, "_", idx), ordinate(PHYLOSEQ_norm, type, idx, k=4, trymax=100))
        }
    }
}

# ~~~~~~~~~~~~~~~~~~~~~ #
# Variance explaination #
# ~~~~~~~~~~~~~~~~~~~~~ #

### using user-specified variable #### 
assign(args[3], get_variable(PHYLOSEQ_norm, args[3]))

for(idx in index) {
    assign(paste0("adonis_result_", idx), adonis2(phyloseq::distance(PHYLOSEQ_norm, idx) ~ eval(parse(text=args[3])), permutations = 999))
    adonis_result = eval(parse(text=paste0("adonis_result_", idx)))
    rownames(adonis_result)[1] <- args[3]
    assign(paste0("adonis_result_", idx), adonis_result)
    write.table(eval(parse(text=paste0("adonis_result_", idx))), paste0(args[3], "_permanova_result_", args[2], ".tsv"), col.names=NA, sep="\t", quote=F)
    assign(paste0("pairwiseadonis_result_", idx), pairwise.adonis(phyloseq::distance(PHYLOSEQ_norm, idx), eval(parse(text=args[3])), p.adjust.m = "bonferroni"))
    write.table(eval(parse(text=paste0("pairwiseadonis_result_", idx))), paste0(args[3], "_pairwiseAdonis_result_", args[2], ".tsv"), col.names=T, row.names=F, sep="\t", quote=F)
}

### using all available variables ###
all_var <- c()
for (var in sample_variables(PHYLOSEQ_norm) ) {
    l <- length(levels(as.factor(get_variable(PHYLOSEQ_norm, var))))
    if(l > 1 && l < nsamples(PHYLOSEQ_norm)){
        all_var <- cbind(all_var, var)
    }
}

variables <- paste(collapse =" + ", all_var )

justify <- function(x, hjust="center", vjust="center", draw=TRUE){
  w <- sum(x$widths)
  h <- sum(x$heights)
  xj <- switch(hjust,
               center = 0.4,
               left = 0.5*w,
               right=unit(1,"npc") - 0.5*w)
  yj <- switch(vjust,
               center = 0.5,
               bottom = 2*h,
               top=unit(1,"npc") - 0.5*h)
  x$vp <- viewport(x=xj, y=yj)
  if(draw) grid.draw(x)
  return(x)
}


for(idx in index) {
    f  <- paste("phyloseq::distance(PHYLOSEQ_norm, idx)"," ~ ", variables)
    permanova_result <- adonis2(as.formula(f), data=data.frame(sample_data(PHYLOSEQ_norm)), permutations=999, by="terms")
    assign(paste0("permanova_result_", idx), permanova_result)
    permanova_result_summary <- permanova_result[1:length(permanova_result$R2)-1,c(3,5)]
    permanova_result_summary$R2 <- round(permanova_result_summary$R2 * 100, 2)
    colnames(permanova_result_summary) <- c("R2", "pvalue")
    permanova_result_summary <- permanova_result_summary %>% mutate(Significance = case_when(pvalue <= 0.001 ~ "***", pvalue > 0.001 & pvalue < 0.01 ~ "**", pvalue > 0.01 & pvalue < 0.05 ~ "*", pvalue > 0.05 ~ "ns"))
    permanova_result_summary <- permanova_result_summary %>% replace(is.na(.), "")
    df_permanova <- tableGrob(permanova_result_summary,
                              theme=ttheme_minimal(core=list(bg_params=list(fill=adjustcolor("#FAE500",alpha.f=0.2)),
                                                             fg_params=list(fontface=3, fontsize=6)),
                                                   colhead=list(bg_params=list(fill="#00609B"),
                                                                fg_params=list(col="white", fontface=4L, fontsize=6)),
                                                   rowhead=list(fg_params=list(col="#00609B", fontface=3L, fontsize=6))))

    plot_df_permanova <- list(df_permanova)
    plot_df_permanova <- lapply(plot_df_permanova, justify, vjust="center", draw=FALSE)
    plot_h <- grid::convertHeight(sum(df_permanova$heights)*1.1, "in", TRUE)
    plot_w <- grid::convertWidth(sum(df_permanova$widths)*2.5, "in", TRUE) 

    ggsave(paste0("permanova_table_", idx, "_", args[2], ".png"), grid.arrange(grobs=plot_df_permanova, ncol=2), height=plot_h, width=plot_w, device="png")

}

# ~~~~~~~~~~~~~~~~~~~~~~~~ #
# Global NMDS construction #
# ~~~~~~~~~~~~~~~~~~~~~~~~ #

## function ####
plot.ordi <- function(physeq, ordination, ordi_type, idx, var, ordi_color, permanova_result, norm) {
  assign("R2_var", eval(parse(text=paste0("permanova_result_", idx)))[var, 3] * 100)
  if(ordi_type == "NMDS") {
      plot_ordination(physeq, ordination, type="samples", color=var) +
        theme_classic() +
        geom_point(size=5) +
        theme(legend.text=element_text(size=13)) +
        theme(legend.title=element_text(size=14)) +
        labs(color=var) +
        theme(axis.text=element_text(size=12,color="black")) +
        scale_fill_manual(values=alpha(ordi_color,0.4)) +
        scale_color_manual(values=ordi_color) +
        stat_ellipse(geom="polygon", alpha=0.1, type="t", aes_string(fill=var)) +
        theme(plot.caption = element_text(size=14)) +
        labs(caption = paste("Stress:",round(ordination$stress,4),
                             "\nPermutation test R2:",round(R2_var,2),
                             paste("\n", var, "variable significance: p-value"),permanova_result[var,]$`Pr(>F)`, sep=" "))
      ggsave(filename=paste0("ordination_", ordi_type, "_", idx, "_", var, "_", norm ,".svg",sep=""), device="svg", width=12, height=10)
      ggsave(filename=paste0("ordination_", ordi_type, "_", idx, "_", var, "_", norm, ".png",sep=""), device="png", width=12, height=10)
  } else {
      plot_ordination(physeq, ordination, type="samples", color=var) +
        theme_classic() +
        geom_point(size=5) +
        theme(legend.text=element_text(size=13)) +
        theme(legend.title=element_text(size=14)) +
        labs(color=var) +
        theme(axis.text=element_text(size=12,color="black")) +
        scale_fill_manual(values=alpha(ordi_color,0.4)) +
        scale_color_manual(values=ordi_color) +
        stat_ellipse(geom="polygon", alpha=0.1, type="t", aes_string(fill=var)) +
        theme(plot.caption = element_text(size=14)) +
        labs(caption = paste("Permutation test R2:",round(R2_var,2), 
                             paste("\n", var, "variable significance: p-value"),permanova_result[var,]$`Pr(>F)`, sep=" "))
      ggsave(filename=paste0("ordination_", ordi_type, "_", idx, "_", var, "_", norm, ".svg",sep=""), device="svg", width=12, height=10)
      ggsave(filename=paste0("ordination_", ordi_type, "_", idx, "_", var, "_", norm, ".png",sep=""), device="png", width=12, height=10)
  }
}

## plotting ####
ordi_color = adjustcolor(c("darkolivegreen3", "cornflowerblue", "darkorange2", "red", "deepskyblue4", "deeppink", "gray55", "khaki1", "mediumpurple4", "peachpuff2"), alpha.f=0.8)

for(dist in index) {
    plot.ordi(PHYLOSEQ_norm, eval(parse(text=paste0("ordination_NMDS_", dist))), "NMDS", dist, args[3], ordi_color, eval(parse(text=paste0("permanova_result_", dist))), args[2])
    plot.ordi(PHYLOSEQ_norm, eval(parse(text=paste0("ordination_PCoA_", dist))), "PCoA", dist, args[3], ordi_color, eval(parse(text=paste0("permanova_result_", dist))), args[2])
}
