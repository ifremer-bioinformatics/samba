
#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Script name: R script for the SAMBA-nextflow workflow                   ####
##                                                                           ##
## Purpose of script: Automated Statistical Analyses of Metabarcoding Data   ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Cyril NOEL and Laure QUINTRIC                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2019-12-13                                                 ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## Emails: cyril.noel@ifremer.fr and laure.quintric@ifremer.fr             ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Notes: This part of the script performs the beta diversity (NMDS, PCoa &  ##
##        Hierachical Clustering) on CSS normalized ASV table based on       ##
##        four distance matrices (Jaccard, Bray-Curtis, UniFrac & Weighted   ##
##        UniFrac) 							     ##   
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Load up the needed packages ####
requiredPackages = c("dplyr","stringr","ggplot2","RColorBrewer","svglite","tidyr","gridExtra","egg","vegan","dendextend","BiocManager","phyloseq","metagenomeSeq")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@ #
#                           #
# CSS normalization process #
#                           #
# @@@@@@@@@@@@@@@@@@@@@@@@@ #

args = commandArgs(trailingOnly=TRUE)
PHYLOSEQ = readRDS(args[1])
final_css_ASV_table_with_taxonomy = args[2]

css = phyloseq_to_metagenomeSeq(PHYLOSEQ)
p = cumNormStatFast(css)
css = cumNorm(css, p=p)
css_norm_factor = normFactors(css)
CSS_TABLE = MRcounts(css, norm = T)
PHYLOSEQ_css = PHYLOSEQ
otu_table(PHYLOSEQ_css) = otu_table(CSS_TABLE,taxa_are_rows=TRUE)
CSS_normalized_table = cbind(as.data.frame(otu_table(PHYLOSEQ_css)),as.data.frame(tax_table(PHYLOSEQ_css)))
write.table(CSS_normalized_table,final_css_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=",")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#										#
# Function to standardize ordination analysis and plot for each distance matrix #
#										#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

betadiversity_css <- function (PHYLOSEQ_css, distance, metadata, criteria, nmds_css, pcoa_css, method_hc, plot_hc) {

    #~~~~~~~~~~~~~~~~~~~#
    # CSS normalization #
    #~~~~~~~~~~~~~~~~~~~#
    
    ## Ordination process ####
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    
    if (distance == "jaccard" | distance == "bray") {    
      ord_css_nmds = ordinate(PHYLOSEQ_css,"NMDS", distance, trymax = 1000)
      ord_css_pcoa = ordinate(PHYLOSEQ_css,"PCoA", distance, trymax = 1000)
    }
    else {
      ord_css_nmds = ordinate(PHYLOSEQ_css,"NMDS", distance)
      ord_css_pcoa = ordinate(PHYLOSEQ_css,"PCoA", distance)
    }

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]
    
    group_css = get_variable(PHYLOSEQ_css, criteria)
    anosim_result_css = anosim(distance(PHYLOSEQ_css,distance),group_css, permutations = 999)
    
    ## Ordination plots ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, anosim result, ordination plot name, distance, width of graph, heigth of graph, graph title
    plot.nmds(PHYLOSEQ_css, ord_css_nmds, criteria, color_samples, anosim_result_css, nmds_css, distance, 12, 10, paste("NMDS on CSS normalized data","based on",distance,"distance",sep=" "))
    plot.pcoa(PHYLOSEQ_css, ord_css_pcoa, criteria, color_samples, anosim_result_css, pcoa_css, distance, 12, 10, paste("MDS-PCoA on CSS normalized data","based on",distance,"distance",sep=" "))

    ## Hierarchical clustering ####    
    dist = distance(PHYLOSEQ_css, distance, type="samples")
    hc = hclust(dist, method = method_hc)
    dendro = as.dendrogram(hc)
    group = data.frame(PHYLOSEQ_css@sam_data[,criteria])[,1]
    n_group = length(unique(group))
    cols = color_vector[1:n_group]
    col_group = cols[group]
    plot.hc(dendro, group, cols, col_group, method_hc, plot_hc, distance, 12, 10)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#                                             #
# Ordination process for each distance matrix #
#                                             #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

main_jaccard <- function(){
    PHYLOSEQ_css = PHYLOSEQ_css
    distance = "jaccard"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_css = args[6]
    pcoa_css = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_css(PHYLOSEQ_css, distance, metadata, criteria, nmds_css, pcoa_css, method_hc, plot_hc)
}

if (!interactive()) {
        main_jaccard()
}

main_bray <- function(){
    PHYLOSEQ_css = PHYLOSEQ_css
    distance = "bray"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_css = args[6]
    pcoa_css = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_css(PHYLOSEQ_css, distance, metadata, criteria, nmds_css, pcoa_css, method_hc, plot_hc)
}

if (!interactive()) {
        main_bray()
}

main_unifrac <- function(){
    PHYLOSEQ_css = PHYLOSEQ_css
    distance = "unifrac"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_css = args[6]
    pcoa_css = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_css(PHYLOSEQ_css, distance, metadata, criteria, nmds_css, pcoa_css, method_hc, plot_hc)
}

if (!interactive()) {
        main_unifrac()
}

main_wunifrac <- function(){
    PHYLOSEQ_css = PHYLOSEQ_css
    distance = "wunifrac"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_css = args[6]
    pcoa_css = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/lib/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_css(PHYLOSEQ_css, distance, metadata, criteria, nmds_css, pcoa_css, method_hc, plot_hc)
}

if (!interactive()) {
        main_wunifrac()
}
