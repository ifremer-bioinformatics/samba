#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Automated Statistical Analyses of Metabarcoding Data   ##
##                                                                           ##
###############################################################################

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

otu_table(PHYLOSEQ) = otu_table(PHYLOSEQ)[,colSums(otu_table(PHYLOSEQ) !=0) > 1]
otu_table(PHYLOSEQ) = otu_table(PHYLOSEQ)[rowSums(otu_table(PHYLOSEQ)) > 0,]
css = phyloseq_to_metagenomeSeq(PHYLOSEQ)
p = cumNormStatFast(css)
css = cumNorm(css, p=p)
css_norm_factor = normFactors(css)
CSS_TABLE = MRcounts(css, norm = T)
PHYLOSEQ_css = PHYLOSEQ
otu_table(PHYLOSEQ_css) = otu_table(CSS_TABLE,taxa_are_rows=TRUE)
CSS_normalized_table = cbind(as.data.frame(otu_table(PHYLOSEQ_css)),as.data.frame(tax_table(PHYLOSEQ_css)))
write.table(CSS_normalized_table,final_css_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=".")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#                                                                               #
# Function to standardize ordination analysis and plot for each distance matrix #
#                                                                               #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

betadiversity_css <- function (PHYLOSEQ_css, distance, metadata, variance_significance_tests_css, criteria, nmds_css, pcoa_css, method_hc, plot_hc, plot_pie) {

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
    
    ## Statistic about variance explained ####
    ### using user-specified variable #### 
    group_css = get_variable(PHYLOSEQ_css, criteria)
    adonis_result_css = adonis(distance(PHYLOSEQ_css,distance) ~ group_css, permutations = 9999)
    
    ### using all available variables ####
    all_var = c()
    for (var in sample_variables(PHYLOSEQ_css) ) {
      l = length(levels(as.factor(get_variable(PHYLOSEQ_css,var))))
      if(l > 1 && l < nsamples(PHYLOSEQ_css)){
        all_var <- cbind(all_var,var)
      }
    }

    variables = paste(collapse =" + ", all_var )
    
    sink(file = paste(variance_significance_tests_css,distance,".txt",sep="") , type = "output")
      f  = paste("distance(PHYLOSEQ_css,distance)"," ~ ", variables)
      cat(sep = "", "###############################################################\n",
                "#Perform Adonis test on multiple variables: ",variables," using the ",distance," distance matrix")
      adonis_all_css=adonis(as.formula(f), data=metadata, perm = 9999)
      print(adonis_all_css)
      cat("\n\n")
    sink()

    ### Explained variance graphs ####
    ExpVar_perc = adonis_all_css$aov.tab$R2[-length(adonis_all_css$aov.tab$R2)]*100
    ExpVar_name = rownames(adonis_all_css$aov.tab)[-length(rownames(adonis_all_css$aov.tab))]
    ExpVar_piedata = data.frame(ExpVar_name,ExpVar_perc)
    ExpVar_piedata = ExpVar_piedata[order(ExpVar_piedata$ExpVar_perc),]
    ExpVar_pielabels = sprintf("%s = %3.1f%s", ExpVar_piedata$ExpVar_name,ExpVar_piedata$ExpVar_perc, "%")

    plot.pie(ExpVar_piedata$ExpVar_perc, ExpVar_pielabels, distance, plot_pie, 12, 10)

    ## Ordination plots ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, adonis result, ordination plot name, distance, width of graph, heigth of graph, graph title
    plot.nmds(PHYLOSEQ_css, ord_css_nmds, criteria, color_samples, adonis_result_css, nmds_css, distance, 12, 10, paste("NMDS on CSS normalized data","based on",distance,"distance",sep=" "))
    plot.pcoa(PHYLOSEQ_css, ord_css_pcoa, criteria, color_samples, adonis_result_css, pcoa_css, distance, 12, 10, paste("MDS-PCoA on CSS normalized data","based on",distance,"distance",sep=" "))

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
    variance_significance_tests_css=args[10]
    plot_pie = args[11]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/bin/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_css(PHYLOSEQ_css, distance, metadata, variance_significance_tests_css, criteria, nmds_css, pcoa_css, method_hc, plot_hc, plot_pie)
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
    variance_significance_tests_css=args[10]
    plot_pie = args[11]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/bin/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_css(PHYLOSEQ_css, distance, metadata, variance_significance_tests_css, criteria, nmds_css, pcoa_css, method_hc, plot_hc, plot_pie)
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
    variance_significance_tests_css=args[10]
    plot_pie = args[11]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/bin/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_css(PHYLOSEQ_css, distance, metadata, variance_significance_tests_css, criteria, nmds_css, pcoa_css, method_hc, plot_hc, plot_pie)
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
    variance_significance_tests_css=args[10]
    plot_pie = args[11]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/bin/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_css(PHYLOSEQ_css, distance, metadata, variance_significance_tests_css, criteria, nmds_css, pcoa_css, method_hc, plot_hc, plot_pie)
}

if (!interactive()) {
        main_wunifrac()
}
