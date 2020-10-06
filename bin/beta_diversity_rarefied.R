#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Automated Statistical Analyses of Metabarcoding Data   ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("dplyr","stringr","ggplot2","RColorBrewer","svglite","tidyr","gridExtra","egg","vegan","dendextend","BiocManager","phyloseq","plotly","htmlwidgets")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

# @@@@@@@@@@@@@@@@@@@ #
#                     #
# Rarefaction process #
#                     #
# @@@@@@@@@@@@@@@@@@@ #

args = commandArgs(trailingOnly=TRUE)
PHYLOSEQ = readRDS(args[1])
final_rarefied_ASV_table_with_taxonomy = args[2]

PHYLOSEQ_rarefied = rarefy_even_depth(PHYLOSEQ, sample.size=min(sample_sums(PHYLOSEQ)),rngseed=1000,replace=TRUE,trimOTUs=TRUE)
rarefied_table = cbind(as.data.frame(otu_table(PHYLOSEQ_rarefied)),as.data.frame(tax_table(PHYLOSEQ_rarefied)))
write.table(rarefied_table,final_rarefied_ASV_table_with_taxonomy,sep="\t",col.names=T,row.names=T,dec=".",quote=F)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#                                                                               #
# Function to standardize ordination analysis and plot for each distance matrix #
#                                                                               #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

betadiversity_rarefied <- function (PHYLOSEQ_rarefied, distance, metadata, variance_significance_tests_rarefied, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc, plot_pie, plotly_js) {
   
    #~~~~~~~~~~~~~~~#
    # Rarefied data #
    #~~~~~~~~~~~~~~~#

    ## Ordination process ####
    metadata = read.table(metadata, row.names=1, h=T, sep="\t", check.names=FALSE)
    
    if (distance == "jaccard" | distance == "bray") {
      ord_rarefied_nmds = ordinate(PHYLOSEQ_rarefied, "NMDS", distance, trymax = 1000)
      ord_rarefied_pcoa = ordinate(PHYLOSEQ_rarefied, "PCoA", distance, trymax = 1000)
    }
    else {
      ord_rarefied_nmds = ordinate(PHYLOSEQ_rarefied, "NMDS", distance)
      ord_rarefied_pcoa = ordinate(PHYLOSEQ_rarefied, "PCoA", distance)
    }

    color_vector = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))
    color_samples = color_vector[1:length(levels(metadata[,criteria]))]
    
    ## Statistic about variance explained ####
    ### using user-specified variable #### 
    group_rarefied = get_variable(PHYLOSEQ_rarefied, criteria)
    adonis_result_rarefied = adonis(distance(PHYLOSEQ_rarefied,distance) ~ group_rarefied, permutations = 9999)

    ### using all available variables ####
    all_var = c()
    for (var in sample_variables(PHYLOSEQ_rarefied) ) {
      l = length(levels(as.factor(get_variable(PHYLOSEQ_rarefied,var))))
      if(l > 1 && l < nsamples(PHYLOSEQ_rarefied)){
        all_var <- cbind(all_var,var)
      }
    }

    variables = paste(collapse =" + ", all_var )
    
    sink(file = paste(variance_significance_tests_rarefied,distance,".txt",sep="") , type = "output")
      f  = paste("distance(PHYLOSEQ_rarefied,distance)"," ~ ", variables)
      cat(sep = "", "###############################################################\n",
                "#Perform Adonis test on multiple variables: ",variables," using the ",distance," distance matrix")
      adonis_all_rarefied=adonis(as.formula(f), data=metadata, perm = 9999)
      print(adonis_all_rarefied)
      cat("\n\n")
    sink()
    
    ### Explained variance graphs ####
    ExpVar_perc = adonis_all_rarefied$aov.tab$R2[-length(adonis_all_rarefied$aov.tab$R2)]*100
    ExpVar_perc[ExpVar_perc < 0] <- 0
    ExpVar_name = rownames(adonis_all_rarefied$aov.tab)[-length(rownames(adonis_all_rarefied$aov.tab))]
    ExpVar_piedata = data.frame(ExpVar_name,ExpVar_perc)
    ExpVar_piedata = ExpVar_piedata[order(ExpVar_piedata$ExpVar_perc),]
    ExpVar_pielabels = sprintf("%s = %3.1f%s", ExpVar_piedata$ExpVar_name,ExpVar_piedata$ExpVar_perc, "%")

    plot.pie(ExpVar_piedata$ExpVar_perc, ExpVar_pielabels, distance, plot_pie, 12, 10)
    plot.pie.interactive(ExpVar_piedata,ExpVar_name,ExpVar_perc,plot_pie,distance, plotly_js)

    ## Ordination plots ####
    ### PHYLOSEQ_OBJ, Ordination, variable to test, colors to use, adonis result, ordination plot name, distance, width of graph, heigth of graph, graph title
    plot.nmds(PHYLOSEQ_rarefied, ord_rarefied_nmds, criteria, color_samples, adonis_result_rarefied, nmds_rarefied, distance, 12, 10, paste("NMDS on rarefied data","based on",distance,"distance",sep=" "))
    plot.nmds.interactive(ord_rarefied_nmds, metadata, criteria, color_samples, adonis_result_rarefied, nmds_rarefied, distance, paste("NMDS on rarefied data","based on",distance,"distance",sep=" "), plotly_js)
    plot.pcoa(PHYLOSEQ_rarefied, ord_rarefied_pcoa, criteria, color_samples, adonis_result_rarefied, pcoa_rarefied, distance, 12, 10, paste("MDS-PCoA on rarefied data","based on",distance,"distance",sep=" "))
    plot.pcoa.interactive(ord_rarefied_pcoa, metadata, criteria, color_samples, adonis_result_rarefied, pcoa_rarefied, distance, paste("MDS-PCoA on rarefied data","based on",distance,"distance",sep=" "), plotly_js)

    ## Hierarchical clustering ####
    dist = distance(PHYLOSEQ_rarefied, distance, type="samples")
    hc = hclust(dist, method = method_hc)
    dendro = as.dendrogram(hc)
    group = data.frame(PHYLOSEQ_rarefied@sam_data[,criteria])[,1]
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
    PHYLOSEQ_rarefied = PHYLOSEQ_rarefied 
    distance = "jaccard"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_rarefied = args[6]
    pcoa_rarefied = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    variance_significance_tests_rarefied=args[10]
    plot_pie = args[11]
    plotly_js = args[12]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/bin/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_rarefied(PHYLOSEQ_rarefied, distance, metadata, variance_significance_tests_rarefied, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc, plot_pie, plotly_js)
}

if (!interactive()) {
        main_jaccard()
}

main_bray <- function(){
    PHYLOSEQ_rarefied = PHYLOSEQ_rarefied
    distance = "bray"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_rarefied = args[6]
    pcoa_rarefied = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    variance_significance_tests_rarefied=args[10]
    plot_pie = args[11]
    plotly_js = args[12]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/bin/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_rarefied(PHYLOSEQ_rarefied, distance, metadata, variance_significance_tests_rarefied, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc, plot_pie, plotly_js)
}

if (!interactive()) {
        main_bray()
}

main_unifrac <- function(){
    PHYLOSEQ_rarefied = PHYLOSEQ_rarefied
    distance = "unifrac"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_rarefied = args[6]
    pcoa_rarefied = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    variance_significance_tests_rarefied=args[10]
    plot_pie = args[11]
    plotly_js = args[12]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/bin/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_rarefied(PHYLOSEQ_rarefied, distance, metadata, variance_significance_tests_rarefied, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc, plot_pie, plotly_js)
}

if (!interactive()) {
        main_unifrac()
}

main_wunifrac <- function(){
    PHYLOSEQ_rarefied = PHYLOSEQ_rarefied
    distance = "wunifrac"
    # Get criteria and replace "-" character by "_"
    criteria = str_replace(args[3], "-", "_")
    metadata = args[4]
    workflow_dir = args[5]
    nmds_rarefied = args[6]
    pcoa_rarefied = args[7]
    method_hc = args[8]
    plot_hc = args[9]
    variance_significance_tests_rarefied=args[10]
    plot_pie = args[11]
    plotly_js = args[12]
    # Check if functions are loaded, if not source them
    if (!exists("plot.nmds", mode="function")) source(gsub(" ", "", paste(workflow_dir,"/bin/beta_diversity_graphs.R")))
    # Beta diversity analyses
    betadiversity_rarefied(PHYLOSEQ_rarefied, distance, metadata, variance_significance_tests_rarefied, criteria, nmds_rarefied, pcoa_rarefied, method_hc, plot_hc, plot_pie, plotly_js)
}

if (!interactive()) {
        main_wunifrac()
}