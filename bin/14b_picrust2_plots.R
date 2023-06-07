#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Statistical representation of functional predictions   ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("dplyr","vegan","ggplot2","svglite","stringr","ggord","compositions")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

functional_predictions <- function(pred,table,metadata,pred_plot,name,var) {

  # Format functional predictions tables ####
  rownames(pred) = pred[,1]
  pred = pred %>% select(-1)
  pred_clr = clr(pred)
  t_pred_clr = as.data.frame(t(pred_clr))

  # Format ASV table
  asv_table = read.table(table, h=T, sep="\t", row.names=1, check.names=F)
  asv_table = asv_table[-length(asv_table)]
  t_asv_table = t(asv_table)

  # Load metatada
  metadata = read.table(metadata, h=T, sep="\t", check.names=F)
  metadata_filtered = metadata[metadata[,1] %in% rownames(t_pred_clr),]

  # RDA analysis
  rda_intercept = rda(t_asv_table ~ 1, data=t_pred_clr)
  rda_all_var_exp = rda(t_asv_table ~ ., data=t_pred_clr)

  # Forward selection of significant predictions
  fwd.sel = ordiR2step(rda_intercept, scope=rda_all_var_exp, R2scope=FALSE, perm.max=1000, trace=FALSE)
  f = str_split(as.character(fwd.sel$call), "~", 1)[[2]]
  rda_signif = rda(as.formula(f), data=t_pred_clr)
  R2 = RsquareAdj(rda_signif)

  # Test significance of model and terms
  anova_model = anova.cca(rda_signif, perm.max=500)
  anova_terms = anova.cca(rda_signif, permu=500, by="term")
  sink(file=paste("significance_",name,".txt",sep="") , type = "output")
  print(anova_terms)
  sink()

  # Plot
  cols = adjustcolor(c("darkolivegreen3", "cornflowerblue", "darkorange2", "red", "deepskyblue4", "deeppink", "gray55", "khaki1", "mediumpurple4", "peachpuff2"), alpha.f=0.8)
  cols_var = cols[1:length(unique(metadata_filtered[, var]))]

  plot = ggord(rda_signif, metadata_filtered[, var], cols=cols_var, txt=4, arrow=0.2, addsize=-5, size=3, hull=FALSE, alpha_el=0.2, repel=F, ext=1, vec_ext=1) +
    theme_light() +
    theme(legend.title=element_blank()) +
    labs(caption=paste("Model significance pvalue: ", anova_model$`Pr(>F)`[1],"\n",
                       "Adjusted R2: ", round(R2$adj.r.squared*100,2),"%", sep=""))
  ggsave(filename=paste(name,pred_plot,"_",var,".svg",sep=""), device="svg", width = 12, height = 10)
  ggsave(filename=paste(name,pred_plot,"_",var,".png",sep=""), device="png", width = 12, height = 10)

}

args = commandArgs(trailingOnly=TRUE)
traits = unlist(strsplit(args[1], ","))

if( "EC" %in% traits) {
    main_ec <- function(){
      pred_ec = read.table("14_PICRUSt2_predictions_output/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz", h=T, sep="\t", check.names=F, quote="")
      pred_ec = data.frame(pred_ec[,2:length(pred_ec)] %>% group_by(description) %>% summarise_all(.funs = sum), check.names=F)
      table = args[2]
      metadata = args[3]
      pred_plot = args[4]
      var = args[5]
      name = "EC"
      functional_predictions(pred_ec,table,metadata,pred_plot,name, var)
    }
    
    if (!interactive()) {
      main_ec()
    }
}

if( "KO" %in% traits) {
    main_ko <- function(){
      pred_ko = read.table("14_PICRUSt2_predictions_output/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz", h=T, sep="\t", check.names=F, quote="")
      pred_ko = data.frame(pred_ko[,2:length(pred_ko)] %>% group_by(description) %>% summarise_all(.funs = sum), check.names=F)
      table = args[2]
      metadata = args[3]
      pred_plot = args[4]
      var = args[5]
      name = "KO"
      functional_predictions(pred_ko,table,metadata,pred_plot,name,var)
    }

    if (!interactive()) {
      main_ko()
    }
}


if( "COG" %in% traits) {
    main_cog <- function(){
      pred_cog = read.table("14_PICRUSt2_predictions_output/COG_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz", h=T, sep="\t", check.names=F, quote="")
      pred_cog = data.frame(pred_cog[,2:length(pred_cog)] %>% group_by(description) %>% summarise_all(.funs = sum), check.names=F)
      table = args[2]
      metadata = args[3]
      pred_plot = args[4]
      var = args[5]
      name = "COG"
      functional_predictions(pred_cog,table,metadata,pred_plot,name,var)
    }

    if (!interactive()) {
      main_cog()
    }
}

if( "PFAM" %in% traits) {
    main_pfam <- function(){
      pred_pfam = read.table("14_PICRUSt2_predictions_output/PFAM_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz", h=T, sep="\t", check.names=F, quote="")
      pred_pfam = data.frame(pred_pfam[,2:length(pred_pfam)] %>% group_by(description) %>% summarise_all(.funs = sum), check.names=F)
      table = args[2]
      metadata = args[3]
      pred_plot = args[4]
      var = args[5]
      name = "PFAM"
      functional_predictions(pred_pfam,table,metadata,pred_plot,name,var)
    }

    if (!interactive()) {
      main_pfam()
    }
}

if( "TIGRFAM" %in% traits) {
    main_tigrfam <- function(){
      pred_tigrfam = read.table("14_PICRUSt2_predictions_output/TIGRFAM_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz", h=T, sep="\t", check.names=F, quote="")
      pred_tigrfam = data.frame(pred_tigrfam[,2:length(pred_tigrfam)] %>% group_by(description) %>% summarise_all(.funs = sum), check.names=F)
      table = args[2]
      metadata = args[3]
      pred_plot = args[4]
      var = args[5]
      name = "TIGRFAM"
      functional_predictions(pred_tigrfam,table,metadata,pred_plot,name,var)
    }

    if (!interactive()) {
      main_tigrfam()
    }
}

main_metacyc <- function(){
  pred_path = read.table("14_PICRUSt2_predictions_output/pathways_out/path_abun_unstrat_descrip.tsv.gz", h=T, sep="\t", check.names=F, quote="")
  pred_path = data.frame(pred_path[,2:length(pred_path)] %>% group_by(description) %>% summarise_all(.funs = sum), check.names=F)
  table = args[2]
  metadata = args[3]
  pred_plot = args[4]
  var = args[5]
  name = "MetaCyc"
  functional_predictions(pred_path,table,metadata,pred_plot,name,var)
}

if (!interactive()) {
  main_metacyc()
}

main_metacyc_sec <- function(){
  pred_path_sec = read.table("14_PICRUSt2_predictions_output/pathways_out/path_abun_unstrat_descrip_secondary_level.tsv.gz", h=T, sep="\t", check.names=F, quote="")
  pred_path_sec = data.frame(pred_path_sec[,2:length(pred_path_sec)] %>% group_by(description) %>% summarise_all(.funs = sum), check.names=F)
  table = args[2]
  metadata = args[3]
  pred_plot = args[4]
  var = args[5]
  name = "MetaCyc_secondary_level"
  functional_predictions(pred_path_sec,table,metadata,pred_plot,name,var)
}

if (!interactive()) {
  main_metacyc_sec()
}
