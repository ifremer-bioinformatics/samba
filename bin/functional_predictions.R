#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Statistical representation of functional predictions   ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("dplyr","ggplot2","RColorBrewer","svglite","vegan","stringr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

functional_predictions <- function(pred,metadata,criteria,pred_plot,name,microDecon,control) {

  # Format functional predictions tables ####
  pred = read.table(pred,h=T,sep="\t",check.names=F,quote="")
  pred = pred[,c(1,3:length(pred))]
  rownames(pred) = pred[,1]
  pred = pred %>% select (-1)
  pred = pred[,colSums(pred) > 0]
  t_pred = as.data.frame(t(pred))
  t_pred = round(t_pred)

  # Load metatada
  metadata = read.table(metadata, h=T, sep="\t", check.names=F)
  metadata = metadata[metadata[,1] %in% rownames(t_pred),]
  if (microDecon == "true") {
     control_list =  unlist(strsplit(control,","))
     metadata = metadata[!metadata[,1] %in% control_list, ]
  }

  # Perform normalization
  min_depth = min(rowSums(t_pred))
  t_pred_norm = as.data.frame(round(rrarefy(t_pred, min_depth)))

  # Calculate the distance matrix
  pred_dist = as.matrix((vegdist(t_pred_norm, "bray")))

  # Perform NMDS
  pred_nmds = metaMDS(pred_dist, k=2, trymax=100)
  pred_nmds_stress = round(pred_nmds$stress,3)

  # Build a data frame with NMDS coordinates and metadata
  pred_nmds1 = pred_nmds$points[,1]
  pred_nmds2 = pred_nmds$points[,2]
  pred_nmds_data = data.frame(MDS1 = pred_nmds1, MDS2 = pred_nmds2, metadata = metadata[,criteria], samples = rownames(t_pred))

  # ADONIS statistic
  pred_adonis = adonis(pred_dist ~ metadata[,criteria])

  # Plot
  ggplot(pred_nmds_data, aes(x=MDS1, y=MDS2, col=metadata)) +
    theme_classic() +
    geom_point(shape=19, size=3) +
    geom_text(data=pred_nmds_data,aes(x=MDS1,y=MDS2,label=samples),size=3,vjust=2) +
    stat_ellipse(geom="polygon",alpha=0.1,type="t",aes(fill=metadata),lty=2) +
    scale_color_brewer(palette="Set1") +
    theme(legend.title = element_blank()) +
    theme(axis.text = element_text(colour = "black", size = 10)) +
    labs(caption = paste("Stress:",pred_nmds_stress,
                         "\nAdonis statistic R:",round(pred_adonis$aov.tab$R2[1]*100,2),
                         paste("\nAdonis based on ", "transect_name",": p-value"),pred_adonis$aov.tab$`Pr(>F)`[1],sep=" "))
  ggsave(filename=paste(name,pred_plot,".svg",sep=""), device="svg", width = 12, height = 10)
  ggsave(filename=paste(name,pred_plot,".png",sep=""), device="png", width = 12, height = 10)
}

main_ec <- function(){
  # Get arguments from RScript command line
  args = commandArgs(trailingOnly=TRUE)
  pred_ec = args[1]
  metadata =   args[4]
  criteria = args[5]
  pred_plot = args[6]
  name = "EC_"
  microDecon = args[7]
  control = args[8]
  functional_predictions(pred_ec,metadata,criteria,pred_plot,name,microDecon,control)
}

if (!interactive()) {
  main_ec()
}

main_ko <- function(){
  # Get arguments from RScript command line
  args = commandArgs(trailingOnly=TRUE)
  pred_ko = args[2]
  metadata = args[4]
  criteria = args[5]
  pred_plot = args[6]
  name = "KO_"
  microDecon = args[7]
  control = args[8]
  functional_predictions(pred_ko,metadata,criteria,pred_plot,name,microDecon,control)
}

if (!interactive()) {
  main_ko()
}

main_metacyc <- function(){
  # Get arguments from RScript command line
  args = commandArgs(trailingOnly=TRUE)
  pred_path = args[3]
  metadata = args[4]
  criteria = args[5]
  pred_plot = args[6]
  name = "MetaCyc_"
  microDecon = args[7]
  control = args[8]
  functional_predictions(pred_path,metadata,criteria,pred_plot,name,microDecon,control)
}

if (!interactive()) {
  main_metacyc()
}
