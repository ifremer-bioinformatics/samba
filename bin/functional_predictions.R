#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Statistical representation of functional predictions   ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("dplyr","ggplot2","RColorBrewer","svglite","vegan","stringr","plotly","htmlwidgets")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

functional_predictions <- function(pred,metadata,criteria,pred_plot,name,microDecon,control, plotly_js) {

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
  pred_nmds_data = data.frame(MDS1 = pred_nmds1, MDS2 = pred_nmds2, Condition = metadata[,criteria], SampleID = rownames(t_pred))

  # ADONIS statistic
  pred_adonis = adonis(pred_dist ~ metadata[,criteria])

  # Plot
  ggplot(pred_nmds_data, aes(x=MDS1, y=MDS2, col=Condition)) +
    theme_classic() +
    geom_point(shape=19, size=3) +
    geom_text(data=pred_nmds_data,aes(x=MDS1,y=MDS2,label=SampleID),size=3,vjust=2) +
    stat_ellipse(geom="polygon",alpha=0.1,type="t",aes(fill=Condition),lty=2) +
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") +
    theme(legend.title = element_blank()) +
    theme(axis.text = element_text(colour = "black", size = 10)) +
    labs(caption = paste("Stress:",pred_nmds_stress,
                         "\nAdonis statistic R:",round(pred_adonis$aov.tab$R2[1]*100,2),
                         paste("\nAdonis based on ", "transect_name",": p-value"),pred_adonis$aov.tab$`Pr(>F)`[1],sep=" "))
  ggsave(filename=paste(name,pred_plot,".svg",sep=""), device="svg", width = 12, height = 10)
  ggsave(filename=paste(name,pred_plot,".png",sep=""), device="png", width = 12, height = 10)

  plotly_mod_dep = function(p,js_file){
    deps <- p$dependencies
    deps_urls <- purrr::map(
      deps,
      ~if(.x$name == "plotly-basic") {
        .x$src = list(file=getwd())
        .x$script = js_file
        .x
      } else {
        .x
      }
    )
    p$dependencies <- deps_urls
    p
  }

  plot_pred_nmds = ggplot(pred_nmds_data, aes(x=MDS1, y=MDS2)) +
    theme_classic() +
    stat_ellipse(geom="polygon",alpha=0.2,type="t",aes(fill=Condition),lty=1) +
    geom_point(shape=19, size=2, alpha=0.8,aes(color=Condition, label = SampleID)) +
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") +
    theme(legend.title = element_blank()) +
    theme(plot.background = element_rect(fill="#fafafa")) +
    theme(axis.text = element_text(colour = "black", size = 10)) +
    labs(caption = paste("Stress:",pred_nmds_stress,
                       "\nAdonis statistic R:",round(pred_adonis$aov.tab$R2[1]*100,2),
                       paste("\nAdonis based on ", "transect_name",": p-value"),pred_adonis$aov.tab$`Pr(>F)`[1],sep=" "))

  plotly_nmds = ggplotly(plot_pred_nmds) %>% partial_bundle(local=FALSE) %>% plotly_mod_dep(js_file=plotly_js) %>% layout(autosize=F,margin=list(r=0,l=0,t=0,b=0,pad=0),legend=list(bgcolor="#fafafa"),paper_bgcolor="#fafafa")
  for (i in c(1:length(plotly_nmds$x$data))) {
    tmp_replace_name = str_remove_all(plotly_nmds$x$data[[i]]$name,"\\(")
    tmp_replace_name = str_remove_all(tmp_replace_name, ",1\\)")
    plotly_nmds$x$data[[i]]$name = tmp_replace_name
    tmp_replace_legendgroup = str_remove_all(plotly_nmds$x$data[[i]]$name,"\\(")
    tmp_replace_legendgroup = str_remove_all(tmp_replace_legendgroup, ",1\\)")
    plotly_nmds$x$data[[i]]$legendgroup = tmp_replace_legendgroup
  }
  htmlwidgets::saveWidget(as_widget(plotly_nmds),file=paste0(name,pred_plot,"_interactive.html"),background="#fafafa",selfcontained=FALSE)
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
  plotly_js = args[9]
  functional_predictions(pred_ec,metadata,criteria,pred_plot,name,microDecon,control, plotly_js)
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
  plotly_js = args[9]
  functional_predictions(pred_ko,metadata,criteria,pred_plot,name,microDecon,control,plotly_js)
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
  plotly_js = args[9]
  functional_predictions(pred_path,metadata,criteria,pred_plot,name,microDecon,control, plotly_js)
}

if (!interactive()) {
  main_metacyc()
}