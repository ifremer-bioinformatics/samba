#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Automated Statistical Analyses of Metabarcoding Data   ##
##                                                                           ##
###############################################################################

plot.nmds <- function(PHYLOSEQ, ord_nmds, criteria, color_samples, adonis_result, nmds, distance, width, height, graph_title) {
    ## Sample ordination - NMDS ####
    plot_ordination(PHYLOSEQ,ord_nmds,type="samples",color=criteria, title=graph_title) +
      theme_classic() +
      geom_point(size=3) +
      geom_text(aes(label=rownames(sample_data(PHYLOSEQ))),col="black",size=2.5,vjust=2,hjust=1) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      labs(color=criteria) +
      theme(axis.text=element_text(size=12,color="black")) +
      scale_fill_manual(values=alpha(color_samples,0.4)) +
      scale_color_manual(values=color_samples) +
      stat_ellipse(geom="polygon",alpha=0.1,type="t",aes_string(fill=criteria)) +
      labs(caption = paste("Stress:",round(ord_nmds$stress,4),
                           "\nAdonis statistic R:",round(adonis_result$aov.tab$R2[1]*100,2),
                     paste("\nAdonis based on ", criteria,": p-value"),adonis_result$aov.tab$`Pr(>F)`[1],sep=" "))
    ggsave(filename=paste(nmds,distance,".svg",sep=""), device="svg", width = width, height = height)
    ggsave(filename=paste(nmds,distance,".png",sep=""), device="png", width = width, height = height)
}

plot.nmds.interactive <- function(ord_nmds, metadata, Variable, color_samples, adonis_result, nmds, distance, graph_title) {
  nmds1 = ord_nmds$points[,1]
  nmds2 = ord_nmds$points[,2]
  nmds_data = data.frame(NMDS1=nmds1, NMDS2=nmds2, Condition=metadata[,Variable], SampleID = rownames(ord_nmds$points))
  
  ## Sample ordination - NMDS ####
  nmds.interactive = ggplot(nmds_data, aes(x=NMDS1,y=NMDS2), title=graph_title) +
    theme_classic() +
    stat_ellipse(geom="polygon",alpha=0.2,type="t",aes(fill=Condition)) +
    geom_point(shape=19, size=2.5, alpha=0.8, aes(fill=Condition,label = SampleID)) +
    theme(legend.text=element_text(size=13)) +
    theme(legend.title=element_text(size=14)) +
    theme(axis.text=element_text(size=12,color="black")) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values=color_samples) +
    scale_color_manual(values=color_samples)
  
  text_annotation = paste(paste("Stress:",round(ord_nmds$stress,4)), paste0("Adonis R2:"," ",round(adonis_result$aov.tab$R2[1]*100,2),"%"),paste0("p-value: ",
adonis_result$aov.tab$`Pr(>F)`[1]), sep="\n")
  plotly_nmds = ggplotly(nmds.interactive) %>% 
    layout(title=graph_title,margin=list(r=200,l=0,pad=4)) %>%
    add_annotations(x=1.08,y=0.5,text=text_annotation,
                    showarrow=F,xref='paper', yref='paper', xanchor='right', yanchor='auto', xshift=0, yshift=0,font=list(size=13, color="black"))
  for (i in c(1:length(plotly_nmds$x$data))) {
    tmp_replace_name = str_remove_all(plotly_nmds$x$data[[i]]$name,"\\(")
    tmp_replace_name = str_remove_all(tmp_replace_name, ",1\\)")
    plotly_nmds$x$data[[i]]$name = tmp_replace_name
    tmp_replace_legendgroup = str_remove_all(plotly_nmds$x$data[[i]]$name,"\\(")
    tmp_replace_legendgroup = str_remove_all(tmp_replace_legendgroup, ",1\\)")
    plotly_nmds$x$data[[i]]$legendgroup = tmp_replace_legendgroup
  }
  htmlwidgets::saveWidget(as_widget(plotly_nmds),file=paste0(nmds,distance,"_interactive.html"),background="#fafafa",selfcontained=FALSE)
}

plot.pcoa <- function(PHYLOSEQ, ord_pcoa, criteria, color_samples, adonis_result, pcoa, distance, width, height, graph_title) {
    ## Sample ordination - MDS-PCoA ####
    plot_ordination(PHYLOSEQ,ord_pcoa,type="samples",color=criteria, title=graph_title) +
      theme_classic() +
      geom_point(size=3) +
      geom_text(aes(label=rownames(sample_data(PHYLOSEQ))),col="black",size=2.5,vjust=2,hjust=1) +
      theme(legend.text=element_text(size=13)) +
      theme(legend.title=element_text(size=14)) +
      labs(color=criteria) +
      theme(axis.text=element_text(size=12,color="black")) +
      scale_fill_manual(values=alpha(color_samples,0.4)) +
      scale_color_manual(values=color_samples) +
      stat_ellipse(geom="polygon",alpha=0.1,type="t",aes_string(fill=criteria)) +
      labs(caption = paste("\nAdonis statistic R:",round(adonis_result$aov.tab$R2[1]*100,2),
                     paste("\nAdonis based on ", criteria,": p-value"),adonis_result$aov.tab$`Pr(>F)`[1],sep=" "))
    ggsave(filename=paste(pcoa,distance,".svg",sep=""), device="svg", width = width, height = height)
    ggsave(filename=paste(pcoa,distance,".png",sep=""), device="png", width = width, height = height)
}

plot.pcoa.interactive <- function(ord_pcoa, metadata, Variable, color_samples, adonis_result, pcoa, distance, graph_title) {
  mds1 = ord_pcoa$vectors[,1]
  mds2 = ord_pcoa$vectors[,2]
  pcoa_data = data.frame(Axis1=mds1, Axis2=mds2, Condition=metadata[,Variable], SampleID = rownames(ord_pcoa$vectors))
  
  ## Sample ordination - NMDS ####
  pcoa.interactive = ggplot(pcoa_data, aes(x=Axis1,y=Axis2), title=graph_title) +
    theme_classic() +
    stat_ellipse(geom="polygon",alpha=0.2,type="t",aes(fill=Condition)) +
    geom_point(shape=19, size=2.5, alpha=0.8, aes(fill=Condition,label = SampleID)) +
    theme(legend.text=element_text(size=13)) +
    theme(legend.title=element_text(size=14)) +
    theme(axis.text=element_text(size=12,color="black")) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values=color_samples) +
    scale_color_manual(values=color_samples)
  
  text_annotation = paste(paste0("Adonis R2:"," ",round(adonis_result$aov.tab$R2[1]*100,2),"%"),paste0("p-value: ",adonis_result$aov.tab$`Pr(>F)`[1]), sep="\n
")
  plotly_pcoa = ggplotly(pcoa.interactive) %>% 
    layout(title=graph_title,margin=list(r=200,l=0,pad=4)) %>%
    add_annotations(x=1.08,y=0.5,text=text_annotation,
                    showarrow=F,xref='paper', yref='paper', xanchor='right', yanchor='auto', xshift=0, yshift=0,font=list(size=13, color="black"))
  for (i in c(1:length(plotly_pcoa$x$data))) {
    tmp_replace_name = str_remove_all(plotly_pcoa$x$data[[i]]$name,"\\(")
    tmp_replace_name = str_remove_all(tmp_replace_name, ",1\\)")
    plotly_pcoa$x$data[[i]]$name = tmp_replace_name
    tmp_replace_legendgroup = str_remove_all(plotly_pcoa$x$data[[i]]$name,"\\(")
    tmp_replace_legendgroup = str_remove_all(tmp_replace_legendgroup, ",1\\)")
    plotly_pcoa$x$data[[i]]$legendgroup = tmp_replace_legendgroup
  }
  htmlwidgets::saveWidget(as_widget(plotly_pcoa),file=paste0(pcoa,distance,"_interactive.html"),background="#fafafa",selfcontained=FALSE)
}

plot.hc <- function(dendro, group, cols, col_group, method_hc, plot_hc, distance, width, height) {
    ## Sort GROUP color palette according to dend ####
    color = col_group[order.dendrogram(dendro)]
    ## Plot dendrogram ####
    svglite(paste(plot_hc,distance,".svg",sep=""), width = width, height = height)
    plot = dendro %>% set("labels_colors", color) %>% plot(main = paste("Hierarchical clustering with the", method_hc, "method", "based on", distance, "distance", sep=" "))
    legend("topright", legend = levels(group), fill = cols, cex = 0.8, horiz=FALSE, border="white",box.lty=0)
    dev.off()
    png(filename=paste(plot_hc,distance,".png",sep=""), res=150, width = 2000, height = 1200)
    plot = dendro %>% set("labels_colors", color) %>% plot(main = paste("Hierarchical clustering with the", method_hc, "method", "based on", distance, "distance", sep=" "))
    legend("topright", legend = levels(group), fill = cols, cex = 0.8, horiz=FALSE, border="white", box.lty=0)
    dev.off()

}

plot.pie <- function(ExpVar_perc, labels, distance, plot_pie, width, height) {
    svglite(paste(plot_pie,distance,".svg",sep=""), width = width, height = height)
    pie(ExpVar_perc,
      labels=labels,
      clockwise=TRUE,
      radius=1,
      col=brewer.pal(length(ExpVar_perc),"Set1"),
      border="white",
      cex=1,
      main=paste("Percentage of variance explained by each variable for",distance,"matrix",sep=" "))
    dev.off()
    png(filename=paste(plot_pie,distance,".png",sep=""), res=150, width = 2000, height = 1200)
    pie(ExpVar_perc,
      labels=labels,
      clockwise=TRUE,
      radius=1,
      col=brewer.pal(length(ExpVar_perc),"Set1"),
      border="white",
      cex=1,
      main=paste("Percentage of variance explained by each variable for",distance,"matrix",sep=" "))
    dev.off()
}

plot.pie.interactive <- function(pie_data,labels,values,plot_pie,distance) {
  pie.interactive = plot_ly(pie_data, labels=~labels,values=~values,type="pie", 
                                 textposition='auto', textinfo='label+percent', hoverinfo = 'text', 
                                 text = ~paste0("Variable: ",labels,"\n","Percentage of explained variance: ",round(values,2),"%")) %>% 
    layout(showlegend = FALSE)
  htmlwidgets::saveWidget(as_widget(pie.interactive),file=paste0(plot_pie,distance,"_interactive.html"),background="#fafafa",selfcontained=FALSE)
}
