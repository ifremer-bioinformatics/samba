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
## Creation Date: 2019-12-13                                               ####
## Modified on: 2019-12-13                                                 ####
##                                                                           ##
## Copyright (c) SeBiMER, december-2019                                    ####
## Emails: cyril.noel@ifremer.fr and laure.quintric@ifremer.fr             ####
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                                                           ##
## Note: This part contains the function for construction of the barplot     ##   
## 									     ##
## Modified from https://raw.githubusercontent.com/mahendra-mariadassou/     ##
## phyloseq-extended/master/R/graphical_methods.R"                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## return data frame for relative abundance plots in ggplot2
## Return relative abundance of top NumberOfTaxa OTUs at the taxaRank2 level
## within taxaSet1 at taxaRank1 level
ggformat <- function(physeq, taxaRank1 = "Phylum", taxaSet1 = "Proteobacteria",
                     taxaRank2 = "Family", numberOfTaxa = 9) {
    ## Args:
    ## - physeq: phyloseq class object
    ## - taxaRank1: taxonomic level in which to do the first subsetting
    ## - taxaSet1: subset of level taxaRank1 to use
    ## - taxaRank2: taxonomic level used to agglomerate
    ## - numberOfTaxa: number of (top) taxa to keep at level taxaRank2
    ##
    ## Returns:
    ## - data frame with taxonomic, sample and otu counts informations melted together
    ## Enforce orientation and transform count to relative abundances
    stopifnot(!is.null(sample_data(physeq, FALSE)),
              !is.null(tax_table(physeq, FALSE)))
    otutab <- otu_table(physeq)
    if ( !taxa_are_rows(otutab) ) {
    otutab = t(otutab)
    }
    otutab <- as(otutab, "matrix")
    otutab <- apply(otutab, 2, function(x) x / sum(x))
    ## Subset to OTUs belonging to taxaSet1 to fasten process
    stopifnot(all(c(taxaRank1, taxaRank2) %in% colnames(tax_table(physeq))))
    otutab <- otutab[tax_table(physeq)[ , taxaRank1] %in% taxaSet1, , drop = FALSE]
    if (nrow(otutab) == 0) {
        stop(paste("No otu belongs to", paste(taxaSet1, collapse = ","), "\n",
                   "at taxonomic level", taxaRank1))
    }
    mdf <- melt(data = otutab, varnames = c("OTU", "Sample"))
    colnames(mdf)[3] <- "Abundance"
    ## mdf <- mdf[mdf$Abundance > 0, ] ## Remove absent taxa
    ## Add taxonomic information and replace NA and unclassified Unknown
    tax <- as(tax_table(physeq), "matrix")
    tax[is.na(tax)] <- "Unknown"
    tax[tax %in% c("", "unclassified", "Unclassified")] <- "Unknown"
    tax <- data.frame(OTU = rownames(tax), tax)
    mdf <- merge(mdf, tax, by.x = "OTU")
    ## Aggregate by taxaRank2
    mdf <- aggregate(as.formula(paste("Abundance ~ Sample +", taxaRank2)), data = mdf, FUN = sum)
    topTaxa <- aggregate(as.formula(paste("Abundance ~ ", taxaRank2)), data = mdf, FUN = sum)
    ## Keep only numberOfTaxa top taxa and aggregate the rest as "Other"
    topTax <- as.character(topTaxa[ order(topTaxa[ , "Abundance"], decreasing = TRUE), taxaRank2])
    topTax <- topTax[topTax != "Unknown"]
    topTax <- topTax[1:min(length(topTax), numberOfTaxa)]
    ## Change to character
    mdf[ , taxaRank2] <- as.character(mdf[ , taxaRank2])
    ii <- (mdf[ , taxaRank2] %in% c(topTax, "Unknown"))
    mdf[!ii , taxaRank2] <- "Other"
    mdf <- aggregate(as.formula(paste("Abundance ~ Sample +", taxaRank2)), data = mdf, FUN = sum)
    mdf[, taxaRank2] <- factor(mdf[, taxaRank2], levels = c(sort(topTax), "Unknown", "Other"))
    ## Add sample data.frame
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- sample_names(physeq)
    mdf <- merge(mdf, sdf, by.x = "Sample")
    ## Sort the entries by abundance to produce nice stacked bars in ggplot
    mdf <- mdf[ order(mdf[ , taxaRank2], mdf$Abundance, decreasing = TRUE), ]
    return(mdf)
}

## ggplot hue color scale
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

## Plot composition at specific level within the studied kingdom at the kingdom rank
## Restricts plot to x tax (define by the user : nbtax)
composition <- function(PHYLOSEQ, taxaRank1, taxaSet1, taxaRank2, nbtax, fill, group, color_bar, barplot) {
  ggdata = ggformat(PHYLOSEQ, taxaRank1, taxaSet1, taxaRank2, nbtax)
  p = ggplot(ggdata, aes_string(x = "Sample", y = "Abundance", fill = taxaRank2, color = fill))
  ## Manually change color scale to assign grey to "Unknown" (if any)
  if (!is.null(fill) && any(c("Unknown", "Other") %in% unique(ggdata[, fill]))) {
      ranks = as.character(unique(ggdata[, fill]))
      ranks = ranks[ ! ranks %in% c("Unknown", "Other")]
      colvals = c(color_bar, "grey45", "black")
      names(colvals) = c(ranks, "Unknown", "Other")
      ## Now add the manually re-scaled layer with Unassigned as grey
      p = p + scale_fill_manual(values=colvals) + scale_color_manual(values = colvals)
  }
  p = p + geom_bar(stat = "identity", position = "stack") + theme(axis.text.x=element_text(angle=90)) + ggtitle(paste("Composition within", taxaSet1, "(", nbtax, "top", taxaRank2, ")")) + theme_classic() + labs(x="Samples",y="Abundance",fill=fill) + scale_y_continuous(expand=c(0,0),labels=c("0","25","50","75","100")) + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black",size=10)) + theme(axis.title.x=element_text(vjust=-1,color="black",size=13)) + theme(axis.text.y=element_text(hjust=0.8,color="black",size=13)) + theme(axis.title.y=element_text(size=11)) + theme(legend.text=element_text(size=16)) + theme(legend.title=element_text(size=16,face="bold")) + facet_wrap(group, scales = "free_x", nrow = 1) + theme(strip.text=element_text(size=13)) + theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid")) +theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"))

 ggsave(filename=paste(barplot,".svg",sep=""), device="svg", width = 20, height = 12)
 ggsave(filename=paste(barplot,".png",sep=""), device="png", width = 20, height = 12)
}
