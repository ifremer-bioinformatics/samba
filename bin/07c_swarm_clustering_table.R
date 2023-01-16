#!/usr/bin/env Rscript
###############################################################################
##                                                                           ##
## Purpose of script: Reconstruct ASV table from ASV clustered by swarm      ##
##                                                                           ##
###############################################################################

## Load up the needed packages ####
requiredPackages = c("reshape2","dplyr","stringr")
for(package in requiredPackages){
  library(package,character.only = TRUE)
}

# @@@@@@@@@@@@@@@@@@@@@@@@ #
#                          #
# Decontamination analysis #
#                          #
# @@@@@@@@@@@@@@@@@@@@@@@@ #
args = commandArgs(trailingOnly=TRUE)

ASV_TABLE = args[1]
SWARM_CLUSTER_LIST =  args[2]
CLUSTER_ID = args[3]
ASV_TABLE_CLUSTERED = args[4]

# Load ASV table
IMPORTED_TABLE = read.table(ASV_TABLE, h=T, sep="\t", check.names=FALSE)

# Load cluster list file
DF_CLUSTER = read.table(SWARM_CLUSTER_LIST, h=F, sep=" ", dec=".", fill=TRUE)

# Reformat cluster dataframe
DF_CLUSTER = data.frame(apply(DF_CLUSTER,2,function(x)gsub('_.*', '',x)))
DF_CLUSTER$clusterID = paste0("cluster_", 1:length(DF_CLUSTER$V1))
ASV_REF_ID = data.frame(DF_CLUSTER$V1, DF_CLUSTER$clusterID)
colnames(ASV_REF_ID) = c("ASV_ID","clusterID")
MELT_DF_CLUSTER = melt(DF_CLUSTER, id.vars="clusterID")
MELT_DF_CLUSTER = MELT_DF_CLUSTER %>% select(-variable)
MELT_DF_CLUSTER = MELT_DF_CLUSTER[!MELT_DF_CLUSTER$value=="", ]
colnames(MELT_DF_CLUSTER) = c("clusterID","ASV_ID")
MELT_DF_CLUSTER = MELT_DF_CLUSTER[order(MELT_DF_CLUSTER$clusterID),]

# Construct the ASV table of clustered ASV
DF_TABLE_CLUSTER = merge(IMPORTED_TABLE, MELT_DF_CLUSTER, by.x="ASV_ID", by.y="ASV_ID", all.x=T)
row.names(DF_TABLE_CLUSTER) = DF_TABLE_CLUSTER$ASV_ID
DF_TABLE_CLUSTER = DF_TABLE_CLUSTER %>% select(-ASV_ID)
DF_TABLE_CLUSTER_SUMMED = aggregate(.~clusterID, DF_TABLE_CLUSTER, sum)

FINAL_TABLE_CLUSTERED = merge(DF_TABLE_CLUSTER_SUMMED, ASV_REF_ID, by.x="clusterID", by.y="clusterID")
FINAL_TABLE_CLUSTERED = FINAL_TABLE_CLUSTERED %>% select(-clusterID) %>% select("ASV_ID", everything())

# Save the new ASV table
write.table(MELT_DF_CLUSTER, file=CLUSTER_ID, row.names=F, col.names=T, quote=F, sep="\t")
write.table(FINAL_TABLE_CLUSTERED, file=ASV_TABLE_CLUSTERED, row.names=F, col.names=T, quote=F, sep="\t")
