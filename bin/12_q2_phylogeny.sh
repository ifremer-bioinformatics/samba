#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Build phylogeny of ASVs                                ##
##                                                                           ##
###############################################################################

args=("$@")

CPUS=${args[0]}
ASV_SEQS=${args[1]}
ASV_PHYLOGENY_OUTPUT=${args[2]}
ASV_PHYLOGENY_NWK=${args[3]}
LOGCMD=${args[4]}

# Run the phylogeny pipeline
# 1. sequence alignment using MAFFT
# 2. reducing alignment ambiguity (masking)
# 3. construct a phylogeny using fasttree
# 4. export output
CMD="qiime phylogeny align-to-tree-mafft-fasttree --i-sequences ${ASV_SEQS} --output-dir ${ASV_PHYLOGENY_OUTPUT}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Export Newick asv tree
CMD="qiime tools export --input-path ${ASV_PHYLOGENY_OUTPUT}/rooted_tree.qza --output-path newick_export ;
cp newick_export/tree.nwk ${ASV_PHYLOGENY_NWK}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
