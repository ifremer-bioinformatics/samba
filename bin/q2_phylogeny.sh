#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Build phylogeny of ASVs                                ##
##                                                                           ##
###############################################################################

args=("$@")
repseq=${args[0]}
aligned=${args[1]}
masked=${args[2]}
tree=${args[3]}
tree_log=${args[4]}
rooted_tree=${args[5]}
export=${args[6]}
tree_export=${args[7]}
logcmd=${args[8]}

# Alignment of representative sequences
cmd="qiime alignment mafft --i-sequences $repseq --o-alignment $aligned"
echo $cmd > $logcmd
eval $cmd

# Reducing alignment ambiguity: masking and reference alignments
cmd="qiime alignment mask --i-alignment $aligned --o-masked-alignment $masked"
echo $cmd >> $logcmd
eval $cmd

# fasttree phylogeny
cmd="qiime phylogeny fasttree --i-alignment $masked --o-tree $tree >& $tree_log 2>&1"
echo $cmd >> $logcmd
eval $cmd

# Root tree
cmd="qiime phylogeny midpoint-root --i-tree $tree --o-rooted-tree $rooted_tree"
echo $cmd >> $logcmd
eval $cmd

# tree export
cmd="qiime tools export --input-path $rooted_tree --output-path $export >& $tree_export 2>&1"
echo $cmd >> $logcmd
eval $cmd


