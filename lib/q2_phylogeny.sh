#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Script name: q2_phylogeny.sh                                            ####
##                                                                           ##
## Purpose of script: Build phylogeny of ASVs                                ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2020-01-31                                                 ####
##                                                                           ##
## Email: samba-sebimer@ifremer.fr                                         ####
## 									     ##
## Copyright (c) SeBiMER, august-2019                                      ####
## This program is free software: you can redistribute it and/or modify it   ##
## under the terms of the GNU Affero General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or         ##
## (at your option) any later version.                                       ## 
##									     ##
## License at https://www.gnu.org/licenses/agpl-3.0.txt 		     ##
## 									     ##
## This program is distributed in the hope that it will be useful, but       ##
## WITHOUT ANY WARRANTY; without even the implied warranty of                ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                      ##
## See the GNU Affero General Public License for more details.               ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Command run by snakemake :
### q2_phylogeny.sh ${repseqs_phylo} aligned_repseq.qza masked-aligned_repseq.qza tree.qza tree.log tree_export_dir tree_export.log completecmd > q2_phylogeny.log 2>&1

# Arguments
args=("$@")
repseq=${args[0]}
aligned=${args[1]}
masked=${args[2]}
tree=${args[3]}
tree_log=${args[4]}
export=${args[5]}
tree_export=${args[6]}
logcmd=${args[7]}

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

# tree export
cmd="qiime tools export --input-path $tree --output-path $export >& $tree_export 2>&1"
echo $cmd >> $logcmd
eval $cmd

