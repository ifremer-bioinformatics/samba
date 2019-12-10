#!/usr/bin/env bash
## Command run by snakemake :
### ${baseDir}/lib/q2_phylogeny.sh ${repseqs_phylo} aligned_repseq.qza masked-aligned_repseq.qza tree.qza tree.log tree_export_dir tree_export.log completecmd > q2_phylogeny.log 2>&1

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

