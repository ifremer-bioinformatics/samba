#!/usr/bin/env bash
## Command run by snakemake :
### ${baseDir}/lib/q2_phylogeny.sh repseqs aligned_repseq masked_aligned tree model ${params.phylogeny.alrt} ${params.phylogeny.bootstrap} tree_bestmodel tree_log completecmd > q2_phylogeny.log 2>&1

# Arguments
args=("$@")
repseq=${args[0]}
aligned=${args[1]}
masked=${args[2]}
tree=${args[3]}
model=${args[4]}
alrt=${args[5]}
bootstrap=${args[6]}
tree_bestmodel=${args[7]}
tree_log=${args[8]}
logcmd=${args[9]}

# Alignment of representative sequences
cmd="qiime alignment mafft --i-sequences $repseq --o-alignment $aligned"
echo $cmd > $logcmd
eval $cmd

# Reducing alignment ambiguity: masking and reference alignments
cmd="qiime alignment mask --i-alignment $aligned --o-masked-alignment $masked"
echo $cmd >> $logcmd
eval $cmd

# IQ-tree phylogeny
cmd="qiime phylogeny iqtree-ultrafast-bootstrap --i-alignment $masked --p-n-cores 0 --o-tree $tree --verbose >& $model 2>&1"
echo $cmd >> $logcmd
eval $cmd

best_model="$(grep 'Best-fit model' $model_log | cut -d ' ' -f3)"

cmd="qiime phylogeny iqtree-ultrafast-bootstrap --i-alignment $masked --p-n-cores 0 --p-alrt $alrt --p-abayes --p-lbp $bootstrap --p-substitution-model $best_model --o-tree $tree_bestmodel --verbose >& $tree_log 2>&1"
echo $cmd >> $logcmd
eval $cmd

