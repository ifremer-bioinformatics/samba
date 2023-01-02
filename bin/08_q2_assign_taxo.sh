#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Taxonomic assignation of ASVs                          ##
##                                                                           ##
###############################################################################

args=("$@")

Q2_TMP_DIR=${args[0]}
CPUS=${args[1]}
CONFIDENCE=${args[2]}
DATABASE=${args[3]}
ASV_SEQS=${args[4]}
ASSIGNED_TAXO_QZA=${args[5]}
ASSIGNED_TAXO_QZV=${args[6]}
ASSIGN_TAXO_OUTPUT=${args[7]}
TAX_ASSIGN_TSV=${args[8]}
ASV_OUTDIR=${args[9]}
ASV_TAX_TABLE_BIOM=${args[10]}
ASV_TAX_TABLE_TSV=${args[11]}
LOGCMD=${args[12]}

mkdir ${Q2_TMP_DIR}
export TMPDIR=${Q2_TMP_DIR}

# Run the Naive Bayesian Classifier for taxonomy assignment
CMD="qiime feature-classifier classify-sklearn --p-n-jobs ${CPUS} --p-confidence ${CONFIDENCE} --i-classifier ${DATABASE} --i-reads ${ASV_SEQS} --o-classification ${ASSIGNED_TAXO_QZA}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Generate visual and tabular view of the taxonomy assignation
CMD="qiime metadata tabulate --m-input-file ${ASSIGNED_TAXO_QZA} --o-visualization ${ASSIGNED_TAXO_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export all data to an QIIME 2 html report
CMD="qiime tools export --input-path ${ASSIGNED_TAXO_QZV} --output-path ${ASSIGN_TAXO_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Rename tabular taxonomy file and modify header
CMD="cp ${ASSIGN_TAXO_OUTPUT}/metadata.tsv ${TAX_ASSIGN_TSV} ;
sed -i '1,2d' ${TAX_ASSIGN_TSV} ;
sed -i '1 i\\#OTUID\ttaxonomy\tconfidence' ${TAX_ASSIGN_TSV} ; 
cp ${TAX_ASSIGN_TSV} ${ASSIGN_TAXO_OUTPUT}/"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Add taxonomy to count table (biom format)
CMD="biom add-metadata -i ${ASV_OUTDIR}/feature-table.biom --observation-metadata-fp ${TAX_ASSIGN_TSV} -o ${ASV_TAX_TABLE_BIOM} --sc-separated taxonomy"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert biom table to tabular
CMD="biom convert -i ${ASV_TAX_TABLE_BIOM} -o ${ASV_TAX_TABLE_TSV} --to-tsv --header-key taxonomy ;
cp ${ASV_TAX_TABLE_TSV} ${ASSIGN_TAXO_OUTPUT}/"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
