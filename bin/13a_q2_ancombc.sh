#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Differential abundance testing with ANCOM-BC           ##
##                                                                           ##
###############################################################################

args=("$@")

ASV_TABLE=${args[0]}
METADATA=${args[1]}
REFERENCE_ENABLED=${args[2]}
REFERENCE_LEVEL=${args[3]}
ADJ_METHOD=${args[4]}
MAX_ITER=${args[5]}
ALPHA=${args[6]}
FORMULA=${args[7]}
ANCOMBC_OUTPUT_QZA=${args[8]}
ANCOMBC_EXPORT=${args[9]}
TAXONOMY=${args[10]}
ASV_TABLE_FAMILY=${args[11]}
ANCOMBC_FAMILY_OUTPUT_QZA=${args[12]}
ANCOMBC_FAMILY_EXPORT=${args[13]}
ASV_TABLE_GENUS=${args[14]}
ANCOMBC_GENUS_OUTPUT_QZA=${args[15]}
ANCOMBC_GENUS_EXPORT=${args[16]}
LOGCMD=${args[17]}

if "${REFERENCE}"
then
    options="--p-reference-levels ${REFERENCE_LEVEL} --p-p-adj-method ${ADJ_METHOD} --p-max-iter ${MAX_ITER} --p-alpha ${ALPHA}"
fi

# Run ANCOM-BC
CMD="qiime composition ancombc --i-table ${ASV_TABLE} --m-metadata-file ${METADATA} --p-formula ${FORMULA} ${options} --o-differentials ${ANCOMBC_OUTPUT_QZA} ;
qiime tools export --input-path ${ANCOMBC_OUTPUT_QZA} --output-path ${ANCOMBC_EXPORT}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Run ANCOM-BC at the family level
CMD="qiime taxa collapse --i-table ${ASV_TABLE} --i-taxonomy ${TAXONOMY} --p-level 5 --o-collapsed-table ${ASV_TABLE_FAMILY} ;
qiime composition ancombc --i-table ${ASV_TABLE_FAMILY} --m-metadata-file ${METADATA} --p-formula ${FORMULA} ${options} --o-differentials ${ANCOMBC_FAMILY_OUTPUT_QZA} ;
qiime tools export --input-path ${ANCOMBC_FAMILY_OUTPUT_QZA} --output-path ${ANCOMBC_FAMILY_EXPORT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Run ANCOM-BC at the genus level
CMD="qiime taxa collapse --i-table ${ASV_TABLE} --i-taxonomy ${TAXONOMY} --p-level 6 --o-collapsed-table ${ASV_TABLE_GENUS} ;
qiime composition ancombc --i-table ${ASV_TABLE_GENUS} --m-metadata-file ${METADATA} --p-formula ${FORMULA} ${options} --o-differentials ${ANCOMBC_GENUS_OUTPUT_QZA} ;
qiime tools export --input-path ${ANCOMBC_GENUS_OUTPUT_QZA} --output-path ${ANCOMBC_GENUS_EXPORT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
