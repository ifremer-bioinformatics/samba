#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Functional predictions using PICRUSt2                  ##
##                                                                           ##
###############################################################################

args=("$@")

ASV_SEQS_FASTA=${args[0]}
ASV_TABLE_BIOM=${args[1]}
PICRUST2_OUTPUT=${args[2]}
CPUS=${args[3]}
TRAITS_DB=${args[4]}
NSTI=${args[5]}
HSP_METHOD=${args[6]}
MIN_READS=${args[7]}
MIN_SAMPLES=${args[8]}
TOP_PATHWAY=${args[9]}
SEC_PATHWAY=${args[10]}
LOGCMD=${args[11]}

# PICRUSt2 pipeline
CMD="picrust2_pipeline.py -s ${ASV_SEQS_FASTA} -i ${ASV_TABLE_BIOM} -o ${PICRUST2_OUTPUT} -p ${CPUS} --in_traits ${TRAITS_DB} --max_nsti ${NSTI} -m ${HSP_METHOD} --min_reads ${MIN_READS} --min_samples ${MIN_SAMPLES} --coverage --remove_intermediate ;
add_descriptions.py -i ${PICRUST2_OUTPUT}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o ${PICRUST2_OUTPUT}/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz ;
add_descriptions.py -i ${PICRUST2_OUTPUT}/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO -o ${PICRUST2_OUTPUT}/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz ;
add_descriptions.py -i ${PICRUST2_OUTPUT}/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o ${PICRUST2_OUTPUT}/pathways_out/path_abun_unstrat_descrip.tsv.gz;
add_descriptions.py -i ${PICRUST2_OUTPUT}/pathways_out/path_abun_unstrat.tsv.gz --custom_map_tabl ${TOP_PATHWAY} -o ${PICRUST2_OUTPUT}/pathways_out/path_abun_unstrat_descrip_top_level.tsv.gz;
add_descriptions.py -i ${PICRUST2_OUTPUT}/pathways_out/path_abun_unstrat.tsv.gz --custom_map_tabl ${SEC_PATHWAY} -o ${PICRUST2_OUTPUT}/pathways_out/path_abun_unstrat_descrip_secondary_level.tsv.gz"
echo ${CMD} > ${LOGCMD}
eval ${CMD}
