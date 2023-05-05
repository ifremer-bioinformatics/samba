#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: All Dada2 processes                                    ##
##                                                                           ##
###############################################################################

args=("$@") 

SINGLE_END=${args[0]}
DADA2_INPUT_DATA=${args[1]}
FIGARO_ENABLE=${args[2]}
RAW_LENGTH=${args[3]}
PRIMER_F=${args[4]}
PRIMER_R=${args[5]}
F_TRIM_LEFT=${args[6]}
R_TRIM_LEFT=${args[7]}
F_TRUNC_LENGTH=${args[8]}
R_TRUNC_LENGTH=${args[9]}
TRUNC_Q=${args[10]}
F_MAX_EE=${args[11]}
R_MAX_EE=${args[12]}
N_READS_LEARN=${args[13]}
POOLING_METHOD=${args[14]}
CHIMERA_METHOD=${args[15]}
CPUS=${args[16]}
REP_SEQS_QZA=${args[17]}
TABLE_QZA=${args[18]}
DADA2_STATS_QZA=${args[19]}
DADA2_STATS_QZV=${args[20]}
TABLE_QZV=${args[21]}
METADATA=${args[22]}
REP_SEQS_QZV=${args[23]}
DADA2_OUTPUT=${args[24]}
LOGCMD=${args[25]}

if ${FIGARO_ENABLE}
then
    LENGTH_PRIMER_F=$(echo ${PRIMER_F} | awk '{print length}')
    LENGTH_PRIMER_R=$(echo ${PRIMER_R} | awk '{print length}')
    F_LENGTH_TEST=$((F_TRUNC_LENGTH + LENGTH_PRIMER_F))
    R_LENGTH_TEST=$((R_TRUNC_LENGTH + LENGTH_PRIMER_R)) 
    if [ "${F_LENGTH_TEST}" -ge "${RAW_LENGTH}" ]
    then
        F_TRUNC_LENGTH=0
    else
        F_TRUNC_LENGTH=${F_TRUNC_LENGTH}
    fi
    if ["${R_LENGTH_TEST}" -ge "${RAW_LENGTH}"]
    then
        R_TRUNC_LENGTH=0
    else
        R_TRUNC_LENGTH=${R_TRUNC_LENGTH}
    fi
fi

if ${SINGLE_END}
then
    options="qiime dada2 denoise-single --p-trim-left ${F_TRIM_LEFT} --p-trunc-len ${F_TRUNC_LENGTH} --p-max-ee ${F_MAX_EE}"
else
    options="qiime dada2 denoise-paired --p-trim-left-f ${F_TRIM_LEFT} --p-trim-left-r ${R_TRIM_LEFT} --p-trunc-len-f ${F_TRUNC_LENGTH} --p-trunc-len-r ${R_TRUNC_LENGTH} --p-max-ee-f ${F_MAX_EE} --p-max-ee-r ${R_MAX_EE}"
fi

# Run DADA2
CMD="${options} --verbose --i-demultiplexed-seqs ${DADA2_INPUT_DATA} --p-trunc-q ${TRUNC_Q} --p-n-reads-learn ${N_READS_LEARN} --p-pooling-method ${POOLING_METHOD} --p-chimera-method ${CHIMERA_METHOD} --p-n-threads ${CPUS} --o-representative-sequences ${REP_SEQS_QZA} --o-table ${TABLE_QZA} --o-denoising-stats ${DADA2_STATS_QZA}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Generate a tabular view of DADA2 process stats
CMD="qiime metadata tabulate --verbose --m-input-file ${DADA2_STATS_QZA} --o-visualization ${DADA2_STATS_QZV}"
echo ${CMD} >> ${LOGCMD} 
eval ${CMD}

# Generate visual and tabular summaries of the DADA2 ASV table
CMD="qiime feature-table summarize --verbose --i-table ${TABLE_QZA} --o-visualization ${TABLE_QZV} --m-sample-metadata-file ${METADATA}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Generate tabular view of ASV sequences
CMD="qiime feature-table tabulate-seqs --verbose --i-data ${REP_SEQS_QZA} --o-visualization ${REP_SEQS_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export all data to an QIIME 2 html report
CMD="qiime tools export --input-path ${REP_SEQS_QZV} --output-path ${DADA2_OUTPUT} ; 
qiime tools export --input-path ${TABLE_QZA} --output-path ${DADA2_OUTPUT} ;
qiime tools export --input-path ${TABLE_QZV} --output-path ${DADA2_OUTPUT} ;
qiime tools export --input-path ${DADA2_STATS_QZV} --output-path ${DADA2_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Convert BIOM table to TSV
CMD="biom convert -i ${DADA2_OUTPUT}/feature-table.biom -o ${DADA2_OUTPUT}/dada2_asv_table.tsv --to-tsv ;
sed -i '1d' ${DADA2_OUTPUT}/dada2_asv_table.tsv ;
sed -i 's/#OTU ID/ASV_ID/g' ${DADA2_OUTPUT}/dada2_asv_table.tsv"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Associate abundance to each ASV sequence
# -> for further header modification of the ASV FASTA file
CMD="sed '1d' ${DADA2_OUTPUT}/feature-frequency-detail.csv > asv_abundance_header.txt ;
sed -i 's/,/_/g' asv_abundance_header.txt ;
sed -i 's/.0$//g' asv_abundance_header.txt ;
cut -d '_' -f1 asv_abundance_header.txt > asv_id.txt"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Modify the ASV sequence header of the ASV FASTA file
CMD="cp ${DADA2_OUTPUT}/sequences.fasta sequences_with_abundance.fasta ;
sed 's/^/sed -i \"s\//g' asv_id.txt > tmp_rename_ASV_sequence_header.sh ;
paste -d '\/' tmp_rename_ASV_sequence_header.sh asv_abundance_header.txt > rename_ASV_sequence_header.sh ;
sed -i 's/$/\/g\" sequences_with_abundance.fasta/g' rename_ASV_sequence_header.sh ;
sh rename_ASV_sequence_header.sh"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
