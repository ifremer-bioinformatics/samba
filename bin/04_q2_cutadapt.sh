#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Remove primers in raw sequences                        ##
##                                                                           ##
###############################################################################


args=("$@")

SINGLE_END=${args[0]}
CPUS=${args[1]}
IMPORTED_DATA=${args[2]}
F_PRIMER=${args[3]}
R_PRIMER=${args[4]}
ERROR_RATE=${args[5]}
TRIMMED_DATA_QZA=${args[6]}
TRIMMED_DATA_QZV=${args[7]}
TRIMMED_OUTPUT=${args[8]}
LOGCMD=${args[9]}

# Identify optimal overlap length to minimize false positive primer identification
LENGTH_PRIMER_F=$(echo ${F_PRIMER} | awk '{print length}')
LENGTH_PRIMER_R=$(echo ${R_PRIMER} | awk '{print length}')

if [[ "${LENGTH_PRIMER_F}" -le "${LENGTH_PRIMER_R}" ]]
then 
    OVERLAP=$(( ${LENGTH_PRIMER_F} - 1 ))
else 
    OVERLAP=$(( ${LENGTH_PRIMER_R} - 1 ))
fi

# Run cutadapt
if ${SINGLE_END}
then
    options="qiime cutadapt trim-single --p-front ${F_PRIMER}"
else
    options="qiime cutadapt trim-paired --p-front-f ${F_PRIMER} --p-front-r ${R_PRIMER}" 
fi

CMD="${options} --verbose --p-cores ${CPUS} --i-demultiplexed-sequences ${IMPORTED_DATA} --p-error-rate ${ERROR_RATE} --p-discard-untrimmed --p-match-read-wildcards --p-overlap ${OVERLAP} --o-trimmed-sequences ${TRIMMED_DATA_QZA}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}

# Summarize counts per sample
CMD="qiime demux summarize --verbose --i-data ${TRIMMED_DATA_QZA} --o-visualization ${TRIMMED_DATA_QZV}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}

# Export html report
CMD="qiime tools export --input-path ${TRIMMED_DATA_QZV} --output-path ${TRIMMED_OUTPUT}"
echo ${CMD} >> ${LOGCMD}
eval ${CMD}
