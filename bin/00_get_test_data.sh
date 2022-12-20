#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: get test data for testing samba                        ##
##                                                                           ##
###############################################################################

args=("$@")
BASEDIR=${args[0]}
READY=${args[1]}
DATATYPE=${args[2]}
XLS=${args[3]}

if [ "${DATATYPE}" == "longreads" ]
then
  # long reads dataset
  DATADIR="${BASEDIR}/training_dataset/${DATATYPE}"
  DB="silva_v138.1_16S_NR99_SEQ_k15.mmi"
  TAX="silva_v138.1_16S_NR99_TAX.txt"
elif [ "${DATATYPE}" == "shortreads" ]
then
  # short reads dataset
  DATADIR="${BASEDIR}/training_dataset/${DATATYPE}"
  DB="silva_v138.1_16S_99_V4_515F-806RB.qza"
else 
  echo "DATATYPE is incorrect"
  exit 1 
fi

if [ ! -d "${DATADIR}" ] || ([ -d "${DATADIR}" ] && [ ! "$(ls -A ${DATADIR})" ])
then 
     mkdir -p ${DATADIR}
     wget -r -nc -l2 -nH --cut-dirs=9 ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/v4/training_dataset/${DATATYPE}/excel_sample_file.xls -P ${DATADIR}
     wget -r -nc -l2 -nH --cut-dirs=9 -A '_001.fastq.gz' ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/v4/training_dataset/${DATATYPE}/dna-sequence-raw -P ${DATADIR}
fi
if [ ! -f "${BASEDIR}/tax.databases.test/${DATATYPE}/${DB}" ]
then
    mkdir -p ${BASEDIR}/tax.databases.test/${DATATYPE}
    wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/v4/training_dataset/${DATATYPE}/${DB} -O ${BASEDIR}/tax.databases.test/${DATATYPE}/${DB}
    if [ "${DATATYPE}" == "longreads" ] 
    then
       wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/v4/training_dataset/${DATATYPE}/${TAX} -O ${BASEDIR}/tax.databases.test/${DATATYPE}/${TAX}
    fi
fi

if ([ -f "${BASEDIR}/tax.databases.test/${DATATYPE}/${DB}" ] && [ -f "${BASEDIR}/training_dataset/${DATATYPE}/excel_sample_file.xls" ])
then
   touch ${READY}
   ln -s "${BASEDIR}/training_dataset/${DATATYPE}/excel_sample_file.xls" "${XLS}"
fi
