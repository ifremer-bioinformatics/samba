#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: get test data for testing samba                        ##
##                                                                           ##
###############################################################################

args=("$@")
BASEDIR=${args[0]}
ready=${args[1]}
longreads=${args[2]}

if [ "$longreads" = true ]
then
  datatype="longreads"
  datadir="$BASEDIR/training_dataset/$datatype"
  DB="silva_16S_99_k15.mmi"
  TAX="silva_taxonomy_16S_99_majority.txt"
else
  # short reads dataset
  datatype="shortreads"
  datadir="$BASEDIR/training_dataset/$datatype"
  DB="DATABASE_silva_v132_99_16S.qza"
  TAX=""
fi

if [ ! -d "$datadir" ] || ([ -d "$datadir" ] && [ ! "$(ls -A $datadir)" ])
then 
     mkdir -p $datadir
     wget -r -nc -l2 -nH --cut-dirs=88888888 -A 'q2*' ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/training_dataset/$datatype -P $datadir
     sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/$datatype/q2_manifest
     wget -r -nc -l2 -nH --cut-dirs=8 -A '_subsampled.fastq.gz' ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/training_dataset/$datatype/dna-sequence-raw -P $datadir
     if [ ! "$longreads" ] 
     then
        sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/$datatype/q2_manifest.single
     fi
fi
if [ ! -f "$BASEDIR/tax.databases.test/$datatype/$DB" ]
then
    mkdir -p $BASEDIR/tax.databases.test/$datatype
    wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/2019.10/SILVA_v132/$datatype/$DB -O $BASEDIR/tax.databases.test/$datatype/$DB
    if test -z "$TAX"
    then
       wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/2019.10/SILVA_v132/$datatype/$TAX -O $BASEDIR/tax.databases.test/$datatype/$TAX
    fi
fi
if ([ -f "$BASEDIR/tax.databases.test/$datatype/$DB" ] && [ -f "$BASEDIR/training_dataset/$datatype/q2_manifest" ])
then
   touch $ready
fi
