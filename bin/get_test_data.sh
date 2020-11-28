#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: get test data for testing samba                        ##
##                                                                           ##
###############################################################################

args=("$@")
BASEDIR=${args[0]}
ready=${args[1]}
datatype=${args[2]}
manifest=${args[3]}
metadata=${args[4]}

if [ "$datatype" == "longreads" ]
then
  datadir="$BASEDIR/training_dataset/$datatype"
  DB="silva_16S_99_k15.mmi"
  TAX="silva_taxonomy_16S_99_majority.txt"
elif [ "$datatype" == "shortreads" ]
then
  # short reads dataset
  datadir="$BASEDIR/training_dataset/$datatype"
  DB="DATABASE_silva_v132_99_16S_V4_515F-806R.qza"
  TAX=""
else 
  echo "Datatype is incorrect"
  exit 1 
fi

if [ ! -d "$datadir" ] || ([ -d "$datadir" ] && [ ! "$(ls -A $datadir)" ])
then 
     mkdir -p $datadir
     wget -r -nc -l2 -nH --cut-dirs=8 -A 'q2*' ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/training_dataset/$datatype -P $datadir
     sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/$datatype/q2_manifest
     wget -r -nc -l2 -nH --cut-dirs=8 -A '_subsampled.fastq.gz' ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/training_dataset/$datatype/dna-sequence-raw -P $datadir
     if [ "$datatype" == "shortreads" ] 
     then
        sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/$datatype/q2_manifest.single
     fi
fi
if [ ! -f "$BASEDIR/tax.databases.test/$datatype/$DB" ]
then
    mkdir -p $BASEDIR/tax.databases.test/$datatype
    wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/2019.10/SILVA_v132/$datatype/$DB -O $BASEDIR/tax.databases.test/$datatype/$DB
    if [ "$datatype" == "longreads" ] 
    then
       wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/2019.10/SILVA_v132/$datatype/$TAX -O $BASEDIR/tax.databases.test/$datatype/$TAX
    fi
fi
if ([ -f "$BASEDIR/tax.databases.test/$datatype/$DB" ] && [ -f "$BASEDIR/training_dataset/$datatype/q2_manifest" ])
then
   touch $ready
   ln -s "$BASEDIR/training_dataset/$datatype/q2_manifest" "$manifest"
   ln -s "$BASEDIR/training_dataset/$datatype/q2_metadata" "$metadata"
fi
