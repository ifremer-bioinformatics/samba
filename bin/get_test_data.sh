#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: get test data for testing samba                        ##
##                                                                           ##
###############################################################################

args=("$@")
BASEDIR=${args[0]}
ready=${args[1]}
if [ ! -d "$BASEDIR/training_dataset" ] || ([ -d "$BASEDIR/training_dataset" ] && [ ! "$(ls -A $BASEDIR/training_dataset)" ])
  then 
     mkdir -p $BASEDIR/training_dataset
     wget -r -nc -l2 -nH --cut-dirs=7 -A 'q2*' ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/training_dataset/ -P $BASEDIR/training_dataset
     sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest
     sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest.single
     wget -r -nc -l2 -nH --cut-dirs=7 -A '*_subsampled.fastq.gz' ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/training_dataset/dna-sequence-raw/ -P $BASEDIR/training_dataset
fi
#download taxonomic database
DB=$BASEDIR/tax.databases.test/DATABASE_silva_v132_99_16S.qza
if [ ! -f "$DB" ]
  then
    mkdir -p $BASEDIR/tax.databases.test
    wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/2019.10/SILVA_v132/DATABASE_silva_v132_99_16S.qza -O $BASEDIR/tax.databases.test/DATABASE_silva_v132_99_16S.qza
fi
if ([ -f "$BASEDIR/tax.databases.test/DATABASE_silva_v132_99_16S.qza" ] && [ -f "$BASEDIR/training_dataset/q2_manifest" ])
then
   touch $ready
fi
