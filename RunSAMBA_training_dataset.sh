#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/conf/conda_envs/nextflow_env.sh

if [ "$1" != "-resume" ]
then
    if [ ! -d "$BASEDIR/training_dataset" ]
    then 
        mkdir -p $BASEDIR/training_dataset
        wget -r -nc -l2 -nH --cut-dirs=7 ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/training_dataset/q2* -P $BASEDIR/training_dataset
        sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest
        sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest.single
        wget -r -nc -l2 -nH --cut-dirs=7 ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/training_dataset/dna-sequence-raw/*_subsampled.fastq.gz -P $BASEDIR/training_dataset
    fi
    #nextflow temp directory
    export NXF_TEMP=$BASEDIR/.nxf_temp
    mkdir -p $NXF_TEMP
    export tax_db_dir=$BASEDIR/tax.databases.test
    #download taxonomic database
    DB=$tax_db_dir/DATABASE_silva_v132_99_16S.qza
  if [ -f "$DB" ]
  then
    echo "$DB exist"
  else
    mkdir -p $tax_db_dir
    wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/2020.02/SILVA_v132/DATABASE_silva_132_99_16S.qza -O $tax_db_dir/DATABASE_silva_v132_99_16S.qza
  fi
else
  echo "Resume of the previous analysis"
fi

#run nextflow nextmb workflow 
nextflow -trace nextflow.executor run main.nf -profile test,conda $@

#deactivate nextflow environment
. $BASEDIR/conf/conda_envs/delenv.sh
