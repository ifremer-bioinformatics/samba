#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/config/conda_envs/nextflow_env.sh

if [ "$1" != "-resume" ]
then
    sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest 
    #nextflow temp directory
    export NXF_TEMP=$BASEDIR/.nxf_temp
    mkdir -p $NXF_TEMP
    export tax_db_dir=$BASEDIR/tax.databases.test
    sed -i "s|/PATH-TO/database.qza|$BASEDIR/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g" config/params.config
    #download taxonomic database
    DB=$tax_db_dir/DATABASE_silva_v132_99_16S.qza
  if [ -f "$DB" ]
  then
    echo "$DB exist"
  else
    mkdir -p $tax_db_dir
    wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/qiime2/2020.02/SILVA_v132/DATABASE_silva_132_99_16S.qza -O $tax_db_dir/DATABASE_silva_v132_99_16S.qza
  fi
else
  echo "Resume of the previous analysis"
fi

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run SAMBA.nf $1

#deactivate nextflow environment
. $BASEDIR/config/conda_envs/delenv.sh
