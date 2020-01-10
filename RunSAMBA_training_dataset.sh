#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/config/conda_envs/nextflow_env.sh

if [ "$1" != "-resume" ]
then
  #set test directories
  if [ ! -z $SCRATCH ] 
  then 
    sed -i 's|/PATH-TO/$projectName|$SCRATCH/SAMBA_results_of_${projectName}|g' config/params.config
    sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest
    #nextflow temp directory
    export NXF_TEMP=$SCRATCH
    export tax_db_dir=$SCRATCH/tax.databases.test
    sed -i 's|/PATH-TO/database.qza|$SCRATCH/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g' config/params.config
  elif [  ! -z $TMP  ] 
  then
    sed -i 's|/PATH-TO/$projectName|$TMP/SAMBA_results_of_${projectName}|g' config/params.config
    sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest 
    #nextflow temp directory
    export NXF_TEMP=$TMP
    export tax_db_dir=$TMP/tax.databases.test
    sed -i 's|/PATH-TO/database.qza|$TMP/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g' config/params.config
  else
    sed -i 's|/PATH-TO/$projectName|${baseDir}/SAMBA_results_of_${projectName}|g' config/params.config
    sed -i "s|/PATH-TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest 
    #nextflow temp directory
    export NXF_TEMP=$BASEDIR
    export tax_db_dir=$BASEDIR/tax.databases.test
    sed -i 's|/PATH-TO/database.qza|$BASEDIR/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g' config/params.config
  fi
  #download taxonomic database
  DB=$tax_db_dir/DATABASE_silva_v132_99_16S.qza
  if [ -f "$DB" ]
  then
    echo "$DB exist"
  else
    mkdir -p $tax_db_dir
    wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/qiime2/2019.10/SILVA_v132/DATABASE_silva_v132_99_16S.qza -O $tax_db_dir/DATABASE_silva_v132_99_16S.qza
  fi
else
  echo "Resume of the previous analysis"
fi

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run SAMBA.nf $1

#deactivate nextflow environment
. $BASEDIR/config/conda_envs/delenv.sh
