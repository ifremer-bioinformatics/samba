#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/config/conda_envs/nextflow_env.sh

#set test directories
if [ ! -z $TMP ] 
then 
  sed -i 's|/PATH-TO/$projectName|$TMP/SAMBA_results_of_${projectName}|g' config/params.config
  export tax_db_dir=$TMP/tax.databases.test
  sed -i "s|/PATH-TO/database.qza|$TMP/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g" config/params.config
  #nextflow temp directory
  export NXF_TEMP=$TMP
elif [ ! -z $SCRATCH ] 
then
  sed -i 's|/PATH-TO/$projectName|$SCRATCH/SAMBA_results_of_${projectName}|g' config/params.config
  export tax_db_dir=$SCRATCH/tax.databases.test
  sed -i "s|/PATH-TO/database.qza|$SCRATCH/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g" config/params.config
  #nextflow temp directory
  export NXF_TEMP=$SCRATCH
else
  sed -i 's|/PATH-TO/$projectName|${baseDir}/SAMBA_results_of_${projectName}|g' config/params.config
  export tax_db_dir=$BASEDIR/tax.databases.test
  sed -i "s|/PATH-TO/database.qza|$BASEDIR/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g" config/params.config
  #nextflow temp directory
  export NXF_TEMP=$BASEDIR
fi

#download taxonomic database
DB=$tax_db_dir/DATABASE_silva_v132_99_16S.qza
if [ -f "$DB" ]; then
    echo "$DB exist"
else
    mkdir -p $tax_db_dir
    wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/qiime2/2019.07/DATABASE_silva_v132_99_16S.qza -O $tax_db_dir/DATABASE_silva_v132_99_16S.qza
fi

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run SAMBA.nf $1

#deactivate nextflow environment
. $BASEDIR/config/conda_envs/delenv.sh
