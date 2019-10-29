#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/config/conda_envs/nextflow_env.sh

#modify path to the training dataset data
sed -i "s|PATH/TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest

#nextflow temp directory
export NXF_TEMP=$SCRATCH

#set test directories
if [ ! -z $TMP ] 
then 
  sed -i 's|/PATH/TO/OUTDIR/$projectName|$TMP/output.test/$projectName|g' config/params.config
  mkdir -p $TMP/tax.databases.test/
  export tax_db_dir=$TMP/tax.databases.test/
  sed -i "s|/PATH/TO/qiime2/2019.07/DATABASE|$TMP/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g" config/params.config
else 
  sed -i 's|/PATH/TO/OUTDIR/$projectName|${basedir}/output.test/$projectName|g' config/params.config
  mkdir -p $BASEDIR/tax.databases.test
  export tax_db_dir=$BASEDIR/tax.databases.test/
  sed -i "s|/PATH/TO/qiime2/2019.07/DATABASE|$BASEDIR/tax.databases.test/DATABASE_silva_v132_99_16S.qza|g" config/params.config
fi

#download taxonomic database
wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/qiime2/2019.07/DATABASE_silva_v132_99_16S.qza -O $tax_db_dir/DATABASE_silva_v132_99_16S.qza 

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run SAMBA.nf $1

#deactivate nextflow environment
. $BASEDIR/config/conda_envs/delenv.sh
