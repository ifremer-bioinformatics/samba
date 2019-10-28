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
  sed "s|outdir = "/PATH/TO/OUTDIR|outdir = "$TMP/output.test|g"
  tax.db.test="$TMP/tax.databases.test"
else 
  sed "s|outdir = "/PATH/TO/OUTDIR|outdir = "$BASEDIR/output.test|g"
  tax.db.test="$BASEDIR/tax.databases.test"
fi

#download taxonomic database
wget ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/qiime2/2019.07/DATABASE_silva_v132_99_16S.qza -O $tax.db.test

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run SAMBA.nf $1

#deactivate nextflow environment
. $BASEDIR/config/conda_envs/delenv.sh
