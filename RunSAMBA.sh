#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/config/conda_envs/nextflow_env.sh

#nextflow temp directory
if [ ! -z $TMP ] 
then 
  sed -i 's|/PATH-TO/$projectName|$TMP/SAMBA_results_of_${projectName}|g' config/params.config
  #nextflow temp directory
  export NXF_TEMP=$TMP
elif [ ! -z $SCRATCH ] 
then
  sed -i 's|/PATH-TO/$projectName|$SCRATCH/SAMBA_results_of_${projectName}|g' config/params.config
  #nextflow temp directory
  export NXF_TEMP=$SCRATCH
else
  sed -i 's|/PATH-TO/$projectName|${baseDir}/SAMBA_results_of_${projectName}|g' config/params.config
  #nextflow temp directory
  export NXF_TEMP=$BASEDIR
fi

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run SAMBA.nf $1

#deactivate nextflow environment
. $BASEDIR/config/conda_envs/delenv.sh
