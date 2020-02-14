#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/config/conda_envs/nextflow_env.sh

if [ "$1" != "-resume" ]
then
    sed -i 's|/PATH-TO/$projectName|${baseDir}/SAMBA_results_of_${projectName}|g' config/params.config
    #nextflow temp directory
    export NXF_TEMP=$BASEDIR/.nfx_temp
    mkdir -p $NXF_TEMP
fi

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run SAMBA.nf $1

#deactivate nextflow environment
. $BASEDIR/config/conda_envs/delenv.sh
