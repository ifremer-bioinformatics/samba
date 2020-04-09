#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/conf/conda_envs/nextflow_env.sh

if [ "$1" != "-resume" ]
then
    #nextflow temp directory
    export NXF_TEMP=$BASEDIR/.nfx_temp
    mkdir -p $NXF_TEMP
fi

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run main.nf $1

#deactivate nextflow environment
. $BASEDIR/conf/conda_envs/delenv.sh
