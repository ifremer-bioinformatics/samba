#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/conf/conda_envs/nextflow_env.sh
export NXF_TEMP=$BASEDIR/.nxf_temp
mkdir -p $NXF_TEMP

#run nextflow nextmb workflow 
nextflow -trace nextflow.executor run main.nf -profile conda,test $@

#deactivate nextflow environment
. $BASEDIR/conf/conda_envs/delenv.sh
