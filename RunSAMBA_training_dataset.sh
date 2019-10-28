#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. $BASEDIR/config/conda_envs/nextflow_env.sh

#nextflow temp directory
export NXF_TEMP=$SCRATCH

#modify path to the training dataset data
sed -i "s|PATH/TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run MB.nf $1

#deactivate nextflow environment
. $BASEDIR/config/conda_envs/delenv.sh
