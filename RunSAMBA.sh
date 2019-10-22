#!/usr/bin/env bash
#activate nextflow environment
. /appli/bioinfo/nextflow/19.07.0/env.sh

#nextflow temp directory
export NXF_TEMP=$SCRATCH

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -trace nextflow.executor run MB.nf $1

#deactivate nextflow environment
. /appli/bioinfo/nextflow/19.07.0/delenv.sh
