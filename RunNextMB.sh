#!/usr/bin/env bash
#activate nextflow environment
. /appli/bioinfo/nextflow/19.07.0/env.sh

TIMESTAMP=$(date +"%Y%m%d-%H%M%S")

#nextflow temp directory
export NXF_TEMP=$SCRATCH

#run nextflow nextmb workflow
#nextflow run MB.nf $1 > $TIMESTAMP-nextflow.log 2>&1
nextflow -trace nextflow.executor run MB.nf $1

#deactivate nextflow environment
. /appli/bioinfo/nextflow/19.07.0/delenv.sh
