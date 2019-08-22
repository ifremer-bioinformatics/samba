#!/usr/bin/env bash
#activate nextflow environment
. /appli/bioinfo/nextflow/19.07.0/env.sh

TIMESTAMP=$(date +"%Y%m%d-%H%M%S")

#run nextflow nextmb workflow
nextflow run MB.nf $1 2>&1 $TIMESTAMP-nextflow.log

#deactivate nextflow environment
. /appli/bioinfo/nextflow/19.07.0/delenv.sh
