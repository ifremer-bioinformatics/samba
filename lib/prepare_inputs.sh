#!/usr/bin/env bash
## Command run by nextflow :
### prepare_input.sh ${manifest} > data_integrity.log 2>&1

# Arguments
args=("$@")
manifest=${args[0]}

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}")" ; cd ..  >/dev/null && pwd )"

# Format the manifest input file
sed -i "s|PATH/TO|$BASEDIR|g" $BASEDIR/training_dataset/q2_manifest

