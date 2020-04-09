#!/usr/bin/env bash
# Deactivate nextflow conda environment
conda deactivate
# Remove condabin from PATH (Modify the path according to your local installation)
PATH=$(echo "$PATH" | sed -e 's@/appli/anaconda/versions/.*/condabin:@@g')
