#!/usr/bin/env bash
conda deactivate
PATH=$(echo "$PATH" | sed -e 's@/appli/anaconda/latest/bin/:@@g')
