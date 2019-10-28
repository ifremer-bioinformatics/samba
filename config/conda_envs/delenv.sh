#!/usr/bin/env bash
conda deactivate
PATH=$(echo "$PATH" | sed -e 's@/appli/anaconda/3.7/bin/:@@g')