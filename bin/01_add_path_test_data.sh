#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Add user local path for test data                      ##
##                                                                           ##
###############################################################################

args=("$@")
MANIFEST_XLS=${args[0]}
BASEDIR=${args[1]}
MANIFEST_FINAL=${args[2]}
LOGCMD=${args[3]}

# Replace PATH-TO by local path
CMD="sed 's|/PATH-TO|${BASEDIR}|g' ${MANIFEST_XLS} > ${MANIFEST_FINAL}"
echo ${CMD} > ${LOGCMD}
eval ${CMD}
