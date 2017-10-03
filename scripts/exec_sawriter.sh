#!/usr/bin/env bash

# get input params-----------------------------------
SEYMOUR_HOME=$1
REF_FILEPATH=$2


# Setting up SMRTpipe environment
echo "Setting up ENV on $(uname -n) for methylation detection"
source $SEYMOUR_HOME/etc/setup.sh

# run sawriter
sawriter $REF_FILEPATH

# exit-----------------------------------------------
rcode=$?
echo "sa operations finished on $(date) with exit code ${rcode}."
exit ${rcode}