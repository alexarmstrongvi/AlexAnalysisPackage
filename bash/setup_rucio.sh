#!/bin/bash

# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. " 
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi

# Make sure grid is setup
source $ANALYSIS_DIR/analysis/AlexAnalysisPackage/bash/setup_grid.sh

# Setup rucio
# RUCIO_HOME is an environment variable that
# should be defined if rucio is setup
if [ -z "$RUCIO_HOME" ]; then
    echo -e "\nSetting up rucio"
    lsetup rucio
else
    echo -e "\nRucio is already setup"
fi
