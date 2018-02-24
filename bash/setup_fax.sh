#!/bin/bash

# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. " 
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi

# Make sure grid is setup
source $ANALYSIS_DIR/analysis/AlexAnalysisPackage/bash/setup_grid.sh

# Setup fax
# STORAGEPREFIX is an environment variable that
# should be defined if fax is setup
if [ -z "$STORAGEPREFIX" ]; then
    echo -e "\nSetting up fax"
    localSetupFAX
else
    echo -e "\nFax is already setup"
fi
