#!/bin/bash

# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. " 
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi

# Setup grid proxy
# X509_USER_PROXY is an environment variable that
# should be defined if grid is already setup
if [ -z "$X509_USER_PROXY" ]; then
    echo -e "\nSetting up grid proxy"
    voms-proxy-init -voms atlas -valid 96:00
else
    echo -e "\nGrid proxy is already setup ($X509_USER_PROXY)"
fi




















































