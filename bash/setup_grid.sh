#!/bin/bash

# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. " 
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi

# Setup grid proxy
# X509_USER_PROXY is an environment variable that
# should be defined if grid proxy has been setup
time_left=`voms-proxy-info | grep timeleft`
hours_left=${time_left:12:2}
#min_left=${time_left:15:2}
#sec_left=${time_left:18:2}

if [ -z "$time_left" ]; then
    echo -e "\nSetting up grid proxy"
    voms-proxy-init -voms atlas -valid 96:00
elif [ "$hours_left" == "00" ]; then
    echo -e "\nRenewing Grid Proxy ($X509_USER_PROXY -> $time_left)"
    voms-proxy-init -voms atlas -valid 96:00
else
    echo -e "\nGrid proxy is already setup ($X509_USER_PROXY -> $time_left)"
fi

return 0
