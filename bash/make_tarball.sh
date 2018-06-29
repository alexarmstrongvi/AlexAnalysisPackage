#!/bin/bash

dir=$PWD
# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. " 
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi

# store tar option settings not as a string
exclude_ops=`eval echo '--exclude=".*" --exclude="*git*"'`

# tar analysis file
tar cfvz $ANALYSIS_DIR/analysis.tgz ${exclude_ops} $ANALYSIS_DIR/analysis/

# create symbolic link if not present
ln -fs $ANALYSIS_DIR/analysis.tgz $ANALYSIS_DIR/area.tgz 

echo -e "\n Tar file created"
echo "-----------------------------"

cd $dir
