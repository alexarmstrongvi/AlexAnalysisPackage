#!/bin/bash

dir=$PWD
# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. "
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi


cd $ANALYSIS_DIR/analysis_run/ntuples
echo "Creating data and MC directories"
for group in "data" "mc"; do
    DATE=`date +%Y_%m_%d`
    DIR="${group}_${DATE}/"
    if [ -d "$DIR" ] && [ "$(ls -A ${DIR})" ]; then
        echo "$DIR is already filled"
        cd $dir
        return 1
    fi
    mkdir -p $DIR
    ln -s $DIR $group 
done

echo "Moving output files"
cd $ANALYSIS_DIR/analysis_run
# Order is important
mv outputs/CENTRAL_physics_Main_??????.root ntuples/data/
mv outputs/CENTRAL_??????.root ntuples/mc/

cd $dir
echo "Ready for plotting"
