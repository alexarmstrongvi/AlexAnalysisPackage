#!/bin/bash

# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. " 
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi
# Setup Grid, Rucio, and FAX
source $ANALYSIS_DIR/analysis/AlexAnalysisPackage/bash/setup_grid.sh
echo -e "\n-------------------------------"
source $ANALYSIS_DIR/analysis/AlexAnalysisPackage/bash/setup_rucio.sh
echo -e "\n-------------------------------"
source $ANALYSIS_DIR/analysis/AlexAnalysisPackage/bash/setup_fax.sh
echo -e "\n-------------------------------"

# Note: No samples are available on rucio so I am unable to test/validate
#       that this code works

## Get available datasets
#cd $ANALYSIS_DIR/analysis/inputs
#if ls $ANALYSIS_DIR/analysis/inputs/*txt &>/dev/null; then
#    echo "Dataset files already exist. Moving on to grouping"
#else
#    python $ANALYSIS_DIR/analysis/susynt-read/python/available_datasets
#fi
#
## Get condor lists
#make_condor_lists="python $ANALYSIS_DIR/analysis/susynt-read/python/make_condor_lists.py"
#condor_sample_groupings="python $ANALYSIS_DIR/analysis/AlexAnalysisPackage/python/condor_sample_groupings.py"
#for txt in $ANALYSIS_DIR/analysis/inputs/*txt 
#do 
#    echo -e "\n-------------------------------"
#    dir=${txt%.txt}
#    mkdir -p $dir
#    eval "$make_condor_lists -i $txt -o $dir" 
#done  
#
## Group sample lists
#for dir in $ANALYSIS_DIR/analysis/inputs/*mc*/
#do
#    eval "$condor_sample_groupings -i $dir"
#done
#
## Get local input file links
#python $ANALYSIS_DIR/analysis/AlexAnalysisPackage/python/make_local_input_lists.py
#
# Check all DSIDs have a cross section defined in SUSYTools


