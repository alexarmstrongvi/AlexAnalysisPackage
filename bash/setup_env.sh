# Setup Atlas software
export ANALYSIS_DIR='X_ANALYSIS_DIR_X'
dir=$PWD
cd $ANALYSIS_DIR
if [ -z "$alrb_AvailableToolsPre" ]; then
    setupATLAS
fi

# Define PATH variables
if [[ $PYTHONPATH != *"AlexAnalysisPackage"* ]]; then
    export PYTHONPATH="$PYTHONPATH:$ANALYSIS_DIR/analysis/AlexAnalysisPackage/python/"
    export PYTHONPATH="$PYTHONPATH:$ANALYSIS_DIR/analysis/AlexAnalysisPackage/plotting/"
fi

# Setup Root and RootCore
cd $ANALYSIS_DIR/analysis/
source $ANALYSIS_DIR/analysis/susynt-read/bash/setup_release.sh

cd $dir
