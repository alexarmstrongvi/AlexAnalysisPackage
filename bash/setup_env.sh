# Setup Atlas software
export ANALYSIS_DIR='/data/uclhc/uci/user/armstro1/SusyNt/SusyNt_n0235_LFV_analysis'
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
rootver=6.04.16-x86_64-slc6-gcc49-opt
echo ""
echo "Setting up ROOT ${rootver}"
lsetup "root ${rootver} --skipConfirm"

# Setup python and desired packages
lsetup python
lsetup "lcgenv -p LCG_93 x86_64-slc6-gcc62-dbg pandas"

# if rootcore is already set up, clean up the env
echo ""
echo "Setting up RootCore"
if [ -d "${ROOTCOREDIR}" ]; then
    source ${ROOTCOREDIR}/scripts/unsetup.sh
fi
source RootCoreBin/local_setup.sh

#rc find_packages
#rc clean
#rc compile

cd $dir
