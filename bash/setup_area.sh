# Check Environment setup 
    # correct directory
if [ ! -d "AlexAnalysisPackage" ]; then
    echo -e "\nERROR :: Must run setup script in parent directory of AlexAnalysisPackage.\n"
    return 1
fi
if [ -z "$ATLAS_LOCAL_ROOT_BASE" ]; then
    echo "ERROR :: ATLAS_LOCAL_ROOT_BASE not defined. Add this to enviornment setup:"
    echo "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase"
    echo "alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'"
    return 1
fi
example_tag="n####"
susynt_tag=$1
if [ -z "$1" ] || [ ${#susynt_tag} != ${#example_tag} ]; then
    echo -e "\nERROR :: No susynt-read tag provided as input\n"
    return 1
fi
if [ -d "analysis" ] || [ -d "analysis_run" ]; then
    echo -e "\nERROR :: dependencies already exist. Only run setup in clean environment\n"
    return 1
fi

# TODOs
    # kerberos check

# Setup environment
echo -e "\n\n------------------------------------------------------------"
echo -e "\nSetting up ATLAS software environment"
export ANALYSIS_DIR=$PWD
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

mkdir -p analysis analysis_run
cd $ANALYSIS_DIR/analysis

if [ -d "susynt-read" ] || [ -d "Superflow" ] || [ -d "RestFrames" ]; then
    echo -e "\nERROR :: dependencies already exist. Only run setup in clean environment\n"
    return 1
fi

# Setup susynt-read
echo -e "\n\n------------------------------------------------------------"
echo -e "\nCloning susynt-read"
git clone git@github.com:susynt/susynt-read.git
cd $ANALYSIS_DIR/analysis/susynt-read
echo -e "\nChecking out proper tag"
git checkout $susynt_tag/AA/devel
git_branch=$(git rev-parse --abbrev-ref HEAD)
if [[ "$git_branch" != "$susynt_tag/AA/devel" ]]; then
    echo "ERROR :: Unable to checkout $susynt_tag/AA/devel branch"
    return 1
fi

echo -e "\nSetting up susynt-read"
./bash/setup_area.sh [--stable] 2>&1 |tee setup_area.log

cd $ANALYSIS_DIR/analysis/susynt-read/SusyNtuple
git checkout $susynt_tag/AA/devel

echo -e "\nSetting up susynt-read environment"
cd $ANALYSIS_DIR/analysis/susynt-read/
source bash/setup_release.sh

cd $ANALYSIS_DIR/analysis/
mv susynt-read/RootCore ./
mv susynt-read/SusyNtuple ./
mv susynt-read/SUSYTools ./

# Setup Superflow and RestFrames
echo -e "\n\n------------------------------------------------------------"
echo -e "\nCloning Superflow"
cd $ANALYSIS_DIR/analysis/
git clone git@github.com:alexarmstrongvi/Superflow.git
cd $ANALYSIS_DIR/analysis/Superflow
git remote add upstream git@github.com:dantrim/Superflow.git

echo -e "\nCloning RestFrames"
cd $ANALYSIS_DIR/analysis/
git clone git@github.com:crogan/RestFrames.git

echo -e "\nConfiguring RestFrames (may take a while)"
cd $ANALYSIS_DIR/analysis/RestFrames
./configure
make 
make check
make install
make installcheck

# Setup directory structure
echo -e "\n\n------------------------------------------------------------"
echo -e "\nSetting up directory structure"
cd $ANALYSIS_DIR/analysis/
mkdir -p inputs

cd $ANALYSIS_DIR/analysis_run/
mkdir -p ntuples outputs logs plots yields lists

cd $ANALYSIS_DIR/analysis_run/ntuples
mkdir data mc

mv $ANALYSIS_DIR/AlexAnalysisPackage $ANALYSIS_DIR/analysis/

# Add useful symbolic links
cd $ANALYSIS_DIR/analysis_run/outputs/
ln -s $ANALYSIS_DIR/analysis/AlexAnalysisPackage/submit/submitMakeMiniNtupleToCondor.py ./submitMakeMiniNtupleToCondor.py
#cd $ANALYSIS_DIR
#ln -s $ANALYSIS_DIR/analysis/AlexAnalysisPackage/bash/setup_env.sh ./setup_env.sh


# Configure setup_env.sh
sed -i'' s:X_ANALYSIS_DIR_X:$ANALYSIS_DIR: $ANALYSIS_DIR/analysis/AlexAnalysisPackage/bash/setup_env.sh

# Final setup and compile
echo -e "\n\n------------------------------------------------------------"
echo -e "\nFinal setup and compile"
source $ANALYSIS_DIR/analysis/AlexAnalysisPackage/bash/setup_env.sh
rc find_packages
rc clean
rc compile


cd $ANALYSIS_DIR
return 0
