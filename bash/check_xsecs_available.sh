#!/bin/bash
cur_dir=$PWD

# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. "
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi

dirs_to_check="$ANALYSIS_DIR/analysis/SUSYTools/data/mc15_13TeV/*"
#dirs_to_check="$dirs_to_check /cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/mc15_13TeV/*"

dsids_searched=0
dsids_found=0
while read dsid; do
    if [[ -z "$dsid" ]] || [[ $dsid == \#* ]]; then continue; fi
    ((dsids_searched++))
    flag=false
    for dir in $dirs_to_check; do
        CMD="grep -r $dsid $dir"
        if [[ ! -z "$($CMD)" ]]; then
            #echo "FOUND: $($CMD)"
            if [[ "$flag" = true ]]; then
                echo "DSID xsec info duplicated in multiple dirs: $dsid"
                echo "$(grep -rl $dsid $dirs_to_check)" 
            else
                ((dsids_found++))
            fi
            flag=true
            #break
        fi  
    done
    if [[ "$flag" = false ]]; then
        echo "No Xsec found for DSID: $dsid"
    fi
done < $1
echo "${dsids_found}/${dsids_searched} DSID cross-sections found" 
echo "=== COMPLETED ==="

echo $cur_dir
