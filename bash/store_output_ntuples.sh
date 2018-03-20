#!/bin/bash


dir=$PWD
# Check environment setup
if [ -z "$ANALYSIS_DIR" ]; then
    printf "ERROR :: ANALYSIS_DIR not defined. "
    printf "Make sure to setup environment with setup_env.sh.\n"
    return 1
fi

function usage()
{
    echo -e "Move flat ntuples from outputs into ntuples directory\n"
    echo "./store_output_ntuples.sh"
    echo -e "\t-h --help"
    echo -e "\t--update : move flat ntuples into current directory instead of new directory\n"
}

update=false
while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            return 1
            ;;
        --update)
            update=true
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            return 1
            ;;
    esac
    shift
done

cd $ANALYSIS_DIR/analysis_run/ntuples
if [ "$update" = false ]; then
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
        unlink $group
        ln -s $DIR $group 
    done
    echo "Moving output files"
else
    echo "Updating data and MC directories"
fi
cd $ANALYSIS_DIR/analysis_run
# Order is important
mv outputs/CENTRAL_physics_Main_??????.root ntuples/data/
mv outputs/CENTRAL_??????.root ntuples/mc/

cd $dir
echo "Ready for plotting"
