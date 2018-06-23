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
    echo -e "\t-a : Add to name of output directories [mc_date -> mc_<input>_data]\n"
    echo -e "\t-n : numerator samples\n"
    echo -e "\t-d : denominator samples\n"
    echo -e "\t--update : move flat ntuples into current directory instead of new directory\n"
}

update=false
name_mod=""
num=false
den=false
while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    echo "$PARAM -> $VALUE"
    case $PARAM in
        -h | --help)
            usage
            return 1
            ;;
        -a)
            name_mod="_${VALUE}"
            ;;
        --update)
            update=true
            ;;
        -n)
            num=true
            ;;
        -d)
            den=true
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
        DIR="${group}${name_mod}_${DATE}/"
        echo "$DIR"
        if [ -d "$DIR" ] && [ "$(ls -A ${DIR})" ]; then
            echo "$DIR is already filled"
            cd $dir
            return 1
        fi
        mkdir -p $DIR
        link_name="${group}"
        if [ "$num" = true ]; then
            link_name="${link_name}_num"
        elif [ "$den" = true ]; then
            link_name="${link_name}_den"
        fi
        unlink $link_name
        ln -s $DIR $link_name 
    done
    echo "Moving output files"
else
    echo "Updating data and MC directories"
fi
cd $ANALYSIS_DIR/analysis_run
# Order is important
if [ "$num" = true ]; then
    suffix="_num"
elif [ "$den" = true ]; then
    suffix="_den"
else
    suffix=""
fi
mv outputs/CENTRAL_physics_Main_*.root ntuples/data${suffix}/

mv outputs/CENTRAL_*.root ntuples/mc${suffix}/
rm outputs/RunCondorSF.sh outputs/submitFile_TEMPLATE.condor

cd $dir
echo "Ready for plotting"
