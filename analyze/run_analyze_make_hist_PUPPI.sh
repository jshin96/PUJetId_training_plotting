#!/bin/bash

DIR=$( dirname "${BASH_SOURCE[0]}" )

cd $DIR

source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh

export PYTHONUSERBASE=`pwd`

PY_VER=`python -c "import sys; print('python{0}.{1}'.format(*sys.version_info))"`

export PYTHONPATH=$PYTHONUSERBASE/lib/$PY_VER/site-packages:$PYTHONPATH

PATH=$PYTHONUSERBASE/bin:$PATH

echo "now installing pip pkgs"

#pip install --user --upgrade pip

pip install --user --upgrade uproot

# for pyroot_cms_scripts
install_version=0.3.1
current_version=`pip show pyroot_cms_scripts | grep Version | cut -d' ' -f 2`

if [ current_version!=install_version ]; then

    if [ ! -d pyroot_cms_scripts ]; then
        git clone https://github.com/singh-ramanpreet/pyroot_cms_scripts.git --quiet
    fi
    (cd pyroot_cms_scripts
    git pull origin $install_version
    git checkout $install_version --quiet
    pip install --user .
    cd ..)
fi

cd -

#era=$1
#jet_type=$2
#data_type=$3
#year=$4
#period=$5
#lumi=$6
#N_mc=$7
#xs=$8
#input_filename=$9
#output_filename=${10}


#./analyze_make_hists_PUPPI.py \
#--era=$era \
#--jet_type=$jet_type \
#--data_type=$data_type \
#--year=$year \
#--period=$period \
#--lumi=$lumi \
#--N_mc=$N_mc \
#--xs=$xs \
#--input_filename=$input_filename \
#--output_filename=$output_filename

export era=106X
export jet_type=puppi
export data_type=mc
export year=2017
export period=full
export lumi=41.54
export N_mc=102941427
export xs=5343.0
export input_filename=dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/shin/PU_Training_2017_with_weights/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NanoTestPost/220907_041726/0000/tree_12.root
export output_filename=dy_amcatnlo_${data_type}_${year}${period}.root

./analyze_make_hists_PUPPI.py --era=$era --jet_type=$jet_type --data_type=$data_type --year=$year --period=$period --lumi=$lumi --N_mc=$N_mc --xs=$xs --input_filename=$input_filename --output_filename=$output_filename

mkdir -p output_$year

#mv $output_filename output_$year/$output_filename
mv $output_filename output_$year/
