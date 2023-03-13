#!/bin/bash

DIR=$( dirname "${BASH_SOURCE[0]}" )

cd $DIR

source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh

export PYTHONUSERBASE=`pwd`

PY_VER=`python -c "import sys; print('python{0}.{1}'.format(*sys.version_info))"`

export PYTHONPATH=$PYTHONUSERBASE/lib/$PY_VER/site-packages:$PYTHONPATH

PATH=$PYTHONUSERBASE/bin:$PATH

echo "now installing pip pkgs"

pip install --user --upgrade pip

pip install --user --upgrade uproot

# for pyroot_cms_scripts
install_version=0.3.1
current_version=`pip show pyroot_cms_scripts | grep Version | cut -d' ' -f 2`

if [ current_version!=install_version ]; then

    if [ ! -d pyroot_cms_scripts ]; then
        git clone https://github.com/jshin96/pyroot_cms_scripts.git --quiet
    fi
    (cd pyroot_cms_scripts
    git pull origin $install_version
    git checkout $install_version --quiet
    pip install --user .
    cd ..)
fi

cd -

echo "Setting Env"

export era=$1
export max_N=$2
export jet_type=$3
export in_dir=$4
export eta_bin=$5
export minJetpt=$6
export maxJetpt=$7
export input_index=$8
export year=$9
export d_name="BDT_${jet_type}_${era}_${year}_NTree500"
mkdir -p output
cd output

mkdir -p $d_name

chmod +x /u/user/shin/scratch/training_code/PUJetId_training_plotting/training/train_bdt_RUN3.py
/u/user/shin/scratch/training_code/PUJetId_training_plotting/training/train_bdt_RUN3.py --era $era --max_N $max_N --jet_type $jet_type --in_dir $in_dir  --d_name  $d_name --eta_bins $eta_bin --minJet_pt $minJetpt --maxJet_pt $maxJetpt --input_index $input_index --year $year

cd ..
