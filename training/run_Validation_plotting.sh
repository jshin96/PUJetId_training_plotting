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


mkdir Validation_plots
cd Validation_plots
mkdir $2
cd $2
mkdir $1
cd $1


chmod +x /u/user/shin/scratch/training_code/PUJetId_training_plotting/training/plot_input_variables_Validation.py
/u/user/shin/scratch/training_code/PUJetId_training_plotting/training/plot_input_variables_Validation.py --jet_type $1 --process_name $2 --eta_s $3 --eta_min $4 --eta_max $5
cd ..
