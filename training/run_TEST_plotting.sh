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

export eta_s=$1
export eta_f=$2
export process=$3
export dz_name=$4

mkdir -p compare_plots
cd compare_plots
mkdir -p ${process}_${dz_name}
cd ${process}_${dz_name}


chmod +x /u/user/shin/scratch/training_code/PUJetId_training_plotting/training/plot_input_variables_TEST_puId_PFCands.py
/u/user/shin/scratch/training_code/PUJetId_training_plotting/training/plot_input_variables_TEST_puId_PFCands.py --eta_s $eta_s --eta_f $eta_f --process $process --dz $dz_name
cd ..
