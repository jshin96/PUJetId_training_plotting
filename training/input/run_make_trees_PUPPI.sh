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
	git clone https://github.com/jshin96/pyroot_cms_scripts.git --quiet
    fi
    (cd pyroot_cms_scripts
    git pull origin $install_version
    git checkout $install_version --quiet
    pip install --user .
    cd ..)
fi

cd -



echo "Starting"

chmod +x make_training_trees_PUPPI.py
./make_training_trees_PUPPI.py --era $1 --year $2 --jet_type $3 --inputFilesList $4 --minJet_pt $5 --maxJet_pt $6 --output_index $7
