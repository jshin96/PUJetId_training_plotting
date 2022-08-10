#!/bin/bash
echo "Setting Env"
source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh

export era=$1
export max_N=$2
export jet_type=$3
export in_dir=$4
export d_name="BDT_${jet_type}_${era}"
export eta_bin=$5
export minJetpt=$6
export maxJetpt=$7


mkdir -p output
cd output

mkdir -p $d_name

chmod +x /u/user/shin/scratch/training_code/PUJetId_training_plotting/training/train_bdt_test.py
/u/user/shin/scratch/training_code/PUJetId_training_plotting/training/train_bdt_test.py --era $era --max_N $max_N --jet_type $jet_type --in_dir $in_dir  --d_name  $d_name --eta_bins $eta_bin --minJet_pt $minJetpt --maxJet_pt $maxJetpt

cd ..
