#!/bin/bash
echo "Starting"
source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh
chmod +x make_training_trees_test.py
./make_training_trees_test.py --era $1 --jet_type $2 --inputFilesList $3 --minJet_pt $4 --maxJet_pt $5 --output_index $6
