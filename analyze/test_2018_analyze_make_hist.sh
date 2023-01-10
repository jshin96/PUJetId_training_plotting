#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh

era="106X"
data_type="mc"
#data_type="data"
year="2018"
period="full"
#period="D"
lumi=59.96
N_mc=96004933
xs=5321.0
input_filename="/u/user/shin/SE_UserHome/PU_Training_Samples/2018/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_0.root"
output_filename="output_hists_${data_type}_${year}${period}.root"

./analyze_make_hists.py \
--era=$era \
--data_type=$data_type \
--year=$year \
--period=$period \
--lumi=$lumi \
--N_mc=$N_mc \
--input_filename $input_filename \
--output_filename $output_filename
