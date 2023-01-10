#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh

era=$1
jet_type=$2
data_type=$3
year=$4
period=$5
lumi=$6
N_mc=$7
xs=$8
input_filename=$9
output_filename=$10

./analyze_make_hists_CHS.py \
--era=$era \
--jet_type=$jet_type \
--data_type=$data_type \
--year=$year \
--period=$period \
--lumi=$lumi \
--N_mc=$N_mc \
--input_filename $input_filename \
--output_filename $output_filename

mkdir -p output
mv $output_filename output/$output_filename
