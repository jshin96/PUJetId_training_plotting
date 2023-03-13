#!/bin/bash

export nFiles=$1
export year=$2
cp submit_make_trees_PUPPI_default.sub submit_make_trees_PUPPI.sub
for i in $(seq 0 $nFiles);do if [ $i == $nFiles ]; then sed -i "s/index.txt/$i.txt\noutput_index=$i/g" submit_make_trees_PUPPI.sub; else sed -i "s/index.txt/$i.txt\noutput_index=$i\nqueue\n\ninput=input_year_index.txt/g" submit_make_trees_PUPPI.sub;fi; done
sed -i "s/year/${year}/g" submit_make_trees_PUPPI.sub


cp submit_make_trees_CHS_default.sub submit_make_trees_CHS.sub
for i in $(seq 0 $nFiles);do if [ $i == $nFiles ]; then sed -i "s/index.txt/$i.txt\noutput_index=$i/g" submit_make_trees_CHS.sub; else sed -i "s/index.txt/$i.txt\noutput_index=$i\nqueue\n\ninput=input_year_index.txt/g" submit_make_trees_CHS.sub;fi; done
sed -i "s/year/${year}/g" submit_make_trees_CHS.sub

cp submit_make_trees_PUPPI_CHS_matching_default.sub submit_make_trees_PUPPI_CHS_matching.sub
for i in $(seq 0 $nFiles);do if [ $i == $nFiles ]; then sed -i "s/index.txt/$i.txt\noutput_index=$i/g" submit_make_trees_PUPPI_CHS_matching.sub; else sed -i "s/index.txt/$i.txt\noutput_index=$i\nqueue\n\ninput=input_year_index.txt/g" submit_make_trees_PUPPI_CHS_matching.sub;fi; done
sed -i "s/year/${year}/g" submit_make_trees_PUPPI_CHS_matching.sub

cp submit_make_trees_PUPPI_PFCands_default.sub submit_make_trees_PUPPI_PFCands.sub
for i in $(seq 0 $nFiles);do if [ $i == $nFiles ]; then sed -i "s/index.txt/$i.txt\noutput_index=$i/g" submit_make_trees_PUPPI_PFCands.sub; else sed -i "s/index.txt/$i.txt\noutput_index=$i\nqueue\n\ninput=input_PFCands_index.txt/g" submit_make_trees_PUPPI_PFCands.sub;fi; done
sed -i "s/year/${year}/g" submit_make_trees_PUPPI_PFCands.sub

cp submit_make_trees_PUPPI_RUN3_default.sub submit_make_trees_PUPPI_RUN3.sub
for i in $(seq 0 $nFiles);do if [ $i == $nFiles ]; then sed -i "s/index.txt/$i.txt\noutput_index=$i/g" submit_make_trees_PUPPI_RUN3.sub; else sed -i "s/index.txt/$i.txt\noutput_index=$i\nqueue\n\ninput=input_RUN3_year_index.txt/g" submit_make_trees_PUPPI_RUN3.sub;fi; done
sed -i "s/year/${year}/g" submit_make_trees_PUPPI_RUN3.sub
