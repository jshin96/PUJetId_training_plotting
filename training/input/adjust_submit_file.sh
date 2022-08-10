#!/bin/bash

export nFiles=$1
cp submit_make_trees_default.sub submit_make_trees.sub
for i in $(seq 0 $nFiles);do if [ $i == $nFiles ]; then sed -i "s/index.txt/$i.txt\noutput_index=$i/g" submit_make_trees.sub; else sed -i "s/index.txt/$i.txt\noutput_index=$i\nqueue\n\ninput=input_index.txt/g" submit_make_trees.sub;fi; done
