#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser(
    description="split input files for condor job"
    )
parser.add_argument(
    "--txtname", type=str, default="input_RUN3_2018", help="name of the input txt file with all the input files without txt"
    )
parser.add_argument(
    "--nFiles", type=int, default="10", help="number of files to go in as input for each batch"
    )


args = parser.parse_args()


inputtxt=args.txtname
nFiles=args.nFiles


inputFilesList_o=open("%s.txt" %inputtxt, "r")
inputFiles = inputFilesList_o.read().splitlines()
inputFilesList_o.close()

index=0
for i, name in enumerate(inputFiles):
   if i % nFiles == 0:
      f=open("%s_%i.txt" %(inputtxt,index), "a")
      index=index+1
   f.write("%s \n" %name)
   if i % nFiles == (nFiles-1):
      f.close()

      
