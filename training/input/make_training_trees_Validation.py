#!/usr/bin/env python3

import ROOT
from array import array
import argparse
from collections import Counter



inputFiles = ["/d0/scratch/shin/training_code/PUJetId_training_plotting/training/tree_jmepfnano_qcd_before.root"]


tChain  = ROOT.TChain("Events")
outFile = ROOT.TFile("Validatopm_tree_before.root", "RECREATE")
for inputFile in inputFiles:
    tChain.Add(inputFile)
    print("File name, %s, is added" %inputFile)


outTrees = []

outTrees.append(ROOT.TTree("all","all"))

def book_float_branch(ttree, branch_name, default_value=-999.0):
    branch_array = array("f", [default_value])
    ttree.Branch(branch_name, branch_array, "%s/F" % branch_name)        
    return branch_array

def book_int_branch(ttree, branch_name, default_value=-999):
    branch_array = array("i", [default_value])
    ttree.Branch(branch_name, branch_array, "%s/I" % branch_name)        
    return branch_array

NTrees = len(outTrees)

PV_npvsGood                   = NTrees*[0]
fixedGridRhoFastjetAll        = NTrees*[0]
Jet_pt                   = NTrees*[0]
Jet_eta                  = NTrees*[0]
Jet_puId_dR2Mean         = NTrees*[0]
Jet_puId_majW            = NTrees*[0]
Jet_puId_minW            = NTrees*[0]
Jet_puId_frac01          = NTrees*[0]
Jet_puId_frac02          = NTrees*[0]
Jet_puId_frac03          = NTrees*[0]
Jet_puId_frac04          = NTrees*[0]
Jet_puId_ptD             = NTrees*[0]
Jet_puId_beta            = NTrees*[0]
Jet_puId_pull            = NTrees*[0]
Jet_puId_jetR            = NTrees*[0]
Jet_puId_jetRchg         = NTrees*[0]
Jet_nConstituents        = NTrees*[0]
Jet_puId_nCharged        = NTrees*[0]


for i, outTree in enumerate(outTrees):
    PV_npvsGood                  [i] = book_int_branch  (outTree, "PV_npvsGood"                  )
#    fixedGridRhoFastjetAll       [i] = book_float_branch(outTree, "fixedGridRhoFastjetAll"       )
    Jet_pt                  [i] = book_float_branch(outTree, "Jet_pt"                  )
    Jet_eta                 [i] = book_float_branch(outTree, "Jet_eta"                 )
    Jet_puId_dR2Mean        [i] = book_float_branch(outTree, "Jet_puId_dR2Mean"        )
    Jet_puId_majW           [i] = book_float_branch(outTree, "Jet_puId_majW"           )
    Jet_puId_minW           [i] = book_float_branch(outTree, "Jet_puId_minW"           )
    Jet_puId_frac01         [i] = book_float_branch(outTree, "Jet_puId_frac01"         )
    Jet_puId_frac02         [i] = book_float_branch(outTree, "Jet_puId_frac02"         )
    Jet_puId_frac03         [i] = book_float_branch(outTree, "Jet_puId_frac03"         )
    Jet_puId_frac04         [i] = book_float_branch(outTree, "Jet_puId_frac04"         )
    Jet_puId_ptD            [i] = book_float_branch(outTree, "Jet_puId_ptD"            )
    Jet_puId_beta           [i] = book_float_branch(outTree, "Jet_puId_beta"           )
    Jet_puId_pull           [i] = book_float_branch(outTree, "Jet_puId_pull"           )
    Jet_puId_jetR           [i] = book_float_branch(outTree, "Jet_puId_jetR"           )
    Jet_puId_jetRchg        [i] = book_float_branch(outTree, "Jet_puId_jetRchg"        )
    Jet_nConstituents       [i] = book_int_branch  (outTree, "Jet_nConstituents"  )
    Jet_puId_nCharged       [i] = book_int_branch  (outTree, "Jet_puId_nCharged"       )


def dphi(phi1, phi2):
    dphi = phi1 - phi2
    Pi = ROOT.TMath.Pi()
    if dphi > Pi:
        dphi -= 2 * Pi
    elif dphi <= -Pi:
        dphi += 2 * Pi
    return dphi
nEvent=0

ROOT.gSystem.Load("libFWCoreFWLite.so")

for ievent, event in enumerate(tChain):
    if ievent % 1000 == 0:
        print("processing %s" % ievent)
#    if ievent == 30000: break
    nEvent+=1
    for k in range(event.nJet):

        key = 0

        outTree_toFill = outTrees[0]
        
        PV_npvsGood                        [key][0] = ord(event.PV_npvsGood)
#        fixedGridRhoFastjetAll             [key][0] = event.fixedGridRhoFastjetAll
        Jet_pt                        [key][0] = event.Jet_pt[k]
        Jet_eta                       [key][0] = event.Jet_eta[k]
        Jet_puId_dR2Mean              [key][0] = event.Jet_puId_dR2Mean[k]
        Jet_puId_majW                 [key][0] = event.Jet_puId_majW[k]
        Jet_puId_minW                 [key][0] = event.Jet_puId_minW[k]
        Jet_puId_frac01               [key][0] = event.Jet_puId_frac01[k]
        Jet_puId_frac02               [key][0] = event.Jet_puId_frac02[k]
        Jet_puId_frac03               [key][0] = event.Jet_puId_frac03[k]
        Jet_puId_frac04               [key][0] = event.Jet_puId_frac04[k]
        Jet_puId_ptD                  [key][0] = event.Jet_puId_ptD[k]
        Jet_puId_beta                 [key][0] = event.Jet_puId_beta[k]
        Jet_puId_pull                 [key][0] = event.Jet_puId_pull[k]
        Jet_puId_jetR                 [key][0] = event.Jet_puId_jetR[k]
        Jet_puId_jetRchg              [key][0] = event.Jet_puId_jetRchg[k]
        Jet_nConstituents             [key][0] = ord(event.Jet_nConstituents[k])
        Jet_puId_nCharged             [key][0] = event.Jet_puId_nCharged[k]



        outTree_toFill.Fill()

for outTree in outTrees:
    outTree.Write("", ROOT.TObject.kOverwrite)




outFile.Write("", ROOT.TObject.kOverwrite)
outFile.Close()
