#!/usr/bin/env python3

import ROOT
from array import array
import argparse

parser = argparse.ArgumentParser(
    description="make TTrees for training"
    )
parser.add_argument(
    "--era", type=str, default="106X", help="MC era, like 94X, 102X"
    )
parser.add_argument(
    "--jet_type", type=str, default="puppi", help="jet type puppi or CHS"
    )
parser.add_argument(
    "--inputFilesList", type=str, default="", help="txt file with the full path to input files in each line"
    )
parser.add_argument(
    "--minJet_pt", type=str, default="10", help="minimum pt of jets"
    )
parser.add_argument(
    "--maxJet_pt", type=str, default="9999", help="maximum pt of jets"
    )
parser.add_argument(
    "--output_index", type=str, default="", help="index for multiple output files"
    )


args = parser.parse_args()

era = args.era
jet_type = args.jet_type
inputFilesList = args.inputFilesList
minJetpt=args.minJet_pt
maxJetpt=args.maxJet_pt
f_minJetpt=float(minJetpt)
f_maxJetpt=float(maxJetpt)
output_index=args.output_index


print("Jet Transverse momentum range is (%s , %s)" %(minJetpt, maxJetpt))

inputFilesList_o=open(inputFilesList, "r")
inputFiles = inputFilesList_o.read().splitlines()
inputFilesList_o.close()




tChain  = ROOT.TChain("Events")
outFile = ROOT.TFile("training_trees_pt"+minJetpt+"To"+maxJetpt+"_%s_%s_%s.root" % (jet_type, era,output_index), "RECREATE")

for inputFile in inputFiles:
    tChain.Add(inputFile)
    print("File name, %s, is added" %inputFile)

eta_bins = [
    "Eta0p0To2p5" ,
    "Eta2p5To2p75",
    "Eta2p75To3p0",
    "Eta3p0To5p0" ,
]

outTrees = []

for eta_bin in eta_bins:
    outTrees.append(ROOT.TTree(eta_bin + "_Prompt", eta_bin + "_Prompt"))
    
for eta_bin in eta_bins:
    outTrees.append(ROOT.TTree(eta_bin + "_Pileup", eta_bin + "_Pileup"))

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
JetPuppi_pt                   = NTrees*[0]
JetPuppi_eta                  = NTrees*[0]
JetPuppi_puId_dR2Mean         = NTrees*[0]
JetPuppi_puId_majW            = NTrees*[0]
JetPuppi_puId_minW            = NTrees*[0]
JetPuppi_puId_frac01          = NTrees*[0]
JetPuppi_puId_frac02          = NTrees*[0]
JetPuppi_puId_frac03          = NTrees*[0]
JetPuppi_puId_frac04          = NTrees*[0]
JetPuppi_puId_ptD             = NTrees*[0]
JetPuppi_puId_beta            = NTrees*[0]
JetPuppi_puId_pull            = NTrees*[0]
JetPuppi_puId_jetR            = NTrees*[0]
JetPuppi_puId_jetRchg         = NTrees*[0]
JetPuppi_nConstituents        = NTrees*[0]
JetPuppi_puId_nCharged        = NTrees*[0]
JetPuppi_dRMatch              = NTrees*[0]
JetPuppi_genJetIdx            = NTrees*[0]
JetPuppi_partonFlavour        = NTrees*[0]


for i, outTree in enumerate(outTrees):
    PV_npvsGood                  [i] = book_int_branch  (outTree, "PV_npvsGood"                  )
    fixedGridRhoFastjetAll       [i] = book_float_branch(outTree, "fixedGridRhoFastjetAll"       )
    JetPuppi_pt                  [i] = book_float_branch(outTree, "JetPuppi_pt"                  )
    JetPuppi_eta                 [i] = book_float_branch(outTree, "JetPuppi_eta"                 )
    JetPuppi_puId_dR2Mean        [i] = book_float_branch(outTree, "JetPuppi_puId_dR2Mean"        )
    JetPuppi_puId_majW           [i] = book_float_branch(outTree, "JetPuppi_puId_majW"           )
    JetPuppi_puId_minW           [i] = book_float_branch(outTree, "JetPuppi_puId_minW"           )
    JetPuppi_puId_frac01         [i] = book_float_branch(outTree, "JetPuppi_puId_frac01"         )
    JetPuppi_puId_frac02         [i] = book_float_branch(outTree, "JetPuppi_puId_frac02"         )
    JetPuppi_puId_frac03         [i] = book_float_branch(outTree, "JetPuppi_puId_frac03"         )
    JetPuppi_puId_frac04         [i] = book_float_branch(outTree, "JetPuppi_puId_frac04"         )
    JetPuppi_puId_ptD            [i] = book_float_branch(outTree, "JetPuppi_puId_ptD"            )
    JetPuppi_puId_beta           [i] = book_float_branch(outTree, "JetPuppi_puId_beta"           )
    JetPuppi_puId_pull           [i] = book_float_branch(outTree, "JetPuppi_puId_pull"           )
    JetPuppi_puId_jetR           [i] = book_float_branch(outTree, "JetPuppi_puId_jetR"           )
    JetPuppi_puId_jetRchg        [i] = book_float_branch(outTree, "JetPuppi_puId_jetRchg"        )
    JetPuppi_nConstituents       [i] = book_int_branch  (outTree, "JetPuppi_nConstituents"  )
    JetPuppi_puId_nCharged       [i] = book_int_branch  (outTree, "JetPuppi_puId_nCharged"       )
    JetPuppi_dRMatch             [i] = book_float_branch(outTree, "JetPuppi_dRMatch"             )
    JetPuppi_genJetIdx           [i] = book_int_branch  (outTree, "JetPuppi_genJetIdx"           )
    JetPuppi_partonFlavour       [i] = book_int_branch  (outTree, "JetPuppi_partonFlavour"       )


#dR of Puppi Jet and Gen Jet
def dRGenMatch(event, index):
   PuppiJet_eta = event.JetPuppi_eta[index]
   PuppiJet_phi = event.JetPuppi_phi[index]
   Gen_index = event.JetPuppi_genJetIdx[index]
   Matched_GenJet_eta = event.GenJet_eta[Gen_index]
   Matched_GenJet_phi = event.GenJet_phi[Gen_index]
   dRMatch = ((PuppiJet_eta-Matched_GenJet_eta)**2+(PuppiJet_phi-Matched_GenJet_phi)**2)**0.5
   return dRMatch
   


for ievent, event in enumerate(tChain):
    if ievent % 10000 == 0:
        print("processing %s" % ievent)
#    if ievent > 40000: break
    nLeptons = event.nElectron + event.nMuon
    #stop if there is no 2 leptons, nor there is no jet
    if nLeptons != 2: continue
    if event.nJetPuppi == 0: continue
    for i in range(event.nJetPuppi):
        PuppiJet_eta = event.JetPuppi_eta[i]
        PuppiJet_phi = event.JetPuppi_phi[i]
        Gen_index = event.JetPuppi_genJetIdx[i]
        if Gen_index != -1:
            Matched_GenJet_eta = event.GenJet_eta[Gen_index]
            Matched_GenJet_phi = event.GenJet_phi[Gen_index]
            dRMatch_ = ((PuppiJet_eta-Matched_GenJet_eta)**2+(PuppiJet_phi-Matched_GenJet_phi)**2)**0.5
        eta_     = event.JetPuppi_eta[i]
        flavor_  = event.JetPuppi_partonFlavour[i]

        #Jet pt upper and lower bound condition
        if event.JetPuppi_pt[i] > f_maxJetpt or event.JetPuppi_pt[i] < f_minJetpt: continue
        
        isPrompt = False
        isPileup = False
        
        if (dRMatch_ <= 0.2):
            isPrompt = True
        
        if (dRMatch_ >= 0.4 and abs(flavor_) == 0 ):
            isPileup = True
            
        key = ""
        
        if (isPrompt and (        abs(eta_) <= 2.5 )) : key = 0
        if (isPrompt and ( 2.5  < abs(eta_) <= 2.75)) : key = 1
        if (isPrompt and ( 2.75 < abs(eta_) <= 3.0 )) : key = 2
        if (isPrompt and ( 3.0  < abs(eta_) <= 5.0 )) : key = 3
        
        if (isPileup and (        abs(eta_) <= 2.5 )) : key = 4
        if (isPileup and ( 2.5  < abs(eta_) <= 2.75)) : key = 5
        if (isPileup and ( 2.75 < abs(eta_) <= 3.0 )) : key = 6
        if (isPileup and ( 3.0  < abs(eta_) <= 5.0 )) : key = 7

        # nothing matched
        if key == "": continue

        outTree_toFill = outTrees[key]
        
        PV_npvsGood                        [key][0] = event.PV_npvsGood
        fixedGridRhoFastjetAll             [key][0] = event.fixedGridRhoFastjetAll
        JetPuppi_pt                        [key][0] = event.JetPuppi_pt[i]
        JetPuppi_eta                       [key][0] = event.JetPuppi_eta[i]
        JetPuppi_puId_dR2Mean              [key][0] = event.JetPuppi_puId_dR2Mean[i]
        JetPuppi_puId_majW                 [key][0] = event.JetPuppi_puId_majW[i]
        JetPuppi_puId_minW                 [key][0] = event.JetPuppi_puId_minW[i]
        JetPuppi_puId_frac01               [key][0] = event.JetPuppi_puId_frac01[i]
        JetPuppi_puId_frac02               [key][0] = event.JetPuppi_puId_frac02[i]
        JetPuppi_puId_frac03               [key][0] = event.JetPuppi_puId_frac03[i]
        JetPuppi_puId_frac04               [key][0] = event.JetPuppi_puId_frac04[i]
        JetPuppi_puId_ptD                  [key][0] = event.JetPuppi_puId_ptD[i]
        JetPuppi_puId_beta                 [key][0] = event.JetPuppi_puId_beta[i]
        JetPuppi_puId_pull                 [key][0] = event.JetPuppi_puId_pull[i]
        JetPuppi_puId_jetR                 [key][0] = event.JetPuppi_puId_jetR[i]
        JetPuppi_puId_jetRchg              [key][0] = event.JetPuppi_puId_jetRchg[i]
        JetPuppi_nConstituents             [key][0] = ord(event.JetPuppi_nConstituents[i])
        JetPuppi_puId_nCharged             [key][0] = event.JetPuppi_puId_nCharged[i]
        JetPuppi_dRMatch                   [key][0] = dRMatch_
        JetPuppi_genJetIdx                 [key][0] = event.JetPuppi_genJetIdx[i]
        JetPuppi_partonFlavour             [key][0] = event.JetPuppi_partonFlavour[i]



        outTree_toFill.Fill()

for outTree in outTrees:
    outTree.Write("", ROOT.TObject.kOverwrite)

outFile.Write("", ROOT.TObject.kOverwrite)
outFile.Close()
