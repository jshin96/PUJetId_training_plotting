#!/usr/bin/env python3

import ROOT
import os
import argparse

ROOT.gROOT.SetBatch(True)

if "HOME" not in os.environ:
    os.environ["HOME"] = "" # needed for condor batch mode

parser = argparse.ArgumentParser(
    description="train BDT PU JeT ID"
    )
parser.add_argument(
    "--era", type=str, default="106X", help="MC era, like 94X, 102X"
    )
parser.add_argument(
    "--max_N", type=int, default=200000, help="max N entries for testing and training"
    )
parser.add_argument(
    "--jet_type", type=str, default="puppi", help="chs or puppi"
    )
parser.add_argument(
    "--in_dir", type=str, default="", help="directory of root file"
    )
parser.add_argument(
    "--d_name", type=str, default="BDT_puppi_106X", help="dataloader name"
    )
parser.add_argument(
    "--eta_bins", type=str, nargs="+", default=["Eta0p0To2p5","Eta2p5To2p75","Eta2p75To3p0","Eta3p0To5p0"], help="eta bins names in input tree"
    )
parser.add_argument(
    "--minJet_pt", type=str, default="10", help="minimum pt of jets"
    )
parser.add_argument(
    "--maxJet_pt", type=str, default="100", help="maximum pt of jets"
    )
parser.add_argument(
    "--input_index", type=str, default="", help="index for multiple input files"
    )
parser.add_argument(
    "--year", type=str, default="", help="year of dataset"
    )

args = parser.parse_args()

era = args.era
jet_type = args.jet_type
max_N = args.max_N
in_dir = args.in_dir
year=args.year
d_name = args.d_name
eta_bins = args.eta_bins
input_index=args.input_index
i_input_index=int(input_index)+1
minJetpt=args.minJet_pt
maxJetpt=args.maxJet_pt
f_minJetpt=float(minJetpt)
f_maxJetpt=float(maxJetpt)


input_files={}
for i in range(i_input_index):
    input_files["input_file_{0}".format(i)] = ROOT.TFile.Open(in_dir + "/training_trees_pt"+minJetpt+"To"+maxJetpt+"_%s_%s_%i_%s.root" % (jet_type, era, i, year))
#    input_files["input_file_{0}".format(i)] = ROOT.TFile.Open(in_dir + "/training_trees_pt10To100_%s_%s_%i_%s.root" % (jet_type, era, i, year))

#input_file = ROOT.TFile.Open(
#    in_dir + "training_trees_pt"+minJetpt+"To"+maxJetpt+"_%s_%s.root" % (jet_type, era)
#)

for eta_bin in eta_bins:
    print("===================")
    print(".")
    print(".")
    print("Starting")
    print("eta bin:" + eta_bin)
    print(".")
    print(".")
    print("===================")
    PromptTrees={}
    PileupTrees={}
    for j in range(i_input_index):
        PromptTrees["PromptTree_{0}".format(j)] = input_files["input_file_%i" %j].Get(eta_bin + "_Prompt")
        PileupTrees["PileupTree_{0}".format(j)] = input_files["input_file_%i" %j].Get(eta_bin + "_Pileup")

    output_file = ROOT.TFile(
        d_name +"/" + "tmva_output_Pt"+minJetpt + "To" + maxJetpt + eta_bin + "_" + jet_type +".root",
#        d_name +"/" + "tmva_output_Pt"+"40" + "To" + maxJetpt + eta_bin + "_" + jet_type +"GenIdx_Method.root",
        "RECREATE"
    )
    N_Prompt=0
    N_Pileup=0
    for l in range(i_input_index):
#        for a in range(PromptTrees["PromptTree_%i" %l].GetEntries()):
#            PromptTrees["PromptTree_%i" %l].GetEntry(a)
#            if getattr(PromptTrees["PromptTree_%i" %l],"JetCHS_pt") >40:
#                N_Prompt+=1
#        for b in range(PileupTrees["PileupTree_%i" %l].GetEntries()):
#            PileupTrees["PileupTree_%i" %l].GetEntry(b)
#            if getattr(PileupTrees["PileupTree_%i" %l],"JetCHS_pt") >40:
#                N_Pileup+=1
        N_Prompt=N_Prompt+PromptTrees["PromptTree_%i" %l].GetEntries()
        N_Pileup=N_Pileup+PileupTrees["PileupTree_%i" %l].GetEntries()
    print("number of Prompt Jet (%s, %s): %i" %(year, eta_bin, int(N_Prompt)))
    print("number of Pileup Jet (%s, %s): %i" %(year, eta_bin, int(N_Pileup)))
    NE = min(N_Prompt, N_Pileup)
    N = min(NE, max_N)
    N = int(N/2)
    
    #-----------------------------------------------------------
    factory = ROOT.TMVA.Factory(
        "pileupJetId_" + era + "_" + eta_bin + "_" + jet_type,
        output_file,
        "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification"
    )

    if jet_type == "chs":
        # --------------SET 1 For ETA < 3----------------------
        var_set1 = [
            ("PV_npvsGood"                     , "I"),
            ("JetCHS_puId_beta"          , "F"),
            ("JetCHS_puId_dR2Mean"   , "F"),
            ("JetCHS_puId_frac01"        , "F"),
            ("JetCHS_puId_frac02"        , "F"),
            ("JetCHS_puId_frac03"        , "F"),
            ("JetCHS_puId_frac04"        , "F"),
            ("JetCHS_puId_majW"          , "F"),
            ("JetCHS_puId_minW"          , "F"),
            ("JetCHS_puId_jetR"          , "F"),
            ("JetCHS_puId_jetRchg"   , "F"),
            ("JetCHS_nConstituents"  , "I"),
            ("JetCHS_puId_nCharged"  , "I"),
            ("JetCHS_puId_ptD"           , "F"),
            ("JetCHS_puId_pull"          , "F"),
            ("fixedGridRhoFastjetAll"          , "F"),
        ]
        
        # --------------SET 2 For ETA > 3----------------------
        var_set2 = [
            ("PV_npvsGood"          , "I"),
            #("JetCHS_puId_beta"         , "F"),
            ("JetCHS_puId_dR2Mean"   , "F"),
            ("JetCHS_puId_frac01"        , "F"),
            ("JetCHS_puId_frac02"        , "F"),
            ("JetCHS_puId_frac03"        , "F"),
            ("JetCHS_puId_frac04"        , "F"),
            ("JetCHS_puId_majW"          , "F"),
            ("JetCHS_puId_minW"          , "F"),
            ("JetCHS_puId_jetR"          , "F"),
            #("JetCHS_puId_jetRchg"  , "F"),
            ("JetCHS_nConstituents"  , "I"),
            #("JetCHS_puId_nCharged" , "I"),
            ("JetCHS_puId_ptD"           , "F"),
            ("JetCHS_puId_pull"          , "F"),
            ("fixedGridRhoFastjetAll"          , "F"),
        ]
        spectator = [
            ("JetCHS_pt" , "F"),
            ("JetCHS_eta", "F"),
        ]
     

    if jet_type=="puppi":
        # --------------SET 1 For ETA < 3----------------------
        var_set1 = [
            ("PV_npvsGood"                     , "I"),
            ("Jet_puId_beta"          , "F"),
            ("Jet_puId_dR2Mean"   , "F"),
            ("Jet_puId_frac01"        , "F"),
            ("Jet_puId_frac02"        , "F"),
            ("Jet_puId_frac03"        , "F"),
            ("Jet_puId_frac04"        , "F"),
            ("Jet_puId_majW"          , "F"),
            ("Jet_puId_minW"          , "F"),
            ("Jet_puId_jetR"          , "F"),
            ("Jet_puId_jetRchg"   , "F"),
            ("Jet_nConstituents"  , "I"),
            ("Jet_puId_nCharged"  , "I"),
            ("Jet_puId_ptD"           , "F"),
            ("Jet_puId_pull"          , "F"),
            ("fixedGridRhoFastjetAll"          , "F"),
        ]
        
        # --------------SET 2 For ETA > 3----------------------
        var_set2 = [
            ("PV_npvsGood"          , "I"),
            #("Jet_puId_beta"         , "F"),
            ("Jet_puId_dR2Mean"   , "F"),
            ("Jet_puId_frac01"        , "F"),
            ("Jet_puId_frac02"        , "F"),
            ("Jet_puId_frac03"        , "F"),
            ("Jet_puId_frac04"        , "F"),
            ("Jet_puId_majW"          , "F"),
            ("Jet_puId_minW"          , "F"),
            ("Jet_puId_jetR"          , "F"),
            #("Jet_puId_jetRchg"  , "F"),
            ("Jet_nConstituents"  , "I"),
            #("Jet_puId_nCharged" , "I"),
            ("Jet_puId_ptD"           , "F"),
            ("Jet_puId_pull"          , "F"),
            ("fixedGridRhoFastjetAll"          , "F"),
        ]
        spectator = [
            ("Jet_pt" , "F"),
            ("Jet_eta", "F"),
        ]
        
   
    variables = var_set1
    
    if eta_bin == "Eta3p0To5p0":
        variables = var_set2
    
    loader = ROOT.TMVA.DataLoader(d_name)
    
    for var, var_type in variables:
        loader.AddVariable(var, var_type)
    
    for spec, spec_type in spectator:
        loader.AddSpectator(spec, spec_type)
    
    for k in range(i_input_index):
        N_Prompt=0
        N_Pileup=0
        N_Prompt=PromptTrees["PromptTree_%i" %k].GetEntries()
        N_Pileup=PileupTrees["PileupTree_%i" %k].GetEntries()
        if N_Prompt !=0:
            loader.AddTree(PromptTrees["PromptTree_%i" %k], "Signal")
        if N_Pileup !=0:
            loader.AddTree(PileupTrees["PileupTree_%i" %k], "Background")
#    cuts=ROOT.TCut("Jet_pt>%s && Jet_pt<%s" %(minJetpt,maxJetpt))
#    cutb=ROOT.TCut("Jet_pt>%s && Jet_pt<%s" %(minJetpt,maxJetpt))
#    loader.PrepareTrainingAndTestTree(cuts,cutb, "SplitMode=Random:NormMode=NumEvents:!V:nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:" %(N, N, N, N))
    loader.PrepareTrainingAndTestTree(ROOT.TCut(''),"SplitMode=Random:NormMode=NumEvents:!V:nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:" %(N, N, N, N))
    
    #-------------------------------------------------------------------
    factory.BookMethod(loader,ROOT.TMVA.Types.kBDT,"BDT","!H:!V:NTrees=500:BoostType=AdaBoost:Shrinkage=0.1:nCuts=20:DoBoostMonitor")
#    factory.BookMethod(loader,ROOT.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=200:BoostType=Grad:Shrinkage=0.1:nCuts=20:DoBoostMonitor")
#    factory.BookMethod(loader,ROOT.TMVA.Types.kBDT,"BDTB","!H:!V:NTrees=200:BoostType=Bagging:Shrinkage=0.1:nCuts=20:DoBoostMonitor")
#    factory.BookMethod(loader,ROOT.TMVA.Types.kDL,"DNN","!H:V:ErrorStrategy=CROSSENTROPY:Layout=TANH|128,TANH|128,TANH|128,LINEAR:TrainingStrategy=LearningRate=1e-1,Momentum=0.9,Optimizer=ADAM,BatchSize=256:Architecture=CPU")
    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()
    
    output_file.cd()
    output_file.Close()
