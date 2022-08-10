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


args = parser.parse_args()

era = args.era
jet_type = args.jet_type
max_N = args.max_N
in_dir = args.in_dir
d_name = args.d_name
eta_bins = args.eta_bins

minJetpt=args.minJet_pt
maxJetpt=args.maxJet_pt
f_minJetpt=float(minJetpt)
f_maxJetpt=float(maxJetpt)


input_file = ROOT.TFile.Open(
    in_dir + "training_trees_pt"+minJetpt+"To"+maxJetpt+"_%s_%s.root" % (jet_type, era)
)

for eta_bin in eta_bins:
    print("===================")
    print(".")
    print(".")
    print("Starting")
    print("eta bin:" + eta_bin)
    print(".")
    print(".")
    print("===================")
    PromptTree = input_file.Get(eta_bin + "_Prompt")
    PileupTree = input_file.Get(eta_bin + "_Pileup")

    output_file = ROOT.TFile(
        d_name + "/" + "tmva_output_Pt"+minJetpt + "To" + maxJetpt + eta_bin + "_" + jet_type + ".root",
        "RECREATE"
    )

    N = min(PromptTree.GetEntries(), PileupTree.GetEntries())
    N = min(N, max_N)
    N = int(N/2)
    
    #-----------------------------------------------------------
    factory = ROOT.TMVA.Factory(
        "pileupJetId_" + era + "_" + eta_bin + "_" + jet_type,
        output_file,
        "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification"
    )

    # --------------SET 1 For ETA < 3----------------------
    var_set1 = [
        ("PV_npvsGood"             , "I"),
        ("JetPuppi_puId_beta"      , "F"),
        ("JetPuppi_puId_dR2Mean"   , "F"),
        ("JetPuppi_puId_frac01"    , "F"),
        ("JetPuppi_puId_frac02"    , "F"),
        ("JetPuppi_puId_frac03"    , "F"),
        ("JetPuppi_puId_frac04"    , "F"),
        ("JetPuppi_puId_majW"      , "F"),
        ("JetPuppi_puId_minW"      , "F"),
        ("JetPuppi_puId_jetR"      , "F"),
        ("JetPuppi_puId_jetRchg"   , "F"),
        ("JetPuppi_nConstituents"  , "I"),
        ("JetPuppi_puId_nCharged"  , "I"),
        ("JetPuppi_puId_ptD"       , "F"),
        ("JetPuppi_puId_pull"      , "F"),
    ]
    
    # --------------SET 2 For ETA > 3----------------------
    var_set2 = [
        ("PV_npvsGood"      , "I"),
        #("JetPuppi_puId_beta"     , "F"),
        ("JetPuppi_puId_dR2Mean"   , "F"),
        ("JetPuppi_puId_frac01"    , "F"),
        ("JetPuppi_puId_frac02"    , "F"),
        ("JetPuppi_puId_frac03"    , "F"),
        ("JetPuppi_puId_frac04"    , "F"),
        ("JetPuppi_puId_majW"      , "F"),
        ("JetPuppi_puId_minW"      , "F"),
        ("JetPuppi_puId_jetR"      , "F"),
        #("JetPuppi_puId_jetRchg"  , "F"),
        ("JetPuppi_nConstituents"  , "I"),
        #("JetPuppi_puId_nCharged" , "I"),
        ("JetPuppi_puId_ptD"       , "F"),
        ("JetPuppi_puId_pull"      , "F"),
    ]
    
    variables = var_set1
    spectator = [
        ("JetPuppi_pt" , "F"),
        ("JetPuppi_eta", "F"),
    ]
    
    if eta_bin == "Eta3p0To5p0":
        variables = var_set2
    
    loader = ROOT.TMVA.DataLoader(d_name)
    
    for var, var_type in variables:
        loader.AddVariable(var, var_type)
    
    for spec, spec_type in spectator:
        loader.AddSpectator(spec, spec_type)
    
    loader.AddTree(PromptTree, "Signal")
    loader.AddTree(PileupTree, "Background")
    
    loader.PrepareTrainingAndTestTree(
        ROOT.TCut(''),
        "SplitMode=Random:NormMode=NumEvents:!V:"
        + "nTrain_Signal=%d:nTest_Signal=%d:" % (N, N)
        + "nTrain_Background=%d:nTest_Background=%d:" % (N, N)
    )
    
    #-------------------------------------------------------------------
    factory.BookMethod(loader,ROOT.TMVA.Types.kBDT,"BDT","!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1:DoBoostMonitor")
    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()
    
    output_file.cd()
    output_file.Close()
