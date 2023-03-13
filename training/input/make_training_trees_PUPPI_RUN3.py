#!/usr/bin/env python3

import numpy as np
from numpy.linalg import eig
import ROOT
from array import array
import argparse
from collections import Counter
parser = argparse.ArgumentParser(
    description="make TTrees for training"
    )
parser.add_argument(
    "--era", type=str, default="106X", help="MC era, like 94X, 102X, 106X, 130X"
    )
parser.add_argument(
    "--year", type=str, default="", help="sample year"
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
year = args.year

print("Jet Transverse momentum range is (%s , %s)" %(minJetpt, maxJetpt))

inputFilesList_o=open(inputFilesList, "r")
inputFiles = inputFilesList_o.read().splitlines()
inputFilesList_o.close()




tChain  = ROOT.TChain("Events")
outFile = ROOT.TFile("EE_Comparison_training_trees_pt"+minJetpt+"To"+maxJetpt+"_%s_%s_%s_%s.root" % (jet_type, era,output_index,year), "RECREATE")
for inputFile in inputFiles:
    tChain.Add(inputFile)
#    event.Add(inputFile)
    print("File name, %s, is added" %inputFile)

eta_bins = [
    "Eta0p0To2p5" ,
    "Eta2p5To2p75",
    "Eta2p75To3p0",
    "Eta3p0To5p0" ,
]

outTrees = []
PFCandTrees = []
for eta_bin in eta_bins:
    outTrees.append(ROOT.TTree(eta_bin + "_Prompt", eta_bin + "_Prompt"))
    PFCandTrees.append(ROOT.TTree(eta_bin + "_Prompt_PFCands", eta_bin + "_Prompt_PFCands"))
    
for eta_bin in eta_bins:
    outTrees.append(ROOT.TTree(eta_bin + "_Pileup", eta_bin + "_Pileup"))
    PFCandTrees.append(ROOT.TTree(eta_bin + "_Pileup_PFCands", eta_bin + "_Pileup_PFCands"))

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
Jet_dRMatch              = NTrees*[0]
Jet_genJetIdx            = NTrees*[0]
Jet_partonFlavour        = NTrees*[0]
Zpt_Jetpt                = NTrees*[0]
dphi_zj                  = NTrees*[0]
Jet_PFCand_dR               = NTrees*[0]
Jet_PFCand_deta               = NTrees*[0]
Jet_PFCand_dphi               = NTrees*[0]
Jet_PFCand_Absptdeta               = NTrees*[0]
Jet_PFCand_Absptdphi               = NTrees*[0]
Jet_PFCand_Frcptdeta               = NTrees*[0]
Jet_PFCand_Frcptdphi               = NTrees*[0]
Jet_PFCand_puppiWeight               = NTrees*[0]
Jet_PFCand_puppiWeightNoLep               = NTrees*[0]
Jet_PFCand_trkChi2               = NTrees*[0]
Jet_PFCand_vtxChi2               = NTrees*[0]

for i, outTree in enumerate(outTrees):
    PV_npvsGood                  [i] = book_int_branch  (outTree, "PV_npvsGood"                  )
    fixedGridRhoFastjetAll       [i] = book_float_branch(outTree, "fixedGridRhoFastjetAll"       )
    Jet_pt                       [i] = book_float_branch(outTree, "Jet_pt"                  )
    Jet_eta                      [i] = book_float_branch(outTree, "Jet_eta"                 )
    Jet_puId_dR2Mean             [i] = book_float_branch(outTree, "Jet_puId_dR2Mean"        )
    Jet_puId_majW                [i] = book_float_branch(outTree, "Jet_puId_majW"           )
    Jet_puId_minW                [i] = book_float_branch(outTree, "Jet_puId_minW"           )
    Jet_puId_frac01              [i] = book_float_branch(outTree, "Jet_puId_frac01"         )
    Jet_puId_frac02              [i] = book_float_branch(outTree, "Jet_puId_frac02"         )
    Jet_puId_frac03              [i] = book_float_branch(outTree, "Jet_puId_frac03"         )
    Jet_puId_frac04              [i] = book_float_branch(outTree, "Jet_puId_frac04"         )
    Jet_puId_ptD                 [i] = book_float_branch(outTree, "Jet_puId_ptD"            )
    Jet_puId_beta                [i] = book_float_branch(outTree, "Jet_puId_beta"           )
    Jet_puId_pull                [i] = book_float_branch(outTree, "Jet_puId_pull"           )
    Jet_puId_jetR                [i] = book_float_branch(outTree, "Jet_puId_jetR"           )
    Jet_puId_jetRchg             [i] = book_float_branch(outTree, "Jet_puId_jetRchg"        )
    Jet_nConstituents            [i] = book_int_branch  (outTree, "Jet_nConstituents"  )
    Jet_puId_nCharged            [i] = book_int_branch  (outTree, "Jet_puId_nCharged"       )
    Jet_dRMatch                  [i] = book_float_branch(outTree, "Jet_dRMatch"             )
    Jet_genJetIdx                [i] = book_int_branch  (outTree, "Jet_genJetIdx"           )
    Jet_partonFlavour            [i] = book_int_branch  (outTree, "Jet_partonFlavour"       )
    Zpt_Jetpt                    [i] = book_float_branch  (outTree, "Zpt_Jetpt"                    )
    dphi_zj                      [i] = book_float_branch  (outTree, "dphi_zj"                    )

for i, PFCandTree in enumerate(PFCandTrees):
    Jet_PFCand_dR                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_dR" )
    Jet_PFCand_puppiWeight                  [i] = book_float_branch  (PFCandTree, "Jet_PFCand_puppiWeight" )
    Jet_PFCand_puppiWeightNoLep             [i] = book_float_branch  (PFCandTree, "Jet_PFCand_puppiWeightNoLep" )
    Jet_PFCand_deta                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_deta" )
    Jet_PFCand_dphi                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_dphi" )
    Jet_PFCand_Absptdeta                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_Absptdeta" )
    Jet_PFCand_Absptdphi                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_Absptdphi" )
    Jet_PFCand_Frcptdeta                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_Frcptdeta" )
    Jet_PFCand_Frcptdphi                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_Frcptdphi" )
    Jet_PFCand_trkChi2                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_trkChi2" )
    Jet_PFCand_vtxChi2                           [i] = book_float_branch  (PFCandTree, "Jet_PFCand_vtxChi2" )
    

def dphi(phi1, phi2):
    dphi = phi1 - phi2
    Pi = ROOT.TMath.Pi()
    if dphi > Pi:
        dphi -= 2 * Pi
    elif dphi <= -Pi:
        dphi += 2 * Pi
    return dphi

pass_trig=0
pass_muon=0
pass_jet=0
pass_Zcut=0
nEvent=0
nPrompt=0
nPileup=0

ROOT.gSystem.Load("libFWCoreFWLite.so")
#AutoLibraryLoader.enable()
#ResJetPar=ROOT.JetCorrectorParameters("/u/user/shin/scratch/training_code/PUJetId_training_plotting/training/input/PUPPI_JEC/Winter22Run3_V1_MC_L2L3Residual_AK4PFPuppi.txt")
#L3JetPar=ROOT.JetCorrectorParameters("/u/user/shin/scratch/training_code/PUJetId_training_plotting/training/input/PUPPI_JEC/Winter22Run3_V1_MC_L3Absolute_AK4PFPuppi.txt")
#L2JetPar=ROOT.JetCorrectorParameters("/u/user/shin/scratch/training_code/PUJetId_training_plotting/training/input/PUPPI_JEC/Winter22Run3_V1_MC_L2Relative_AK4PFPuppi.txt")
#L1JetPar=ROOT.JetCorrectorParameters("/u/user/shin/scratch/training_code/PUJetId_training_plotting/training/input/PUPPI_JEC/Winter22Run3_V1_MC_L1FastJet_AK4PFPuppi.txt")

#vPar=[L1JetPar,L2JetPar,L3JetPar,ResJetPar]
#JetCorrector=ROOT.FactorizedJetCorrector(vPar)


#skip_entry=[]
#try:
#   for ievent, event in enumerate(tChain):
#      if ievent % 1000 == 0:
#         print("testing %s" % ievent)
#except RuntimeError:
   
#nEvents=tChain.GetEntries()
 

#nEvents = event.GetEntries()

#print("number of entries: %i" %nEvents)
for ievent, event in enumerate(tChain):
#for ievent in range(nEvents):
#    print(event)
#    print(type(event))
    if ievent % 1000 == 0:
        print("processing %s" % ievent)
#    if ievent == 1000: break
    nEvent+=1
#    event.GetEntry(ievent)    
    try:
        if event.nGenJet==0: 
            continue
    except AttributeError:
        print("Error in nGenJet")
        continue
############################################## add selection rules####################################################
####### first the event passes the trigger########3
#    if year=="2016":
#        trigger = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
#    if year =="2016APV":
#        trigger = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
#    if year =="2017" or year=="2018":
#        trigger = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8
#    if trigger != 1: continue


    nLeptons = event.nMuon
######require exactly 2 leptons and at least 1 jet to start ################
    if nLeptons != 2 or event.nJet<1: continue

    good_lep = ROOT.TLorentzVector()
    good_leps = 0
    good_jets = 0
    base_jet = 0
    good_jet=ROOT.TLorentzVector()
    good_jet_index=[]
    good_lep_index=[]
    clean_jet=True
    lep1=ROOT.TLorentzVector()
    lep2=ROOT.TLorentzVector()
    pass_trig += 1
    temp_lep=ROOT.TLorentzVector()
    for i in range(nLeptons):
############# muon pt cut, muon eta cut, tight muon ID, PfIsoId tight############ 
        if event.Muon_pt[i]>20 and abs(event.Muon_eta[i])<2.4 and event.Muon_tightId[i]==1 and ord(event.Muon_pfIsoId[i])>=4:
            good_lep.SetPtEtaPhiM(event.Muon_pt[i], event.Muon_eta[i],event.Muon_phi[i],event.Muon_mass[i])
            good_leps=good_leps+1
            good_lep_index.append(i)
####### at least two leptons that pass selection#######
    if good_leps<2: continue
    pass_muon += 1
###############select leading and subleading good muons to reconstruct Z candidate
    Z_cand=ROOT.TLorentzVector()
    lep1.SetPtEtaPhiM(event.Muon_pt[good_lep_index[0]], event.Muon_eta[good_lep_index[0]],event.Muon_phi[good_lep_index[0]],event.Muon_mass[good_lep_index[0]])
    lep2.SetPtEtaPhiM(event.Muon_pt[good_lep_index[1]], event.Muon_eta[good_lep_index[1]],event.Muon_phi[good_lep_index[1]],event.Muon_mass[good_lep_index[1]])
    Z_cand = lep1+lep2
    if Z_cand.M() < 70 or Z_cand.M() > 110: continue
    pass_Zcut += 1

#################################jet selection#########################    
    for j in range(event.nJet):
        if event.Jet_eta[j]>2.2 or event.Jet_eta[j]<1.4 or event.Jet_phi[j]>2.8 or event.Jet_phi[j]<1.2: continue 
        clean_jet = True
        good_jet = ROOT.TLorentzVector()

###############Jet Correction######################
#        JetCorrector.setJetEta()
#        JetCorrector.setJetPt()
#        JetCorrector.setJetA()
#        JetCorrector.setJetRho()



#        print("Jet number %i" %j)
#############tight jet ID, jet pt cut, jet eta cut ########################
        if event.Jet_jetId[j]!=2 and event.Jet_pt[j] <= 10 and abs(event.Jet_eta[j]) > 5: continue
        good_jet.SetPtEtaPhiM(event.Jet_pt[j], event.Jet_eta[j],event.Jet_phi[j],event.Jet_mass[j])
        base_jet += 1
#############Cross Clean Jet and make sure that muon is not near the jet###############
        for lep_index in range(nLeptons): 
            if event.Muon_looseId[lep_index]==0: continue
            temp_lep.SetPtEtaPhiM(event.Muon_pt[lep_index], event.Muon_eta[lep_index],event.Muon_phi[lep_index],event.Muon_mass[lep_index])
            if temp_lep.DeltaR(good_jet) < 0.4:
                clean_jet=False
        if clean_jet:
            good_jets=good_jets+1
            good_jet_index.append(j)
######at least one clean jet that pass the selection#####
    if good_jets<1: continue
    pass_jet += 1
################################################### Selection complete ########################################################
    matched_genjet_idx=[]
    Pileup_GenJet_idx=[]
#    for k in good_jet_index:
#        if event.Jet_genJetIdx[k] >= 0 and event.Jet_genJetIdx[k] < event.nGenJet:
#            matched_genjet_idx.append(event.Jet_genJetIdx[k])


#    for m in range(event.nGenJet):
#        if m in matched_genjet_idx: continue
#        Pileup_GenJet_idx.append(m)
#    for k in good_jet_index:
#        PU_dR=10000
#        if event.Jet_genJetIdx[k] > 0: continue
#        for m in Pileup_GenJet_idx:
#            temp_PU_dR=((event.Jet_eta[k]-event.GenJet_eta[m])**2+dphi(event.Jet_phi[k],event.GenJet_phi[m])**2)**0.5
            
    
        
    

    for k in good_jet_index:
        if k >= event.nJet: continue
        try:
            GenJetIdx= event.Jet_genJetIdx[k]
        except IndexError:
            print("Error in GenJetIdx")
            continue

        PFCand_puppiWeight_list = []
        PFCand_puppiWeightNoLep_list = []
        PFCand_dR_list = []        
        PFCand_Abspt_list = []        
        PFCand_Frcpt_list = []        
        PFCand_deta_list = []        
        PFCand_dphi_list = []        
        PFCand_Absptdeta_list = []        
        PFCand_Absptdphi_list = []        
        PFCand_Frcptdeta_list = []        
        PFCand_Frcptdphi_list = []        
        PFCand_trkChi2_list = []        
        PFCand_vtxChi2_list = []        

        pt_      = event.Jet_pt[k]
        eta_     = event.Jet_eta[k]
        phi_     = event.Jet_phi[k]
        flavor_  = event.Jet_partonFlavour[k]

        for m,p in enumerate(event.JetPFCand_jetIdx):
            if p == k:
                temp_PFCandIdx=event.JetPFCand_pfCandIdx[m]

                PFCand_puppiWeight_list.append(event.PFCand_puppiWeight[temp_PFCandIdx])
                PFCand_puppiWeightNoLep_list.append(event.PFCand_puppiWeightNoLep[temp_PFCandIdx])

                temp_PFCand_dR=-1
                temp_PFCand_dR = ((eta_-event.PFCand_eta[temp_PFCandIdx])**2+dphi(phi_,event.PFCand_phi[temp_PFCandIdx])**2)**0.5
                PFCand_dR_list.append(temp_PFCand_dR)
                temp_deta=eta_-event.PFCand_eta[temp_PFCandIdx]
                temp_dphi=dphi(phi_,event.PFCand_phi[temp_PFCandIdx])
                temp_Absptdeta=temp_deta*event.PFCand_pt[temp_PFCandIdx]
                temp_Absptdphi=temp_dphi*event.PFCand_pt[temp_PFCandIdx]
                temp_Frcptdeta=temp_Absptdeta/pt_
                temp_Frcptdphi=temp_Absptdphi/pt_


                PFCand_deta_list.append(temp_deta)
                PFCand_dphi_list.append(temp_dphi)
                PFCand_Absptdeta_list.append(temp_Absptdeta)
                PFCand_Absptdphi_list.append(temp_Absptdphi)
                PFCand_Frcptdeta_list.append(temp_Frcptdeta)
                PFCand_Frcptdphi_list.append(temp_Frcptdphi)
                PFCand_Abspt_list.append(event.PFCand_pt[temp_PFCandIdx])
                PFCand_Frcpt_list.append(event.PFCand_pt[temp_PFCandIdx]/pt_)

                PFCand_trkChi2_list.append(event.PFCand_trkChi2[temp_PFCandIdx])       
                PFCand_vtxChi2_list.append(event.PFCand_vtxChi2[temp_PFCandIdx])        



        dRMatch_ = 10
        temp_dphi_zj = dphi(phi_,Z_cand.Phi())
        if GenJetIdx < event.nGenJet and GenJetIdx >= 0:
            try:
                dRMatch_ = ((eta_-event.GenJet_eta[GenJetIdx])**2+dphi(phi_,event.GenJet_phi[GenJetIdx])**2)**0.5
            except IndexError:
                print("Error in dRMatch_")
                continue
#        if GenJetIdx < 0:
            






        #Jet pt upper and lower bound condition
#        if event.Jet_pt[k] > f_maxJetpt or event.Jet_pt[k] < f_minJetpt or Gen_index == -1: continue
        if event.Jet_pt[k] > f_maxJetpt or event.Jet_pt[k] < f_minJetpt: continue
        Z_jet_pt_ratio=Z_cand.Pt()/pt_
        isPrompt = False
        isPileup = False
                


        if (dRMatch_ <= 0.2 and GenJetIdx < event.nGenJet and GenJetIdx >= 0):
            isPrompt = True
            nPrompt += 1


        if (abs(flavor_) == 0 and GenJetIdx < 0):
#        if (GenJetIdx < 0):
            isPileup = True
            nPileup += 1
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
        fixedGridRhoFastjetAll             [key][0] = event.Rho_fixedGridRhoFastjetAll
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
#        Jet_dRMatch                   [key][0] = dRMatch_
        Jet_genJetIdx                 [key][0] = event.Jet_genJetIdx[k]
        Jet_partonFlavour             [key][0] = event.Jet_partonFlavour[k]
        Zpt_Jetpt                          [key][0] = Z_jet_pt_ratio
        dphi_zj                            [key][0] = temp_dphi_zj

        PFCandTree_toFill = PFCandTrees[key]
        for aa in range(len(PFCand_puppiWeight_list)):
            Jet_PFCand_dR                      [key][0] = PFCand_dR_list[aa]
            Jet_PFCand_deta                      [key][0] = PFCand_deta_list[aa]
            Jet_PFCand_dphi                      [key][0] = PFCand_dphi_list[aa]
            Jet_PFCand_puppiWeight             [key][0] = PFCand_puppiWeight_list[aa]
            Jet_PFCand_puppiWeightNoLep        [key][0] = PFCand_puppiWeightNoLep_list[aa]
            Jet_PFCand_Absptdeta                      [key][0] = PFCand_Absptdeta_list[aa]
            Jet_PFCand_Frcptdeta                      [key][0] = PFCand_Frcptdeta_list[aa]
            Jet_PFCand_Absptdphi                      [key][0] = PFCand_Absptdphi_list[aa]
            Jet_PFCand_Frcptdphi                      [key][0] = PFCand_Frcptdphi_list[aa]
            Jet_PFCand_trkChi2                      [key][0] = PFCand_trkChi2_list[aa]
            Jet_PFCand_vtxChi2                      [key][0] = PFCand_vtxChi2_list[aa]
            PFCandTree_toFill.Fill()
        outTree_toFill.Fill()

for outTree in outTrees:
    outTree.Write("", ROOT.TObject.kOverwrite)




print("trig cut:%f" %(pass_trig/nEvent))
print("muon cut:%f" %(pass_muon/pass_trig))
print("Z cut:%f" %(pass_Zcut/pass_muon))
print("jet cut:%f" %(pass_jet/pass_Zcut))
print("after selection:%f" %(pass_jet/nEvent))
print("number of Prompt Jet:%f" %nPrompt)
print("number of Pileup Jet:%f" %nPileup)
outFile.Write("", ROOT.TObject.kOverwrite)
outFile.Close()
