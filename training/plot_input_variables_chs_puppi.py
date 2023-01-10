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
   "--eta_s", type=str, default="Eta0p0To2p5", help="string of eta region Eta0p0To2p5, Eta2p5To2p75, Eta2p75To3p0, Eta3p0To5p0"
   )

parser.add_argument(
   "--eta_f", type=str, default="[0.0,2.5]", help="bracket range of eta region, [0.0,2.5], [2.5,2.75], [2.75,3.0], [3.0,5.0]"
   )
args = parser.parse_args()

eta_s=args.eta_s
eta_f=args.eta_f

eta_bins=[eta_s]
eta_bin_range=[eta_f]
chs_files=[]
puppi_files=[]
for idx in range(17):
   inputfile="/d0/scratch/shin/training_code/PUJetId_training_plotting/training/input/chs_2018_GenJetIdx/training_trees_pt10To100_chs_106X_%i_2018.root" %idx
   chs_files.append(inputfile)
   inputfile="/d0/scratch/shin/training_code/PUJetId_training_plotting/training/input/puppi_2018_GenJetIdx/training_trees_pt10To100_puppi_106X_%i_2018.root" %idx
   puppi_files.append(inputfile)


for x,eta_bin in enumerate(eta_bins):
   c_nsig=0
   c_nbkg=0
   p_nsig=0
   p_nbkg=0
   ps_npvs=ROOT.TH1D("npvs","npvsGood #eta #in %s" %eta_bin_range[x],70,0,70)
   ps_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_dR2Mean=ROOT.TH1D("dR2Mean","dR2Mean #eta #in %s" %eta_bin_range[x],30,0,0.1)
   ps_frac01=ROOT.TH1D("frac_01","frac_01 #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_frac02=ROOT.TH1D("frac_02","frac_02 #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_frac03=ROOT.TH1D("frac_03","frac_03 #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_frac04=ROOT.TH1D("frac_04","frac_04 #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_majW=ROOT.TH1D("majW","majW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   ps_minW=ROOT.TH1D("minW","minW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   ps_jetR=ROOT.TH1D("jetR","jetR #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_nConstituents=ROOT.TH1D("nConstituents","nConstituents #eta #in %s" %eta_bin_range[x],40,0,40)
   ps_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
   ps_ptD=ROOT.TH1D("ptD","ptD #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_pull=ROOT.TH1D("pull","pull #eta #in %s" %eta_bin_range[x],30,0,0.05)
   ps_Rho=ROOT.TH1D("Rho","Rho #eta #in %s" %eta_bin_range[x],70,0,70)

   cs_npvs=ROOT.TH1D("npvs","npvsGood #eta #in %s" %eta_bin_range[x],70,0,70)
   cs_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_dR2Mean=ROOT.TH1D("dR2Mean","dR2Mean #eta #in %s" %eta_bin_range[x],30,0,0.1)
   cs_frac01=ROOT.TH1D("frac_01","frac_01 #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_frac02=ROOT.TH1D("frac_02","frac_02 #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_frac03=ROOT.TH1D("frac_03","frac_03 #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_frac04=ROOT.TH1D("frac_04","frac_04 #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_majW=ROOT.TH1D("majW","majW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   cs_minW=ROOT.TH1D("minW","minW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   cs_jetR=ROOT.TH1D("jetR","jetR #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_nConstituents=ROOT.TH1D("nConstituents","nConstituents #eta #in %s" %eta_bin_range[x],40,0,40)
   cs_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
   cs_ptD=ROOT.TH1D("ptD","ptD #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_pull=ROOT.TH1D("pull","pull #eta #in %s" %eta_bin_range[x],30,0,0.05)
   cs_Rho=ROOT.TH1D("Rho","Rho #eta #in %s" %eta_bin_range[x],70,0,70)
 
   pb_npvs=ROOT.TH1D("npvs","npvsGood #eta #in %s" %eta_bin_range[x],70,0,70)
   pb_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_dR2Mean=ROOT.TH1D("dR2Mean","dR2Mean #eta #in %s" %eta_bin_range[x],30,0,0.1)
   pb_frac01=ROOT.TH1D("frac_01","frac_01 #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_frac02=ROOT.TH1D("frac_02","frac_02 #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_frac03=ROOT.TH1D("frac_03","frac_03 #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_frac04=ROOT.TH1D("frac_04","frac_04 #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_majW=ROOT.TH1D("majW","majW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   pb_minW=ROOT.TH1D("minW","minW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   pb_jetR=ROOT.TH1D("jetR","jetR #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_nConstituents=ROOT.TH1D("nConstituents","nConstituents #eta #in %s" %eta_bin_range[x],40,0,40)
   pb_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
   pb_ptD=ROOT.TH1D("ptD","ptD #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_pull=ROOT.TH1D("pull","pull #eta #in %s" %eta_bin_range[x],30,0,0.05)
   pb_Rho=ROOT.TH1D("Rho","Rho #eta #in %s" %eta_bin_range[x],70,0,70)

   cb_npvs=ROOT.TH1D("npvs","npvsGood #eta #in %s" %eta_bin_range[x],70,0,70)
   cb_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_dR2Mean=ROOT.TH1D("dR2Mean","dR2Mean #eta #in %s" %eta_bin_range[x],30,0,0.1)
   cb_frac01=ROOT.TH1D("frac_01","frac_01 #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_frac02=ROOT.TH1D("frac_02","frac_02 #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_frac03=ROOT.TH1D("frac_03","frac_03 #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_frac04=ROOT.TH1D("frac_04","frac_04 #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_majW=ROOT.TH1D("majW","majW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   cb_minW=ROOT.TH1D("minW","minW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   cb_jetR=ROOT.TH1D("jetR","jetR #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_nConstituents=ROOT.TH1D("nConstituents","nConstituents #eta #in %s" %eta_bin_range[x],40,0,40)
   cb_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
   cb_ptD=ROOT.TH1D("ptD","ptD #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_pull=ROOT.TH1D("pull","pull #eta #in %s" %eta_bin_range[x],30,0,0.05)
   cb_Rho=ROOT.TH1D("Rho","Rho #eta #in %s" %eta_bin_range[x],70,0,70)
 

   for idx in range(17):
      chs_file=ROOT.TFile.Open(chs_files[idx],"READ")
      puppi_file=ROOT.TFile.Open(puppi_files[idx],"READ")
      chs_sig_tree=chs_file.Get("%s_Prompt" %eta_bin) 
      chs_bkg_tree=chs_file.Get("%s_Pileup" %eta_bin)
      puppi_sig_tree=puppi_file.Get("%s_Prompt" %eta_bin) 
      puppi_bkg_tree=puppi_file.Get("%s_Pileup" %eta_bin)
      for entry in range(chs_sig_tree.GetEntries()): 
         c_nsig+=1
         if c_nsig>250000: break
         chs_sig_tree.GetEntry(entry)
         Prompt_npvs=getattr(chs_sig_tree,"PV_npvsGood")
         Prompt_dR2Mean=getattr(chs_sig_tree,"Jet_puId_dR2Mean")
         Prompt_frac01=getattr(chs_sig_tree,"Jet_puId_frac01")
         Prompt_frac02=getattr(chs_sig_tree,"Jet_puId_frac02")
         Prompt_frac03=getattr(chs_sig_tree,"Jet_puId_frac03")
         Prompt_frac04=getattr(chs_sig_tree,"Jet_puId_frac04")
         Prompt_majW=getattr(chs_sig_tree,"Jet_puId_majW")
         Prompt_minW=getattr(chs_sig_tree,"Jet_puId_minW")
         Prompt_jetR=getattr(chs_sig_tree,"Jet_puId_jetR")
         Prompt_nConstituents=getattr(chs_sig_tree,"Jet_nConstituents")
         Prompt_ptD=getattr(chs_sig_tree,"Jet_puId_ptD")
         Prompt_pull=getattr(chs_sig_tree,"Jet_puId_pull")
         Prompt_Rho=getattr(chs_sig_tree,"fixedGridRhoFastjetAll")
         if eta_bin!="Eta3p0To5p0":
            Prompt_Charged=getattr(chs_sig_tree,"Jet_puId_nCharged")
            Prompt_jetRchg=getattr(chs_sig_tree,"Jet_puId_jetRchg")
            Prompt_beta=getattr(chs_sig_tree,"Jet_puId_beta")
         cs_npvs.Fill(Prompt_npvs)
         cs_dR2Mean.Fill(Prompt_dR2Mean)
         cs_frac01.Fill(Prompt_frac01)
         cs_frac02.Fill(Prompt_frac02)
         cs_frac03.Fill(Prompt_frac03)
         cs_frac04.Fill(Prompt_frac04)
         cs_majW.Fill(Prompt_majW)
         cs_minW.Fill(Prompt_minW)
         cs_jetR.Fill(Prompt_jetR)
         cs_nConstituents.Fill(Prompt_nConstituents)
         cs_ptD.Fill(Prompt_ptD)
         cs_pull.Fill(Prompt_pull)
         cs_Rho.Fill(Prompt_Rho)
         if eta_bin!="Eta3p0To5p0":
            cs_beta.Fill(Prompt_beta)
            cs_jetRchg.Fill(Prompt_jetRchg)
            cs_Charged.Fill(Prompt_Charged)

      for entry in range(chs_bkg_tree.GetEntries()): 
         c_nbkg+=1
         if c_nbkg>250000: break
         chs_bkg_tree.GetEntry(entry)
         Pileup_npvs=getattr(chs_bkg_tree,"PV_npvsGood")
         Pileup_dR2Mean=getattr(chs_bkg_tree,"Jet_puId_dR2Mean")
         Pileup_frac01=getattr(chs_bkg_tree,"Jet_puId_frac01")
         Pileup_frac02=getattr(chs_bkg_tree,"Jet_puId_frac02")
         Pileup_frac03=getattr(chs_bkg_tree,"Jet_puId_frac03")
         Pileup_frac04=getattr(chs_bkg_tree,"Jet_puId_frac04")
         Pileup_majW=getattr(chs_bkg_tree,"Jet_puId_majW")
         Pileup_minW=getattr(chs_bkg_tree,"Jet_puId_minW")
         Pileup_jetR=getattr(chs_bkg_tree,"Jet_puId_jetR")
         Pileup_nConstituents=getattr(chs_bkg_tree,"Jet_nConstituents")
         Pileup_ptD=getattr(chs_bkg_tree,"Jet_puId_ptD")
         Pileup_pull=getattr(chs_bkg_tree,"Jet_puId_pull")
         Pileup_Rho=getattr(chs_bkg_tree,"fixedGridRhoFastjetAll")
         if eta_bin!="Eta3p0To5p0":
            Pileup_Charged=getattr(chs_bkg_tree,"Jet_puId_nCharged")
            Pileup_jetRchg=getattr(chs_bkg_tree,"Jet_puId_jetRchg")
            Pileup_beta=getattr(chs_bkg_tree,"Jet_puId_beta")
         cb_npvs.Fill(Pileup_npvs)
         cb_dR2Mean.Fill(Pileup_dR2Mean)
         cb_frac01.Fill(Pileup_frac01)
         cb_frac02.Fill(Pileup_frac02)
         cb_frac03.Fill(Pileup_frac03)
         cb_frac04.Fill(Pileup_frac04)
         cb_majW.Fill(Pileup_majW)
         cb_minW.Fill(Pileup_minW)
         cb_jetR.Fill(Pileup_jetR)
         cb_nConstituents.Fill(Pileup_nConstituents)
         cb_ptD.Fill(Pileup_ptD)
         cb_pull.Fill(Pileup_pull)
         cb_Rho.Fill(Pileup_Rho)
         if eta_bin!="Eta3p0To5p0":
            cb_beta.Fill(Pileup_beta)
            cb_jetRchg.Fill(Pileup_jetRchg)
            cb_Charged.Fill(Pileup_Charged)

      for entry in range(puppi_sig_tree.GetEntries()): 
         p_nsig+=1
         if p_nsig>250000: break
         puppi_sig_tree.GetEntry(entry)
         Prompt_npvs=getattr(puppi_sig_tree,"PV_npvsGood")
         Prompt_dR2Mean=getattr(puppi_sig_tree,"JetPuppi_puId_dR2Mean")
         Prompt_frac01=getattr(puppi_sig_tree,"JetPuppi_puId_frac01")
         Prompt_frac02=getattr(puppi_sig_tree,"JetPuppi_puId_frac02")
         Prompt_frac03=getattr(puppi_sig_tree,"JetPuppi_puId_frac03")
         Prompt_frac04=getattr(puppi_sig_tree,"JetPuppi_puId_frac04")
         Prompt_majW=getattr(puppi_sig_tree,"JetPuppi_puId_majW")
         Prompt_minW=getattr(puppi_sig_tree,"JetPuppi_puId_minW")
         Prompt_jetR=getattr(puppi_sig_tree,"JetPuppi_puId_jetR")
         Prompt_nConstituents=getattr(puppi_sig_tree,"JetPuppi_nConstituents")
         Prompt_ptD=getattr(puppi_sig_tree,"JetPuppi_puId_ptD")
         Prompt_pull=getattr(puppi_sig_tree,"JetPuppi_puId_pull")
         Prompt_Rho=getattr(puppi_sig_tree,"fixedGridRhoFastjetAll")
         if eta_bin!="Eta3p0To5p0":
            Prompt_Charged=getattr(puppi_sig_tree,"JetPuppi_puId_nCharged")
            Prompt_jetRchg=getattr(puppi_sig_tree,"JetPuppi_puId_jetRchg")
            Prompt_beta=getattr(puppi_sig_tree,"JetPuppi_puId_beta")
         ps_npvs.Fill(Prompt_npvs)
         ps_dR2Mean.Fill(Prompt_dR2Mean)
         ps_frac01.Fill(Prompt_frac01)
         ps_frac02.Fill(Prompt_frac02)
         ps_frac03.Fill(Prompt_frac03)
         ps_frac04.Fill(Prompt_frac04)
         ps_majW.Fill(Prompt_majW)
         ps_minW.Fill(Prompt_minW)
         ps_jetR.Fill(Prompt_jetR)
         ps_nConstituents.Fill(Prompt_nConstituents)
         ps_ptD.Fill(Prompt_ptD)
         ps_pull.Fill(Prompt_pull)
         ps_Rho.Fill(Prompt_Rho)
         if eta_bin!="Eta3p0To5p0":
            ps_beta.Fill(Prompt_beta)
            ps_jetRchg.Fill(Prompt_jetRchg)
            ps_Charged.Fill(Prompt_Charged)
 

      for entry in range(puppi_bkg_tree.GetEntries()): 
         p_nbkg+=1
         if p_nbkg>250000: break
         puppi_bkg_tree.GetEntry(entry)
         Pileup_npvs=getattr(puppi_bkg_tree,"PV_npvsGood")
         Pileup_dR2Mean=getattr(puppi_bkg_tree,"JetPuppi_puId_dR2Mean")
         Pileup_frac01=getattr(puppi_bkg_tree,"JetPuppi_puId_frac01")
         Pileup_frac02=getattr(puppi_bkg_tree,"JetPuppi_puId_frac02")
         Pileup_frac03=getattr(puppi_bkg_tree,"JetPuppi_puId_frac03")
         Pileup_frac04=getattr(puppi_bkg_tree,"JetPuppi_puId_frac04")
         Pileup_majW=getattr(puppi_bkg_tree,"JetPuppi_puId_majW")
         Pileup_minW=getattr(puppi_bkg_tree,"JetPuppi_puId_minW")
         Pileup_jetR=getattr(puppi_bkg_tree,"JetPuppi_puId_jetR")
         Pileup_nConstituents=getattr(puppi_bkg_tree,"JetPuppi_nConstituents")
         Pileup_ptD=getattr(puppi_bkg_tree,"JetPuppi_puId_ptD")
         Pileup_pull=getattr(puppi_bkg_tree,"JetPuppi_puId_pull")
         Pileup_Rho=getattr(puppi_bkg_tree,"fixedGridRhoFastjetAll")
         if eta_bin!="Eta3p0To5p0":
            Pileup_Charged=getattr(puppi_bkg_tree,"JetPuppi_puId_nCharged")
            Pileup_jetRchg=getattr(puppi_bkg_tree,"JetPuppi_puId_jetRchg")
            Pileup_beta=getattr(puppi_bkg_tree,"JetPuppi_puId_beta")
         pb_npvs.Fill(Pileup_npvs)
         pb_dR2Mean.Fill(Pileup_dR2Mean)
         pb_frac01.Fill(Pileup_frac01)
         pb_frac02.Fill(Pileup_frac02)
         pb_frac03.Fill(Pileup_frac03)
         pb_frac04.Fill(Pileup_frac04)
         pb_majW.Fill(Pileup_majW)
         pb_minW.Fill(Pileup_minW)
         pb_jetR.Fill(Pileup_jetR)
         pb_nConstituents.Fill(Pileup_nConstituents)
         pb_ptD.Fill(Pileup_ptD)
         pb_pull.Fill(Pileup_pull)
         pb_Rho.Fill(Pileup_Rho)
         if eta_bin!="Eta3p0To5p0":
            pb_beta.Fill(Pileup_beta)
            pb_jetRchg.Fill(Pileup_jetRchg)
            pb_Charged.Fill(Pileup_Charged)
      chs_file.Close()
      puppi_file.Close()





         






         


   c_npvs=cs_npvs
   c_dR2Mean=cs_dR2Mean
   c_frac01=cs_frac01
   c_frac02=cs_frac02
   c_frac03=cs_frac03
   c_frac04=cs_frac04
   c_majW=cs_majW
   c_minW=cs_minW
   c_jetR=cs_jetR
   c_nConstituents=cs_nConstituents
   c_ptD=cs_ptD
   c_pull=cs_pull
   c_Rho=cs_Rho
   c_beta=cs_beta
   c_jetRchg=cs_jetRchg
   c_Charged=cs_Charged

   c_npvs.Add(cb_npvs)
   c_dR2Mean.Add(cb_dR2Mean)
   c_frac01.Add(cb_frac01)
   c_frac02.Add(cb_frac02)
   c_frac03.Add(cb_frac03)
   c_frac04.Add(cb_frac04)
   c_majW.Add(cb_majW)
   c_minW.Add(cb_minW)
   c_jetR.Add(cb_jetR)
   c_nConstituents.Add(cb_nConstituents)
   c_ptD.Add(cb_ptD)
   c_pull.Add(cb_pull)
   c_Rho.Add(cb_Rho)
   c_beta.Add(cb_beta)
   c_jetRchg.Add(cb_jetRchg)
   c_Charged.Add(cb_Charged)


   p_npvs=ps_npvs
   p_dR2Mean=ps_dR2Mean
   p_frac01=ps_frac01
   p_frac02=ps_frac02
   p_frac03=ps_frac03
   p_frac04=ps_frac04
   p_majW=ps_majW
   p_minW=ps_minW
   p_jetR=ps_jetR
   p_nConstituents=ps_nConstituents
   p_ptD=ps_ptD
   p_pull=ps_pull
   p_Rho=ps_Rho
   p_beta=ps_beta
   p_jetRchg=ps_jetRchg
   p_Charged=ps_Charged

   p_npvs.Add(pb_npvs)
   p_dR2Mean.Add(pb_dR2Mean)
   p_frac01.Add(pb_frac01)
   p_frac02.Add(pb_frac02)
   p_frac03.Add(pb_frac03)
   p_frac04.Add(pb_frac04)
   p_majW.Add(pb_majW)
   p_minW.Add(pb_minW)
   p_jetR.Add(pb_jetR)
   p_nConstituents.Add(pb_nConstituents)
   p_ptD.Add(pb_ptD)
   p_pull.Add(pb_pull)
   p_Rho.Add(pb_Rho)
   p_beta.Add(pb_beta)
   p_jetRchg.Add(pb_jetRchg)
   p_Charged.Add(pb_Charged)



   if eta_bin!="Eta3p0To5p0":
      chs_signal=[cs_npvs,cs_dR2Mean,cs_frac01,cs_frac02,cs_frac03,cs_frac04,cs_majW,cs_minW,cs_jetR,cs_nConstituents,cs_ptD,cs_pull,cs_Rho,cs_beta,cs_jetRchg,cs_Charged]
      chs_background=[cb_npvs,cb_dR2Mean,cb_frac01,cb_frac02,cb_frac03,cb_frac04,cb_majW,cb_minW,cb_jetR,cb_nConstituents,cb_ptD,cb_pull,cb_Rho,cb_beta,cb_jetRchg,cb_Charged]
      chs_combine=[c_npvs,c_dR2Mean,c_frac01,c_frac02,c_frac03,c_frac04,c_majW,c_minW,c_jetR,c_nConstituents,c_ptD,c_pull,c_Rho,c_beta,c_jetRchg,c_Charged]
   if eta_bin=="Eta3p0To5p0":
      chs_signal=[cs_npvs,cs_dR2Mean,cs_frac01,cs_frac02,cs_frac03,cs_frac04,cs_majW,cs_minW,cs_jetR,cs_nConstituents,cs_ptD,cs_pull,cs_Rho]
      chs_background=[cb_npvs,cb_dR2Mean,cb_frac01,cb_frac02,cb_frac03,cb_frac04,cb_majW,cb_minW,cb_jetR,cb_nConstituents,cb_ptD,cb_pull,cb_Rho]
      chs_combine=[c_npvs,c_dR2Mean,c_frac01,c_frac02,c_frac03,c_frac04,c_majW,c_minW,c_jetR,c_nConstituents,c_ptD,c_pull,c_Rho]
   if eta_bin!="Eta3p0To5p0":
      puppi_signal=[ps_npvs,ps_dR2Mean,ps_frac01,ps_frac02,ps_frac03,ps_frac04,ps_majW,ps_minW,ps_jetR,ps_nConstituents,ps_ptD,ps_pull,ps_Rho,ps_beta,ps_jetRchg,ps_Charged]
      puppi_background=[pb_npvs,pb_dR2Mean,pb_frac01,pb_frac02,pb_frac03,pb_frac04,pb_majW,pb_minW,pb_jetR,pb_nConstituents,pb_ptD,pb_pull,pb_Rho,pb_beta,pb_jetRchg,pb_Charged]
      puppi_combine=[p_npvs,p_dR2Mean,p_frac01,p_frac02,p_frac03,p_frac04,p_majW,p_minW,p_jetR,p_nConstituents,p_ptD,p_pull,p_Rho,p_beta,p_jetRchg,p_Charged]
   if eta_bin=="Eta3p0To5p0":
      puppi_signal=[ps_npvs,ps_dR2Mean,ps_frac01,ps_frac02,ps_frac03,ps_frac04,ps_majW,ps_minW,ps_jetR,ps_nConstituents,ps_ptD,ps_pull,ps_Rho]
      puppi_background=[pb_npvs,pb_dR2Mean,pb_frac01,pb_frac02,pb_frac03,pb_frac04,pb_majW,pb_minW,pb_jetR,pb_nConstituents,pb_ptD,pb_pull,pb_Rho]
      puppi_combine=[p_npvs,p_dR2Mean,p_frac01,p_frac02,p_frac03,p_frac04,p_majW,p_minW,p_jetR,p_nConstituents,p_ptD,p_pull,p_Rho]
   
   if eta_bin!="Eta3p0To5p0":
      var=["npvs","dR2Mean","frac_01","frac_02","frac_03","frac_04","majW","minW","jetR","nConstitunents","ptD","pull","Rho","beta","jetRchg","Charged"]
   if eta_bin=="Eta3p0To5p0":
      var=["npvs","dR2Mean","frac_01","frac_02","frac_03","frac_04","majW","minW","jetR","nConstitunents","ptD","pull","Rho"]

   for i,name in enumerate(var):
      ROOT.gStyle.SetHatchesLineWidth(2)
      chs_signal[i].SetLineColor(ROOT.kBlue)
      chs_signal[i].SetFillColorAlpha(ROOT.kBlue-7,0.35)
      chs_signal[i].SetStats(0)
      chs_signal[i].SetFillStyle(1001)
      puppi_signal[i].SetLineColor(ROOT.kRed)
      puppi_signal[i].SetFillColor(ROOT.kRed+1)
      puppi_signal[i].SetStats(0)
      puppi_signal[i].SetLineWidth(2)
      puppi_signal[i].SetFillStyle(3554)

      chs_signal[i].GetXaxis().SetTitle(name)
      puppi_signal[i].GetXaxis().SetTitle(name)
      chs_signal[i].GetYaxis().SetTitle("Events")
      puppi_signal[i].GetYaxis().SetTitle("Events")


      sig_max_bin=chs_signal[i].GetMaximumBin()
      sig_max=chs_signal[i].GetBinContent(sig_max_bin)
      bkg_max_bin=puppi_signal[i].GetMaximumBin()
      bkg_max=puppi_signal[i].GetBinContent(bkg_max_bin)

      lgd = ROOT.TLegend(0.63, 0.9, 0.88, 0.75)
      lgd.SetBorderSize(1)
      lgd.SetFillColor(0)
      lgd.SetFillStyle(0)
      lgd.SetTextFont(42)
      lgd.SetTextSize(0.030)
      lgd.AddEntry(chs_signal[i],"CHS")
      lgd.AddEntry(puppi_signal[i],"PUPPI")
      c=ROOT.TCanvas("canvas","canvas",500,500)
      c.cd()
      if sig_max>bkg_max:
         chs_signal[i].Draw("HIST")
         puppi_signal[i].Draw("HIST SAME")
      if sig_max <= bkg_max:
         puppi_signal[i].Draw("HIST")
         chs_signal[i].Draw("HIST SAME")
         puppi_signal[i].Draw("HIST SAME")
      lgd.Draw("SAME")
      c.SaveAs("GenIdx_plot_%s_%s_Prompt.png" %(name,eta_bin))

   for i,name in enumerate(var):
      ROOT.gStyle.SetHatchesLineWidth(2)
      chs_background[i].SetLineColor(ROOT.kBlue)
      chs_background[i].SetFillColorAlpha(ROOT.kBlue-7,0.35)
      chs_background[i].SetStats(0)
      chs_background[i].SetFillStyle(1001)
      puppi_background[i].SetLineColor(ROOT.kRed)
      puppi_background[i].SetFillColor(ROOT.kRed+1)
      puppi_background[i].SetStats(0)
      puppi_background[i].SetLineWidth(2)
      puppi_background[i].SetFillStyle(3554)

      chs_background[i].GetXaxis().SetTitle(name)
      puppi_background[i].GetXaxis().SetTitle(name)
      chs_background[i].GetYaxis().SetTitle("Events")
      puppi_background[i].GetYaxis().SetTitle("Events")


      sig_max_bin=chs_background[i].GetMaximumBin()
      sig_max=chs_background[i].GetBinContent(sig_max_bin)
      bkg_max_bin=puppi_background[i].GetMaximumBin()
      bkg_max=puppi_background[i].GetBinContent(bkg_max_bin)

      lgd = ROOT.TLegend(0.63, 0.9, 0.88, 0.75)
      lgd.SetBorderSize(1)
      lgd.SetFillColor(0)
      lgd.SetFillStyle(0)
      lgd.SetTextFont(42)
      lgd.SetTextSize(0.030)
      lgd.AddEntry(chs_background[i],"CHS")
      lgd.AddEntry(puppi_background[i],"PUPPI")
      c=ROOT.TCanvas("canvas","canvas",500,500)
      c.cd()
      if sig_max>bkg_max:
         chs_background[i].Draw("HIST")
         puppi_background[i].Draw("HIST SAME")
      if sig_max <= bkg_max:
         puppi_background[i].Draw("HIST")
         chs_background[i].Draw("HIST SAME")
         puppi_background[i].Draw("HIST SAME")
      lgd.Draw("SAME")
      c.SaveAs("GenIdx_plot_%s_%s_Pileup.png" %(name,eta_bin))

   for i,name in enumerate(var):
      ROOT.gStyle.SetHatchesLineWidth(2)
      chs_combine[i].SetLineColor(ROOT.kBlue)
      chs_combine[i].SetFillColorAlpha(ROOT.kBlue-7,0.35)
      chs_combine[i].SetStats(0)
      chs_combine[i].SetFillStyle(1001)
      puppi_combine[i].SetLineColor(ROOT.kRed)
      puppi_combine[i].SetFillColor(ROOT.kRed+1)
      puppi_combine[i].SetStats(0)
      puppi_combine[i].SetLineWidth(2)
      puppi_combine[i].SetFillStyle(3554)

      chs_combine[i].GetXaxis().SetTitle(name)
      puppi_combine[i].GetXaxis().SetTitle(name)
      chs_combine[i].GetYaxis().SetTitle("Events")
      puppi_combine[i].GetYaxis().SetTitle("Events")


      sig_max_bin=chs_combine[i].GetMaximumBin()
      sig_max=chs_combine[i].GetBinContent(sig_max_bin)
      bkg_max_bin=puppi_combine[i].GetMaximumBin()
      bkg_max=puppi_combine[i].GetBinContent(bkg_max_bin)

      lgd = ROOT.TLegend(0.63, 0.9, 0.88, 0.75)
      lgd.SetBorderSize(1)
      lgd.SetFillColor(0)
      lgd.SetFillStyle(0)
      lgd.SetTextFont(42)
      lgd.SetTextSize(0.030)
      lgd.AddEntry(chs_combine[i],"CHS")
      lgd.AddEntry(puppi_combine[i],"PUPPI")
      c=ROOT.TCanvas("canvas","canvas",500,500)
      c.cd()
      if sig_max>bkg_max:
         chs_combine[i].Draw("HIST")
         puppi_combine[i].Draw("HIST SAME")
      if sig_max <= bkg_max:
         puppi_combine[i].Draw("HIST")
         chs_combine[i].Draw("HIST SAME")
         puppi_combine[i].Draw("HIST SAME")
      lgd.Draw("SAME")
      c.SaveAs("GenIdx_plot_%s_%s_Combined.png" %(name,eta_bin))





