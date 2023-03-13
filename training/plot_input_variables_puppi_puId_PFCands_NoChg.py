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
for idx in range(27):
   inputfile="/d0/scratch/shin/training_code/PUJetId_training_plotting/training/input/puppi_2018_GenJetIdx_PFCands_input/training_trees_pt10To100_puppi_PFCands_106X_%i_2018.root" %idx
   chs_files.append(inputfile)
   inputfile="/d0/scratch/shin/training_code/PUJetId_training_plotting/training/input/puppi_2018_GenJetIdx_PFCands_input/training_trees_pt10To100_puppi_PFCands_106X_%i_2018.root" %idx
   puppi_files.append(inputfile)


for x,eta_bin in enumerate(eta_bins):
   c_nsig=0
   c_nbkg=0
   p_nsig=0
   p_nbkg=0
   ps_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   ps_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)

   cs_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   cs_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
 
   pb_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   pb_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)

   cb_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   cb_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
 
   c_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   c_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   c_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
 
   p_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   p_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   p_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)

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
         if eta_bin!="Eta3p0To5p0":
            Prompt_Charged=getattr(chs_sig_tree,"JetPuppi_PFCands_nCharged")
            Prompt_jetRchg=getattr(chs_sig_tree,"JetPuppi_PFCands_jetRchg")
            Prompt_beta=getattr(chs_sig_tree,"JetPuppi_PFCands_beta")
         if eta_bin!="Eta3p0To5p0" and Prompt_Charged==0:
            cs_beta.Fill(Prompt_beta)
            cs_jetRchg.Fill(Prompt_jetRchg)
            cs_Charged.Fill(Prompt_Charged)

         if eta_bin!="Eta3p0To5p0" and Prompt_Charged==0:
            c_beta.Fill(Prompt_beta)
            c_jetRchg.Fill(Prompt_jetRchg)
            c_Charged.Fill(Prompt_Charged)



      for entry in range(chs_bkg_tree.GetEntries()): 
         c_nbkg+=1
         if c_nbkg>250000: break
         chs_bkg_tree.GetEntry(entry)
         if eta_bin!="Eta3p0To5p0":
            Pileup_Charged=getattr(chs_bkg_tree,"JetPuppi_PFCands_nCharged")
            Pileup_jetRchg=getattr(chs_bkg_tree,"JetPuppi_PFCands_jetRchg")
            Pileup_beta=getattr(chs_bkg_tree,"JetPuppi_PFCands_beta")
         if eta_bin!="Eta3p0To5p0" and Pileup_Charged==0:
            cb_beta.Fill(Pileup_beta)
            cb_jetRchg.Fill(Pileup_jetRchg)
            cb_Charged.Fill(Pileup_Charged)
         if eta_bin!="Eta3p0To5p0" and Pileup_Charged==0:
            c_beta.Fill(Pileup_beta)
            c_jetRchg.Fill(Pileup_jetRchg)
            c_Charged.Fill(Pileup_Charged)

      for entry in range(puppi_sig_tree.GetEntries()): 
         p_nsig+=1
         if p_nsig>250000: break
         puppi_sig_tree.GetEntry(entry)
         if eta_bin!="Eta3p0To5p0":
            Prompt_Charged=getattr(puppi_sig_tree,"JetPuppi_puId_nCharged")
            Prompt_jetRchg=getattr(puppi_sig_tree,"JetPuppi_puId_jetRchg")
            Prompt_beta=getattr(puppi_sig_tree,"JetPuppi_puId_beta")
         if eta_bin!="Eta3p0To5p0"  and Prompt_Charged==0:
            ps_beta.Fill(Prompt_beta)
            ps_jetRchg.Fill(Prompt_jetRchg)
            ps_Charged.Fill(Prompt_Charged)
         if eta_bin!="Eta3p0To5p0" and Prompt_Charged==0:
            p_beta.Fill(Prompt_beta)
            p_jetRchg.Fill(Prompt_jetRchg)
            p_Charged.Fill(Prompt_Charged)
 

      for entry in range(puppi_bkg_tree.GetEntries()): 
         p_nbkg+=1
         if p_nbkg>250000: break
         puppi_bkg_tree.GetEntry(entry)
         if eta_bin!="Eta3p0To5p0":
            Pileup_Charged=getattr(puppi_bkg_tree,"JetPuppi_puId_nCharged")
            Pileup_jetRchg=getattr(puppi_bkg_tree,"JetPuppi_puId_jetRchg")
            Pileup_beta=getattr(puppi_bkg_tree,"JetPuppi_puId_beta")
         if eta_bin!="Eta3p0To5p0"  and Pileup_Charged==0:
            pb_beta.Fill(Pileup_beta)
            pb_jetRchg.Fill(Pileup_jetRchg)
            pb_Charged.Fill(Pileup_Charged)
         if eta_bin!="Eta3p0To5p0" and Pileup_Charged==0:
            p_beta.Fill(Pileup_beta)
            p_jetRchg.Fill(Pileup_jetRchg)
            p_Charged.Fill(Pileup_Charged)
      chs_file.Close()
      puppi_file.Close()





         






         





   if eta_bin!="Eta3p0To5p0":
      chs_signal=[cs_beta,cs_jetRchg,cs_Charged]
      chs_background=[cb_beta,cb_jetRchg,cb_Charged]
      chs_combine=[c_beta,c_jetRchg,c_Charged]
   if eta_bin!="Eta3p0To5p0":
      puppi_signal=[ps_beta,ps_jetRchg,ps_Charged]
      puppi_background=[pb_beta,pb_jetRchg,pb_Charged]
      puppi_combine=[p_beta,p_jetRchg,p_Charged]
   
   if eta_bin!="Eta3p0To5p0":
      var=["beta","jetRchg","Charged"]

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
      lgd.AddEntry(chs_signal[i],"PFCands")
      lgd.AddEntry(puppi_signal[i],"puId")
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
      c.SaveAs("puppi_PFCands_puId_plot_%s_%s_Prompt.png" %(name,eta_bin))

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
      lgd.AddEntry(chs_background[i],"PFCands")
      lgd.AddEntry(puppi_background[i],"puId")
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
      c.SaveAs("puppi_PFCands_puId_plot_%s_%s_Pileup.png" %(name,eta_bin))

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
      lgd.AddEntry(chs_combine[i],"PFCands")
      lgd.AddEntry(puppi_combine[i],"puId")
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
      c.SaveAs("puppi_PFCands_puId_plot_%s_%s_Combined.png" %(name,eta_bin))





