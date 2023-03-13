#!/usr/bin/env python3


import ROOT
import os
import argparse

ROOT.gROOT.SetBatch(True)

if "HOME" not in os.environ:
   os.environ["HOME"] = "" # needed for condor batch mode

parser = argparse.ArgumentParser(
   description="Validation plot"
   )
parser.add_argument(
   "--jet_type", type=str, default="PUPPI", help=""
   )
parser.add_argument(
   "--process_name", type=str, default="QCD", help=""
   )
parser.add_argument(
   "--eta_s", type=str, default="Eta0p0To2p5", help="string of eta region Eta0p0To2p5, Eta2p5To2p75, Eta2p75To3p0, Eta3p0To5p0"
   )
parser.add_argument(
   "--eta_min", type=float, default=0.0, help="min eta"
   )
parser.add_argument(
   "--eta_max", type=float, default=2.5, help="max eta"
   )

args = parser.parse_args()

jet_type=args.jet_type
process_name = args.process_name
eta_min=args.eta_min
eta_max=args.eta_max
eta_s=args.eta_s

b_npvs=ROOT.TH1D("npvs","npvsGood" ,70,0,70)
b_beta=ROOT.TH1D("beta","beta", 30,0,0.1)
b_dR2Mean=ROOT.TH1D("dR2Mean","dR2Mean",30,0,0.1)
b_frac01=ROOT.TH1D("frac_01","frac_01",30,0,1)
b_frac02=ROOT.TH1D("frac_02","frac_02",30,0,1)
b_frac03=ROOT.TH1D("frac_03","frac_03",30,0,1)
b_frac04=ROOT.TH1D("frac_04","frac_04",30,0,1)
b_majW=ROOT.TH1D("majW","majW",30,0,0.4)
b_minW=ROOT.TH1D("minW","minW",30,0,0.4)
b_jetR=ROOT.TH1D("jetR","jetR",30,0,1)
b_jetRchg=ROOT.TH1D("jetRchg","jetRchg",30,0,1)
#b_nConstituents=ROOT.TH1D("nConstituents","nConstituents",40,0,40)
b_nCharged=ROOT.TH1D("nCharged","nCharged",30,0,30)
b_ptD=ROOT.TH1D("ptD","ptD",30,0,1)
b_pull=ROOT.TH1D("pull","pull",30,0,0.05)


a_npvs=ROOT.TH1D("npvs","npvsGood" ,70,0,70)
a_beta=ROOT.TH1D("beta","beta", 30,0,0.1)
a_dR2Mean=ROOT.TH1D("dR2Mean","dR2Mean",30,0,0.1)
a_frac01=ROOT.TH1D("frac_01","frac_01",30,0,1)
a_frac02=ROOT.TH1D("frac_02","frac_02",30,0,1)
a_frac03=ROOT.TH1D("frac_03","frac_03",30,0,1)
a_frac04=ROOT.TH1D("frac_04","frac_04",30,0,1)
a_majW=ROOT.TH1D("majW","majW",30,0,0.4)
a_minW=ROOT.TH1D("minW","minW",30,0,0.4)
a_jetR=ROOT.TH1D("jetR","jetR",30,0,1)
a_jetRchg=ROOT.TH1D("jetRchg","jetRchg",30,0,1)
#a_nConstituents=ROOT.TH1D("nConstituents","nConstituents",40,0,40)
a_nCharged=ROOT.TH1D("nCharged","nCharged",30,0,30)
a_ptD=ROOT.TH1D("ptD","ptD",30,0,1)
a_pull=ROOT.TH1D("pull","pull",30,0,0.05)


#before_tree = before_file.Get("Events")
#after_tree = after_file.Get("Events")

before_Chain = ROOT.TChain("Events")
after_Chain = ROOT.TChain("Events")

before_Chain.Add("/u/user/shin/Validataion_file/Before_fix_tree_jmepfnano_%s.root" %process_name.lower())
after_Chain.Add("/u/user/shin/Validataion_file/After_fix_tree_jmepfnano_%s.root" %process_name.lower())

for event in before_Chain:
  Before_npvs= event.PV_npvsGood
  if (jet_type == "PUPPI"):
    for k in range(event.nJet):
      if event.Jet_pt[k]>20 and event.Jet_eta[k]>eta_min and event.Jet_eta[k]<=eta_max:
        Before_dR2Mean=event.Jet_puId_dR2Mean[k]
        Before_frac01=event.Jet_puId_frac01[k]
        Before_frac02=event.Jet_puId_frac02[k]
        Before_frac03=event.Jet_puId_frac03[k]
        Before_frac04=event.Jet_puId_frac04[k]
        Before_majW=event.Jet_puId_majW[k]
        Before_minW=event.Jet_puId_minW[k]
        Before_jetR=event.Jet_puId_jetR[k]
#      if ("0x" in event.Jet_nConstituents[k].encode('utf-8').strip()):
#         Before_nConstituents=int(event.Jet_nConstituents[k],0)
#      else:
#         Before_nConstituents=ord(event.Jet_nConstituents[k])
        Before_ptD=event.Jet_puId_ptD[k]
        Before_pull=event.Jet_puId_pull[k]
        Before_nCharged=event.Jet_puId_nCharged[k]
        Before_jetRchg=event.Jet_puId_jetRchg[k]
        Before_beta=event.Jet_puId_beta[k]

        b_dR2Mean.Fill(Before_dR2Mean)
        b_frac01.Fill(Before_frac02)
        b_frac02.Fill(Before_frac02)
        b_dR2Mean.Fill(Before_dR2Mean)
        b_frac01.Fill(Before_frac02)
        b_frac02.Fill(Before_frac02)
        b_frac03.Fill(Before_frac03)
        b_frac04.Fill(Before_frac04)
        b_majW.Fill(Before_majW)
        b_minW.Fill(Before_minW)
        b_jetR.Fill(Before_jetR)
 #    b_nConstituents.Fill(Before_nConstituents,0)
        b_ptD.Fill(Before_ptD)
        b_pull.Fill(Before_pull)
        b_nCharged.Fill(Before_nCharged)
        b_jetRchg.Fill(Before_jetRchg)
        b_beta.Fill(Before_beta)
  if (jet_type == "CHS"):
    for k in range(event.nJetCHS):
      if event.JetCHS_pt[k]>20 and event.JetCHS_eta[k]>eta_min and event.JetCHS_eta[k]<=eta_max:
        Before_dR2Mean=event.JetCHS_puId_dR2Mean[k]
        Before_frac01=event.JetCHS_puId_frac01[k]
        Before_frac02=event.JetCHS_puId_frac02[k]
        Before_frac03=event.JetCHS_puId_frac03[k]
        Before_frac04=event.JetCHS_puId_frac04[k]
        Before_majW=event.JetCHS_puId_majW[k]
        Before_minW=event.JetCHS_puId_minW[k]
        Before_jetR=event.JetCHS_puId_jetR[k]
#       if ("0x" in event.Jet_nConstituents[k]):
#         Before_nConstituents=int(event.JetCHS_nConstituents[k],0)
#       else:
#          Before_nConstituents=ord(event.JetCHS_nConstituents[k])
        Before_ptD=event.JetCHS_puId_ptD[k]
        Before_pull=event.JetCHS_puId_pull[k]
        Before_nCharged=event.JetCHS_puId_nCharged[k]
        Before_jetRchg=event.JetCHS_puId_jetRchg[k]
        Before_beta=event.JetCHS_puId_beta[k]
 
        b_dR2Mean.Fill(Before_dR2Mean)
        b_frac01.Fill(Before_frac02)
        b_frac02.Fill(Before_frac02)
        b_dR2Mean.Fill(Before_dR2Mean)
        b_frac01.Fill(Before_frac02)
        b_frac02.Fill(Before_frac02)
        b_frac03.Fill(Before_frac03)
        b_frac04.Fill(Before_frac04)
        b_majW.Fill(Before_majW)
        b_minW.Fill(Before_minW)
        b_jetR.Fill(Before_jetR)
  #     b_nConstituents.Fill(Before_nConstituents,0)
        b_ptD.Fill(Before_ptD)
        b_pull.Fill(Before_pull)
        b_nCharged.Fill(Before_nCharged)
        b_jetRchg.Fill(Before_jetRchg)
        b_beta.Fill(Before_beta)
  b_npvs.Fill(ord(Before_npvs))

for event in after_Chain:
  After_npvs=event.PV_npvsGood
  if (jet_type == "PUPPI"):
    for k in range(event.nJet):
      if event.Jet_pt[k]>20 and event.Jet_eta[k]>eta_min and event.Jet_eta[k]<=eta_max:
        After_dR2Mean=event.Jet_puId_dR2Mean[k]
        After_frac01=event.Jet_puId_frac01[k]
        After_frac02=event.Jet_puId_frac02[k]
        After_frac03=event.Jet_puId_frac03[k]
        After_frac04=event.Jet_puId_frac04[k]
        After_majW=event.Jet_puId_majW[k]
        After_minW=event.Jet_puId_minW[k]
        After_jetR=event.Jet_puId_jetR[k]
 #      if ("0x" in event.Jet_nConstituents[k]):
  #        After_nConstituents= int(event.Jet_nConstituents[k],0)
  #      else:
  #        After_nConstituents= ord(event.Jet_nConstituents[k])
        After_ptD=event.Jet_puId_ptD[k]
        After_pull=event.Jet_puId_pull[k]
        After_nCharged=event.Jet_puId_nCharged[k]
        After_jetRchg=event.Jet_puId_jetRchg[k]
        After_beta=event.Jet_puId_beta[k]
 
        a_dR2Mean.Fill(After_dR2Mean)
        a_frac01.Fill(After_frac02)
        a_frac02.Fill(After_frac02)
        a_dR2Mean.Fill(After_dR2Mean)
        a_frac01.Fill(After_frac02)
        a_frac02.Fill(After_frac02)
        a_frac03.Fill(After_frac03)
        a_frac04.Fill(After_frac04)
        a_majW.Fill(After_majW)
        a_minW.Fill(After_minW)
        a_jetR.Fill(After_jetR)
  #     a_nConstituents.Fill(After_nConstituents)
        a_ptD.Fill(After_ptD)
        a_pull.Fill(After_pull)
        a_nCharged.Fill(After_nCharged)
        a_jetRchg.Fill(After_jetRchg)
        a_beta.Fill(After_beta)
  if (jet_type == "CHS"):
    for k in range(event.nJetCHS):
      if event.JetCHS_pt[k]>20 and event.JetCHS_eta[k]>eta_min and event.JetCHS_eta[k]<=eta_max:
        After_dR2Mean=event.JetCHS_puId_dR2Mean[k]
        After_frac01=event.JetCHS_puId_frac01[k]
        After_frac02=event.JetCHS_puId_frac02[k]
        After_frac03=event.JetCHS_puId_frac03[k]
        After_frac04=event.JetCHS_puId_frac04[k]
        After_majW=event.JetCHS_puId_majW[k]
        After_minW=event.JetCHS_puId_minW[k]
        After_jetR=event.JetCHS_puId_jetR[k]
#       if ("0x" in event.JetCHS_nConstituents[k]):
#         After_nConstituents= int(event.JetCHS_nConstituents[k],0)
#       else:
#         After_nConstituents= ord(event.JetCHS_nConstituents[k])
        After_ptD=event.JetCHS_puId_ptD[k]
        After_pull=event.JetCHS_puId_pull[k]
        After_nCharged=event.JetCHS_puId_nCharged[k]
        After_jetRchg=event.JetCHS_puId_jetRchg[k]
        After_beta=event.JetCHS_puId_beta[k]
 
        a_dR2Mean.Fill(After_dR2Mean)
        a_frac01.Fill(After_frac02)
        a_frac02.Fill(After_frac02)
        a_dR2Mean.Fill(After_dR2Mean)
        a_frac01.Fill(After_frac02)
        a_frac02.Fill(After_frac02)
        a_frac03.Fill(After_frac03)
        a_frac04.Fill(After_frac04)
        a_majW.Fill(After_majW)
        a_minW.Fill(After_minW)
        a_jetR.Fill(After_jetR)
#       a_nConstituents.Fill(After_nConstituents)
        a_ptD.Fill(After_ptD)
        a_pull.Fill(After_pull)
        a_nCharged.Fill(After_nCharged)
        a_jetRchg.Fill(After_jetRchg)
        a_beta.Fill(After_beta)
  a_npvs.Fill(ord(After_npvs))
 #after_file.Close()
#before_file.Close()

#a_var_list=[a_dR2Mean,a_frac01,a_frac02,a_frac03,a_frac04,a_majW,a_minW,a_jetR,a_ptD,a_pull,a_beta,a_jetRchg,a_nCharged,a_nConstituents, a_npvs]
#b_var_list=[b_dR2Mean,b_frac01,b_frac02,b_frac03,b_frac04,b_majW,b_minW,b_jetR,b_ptD,b_pull,b_beta,b_jetRchg,b_nCharged,b_nConstituents, b_npvs]

a_var_list=[a_dR2Mean,a_frac01,a_frac02,a_frac03,a_frac04,a_majW,a_minW,a_jetR,a_ptD,a_pull,a_beta,a_jetRchg,a_nCharged, a_npvs]
b_var_list=[b_dR2Mean,b_frac01,b_frac02,b_frac03,b_frac04,b_majW,b_minW,b_jetR,b_ptD,b_pull,b_beta,b_jetRchg,b_nCharged, b_npvs]

#var=["dR2Mean","frac_01","frac_02","frac_03","frac_04","majW","minW","jetR","ptD","pull","beta","jetRchg","Charged", "nConstituents","npvs"]
var=["dR2Mean","frac_01","frac_02","frac_03","frac_04","majW","minW","jetR","ptD","pull","beta","jetRchg","Charged","npvs"]


r_dR2Mean=a_dR2Mean.Clone()
r_frac01=a_frac01.Clone()
r_frac02=a_frac02.Clone()
r_frac03=a_frac03.Clone()
r_frac04=a_frac04.Clone()
r_majW=a_majW.Clone()
r_minW=a_minW.Clone()
r_jetR=a_jetR.Clone()
r_ptD=a_ptD.Clone()
r_pull=a_pull.Clone()
r_beta=a_beta.Clone()
r_jetRchg=a_jetRchg.Clone()
r_nCharged=a_nCharged.Clone()
#r_nConstituents= a_nConstituents.Clone()
r_npvs=a_npvs.Clone()

r_dR2Mean.Divide(b_dR2Mean)
r_frac01.Divide(b_frac01)
r_frac02.Divide(b_frac02)
r_frac03.Divide(b_frac03)
r_frac04.Divide(b_frac04)
r_majW.Divide(b_majW)
r_minW.Divide(b_minW)
r_jetR.Divide(b_jetR)
r_ptD.Divide(b_ptD)
r_pull.Divide(b_pull)
r_beta.Divide(b_beta)
r_jetRchg.Divide(b_jetRchg)
r_nCharged.Divide(b_nCharged)
#r_nConstituents.Divide(b_nConstituents)
r_npvs.Divide(b_npvs)



#r_var_list=[r_dR2Mean,r_frac01,r_frac02,r_frac03,r_frac04,r_majW,r_minW,r_jetR,r_ptD,r_pull,r_beta,r_jetRchg,r_nCharged,r_nConstituents, r_npvs]
r_var_list=[r_dR2Mean,r_frac01,r_frac02,r_frac03,r_frac04,r_majW,r_minW,r_jetR,r_ptD,r_pull,r_beta,r_jetRchg,r_nCharged, r_npvs]





if True:
   for i,name in enumerate(var):
      ROOT.gStyle.SetHatchesLineWidth(2)
      b_var_list[i].SetLineColor(ROOT.kBlue)
#      b_var_list[i].SetFillColorAlpha(ROOT.kBlue-7,0.35)
      b_var_list[i].SetStats(0)
      b_var_list[i].SetFillStyle(0)
      b_var_list[i].SetLineWidth(2)
      a_var_list[i].SetLineColor(ROOT.kRed)
#      a_var_list[i].SetFillColorAlpha(ROOT.kRed+1,0)
      a_var_list[i].SetStats(0)
      a_var_list[i].SetLineWidth(2)
      a_var_list[i].SetFillStyle(0)
      r_var_list[i].SetLineColor(ROOT.kRed)
      r_var_list[i].SetLineStyle(2)
      r_var_list[i].SetLineWidth(1)
      r_var_list[i].SetStats(0)

      b_var_list[i].GetXaxis().SetTitle(name)
      a_var_list[i].GetXaxis().SetTitle(name)
      b_var_list[i].GetYaxis().SetTitle("Events")
      a_var_list[i].GetYaxis().SetTitle("Events")
      r_var_list[i].GetYaxis().SetTitle("ratio")


      a_max_bin=a_var_list[i].GetMaximumBin()
      a_max=a_var_list[i].GetBinContent(a_max_bin)
      b_max_bin=b_var_list[i].GetMaximumBin()
      b_max=b_var_list[i].GetBinContent(b_max_bin)

      lgd = ROOT.TLegend(0.63, 0.9, 0.88, 0.75)
      lgd.SetBorderSize(1)
      lgd.SetFillColor(0)
      lgd.SetFillStyle(0)
      lgd.SetTextFont(42)
      lgd.SetTextSize(0.030)
      lgd.AddEntry(b_var_list[i],"Before")
      lgd.AddEntry(a_var_list[i],"After")
      c=ROOT.TCanvas("canvas","canvas",500,500)
      c.cd()
      pad1 = ROOT.TPad("pad1","pad1",0,0.2,1,1) 
      pad1.Draw("H")
      pad1.SetBottomMargin(0)
      pad1.cd()     
      if a_max>b_max:
         a_var_list[i].Draw("HIST")
         b_var_list[i].Draw("HIST SAME")
         a_var_list[i].Draw("HIST SAME")
      if a_max <= b_max:
         b_var_list[i].Draw("HIST")
         a_var_list[i].Draw("HIST SAME")
      lgd.Draw("SAME")
      a_var_list[i].GetXaxis().SetLabelSize(0)
      b_var_list[i].GetXaxis().SetLabelSize(0)
      a_var_list[i].GetXaxis().SetTitleSize(0)
      b_var_list[i].GetXaxis().SetTitleSize(0)


      c.cd()
      pad2 = ROOT.TPad("pad2","pad2",0,0.02,1,0.2) 
      pad2.Draw()
      pad2.cd()
      pad2.SetTopMargin(0)
      pad2.SetBottomMargin(0.25)
      r_var_list[i].SetTitle("")
      r_var_list[i].GetYaxis().SetRangeUser(0.5,1.5)
      r_var_list[i].GetYaxis().SetLabelSize(0.13)
      r_var_list[i].GetYaxis().SetTitleSize(0.2)
      r_var_list[i].GetYaxis().SetTitle("After/Before")
      r_var_list[i].GetYaxis().SetNdivisions(207)
      r_var_list[i].GetXaxis().SetLabelSize(0.15)
      r_var_list[i].Draw("Hist")
      c.SaveAs("%s_Validation_plot_%s_%s_%s.png" %(process_name.upper(),name,jet_type.upper(),eta_s))
      c.Clear()



