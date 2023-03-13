import ROOT


eta_bins=["Eta0p0To2p5","Eta2p5To2p75","Eta2p75To3p0","Eta3p0To5p0"]
eta_bin_range=["[0.0,2.5]","[2.5,2.75]","[2.75,3.0]","[3.0,5.0]"]
input_files=[]
for idx in range(65):
   inputfile="input/PUPPI_RUN3_2022/training_trees_pt10To100_puppi_130X_%i_2022.root" %idx
   input_files.append(inputfile)

for x, eta_bin in enumerate(eta_bins):
   nsig=0
   nbkg=0
 
   s_npvs=ROOT.TH1D("npvs","npvsGood #eta #in %s" %eta_bin_range[x],70,0,70)
   s_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   s_dR2Mean=ROOT.TH1D("dR2Mean","dR2Mean #eta #in %s" %eta_bin_range[x],30,0,0.1)
   s_frac01=ROOT.TH1D("frac_01","frac_01 #eta #in %s" %eta_bin_range[x],30,0,1)
   s_frac02=ROOT.TH1D("frac_02","frac_02 #eta #in %s" %eta_bin_range[x],30,0,1)
   s_frac03=ROOT.TH1D("frac_03","frac_03 #eta #in %s" %eta_bin_range[x],30,0,1)
   s_frac04=ROOT.TH1D("frac_04","frac_04 #eta #in %s" %eta_bin_range[x],30,0,1)
   s_majW=ROOT.TH1D("majW","majW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   s_minW=ROOT.TH1D("minW","minW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   s_jetR=ROOT.TH1D("jetR","jetR #eta #in %s" %eta_bin_range[x],30,0,1)
   s_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   s_nConstituents=ROOT.TH1D("nConstituents","nConstituents #eta #in %s" %eta_bin_range[x],40,0,40)
   s_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
   s_ptD=ROOT.TH1D("ptD","ptD #eta #in %s" %eta_bin_range[x],30,0,1)
   s_pull=ROOT.TH1D("pull","pull #eta #in %s" %eta_bin_range[x],30,0,0.05)
   s_Rho=ROOT.TH1D("Rho","Rho #eta #in %s" %eta_bin_range[x],70,0,70)
   s_pt=ROOT.TH1D("Jet_Pt","Jet p_{T} #eta #in %s" %eta_bin_range[x],100,0,100)
   s_eta=ROOT.TH1D("Jet_Eta","Jet #eta #eta #in %s" %eta_bin_range[x],100,-5,5)

   b_npvs=ROOT.TH1D("npvs","npvsGood #eta #in %s" %eta_bin_range[x],70,0,70)
   b_beta=ROOT.TH1D("beta","beta #eta #in %s" %eta_bin_range[x],30,0,1)
   b_dR2Mean=ROOT.TH1D("dR2Mean","dR2Mean #eta #in %s" %eta_bin_range[x],30,0,0.1)
   b_frac01=ROOT.TH1D("frac_01","frac_01 #eta #in %s" %eta_bin_range[x],30,0,1)
   b_frac02=ROOT.TH1D("frac_02","frac_02 #eta #in %s" %eta_bin_range[x],30,0,1)
   b_frac03=ROOT.TH1D("frac_03","frac_03 #eta #in %s" %eta_bin_range[x],30,0,1)
   b_frac04=ROOT.TH1D("frac_04","frac_04 #eta #in %s" %eta_bin_range[x],30,0,1)
   b_majW=ROOT.TH1D("majW","majW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   b_minW=ROOT.TH1D("minW","minW #eta #in %s" %eta_bin_range[x],30,0,0.4)
   b_jetR=ROOT.TH1D("jetR","jetR #eta #in %s" %eta_bin_range[x],30,0,1)
   b_jetRchg=ROOT.TH1D("jetRchg","jetRchg #eta #in %s" %eta_bin_range[x],30,0,1)
   b_nConstituents=ROOT.TH1D("nConstituents","nConstituents #eta #in %s" %eta_bin_range[x],40,0,40)
   b_Charged=ROOT.TH1D("nCharged","nCharged #eta #in %s" %eta_bin_range[x],30,0,30)
   b_ptD=ROOT.TH1D("ptD","ptD #eta #in %s" %eta_bin_range[x],30,0,1)
   b_pull=ROOT.TH1D("pull","pull #eta #in %s" %eta_bin_range[x],30,0,0.05)
   b_Rho=ROOT.TH1D("Rho","Rho #eta #in %s" %eta_bin_range[x],70,0,70)
   b_pt=ROOT.TH1D("Jet_Pt","Jet p_{T} #eta #in %s" %eta_bin_range[x],100,0,100)
   b_eta=ROOT.TH1D("Jet_Eta","Jet #eta #eta #in %s" %eta_bin_range[x],100,-5,5)

   for input_file_name in pre_input_files:
      input_file=ROOT.TFile.Open(input_file_name,"READ")
      sig_tree=input_file.Get("%s_Prompt" %eta_bin) 
      for entry in range(sig_tree.GetEntries()): 
         nsig+=1
         if nsig>250000: break
         sig_tree.GetEntry(entry)
         Prompt_npvs=getattr(sig_tree,"PV_npvsGood")
         Prompt_dR2Mean=getattr(sig_tree,"Jet_puId_dR2Mean")
         Prompt_frac01=getattr(sig_tree,"Jet_puId_frac01")
         Prompt_frac02=getattr(sig_tree,"Jet_puId_frac02")
         Prompt_frac03=getattr(sig_tree,"Jet_puId_frac03")
         Prompt_frac04=getattr(sig_tree,"Jet_puId_frac04")
         Prompt_majW=getattr(sig_tree,"Jet_puId_majW")
         Prompt_minW=getattr(sig_tree,"Jet_puId_minW")
         Prompt_jetR=getattr(sig_tree,"Jet_puId_jetR")
         Prompt_nConstituents=getattr(sig_tree,"Jet_nConstituents")
         Prompt_ptD=getattr(sig_tree,"Jet_puId_ptD")
         Prompt_pull=getattr(sig_tree,"Jet_puId_pull")
         Prompt_Rho=getattr(sig_tree,"fixedGridRhoFastjetAll")
         Prompt_Pt=getattr(sig_tree,"Jet_pt")
         Prompt_Eta=getattr(sig_tree,"Jet_eta")
         if eta_bin!="Eta3p0To5p0":
            Prompt_Charged=getattr(sig_tree,"Jet_puId_nCharged")
            Prompt_jetRchg=getattr(sig_tree,"Jet_puId_jetRchg")
            Prompt_beta=getattr(sig_tree,"Jet_puId_beta")
         s_npvs.Fill(Prompt_npvs)
         s_dR2Mean.Fill(Prompt_dR2Mean)
         s_frac01.Fill(Prompt_frac01)
         s_frac02.Fill(Prompt_frac02)
         s_frac03.Fill(Prompt_frac03)
         s_frac04.Fill(Prompt_frac04)
         s_majW.Fill(Prompt_majW)
         s_minW.Fill(Prompt_minW)
         s_jetR.Fill(Prompt_jetR)
         s_nConstituents.Fill(Prompt_nConstituents)
         s_ptD.Fill(Prompt_ptD)
         s_pull.Fill(Prompt_pull)
         s_Rho.Fill(Prompt_Rho)
         s_pt.Fill(Prompt_Pt)
         s_eta.Fill(Prompt_Eta)
         if eta_bin!="Eta3p0To5p0":
            s_beta.Fill(Prompt_beta)
            s_jetRchg.Fill(Prompt_jetRchg)
            s_Charged.Fill(Prompt_Charged)
 
      bkg_tree=input_file.Get("%s_Pileup" %eta_bin)
      for entry in range(bkg_tree.GetEntries()): 
         nbkg+=1
         if nbkg>250000: break
         bkg_tree.GetEntry(entry)
         Pileup_npvs=getattr(bkg_tree,"PV_npvsGood")
         Pileup_dR2Mean=getattr(bkg_tree,"Jet_puId_dR2Mean")
         Pileup_frac01=getattr(bkg_tree,"Jet_puId_frac01")
         Pileup_frac02=getattr(bkg_tree,"Jet_puId_frac02")
         Pileup_frac03=getattr(bkg_tree,"Jet_puId_frac03")
         Pileup_frac04=getattr(bkg_tree,"Jet_puId_frac04")
         Pileup_majW=getattr(bkg_tree,"Jet_puId_majW")
         Pileup_minW=getattr(bkg_tree,"Jet_puId_minW")
         Pileup_jetR=getattr(bkg_tree,"Jet_puId_jetR")
         Pileup_nConstituents=getattr(bkg_tree,"Jet_nConstituents")
         Pileup_ptD=getattr(bkg_tree,"Jet_puId_ptD")
         Pileup_pull=getattr(bkg_tree,"Jet_puId_pull")
         Pileup_Rho=getattr(bkg_tree,"fixedGridRhoFastjetAll")
         Pileup_Pt=getattr(bkg_tree,"Jet_pt")
         Pileup_Eta=getattr(bkg_tree,"Jet_eta")
         if eta_bin!="Eta3p0To5p0":
            Pileup_Charged=getattr(bkg_tree,"Jet_puId_nCharged")
            Pileup_jetRchg=getattr(bkg_tree,"Jet_puId_jetRchg")
            Pileup_beta=getattr(bkg_tree,"Jet_puId_beta")
         b_npvs.Fill(Pileup_npvs)
         b_dR2Mean.Fill(Pileup_dR2Mean)
         b_frac01.Fill(Pileup_frac01)
         b_frac02.Fill(Pileup_frac02)
         b_frac03.Fill(Pileup_frac03)
         b_frac04.Fill(Pileup_frac04)
         b_majW.Fill(Pileup_majW)
         b_minW.Fill(Pileup_minW)
         b_jetR.Fill(Pileup_jetR)
         b_nConstituents.Fill(Pileup_nConstituents)
         b_ptD.Fill(Pileup_ptD)
         b_pull.Fill(Pileup_pull)
         b_Rho.Fill(Pileup_Rho)
         b_pt.Fill(Pileup_Pt)
         b_eta.Fill(Pileup_Eta)
         if eta_bin!="Eta3p0To5p0":
            b_beta.Fill(Pileup_beta)
            b_jetRchg.Fill(Pileup_jetRchg)
            b_Charged.Fill(Pileup_Charged)
      input_file.Close()





         


   
   if eta_bin!="Eta3p0To5p0":
      signal=[s_npvs,s_dR2Mean,s_frac01,s_frac02,s_frac03,s_frac04,s_majW,s_minW,s_jetR,s_nConstituents,s_ptD,s_pull,s_Rho,s_beta,s_jetRchg,s_Charged,s_pt, s_eta]
      background=[b_npvs,b_dR2Mean,b_frac01,b_frac02,b_frac03,b_frac04,b_majW,b_minW,b_jetR,b_nConstituents,b_ptD,b_pull,b_Rho,b_beta,b_jetRchg,b_Charged,b_pt,b_eta]
   if eta_bin=="Eta3p0To5p0":
      signal=[s_npvs,s_dR2Mean,s_frac01,s_frac02,s_frac03,s_frac04,s_majW,s_minW,s_jetR,s_nConstituents,s_ptD,s_pull,s_Rho, s_pt, s_eta]
      background=[b_npvs,b_dR2Mean,b_frac01,b_frac02,b_frac03,b_frac04,b_majW,b_minW,b_jetR,b_nConstituents,b_ptD,b_pull,b_Rho, b_pt, b_eta]
   
   if eta_bin!="Eta3p0To5p0":
      var=["npvs","dR2Mean","frac_01","frac_02","frac_03","frac_04","majW","minW","jetR","nConstitunents","ptD","pull","Rho","beta","jetRchg","nCharged", "Jet p_{T}", "Jet #eta"]
   if eta_bin=="Eta3p0To5p0":
      var=["npvs","dR2Mean","frac_01","frac_02","frac_03","frac_04","majW","minW","jetR","nConstitunents","ptD","pull","Rho","Jet p_{T}","Jet #eta"]

   for i,name in enumerate(var):
      ROOT.gStyle.SetHatchesLineWidth(2)
      signal[i].SetLineColor(ROOT.kBlue)
      signal[i].SetFillColorAlpha(ROOT.kBlue-7,0.35)
      signal[i].SetStats(0)
      signal[i].SetFillStyle(1001)
      background[i].SetLineColor(ROOT.kRed)
      background[i].SetFillColor(ROOT.kRed+1)
      background[i].SetStats(0)
      background[i].SetLineWidth(2)
      background[i].SetFillStyle(3554)

      signal[i].GetXaxis().SetTitle(name)
      background[i].GetXaxis().SetTitle(name)
      signal[i].GetYaxis().SetTitle("Events")
      background[i].GetYaxis().SetTitle("Events")


      sig_max_bin=signal[i].GetMaximumBin()
      sig_max=signal[i].GetBinContent(sig_max_bin)
      bkg_max_bin=background[i].GetMaximumBin()
      bkg_max=background[i].GetBinContent(bkg_max_bin)

      lgd = ROOT.TLegend(0.63, 0.9, 0.88, 0.75)
      lgd.SetBorderSize(1)
      lgd.SetFillColor(0)
      lgd.SetFillStyle(0)
      lgd.SetTextFont(42)
      lgd.SetTextSize(0.030)
      lgd.AddEntry(signal[i],"Prompt")
      lgd.AddEntry(background[i],"Pileup")
      c=ROOT.TCanvas("canvas","canvas",500,500)
      c.cd()
      if sig_max>bkg_max:
         signal[i].Draw("HIST")
         background[i].Draw("HIST SAME")
      if sig_max <= bkg_max:
         background[i].Draw("HIST")
         signal[i].Draw("HIST SAME")
         background[i].Draw("HIST SAME")
      lgd.Draw("SAME")
      c.SaveAs("RUN3_plot_%s_%s.png" %(name,eta_bin))


