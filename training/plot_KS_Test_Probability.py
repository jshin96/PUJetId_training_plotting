import ROOT
from array import array
eta_bins=["0p0To2p5","2p5To2p75","2p75To3p0","3p0To5p0"]
for eta_bin in eta_bins:
   x=array('d')
   sig_y=array('d')
   bkg_y=array('d')
   for NTree in range(50, 1001, 50):
      inputfile=ROOT.TFile.Open("output/BDT_puppi_106X_2018/tmva_output_Pt10To100Eta%s_puppi_NTrees_%i.root" %(eta_bin,NTree), "READ")

      sig_hist_test=inputfile.Get("BDT_puppi_106X_2018/Method_BDT/BDT/MVA_BDT_S")
      sig_hist_train=inputfile.Get("BDT_puppi_106X_2018/Method_BDT/BDT/MVA_BDT_Train_S")
      
      bkg_hist_test=inputfile.Get("BDT_puppi_106X_2018/Method_BDT/BDT/MVA_BDT_B")
      bkg_hist_train=inputfile.Get("BDT_puppi_106X_2018/Method_BDT/BDT/MVA_BDT_Train_B")


      sig_SK_score=sig_hist_train.KolmogorovTest(sig_hist_test,"X")
      bkg_SK_score=bkg_hist_train.KolmogorovTest(bkg_hist_test,"X")

      x.append(NTree)
      bkg_y.append(sig_SK_score)
      sig_y.append(bkg_SK_score)


      inputfile.Close()

   sig_graph=ROOT.TGraph(20,x,sig_y)
   bkg_graph=ROOT.TGraph(20,x,bkg_y)
   c=ROOT.TCanvas("canvas","canvas",1000,500)
   c.Divide(2,1)
   c.cd(1)
   sig_graph.SetTitle('signal SK score Eta%s' %eta_bin)
   sig_graph.Draw("AL*")
   c.cd(2)
   bkg_graph.SetTitle('background SK score Eta%s' %eta_bin)
   bkg_graph.Draw("AL*")
   c.SaveAs("SK_Test_Probability_Eta%s.png" %eta_bin)

