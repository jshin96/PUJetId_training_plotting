#!/usr/bin/env python3

import ROOT
import os
import sys
from array import array
import argparse

ROOT.gROOT.SetBatch(True)

if "HOME" not in os.environ:
    os.environ["HOME"] = "" # needed for condor batch mode

parser = argparse.ArgumentParser(
    description="analyze and make hists from JME flat ntuple"
    )
parser.add_argument(
    "--era", type=str, default="106X", help="data Era, default=%(default)s"
    )
parser.add_argument(
    "--jet_type", type=str, default="puppi", help="jet type, default=%(default)s"
    )
parser.add_argument(
    "--data_type", type=str, default="mc", help="data or mc, default=%(default)s"
    )
parser.add_argument(
    "--year", type=str, default="2018",
    help="Data year: 2017 or 2018, default=%(default)s"
    )
parser.add_argument(
    "--period", type=str, default="F",
    help="Data period: A,B,C,D... , default=%(default)s"
    )
parser.add_argument(
    "--lumi", type=float, default=63.67,
    help="Data integrated lumi (1/fb) , default=%(default)s"
    )
parser.add_argument(
    "--N_mc", type=int, default=48675378,
    help="# of generated evts of MC sample, default=%(default)s"
    )
parser.add_argument(
    "--xs", type=float, default=5343.0,
    help="Cross Section (pb) of MC sample, default=%(default)s"
    )
parser.add_argument(
    "--input_filename", type=str, help="input root file location"
    )
parser.add_argument(
    "--output_filename", type=str, help="output root file name"
    )

args = parser.parse_args()

era = args.era
jet_type = args.jet_type
data_type = args.data_type
year = args.year
period = args.period
lumi = args.lumi
N_mc = args.N_mc
xs = args.xs
input_filename = args.input_filename
output_filename = args.output_filename

xs_weight = lumi * xs * 1000 / N_mc

print("====================")
print("DataSet:\t", data_type)
print("Era, Year, Period:\t", era, year, period)
print("====================")
print("XS Weight: %s (1/fb) * 1000 * %s (pb) / (%s) = %s" % (lumi, xs, N_mc, xs_weight))
print("====================")
print("Input File:\t", input_filename)
print("Output File:\t", output_filename)

input_tree = "Events"
#tChain = ROOT.TChain("Events")
#tChain.Add(input_filename)
input_file = ROOT.TFile.Open(input_filename)
#for i in tChain:
#    events = i
events = input_file.Get("Events")

eta_bins = [
    "Eta0p0To2p5",
    "Eta2p5To2p75",
    "Eta2p75To3p0",
    "Eta3p0To5p0"
]

f_eta_bins = [
    tuple(float(x) for x in eta_bin.replace("Eta", "").replace("p", ".").split("To")) for eta_bin in eta_bins
]

pt_bins = [
    "Pt10To20",
    "Pt20To30",
    "Pt30To40",
    "Pt40To50",
    "Pt50To100",
]

f_pt_bins = [
    tuple(float(x) for x in pt_bin.replace("Pt", "").split("To")) for pt_bin in pt_bins
]

tmva_weight_filenames = [
    "tmva_weights_%s/pileupJetId_%s_%s_%s_BDT.weights.xml" % (year, era, i, jet_type) for i in eta_bins
]

# uncomment below for testing old training,
# and comment above ones

#tmva_weight_filenames = [
#    "old_training_files/pileupJetId_80XvarFix_Eta0to2p5_BDT.weights.xml",
#    "old_training_files/pileupJetId_80XvarFix_Eta2p5to2p75_BDT.weights.xml",
#    "old_training_files/pileupJetId_80XvarFix_Eta2p75to3_BDT.weights.xml",
#    "old_training_files/pileupJetId_80XvarFix_Eta3to5_BDT.weights.xml",
#]

tmva_s_jetPt = array("f", [-999])
tmva_s_jetEta = array("f", [-999])

tmva_spectators = [
    ("JetPuppi_pt", tmva_s_jetPt),
    ("JetPuppi_eta", tmva_s_jetEta),
]

tmva_v_nvtx = array("f", [-999])
tmva_v_beta = array("f", [-999])
tmva_v_dR2Mean = array("f", [-999])
tmva_v_frac01 = array("f", [-999])
tmva_v_frac02 = array("f", [-999])
tmva_v_frac03 = array("f", [-999])
tmva_v_frac04 = array("f", [-999])
tmva_v_majW = array("f", [-999])
tmva_v_minW = array("f", [-999])
tmva_v_jetR = array("f", [-999])
tmva_v_jetRchg = array("f", [-999])
tmva_v_nConstituents = array("f", [-999])
tmva_v_nCharged = array("f", [-999])
tmva_v_ptD = array("f", [-999])
tmva_v_pull = array("f", [-999])
tmva_v_fixedGridRhoFastjetAll = array("f", [-999])


tmva_variables = [
    ("PV_npvsGood" , tmva_v_nvtx),
    ("JetPuppi_puId_beta", tmva_v_beta),
    ("JetPuppi_puId_dR2Mean", tmva_v_dR2Mean),
    ("JetPuppi_puId_frac01", tmva_v_frac01),
    ("JetPuppi_puId_frac02", tmva_v_frac02),
    ("JetPuppi_puId_frac03", tmva_v_frac03),
    ("JetPuppi_puId_frac04", tmva_v_frac04),
    ("JetPuppi_puId_majW", tmva_v_majW),
    ("JetPuppi_puId_minW", tmva_v_minW),
    ("JetPuppi_puId_jetR", tmva_v_jetR),
    ("JetPuppi_puId_jetRchg", tmva_v_jetRchg),
    ("JetPuppi_nConstituents", tmva_v_nConstituents),
    ("JetPuppi_puId_nCharged", tmva_v_nCharged),
    ("JetPuppi_puId_ptD", tmva_v_ptD),
    ("JetPuppi_puId_pull", tmva_v_pull),
    ("fixedGridRhoFastjetAll", tmva_v_fixedGridRhoFastjetAll),
]

# uncomment below for testing old training,
# and comment above ones

#tmva_variables = [
#    ("nvtx" , tmva_v_nvtx),
#    ("dR2Mean", tmva_v_dR2Mean),
#    ("nConstituents", tmva_v_nConstituents),
#    ("nCharged", tmva_v_nCharged),
#    ("majW", tmva_v_majW),
#    ("minW", tmva_v_minW),
#    ("frac01", tmva_v_frac01),
#    ("frac02", tmva_v_frac02),
#    ("frac03", tmva_v_frac03),
#    ("frac04", tmva_v_frac04),
#    ("ptD", tmva_v_ptD),
#    ("beta", tmva_v_beta),
#    ("pull", tmva_v_pull),
#    ("jetR", tmva_v_jetR),
#    ("jetRchg", tmva_v_jetRchg),
#]

tmva_readers = []

for i, eta_bin in enumerate(eta_bins):
    
    tmva_readers.append(ROOT.TMVA.Reader())
    
    for spec_name, spec_address in tmva_spectators:
        
        tmva_readers[i].AddSpectator(spec_name, spec_address)
    
    for var_name, var_address in tmva_variables:
        
        if eta_bin == "Eta3p0To5p0":
            if var_name == "JetPuppi_puId_beta": continue
            if var_name == "JetPuppi_puId_jetRchg": continue
            if var_name == "JetPuppi_puId_nCharged": continue
        tmva_readers[i].AddVariable(var_name, var_address)
    
    tmva_readers[i].BookMVA("BDT", tmva_weight_filenames[i])

# some helping functions
def dphi(phi1, phi2):
    dphi = phi1 - phi2
    Pi = ROOT.TMath.Pi()
    if dphi > Pi:
        dphi -= 2 * Pi
    elif dphi <= -Pi:
        dphi += 2 * Pi
    return dphi

# book histograms
class book_hist_dict:
    def __init__(self, xbins, xlow, xup, titleX, units="", name="", 
                 keys=["main"], keys_sub=["sub_category"]):
        self.xbins = xbins
        self.xlow = xlow
        self.xup = xup
        self.titleX = titleX
        self.units = units
        self.name = name
        self.keys = keys
        self.keys_sub = keys_sub
    
    def hist_1D(self):
        variable = ROOT.TH1F("", "", self.xbins, self.xlow, self.xup)
        bw = variable.GetBinWidth(1)
        
        titleX = self.titleX
        titleY = "Enteries/%s" % bw
        
        if self.units != "":
            titleX = self.titleX + " [" + self.units + "]"
            titleY = titleY + " " + self.units
        
        variable.SetTitle("%s;%s;%s" % (titleX, titleX, titleY))
        
        if self.name != "":
            variable.SetName(self.name)
        else:
            variable.SetName(self.titleX)
        
        return variable
    
    def clone(self):
        hist_dict = {}
        
        for key in self.keys:
            name_ = key + "_" + self.hist_1D().GetName()
            hist_dict[key] = self.hist_1D().Clone()
            hist_dict[key].SetName(name_)

            for key_sub in self.keys_sub:
                name = name_ + "_" + key_sub
                hist_dict[key + "_" + key_sub] = self.hist_1D().Clone()
                hist_dict[key + "_" + key_sub].SetName(name)
                
        return hist_dict

# keys for per event histograms
# should be just one element which is "mc" or "data", hardcoded
h_keys1 = [data_type]

h_nvtx = book_hist_dict(100, 0, 100, "nvtx", keys=h_keys1).clone()

h_z_mass = book_hist_dict(40, 70, 110, "Z_mass", keys=h_keys1).clone()
h_z_pt = book_hist_dict(60, 0, 300, "Z_pt", keys=h_keys1).clone()
h_z_phi = book_hist_dict(32, -3.2, 3.2, "Z_phi", keys=h_keys1).clone()
h_z_rapidity = book_hist_dict(30, -3.0, 3.0, "Z_rapidity", keys=h_keys1).clone()

h_lept_pt1 = book_hist_dict(50, 0, 250, "lept_pt1", keys=h_keys1).clone()
h_lept_pt2 = book_hist_dict(50, 0, 250, "lept_pt2", keys=h_keys1).clone()
h_lept_eta1 = book_hist_dict(60, -3.0, 3.0, "lept_eta1", keys=h_keys1).clone()
h_lept_eta2 = book_hist_dict(60, -3.0, 3.0, "lept_eta2", keys=h_keys1).clone()
h_lept_phi1 = book_hist_dict(64, -3.2, 3.2, "lept_phi1", keys=h_keys1).clone()
h_lept_phi2 = book_hist_dict(64, -3.2, 3.2, "lept_phi2", keys=h_keys1).clone()

h_dphi_zj = book_hist_dict(64, -3.2, 3.2, "dphi_zj", keys=h_keys1).clone()
h_n_jets = book_hist_dict(40, 0, 40, "n_jets", keys=h_keys1).clone()
    
# keys for per jet histograms
h_keys2 = eta_bins

h_keys3 = [
    eta_bin + "_" + pt_bin for eta_bin in eta_bins for pt_bin in pt_bins
]

if data_type == "mc":

    jet_tags = ["all", "prompt", "pileup", "rest"]

    h_keys2 = [
        i + "_" + j for i in h_keys2 for j in jet_tags
    ]

    h_keys3 = [
        i + "_" + j for i in h_keys3 for j in jet_tags
    ]

h_keys4 = h_keys2 + h_keys3

h_jet_pt = book_hist_dict(100, 0.0, 100.0, "jet_pt", units="GeV", keys=h_keys1, keys_sub=h_keys4).clone()
h_jet_eta = book_hist_dict(50, -5.0, 5.0, "jet_eta", keys=h_keys1, keys_sub=h_keys4).clone()
h_jet_energy = book_hist_dict(100, 0.0, 500.0, "jet_energy", units="GeV", keys=h_keys1, keys_sub=h_keys4).clone()
h_jet_phi = book_hist_dict(64, -3.2, 3.2, "jet_phi", keys=h_keys1, keys_sub=h_keys4).clone()
# min dR jet and genjet
if data_type == "mc":
    h_jet_dR = book_hist_dict(250, 0.0, 5.0, "jet_dR", keys=h_keys1, keys_sub=h_keys4).clone()
#chargedHadronEnergyFraction
h_jet_chf= book_hist_dict(52, -0.02, 1.02, "jet_chf", keys=h_keys1, keys_sub=h_keys4).clone()
#chargedEmEnergyFraction
h_jet_cemf = book_hist_dict(52, -0.02, 1.02, "jet_cemf", keys=h_keys1, keys_sub=h_keys4).clone()
#neutralHadronEnergyFraction
h_jet_nhf = book_hist_dict(52, -0.02, 1.02, "jet_nhf", keys=h_keys1, keys_sub=h_keys4).clone()
#neutralEmEnergyFraction
h_jet_nemf = book_hist_dict(52, -0.02, 1.02, "jet_nemf", keys=h_keys1, keys_sub=h_keys4).clone()
#photonEnergyFraction
h_jet_phf = book_hist_dict(52, -0.02, 1.02, "jet_phf", keys=h_keys1, keys_sub=h_keys4).clone()
#muonEnergyFraction
h_jet_muf = book_hist_dict(52, -0.02, 1.02, "jet_muf", keys=h_keys1, keys_sub=h_keys4).clone()
#electronEnergyFraction
h_jet_elf = book_hist_dict(52, -0.02, 1.02, "jet_elf", keys=h_keys1, keys_sub=h_keys4).clone()
#neutralMultiplicity
h_jet_npr = book_hist_dict(60, 0, 60, "jet_npr", keys=h_keys1, keys_sub=h_keys4).clone()
#chargedHadronMultiplicity
h_jet_chm = book_hist_dict(40, 0, 40, "jet_chm", keys=h_keys1, keys_sub=h_keys4).clone()
#neutralHadronMultiplicity
h_jet_nhm = book_hist_dict(10, 0, 10, "jet_nhm", keys=h_keys1, keys_sub=h_keys4).clone()
#photonMultiplicity
h_jet_phm = book_hist_dict(40, 0, 40, "jet_phm", keys=h_keys1, keys_sub=h_keys4).clone()
#muonMultiplicity
h_jet_mum = book_hist_dict(5, 0, 5, "jet_mum", keys=h_keys1, keys_sub=h_keys4).clone()
#electronMultiplicity
h_jet_elm = book_hist_dict(5, 0, 5, "jet_elm", keys=h_keys1, keys_sub=h_keys4).clone()
# ------------
# BDT variables
h_jet_beta = book_hist_dict(52, -0.02, 1.02, "beta", keys=h_keys1, keys_sub=h_keys4).clone()
h_dR2Mean = book_hist_dict(60, -0.01, 0.11, "dR2Mean", keys=h_keys1, keys_sub=h_keys4).clone()
h_frac01 = book_hist_dict(52, -0.02, 1.02, "frac01", keys=h_keys1, keys_sub=h_keys4).clone()
h_frac02 = book_hist_dict(52, -0.02, 1.02, "frac02", keys=h_keys1, keys_sub=h_keys4).clone()
h_frac03 = book_hist_dict(52, -0.02, 1.02, "frac03", keys=h_keys1, keys_sub=h_keys4).clone()
h_frac04 = book_hist_dict(52, -0.02, 1.02, "frac04", keys=h_keys1, keys_sub=h_keys4).clone()
h_jetRchg = book_hist_dict(52, -0.02, 1.02, "jetRchg", keys=h_keys1, keys_sub=h_keys4).clone()
h_jetR = book_hist_dict(52, -0.02, 1.02, "jetR", keys=h_keys1, keys_sub=h_keys4).clone()
h_majW = book_hist_dict(60, 0.0, 0.3, "majW", keys=h_keys1, keys_sub=h_keys4).clone()
h_minW = book_hist_dict(40, 0.0, 0.2, "minW", keys=h_keys1, keys_sub=h_keys4).clone()
h_nCharged = book_hist_dict(40, 0, 40, "nCharged", keys=h_keys1, keys_sub=h_keys4).clone()
h_nConstituents = book_hist_dict(60, 0, 60, "nConstituents", keys=h_keys1, keys_sub=h_keys4).clone()
h_ptD = book_hist_dict(50, 0.0, 1.0, "ptD", keys=h_keys1, keys_sub=h_keys4).clone()
h_pull = book_hist_dict(50, 0.0, 0.025, "pull", keys=h_keys1, keys_sub=h_keys4).clone()
h_fixedGridRhoFastjetAll = book_hist_dict(100, 0.0, 60.0, "fixedGridRhoFastjetAll", keys=h_keys1, keys_sub=h_keys4).clone()

#h_old_mva = book_hist_dict(200, -1.0, 1.0, "old_mva", keys=h_keys1, keys_sub=h_keys4).clone()
h_new_mva = book_hist_dict(200, -1.0, 1.0, "new_mva", keys=h_keys1, keys_sub=h_keys4).clone()


control_wp_keys = [
    "pass_wp1", "fail_wp1", 
    "pass_wp2", "fail_wp2",
    "pass_wp3", "fail_wp3",
]

h_control_jet_pt = book_hist_dict(80, 20, 100, "control_jet_pt", keys=["tight", "loose"], keys_sub=control_wp_keys).clone()
h_control_jet_eta = book_hist_dict(50, -5.0, 5.0, "control_jet_eta", keys=["tight", "loose"], keys_sub=control_wp_keys).clone()

h_dphi_zj_ptj_z_ = ROOT.TH2F("dphi_zj_ptj_z", "dphi_zj_ptj_z", 64, -3.2, 3.2, 100, 0, 10)

if data_type == "mc":
    h_dphi_zj_ptj_z = {key: h_dphi_zj_ptj_z_.Clone(f"{key}_dphi_zj_ptj_z") for key in ["prompt", "pileup"]}

else:
    h_dphi_zj_ptj_z = {key: h_dphi_zj_ptj_z_.Clone(f"{key}_dphi_zj_ptj_z") for key in ["data"]}

    
def write_hists(k=""):
    # per event
    if k in h_nvtx: h_nvtx[k].Write()
    if k in h_z_mass: h_z_mass[k].Write()
    if k in h_z_pt: h_z_pt[k].Write()
    if k in h_z_phi: h_z_phi[k].Write()
    if k in h_z_rapidity: h_z_rapidity[k].Write()
    if k in h_lept_pt1: h_lept_pt1[k].Write()
    if k in h_lept_pt2: h_lept_pt2[k].Write()
    if k in h_lept_eta1: h_lept_eta1[k].Write()
    if k in h_lept_eta2: h_lept_eta2[k].Write()
    if k in h_lept_phi1: h_lept_phi1[k].Write()
    if k in h_lept_phi2: h_lept_phi2[k].Write()
    if k in h_dphi_zj: h_dphi_zj[k].Write()
    if k in h_n_jets: h_n_jets[k].Write()
    
    # per jet
    if k in h_jet_pt: h_jet_pt[k].Write()
    if k in h_jet_eta: h_jet_eta[k].Write()
    if k in h_jet_phi: h_jet_phi[k].Write()
    if k in h_jet_energy: h_jet_energy[k].Write()
    if data_type == "mc":
        if k in h_jet_dR: h_jet_dR[k].Write()
    if k in h_jet_chf: h_jet_chf[k].Write()
    if k in h_jet_cemf: h_jet_cemf[k].Write()
    if k in h_jet_nhf: h_jet_nhf[k].Write()
    if k in h_jet_nemf: h_jet_nemf[k].Write()
    if k in h_jet_phf: h_jet_phf[k].Write()
    if k in h_jet_muf: h_jet_muf[k].Write()
    if k in h_jet_elf: h_jet_elf[k].Write()
    if k in h_jet_npr: h_jet_npr[k].Write()
    if k in h_jet_chm: h_jet_chm[k].Write()
    if k in h_jet_nhm: h_jet_nhm[k].Write()
    if k in h_jet_phm: h_jet_phm[k].Write()
    if k in h_jet_mum: h_jet_mum[k].Write()
    if k in h_jet_elm: h_jet_elm[k].Write()
    if k in h_jet_beta: h_jet_beta[k].Write()
    if k in h_dR2Mean: h_dR2Mean[k].Write()
    if k in h_frac01: h_frac01[k].Write()
    if k in h_frac02: h_frac02[k].Write()
    if k in h_frac03: h_frac03[k].Write()
    if k in h_frac04: h_frac04[k].Write()
    if k in h_jetRchg: h_jetRchg[k].Write()
    if k in h_jetR: h_jetR[k].Write()
    if k in h_majW: h_majW[k].Write()
    if k in h_minW: h_minW[k].Write()
    if k in h_nCharged: h_nCharged[k].Write()
    if k in h_nConstituents: h_nConstituents[k].Write()
    if k in h_ptD: h_ptD[k].Write()
    if k in h_pull: h_pull[k].Write()
    if k in h_fixedGridRhoFastjetAll: h_fixedGridRhoFastjetAll[k].Write()
#    if k in h_old_mva: h_old_mva[k].Write()
    if k in h_new_mva: h_new_mva[k].Write()

    return 0


# pileup weight

#f_pu_weights = ROOT.TFile.Open("pu_weights_${year}.root","READ")
#print("opened pu_weights_${year}.root")
#weights_tree = f_pu_weights.Get("Events")

def process_event(entry,e):
    
    weight = 1.0
    if data_type == "mc":
        weight = weight * xs_weight

    gen_weight = 1.0
    pu_weight = 1.0

    if data_type == "mc":

        if e.Generator_weight < 0:
            gen_weight = -1.0
        pu_weight = e.weights

        
    weight = weight * gen_weight * pu_weight
    #number of good primary vertex 
    nvtx = e.PV_npvsGood
    #select two same flavour leptons and check that mass around Z Boson mass to select DY leptons (Muon channel only for now)

    Lep1=ROOT.TLorentzVector()
    Lep2=ROOT.TLorentzVector()
    mu_Lep1 = ROOT.TLorentzVector()
    mu_Lep2 = ROOT.TLorentzVector()
    mu_Z_Boson = mu_Lep1+mu_Lep2
    el_Lep1 = ROOT.TLorentzVector()
    el_Lep2 = ROOT.TLorentzVector()
    el_Z_Boson = el_Lep1+el_Lep2
    Z_Boson = ROOT.TLorentzVector()
    if e.nMuon==2:
        mu_Lep1.SetPtEtaPhiM(e.Muon_pt[0], e.Muon_eta[0],e.Muon_phi[0],e.Muon_mass[0])
        mu_Lep2.SetPtEtaPhiM(e.Muon_pt[1], e.Muon_eta[1],e.Muon_phi[1],e.Muon_mass[1])
        mu_Z_Boson = mu_Lep1+mu_Lep2
    else:
        return

#    if e.nElectron==2:
#        el_Lep1.SetPtEtaPhiM(e.Electron_pt[0], e.Electron_eta[0],e.Electron_phi[0],e.Electron_mass[0])
#        el_Lep2.SetPtEtaPhiM(e.Electron_pt[1], e.Electron_eta[1],e.Electron_phi[1],e.Electron_mass[1])
#        el_Z_Boson = el_Lep1+el_Lep2

    if mu_Z_Boson.M() > 70 and mu_Z_Boson.M() < 110:
         Lep1 = mu_Lep1
         Lep2 = mu_Lep2
         Z_Boson = mu_Z_Boson

    else:
        return

#    elif el_Z_Boson.M() > 70 and el_Z_Boson.M() < 110:
#         Lep1 = el_Lep1
#         Lep2 = el_Lep2
#         Z_Boson = el_Z_Boson

            

    lept_pt1 = Lep1.Pt()
    lept_eta1 = Lep1.Eta()
    lept_phi1 = Lep1.Phi()

    lept_pt2 = Lep2.Pt()
    lept_eta2 = Lep2.Eta()
    lept_phi2 = Lep2.Phi()


    z_mass = Z_Boson.M()
    z_pt = Z_Boson.Pt()
    z_phi = Z_Boson.Phi()
    z_rapidity = Z_Boson.Eta()


    #Phi difference between Z Boson and leading jet
    #Default for Puppi Jet, needs to be adjusted to Jet instead of JetPuppi for CHS jet
    n_jets = len(e.JetPuppi_pt)
    assert (n_jets == e.nJetPuppi)

    def fill_hist_per_event(k=""):
        h_nvtx[k].Fill(nvtx, weight)

        h_z_mass[k].Fill(z_mass, weight)
        h_z_pt[k].Fill(z_pt, weight)
        h_z_phi[k].Fill(z_phi, weight) 
        h_z_rapidity[k].Fill(z_rapidity, weight)

        h_lept_pt1[k].Fill(lept_pt1, weight)
        h_lept_pt2[k].Fill(lept_pt2, weight)
        h_lept_eta1[k].Fill(lept_eta1, weight)
        h_lept_eta2[k].Fill(lept_eta2, weight)
        h_lept_phi1[k].Fill(lept_phi1, weight)
        h_lept_phi2[k].Fill(lept_phi2, weight)
        if dphi_zj != "something":
            h_dphi_zj[k].Fill(dphi_zj, weight)
        h_n_jets[k].Fill(n_jets, weight)
        h_fixedGridRhoFastjetAll[k].Fill(tmva_v_fixedGridRhoFastjetAll[0], weight)

        return 0
    lead_jet_index=0
    if n_jets > 0:
        #identify index of jet for each photons and electrons in jet (-1 if it is not in a jet)
        photon_in_jet_index=[]
        electron_in_jet_index=[]

        if e.nPhoton != 0:    
            for j in range(e.nPhoton):
                photon_in_jet_index.append(j)

        if e.nElectron != 0:
            for k in range(e.nElectron):
                electron_in_jet_index.append(k)
        #identify leading jet index
        for i in range(n_jets):
            if e.JetPuppi_pt[i]>e.JetPuppi_pt[lead_jet_index]:
                lead_jet_index = i
        #dphi of leading jet and Z boson
        dphi_zj = "something"
        dphi_zj = dphi(z_phi, e.JetPuppi_phi[lead_jet_index])
    

        e_new_mva = []


    #Default for Puppi Jet, needs to b e adjusted to Jet instead of JetPuppi for CHS jet
        for i in range(e.nJetPuppi):

            jet_pt = e.JetPuppi_pt[i]
            jet_eta = e.JetPuppi_eta[i]
            jet_phi = e.JetPuppi_phi[i]
            jet_mass = e.JetPuppi_mass[i]
            temp_jet=ROOT.TLorentzVector()
            temp_jet.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass)
            jet_energy = temp_jet.E()

            jet_chf = e.JetPuppi_chHEF[i]
            jet_cemf = e.JetPuppi_chEmEF[i]
            jet_nemf = e.JetPuppi_neEmEF[i]
            jet_nhf = e.JetPuppi_neHEF[i]

        #select photons and electrons associated with the current jet to calculate the energy fraction.
            photon_in_jet=ROOT.TLorentzVector()
            if e.nPhoton !=0 and e.JetPuppi_nConstPhotons != 0:
                for j in range(e.nPhoton):
                    temp_photon_in_jet=ROOT.TLorentzVector()
                    if photon_in_jet_index[j]==i:
                        temp_photon_in_jet.SetPtEtaPhiM(e.Photon_pt[j], e.Photon_eta[j], e.Photon_phi[j], e.Photon_mass[j])
                        photon_in_jet = photon_in_jet+temp_photon_in_jet
                jet_phf = photon_in_jet.E()/jet_energy

            electron_in_jet=ROOT.TLorentzVector()
            if e.nElectron != 0 and e.JetPuppi_nElectrons != 0: 
                for j in range(e.nElectron):
                    temp_electron_in_jet=ROOT.TLorentzVector()
                    if electron_in_jet_index[j]==i:
                        temp_electron_in_jet.SetPtEtaPhiM(e.Electron_pt[j], e.Electron_eta[j], e.Electron_phi[j], e.Electron_mass[j])
                        electron_in_jet = electron_in_jet+temp_electron_in_jet
                jet_elf = electron_in_jet.E()/jet_energy

            jet_muf = e.JetPuppi_muEF[i]

            jet_chm = e.JetPuppi_nConstChHads[i]
            jet_nhm = e.JetPuppi_nConstNeuHads[i]
            jet_phm = e.JetPuppi_nConstPhotons[i]
            jet_mum = e.JetPuppi_nConstMuons[i]
            jet_elm = e.JetPuppi_nConstElecs[i]
            jet_npr = ord(e.JetPuppi_nConstituents[i]) - int(e.JetPuppi_puId_nCharged[i])

            if data_type == "mc":
                jet_dR = -1
                Gen_index = event.JetPuppi_genJetIdx[i]
                if Gen_index != -1:
                    Matched_GenJet_eta = event.GenJet_eta[Gen_index]
                    Matched_GenJet_phi = event.GenJet_phi[Gen_index]
                    jet_dR = ((jet_eta-Matched_GenJet_eta)**2+(jet_phi-Matched_GenJet_phi)**2)**0.5
                jet_flavor = e.JetPuppi_partonFlavour[i]

            tmva_s_jetPt[0] = jet_pt
            tmva_s_jetEta[0] = jet_eta

            tmva_v_nvtx[0] = nvtx
            tmva_v_beta[0] = e.JetPuppi_puId_beta[i]
            tmva_v_dR2Mean[0] = e.JetPuppi_puId_dR2Mean[i]
            tmva_v_frac01[0] = e.JetPuppi_puId_frac01[i]
            tmva_v_frac02[0] = e.JetPuppi_puId_frac02[i]
            tmva_v_frac03[0] = e.JetPuppi_puId_frac03[i]
            tmva_v_frac04[0] = e.JetPuppi_puId_frac04[i]
            tmva_v_majW[0] = e.JetPuppi_puId_majW[i]
            tmva_v_minW[0] = e.JetPuppi_puId_minW[i]
            tmva_v_jetR[0] = e.JetPuppi_puId_jetR[i]
            tmva_v_jetRchg[0] = e.JetPuppi_puId_jetRchg[i]
            tmva_v_nConstituents[0] = ord(e.JetPuppi_nConstituents[i])
            tmva_v_nCharged[0] = e.JetPuppi_puId_nCharged[i]
            tmva_v_ptD[0] = e.JetPuppi_puId_ptD[i]
            tmva_v_pull[0] = e.JetPuppi_puId_pull[i]         
            tmva_v_fixedGridRhoFastjetAll[0] = e.fixedGridRhoFastjetAll
#            old_mva = e.pumva[i]
            new_mva = 0

            for i, f_eta in enumerate(f_eta_bins):
                #this -1 is added to avoid the new_mva score array not filling up for jet that does not belong in any of the eta bins
                new_mva=-1
                if abs(jet_eta) > f_eta[0] and abs(jet_eta) <= f_eta[1]: 
                    new_mva = tmva_readers[i].EvaluateMVA("BDT")
                e_new_mva.append(new_mva)

            def fill_hist_per_jet(k=""):
                h_jet_pt[k].Fill(jet_pt, weight)
                h_jet_eta[k].Fill(jet_eta, weight)
                h_jet_phi[k].Fill(jet_phi, weight)
                h_jet_energy[k].Fill(jet_energy, weight)
 
                if data_type == "mc": h_jet_dR[k].Fill(jet_dR, weight)
 
                h_jet_chf[k].Fill(jet_chf, weight)
                h_jet_cemf[k].Fill(jet_cemf, weight)
                h_jet_nhf[k].Fill(jet_nhf, weight)
                h_jet_nemf[k].Fill(jet_nemf, weight)
                if e.nPhoton !=0 and e.JetPuppi_nConstPhotons != 0:
                    h_jet_phf[k].Fill(jet_phf, weight)
                h_jet_muf[k].Fill(jet_muf, weight)
                if e.nElectron != 0 and e.JetPuppi_nElectrons != 0: 
                    h_jet_elf[k].Fill(jet_elf, weight)
                h_jet_npr[k].Fill(jet_npr, weight)
                h_jet_chm[k].Fill(jet_chm, weight)
                h_jet_nhm[k].Fill(jet_nhm, weight)
                h_jet_phm[k].Fill(jet_phm, weight)
                h_jet_mum[k].Fill(jet_mum, weight)
                h_jet_elm[k].Fill(jet_elm, weight)
 
                h_jet_beta[k].Fill(tmva_v_beta[0], weight)
                h_dR2Mean[k].Fill(tmva_v_dR2Mean[0], weight)
                h_frac01[k].Fill(tmva_v_frac01[0], weight)
                h_frac02[k].Fill(tmva_v_frac02[0], weight)
                h_frac03[k].Fill(tmva_v_frac03[0], weight)
                h_frac04[k].Fill(tmva_v_frac04[0], weight)
                h_jetRchg[k].Fill(tmva_v_jetRchg[0], weight)
                h_jetR[k].Fill(tmva_v_jetR[0], weight)
                h_majW[k].Fill(tmva_v_majW[0], weight)
                h_minW[k].Fill(tmva_v_minW[0], weight)
                h_nCharged[k].Fill(tmva_v_nCharged[0], weight)
                h_nConstituents[k].Fill(tmva_v_nConstituents[0], weight)
                h_ptD[k].Fill(tmva_v_ptD[0], weight)
                h_pull[k].Fill(tmva_v_pull[0], weight)

#                h_old_mva[k].Fill(old_mva, weight)
                h_new_mva[k].Fill(new_mva, weight)

                return 0
 
            fill_hist_per_jet(k=h_keys1[0])

            if data_type =="mc":
                if jet_dR < 0.2 and jet_dR >= 0:
                    jet_type = "prompt"

                elif jet_dR > 0.4 and jet_dR >= 0 and abs(jet_flavor) == 0:
                   jet_type = "pileup"

                else:
                   jet_type = "rest"

            for i, f_eta in enumerate(f_eta_bins):
    
                if abs(jet_eta) > f_eta[0] and abs(jet_eta) <= f_eta[1]:
                    key_eta = eta_bins[i]
    
                    if data_type == "mc":
                        k1 = "%s_%s_%s" % (h_keys1[0], key_eta, "all")
                        k2 = "%s_%s_%s" % (h_keys1[0], key_eta, jet_type)
                        fill_hist_per_jet(k1)
                        fill_hist_per_jet(k2)
 
                    if data_type == "data":
                        k1 = "%s_%s" % (h_keys1[0], key_eta)
                        fill_hist_per_jet(k1)
       

                    for j, f_pt in enumerate(f_pt_bins):
 
                        if jet_pt > f_pt[0] and jet_pt <= f_pt[1]:
                            key_pt = pt_bins[j] 
                            if data_type == "mc":
                                k1 = "%s_%s_%s_%s" % (h_keys1[0], key_eta, key_pt, "all")
                                k2 = "%s_%s_%s_%s" % (h_keys1[0], key_eta, key_pt, jet_type)
                                fill_hist_per_jet(k1)
                                fill_hist_per_jet(k2)
 
                            if data_type == "data":
                                k1 = "%s_%s_%s" % (h_keys1[0], key_eta, key_pt)
                                fill_hist_per_jet(k1)
     
  
        lead_jet_pt=0
        lead_jet_index=0
        for n in range(0,e.nJetPuppi):
            if lead_jet_pt < e.JetPuppi_pt[n]:
                lead_jet_pt = e.JetPuppi_pt[n]
                lead_jet_index = n
        if data_type == "mc" and len(e.JetPuppi_genJetIdx)>0:
            lead_jet_dR = -1
            lead_jet_genJet_index = e.JetPuppi_genJetIdx[lead_jet_index]
            if lead_jet_genJet_index != -1:
                Matched_GenJet_eta = e.GenJet_eta[lead_jet_genJet_index]
                Matched_GenJet_phi = e.GenJet_phi[lead_jet_genJet_index]
                lead_jet_dR = ((e.JetPuppi_eta[lead_jet_index]-Matched_GenJet_eta)**2+(e.JetPuppi_phi[lead_jet_index]-Matched_GenJet_phi)**2)**0.5
            if lead_jet_dR < 0.2 and lead_jet_dR >= 0 and z_pt>0:
                 h_dphi_zj_ptj_z["prompt"].Fill(dphi_zj, float(lead_jet_pt/z_pt), weight)

            elif lead_jet_dR > 0.4 and lead_jet_dR >= 0 and z_pt>0 and abs(e.JetPuppi_hadronFlavour[lead_jet_index]) == 0:
                 h_dphi_zj_ptj_z["pileup"].Fill(dphi_zj, float(lead_jet_pt/z_pt), weight)
    
        elif data_type == "data" and z_pt>0:
            h_dphi_zj_ptj_z["data"].Fill(dphi_zj, float(lead_jet_pt/z_pt), weight)

        fill_hist_per_event(k=h_keys1[0])

        for i in range(n_jets):
           
            pt_zj_ratio = 0
            if z_pt == 0: continue
            else:
                pt_zj_ratio = lead_jet_pt/z_pt

            if (abs(dphi_zj) > 2.5) and (pt_zj_ratio > 0.5) and (pt_zj_ratio < 1.5):

                key = "tight"
                if len(e_new_mva)>0:
                    if float(e_new_mva[i]) > mva_wp1:
                        h_control_jet_pt[f"{key}_pass_wp1"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_pass_wp1"].Fill(e.JetPuppi_eta[i], weight)
  
                    else:
                        h_control_jet_pt[f"{key}_fail_wp1"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_fail_wp1"].Fill(e.JetPuppi_eta[i], weight)
  
                    if float(e_new_mva[i]) > mva_wp2:
                        h_control_jet_pt[f"{key}_pass_wp2"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_pass_wp2"].Fill(e.JetPuppi_eta[i], weight)
  
                    else:
                        h_control_jet_pt[f"{key}_fail_wp2"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_fail_wp2"].Fill(e.JetPuppi_eta[i], weight)
  
                    if float(e_new_mva[i]) > mva_wp3:
                        h_control_jet_pt[f"{key}_pass_wp3"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_pass_wp3"].Fill(e.JetPuppi_eta[i], weight)
  
                    else:
                        h_control_jet_pt[f"{key}_fail_wp3"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_fail_wp3"].Fill(e.JetPuppi_eta[i], weight)
 
            if abs(dphi_zj) < 1.5:

                key = "loose"
                if len(e_new_mva)>0:
 
                    if float(e_new_mva[i]) > mva_wp1:
                        h_control_jet_pt[f"{key}_pass_wp1"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_pass_wp1"].Fill(e.JetPuppi_eta[i], weight)
 
                    else:
                        h_control_jet_pt[f"{key}_fail_wp1"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_fail_wp1"].Fill(e.JetPuppi_eta[i], weight)

                    if float(e_new_mva[i]) > mva_wp2:
                        h_control_jet_pt[f"{key}_pass_wp2"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_pass_wp2"].Fill(e.JetPuppi_eta[i], weight)
 
                    else:
                        h_control_jet_pt[f"{key}_fail_wp2"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_fail_wp2"].Fill(e.JetPuppi_eta[i], weight)
 
                    if float(e_new_mva[i]) > mva_wp3:
                        h_control_jet_pt[f"{key}_pass_wp3"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_pass_wp3"].Fill(e.JetPuppi_eta[i], weight)
 
                    else:
                        h_control_jet_pt[f"{key}_fail_wp3"].Fill(e.JetPuppi_pt[i], weight)
                        h_control_jet_eta[f"{key}_fail_wp3"].Fill(e.JetPuppi_eta[i], weight)


total_events = events.GetEntries()

mva_wp1 = -0.21
mva_wp2 = 0.65
mva_wp3 = 0.98

for i, event in enumerate(events):

    if i%1000 == 0:
        print("processing event:   %s/%s" % (i, total_events))

    process_event(i, event)
    
# output file for histograms
output_file = ROOT.TFile(output_filename, "RECREATE")

write_hists(h_keys1[0])

for k in h_keys4:
    write_hists(h_keys1[0] + "_" + k)

for k in ["data", "prompt", "pileup"]:
    if k in h_dphi_zj_ptj_z:
        h_dphi_zj_ptj_z[k].Write()
        
for k in ["tight", "loose"]:
    
    for j in control_wp_keys:
        k_ = k + "_" + j
        
        if k_ in h_control_jet_pt:
            h_control_jet_pt[k_].Write()
        if k_ in h_control_jet_eta:
            h_control_jet_eta[k_].Write()


output_file.cd()
output_file.Write()
output_file.Close()
