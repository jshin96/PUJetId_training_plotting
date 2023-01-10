import ROOT
import sys

variable = sys.args[1]

input_file = ROOT.TFile.Open("","READ")
eta_bins = [
    "Eta0p0To2p5",
    "Eta2p5To2p75",
    "Eta2p75To3p0",
    "Eta3p0To5p0"
]

pileup_hist=input_file.Get("")
prompt_hist=input_file.Get("")
rest_hist=input_file.Get("")
all_hist=input_file.Get("")



