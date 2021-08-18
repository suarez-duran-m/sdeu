import ROOT
import sys
import json
import numpy as np


# ====================
# =======      =======
# *** *** MAIN *** ***
# ====================

elecver = sys.argv[1]
pmtid = sys.argv[2]
stid = sys.argv[3]
month = sys.argv[4]

basename= elecver+"AoPtimePMT"+pmtid+"St"+stid+"Mth"+month+"chpk"
nameElec = basename+".root"
outname = basename+".json"

data = {}
data[month] = []
    
inFile = ROOT.TFile.Open(nameElec, "READ")
print("\nDoing for file: ", nameElec)
histos = inFile.Get("Histograms")

for day in range(0, inFile.apHbase.GetXaxis().GetNbins()):
    data[month].append({
        'Date' : inFile.timestmp.GetBinContent(day),
        'AoP' : inFile.apHbase.GetBinContent(day),
        'rms' : inFile.rmsHbase.GetBinContent(day),
        'pk' : inFile.pkHbase.GetBinContent(day),
        'pkrms' : inFile.rmspkHbase.GetBinContent(day),
        'ch' : inFile.chHbase.GetBinContent(day),
        'chrms' : inFile.rmschHbase.GetBinContent(day)
        })

json_object = json.dumps(data, indent = 2) 
with open(outname, "w") as outfile:
    outfile.write(json_object)
