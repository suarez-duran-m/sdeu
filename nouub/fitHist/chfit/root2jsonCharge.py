import ROOT
import sys
import json
import numpy as np

# ====================
# *** *** MAIN *** ***
# ====================

elecver = sys.argv[1]
pmtid = sys.argv[2]
month = sys.argv[3]
lstSt = sys.argv[4]

stList = np.loadtxt(lstSt, int)

basename= "smooth"+elecver+"ChargePMT"+pmtid+month
nameElec = basename+".root"
outname = basename+".json"

data = {}
for st in stList:
    data[str(st)] = []

inFile = ROOT.TFile.Open(nameElec, "READ")
print("\nDoing for file: ", nameElec)

histos = inFile.Get("Histograms")

for day in range(0, histos.GetEntries()):
    histos.GetEntry(day)
    for st in range(len(stList)):
        data[str(stList[st])].append({
            'Date' : histos.timestmp.GetBinContent(st),
            'TotEvt' : histos.eventStat.GetBinContent(st),
            'VemChHb' : histos.vemchHb.GetBinContent(st),
            'CntVemChHb' : histos.cntvemchHb.GetBinContent(st),
            'ChiHb' : histos.chiHb.GetBinContent(st),
            'VemChrmsHb' : histos.vemchHbrms.GetBinContent(st),
            'CntVemChrmsHb' : histos.cntvemchHbrms.GetBinContent(st),
            'ChirmsHb' : histos.chiHbrms.GetBinContent(st)
            })

json_object = json.dumps(data, indent = 2) 
with open(outname, "w") as outfile:
    outfile.write(json_object)
