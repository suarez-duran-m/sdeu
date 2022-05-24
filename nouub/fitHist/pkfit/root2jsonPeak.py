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

basename= elecver+"PeakPMT"+pmtid+month
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
            'VemPkHb' : histos.vempkHb.GetBinContent(st),
            'CntVemPkHb' : histos.cntvempkHb.GetBinContent(st),
            'ChiHb' : histos.chiHb.GetBinContent(st),
            'VemPkrmsHb' : histos.vempkHbrms.GetBinContent(st),
            'CntVemPkrmsHb' : histos.cntvempkHbrms.GetBinContent(st),
            'ChirmsHb' : histos.chiHbrms.GetBinContent(st),
            'VemPkCa' : histos.vempkCa.GetBinContent(st),
            'CntVemPkCa' : histos.cntvempkCa.GetBinContent(st),
            'ChiCa' : histos.chiCa.GetBinContent(st),
            'VemPkrmsCa' : histos.vempkCarms.GetBinContent(st),
            'CntVemPkrmsCa' : histos.cntvempkCarms.GetBinContent(st),
            'ChirmsCa' : histos.chiCarms.GetBinContent(st)
            })

json_object = json.dumps(data, indent = 2) 
with open(outname, "w") as outfile:
    outfile.write(json_object)
