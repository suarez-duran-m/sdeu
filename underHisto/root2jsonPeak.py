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

basename= elecver+"AoPPMT"+pmtid+"St863Mth"+month
nameElec = basename+".root"
outname = basename+".json"

data = {}
data["863"] = []

inFile = ROOT.TFile.Open(nameElec, "READ")
print("\nDoing for file: ", nameElec)

peakInfo = inFile.Get("PeakData")

tmpchi = ROOT.TH1F('tmpchi', '', 1000, 0, 100)
tmpndf = ROOT.TH1F('tmpndf', '', 1000, 0, 100)
tmpchindf = ROOT.TH1F('tmpchindf', '', 100, 0, 10)

for ntry in range(0, peakInfo.GetEntries()):
  peakInfo.GetEntry(ntry)
  data["863"].append({
    'UtcTime' : peakInfo.timeEvnt,
    'EventId' : peakInfo.eventId,
    'peakVal' : peakInfo.peakVal,
    'chi2' : peakInfo.pkChi2,
    'ndf' : peakInfo.pkNdf
    })
  tmpchi.Fill( peakInfo.pkChi2 )
  tmpndf.Fill( peakInfo.pkNdf )
  tmpchindf.Fill( peakInfo.pkChi2 / peakInfo.pkNdf )

json_object = json.dumps(data, indent = 2) 
with open(outname, "w") as outfile:
  outfile.write(json_object)

c1 = ROOT.TCanvas('c1','', 1600, 900)
tmpchindf.Draw()
c1.Print("kk.pdf")
 
