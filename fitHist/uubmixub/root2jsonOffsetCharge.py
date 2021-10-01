import ROOT
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from extractor import extractor as extrtor


# ====================
# =======      =======
# *** *** MAIN *** ***
# ====================

whichElec = sys.argv[1]
pmtId = "PMT" + sys.argv[2]
pmt = int(sys.argv[2])

nameElec = ""
stLabel = np.loadtxt(sys.argv[3],  int)

if whichElec=='uub':
    nameElec = "../mltpst/uubCalibHist"+pmtId+".root"
    inFile = ROOT.TFile.Open(nameElec, "READ")
    outname = "uubAllAveOffsetCh"+pmtId+".json"
else:
    nameElec = "../../nouub/fitHist/mltpst/ubCalibHist"+pmtId+".root"
    inFile = ROOT.TFile.Open(nameElec, "READ")
    outname = "ubAllAveOffsetCh"+pmtId+".json"

print("Reading files for:\n", nameElec)

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

histos = inFile.Get("Histograms")

extPerStat = extrtor(histos)

extPerStat.getOffsetDistCh(stLabel)
extPerStat.getAveOffsetAveRmsCh()

data = {}
data[whichElec] = []

for i in range(0, len(stLabel)):
    data[whichElec].append({
        'stationId': stLabel[i],
        'extAveOffset': extPerStat.avestoffsetCh[i],
        'extRmsOffset': extPerStat.rmsstoffsetCh[i]
        })

json_object = json.dumps(data, indent = 2)
with open(outname, "w") as outfile:
    outfile.write(json_object)



