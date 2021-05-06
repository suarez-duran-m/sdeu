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
    outname = "uubAllAveHbaseCalib"+pmtId+".json"
else:
    nameElec = "../../nouub/fitHist/mltpst/ubCalibHist"+pmtId+".root"
    inFile = ROOT.TFile.Open(nameElec, "READ")
    outname = "ubAllAveHbaseCalib"+pmtId+".json"

print("Reading files for:\n", nameElec)

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

histos = inFile.Get("Histograms")

extPerStat = extrtor(histos)

extPerStat.getHbaseDist(stLabel)
extPerStat.getCalibDist(stLabel)
extPerStat.getHbaseAveRms()
extPerStat.getCalibAveRms()
extPerStat.getBinCenterDistHbase(stLabel)
extPerStat.getBinCenterDistCalib(stLabel)
extPerStat.getBinCenterAveRmsHbase()
extPerStat.getBinCenterAveRmsCalib()

data = {}
data[whichElec] = []

for i in range(0, len(stLabel)):
    data[whichElec].append({
        'stationId': stLabel[i],
        'extAveHbase': extPerStat.aveHbase[i],
        'extRmsHbase': extPerStat.rmsHbase[i],
        'extAveCalib': extPerStat.aveCalib[i],
        'extRmsCalib': extPerStat.rmsCalib[i],
        'extAveFrstBinCntrHbase': extPerStat.aveBinCntrHbase[i],
        'extAveFrstBinCntrCalib': extPerStat.aveBinCntrCalib[i],
        'extRmsFrstBinCntrHbase': extPerStat.rmsBinCntrHbase[i],
        'extRmsFrstBinCntrCalib': extPerStat.rmsBinCntrCalib[i]
        })

json_object = json.dumps(data, indent = 2) 
with open(outname, "w") as outfile:
    outfile.write(json_object)

