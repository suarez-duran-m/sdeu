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
    outname = "uubAllAveFstBinOffset"+pmtId+".json"
else:
    nameElec = "../../nouub/fitHist/mltpst/ubCalibHist"+pmtId+".root"
    inFile = ROOT.TFile.Open(nameElec, "READ")
    outname = "ubAllAveFstBinOffset"+pmtId+".json"

print("Reading files for:\n", nameElec)

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

histos = inFile.Get("Histograms")

extPerStat = extrtor(histos)

extPerStat.getOffsetDist(stLabel)
extPerStat.getAveOffsetAveRms()

extPerStat.getFrstBinDist.(stLabel)
extPerStat.getFrstBinAveRms()

data = {}
data[whichElec] = []

for i in range(0, len(stLabel)):
    data[whichElec].append({
        'stationId': stLabel[i],
        'extAveOffset': extPerStat.aveStOffset[i],
        'extRmsOffset': extPerStat.rmsStOffset[i],
        'extAveFrstBin': extPerStat.aveFrstBin[i],
        'extRmsFrstBin': extPerStat.rmsFrstBin[i]
        })

json_object = json.dumps(data, indent = 4) 
#outputName = "ubAveOffset"+sys.argv[2]+str(runSt1-1)+".json"
with open(outname, "w") as outfile:
    outfile.write(json_object)


# =======================================
# *** For counts in firt bin for Peak ***

# *** Average ***
'''
aveStatUub = getAveFstBinCnt(histosUub)
print("Done for UUB")
aveStatUb = getAveFstBinCnt(histosUb)
print("Done for UB")

plt.figure(figsize=(16,9))
plt.plot(aveStatUub, 'o', markersize=12, label="UUB")
plt.plot(aveStatUb, 'x', markersize=12, label="UB")
plt.legend()
plt.title("Average of counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation=45)
plt.ylabel("Average Counts (au)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../plots/bothAllAveFrstBinCnt"+pmtId+".png", dpi=100)
plt.show()
plt.clf()
print("bothAllAveFrstBinCnt OK")
'''
