import ROOT
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from extractor import extractor as extrtor

def getAveFstBinCnt(hist):
    lstentry = 32000 #hist.GetEntries()
    totetrySt = 0
    averageStat = []
    tmpave = 0
    tmpentry = 0

    for st in range(0, 19):
        tmpave = 0
        totetrySt = 0
        if st == 16:
            averageStat.append(0)
            continue
        for evt in range(0, lstentry):
            hist.GetEntry(evt)
            if tmpentry != hist.eventStat[st+1]:
                tmpentry = hist.eventStat[st+1]
                totetrySt = tmpentry
                tmpave += hist.firstBinCntPk.GetBinContent(st)
        averageStat.append( tmpave / totetrySt )
    return averageStat

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
    outname = "uubAllAveOffset"+pmtId+".json"
else:
    nameElec = "../../nouub/fitHist/mltpst/ubCalibHist"+pmtId+".root"
    inFile = ROOT.TFile.Open(nameElec, "READ")
    outname = "ubAllAveOffset"+pmtId+".json"

print("Reading files for:\n", nameElec)

stLabel = []
histos = inFile.Get("Histograms")

extAveOffset = extrtor(histos)
extAveOffset.getOffsetDist(stLabel)
extAveOffset.getAveOffsetAveRms()

data = {}
data[whichElec] = []

for i in range(0, len(stLabel)):
    data[whichElec].append({
        'stationId': 
        'extAveOffset': extAveOffset.aveStOffset[i],
        'extRmsOffset': extAveOffset.rmsStOffset[i]
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

# =======================
# *** For Peak Offset ***

# *** Average ***


'''
aveStatUub = getAveOffset(histosUub)
print("Done for UUB")
aveStatUb = getAveOffset(histosUb)
print("Done for UB")
'''

'''
plt.figure(figsize=(16,9))
plt.plot(aveStatUub, 'o', markersize=12, label="UUB")
plt.plot(aveStatUb, 'x', markersize=12, label="UB")
plt.legend()
plt.title("Average of Offset for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Offset (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(50,)
#plt.savefig("../../plots/bothAllAveOffset"+pmtId+".png", dpi=100)
plt.clf()
print("bothAllAveOffset OK")
'''

