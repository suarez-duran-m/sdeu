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

pmtId = "PMT" + sys.argv[1]
pmt = int(sys.argv[1])

nameUub = "../mltpst/uubCalibHist"+pmtId+".root"
nameUb = "../../nouub/fitHist/mltpst/ubCalibHist"+pmtId+".root"

inFileUub = ROOT.TFile.Open(nameUub, "READ")
#inFileUb = ROOT.TFile.Open(nameUb, "READ")

print("Reading files:\n", nameUub, "\nand\n", nameUb)

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

histosUub = inFileUub.Get("Histograms")
#histosUb = inFileUb.Get("Histograms")

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

extrUub = extrtor(histosUub)

extrUub.getAveOffset(0, 1)

print(extrUub.aveStOffset)

'''
data = {}
data['uub'] = []
data['ub'] = []

for i in range(0, 2):
    data['uub'].append({
    'station': i,
    'aveOffset': aveStatUub[i]
    })
    data['ub'].append({
    'station': i,
    'aveOffset': aveStatUb[i]
    })

json_object = json.dumps(data, indent = 3)
with open("sample.json", "w") as outfile:
    outfile.write(json_object)
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

