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
stLabel = np.loadtxt(sys.argv[3],  int)

nameElec = "ubFitaop"+pmtId
outname = "ubFitaop"+pmtId+".json"

months = ['aug', 'sep', 'oct']
stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

data = {}
data[whichElec] = []
lenStati = len(stLabel)
lenMonth = len(months)
tmpaveAoPHb = np.zeros( ( lenStati, lenMonth ) )
tmprmsAoPHb = np.zeros( ( lenStati, lenMonth ) )
tmpaveAoPCa = np.zeros( ( lenStati, lenMonth ) )
tmprmsAoPCa = np.zeros( ( lenStati, lenMonth ) )

for mth in range(0, lenMonth):
    nameElec = "ubFitaop"+pmtId+months[mth]+".root"
    inFile = ROOT.TFile.Open(nameElec, "READ")
    print("\nDoing for file: ", nameElec)
    histos = inFile.Get("Histograms")
    
    extPerStat = extrtor(histos)
    extPerStat.getAoPDistHb(stLabel)
    extPerStat.getAoPDistCa(stLabel)
    extPerStat.getAoPAveRmsHb()
    extPerStat.getAoPAveRmsCa()

    for st in range(0, lenStati):
        tmpaveAoPHb[st][mth] = extPerStat.aveaopHb[st]
        tmprmsAoPHb[st][mth] = extPerStat.rmsaopHb[st]
        tmpaveAoPCa[st][mth] = extPerStat.aveaopCa[st]
        tmprmsAoPCa[st][mth] = extPerStat.rmsaopCa[st]


for st in range(0, lenStati):
    data[whichElec].append({
        'stationId': stLabel[st],
        'extAveAoPHb': np.average(tmpaveAoPHb[st]),
        'extRmsAoPHb': np.average(tmprmsAoPHb[st]),
        'extAveAoPCa': np.average(tmpaveAoPCa[st]),
        'extRmsAoPCa': np.average(tmprmsAoPHb[st])
        })

json_object = json.dumps(data, indent = 2) 
with open(outname, "w") as outfile:
    outfile.write(json_object)
