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

nameElec = "ubFitXi"+pmtId
outname = "kkubFitXi"+pmtId+".json"

months = ['aug', 'sep', 'oct']
stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

data = {}
data[whichElec] = []
lenStati = len(stLabel)
lenMonth = len(months)
tmpaveXiPkHb = np.zeros( ( lenStati, lenMonth ) )
tmprmsXiPkHb = np.zeros( ( lenStati, lenMonth ) )
tmpaveXiPkCa = np.zeros( ( lenStati, lenMonth ) )
tmprmsXiPkCa = np.zeros( ( lenStati, lenMonth ) )

tmpaveXiChHb = np.zeros( ( lenStati, lenMonth ) )
tmprmsXiChHb = np.zeros( ( lenStati, lenMonth ) )
tmpaveXiChCa = np.zeros( ( lenStati, lenMonth ) )
tmprmsXiChCa = np.zeros( ( lenStati, lenMonth ) )

totXiPkHb = np.zeros( lenStati )
totXiPkCa = np.zeros( lenStati )
totXiChHb = np.zeros( lenStati )
totXiChCa = np.zeros( lenStati )

totEvt = np.zeros( lenStati )

for mth in range(0, lenMonth):
    nameElec = "ubFitaop"+pmtId+months[mth]+".root"
    inFile = ROOT.TFile.Open(nameElec, "READ")
    print("\nDoing for file:", nameElec)
    histos = inFile.Get("Histograms")

    extPerStat = extrtor(histos)
    extPerStat.getTotEvents(stLabel)
    
    extPerStat.getXiPkDistHb(stLabel)
    extPerStat.getXiPkveRmsHb()
    extPerStat.getXiPkDistCa(stLabel)
    extPerStat.getXiPkveRmsCa()

    extPerStat.getXiChDistHb(stLabel)
    extPerStat.getXiChveRmsHb()
    extPerStat.getXiChDistCa(stLabel)
    extPerStat.getXiChveRmsCa()

    for st in range(0, lenStati):
        tmpaveXiPkHb[st][mth] = extPerStat.aveXiPkHb[st]
        tmprmsXiPkHb[st][mth] = extPerStat.rmsXiPkHb[st]
        tmpaveXiPkCa[st][mth] = extPerStat.aveXiPkCa[st]
        tmprmsXiPkCa[st][mth] = extPerStat.rmsXiPkCa[st]

        tmpaveXiChHb[st][mth] = extPerStat.aveXiChHb[st]
        tmprmsXiChHb[st][mth] = extPerStat.rmsXiChHb[st]
        tmpaveXiChCa[st][mth] = extPerStat.aveXiChCa[st]
        tmprmsXiChCa[st][mth] = extPerStat.rmsXiChCa[st]

        if type(extPerStat.distXiPkHb[st]) is not int:
            totXiPkHb[st] += len(extPerStat.distXiPkHb[st])
            totXiPkCa[st] += len(extPerStat.distXiPkCa[st])
            totXiChHb[st] += len(extPerStat.distXiChHb[st])
            totXiChCa[st] += len(extPerStat.distXiChCa[st])
        totEvt[st] += extPerStat.totEvts[st]

for st in range(0, lenStati):
    data[whichElec].append({
        'stationId': stLabel[st],
        'extAveXiPkHb': np.average(tmpaveXiPkHb[st]),
        'extRmsXiPkHb': np.average(tmprmsXiPkHb[st]),
        'extAveXiPkCa': np.average(tmpaveXiPkCa[st]),
        'extRmsXiPkCa': np.average(tmprmsXiPkCa[st]),
        'extAveXiChHb': np.average(tmpaveXiChHb[st]),
        'extRmsXiChHb': np.average(tmprmsXiChHb[st]),
        'extAveXiChCa': np.average(tmpaveXiChCa[st]),
        'extRmsXiChCa': np.average(tmprmsXiChCa[st]),
        'extTotXiPkHb': totXiPkHb[st],
        'extTotXiPkCa': totXiPkCa[st],
        'extTotXiChHb': totXiChHb[st],
        'extTotXiChCa': totXiChCa[st],
        'extTotEvet': totEvt[st]
        })

json_object = json.dumps(data, indent = 2) 
with open(outname, "w") as outfile:
    outfile.write(json_object)
