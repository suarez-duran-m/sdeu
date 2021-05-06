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

whichElec = 'uub'
pmtId = "PMT" + sys.argv[1]
pmt = int(sys.argv[1])
stLbl = sys.argv[2]
stId = int(stLbl)

nameElec = "kkuubAoPsgl"+pmtId+"St"+stLbl+".root"
outname = "kkuubAoPsgl"+pmtId+"St"+stLbl+".json"

data = {}
data[whichElec] = []

inFile = ROOT.TFile.Open(nameElec, "READ")
print("\nDoing for file: ", nameElec)
histos = inFile.Get("Histograms")
    
extPerStat = extrtor(histos)

extPerStat.getAoPDistHb()
extPerStat.getAoPDistCa()
extPerStat.getAoPAveRmsHb()
extPerStat.getAoPAveRmsCa()

data[whichElec].append({
    'stationId': stId,
    'extAveAoPHb': np.average(extPerStat.aveaopHb),
    'extRmsAoPHb': np.average(extPerStat.rmsaopHb),
    'extAveAoPCa': np.average(extPerStat.aveaopCa),
    'extRmsAoPCa': np.average(extPerStat.rmsaopCa),
    'extFitOkAoPHb': extPerStat.getFitokHb(),
    'extFitOkAoPCa': extPerStat.getFitokCa(),
    'extTotEvt': extPerStat.getTotEvents()
    })

json_object = json.dumps(data, indent = 2) 
with open(outname, "w") as outfile:
    outfile.write(json_object)
