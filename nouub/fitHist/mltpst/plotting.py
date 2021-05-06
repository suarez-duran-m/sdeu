import ROOT
import sys
import datetime as dt
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from mpl_toolkits.axes_grid1 import make_axes_locatable

pmtId = "PMT" + sys.argv[1]
pmt = int(sys.argv[1])

inFile = ROOT.TFile.Open("ubCalibHist"+pmtId+"AoP.root", "READ")
#inFile = ROOT.TFile.Open("ubCalibHist"+pmtId+"Aug.root", "READ")
#inFile = ROOT.TFile.Open("ubCalibHist"+pmtId+".root", "READ")

#print("Reading file:", "ubCalibHist"+pmtId+"AoP.root")
#print("Reading file:", "ubCalibHist"+pmtId+"Aug.root")
print("Reading file:", "ubCalibHist"+pmtId+".root")

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

histos = inFile.Get("Histograms")

'''
xdatesStation = []
xtmp = []
for st in range(0, 19):
    xtmp = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xtmp.append( dt.datetime.fromtimestamp( histos.evtTime ) )
    xdatesStation.append( xtmp )
print("Dates builded")
'''

# =====================
# *** For Area/Peak ***

# =====================
# ***  For Averages ***

aopAveSt = []
aopDisSt = []
tmpaop = 0.
tmpaopDist = []
cntEvtSt = 0
crrEvtSt = 0
nevts = histos.GetEntries()

for st in range(0, 19):
    tmpaop = 0.
    cntEvtSt = 0
    tmpaopDist = []
    for evt in range(0, nevts):
        histos.GetEntry( evt )
        if crrEvtSt != histos.eventStat[st+1] and histos.apHbase.GetBinContent(st) > 0:
            crrEvtSt = histos.eventStat[st+1]
            if histos.apHbase.GetBinContent(st) > 0:
                tmpaop += histos.apHbase.GetBinContent(st)
                cntEvtSt += 1
                tmpaopDist.append( histos.apHbase.GetBinContent(st) )

    if len (tmpaopDist) > 0:
        aopAveSt.append( tmpaop / cntEvtSt )
        aopDisSt.append( tmpaopDist )
    else:
        aopAveSt.append(0)
        aopDisSt.append(0)
        print("Station", stLabel[st], "wrong", tmpaop)
    print("AoP Average done for", stLabel[st])

plt.figure(figsize=(16,9))
plt.plot(np.arange(1,20), aopAveSt, 'or')
plt.violinplot(aopDisSt, showmeans=True)
plt.title("UB Area over Peak "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xlim(0,20)
plt.ylim(2,)
plt.xticks(np.arange(1,20), stLabel, rotation=45)
plt.ylabel("AoP [25 ns]", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.show()
plt.savefig("../../../plots/ubAllAoPave"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllAoPave OK")




# *** For August ***
'''
# *** HBase ***
xbl = []
ybl = []

plt.figure(figsize=(16,9))
for st in range(0, 19):
    xbl = []
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( xdatesStation[st][evt] )
        ybl.append( histos.apHbase.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UB Area/Peak HBase "+pmtId, fontsize=24)
plt.xlabel("Since August 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / 25 ns", fontsize=22)
#plt.ylim(2.5,4.5)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllapHbaseNoChcorr"+pmtId+".png", dpi=100)

plt.clf()
print("ubAllapHbase OK")
'''

'''
# *** HBase ***
xbl = []
ybl = []

plt.figure(figsize=(16,9))
for st in range(0, 19):
    xbl = []
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( xdatesStation[st][evt] )
        ybl.append( histos.apHbase.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UB Area/Peak HBase "+pmtId, fontsize=24)
plt.xlabel("Since August 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / 25 ns", fontsize=22)
plt.ylim(2.5,4.5)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllapHbase"+pmtId+".png", dpi=100)

plt.clf()
print("ubAllapHbase OK")
'''

'''
# *** Calib ***
xbl = []
ybl = []

plt.figure(figsize=(16,9))
for st in range(0, 19):
    xbl = []
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( xdatesStation[st][evt] )
        ybl.append( histos.apCalib.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UB Area/Peak Calib "+pmtId, fontsize=24)
plt.xlabel("Since August 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / 25 ns", fontsize=22)
plt.ylim(2.5,6)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllapCalib"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllapCalib OK")
'''
