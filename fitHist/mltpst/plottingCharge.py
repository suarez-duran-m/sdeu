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

inFile = ROOT.TFile.Open("uubCalibHist"+pmtId+".root", "READ")
print("Reading file:", "uubCalibHist"+pmtId+".root")

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

histos = inFile.Get("Histograms")

'''
xdatesStation = []
xtmp = []
for st in range(0, 19):
    xtmp = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xtmp.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
    xdatesStation.append( xtmp )
print("Dates builded")
'''

# =========================
# *** For Charge Offset ***

# *** Average ***
lstentry = histos.GetEntries()
totetrySt = 0
averageStat = []
tmpave = 0
tmpentry = 0

for st in range(0, 19):
    tmpave = 0
    totetrySt = 0
    for evt in range(0, lstentry):
        histos.GetEntry(evt)
        if tmpentry != histos.eventStat[st+1]:
            tmpentry = histos.eventStat[st+1]
            totetrySt = tmpentry
            tmpave += histos.offSetCh.GetBinContent(st)
    averageStat.append( tmpave / totetrySt )

plt.figure(figsize=(16,9))
plt.plot(averageStat, 'o')
plt.title("UUB Average of Offset for Charge Histograms "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Offset (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(50,)
plt.savefig("../../plots/uubAllAveOffsetCh"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllAveOffsetCh OK")


'''
xbl = []
ybl = []

plt.figure(figsize=(16,9))
for st in range(0, 19):
    xbl = []
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.offSetCh.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Offset for Charge histograms "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(2.2e2,-1)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllOffsetCh"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllOffsetCh OK")
'''

# ========================================
# *** For GetBinCenter Charge histograms ***
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
        xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.bincentHbaseCh.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Charge histogram First Bin Center HBase "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllBinCenterHbaseCh"+pmtId+".png", dpi=100)
#plt.ylim(-7.5,-5.5)
#plt.savefig("../../plots/uubAllBinCenterHbaseZoom"+pmtId+".png", dpi=100)

plt.clf()
print("uubAllBinCenterHbase OK")

# *** Calib ***
xbl = []
ybl = []

plt.figure(figsize=(16,9))
for st in range(0, 19):
    xbl = []
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.bincentCalibCh.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Charge histogram First Bin Center Calib "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllBinCenterCalibCh"+pmtId+".png", dpi=100)
#plt.ylim(-15,55)
#plt.savefig("../../plots/uubAllBinCenterCalibZoom"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllBinCenterCalib OK")
'''

# =====================
# *** For Area/Peak ***

# =============================
# *** For Chis distribution ***
'''
# *** HBase ***
xbl = []
ybl = []
distChis = []
for st in range(0, 19):
    tmp = 0
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        tmp = histos.chisHbasePk.GetBinContent(st)
        if tmp > 0:
            ybl.append( tmp )
    if len(ybl) > 0:
        distChis.append( ybl )

#plt.boxplot(distChis, stLabel, patch_artist=True)
plt.violinplot(distChis, showmeans=True)
plt.show()
'''

'''
xbl = []
ybl = []

plt.figure(figsize=(16,9))
for st in range(0, 19):
    xbl = []
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( xdatesStation[st][evt] )
        ybl.append( histos.chisHbasePk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Area/Peak HBase "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / 8.33 ns", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllChisChHbase"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllChisChHbase OK")

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
        ybl.append( histos.chisCalibPk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Area/Peak Calib "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / 8.33 ns", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllChisChCalib"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllChisChCalib OK")
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
plt.title("UUB Area/Peak HBase "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / 8.33 ns", fontsize=22)
plt.ylim(5,11)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllapHbase"+pmtId+".png", dpi=100)

plt.clf()
print("uubAllapHbase OK")

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
        #xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.apCalib.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Area/Peak Calib "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / 8.33 ns", fontsize=22)
plt.ylim(4,11)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllapCalib"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllapCalib OK")
'''


# =============================
# *** For Chis distribution ***
'''
# *** HBase ***
xbl = []
ybl = []
distChis = []
for st in range(0, 19):
    tmp = 0
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        tmp = histos.chisHbaseCh.GetBinContent(st)
        if tmp > 0:
            ybl.append( tmp )
    if len(ybl) > 0:
        distChis.append( ybl )

#plt.boxplot(distChis, stLabel, patch_artist=True)
plt.violinplot(distChis, showmeans=True)
plt.savefig("../../plots/uubAllChisChHbaseVio"+pmtId+".png", dpi=100)

# *** Calib ***
ybl = []
distChis = []
for st in range(0, 19):
    tmp = 0
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        tmp = histos.chisCalibCh.GetBinContent(st)
        if tmp > 0:
            ybl.append( tmp )
    if len(ybl) > 0:
        distChis.append( ybl )

#plt.boxplot(distChis, stLabel, patch_artist=True)
plt.violinplot(distChis, showmeans=True)
plt.savefig("../../plots/uubAllChisChHbaseVio"+pmtId+".png", dpi=100)
'''
