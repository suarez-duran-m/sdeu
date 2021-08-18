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
datesOk = int(sys.argv[2])

inFile = ROOT.TFile.Open("uubCalibHist"+pmtId+".root", "READ")
print("Reading file:", "uubCalibHist"+pmtId+".root")

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

histos = inFile.Get("Histograms")

# =================
# *** For Dates ***
'''
xdatesStation = []
xtmp = []
datetmp = 0

if datesOk==1:
    outF = open("xdatesStation.txt", "w")
    xdatesStation = []
    xtmp = []
    for st in range(0, 19):
        xtmp = []
        for evt in range(0, histos.GetEntries()):
            histos.GetEntry(evt)
            datetmp = dt.datetime.fromtimestamp( histos.evtTime )
            xtmp.append( datetmp )
            outF.write("%d "% datetmp)
        outF.write("\n")
        xdatesStation.append( xtmp )
    outF.close
    print("Dates builded and written")
if datesOk==2:
    indates = np.loadtxt("xdatesStation.txt")
'''

# =======================================
# *** For counts in firt bin for Peak ***
'''
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
            tmpave += histos.firstBinCntPk.GetBinContent(st)
    averageStat.append( tmpave / totetrySt )

plt.figure(figsize=(16,9))
plt.plot(averageStat, 'o')
plt.title("UUB Average of counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Counts (au)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../plots/uubAllAveFrstBinCnt"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllAveFrstBinCnt OK")
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
        xbl.append( indates[st][evt] )
        ybl.append( histos.firstBinCntPk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / au", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(50,)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllFrstBinCnt"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllFrstBinCnt OK")

for st in range(0, 19):
    xbl = []
    ybl = []
    if pmt==1:
        if st==2 or st==3 or st==5 or st==6 or st==16 or st==17:
            continue
    elif pmt==2:
        if st==2 or st==3 or st==5 or st==6 or st==9 or st==16 or st==17:
            continue
    elif pmt==3:
        if st==2 or st==3 or st==5 or st==6 or st==9 or st==12 or st==16 or st==17:
            continue
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.firstBinCntPk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / au", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
if pmt==1:
    plt.ylim(50,)
elif pmt==2:
    plt.ylim(50,1200)
elif pmt==3:
     plt.ylim(50,700)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllFrstBinCntFiltUp"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllFrstBinCntFiltUp OK")

for st in range(0, 19):
    xbl = []
    ybl = []
    if pmt==1:
        if st!=2 and st!=3 and st!=5 and st!=6 and st!=16 and st!=17:
            continue
    elif pmt==2:
        if st!=2 and st!=3 and st!=5 and st!=6 and st!=9 and st!=16 and st!=17:
            continue
    elif pmt==3:
        if st!=2 and st!=3 and st!=5 and st!=6 and st!=9 and st!=12 and st!=16 and st!=17:
            continue
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.firstBinCntPk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / au", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(50,)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllFrstBinCntFiltDown"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllFrstBinCntFiltDown OK")

xbl = []
ybl = []
for evt in range(0, histos.GetEntries()):
    histos.GetEntry(evt)
    xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
    ybl.append( histos.firstBinCntPk.GetBinContent(18) )

if pmt==1 or pmt==3:
    plt.scatter(xbl, ybl, color='tab:green', label="St. "+stLabel[18])
elif pmt==2:
    plt.scatter(xbl, ybl, color='tab:orange', label="St. "+stLabel[18])
ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / au", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(50,)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllFrstBinCnt1851"+pmtId+".png", dpi=100)
print("uubAllFrstBinCnt1851 OK")
'''

# =======================
# *** For Peak Offset ***
'''
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
            tmpave += histos.offSetPk.GetBinContent(st)
            histos.GetEntry(lstentry)
    averageStat.append( tmpave / totetrySt )

plt.figure(figsize=(16,9))
plt.plot(averageStat, 'o')
plt.title("UUB Average of Offset for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Offset (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(50,)
plt.savefig("../../plots/uubAllAveOffset"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllAveOffset OK")
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
        xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.offSetPk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Offset for Peak histograms "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
if pmt==1:
    plt.ylim(2.2e2,3e2)
elif pmt==2:
    plt.ylim(2.1e2,2.8e2)
elif pmt==3:
    plt.ylim(2.2e2,2.7e2)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllOffsetPk"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllOffsetPk OK")
'''

# ====================
# *** For Baseline ***
'''
# *** Average ***
lstentry = histos.GetEntries()
averageStat = []
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
            tmpave += histos.baselineHbase.GetBinContent(st)
    averageStat.append( tmpave / totetrySt )

plt.figure(figsize=(16,9))
plt.plot(averageStat, 'o')
plt.title("UUB Average Baseline from mean of HBase "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Mean HBase (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../plots/uubAllAveHBase"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllAveHBase OK")

averageStat = []
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
            tmpave += histos.baselineCalib.GetBinContent(st)
    averageStat.append( tmpave / totetrySt )

plt.figure(figsize=(16,9))
plt.plot(averageStat, 'o')
plt.title("UUB Average Baseline from Calib.Base "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Calib.Base (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../plots/uubAllAveCalib"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllAveCalib OK")
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
        xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.baselineHbase.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Mean Baseline HBase "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
if pmt==1:
    plt.ylim(2.2e2,2.9e2)
elif pmt==2:
    plt.ylim(2.1e2,2.9e2)
elif pmt==3:
    plt.ylim(2.2e2,2.7e2)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllBaselineHbase"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllBaselineHbase OK")

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
        ybl.append( histos.baselineCalib.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Baseline Calib "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
if pmt==1:
    plt.ylim(2.2e2,2.9e2)
elif pmt==2:
    plt.ylim(2.1e2,2.9e2)
elif pmt==3:
    plt.ylim(2.2e2,2.7e2)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllBaselineCalib"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllBaselineCalib OK")
'''

# ========================================
# *** For GetBinCenter Peak histograms ***

# *** Average ***

lstentry = histos.GetEntries()
averageStat = []
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
            tmpave += histos.bincentHbasePk.GetBinContent(st)
    averageStat.append( tmpave / totetrySt )

plt.figure(figsize=(16,9))
plt.plot(averageStat, 'o')
plt.title("UUB Average Peak histogram first bin center HBase "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Baseline from Mean HBase (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../plots/uubAllAveBinCenterHBase"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllAveBinCenterHBase OK")

averageStat = []
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
            tmpave += histos.bincentCalibPk.GetBinContent(st)
    averageStat.append( tmpave / totetrySt )

plt.figure(figsize=(16,9))
plt.plot(averageStat, 'o')
plt.title("UUB Average Peak histogram first bin center Calib "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Baseline form Calib (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../plots/uubAllAveBinCenterCalib"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllAveBinCenterCalib OK")


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
        ybl.append( histos.bincentHbasePk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Peak histogram First Bin Center HBase "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllBinCenterHbase"+pmtId+".png", dpi=100)

plt.ylim(-7.5,-5.5)
plt.savefig("../../plots/uubAllBinCenterHbaseZoom"+pmtId+".png", dpi=100)

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
        ybl.append( histos.bincentCalibPk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UUB Peak histogram First Bin Center Calib "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../plots/uubAllBinCenterCalib"+pmtId+".png", dpi=100)
plt.ylim(-15,55)
plt.savefig("../../plots/uubAllBinCenterCalibZoom"+pmtId+".png", dpi=100)
plt.clf()
print("uubAllBinCenterCalib OK")
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
        tmp = histos.chisHbasePk.GetBinContent(st)
        if tmp > 0:
            ybl.append( tmp )
    if len(ybl) > 0:
        distChis.append( ybl )

#plt.boxplot(distChis, stLabel, patch_artist=True)
plt.violinplot(distChis, showmeans=True)
plt.savefig("../../plots/uubAllChisPkHbaseVio"+pmtId+".png", dpi=100)

# *** Calib ***
ybl = []
distChis = []
for st in range(0, 19):
    tmp = 0
    ybl = []
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        tmp = histos.chisCalibPk.GetBinContent(st)
        if tmp > 0:
            ybl.append( tmp )
    if len(ybl) > 0:
        distChis.append( ybl )

#plt.boxplot(distChis, stLabel, patch_artist=True)
plt.violinplot(distChis, showmeans=True)
plt.savefig("../../plots/uubAllChisPkHbaseVio"+pmtId+".png", dpi=100)
'''
