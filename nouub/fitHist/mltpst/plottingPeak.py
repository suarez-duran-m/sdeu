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

inFile = ROOT.TFile.Open("ubCalibHist"+pmtId+".root", "READ")
print("Reading file:", "ubCalibHist"+pmtId+".root")

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
lstentry = 32000 #histos.GetEntries()
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
plt.title("UB Average of counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Counts (au)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../../plots/ubAllAveFrstBinCnt"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllAveFrstBinCnt OK")
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
        ybl.append( histos.firstBinCntPk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UB Counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Since August 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / au", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(50,)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllFrstBinCnt"+pmtId+".png", dpi=100)
plt.clf()

print("ubAllFrstBinCnt OK")

for st in range(0, 19):
    xbl = []
    ybl = []
    if pmt==1:
        if st==16: 
            continue
    elif pmt==3:
        if st==18:
            continue
    for evt in range(0, histos.GetEntries()):
        histos.GetEntry(evt)
        xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
        ybl.append( histos.firstBinCntPk.GetBinContent(st) )
    plt.scatter(xbl, ybl, label="St. "+stLabel[st])

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UB Counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Since August 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / au", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(-1,10)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllFrstBinCntFiltUp"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllFrstBinCntFiltUp OK")

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
plt.title("UB Counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Since August 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / au", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(50,)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllFrstBinCntFiltDown"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllFrstBinCntFiltDown OK")

xbl = []
ybl = []
for evt in range(0, histos.GetEntries()):
    histos.GetEntry(evt)
    xbl.append( dt.datetime.fromtimestamp( histos.entryEvt ) )
    ybl.append( histos.firstBinCntPk.GetBinContent(0) )

if pmt==1 or pmt==3:
    plt.scatter(xbl, ybl, color='tab:green', label="St. "+stLabel[18])
elif pmt==2:
    plt.scatter(xbl, ybl, color='tab:orange', label="St. "+stLabel[18])
ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title("UB Counts in first bin for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Since August 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / au", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(-1,20)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllFrstBinCnt1851"+pmtId+".png", dpi=100)
print("ubAllFrstBinCnt1851 OK")
'''


# =======================
# *** For Peak Offset ***
'''
# *** Average ***
lstentry = 32000 #histos.GetEntries()
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
    averageStat.append( tmpave / totetrySt )

plt.figure(figsize=(16,9))
plt.plot(averageStat, 'o')
plt.title("UB Average of Offset for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Offset (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(50,)
plt.savefig("../../../plots/ubAllAveOffset"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllAveOffset OK")
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
plt.title("UB Offset for Peak histograms "+pmtId, fontsize=24)
plt.xlabel("Since August 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
if pmt==1:
    plt.ylim(35,75)
elif pmt==2:
    plt.ylim(30,70)
elif pmt==3:
    plt.ylim(28,68)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllOffsetPk"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllFrstBinCnt OK")
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
plt.title("UB Average Baseline from mean of HBase "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Mean HBase (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../../plots/ubAllAveHBase"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllAveHBase OK")

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
plt.title("UB Average Baseline from Calib.Base "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Calib.Base (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../../plots/ubAllAveCalib"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllAveCalib OK")
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
plt.title("UB Mean Baseline HBase "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

if pmt==1:
    plt.ylim(35,75)
elif pmt==2:
    plt.ylim(30,70)
elif pmt==3:
    plt.ylim(25,70)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllBaselineHbase"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllBaselineHbase OK")

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
plt.title("UB Baseline Calib "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

if pmt==1:
    plt.ylim(35,75)
elif pmt==2:
    plt.ylim(30,70)
elif pmt==3:
    plt.ylim(25,70)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllBaselineCalib"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllBaselineCalib OK")
'''

# ========================================
# *** For GetBinCenter Peak histograms ***

# *** Average ***

lstentry = 32000 #histos.GetEntries()
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
plt.title("UB Average Peak histogram first bin center HBase "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Baseline from Mean HBase (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../../plots/ubAllAveBinCenterHBase"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllAveBinCenterHBase OK")

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
plt.title("UB Average Peak histogram first bin center Calib "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(np.arange(19), stLabel, rotation='vertical')
plt.ylabel("Average Baseline from Calib (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../../plots/ubAllAveBinCenterCalib"+pmtId+".png", dpi=100)
plt.clf()
print("ubAllAveBinCenterCalib OK")

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
plt.title("UB Peak histogram First Bin Center HBase "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.ylim(-2,2)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllBinCenterHbase"+pmtId+".png", dpi=100)

plt.clf()
print("ubAllBinCenterHbase OK")

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
plt.title("UB Peak histogram First Bin Center Calib "+pmtId, fontsize=24)
plt.xlabel("Since December 1st, 2020 (month/day)", fontsize=22)
plt.ylabel("1. / FADC", fontsize=22)
plt.ylim(-2,2)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18, loc='upper right')
plt.savefig("../../../plots/ubAllBinCenterCalib"+pmtId+".png", dpi=100)

plt.clf()
print("ubAllBinCenterCalib OK")
'''

'''
# =============================
# *** For Chis distribution ***

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
plt.savefig("../../plots/ubAllChisPkHbaseVio"+pmtId+".png", dpi=100)

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
plt.savefig("../../plots/ubAllChisPkHbaseVio"+pmtId+".png", dpi=100)
'''
