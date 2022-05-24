import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# *** Reading JSON files ***

elecver = sys.argv[1]

basename = elecver+"FitXiPMT"
outnameXiPkHb = "../../../plots/"+elecver+"XiPkHbPMTs.png"
outnameXiPkCa = "../../../plots/"+elecver+"XiPkCaPMTs.png"
outnameRateXiPkHb = "../../../plots/"+elecver+"XiRateXiPkHbPMTs.png"
outnameRateXiPkCa = "../../../plots/"+elecver+"XiRateXiPkCaPMTs.png"
outnameRateXi = "../../../plots/"+elecver+"XiRateXiPkHbPMTs.png"

npmts = 3

stLabel = []
xst = []
yXiPkHbpmt = []
yXiPkHberrpmt = []
yXiPkCapmt = []
yXiPkCaerrpmt = []

yXiChHbpmt = []
yXiChHberrpmt = []
yXiChCapmt = []
yXiChCaerrpmt = []

ytotXiPkHb = []
ytotXiPkCa = []
ytotXiChHb = []
ytotXiChCa = []
totEvt = []

# ========================
# *** Reading for PMTs ***

for pmt in range(0, npmts):
    tmpfile = open(basename+str(pmt+1)+".json", 'r')
    tmpdata = json.load(tmpfile)
    tmpaveXiPkHb = []
    tmperrXiPkHb = []
    tmpaveXiPkCa = []
    tmperrXiPkCa = []

    tmpaveXiChHb = []
    tmperrXiChHb = []
    tmpaveXiChCa = []
    tmperrXiChCa = []

    tmptotXiPkHb = []
    tmptotXiPkCa = []
    tmptotXiChHb = []
    tmptotXiChCa = []
    tmptotEvt = []

    for info in tmpdata[elecver]:
        tmpaveXiPkHb.append( info['extAveXiPkHb'] )
        tmperrXiPkHb.append( info['extRmsXiPkHb'] )
        tmpaveXiPkCa.append( info['extAveXiPkCa'] )
        tmperrXiPkCa.append( info['extRmsXiPkCa'] )

        tmpaveXiChHb.append( info['extAveXiChHb'] )
        tmperrXiChHb.append( info['extRmsXiChHb'] )
        tmpaveXiChCa.append( info['extAveXiChCa'] )
        tmperrXiChCa.append( info['extRmsXiChCa'] )

        tmptotXiPkHb.append( info['extTotXiPkHb'] )
        tmptotXiPkCa.append( info['extTotXiPkCa'] )
        tmptotXiChHb.append( info['extTotXiChHb'] )
        tmptotXiChCa.append( info['extTotXiChCa'] )
        tmptotEvt.append( info['extTotEvet'] )

        if ( pmt==npmts-1 ):
            xst.append( info['stationId'] )
    yXiPkHbpmt.append( tmpaveXiPkHb )
    yXiPkHberrpmt.append( tmperrXiPkHb )
    yXiPkCapmt.append( tmpaveXiPkCa )
    yXiPkCaerrpmt.append( tmperrXiPkCa )
    yXiChHbpmt.append( tmpaveXiChHb )
    yXiChHberrpmt.append( tmperrXiChHb )
    yXiChCapmt.append( tmpaveXiChCa )
    yXiChCaerrpmt.append( tmperrXiChCa )

    ytotXiPkHb.append( tmptotXiPkHb )
    ytotXiPkCa.append( tmptotXiPkCa )
    ytotXiChHb.append( tmptotXiChHb )
    ytotXiChCa.append( tmptotXiChCa )
    totEvt.append( tmptotEvt )

# ==================================
# *** For Xi-Peak HBase Plotting ***

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.errorbar(xst, yXiPkHbpmt[i], yerr=yXiPkHberrpmt[i], fmt='o', markersize=12, label=elecver.upper()+" LPMT"+str(i+1))
plt.legend(fontsize=18)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(rotation=60)
plt.ylabel("Average Xi Peak-HBase (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(1,5)
plt.tight_layout()
plt.savefig(outnameXiPkHb, dpi=100)
print("Average Xi Peak-HBase OK")


# ==================================
# *** For Xi-Peak Calib Plotting ***

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.errorbar(xst, yXiPkCapmt[i], yerr=yXiPkCaerrpmt[i], fmt='o', markersize=12, label=elecver.upper()+" LPMT"+str(i+1))
plt.legend(fontsize=18)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(rotation=60)
plt.ylabel("Average Xi Peak-Calib (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(1,5)
plt.tight_layout()
plt.savefig(outnameXiPkCa, dpi=100)
print("Average Xi Peak-Calib OK")


# =========================================
# *** For Xi Diff Calib/HBase Plotting ***
yrateXi = []
tmpdiff = []
tmpCa = 0.
tmpHb = 0.

for ipmt in range(0, npmts):
    tmpdiff = []
    for st in range(0, len(xst)):
        tmpHb = ytotXiPkHb[ipmt][st]
        tmpCa = ytotXiPkCa[ipmt][st]
        tmpdiff.append( tmpHb/totEvt[ipmt][st] )
    yrateXi.append( tmpdiff )

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.plot(xst, yrateXi[i], 'o', markersize=12, label=elecver.upper()+" LPMT"+str(i+1))
plt.legend(fontsize=18)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(rotation=60)
plt.ylabel("Rate of Succesfull Fit For Peak-HBase (au)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(-0.5,0.5)
#plt.grid()
plt.tight_layout()
plt.savefig(outnameRateXiPkHb, dpi=100)
print("Relative Diff. ( 1 - XiCalib / XiHBase) OK")
