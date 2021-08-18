import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# *** Reading JSON files ***

elecver = sys.argv[1]

basename = elecver+"FitaopPMT"
outnameHbase = elecver+"AoPHbasePMTs.png"
outnameCalib = elecver+"AoPCalibPMTs.png"
outnameDiffCaHb = elecver+"AoPDiffCaHbPMTs.png"

npmts = 3

stLabel = []
xst = []
yaopHbpmt = []
yaopHberrpmt = []
yaopCapmt = []
yaopCaerrpmt = []

ydiffHbCa = []
# ========================
# *** Reading for PMTs ***

for pmt in range(0, npmts):
    tmpfile = open(basename+str(pmt+1)+".json", 'r')
    tmpdata = json.load(tmpfile)
    tmpaveAoPHb = []
    tmperrAoPHb = []
    tmpaveAoPCa = []
    tmperrAoPCa = []
    for info in tmpdata[elecver]:
        tmpaveAoPHb.append( info['extAveAoPHb'] )
        tmperrAoPHb.append( info['extRmsAoPHb'] )
        tmpaveAoPCa.append( info['extAveAoPCa'] )
        tmperrAoPCa.append( info['extRmsAoPCa'] )

        if ( pmt==npmts-1 ):
            xst.append( info['stationId'] )
    yaopHbpmt.append( tmpaveAoPHb )
    yaopHberrpmt.append( tmperrAoPHb )
    yaopCapmt.append( tmpaveAoPCa )
    yaopCaerrpmt.append( tmperrAoPCa )

# ==============================
# *** For AoP HBase Plotting ***
markers = ['o', 's', '^']

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.errorbar(xst, yaopHbpmt[i], yerr=yaopHberrpmt[i], fmt=markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
plt.legend(fontsize=22)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Average AoP Peak-HBase [8.33ns]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(-1,12)
plt.tight_layout()
plt.savefig(outnameHbase, dpi=100)
print("Average AoP Peak-HBase OK")


# ==============================
# *** For AoP Calib Plotting ***

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.errorbar(xst, yaopCapmt[i], yerr=yaopCaerrpmt[i], fmt=markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
plt.legend(fontsize=22)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Average AoP Peak-Calib [8.33ns]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(-1,12)
plt.tight_layout()
plt.savefig(outnameCalib, dpi=100)
print("Average AoP Peak-Calib OK")


# =========================================
# *** For AoP Diff Calib/HBase Plotting ***
tmpdiff = []
tmpCa = 0.
tmpHb = 0.

for ipmt in range(0, npmts):
    tmpdiff = []
    for st in range(0, len(xst)):
        tmpCa = yaopCapmt[ipmt][st]
        tmpHb = yaopHbpmt[ipmt][st]
        if tmpHb > 2 and tmpCa > 2:
            tmpdiff.append( 1 - tmpCa/tmpHb )
        else:
             tmpdiff.append( -1 )
    ydiffHbCa.append( tmpdiff )

yerr = []
for pmt in range(0, npmts):
    tmp = []
    for st in range( len(yaopHbpmt[0]) ):
        if yaopHbpmt[pmt][st] > 0:
            tmp.append( ( yaopCaerrpmt[pmt][st]/yaopHbpmt[pmt][st] )**2 + ( yaopHberrpmt[pmt][st]*yaopCapmt[pmt][st]/yaopHbpmt[pmt][st]**2 )**2 )
        else:
            tmp.append( 0 )
    yerr.append( np.sqrt( tmp ) )

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.errorbar(xst, ydiffHbCa[i], yerr=yerr[i], fmt=markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
    #plt.plot(xst, ydiffHbCa[i], markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
plt.legend(fontsize=22, loc=4)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Relative Diff. ( 1 - AoPCalib / AoPHBase) [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(-0.4,0.4)
plt.grid(axis='y', color='0.7')
plt.tight_layout()
plt.savefig(outnameDiffCaHb, dpi=100)
print("Relative Diff. ( 1 - AoPCalib / AoPHBase) OK")

