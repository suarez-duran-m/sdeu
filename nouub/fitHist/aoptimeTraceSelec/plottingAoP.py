import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

monthUb = sys.argv[1]
monthUub = sys.argv[2]

basenameUub = "uubAoPtimePMT1St863Mth"+monthUub+".json"
basenameUb = "../../nouub/fitHist/aoptime/ubAoPtimePMT1St863Mth"+monthUb+".json"
#outnameHbase = "../../plots/"+elecver+"AoPHbasePMTs.png"
#outnameCalib = "../../plots/"+elecver+"AoPCalibPMTs.png"
#outnameDiffCaHb = "../../plots/"+elecver+"AoPDiffCaHbPMTs.png"

xtime = [[], []]
yaopHb = [[], []]
yaopHberr = [[], []]
yaopCa = [[], []]
yaopCaerr = [[], []]

# ========================
# *** Reading for PMTs ***

factor = 8.33/25.

tmpfile = open(basenameUub, 'r')
tmpdata = json.load(tmpfile)

for info in tmpdata[monthUub]:
    xtime[0].append( datetime.datetime.fromtimestamp( info['Date'] ) )
    yaopHb[0].append( info['AoP']*factor )
    yaopHberr[0].append( info['rms'] )

tmpfile = open(basenameUb, 'r')
tmpdata = json.load(tmpfile)
for info in tmpdata[monthUb]:
    xtime[1].append( datetime.datetime.fromtimestamp( info['Date'] ) )
    yaopHb[1].append( info['AoP'] )
    yaopHberr[1].append( info['rms'] )

# ==============================
# *** For AoP HBase Plotting ***
markers = ['o', 's', '^']

plt.figure(figsize=(16,9))
plt.errorbar(xtime[1], yaopHb[1], yerr=yaopHberr[1], fmt=markers[1], markersize=12, label="UB")
plt.errorbar(xtime[0], yaopHb[0], yerr=yaopHberr[0], fmt=markers[0], markersize=12, label="UUB * 8.33/25.")

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)

plt.legend(fontsize=22)
plt.xlabel("Days since November 1st, 2020 (month/day)", fontsize=24)
#plt.xticks(rotation=60)
plt.ylabel("Average AoP per day Peak-HBase", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.tight_layout()
plt.show()
#plt.savefig(outnameHbase, dpi=100)
print("Average AoP Peak-HBase OK")


'''
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
if elecver=="uub":
    plt.ylim(-1,12)
#else:
#    plt.ylim(30,80)
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
        if tmpHb > 4 and tmpCa > 4:
            tmpdiff.append( 1 - tmpCa/tmpHb )
        else:
             tmpdiff.append( -1 )
    ydiffHbCa.append( tmpdiff )

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.plot(xst, ydiffHbCa[i], markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
plt.legend(fontsize=22, loc=4)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Relative Diff. ( 1 - AoPCalib / AoPHBase) [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(-0.4,0.4)
plt.grid(axis='y', color='0.75')
plt.tight_layout()
plt.savefig(outnameDiffCaHb, dpi=100)
print("Relative Diff. ( 1 - AoPCalib / AoPHBase) OK")
'''
