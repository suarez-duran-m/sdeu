import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# *** Reading JSON files ***

elecver = sys.argv[1]

basename = elecver+"AllAveHbaseCalibPMT"
outnameHbase = "../../plots/"+elecver+"AllAveHbasePMTs.png"
outnameCalib = "../../plots/"+elecver+"AllAveCalibPMTs.png"
outnameBnCntrHbase = "../../plots/"+elecver+"AllAveBinCntrHbasePMTs.png"
outnameBnCntrCalib = "../../plots/"+elecver+"AllAveBinCntrCalibPMTs.png"

stLabel = []
xst = []
yHbasepmt = []
yHbaseerrpmt = []
yCalibpmt = []
yCaliberrpmt = []

yBnCenterHbasepmt = []
yBnCenterHbaseerrpmt = []
yBnCenterCalibpmt = []
yBnCenterCaliberrpmt = []

# ========================
# *** Reading for PMTs ***

for pmt in range(0, 3):
    tmpfile = open(basename+str(pmt+1)+".json", 'r')
    tmpdata = json.load(tmpfile)
    tmpaveHbase = []
    tmperrHbase = []
    tmpaveCalib = []
    tmperrCalib = []
    tmpaveBnCntrHbase = []
    tmperrBnCntrHbase = []
    tmpaveBnCntrCalib = []
    tmperrBnCntrCalib = []
    for info in tmpdata[elecver]:
        tmpaveHbase.append( info['extAveHbase'] )
        tmperrHbase.append( info['extRmsHbase'] )
        tmpaveCalib.append( info['extAveCalib'] )
        tmperrCalib.append( info['extRmsCalib'] )
        tmpaveBnCntrHbase.append( info['extAveFrstBinCntrHbase'] )
        tmperrBnCntrHbase.append( info['extRmsFrstBinCntrHbase'] )
        tmpaveBnCntrCalib.append( info['extAveFrstBinCntrCalib'] )
        tmperrBnCntrCalib.append( info['extRmsFrstBinCntrCalib'] )
        if ( pmt==2 ):
            xst.append( info['stationId'] )
    yHbasepmt.append( tmpaveHbase )
    yHbaseerrpmt.append( tmperrHbase )
    yCalibpmt.append( tmpaveCalib )
    yCaliberrpmt.append( tmperrCalib )
    yBnCenterHbasepmt.append( tmpaveBnCntrHbase )
    yBnCenterHbaseerrpmt.append( tmperrBnCntrHbase )
    yBnCenterCalibpmt.append( tmpaveBnCntrCalib )
    yBnCenterCaliberrpmt.append( tmperrBnCntrCalib )

# ==========================
# *** For HBase Plotting ***
markers = ['o', 's', '^']

plt.figure(figsize=(16,9))
for i in range(0, 3):
    plt.errorbar(xst, yHbasepmt[i], yerr=yHbaseerrpmt[i], fmt=markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
plt.legend(fontsize=20)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Average Baseline from Mean-HBase [FADC]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
if elecver=="uub":
    plt.ylim(150,320)
else:
    plt.ylim(30,80)
plt.tight_layout()
plt.savefig(outnameHbase, dpi=100)
print("bothAllAveHbase OK")

plt.figure(figsize=(16,9))
for i in range(0, 3):
    plt.errorbar(xst, yBnCenterHbasepmt[i], yerr=yBnCenterHbaseerrpmt[i], fmt=markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
plt.legend(fontsize=20)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Ave. first bin center Mean-HBase [FADC]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#if elecver=="uub":
#    plt.ylim(150,320)
plt.tight_layout()
plt.savefig(outnameBnCntrHbase, dpi=100)
print("bothAllAveBnCntrHbase OK")


# ==========================
# *** For Calib Plotting ***

plt.figure(figsize=(16,9))
for i in range(0, 3):
    plt.errorbar(xst, yCalibpmt[i], yerr=yCaliberrpmt[i], fmt=markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
plt.legend(fontsize=22)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Average Baseline from Calib [FADC]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
if elecver=="uub":
    plt.ylim(150,320)
else:
    plt.ylim(30,80)
plt.tight_layout()
plt.savefig(outnameCalib, dpi=100)
print("bothAllAveCalib OK")


shifted = []
if elecver == 'ub':
    for pmt in range(0, 3):
        tmp = []
        for fstbin in yBnCenterCalibpmt[pmt]:
            tmp.append( fstbin + pmt/10.)
        shifted.append( tmp )
else:
    shifted = yBnCenterCalibpmt

plt.figure(figsize=(16,9))
for i in range(0, 3):
    if elecver == 'ub':
        labels = elecver.upper()+" PMT"+str(i+1) + " + " + str(i/10.)
    else:
        labels = elecver.upper()+" PMT"+str(i+1)
    plt.errorbar(xst, shifted[i], yerr=yBnCenterCaliberrpmt[i], fmt=markers[i], markersize=12, label=labels)
plt.legend(fontsize=22)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Ave. first bin center Calib [FADC]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.tight_layout()
plt.savefig(outnameBnCntrCalib, dpi=100)
print("bothAllAveBnCntrCalib OK")
