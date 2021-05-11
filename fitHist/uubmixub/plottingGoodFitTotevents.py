import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

listst = sys.argv[1]
stList = np.loadtxt(listst, int)

monthUb = ['aug', 'sep', 'oct'] #, 'nov']
monthUub = ['dec', 'jan', 'feb', 'mar', 'abr']

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

basenamePkUub = "../pkfit/uubPeakPMT"
basenamePkUb = "../../nouub/fitHist/pkfit/ubPeakPMT"
outnamePkUub = "../../plots/uubGoodFitEvtnPk.png"
outnamePkUb = "../../plots/ubGoodFitEvtnPk.png"

basenameChUub = "../chfit/uubChargePMT"
basenameChUb = "../../nouub/fitHist/chfit/ubChargePMT"
outnameChUub = "../../plots/uubGoodFitEvtnChPMT.png"
outnameChUb = "../../plots/ubGoodFitEvtnChPMT.png"


def getRate(pmtid, ifpeak):
  totEvt= [[], []]
  pmtid = str(pmtid)
  for elec in range(2):
    for i in range( len(stList) ):
      totEvt[elec].append( [] )

  if ifpeak:
    basename = basenamePkUb
  else:
    basename = basenameChUb   

  for st in range(len(stList)):
    tmp = []
    for month in monthUb:
      tmpfile = open(basename+pmtid+month+".json", 'r')
      tmpdata = json.load(tmpfile)
      for info in tmpdata[str(stList[st])]:
        tmp.append( info['TotEvt'] )
    totEvt[0][st].append( np.average(np.array(tmp) ) )

  if ifpeak:
    basename = basenamePkUub
  else:
    basename = basenameChUub
  for st in range(len(stList)):
    tmp = []
    for month in monthUub:
      tmpfile = open(basename+pmtid+month+".json", 'r')
      tmpdata = json.load(tmpfile)
      for info in tmpdata[str(stList[st])]:
        tmp.append( info['TotEvt'] )
    totEvt[1][st].append( np.average(np.array(tmp) ) )
 
  return totEvt


# ==========================
# *** Reading for Months ***


totEvtpmt1 = getRate(1, True)
totEvtpmt2 = getRate(2, True)
totEvtpmt3 = getRate(3, True)

totEvtpmt1 = np.array(totEvtpmt1)
totEvtpmt2 = np.array(totEvtpmt2)
totEvtpmt3 = np.array(totEvtpmt3)

# =======================
# *** For Peak Histos ***

markers = ['o', 's', '^']
plt.figure(figsize=(16,9))
plt.plot(stLabel, totEvtpmt1[0], markers[0], markersize=12, label="UB PMT1")
plt.plot(stLabel, totEvtpmt2[0], markers[1], markersize=12, label="UB PMT2")
plt.plot(stLabel, totEvtpmt3[0], markers[2], markersize=12, label="UB PMT3")

plt.legend(title="For Peak histograms", fontsize=22, title_fontsize=22)
plt.xlabel("Station Id", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Events-Fitted / Events-Total [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0, 1)
plt.tight_layout()
#plt.show()
plt.savefig(outnamePkUb, dpi=100)
print("UB Peak Events-Fitted / Events-Total OK")



markers = ['o', 's', '^']
plt.figure(figsize=(16,9))
plt.plot(stLabel, totEvtpmt1[1], markers[0], markersize=12, label="UUB PMT1")
plt.plot(stLabel, totEvtpmt2[1]-.02, markers[1], markersize=12, label="UUB PMT2 - 0.02")
plt.plot(stLabel, totEvtpmt3[1]-.04, markers[2], markersize=12, label="UUB PMT3 - 0.04")

plt.legend(title="For Peak histograms", fontsize=22, title_fontsize=22)
plt.xlabel("Station Id", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Events-Fitted / Events-Total [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0, 1)
plt.tight_layout()
#plt.show()
plt.savefig(outnamePkUub, dpi=100)
print("UUB Peak Events-Fitted / Events-Total OK")


# =========================
# *** For Charge Histos ***


# ==========================
# *** Reading for Months ***


totEvtpmt1 = getRate(1, False)
totEvtpmt2 = getRate(2, False)
totEvtpmt3 = getRate(3, False)

totEvtpmt1 = np.array(totEvtpmt1)
totEvtpmt2 = np.array(totEvtpmt2)
totEvtpmt3 = np.array(totEvtpmt3)

plt.figure(figsize=(16,9))
plt.plot(stLabel, totEvtpmt1[0], markers[0], markersize=12, label="UB PMT1")
plt.plot(stLabel, totEvtpmt2[0], markers[1], markersize=12, label="UB PMT2")
plt.plot(stLabel, totEvtpmt3[0], markers[2], markersize=12, label="UB PMT3")

plt.legend(title="For Charge histograms", fontsize=22, title_fontsize=22)
plt.xlabel("Station Id", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Events-Fitted / Events-Total [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0, 1)
plt.tight_layout()
#plt.show()
plt.savefig(outnameChUb, dpi=100)
print("UB Charge Events-Fitted / Events-Total OK")


plt.figure(figsize=(16,9))
plt.plot(stLabel, totEvtpmt1[1], markers[0], markersize=12, label="UUB PMT1")
plt.plot(stLabel, totEvtpmt2[1]-.02, markers[1], markersize=12, label="UUB PMT2 - 0.02")
plt.plot(stLabel, totEvtpmt3[1]-.04, markers[2], markersize=12, label="UUB PMT3 - 0.04")

plt.legend(title="For Charge histograms", fontsize=22, title_fontsize=22)
plt.xlabel("Station Id", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Events-Fitted / Events-Total [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0, 1)
plt.tight_layout()
#plt.show()
plt.savefig(outnameChUub, dpi=100)
print("UUB Charge Events-Fitted / Events-Total OK")

