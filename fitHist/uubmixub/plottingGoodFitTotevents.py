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
outnamePkChiUb = "../../plots/ubGoodFitEvtnPkChis.png"
outnamePkChiUub = "../../plots/uubGoodFitEvtnPkChis.png"

basenameChUub = "../chfit/uubChargePMT"
basenameChUb = "../../nouub/fitHist/chfit/ubChargePMT"
outnameChUub = "../../plots/uubGoodFitEvtnChPMT.png"
outnameChUb = "../../plots/ubGoodFitEvtnChPMT.png"
outnameChChiUb = "../../plots/ubGoodFitEvtnChChis.png"
outnameChChiUub = "../../plots/uubGoodFitEvtnChChis.png"

def getRate(pmtid, ifpeak):
  totEvt = [[], []]
  chis = [[], []]
  chisrms = [[], []]
  pmtid = str(pmtid)
  for elec in range(2):
    for i in range( len(stList) ):
      totEvt[elec].append( [] )
      chis[elec].append( [] )
      chisrms[elec].append( [] )

  if ifpeak:
    basename = basenamePkUb
  else:
    basename = basenameChUb   

  for st in range(len(stList)):
    tmp = []
    tmpchis = []
    for month in monthUb:
      tmpfile = open(basename+pmtid+month+".json", 'r')
      tmpdata = json.load(tmpfile)
      for info in tmpdata[str(stList[st])]:
        tmp.append( info['TotEvt'] )
        tmpchis.append( info['ChiHb'] )
    totEvt[0][st].append( np.average(np.array(tmp) ) )
    tmpave = np.average( np.array(tmpchis) )
    tmprms = 0
    for i in tmpchis:
      tmprms += (tmpave - i)*(tmpave - i)
    chis[0][st].append( tmpave )
    chisrms[0][st].append( np.sqrt( tmprms/len(tmpchis) ) )

  if ifpeak:
    basename = basenamePkUub
  else:
    basename = basenameChUub
  for st in range(len(stList)):
    tmp = []
    tmpchis = []
    for month in monthUub:
      tmpfile = open(basename+pmtid+month+".json", 'r')
      tmpdata = json.load(tmpfile)
      for info in tmpdata[str(stList[st])]:
        tmp.append( info['TotEvt'] )
        tmpchis.append( info['ChiHb'] )
    totEvt[1][st].append( np.average(np.array(tmp) ) )
    tmpave = np.average( np.array(tmpchis) )
    tmprms = 0
    for i in tmpchis:
      tmprms += (tmpave - i)*(tmpave - i)
    chis[1][st].append( np.average(np.array(tmpchis) ) )
    chisrms[1][st].append( np.sqrt( tmprms/len(tmpchis) ) )

  return totEvt, chis, chisrms


# ==========================
# *** Reading for Months ***

totEvtpmt1, chispmt1, chisrmspmt1 = getRate(1, True)
totEvtpmt2, chispmt2, chisrmspmt2 = getRate(2, True)
totEvtpmt3, chispmt3, chisrmspmt3 = getRate(3, True)

totEvtpmt1 = np.array(totEvtpmt1)
totEvtpmt2 = np.array(totEvtpmt2)
totEvtpmt3 = np.array(totEvtpmt3)

chispmt1 = np.array(chispmt1)
chispmt2 = np.array(chispmt2)
chispmt3 = np.array(chispmt3)

chisrmspmt1 = np.array(chisrmspmt1)
chisrmspmt2 = np.array(chisrmspmt2)
chisrmspmt3 = np.array(chisrmspmt3)

markers = ['o', 's', '^']

# =======================
# *** For Peak Histos ***

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

# ================
# *** For Chis ***

kk = [[], [], []]
for i in range( len(chisrmspmt1[0]) ):
  kk[0].append( float( chisrmspmt1[0][i] ) )
  kk[1].append( float( chisrmspmt2[0][i] ) )
  kk[2].append( float( chisrmspmt3[0][i] ) )

plt.figure(figsize=(16,9))
plt.errorbar(stLabel, chispmt1[0], yerr=kk[0], fmt=markers[0], markersize=12, label="UB PMT1")
plt.errorbar(stLabel, chispmt2[0], yerr=kk[1], fmt=markers[1], markersize=12, label="UB PMT2")
plt.errorbar(stLabel, chispmt3[0], yerr=kk[2], fmt=markers[2], markersize=12, label="UB PMT3")

plt.legend(title="$\chi^2/$NDF For Peak histograms", fontsize=22, title_fontsize=22)
plt.xlabel("Station Id", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("$\chi^2/$NDF [FADC]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim(0, 1)
plt.tight_layout()
#plt.show()
plt.savefig(outnamePkChiUb, dpi=100)
print("UB Chis Peak OK")

kk = [[], [], []]
for i in range( len(chisrmspmt1[1]) ):
  kk[0].append( float( chisrmspmt1[1][i] ) )
  kk[1].append( float( chisrmspmt2[1][i] ) )
  kk[2].append( float( chisrmspmt3[1][i] ) )

plt.figure(figsize=(16,9))
plt.errorbar(stLabel, chispmt1[0], yerr=kk[0], fmt=markers[0], markersize=12, label="UUB PMT1")
plt.errorbar(stLabel, chispmt2[0], yerr=kk[1], fmt=markers[1], markersize=12, label="UUB PMT2")
plt.errorbar(stLabel, chispmt3[0], yerr=kk[2], fmt=markers[2], markersize=12, label="UUB PMT3")

plt.legend(title="$\chi^2/$NDF For Peak histograms", fontsize=22, title_fontsize=22)
plt.xlabel("Station Id", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("$\chi^2/$NDF [FADC]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim(0, 1)
plt.tight_layout()
#plt.show()
plt.savefig(outnamePkChiUub, dpi=100)
print("UUB Chis Peak OK")



# =========================
# *** For Charge Histos ***

# ==========================
# *** Reading for Months ***

totEvtpmt1, chispmt1, chisrmspmt1 = getRate(1, False)
totEvtpmt2, chispmt2, chisrmspmt2 = getRate(2, False)
totEvtpmt3, chispmt3, chisrmspmt3 = getRate(3, False)

totEvtpmt1 = np.array(totEvtpmt1)
totEvtpmt2 = np.array(totEvtpmt2)
totEvtpmt3 = np.array(totEvtpmt3)

chispmt1 = np.array(chispmt1)
chispmt2 = np.array(chispmt2)
chispmt3 = np.array(chispmt3)

chisrmspmt1 = np.array(chisrmspmt1)
chisrmspmt2 = np.array(chisrmspmt2)
chisrmspmt3 = np.array(chisrmspmt3)

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


# ================
# *** For Chis ***

kk = [[], [], []]
for i in range( len(chisrmspmt1[0]) ):
  kk[0].append( float( chisrmspmt1[0][i] ) )
  kk[1].append( float( chisrmspmt2[0][i] ) )
  kk[2].append( float( chisrmspmt3[0][i] ) )

plt.figure(figsize=(16,9))
plt.errorbar(stLabel, chispmt1[0], yerr=kk[0], fmt=markers[0], markersize=12, label="UB PMT1")
plt.errorbar(stLabel, chispmt2[0], yerr=kk[1], fmt=markers[1], markersize=12, label="UB PMT2")
plt.errorbar(stLabel, chispmt3[0], yerr=kk[2], fmt=markers[2], markersize=12, label="UB PMT3")

plt.legend(title="For Charge histograms", fontsize=22, title_fontsize=22)
plt.xlabel("Station Id", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("$\chi^2/$NDF [FADC]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim(0, 1)
plt.tight_layout()
#plt.show()
plt.savefig(outnamePkChiUb, dpi=100)
print("UB Chis for Charge OK")

kk = [[], [], []]
for i in range( len(chisrmspmt1[1]) ):
  kk[0].append( float( chisrmspmt1[1][i] ) )
  kk[1].append( float( chisrmspmt2[1][i] ) )
  kk[2].append( float( chisrmspmt3[1][i] ) )

plt.figure(figsize=(16,9))
plt.errorbar(stLabel, chispmt1[1], yerr=kk[0], fmt=markers[0], markersize=12, label="UUB PMT1")
plt.errorbar(stLabel, chispmt2[1], yerr=kk[1], fmt=markers[1], markersize=12, label="UUB PMT2")
plt.errorbar(stLabel, chispmt3[1], yerr=kk[2], fmt=markers[2], markersize=12, label="UUB PMT3")

plt.legend(title="$\chi^2/$NDF for Charge histograms", fontsize=22, title_fontsize=22)
plt.xlabel("Station Id", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Events-Fitted / Events-Total [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim(0, 1)
plt.tight_layout()
#plt.show()
plt.savefig(outnamePkChiUub, dpi=100)
print("UUB Chis for Charge OK")
