import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

stchoo = int(sys.argv[1])
listst = sys.argv[2]
stList = np.loadtxt(listst, int)
stList = np.sort( stList )

monthUb = ['aug', 'sep', 'oct', 'nov']
monthUub = ['dec', 'jan', 'feb', 'mar', 'abr']

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

basenameUub = "../chfit/uubChargePMT"
basenameUb = "../../nouub/fitHist/chfit/ubChargePMT"
outnameHbase = "../../plots/uububChargetimeHbSt"+str(stList[stchoo])+"PMTs.png"
outnameCalib = "../../plots/uububChargetimeCaPMT.png"


def getRate(pmtid):
  pmtid = str(pmtid)
  xtime = [[], []]
  vemchhb = [[], []]
  vemchrmshb = [[], []]
  cntvemch = [[], []]
  chihb = [[], []]

  for elec in range(2):
    for i in range( len(stList) ):
      xtime[elec].append( [] )
      vemchhb[elec].append( [] )
      vemchrmshb[elec].append( [] )
      cntvemch[elec].append( [] )
      chihb[elec].append( [] )

  for month in monthUb:
    tmpfile = open(basenameUb+pmtid+month+".json", 'r')
    tmpdata = json.load(tmpfile)
    for st in range(len(stList)):
      for info in tmpdata[str(stList[st])]:
        xtime[0][st].append( datetime.datetime.fromtimestamp( info['Date'] ) )
        vemchhb[0][st].append( info['VemChHb'] )
        vemchrmshb[0][st].append( info['VemChrmsHb'] )
        cntvemch[0][st].append( info['CntVemChHb'] )
        chihb[0][st].append( info['ChiHb'] )

  for month in monthUub:
    tmpfile = open(basenameUub+pmtid+month+".json", 'r')
    tmpdata = json.load(tmpfile)
    for st in range(len(stList)):
      for info in tmpdata[str(stList[st])]:
        xtime[1][st].append( datetime.datetime.fromtimestamp( info['Date'] ) )
        vemchhb[1][st].append( info['VemChHb'] )
        vemchrmshb[1][st].append( info['VemChrmsHb'] )
        cntvemch[1][st].append( info['CntVemChHb'] )
        chihb[1][st].append( info['ChiHb'] )

  return xtime, vemchhb, vemchrmshb


# =================================
# *** For Charge HBase Plotting ***

chpmt1 = getRate(1)
chpmt2 = getRate(2)
chpmt3 = getRate(3)

markers = ['o', 's', '^']
xfmt = md.DateFormatter('%m/%d')
mkrSize = 6
ymin = 0
ymax = 2e3

fig = plt.figure(figsize=(16,9))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)


ax1.errorbar(chpmt1[0][0][stchoo], chpmt1[1][0][stchoo], yerr=chpmt1[2][0][stchoo], fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax1.errorbar(chpmt1[0][1][stchoo], chpmt1[1][1][stchoo], yerr=chpmt1[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax1.legend(loc=2, title='PMT1')
ax1.xaxis.set_major_formatter(xfmt)
ax1.set_ylim(ymin, ymax)
ax1.set_ylabel("Vem-Charge [FADC]", fontsize=12)
ax1.tick_params(axis='y', labelsize=12, labelright=True)
ax1.tick_params(axis='x', labelsize=12)


ax2.errorbar(chpmt2[0][0][stchoo], chpmt2[1][0][stchoo], yerr=chpmt2[2][0][stchoo], fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax2.errorbar(chpmt2[0][1][stchoo], chpmt2[1][1][stchoo], yerr=chpmt2[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax2.legend(loc=2, title='PMT2')
ax2.xaxis.set_major_formatter(xfmt)
ax2.set_ylim(ymin, ymax)
ax2.set_ylabel("Vem-Charge [FADC]", fontsize=12)
ax2.tick_params(axis='both', labelsize=12)
ax2.tick_params(axis='y', labelsize=12, labelright=True)


ax3.errorbar(chpmt3[0][0][stchoo], chpmt3[1][0][stchoo], yerr=chpmt3[2][0][stchoo], fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax3.errorbar(chpmt3[0][1][stchoo], chpmt3[1][1][stchoo], yerr=chpmt3[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax3.legend(loc=2, title='PMT3')
ax3.xaxis.set_major_formatter(xfmt)
ax3.set_xlabel("Days since August 1st, 2020 (month/day)", fontsize=12)
ax3.set_ylim(ymin, ymax)
ax3.set_ylabel("Vem-Charge [FADC]", fontsize=12)
ax3.tick_params(axis='both', labelsize=12)
ax3.tick_params(axis='y', labelsize=12, labelright=True)


fig.suptitle("Station "+str(stList[stchoo]), fontsize=22)
plt.tight_layout()
#plt.show()
plt.savefig(outnameHbase, dpi=100)
print("Done for Charge over time, ST. "+str(stList[stchoo]))
