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

basenameUub = "../pkfit/uubPeakPMT"
basenameUb = "../../nouub/fitHist/pkfit/ubPeakPMT"
outnameHbase = "../../plots/uubPeaktimeHbSt"+str(stList[stchoo])+"PMTs.png"
outnameCalib = "../../plots/uububPeaktimeCaPMT.png"


def getRate(pmtid):
  pmtid = str(pmtid)
  xtime = [[], []]
  vempkhb = [[], []]
  vempkrmshb = [[], []]
  cntvempk = [[], []]
  chihb = [[], []]

  for elec in range(2):
    for i in range( len(stList) ):
      xtime[elec].append( [] )
      vempkhb[elec].append( [] )
      vempkrmshb[elec].append( [] )
      cntvempk[elec].append( [] )
      chihb[elec].append( [] )

  for month in monthUb:
    tmpfile = open(basenameUb+pmtid+month+".json", 'r')
    tmpdata = json.load(tmpfile)
    for st in range(len(stList)):
      for info in tmpdata[str(stList[st])]:
        xtime[0][st].append( datetime.datetime.fromtimestamp( info['Date'] ) )
        vempkhb[0][st].append( info['VemPkHb'] )
        vempkrmshb[0][st].append( info['VemPkrmsHb'] )
        cntvempk[0][st].append( info['CntVemPkHb'] )
        chihb[0][st].append( info['ChiHb'] )

  for month in monthUub:
    tmpfile = open(basenameUub+pmtid+month+".json", 'r')
    tmpdata = json.load(tmpfile)
    for st in range(len(stList)):
      for info in tmpdata[str(stList[st])]:
        xtime[1][st].append( datetime.datetime.fromtimestamp( info['Date'] ) )
        vempkhb[1][st].append( info['VemPkHb'] )
        vempkrmshb[1][st].append( info['VemPkrmsHb'] )
        cntvempk[1][st].append( info['CntVemPkHb'] )
        chihb[1][st].append( info['ChiHb'] )

  return xtime, vempkhb, vempkrmshb


# ===============================
# *** For Peak HBase Plotting ***

pkpmt1 = getRate(1)
pkpmt2 = getRate(2)
pkpmt3 = getRate(3)

markers = ['o', 's', '^']
xfmt = md.DateFormatter('%m/%d')
mkrSize = 6
yminuub = 1e2
ymaxuub = 2.2e2
yminub = 0
ymaxub = 1e2


fig = plt.figure(figsize=(16,9))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

ax11 = ax1.twinx()

ax1.errorbar(pkpmt1[0][0][stchoo], pkpmt1[1][0][stchoo], yerr=pkpmt1[2][0][stchoo], fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax11.errorbar(pkpmt1[0][1][stchoo], pkpmt1[1][1][stchoo], yerr=pkpmt1[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax1.legend(loc=2, title='PMT1')
ax11.legend(loc=0, title='PMT1')
ax1.xaxis.set_major_formatter(xfmt)
ax1.set_ylim(yminub, ymaxub)
ax1.set_ylabel("Vem-Peak [FADC]", fontsize=12)
ax1.tick_params(axis='y', labelsize=12, labelcolor='tab:blue')
ax1.tick_params(axis='x', labelsize=12)
ax11.set_ylabel("Vem-Peak [FADC]", fontsize=12)
ax11.tick_params(axis='both', labelsize=12, labelcolor='tab:orange')
ax11.set_ylim(yminuub, ymaxuub)


ax21 = ax2.twinx()

ax2.errorbar(pkpmt2[0][0][stchoo], pkpmt2[1][0][stchoo], yerr=pkpmt2[2][0][stchoo], fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax21.errorbar(pkpmt2[0][1][stchoo], pkpmt2[1][1][stchoo], yerr=pkpmt2[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax2.legend(loc=2, title='PMT2')
ax21.legend(loc=1, title='PMT2')
ax2.xaxis.set_major_formatter(xfmt)
ax2.set_ylim(yminub, ymaxub)
ax2.set_ylabel("Vem-Peak [FADC]", fontsize=12)
ax2.tick_params(axis='both', labelsize=12)
ax2.tick_params(axis='y', labelsize=12, labelcolor='tab:blue')
ax21.set_ylabel("Vem-Peak [FADC]", fontsize=12)
ax21.tick_params(axis='both', labelsize=12)
ax21.tick_params(axis='y', labelsize=12, labelcolor='tab:orange')
ax21.set_ylim(yminuub, ymaxuub)


ax31 = ax3.twinx()

ax3.errorbar(pkpmt3[0][0][stchoo], pkpmt3[1][0][stchoo], yerr=pkpmt3[2][0][stchoo], fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax31.errorbar(pkpmt3[0][1][stchoo], pkpmt3[1][1][stchoo], yerr=pkpmt3[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax3.legend(loc=2, title='PMT3')
ax31.legend(loc=0, title='PMT3')
ax3.xaxis.set_major_formatter(xfmt)
ax3.set_xlabel("Days since August 1st, 2020 (month/day)", fontsize=12)
ax3.set_ylim(yminub, ymaxub)
ax3.set_ylabel("Vem-Peak [FADC]", fontsize=12)
ax3.tick_params(axis='both', labelsize=12)
ax3.tick_params(axis='y', labelsize=12, labelcolor='tab:blue')
ax31.set_ylabel("Vem-Peak [FADC]", fontsize=12)
ax31.tick_params(axis='both', labelsize=12)
ax31.tick_params(axis='y', labelsize=12, labelcolor='tab:orange')
ax31.set_ylim(yminuub, ymaxuub)


fig.suptitle("Station "+str(stList[stchoo]), fontsize=22)
plt.tight_layout()
#plt.show()
plt.savefig(outnameHbase, dpi=100)
print("Done for Peak over time, ST. "+str(stList[stchoo]))
