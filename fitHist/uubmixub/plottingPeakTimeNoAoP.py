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

fadc2mvub = 1.95;
fadc2mvuub = 0.49;

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
        vempkhb[0][st].append( info['VemPkHb']*fadc2mvub )
        vempkrmshb[0][st].append( info['VemPkrmsHb'] )
        cntvempk[0][st].append( info['CntVemPkHb'] )
        chihb[0][st].append( info['ChiHb'] )

  for month in monthUub:
    tmpfile = open(basenameUub+pmtid+month+".json", 'r')
    tmpdata = json.load(tmpfile)
    for st in range(len(stList)):
      for info in tmpdata[str(stList[st])]:
        xtime[1][st].append( datetime.datetime.fromtimestamp( info['Date'] ) )
        vempkhb[1][st].append( info['VemPkHb']*fadc2mvuub )
        vempkrmshb[1][st].append( info['VemPkrmsHb'] )
        cntvempk[1][st].append( info['CntVemPkHb'] )
        chihb[1][st].append( info['ChiHb'] )

  return xtime, vempkhb, vempkrmshb

# ===============================
# *** For Peak HBase Plotting ***

pkpmt1 = getRate(1)
pkpmt2 = getRate(2)
pkpmt3 = getRate(3)

print( np.average( np.array( pkpmt1[1][0][stchoo] ) ) )
print( np.average( np.array( pkpmt2[1][0][stchoo] ) ) )
print( np.average( np.array( pkpmt3[1][0][stchoo] ) ) )
print("")
print( np.average( np.array( pkpmt1[1][1][stchoo] ) ) )
print( np.average( np.array( pkpmt2[1][1][stchoo] ) ) )
print( np.average( np.array( pkpmt3[1][1][stchoo] ) ) )


markers = ['o', 's', '^']
xfmt = md.DateFormatter('%m/%d')
mkrSize = 6
fs = 20
ls = 18

if 
ymin = 40
ymax = 100


index = [i for i in range(87,106)]
a = np.array(pkpmt1[0][0][stchoo])
b = np.delete(a, index)
a = np.array(pkpmt1[1][0][stchoo])
c = np.delete(a, index)
a = np.array(pkpmt1[2][0][stchoo])
d = np.delete(a, index)

fig = plt.figure(figsize=(16,9))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

ax1.errorbar(pkpmt1[0][0][stchoo], pkpmt1[1][0][stchoo], yerr=pkpmt1[2][0][stchoo], fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax1.errorbar(pkpmt1[0][1][stchoo], pkpmt1[1][1][stchoo], yerr=pkpmt1[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax1.legend(loc=2, title='PMT1', title_fontsize=16, fontsize=16)
ax1.xaxis.set_major_formatter(xfmt)
ax1.set_ylim(ymin, ymax)
ax1.set_ylabel("Amplitude [mV]", fontsize=fs)
ax1.tick_params(axis='y', labelsize=ls, labelright=True)
ax1.tick_params(axis='x', labelsize=ls)

print("\n======================\n")
print("=== PMT1 ===\n")
print("UB:", np.mean(c))
a = np.array(pkpmt1[1][1][stchoo])
print("UUB:", a[a>0].mean())


ymin = 40
ymax = 100 #1e2 #2.2e2

a = np.array(pkpmt2[0][0][stchoo])
b = np.delete(a, index)
a = np.array(pkpmt2[1][0][stchoo])
c = np.delete(a, index)
a = np.array(pkpmt2[2][0][stchoo])
d = np.delete(a, index)

ax2.errorbar(b, c, yerr=d, fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax2.errorbar(pkpmt2[0][1][stchoo], pkpmt2[1][1][stchoo], yerr=pkpmt2[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax2.legend(loc=0, title='PMT2', title_fontsize=16, fontsize=16)
ax2.xaxis.set_major_formatter(xfmt)
ax2.set_ylim(ymin, ymax)
ax2.set_ylabel("Amplitude [mV]", fontsize=fs)
ax2.tick_params(axis='both', labelsize=ls)
ax2.tick_params(axis='y', labelsize=ls, labelright=True)

print("\n======================\n")
print("=== PMT2 ===\n")
print("UB:", np.mean(c))
a = np.array(pkpmt2[1][1][stchoo])
print("UUB:", a[a>80].mean())

ymin = 40
ymax = 100

index = [i for i in range(88,105)]
a = np.array(pkpmt2[0][0][stchoo])
b = np.delete(a, index)
a = np.array(pkpmt2[1][0][stchoo])
c = np.delete(a, index)
a = np.array(pkpmt2[2][0][stchoo])
d = np.delete(a, index)

ax3.errorbar(b, c, yerr=d, fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax3.errorbar(pkpmt3[0][1][stchoo], pkpmt3[1][1][stchoo], yerr=pkpmt3[2][1][stchoo], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax3.legend(loc=2, title='PMT3', title_fontsize=16, fontsize=16)
ax3.xaxis.set_major_formatter(xfmt)
ax3.set_xlabel("Days since August 1st, 2020 (month/day)", fontsize=fs)
ax3.set_ylim(ymin, ymax)
ax3.set_ylabel("Amplitude [mV]", fontsize=fs)
ax3.tick_params(axis='both', labelsize=ls)
ax3.tick_params(axis='y', labelsize=ls, labelright=True)

fig.suptitle("Station "+str(stList[stchoo]), fontsize=22)
plt.tight_layout()
#plt.show()
plt.savefig(outnameHbase, dpi=100)

print("\n======================\n")
print("=== PMT3 ===\n")
print("UB:", np.mean(c))
a = np.array(pkpmt3[1][1][stchoo])
print("UUB:", a[a>70].mean(),"\n")

print("Done for Peak over time, ST. "+str(stList[stchoo]))
