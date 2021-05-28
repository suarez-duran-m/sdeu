import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

stid = sys.argv[1]

monthUb = ['aug', 'sep', 'oct', 'nov']
monthUub = ['dec', 'jan', 'feb', 'mar', 'abr']

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

basenameUub = "../aoptimeTraceSelec/uubAoPtimePMT"
basenameUb = "../../nouub/fitHist/aoptimeTraceSelec/ubAoPtimePMT"
outnameHbase = "../../plots/uubChargetimeHbSt"+stid+"PMTs.png"

fadc2pcub = (25./50.)*(2000./1024) # From FADC to pC
fadc2pcuub = (8.33/50.)*(2000./4096) # From FADC to pC

def getRate(pmtid):
  pmtid = str(pmtid)
  xtime = [[], []]
  vemchhb = [[], []]
  vemchrmshb = [[], []]

  for month in monthUb:
    tmpfile = open(basenameUb+pmtid+"St"+stid+"Mth"+month+".json", 'r')
    tmpdata = json.load(tmpfile)
    for info in tmpdata[month]:
      xtime[0].append( datetime.datetime.fromtimestamp( info['Date'] ) )
      vemchhb[0].append( info['Ch']*fadc2pcub )
      vemchrmshb[0].append( info['Chrms'] )

  for month in monthUub:
    tmpfile = open(basenameUub+pmtid+"St"+stid+"Mth"+month+"chpk.json", 'r')
    tmpdata = json.load(tmpfile)
    for info in tmpdata[month]:
      if info['Date'] > 0:
        xtime[1].append( datetime.datetime.fromtimestamp( info['Date'] ) )
        vemchhb[1].append( info['Ch']*fadc2pcuub )
        vemchrmshb[1].append( info['Chrms'] )      

  return xtime, vemchhb, vemchrmshb


# =================================
# *** For Charge HBase Plotting ***

chpmt1 = getRate(1)
chpmt2 = getRate(2)
chpmt3 = getRate(3)

markers = ['o', 's', '^']
xfmt = md.DateFormatter('%m/%d')
mkrSize = 6
ymin = 40
ymax = 200
fs = 20
ls = 18

fig = plt.figure(figsize=(16,9))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

index = [i for i in range(87,106)]
a = np.array(chpmt1[0][0])
b = np.delete(a, index)
a = np.array(chpmt1[1][0])
c = np.delete(a, index)
a = np.array(chpmt1[2][0])
d = np.delete(a, index)

ax1.errorbar(b, c, yerr=d, fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax1.errorbar(chpmt1[0][1], chpmt1[1][1], yerr=chpmt1[2][1], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax1.legend(loc=2, title='PMT1', title_fontsize=16, fontsize=16)
ax1.xaxis.set_major_formatter(xfmt)
ax1.set_ylim(ymin, ymax)
ax1.set_ylabel("Charge [pC]", fontsize=fs)
ax1.tick_params(axis='y', labelsize=ls, labelright=True)
ax1.tick_params(axis='x', labelsize=ls)

print("\n======================\n")
print("=== PMT1 ===\n")
print("UB:", c[c>70].mean())
a = np.array(chpmt1[1][1])
print("UUB:", a[a>10].mean())

a = np.array(chpmt2[0][0])
b = np.delete(a, index)
a = np.array(chpmt2[1][0])
c = np.delete(a, index)
a = np.array(chpmt2[2][0])
d = np.delete(a, index)

ax2.errorbar(b, c, yerr=d, fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax2.errorbar(chpmt2[0][1], chpmt2[1][1], yerr=chpmt2[2][1], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax2.legend(loc=2, title='PMT2', title_fontsize=16, fontsize=16)
ax2.xaxis.set_major_formatter(xfmt)
ax2.set_ylim(ymin, ymax)
ax2.set_ylabel("Charge [pC]", fontsize=fs)
ax2.tick_params(axis='both', labelsize=ls)
ax2.tick_params(axis='y', labelsize=ls, labelright=True)

print("=== PMT2 ===\n")
print("UB:", c[c>60].mean())
a = np.array(chpmt2[1][1])
print("UUB:", a[a>100].mean())

a = np.array(chpmt3[0][0])
b = np.delete(a, index)
a = np.array(chpmt3[1][0])
c = np.delete(a, index)
a = np.array(chpmt3[2][0])
d = np.delete(a, index)

ax3.errorbar(b, c, yerr=d, fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax3.errorbar(chpmt3[0][1], chpmt3[1][1], yerr=chpmt3[2][1], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax3.legend(loc=2, title='PMT3', title_fontsize=16, fontsize=16)
ax3.xaxis.set_major_formatter(xfmt)
ax3.set_xlabel("Days since August 1st, 2020 (month/day)", fontsize=fs)
ax3.set_ylim(ymin, ymax)
ax3.set_ylabel("Charge [pC]", fontsize=fs)
ax3.tick_params(axis='both', labelsize=ls)
ax3.tick_params(axis='y', labelsize=ls, labelright=True)

print("=== PMT3 ===\n")
print("UB:", c[c>70].mean())
a = np.array(chpmt3[1][1])
print("UUB:", a[a>70].mean())

fig.suptitle("Station "+stid, fontsize=22)
plt.tight_layout()
#plt.show()
plt.savefig(outnameHbase, dpi=100)
print("Done for Charge over time, ST. "+stid)
