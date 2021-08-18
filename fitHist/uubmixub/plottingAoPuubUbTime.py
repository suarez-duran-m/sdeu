import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

stid = sys.argv[1]
pmtid = sys.argv[2]

monthUb = ['aug', 'sep', 'oct', 'nov']
monthUub = ['dec', 'jan', 'feb', 'mar', 'abr']

basenameUub = "../aoptimeTraceSelec/uubAoPtimePMT"+pmtid+"St"+stid+"Mth" #+monthUub+".json"
basenameUb = "../../nouub/fitHist/aoptimeTraceSelec/ubAoPtimePMT"+pmtid+"St"+stid+"Mth" #+monthUb+".json"
outnameHbase = "../../plots/uububAoPtimeHbPMT"+pmtid+"St"+stid+".png"

xtime = [[], []]
yaopHb = [[], []]
yaopHberr = [[], []]
yaopCa = [[], []]
yaopCaerr = [[], []]

# ========================
# *** Reading for PMTs ***

chfactorUb = 25./50.
chfactorUub = 8.33/50.

for month in monthUb:
  tmpfile = open(basenameUb+month+".json", 'r')
  tmpdata = json.load(tmpfile)
  for info in tmpdata[month]:
    xtime[0].append( datetime.datetime.fromtimestamp( info['Date'] ) )
    yaopHb[0].append( chfactorUb*info['AoP'] ) # In nF
    yaopHberr[0].append( chfactorUb*info['rms'] ) # In nF

for month in monthUub:
  tmpfile = open(basenameUub+month+"chpk.json", 'r')
  tmpdata = json.load(tmpfile)
  for info in tmpdata[month]:
    if info['Date'] < 1606867200: # December 2nd
      continue
    xtime[1].append( datetime.datetime.fromtimestamp( info['Date'] ) )
    yaopHb[1].append( chfactorUub*info['AoP'] ) # In nF
    yaopHberr[1].append( chfactorUub*info['rms'] ) # In nF

# ==============================
# *** For AoP HBase Plotting ***
markers = ['o', 's', '^']

index = [i for i in range(87,106)]

a = np.array(xtime[0])
b = np.delete(a, index)
a = np.array(yaopHb[0])
c = np.delete(a, index)
a = np.array(yaopHberr[0])
d = np.delete(a, index)

plt.figure(figsize=(16,9))
plt.errorbar(b, c, yerr=d, fmt=markers[0], markersize=12, label="UB")
plt.errorbar(xtime[1], yaopHb[1], yerr=yaopHberr[1], fmt=markers[1], markersize=12, label="UUB")

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)

plt.legend(title="St. "+stid+" PMT"+pmtid, fontsize=22, title_fontsize=22)
plt.xlabel("Days since August 1st, 2020 (month/day)", fontsize=24)
plt.ylabel("Average AoP per day Peak-HBase [nF]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
if stid == '863' and pmtid=="1":
  plt.ylim(1.2, 2.0) 
elif stid == '863' and pmtid=="2":
  plt.ylim(1.2, 2.2)
elif stid == '863' and pmtid=="3":
  plt.ylim(1.2, 1.9)
elif stid == '1740':
  plt.ylim(1.2, 2.0)

plt.tight_layout()
#plt.show()
plt.savefig(outnameHbase, dpi=100)
print("Average AoP Peak-HBase OK")


