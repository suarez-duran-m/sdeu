import ROOT
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

pmtid = sys.argv[1]
ListSt = sys.argv[2]
stlist = np.sort(np.loadtxt(ListSt, int))

monthUub = ['dec', 'jan', 'feb', 'mar', 'abr']

basenameUub = "../aoptimeTraceSelec/uubAoPtimePMT"
outname0 = "../../plots/uubPeaktimeHbAllSt0PMT"+pmtid+".png"
outname1 = "../../plots/uubPeaktimeHbAllSt1PMT"+pmtid+".png"

fadc2mvuub = 2000./4096 # mV

def getRate(st):
  st = str(st)
  xtime = []
  vempkhb = []
  vempkrmshb = []

  for month in monthUub:
    tmpfile = open(basenameUub+pmtid+"St"+st+"Mth"+month+"chpk.json", 'r')
    tmpdata = json.load(tmpfile)
    for info in tmpdata[month]:
      if info['Date'] > 0:
        xtime.append( datetime.datetime.fromtimestamp( info['Date'] ) )
        vempkhb.append( info['Pk']*fadc2mvuub )
        vempkrmshb.append( info['Pkrms']*fadc2mvuub )

  return xtime, vempkhb, vempkrmshb


# ===============================
# *** For Peak HBase Plotting ***

markers = ['o', 's', '^']
xfmt = md.DateFormatter('%m/%d')
mkrSize = 6
fs = 12
ls = 10

ymin = 50
ymax = 120

fig = plt.figure(figsize=(16,9))

pkStations = []
axs = []
tmp = 911
for st in range(0, 9):
  pkStations.append( getRate(stlist[st]) )
  axs.append( fig.add_subplot(tmp) )

  axs[st].errorbar(pkStations[st][0], pkStations[st][1], yerr=pkStations[st][2], fmt=markers[0], markersize=mkrSize, color='tab:orange')
  #axs[st].legend(loc=2, title='PMT'+pmtid, title_fontsize=16, fontsize=16)
  axs[st].xaxis.set_major_formatter(xfmt)
  axs[st].set_ylim(ymin, ymax)
  axs[st].set_ylabel(str(stlist[st])+"\n[mV]", fontsize=fs)
  axs[st].tick_params(axis='y', labelsize=ls, labelright=True)
  axs[st].tick_params(axis='x', labelsize=ls)
  axs[st].grid()
  tmp += 1


fig.suptitle("PMT"+pmtid+" Peak", fontsize=22)
plt.tight_layout()
#plt.show()
plt.savefig(outname0, dpi=100)

fig = plt.figure(figsize=(16,9))
pkStations = []
axs = []
tmp = 911
for st in range(9, 18):
  pkStations.append( getRate(stlist[st]) )
  axs.append( fig.add_subplot(tmp) )

  axs[st-9].errorbar(pkStations[st-9][0], pkStations[st-9][1], yerr=pkStations[st-9][2], fmt=markers[0], markersize=mkrSize, color='tab:orange')
  #axs[st-9].legend(loc=2, title='PMT'+pmtid, title_fontsize=16, fontsize=16)
  axs[st-9].xaxis.set_major_formatter(xfmt)
  axs[st-9].set_ylim(ymin, ymax)
  axs[st-9].set_ylabel(str(stlist[st])+"\n[mV]", fontsize=fs)
  axs[st-9].tick_params(axis='y', labelsize=ls, labelright=True)
  axs[st-9].tick_params(axis='x', labelsize=ls)
  axs[st-9].grid()
  tmp += 1


fig.suptitle("PMT"+pmtid+" Peak", fontsize=22)
plt.tight_layout()
#plt.show()
plt.savefig(outname1, dpi=100)

