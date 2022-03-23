import ROOT
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

pmtId = sys.argv[1]
pmt = int(pmtId)

inFile = ROOT.TFile.Open ( "uubCalibHistPMT"+pmtId+".root" ," READ ")

ch = inFile.Get("hCh")
pk = inFile.Get("hPk")
ap = inFile.Get("hap")

arr0 = []
chH = []
pkH = []
apH = []

for d in range( 1, ch.GetNbinsX()+1 ):
    arr0 = []
    for st in range( 1, ch.GetNbinsY()+1 ):
        if ch.GetBinContent(d, st) != 0:
            arr0.append( ch.GetBinContent(d, st) )
        else:
            arr0.append( np.nan );
    chH.append( arr0 )

chH = np.array( chH )

cmin = 1.3e3
cmax = 1.9e3

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]


fig1, ax1 = plt.subplots(figsize=(16, 9))
chH = np.ma.masked_invalid( chH )

heatmap = ax1.imshow(chH.T, cmap='RdYlBu', vmin=cmin, vmax=cmax, aspect='auto')
ax1.patch.set(hatch='///', edgecolor='black', fill=False, snap=False)

plt.title( ch.GetTitle(), size=24. )
plt.xlabel("Days since December 1st 2020", size=24.)
plt.ylabel("Station ID", size=24.)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

ax1.set_yticks(np.arange(chH.shape[1]), minor=False)
ax1.set_yticklabels(stLabel, minor=False)

divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar = fig1.colorbar(heatmap, cax=cax)
cbar.ax.tick_params(labelsize=18)
cbar.set_label("Position of VEM-Charge / FADC", size=18) 

plt.savefig("../plots/uuBchargePMT"+pmtId+"Hg.pdf", dpi=300)

print("Charge for PMt"+pmtId+" done")

for d in range( 1, pk.GetNbinsX()+1 ):
    arr0 = []
    for st in range( 1, pk.GetNbinsY()+1 ):
        if pk.GetBinContent(d, st) != 0:
            arr0.append( pk.GetBinContent(d, st) )
        else:
            arr0.append( np.nan );
    pkH.append( arr0 )

pkH = np.array( pkH )

cmin = 1.5e2
cmax = 2.1e2

fig2, ax2 = plt.subplots(figsize=(16, 9))
heatmap = ax2.imshow(pkH.T, cmap='RdYlBu', vmin=cmin, vmax=cmax, aspect='auto')
ax2.patch.set(hatch='///', edgecolor='black', fill=False, snap=False)

plt.title( pk.GetTitle(), size=24. )
plt.xlabel("Days since December 1st 2020", size=24.)
plt.ylabel("Station ID", size=24.)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

ax2.set_yticks(np.arange(pkH.shape[1]), minor=False)
ax2.set_yticklabels(stLabel, minor=False)

divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar = fig2.colorbar(heatmap, cax=cax)
cbar.ax.tick_params(labelsize=18)
cbar.set_label("Position of VEM-Peak / FADC", size=18) 

plt.savefig("../plots/uuBpeakPMT"+pmtId+"Hg.pdf", dpi=300)

print("Peak for PMt"+pmtId+" done")

for d in range( 1, ap.GetNbinsX()+1 ):
    arr0 = []
    for st in range( 1, ap.GetNbinsY()+1 ):
        if ap.GetBinContent(d, st) != 0:
            arr0.append( ap.GetBinContent(d, st) )
        else:
            arr0.append( np.nan );
    apH.append( arr0 )

apH = np.array( apH )

cmin = 8.0
cmax = 9.4

fig3, ax3 = plt.subplots(figsize=(16, 9))
heatmap = ax3.imshow(apH.T, cmap='RdYlBu', vmin=cmin, vmax=cmax, aspect='auto')
ax3.patch.set(hatch='///', edgecolor='black', fill=False, snap=False)

plt.title( ap.GetTitle(), size=24. )
plt.xlabel("Days since December 1st 2020", size=24.)
plt.ylabel("Station ID", size=24.)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

ax3.set_yticks(np.arange(pkH.shape[1]), minor=False)
ax3.set_yticklabels(stLabel, minor=False)

divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar = fig3.colorbar(heatmap, cax=cax)
cbar.ax.tick_params(labelsize=18)
cbar.set_label("A/P / FADC", size=18) 

plt.savefig("../plots/uuBapPMT"+pmtId+"Hg.pdf", dpi=300)

print("A/P for PMt"+pmtId+" done")
