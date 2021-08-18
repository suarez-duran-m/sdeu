import ROOT
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

pmtId = sys.argv[1]

inFile = ROOT.TFile.Open ( "zooTraces100binsPMT"+pmtId+".root" ," READ ")

hist = inFile.Get ("okDiffRmsh")
arr0 = []
arr1 = []

for d in range( 1, hist.GetNbinsX()+1 ):
    arr0 = []
    for st in range( 1, hist.GetNbinsY()+1 ):
        if hist.GetBinContent(d, st) != 0:
            arr0.append( hist.GetBinContent(d, st) )
        else:
            arr0.append( np.nan );
    arr1.append( arr0 )

arr1 = np.array( arr1 )

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]


fig, ax = plt.subplots(figsize=(16, 9))
arr1 = np.ma.masked_invalid(arr1)

heatmap = ax.imshow(arr1.T, cmap='RdYlBu', vmin=-2, vmax=2, aspect='auto')
ax.patch.set(hatch='///', edgecolor='black', fill=False, snap=False)

plt.title( hist.GetTitle(), size=24. )
plt.xlabel("Days since December 1st 2020", size=24.)
plt.ylabel("Station ID", size=24.)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

ax.set_yticks(np.arange(arr1.shape[1]), minor=False)
ax.set_yticklabels(stLabel, minor=False)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar = fig.colorbar(heatmap, cax=cax)
cbar.ax.tick_params(labelsize=18)
cbar.set_label("Diff. RMS / FADC", size=18) 

plt.savefig("trOkPMT"+pmtId+"Hg.pdf", dpi=300)
#plt.show()
