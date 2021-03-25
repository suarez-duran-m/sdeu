import ROOT
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

pmtId = sys.argv[1]
stId = sys.argv[2]
pmt = int(pmtId)
uubFile = ROOT.TFile.Open ("../../getCalib/uubCalibHistPMT"+pmtId+"St"+stId+".root", " READ ")
ubFile = ROOT.TFile.Open ("../getCalibNouub/ubCalibHistPMT"+pmtId+"St"+stId+".root" ," READ ")

hUubAP = uubFile.Get("sglSt")
hUbAP = ubFile.Get("sglSt")

uubAp = []
ubAp = []
bins = 0
x = []

for d in range( 1, hUbAP.GetNbinsX()+1 ):
    ubAp.append( hUbAP.GetBinContent(d) )
    x.append(bins)
    bins += 1
ubAp = np.array(ubAp)
ubBins = bins

for d in range( 1, hUubAP.GetNbinsX()+1 ):
    uubAp.append( hUubAP.GetBinContent(d) )
    x.append(bins)
    bins += 1
uubAp = np.array(uubAp)
uubBins = bins

x = np.array(x)

fig1, ax1 = plt.subplots(figsize=(16, 9))
plt.plot(x[:ubBins],ubAp,'r*', label="UB")
plt.plot(x[ubBins:],uubAp,'b*', label="UUB")
plt.ylim(2,11)
plt.xlim(0,)
plt.title("VEM-Ch/VEM-Pk over time for Station "+r"$\bf{"+stId+"}$ "+r"$\bf{LPMT"+pmtId+"}$", size=22.)
plt.xlabel("Time since first event of September 1st, 2020", size=24.)
plt.ylabel("VEM-Ch/VEM-Pk", size=24.)
plt.legend(fontsize=20.)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("../../plots/areaPeakTimePMT"+pmtId+"St"+stId+".pdf", dpi=300)
