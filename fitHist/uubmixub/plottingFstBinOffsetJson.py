import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# *** Reading JSON files ***

elecver = sys.argv[1]

basename = elecver+"AllAveFstBinOffsetPMT"
outnameFrstBin = "../../plots/"+elecver+"AllAveFrstBinPMTs.png"
outnameOffset = "../../plots/"+elecver+"AllAveOffsetPMTs.png"

stLabel = []
xst = []
yOffsetpmt = []
yOffseterrpmt = []
yFrstBinpmt = []
yFrstBinerrpmt = []

# ========================
# *** Reading for PMTs ***

for pmt in range(0, 3):
    tmpfile = open(basename+str(pmt+1)+".json", 'r')
    tmpdata = json.load(tmpfile)
    tmpaveOffset = []
    tmperrOffset = []
    tmpaveFrstBin = []
    tmperrFrstBin = []
    for info in tmpdata[elecver]: # Filling Merge-dictionary
        tmpaveOffset.append( info['extAveOffset'] )
        tmperrOffset.append( info['extRmsOffset'] )
        tmpaveFrstBin.append( info['extAveFrstBin'] )
        tmperrFrstBin.append( info['extRmsFrstBin'] )
        if ( pmt==2 ):
            xst.append( info['stationId'] )
    yOffsetpmt.append( tmpaveOffset )
    yOffseterrpmt.append( tmperrOffset )
    yFrstBinpmt.append( tmpaveFrstBin )
    yFrstBinerrpmt.append( tmperrFrstBin )


# ==============================
# *** For First Bin Plotting ***

plt.figure(figsize=(16,9))
markers = ['o', 's', '^']
shifted = []
if elecver == 'ub':
    for pmt in range(0, 3):
        tmp = []
        for fstbin in yFrstBinpmt[pmt]:
            tmp.append( fstbin + pmt/10.)
        shifted.append( tmp )
else:
    shifted = yFrstBinpmt

for i in range(0, 3):
    if elecver == 'ub':
        labels = elecver.upper()+" PMT"+str(i+1) + " + " + str(i/10.)
    else:
        labels = elecver.upper()+ " PMT"+str(i+1)
    plt.errorbar(xst, shifted[i], yerr=yFrstBinerrpmt[i], fmt=markers[i], markersize=12, label=labels)
plt.legend(fontsize=20)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Average counts in first bin Peak Histo. [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
if elecver=="uub":
    plt.ylim(-10,1e3)
else:
    plt.ylim(-1,1)
plt.tight_layout()
plt.savefig(outnameFrstBin, dpi=100)
print("bothAllAveFirstBin OK")


# ===========================
# *** For Offset Plotting ***

plt.figure(figsize=(16,9))
for i in range(0, 3):
    plt.errorbar(xst, yOffsetpmt[i], yerr=yOffseterrpmt[i], fmt=markers[i], markersize=12, label=elecver.upper()+" PMT"+str(i+1))
plt.legend(fontsize=20)
#plt.title("Average of Offset for Peak Histograms", fontsize=24)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Average Offset for Peak Histo. [FADC]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.tight_layout()
plt.savefig(outnameOffset, dpi=100)
print("bothAllAveOffset OK")


# ========================
# *** Testing for Bars ***
'''
plt.figure(figsize=(16,9))
width = 0.8 / len(yFrstBinpmt)
Pos = np.array( range( len(xst) ) )
colors = ['tab:blue', 'tab:orange', 'tab:green']

plt.bar(Pos + (0.) * width, yFrstBinpmt[0], yerr=yFrstBinerrpmt[0], width = width, color=colors[0], label=elecver.upper()+" PMT1")
plt.bar(Pos + (1.) * width, yFrstBinpmt[1], yerr=yFrstBinerrpmt[1], width = width, color=colors[1], label=elecver.upper()+" PMT2")
plt.bar(Pos + (2.) * width, yFrstBinpmt[2], yerr=yFrstBinerrpmt[2], width = width, color=colors[2], label=elecver.upper()+" PMT3")

plt.legend(fontsize=22)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks([r + width for r in range(len(xst))], xst, rotation=60)
plt.ylabel("Average counts in first bin Peak Histo. [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0,1e3)
plt.tight_layout()
#plt.grid()
#plt.show()
outnameFrstBin = "../../plots/"+elecver+"AllAveFrstBinPMTsBars.png"
plt.savefig(outnameFrstBin, dpi=100)
print("Relative Diff. For Peak( 1 - XiCalib / XiHBase) OK")
'''
