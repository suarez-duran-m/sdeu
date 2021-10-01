import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# *** Reading JSON files ***

elecver = sys.argv[1]

basename = elecver+"AllAveOffsetChPMT"
outnameOffset = "../../plots/"+elecver+"AllAveOffsetChPMTs.png"
outnameOffset3 = "../../plots/"+elecver+"AllAveOffsetChPMT3s.png"

stLabel = []
xst = []
yOffsetpmt = []
yOffseterrpmt = []

# ========================
# *** Reading for PMTs ***

for pmt in range(0, 3):
    tmpfile = open(basename+str(pmt+1)+".json", 'r')
    tmpdata = json.load(tmpfile)
    tmpaveOffset = []
    tmperrOffset = []
    for info in tmpdata[elecver]: # Filling Merge-dictionary
        tmpaveOffset.append( info['extAveOffset'] )
        tmperrOffset.append( info['extRmsOffset'] )
        if ( pmt==2 ):
            xst.append( info['stationId'] )
    yOffsetpmt.append( tmpaveOffset )
    yOffseterrpmt.append( tmperrOffset )

# ===========================
# *** For Offset Plotting ***
npmts = 3
markers = ['o', 's', '^']
shifted = []
if elecver=='uub':
    npmts = 2
    for pmt in range(0, npmts):
        tmp = []
        for fstbin in yOffsetpmt[pmt]:
            tmp.append( fstbin + pmt/10.)
        shifted.append( tmp )
else:
    shifted = yOffsetpmt

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    if elecver == 'uub':
        labels = elecver.upper()+" PMT"+str(i+1) + " + " + str(i/10.)
    else:
        labels = elecver.upper()+" PMT"+str(i+1)
    plt.errorbar(xst, shifted[i], yerr=yOffseterrpmt[i], fmt=markers[i], markersize=12, label=labels)
plt.legend(fontsize=22)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("Average Offset for Charge Histo. (FADC)", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
if elecver == 'uub':
    plt.ylim(-0.05, 0.15)
plt.tight_layout()
plt.savefig(outnameOffset, dpi=100)
print("bothAllAveOffsetCharge OK")

if elecver=='uub':
    plt.figure(figsize=(16,9))
    plt.errorbar(xst, yOffsetpmt[2], yerr=yOffseterrpmt[2], fmt=markers[2], markersize=12, c='green', label=elecver.upper()+" PMT3")
    plt.legend(fontsize=22)
    plt.xlabel("Station Id.", fontsize=24)
    plt.xticks(rotation=60)
    plt.ylabel("Average Offset for Charge Histo. (FADC)", fontsize=24)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.tight_layout()
    plt.savefig(outnameOffset3, dpi=100)
    print("bothAllAveOffsetCharge3 OK")
