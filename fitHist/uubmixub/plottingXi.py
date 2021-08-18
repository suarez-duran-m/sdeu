import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# *** Reading JSON files ***

def getTotals(npmts, basename, ever):
    loctotXiPkHb = []
    loctotXiPkCa = []
    loctotXiChHb = []
    loctotXiChCa = []

    locTotEvt = []
    locStLabel = []

    for pmt in range(0, npmts):
        tmpfile = open(basename+str(pmt+1)+".json", 'r')
        tmpdata = json.load(tmpfile)

        tmptotXiPkHb = []
        tmptotXiPkCa = []
        tmptotXiChHb = []
        tmptotXiChCa = []

        tmptotEvt = []

        for info in tmpdata[ever]:
            tmptotXiPkHb.append( info['extTotXiPkHb'] )
            tmptotXiPkCa.append( info['extTotXiPkCa'] )
            tmptotXiChHb.append( info['extTotXiChHb'] )
            tmptotXiChCa.append( info['extTotXiChCa'] )

            tmptotEvt.append( info['extTotEvet'] )

            if ( pmt==npmts-1 ):
                locStLabel.append( info['stationId'] )

        loctotXiPkHb.append( tmptotXiPkHb )
        loctotXiPkCa.append( tmptotXiPkCa )
        loctotXiChHb.append( tmptotXiChHb )
        loctotXiChCa.append( tmptotXiChCa )

        locTotEvt.append( tmptotEvt )

    return loctotXiPkHb, loctotXiPkCa, loctotXiChHb, loctotXiChCa, locTotEvt, locStLabel


def getRates(npmts, total, ifpk):
    locHbRteXi = []
    locCaRteXi = []
    
    tmprteHb = []
    tmprteCa = []

    tmpHb = 0.
    tmpCa = 0.

    for ipmt in range(0, npmts):
        tmprteHb = []
        tmprteCa = []
        for st in range(0, len(totUub[5])):
            if ifpk:
                tmpHb = total[0][ipmt][st]
                tmpCa = total[1][ipmt][st]
            else:
                tmpHb = total[2][ipmt][st]
                tmpCa = total[3][ipmt][st]
            tmprteHb.append( tmpHb/total[4][ipmt][st] )
            tmprteCa.append( tmpCa/total[4][ipmt][st] )
        locHbRteXi.append( tmprteHb )
        locCaRteXi.append( tmprteCa )

    return locHbRteXi, locCaRteXi


basenameuub = "../fitaopmlt/uubFitXiPMT"
basenameub = "../../nouub/fitHist/fitaopmlt/ubFitXiPMT"

outnmRteXiHbCaUub = "../../plots/uubXiRteXiHbCaPMTs.png"
outnmRteXiHbCaUb = "../../plots/ubXiRteXiHbCaPMTs.png"

npmts = 3

totUub = getTotals(npmts, basenameuub, 'uub')
totUb = getTotals(npmts, basenameub, 'ub')

# =========================================
# *** For Xi Rate Calib/HBase Plotting ***


RteXiPkElec = getRates(npmts, totUub, True)
plotlabel = "UUB LPMT1 "

barWidth = .25

perStPerPmt = []
tmp = []
for pmt in range(0, npmts):
    tmp.append( RteXiPkElec[0][pmt] )
    tmp.append( RteXiPkElec[1][pmt] )
    perStPerPmt.append( tmp )
    tmp = []


plt.figure(figsize=(16,9))

#print(len(perStPerPmt), len(perStPerPmt[0]), len(perStPerPmt[0][0]))
# 3 2 19
width= 2. / len(perStPerPmt[0][0])
Pos = np.array( range( len(perStPerPmt[0][0]) ) )
colors = ['tab:blue', 'tab:orange']
labels = ['HBase', 'Calib']
for bl in range(len(perStPerPmt[0])):
    plt.bar(Pos - (bl+1.5) * width, perStPerPmt[0][bl], width = width, color=colors[bl])
    plt.bar(Pos + bl * width, perStPerPmt[1][bl], width = width, color=colors[bl])
    plt.bar(Pos + (bl+2.5) * width, perStPerPmt[2][bl], width = width, label=labels[bl], color=colors[bl])

plt.legend(fontsize=22)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks([r + width for r in range(len(totUub[5]))], totUub[5], rotation=60)
plt.ylabel("UUB Rate of Succesfull Fit [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim(.15,)
plt.tight_layout()
#plt.grid()
#plt.show()
plt.savefig(outnmRteXiHbCaUub, dpi=100)
print("Relative Diff. For Peak( 1 - XiCalib / XiHBase) OK")

RteXiPkElec = getRates(npmts, totUb, True)
plotlabel = "UB LPMT1 " 

perStPerPmt = []
tmp = []
for pmt in range(0, npmts):
    tmp.append( RteXiPkElec[0][pmt] )
    tmp.append( RteXiPkElec[1][pmt] )
    perStPerPmt.append( tmp )
    tmp = []

plt.figure(figsize=(16,9))
width = 2. / len(perStPerPmt[0][0])
Pos = np.array( range( len(perStPerPmt[0][0]) ) )
for bl in range(len(perStPerPmt[0])):
    plt.bar(Pos - (bl+1.5) * width, perStPerPmt[0][bl], width = width, color=colors[bl])
    plt.bar(Pos + bl * width, perStPerPmt[1][bl], width = width, color=colors[bl])
    plt.bar(Pos + (bl+2.5) * width, perStPerPmt[2][bl], width = width, label=labels[bl], color=colors[bl])

plt.legend(fontsize=22)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks([r + width for r in range(len(totUub[5]))], totUub[5], rotation=60)
plt.ylabel("UB Rate of Succesfull Fit [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim(.15,)
plt.tight_layout()
#plt.grid()
#plt.show()
plt.savefig(outnmRteXiHbCaUb, dpi=100)
print("Relative Diff. For Peak( 1 - XiCalib / XiHBase) OK")

