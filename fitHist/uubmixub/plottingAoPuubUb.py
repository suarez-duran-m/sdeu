import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# *** Reading JSON files ***

# *** Reading for PMTs ***

def getAoPJson(filename, elecver):
    xst = []
    yaopHbpmt = []
    yaopHberrpmt = []
    yaopCapmt = []
    yaopCaerrpmt = []

    ydiffHbCa = []
    npmts = 3

    for pmt in range(0, npmts):
        tmpfile = open(filename+str(pmt+1)+".json", 'r')
        print("opened:", filename+str(pmt+1)+".json")
        tmpdata = json.load(tmpfile)
        tmpaveAoPHb = []
        tmperrAoPHb = []
        tmpaveAoPCa = []
        tmperrAoPCa = []
        for info in tmpdata[elecver]:
            tmpaveAoPHb.append( info['extAveAoPHb'] )
            tmperrAoPHb.append( info['extRmsAoPHb'] )
            tmpaveAoPCa.append( info['extAveAoPCa'] )
            tmperrAoPCa.append( info['extRmsAoPCa'] )

            if ( pmt==npmts-1 ):
                xst.append( info['stationId'] )
        yaopHbpmt.append( tmpaveAoPHb )
        yaopHberrpmt.append( tmperrAoPHb )
        yaopCapmt.append( tmpaveAoPCa )
        yaopCaerrpmt.append( tmperrAoPCa )

    return xst, yaopHbpmt, yaopHberrpmt, yaopCapmt, yaopCaerrpmt


basenameUub = "../fitaopmlt/uubFitaopPMT"
basenameUb = "../../nouub/fitHist/fitaopmlt/ubFitaopPMT"

outnameHbase = "../../plots/aopHbaseUubUbPMTs.png"
outnameCalib = "../../plots/aoCalibUubUbPMTs.png"

aopUub = getAoPJson(basenameUub, "uub")
aopUb = getAoPJson(basenameUb, "ub")
stLabel = []
npmts = 3

ubouub = 8.33/25.
markers = ['o', 's', '^']

# ==============================
# *** For AoP HBase Plotting ***
diffHb = []
tmpaopUub = np.array( aopUub[1] )
tmpaopUb = np.array( aopUb[1] )

for pmt in range(0, npmts):
    tmp = []
    for st in range(len(aopUub[0])):
        if tmpaopUub[pmt][st] > 4 and tmpaopUb[pmt][st] > 2:
            tmp.append( 1 - ( ubouub*tmpaopUub[pmt][st]/tmpaopUb[pmt][st] ) )
        else:
            tmp.append(1)
    diffHb.append(tmp)

yerr = []
for pmt in range(0, npmts):
    tmp = []
    for st in range( len(aopUub[0]) ):
        if aopUb[1][pmt][st] > 0:
            tmp.append( ( aopUub[2][pmt][st]/aopUb[1][pmt][st] )**2 + ( aopUb[2][pmt][st]*aopUub[1][pmt][st]/aopUb[1][pmt][st]**2 )**2 )
        else:
            tmp.append( 0 )
    yerr.append( np.sqrt( tmp ) )

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.errorbar(aopUub[0], diffHb[i], yerr=yerr[i], fmt=markers[i], markersize=12, label=" PMT"+str(i+1))

plt.legend(fontsize=22, loc=1)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("(1 - (0,33*AoPuub)/AopUb) HBase [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0, 0.6)
#plt.grid(axis='y', color='0.75')
#plt.show()
plt.tight_layout()
plt.savefig(outnameHbase, dpi=100)
print("Average AoP Peak-HBase OK")

tmpsum = []
tmpsum.append( [np.average(diffHb[pmt]) for pmt in range(3)] )
print(np.average(np.array(tmpsum)))

# ==============================
# *** For AoP Calib Plotting ***

diffHb = []
tmpaopUub = np.array( aopUub[3] )
tmpaopUb = np.array( aopUb[3] )

for pmt in range(0, npmts):
    tmp = []
    for st in range(len(aopUub[0])):
        if tmpaopUub[pmt][st] > 4 and tmpaopUb[pmt][st] > 2:
            tmp.append( 1 - ( ubouub*tmpaopUub[pmt][st]/tmpaopUb[pmt][st] ) )
        else:
            tmp.append(1)
    diffHb.append(tmp)

yerr = []
for pmt in range(0, npmts):
    tmp = []
    for st in range( len(aopUub[0]) ):
        if aopUb[1][pmt][st] > 0:
            tmp.append( ( aopUub[4][pmt][st]/aopUb[3][pmt][st] )**2 + ( aopUb[4][pmt][st]*aopUub[3][pmt][st]/aopUb[3][pmt][st]**2 )**2 )
        else:
            tmp.append( 0 )
    yerr.append( np.sqrt( tmp ) )

plt.figure(figsize=(16,9))
for i in range(0, npmts):
    plt.errorbar(aopUub[0], diffHb[i], yerr=yerr[i], fmt=markers[i], markersize=12, label=" PMT"+str(i+1))
    #plt.plot(aopUub[0], diffHb[i], markers[i], markersize=12, label=" PMT"+str(i+1))
plt.legend(fontsize=22, loc=1)
plt.xlabel("Station Id.", fontsize=24)
plt.xticks(rotation=60)
plt.ylabel("(1 - (0,33*AoPuub)/AopUb) Calib [au]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0,0.6)
#plt.grid(axis='y', color='0.75')
#plt.show()
plt.tight_layout()
plt.savefig(outnameCalib, dpi=100)
print("Average AoP Peak-Calib OK")
