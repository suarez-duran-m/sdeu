import ROOT
import sys
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable

def getMean( lst ):
    mean = 0.
    for i in lst:
        if i < 10:
            mean += i
    return mean/(1.*len(lst))

def getRms( lst, mn ):
    rms = 0.
    for i in lst:
        if i < 10:
            rms += (i - mn)*(i - mn)
    return np.sqrt( rms/(1.*len(lst)) )


pmtId = sys.argv[1]
pmt = int(pmtId)

inFile = ROOT.TFile.Open ( "uubCalibHistPMT"+pmtId+".root" ," READ ")

ap = inFile.Get("apDist")
tmp = 0.

stLabel = ["863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"]

apMea = np.zeros( len(stLabel) ) # Mean per station
apRms = np.zeros( len(stLabel) ) # RMS per station
arr0 = []

for st in range( 1, ap.GetNbinsX()+1 ):
    arr0 = []
    for evt in range( 1, ap.GetNbinsY()+1 ):
        if ap.GetBinContent(st, evt) != 0:
            arr0.append( ap.GetBinContent(st, evt) )
    if len(arr0) > 0:
        tmp = getMean( arr0 )
        apMea[st-1] = tmp
        apRms[st-1] = getRms(arr0, tmp)
    else:
        apMea[st-1] = 0.
        apRms[st-1] = 0.

print(apMea)

x = []
y = []
y1 = []

for i in range( 1, ap.GetNbinsY()+1 ):
    x.append(i)
    y.append( ap.GetBinContent(1, i) )
    y1.append( ap.GetBinContent(18, i) )

cmin = 7.5
cmax = 9.5
n=0
for i in stLabel:
    print(n,i)
    n+=1
print(apMea[0], apRms[0])
print(apMea[17], apRms[17])
print(np.amax(y), np.amax(y1))

fig1 = plt.subplots(figsize=(16, 9))
plt.errorbar(stLabel, apMea, yerr=apRms, marker='o', linestyle='none')
#plt.plot(x, y, 'o')
#plt.plot(x, y1, '*')
#plt.axhline(y=apMea[0], c='r', ls='--')


#plt.title( ap.GetTitle(), size=24. )
plt.xticks(rotation=45)
plt.xlabel("Stations", size=24.)
plt.ylabel("A/P", size=24.)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(cmin, cmax)

plt.show()

#plt.savefig("../plots/uuBapPMT"+pmtId+"Hg.pdf", dpi=300)
