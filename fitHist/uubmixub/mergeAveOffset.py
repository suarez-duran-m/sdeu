import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# *** Reading JSON files ***

pmtId = "PMT" + sys.argv[1]

fuub = open("uubAllAveOffset"+pmtId+".json", 'r')
fub = open("ubAllAveOffset"+pmtId+".json", 'r')
dataUub = json.load(fuub)
dataUb = json.load(fub)

stLabel = []
datamix = {} # Creating dictionary for merge
datamix['uub'] = []
datamix['ub'] = []

xst = []
yuub = []
yuuberr = []
yub = []
yuberr = []

# =======================
# *** Reading for UUB ***
tmpst = 0
tmpave = 0.
tmprms = 0.
for info in dataUub['uub']: # Filling Merge-dictionary
    tmpst = info['stationId']
    tmpave = info['extAveOffset']
    tmprms = info['extRmsOffset']
    datamix['uub'].append({
        'station': tmpst,
        'extAveOffset': tmpave,
        'extRmsOffset': tmprms
        })
    xst.append( tmpst )
    yuub.append( tmpave )
    yuuberr.append( tmprms )


# ======================
# *** Reading for UB ***
tmpst = 0
tmpave = 0.
tmprms = 0.
for info in dataUb['ub']: # Filling Merge-dictionary
    tmpst = info['stationId']
    tmpave = info['extAveOffset']
    tmprms = info['extRmsOffset']
    datamix['ub'].append({
        'station': tmpst,
        'extAveOffset': tmpave,
        'extRmsOffset': tmprms
        })
    yub.append( tmpave )
    yuberr.append( tmprms )

    
json_object = json.dumps(datamix, indent = 4)
with open("dataMix.json", "w") as outfile:
    outfile.write(json_object)



# ====================
# *** For Plotting ***


pmtId = str(1)
plt.figure(figsize=(16,9))
plt.errorbar(xst, yuub, yerr=yuuberr, fmt='o', label="UUB")
plt.errorbar(xst, yub, yerr=yuberr, fmt='x', label="UB")
plt.legend()
plt.title("Average of Offset for Peak Histograms "+pmtId, fontsize=24)
plt.xlabel("Station Id.", fontsize=22)
plt.xticks(rotation=60)
plt.ylabel("Average Offset (FADC)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(225,300)
plt.show()
#plt.savefig("../../plots/bothAllAveOffset"+pmtId+".png", dpi=100)
plt.clf()
print("bothAllAveOffset OK")
