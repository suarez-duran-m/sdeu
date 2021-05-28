import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

stchoo = int(sys.argv[1])

#monthUb = ['aug', 'sep', 'oct', 'nov']
monthUub = ['aug', 'sep', 'oct', 'nov', 'dec', 'jan', 'feb', 'mar', 'abr']

basenameUub = "../../readMonit/hv"+str(stchoo)
outnameUub = "../../plots/uubHVSt"+str(stchoo)+".png"

# ===============================
# *** For Peak HBase Plotting ***

hvdataUub = [[], [], [], []]# time, pmt1, pmt2, pmt3
for mth in monthUub:
  tmpdata = np.loadtxt( basenameUub+str(mth)+".dat" )
  for info in tmpdata:
    hvdataUub[0].append( datetime.datetime.fromtimestamp( info[0] ) )
    hvdataUub[1].append( info[1] )
    hvdataUub[2].append( info[2] )
    hvdataUub[3].append( info[3] )

fig = plt.figure(figsize=(16,9))

markers = ['o', 's', '^']
ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)
mkrSize = 6

ymin = 9e2
ymax = 2.75e3

plt.plot(hvdataUub[0], hvdataUub[1], markers[0], markersize=mkrSize, label="UUB HV PMT1") #, color='tab:blue'
plt.plot(hvdataUub[0], hvdataUub[2], markers[1], markersize=mkrSize, label="UUB HV PMT2") #, color='tab:blue'
plt.plot(hvdataUub[0], hvdataUub[3], markers[2], markersize=mkrSize, label="UUB HV PMT3") #, color='tab:blue'

plt.legend(fontsize=22)
#plt.ylim(ymin, ymax)
plt.ylabel("Station "+str(stchoo)+" HV [V]", fontsize=22)
plt.xlabel("Time since August 1st, 2020 [month/day]", fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)


plt.tight_layout()
#plt.show()
plt.savefig(outnameUub, dpi=100)
print("Done for HV over time, ST. "+str(stchoo))
