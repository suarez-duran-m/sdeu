import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

stid = sys.argv[1]
pmtid = sys.argv[2]

monthUb = ['aug', 'sep', 'oct', 'nov']
monthUub = ['dec', 'jan', 'feb', 'mar', 'abr']

basename = "../../getCalib/hv"+stid+"pmt"+pmtid+"web.dat"
outname = "../../plots/uububHv"+stid+"pmt"+pmtid+"web.png"

xtime = [[], []]
yHv = [[], []]

# ========================
# *** Reading for PMTs ***

tmpfile = np.loadtxt(basename)
for i in tmpfile:
  if i[0] < 1607990400: # December 15th
    xtime[0].append(  datetime.datetime.fromtimestamp( i[0] ) )
    yHv[0].append( i[1] )
  else:
    xtime[1].append( datetime.datetime.fromtimestamp( i[0] ) )
    yHv[1].append( i[1] )

markers = ['o', 's', '^']

plt.figure(figsize=(16,9))
plt.plot(xtime[0], yHv[0], markers[0], markersize=12, label="UB")
plt.plot(xtime[1], yHv[1], markers[1], markersize=12, label="UUB")

ax=plt.gca()
xfmt = md.DateFormatter('%m/%d')
ax.xaxis.set_major_formatter(xfmt)

plt.legend(title="St. "+stid+" PMT"+pmtid, fontsize=22, title_fontsize=22)
plt.xlabel("Days since August 1st, 2020 (month/day)", fontsize=24)
plt.ylabel("High Voltage [V]", fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
if stid == '863' and pmtid == "1":
  plt.ylim(1202,1210)
elif stid == '863' and pmtid == "2":
  plt.ylim(1130, 1225)
elif stid == '863' and pmtid == "3":
  plt.ylim(1010, 1030)
elif stid == '1740' and pmtid == "1":
  plt.ylim(1070, 1140)
elif stid == '1740' and pmtid == "2":
  plt.ylim(1190, 1270)
elif stid == '1740' and pmtid == "3":
   plt.ylim(1020, 1060)
plt.tight_layout()
#plt.show()
plt.savefig(outname, dpi=100)
print("Average AoP Peak-HBase OK")

