import numpy as np
from scipy.stats import moyal
import matplotlib.pyplot as plt

x25 = np.linspace(25, 400, int(400/25))
x833 = np.linspace(8.33, 400, int(400/8.33))
x1 = np.linspace(1, 400, 400)

def searchPk(x, y):
    tmp = 0
    pk = 0
    for i in range(len(x)):
        tmp = y[i]
        if tmp > pk:
            pk = tmp
    return pk

def getArea(x, y):
    summ = 0
    for i in range(len(x)):
        summ += y[i]
    return summ

pulses = np.loadtxt("rmdpulses.dat", float)

xtime = pulses.T[0]
yampl = pulses.T[1]

# =======================
# *** Plotting pulses ***
'''
plt.figure(figsize=(16,9))

plt.plot(xtime, yampl, 'o')

plt.ylabel("Amplitude [mV]", fontsize=24)
plt.yticks(fontsize=22)
plt.xlabel("time [ns]", fontsize=24)
plt.xticks(fontsize=22)

plt.tight_layout()
#plt.savefig("../plots/artifitialPulses.png", dpi=100)
plt.show()
'''

sglpulse210 = []
sglpulse212 = []
mv2fadc210 = 1.95 #mV
mv2fadc212 = 0.49 #mV

bintime = [i for i in range(40000)]

tmpbin25 = 0.
tmpcnt = 0

for np in range( int(yampl.size/4e4) ):
  tmp = []
  tmpbin25 = 25.
  for t in bintime:
    if t == int(tmpbin25*1e2):
      tmp.append( int(yampl[tmpcnt]/mv2fadc210) )
      tmpbin25 += 25.
    tmpcnt += 1
  sglpulse210.append( tmp )


tmpbin833 = 0.
tmpcnt = 0
for np in range( int(yampl.size/4e4) ):
  tmp = []
  tmpbin833 = 8.33
  for t in bintime:
    if t == int(1e2*tmpbin833):
      tmp.append( int(yampl[tmpcnt]/mv2fadc212) )
      tmpbin833 += 8.33
    tmpcnt += 1
  sglpulse212.append( tmp )

import numpy as np
pkhist25 = np.zeros(200, int)
pkhist833 = np.zeros(200, int)

chhist25 = np.zeros(400, int)
chhist833 = np.zeros(400, int)

factorCh25 = 25./50.
factorCh833 = 8.33/50.

for i in range(len(sglpulse210)):
  pkhist25[ int(mv2fadc210*searchPk(x25[:-1], sglpulse210[i])) ] += 1
  chhist25[ int((mv2fadc210*factorCh25)*getArea(x25[:-1], sglpulse210[i])) ] += 1

for i in range(len(sglpulse212)):
  pkhist833[ int(mv2fadc212*searchPk(x833[:-1], sglpulse212[i])) ] += 1
  chhist833[ int((mv2fadc212*factorCh833)*getArea(x833[:-1], sglpulse212[i])) ] += 1

plt.figure(figsize=(16,9))

plt.plot(pkhist25, 's', ms=12, label="Peak for Sample/25 ns")
plt.plot(pkhist833, '^', ms=12, label="Peak for Sample/8.33 ns")

plt.legend(fontsize=24)
plt.ylabel("Counts [au]", fontsize=24)
plt.yticks(fontsize=22)
plt.xlabel("Amplitude [mV]", fontsize=24)
plt.xticks(fontsize=22)
plt.xlim(50,80)

plt.tight_layout()
plt.savefig("../plots/samplingDigiPkHistos.png", dpi=100)
#plt.show()


plt.figure(figsize=(16,9))

plt.plot(chhist25, 's', ms=12, label="Charge for Sample/25 ns")
plt.plot(chhist833, '^', ms=12, label="Charge for Sample/8.33 ns")

plt.legend(fontsize=24)
plt.ylabel("Counts [au]", fontsize=24)
plt.yticks(fontsize=22)
plt.xlabel("Charge [pC]", fontsize=24)
plt.xticks(fontsize=22)
plt.xlim(80,140)

plt.tight_layout()
plt.savefig("../plots/samplingDigiChHistos.png", dpi=100)
#plt.show()
