import ROOT
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime


# ==========================
# *** Reading JSON files ***

stid = sys.argv[1]

monthUb = ['aug', 'sep', 'oct', 'nov']
monthUub = ['dec'] #, 'jan', 'feb', 'mar', 'abr']

basenameUub = "../aoptimeTraceSelec/uubAoPtimePMT" #+pmtid+"St"+stid+"Mth"
basenameUb = "../../nouub/fitHist/aoptimeTraceSelec/ubAoPtimePMT" #+pmtid+"St"+stid+"Mth"
outnameHbase = "../../plots/uubPeaktimeHbSt"+stid+"PMTs.png"
#outchisname = "../../plots/ubPeakChisHbSt"+stid+"PMTs.pdf"
outchisname = "../../plots/uubPeakChisHbSt"+stid+"PMTs.pdf"

fadc2mvub = 2000./1024. # mV
fadc2mvuub = 2000./4096 # mV

def getRate(pmtid):
  pmtid = str(pmtid)
  xtime = [[], []]
  vempkhb = [[], []]
  vempkrmshb = [[], []]

  for month in monthUb:
    tmpfile = open(basenameUb+pmtid+"St"+stid+"Mth"+month+".json", 'r')
    tmpdata = json.load(tmpfile)
    for info in tmpdata[month]:
      xtime[0].append( datetime.datetime.fromtimestamp( info['Date'] ) )
      vempkhb[0].append( info['Pk']*fadc2mvub )
      vempkrmshb[0].append( info['Pkrms']*fadc2mvub )

  for month in monthUub:
    tmpfile = open(basenameUub+pmtid+"St"+stid+"Mth"+month+"chpk.json", 'r')
    tmpdata = json.load(tmpfile)
    for info in tmpdata[month]:
      if info['Date'] > 0:
        xtime[1].append( datetime.datetime.fromtimestamp( info['Date'] ) )
        vempkhb[1].append( info['Pk']*fadc2mvuub )
        vempkrmshb[1].append( info['Pkrms']*fadc2mvuub )

  return xtime, vempkhb, vempkrmshb


# ===============================
# *** For Peak HBase Plotting ***

pkpmt1 = getRate(1)
pkpmt2 = getRate(2)
pkpmt3 = getRate(3)

markers = ['o', 's', '^']
xfmt = md.DateFormatter('%m/%d')
mkrSize = 6
fs = 20
ls = 18

c0 = ROOT.TCanvas("c0","c0",1600, 900) 
c0.cd()
pad1 = ROOT.TPad("pad1","The pad with the function", 0 ,0.65 ,1 ,0.95)
pad1.SetFillColor(0)
pad1.SetBottomMargin(0.2)
pad1.Draw()

c0.cd()
pad2 = ROOT.TPad("pad2","The pad with the histogram", 0, 0.35, 1, 0.65)
pad2.SetFillColor(0)
pad2.SetBottomMargin(0.2)
pad2.Draw()

c0.cd()
pad3 = ROOT.TPad("pad3","The pad with the histogram", 0, 0.01, 1, 0.35)
pad3.SetFillColor(0)
pad3.SetBottomMargin(0.2)
pad3.Draw()

chispkpmt1 = ROOT.TH1F("chispkpmt1", "", 500, 0, 50)
chischpmt1 = ROOT.TH1F("chischpmt1", "", 500, 0, 50)

chispkpmt2 = ROOT.TH1F("chispkpmt2", "", 500, 0, 50)
chischpmt2 = ROOT.TH1F("chischpmt2", "", 500, 0, 50)

chispkpmt3 = ROOT.TH1F("chispkpmt3", "", 500, 0, 50)
chischpmt3 = ROOT.TH1F("chischpmt3", "", 500, 0, 50)

pads = [pad1, pad2, pad3]
chiPkdis = [chispkpmt1, chispkpmt2, chispkpmt3]
chiChdis = [chischpmt1, chischpmt2, chischpmt3]

chiDistPk = []
chiDistCh = []

for pmtid in range(0, 3):
  
  for month in monthUub:
  #for month in monthUb:
    #nameElec = "../../nouub/fitHist/aoptimeTraceSelec/ubAoPtimePMT"+str(pmtid+1)+"St"+stid+"Mth"+month+".root"
    nameElec = "../aoptimeTraceSelec/uubAoPtimePMT"+str(pmtid+1)+"St"+stid+"Mth"+month+"chpk.root"
    inFile = ROOT.TFile.Open(nameElec, "READ")
  
    for bn in range(1, inFile.pkChis.GetXaxis().GetNbins() ):
      chiPkdis[pmtid].Fill( bn/10., inFile.pkChis.GetBinContent( bn ) )
    for bn in range(1, inFile.chChis.GetXaxis().GetNbins()):
      chiChdis[pmtid].Fill( bn/10., inFile.chChis.GetBinContent( bn ) )

  if pmtid==0:
    tree = inFile.Get("HistForChi2")
    for evt in range(0, tree.GetEntries()):
      tree.GetEntry(evt)
      chiDistPk.append( tree.pkChi2 )
      chiDistCh.append( tree.chChi2 )
      #print( evt, tree.pkChi2 )

  pads[pmtid].cd()
  
  chiPkdis[pmtid].SetLineColor(ROOT.kBlue)
  chiPkdis[pmtid].SetMarkerStyle(20)
  chiPkdis[pmtid].SetLineWidth(1)
  chiPkdis[pmtid].SetFillColor(ROOT.kBlue)
  chiPkdis[pmtid].SetFillStyle(3001)
  chiPkdis[pmtid].GetXaxis().SetTitle("#chi^{2}/ndf")
  #chiPkdis[pmtid].GetXaxis().SetRangeUser(0, 5)
  chiPkdis[pmtid].GetXaxis().SetLabelSize(0.06)
  chiPkdis[pmtid].GetXaxis().SetTitleSize(0.06)
  chiPkdis[pmtid].GetYaxis().SetTitle("Counts [au]")
  chiPkdis[pmtid].GetYaxis().SetTitleOffset(0.3)
  chiPkdis[pmtid].GetYaxis().SetRangeUser(0, 370)
  chiPkdis[pmtid].GetYaxis().SetLabelSize (0.06)
  chiPkdis[pmtid].GetYaxis().SetLabelSize(0.06)
  chiPkdis[pmtid].GetYaxis().SetTitleSize(0.06)

  chiPkdis[pmtid].Draw("SAMES HIST")

  chiChdis[pmtid].SetLineColor(810)
  chiChdis[pmtid].SetLineWidth(2)
  chiChdis[pmtid].SetFillColor(810)
  chiChdis[pmtid].SetFillStyle(3001)
  chiChdis[pmtid].Draw("HIST SAMES")

  ptstats = ROOT.TPaveStats(0.6,0.4, 0.70,0.8,"brNDC")
  ptstats.SetTextColor(ROOT.kBlue)
  chiPkdis[pmtid].SetName("PMT "+str(pmtid+1)+" Peak Fit")
  chiPkdis[pmtid].GetListOfFunctions().Add(ptstats)

  ptstats = ROOT.TPaveStats(0.73,0.4,0.83,0.8,"brNDC")
  ptstats.SetTextColor(810)
  chiChdis[pmtid].SetName("PMT "+str(pmtid+1)+" Charge Fit")
  chiChdis[pmtid].GetListOfFunctions().Add(ptstats)
  
  pads[pmtid].Update()

c0.Update()
c0.Print(outchisname)


chiDistPk = np.array(chiDistPk)
chiDistCh = np.array(chiDistCh)

tmp = chiDistPk.mean()

c1 = ROOT.TCanvas("c1","c1",1600, 900)
nameElec = "../aoptimeTraceSelec/uubAoPtimePMT1St"+stid+"Mth"+month+"chpk.root"
inFile = ROOT.TFile.Open(nameElec, "READ")
tmptree = inFile.Get("HistForChi2")

ok = 1
ok1 = 1

c1.cd()
tmptree.GetEntry(22)
tmptree.pkHistFit.SetLineColor(ROOT.kBlue)
tmptree.pkHistFit.Draw("SAME AL*")
tmptree.GetEntry(4)
tmptree.pkHistFit.Draw("SAME AL*")

'''
for i in range( len(chiDistPk) ):
  if chiDistPk[i] < tmp*1.2 and chiDistPk[i] > tmp*0.8 and ok:
    c1.Clear()
    c1.cd()
    tmptree.GetEntry(i)
    tmptree.pkHistFit.Draw()
    ok = 0
    print(i)
    #c1.Print("kk0.pdf")

  if chiDistPk[i] > 20. and chiDistPk[i] < 40 and ok1:
    tmptree.GetEntry(i)
    #c1.Clear()
    #c1.cd()
    tmptree.GetEntry(i)
    tmptree.pkHistFit.SetLineColor(ROOT.kBlue)
    tmptree.pkHistFit.Draw("SAMES")
    ok1 = 0
    print(i)
    #c1.Print("kk1.pdf")
'''
c1.Print("kk1.pdf")
ymin = 40
ymax = 100

index = [i for i in range(87,106)]
a = np.array(pkpmt1[0][0])
b = np.delete(a, index)
a = np.array(pkpmt1[1][0])
c = np.delete(a, index)
a = np.array(pkpmt1[2][0])
d = np.delete(a, index)

fig = plt.figure(figsize=(16,9))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

ax1.errorbar(b, c, yerr=d, fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax1.errorbar(pkpmt1[0][1], pkpmt1[1][1], yerr=pkpmt1[2][1], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax1.legend(loc=2, title='PMT1', title_fontsize=16, fontsize=16)
ax1.xaxis.set_major_formatter(xfmt)
ax1.set_ylim(ymin, ymax)
ax1.set_ylabel("Peak [mV]", fontsize=fs)
ax1.tick_params(axis='y', labelsize=ls, labelright=True)
ax1.tick_params(axis='x', labelsize=ls)

print("\n======================\n")
print("=== PMT1 ===\n")
print("UB:", np.mean(c))
a = np.array(pkpmt1[1][1])
print("UUB:", a[a>0].mean())


a = np.array(pkpmt2[0][0])
b = np.delete(a, index)
a = np.array(pkpmt2[1][0])
c = np.delete(a, index)
a = np.array(pkpmt2[2][0])
d = np.delete(a, index)

ax2.errorbar(b, c, yerr=d, fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax2.errorbar(pkpmt2[0][1], pkpmt2[1][1], yerr=pkpmt2[2][1], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax2.legend(loc=0, title='PMT2', title_fontsize=16, fontsize=16)
ax2.xaxis.set_major_formatter(xfmt)
ax2.set_ylim(ymin, ymax)
ax2.set_ylabel("Peak [mV]", fontsize=fs)
ax2.tick_params(axis='both', labelsize=ls)
ax2.tick_params(axis='y', labelsize=ls, labelright=True)

print("\n======================\n")
print("=== PMT2 ===\n")
print("UB:", np.mean(c))
a = np.array(pkpmt2[1][1])
print("UUB:", a[a>80].mean())

a = np.array(pkpmt2[0][0])
b = np.delete(a, index)
a = np.array(pkpmt2[1][0])
c = np.delete(a, index)
a = np.array(pkpmt2[2][0])
d = np.delete(a, index)

ax3.errorbar(b, c, yerr=d, fmt=markers[0], markersize=mkrSize, color='tab:blue', label="UB")
ax3.errorbar(pkpmt3[0][1], pkpmt3[1][1], yerr=pkpmt3[2][1], fmt=markers[1], markersize=mkrSize, color='tab:orange', label="UUB")
ax3.legend(loc=2, title='PMT3', title_fontsize=16, fontsize=16)
ax3.xaxis.set_major_formatter(xfmt)
ax3.set_xlabel("Days since August 1st, 2020 (month/day)", fontsize=fs)
ax3.set_ylim(ymin, ymax)
ax3.set_ylabel("Peak [mV]", fontsize=fs)
ax3.tick_params(axis='both', labelsize=ls)
ax3.tick_params(axis='y', labelsize=ls, labelright=True)

fig.suptitle("Station "+stid, fontsize=22)
plt.tight_layout()
#plt.show()
plt.savefig(outnameHbase, dpi=100)

print("\n======================\n")
print("=== PMT3 ===\n")
print("UB:", np.mean(c))
a = np.array(pkpmt3[1][1])
print("UUB:", a[a>70].mean(),"\n")

print("Done for Peak over time, ST. "+stid)
