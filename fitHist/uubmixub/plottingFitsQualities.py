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



def fitFunction(x0, par):
  x = 1./x0[0]
  lognormal = np.exp(par[0])*x*np.exp( -0.5*pow( ((np.log(x)+np.log(par[1]))*par[2]),2 ) )
  expo = np.exp( par[3] - par[4]*x )
  return lognormal+expo

monthUub = ['dec'] #, 'jan', 'feb', 'mar', 'abr']

basenameUub = "../aoptimeTraceSelec/uubAoPtimePMT" #+pmtid+"St"+stid+"Mth"
outnameExampResi = "uubPeakExampleResiduals.pdf" #"../../plots/uubPeakChisHbSt"+stid+"PMTs.pdf"


markers = ['o', 's', '^']
xfmt = md.DateFormatter('%m/%d')
mkrSize = 6
fs = 20
ls = 18

chiDistPk = []
chiDistCh = []

month = monthUub[0]

for pmtid in range(1, 2):
  nameElec = "../aoptimeTraceSelec/uubAoPtimePMT"+str(pmtid)+"St"+stid+"Mth"+month+"chpk.root"
  inFile = ROOT.TFile.Open(nameElec, "READ")
  print(nameElec)

  tree = inFile.Get("HistForChi2")
  for evt in range(0, tree.GetEntries()):
    tree.GetEntry(evt)
    chiDistPk.append( tree.pkChi2 )
    chiDistCh.append( tree.chChi2 )

chiDistPk = np.array(chiDistPk)
chiDistCh = np.array(chiDistCh)

tmp = chiDistPk.mean()
pmtid = 1

nameElec = basenameUub+str(pmtid)+"St"+stid+"Mth"+month+"chpk.root"
inFile = ROOT.TFile.Open(nameElec, "READ")
tmptree = inFile.Get("HistForChi2")

# =================================
# *** Plotting Exemple Good Fit ***

c1 = ROOT.TCanvas("c1","c1",1600, 900)
c1.cd()
pad1 = ROOT.TPad("pad1","The pad with the function", 0 ,0.5 ,1 ,1)
pad1.SetFillColor(0)
pad1.SetBottomMargin(0.2)
pad1.Draw()

c1.cd()
pad2 = ROOT.TPad("pad2","The pad with the histogram", 0, 0.05, 1, 0.5)
pad2.SetFillColor(0)
pad2.SetBottomMargin(0.2)
pad2.Draw()

for i in range( len(chiDistPk) ):
  if chiDistPk[i] < tmp*1.2 and chiDistPk[i] > tmp*0.8:
    pad1.cd()
    tmptree.GetEntry(i)
    ROOT.gStyle.SetOptTitle(0)
    tmptree.SetTitle("")
    tmptree.pkHistFit.Draw()
    pad1.Update()

    nPoints = tmptree.pkHistFit.GetN()
    x = 0.
    y = 0.
    #for i in range(0, nPoints):
      #tmptree.pkHistFit.
      #tmptree.pkHistFit.GetPoint(i, x, y)
      #print(x, y)

    c1.Update()
    c1.Print("kk.pdf")
    #c1.Print(outnameExampResi)
    break
