import ROOT
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ====================
# =======      =======
# *** *** MAIN *** ***
# ====================

pmtId = "PMT" + sys.argv[1]
nameElec = "uubAoPsgl"+pmtId+".root"
inFile = ROOT.TFile.Open(nameElec, "READ")
print("File:", nameElec, "opened.")

histos = inFile.Get("Histograms")

histos.GetEntry(0)

canvas = ROOT.TCanvas("canvas")
canvas.cd()
'''
# ============================
# *** Plot for Smooth-Peak ***
histos.recePk.SetTitle("")
histos.recePk.SetStats(0)
histos.recePk.Draw("SAME")
histos.recePk.GetYaxis().SetTitle("Counts [au]")
histos.recePk.GetYaxis().SetTitleOffset(1.2)
histos.recePk.GetXaxis().SetTitle("Bin [FADC]")

histos.smoothPk.SetLineColor(2)
histos.smoothPk.SetTitle("")
histos.smoothPk.SetStats(0)
histos.smoothPk.Draw("SAME")

legend = ROOT.TLegend(0.52,0.6,0.88,0.8)
legend.AddEntry(histos.recePk, "Original Peak histogram", "l")
legend.AddEntry(histos.smoothPk, "Smoothed Peak histogram", "l")
ROOT.gStyle.SetLegendFont(32)
ROOT.gStyle.SetLegendTextSize(0.035)
legend.Draw()
canvas.Print("../../plots/uubSmoothOriginalPk.pdf")
'''

'''
# ==============================
# *** Plot for Smooth-Charge ***
histos.receCh.SetTitle("")
histos.receCh.SetStats(0)
histos.receCh.Draw("SAME")
histos.receCh.GetYaxis().SetTitle("Counts [au]")
histos.receCh.GetYaxis().SetTitleOffset(1.3)
histos.receCh.GetXaxis().SetTitle("Bin [FADC]")

histos.smoothCh.SetLineColor(2)
histos.smoothCh.SetTitle("")
histos.smoothCh.SetStats(0)
histos.smoothCh.Draw("SAME")

legend = ROOT.TLegend(0.52,0.65,0.89,0.85)
legend.AddEntry(histos.receCh, "Original Charge histogram", "l")
legend.AddEntry(histos.smoothCh, "Smoothed Charge histogram", "l")
ROOT.gStyle.SetLegendFont(32)
ROOT.gStyle.SetLegendTextSize(0.035)
ROOT.gPad.SetLogy()
legend.Draw()
canvas.Print("../../plots/uubSmoothOriginalCh.pdf")
'''


