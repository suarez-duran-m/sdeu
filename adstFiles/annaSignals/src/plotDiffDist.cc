#include "plotDiffDist.h"
#include <iostream>

#include "TF1.h"
using namespace std;

plotDiffDist::plotDiffDist(TString stId, TString printPath, TH1D *pmt12, TH1D *pmt13, 
    TH1D *pmt23, TH1D *totSignalBef, TH1D *totSignalAft) {

  stName = stId;
  outputPath = printPath;

  outputFile = new TFile(outputPath+stId+".root", "RECREATE");

  diff12 = pmt12;
  diff13 = pmt13;
  diff23 = pmt23;

  totSglBef = totSignalBef;
  totSglAft = totSignalAft;

  canvasDiff = doCanvas(stId+"Diff");
  canvasSignal = doCanvas(stId+"Signal");
  canvasTaus = doCanvas(stId+"Taus");
  doDiffDistPlot();
  doDisTotSglPlot();
}

void plotDiffDist::plotTausMuon(TH1D *tauBefPmt1, TH1D *tauBefPmt2, 
    TH1D *tauBefPmt3, TH1D *tauAftPmt1, TH1D *tauAftPmt2, TH1D *tauAftPmt3) {

  muonTauBefPmt1 = tauBefPmt1;
  muonTauBefPmt2 = tauBefPmt2;
  muonTauBefPmt3 = tauBefPmt3;

  muonTauAftPmt1 = tauAftPmt1;
  muonTauAftPmt2 = tauAftPmt2;
  muonTauAftPmt3 = tauAftPmt3;
  TString pmtId;

  doDisTausPlot(pmtId="PMT1", muonTauBefPmt1, muonTauAftPmt1);
  doDisTausPlot(pmtId="PMT2", muonTauBefPmt2, muonTauAftPmt2);
  doDisTausPlot(pmtId="PMT3", muonTauBefPmt3, muonTauAftPmt3);
}

TCanvas *plotDiffDist::doCanvas(TString canvasName) {
  TCanvas *canvas = new TCanvas(canvasName, canvasName, 102, 76, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetRightMargin(0.017);
  canvas->SetLeftMargin(0.074);
  canvas->SetTopMargin(0.014);
  canvas->SetBottomMargin(0.1);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderMode(0);

  return canvas;
}

void plotDiffDist::doDiffDistPlot() {
  canvasDiff->cd();

  int max = diff12->GetMaximum() > diff13->GetMaximum() 
    ? diff12->GetMaximum() : diff13->GetMaximum();
  max = diff23->GetMaximum() > max ? diff23->GetMaximum() : max;
 
  diff12->SetStats(kFALSE); 
  diff12->GetXaxis()->SetTitle("(S_{pmt_i} - S_{pmt_j}) / ErrS^{Tot}");
  diff12->GetXaxis()->SetTitleOffset(1.3);
  diff12->GetXaxis()->SetLabelSize(0.05);
  diff12->GetXaxis()->SetTitleSize(0.05);
  diff12->GetXaxis()->SetTitleOffset(1.);
  diff12->GetYaxis()->SetRangeUser(0, 1.1*max);
  diff12->GetYaxis()->SetTitle("Counts [au]");
  diff12->GetYaxis()->SetTitleSize(0.05);
  diff12->GetYaxis()->SetLabelSize(0.7);
  diff12->SetLineColor(kRed);
  diff12->Draw();
  diff13->SetLineColor(kBlue);
  diff13->Draw("same");
  diff23->SetLineColor(kGreen+3);
  diff23->Draw("same");

  doLegendDiff();

  canvasDiff->Print(outputPath+"diffDist"+stName+".pdf");
}

void plotDiffDist::doDisTausPlot(TString pmtId, TH1D *tausBef, TH1D *tausAft) {
  canvasTaus->cd();
  canvasTaus->SetLogy();

  int max = tausBef->GetMaximum() > tausAft->GetMaximum() 
    ? tausBef->GetMaximum() : tausAft->GetMaximum();
 
  tausBef->SetStats(kFALSE); 
  tausBef->GetXaxis()->SetTitle("#tau [ns]");
  tausBef->GetXaxis()->SetLabelSize(0.05);
  tausBef->GetXaxis()->SetTitleSize(0.05);
  tausBef->GetXaxis()->SetTitleOffset(1.);
  tausBef->GetYaxis()->SetRangeUser(1e-1, 1.1*max);
  tausBef->GetYaxis()->SetTitle("Counts [au]");
  tausBef->GetYaxis()->SetLabelSize(0.05);
  tausBef->GetYaxis()->SetTitleSize(0.05);
  tausBef->GetYaxis()->SetTitleOffset(.7);
  tausBef->SetLineColor(kRed);
  tausBef->Draw();

  tausAft->SetLineColor(kBlue);
  tausAft->Draw("same");

  doLegendTaus(pmtId, tausBef, tausAft);

  canvasTaus->Print(outputPath+"distTaus"+stName+pmtId+".pdf");
}

void plotDiffDist::doDisTotSglPlot() {
  canvasSignal->cd();
  canvasSignal->SetLeftMargin(0.094);

  TF1 *powerLow = new TF1 ("powerLow", "[0]*x^(-[1])",4, 18); 

  totSglBef->Sumw2(true);
  totSglBef->Scale(1./totSglBef->Integral());
  totSglBef->Fit("powerLow", "QR");
  totSglAft->Sumw2(true);
  totSglAft->Scale(1./totSglAft->Integral());
  totSglAft->Fit("powerLow", "QR");

  totSglBef->GetFunction("powerLow")->SetLineColor(kBlue);
  totSglBef->GetFunction("powerLow")->SetLineWidth(1);
  totSglBef->SetStats(kFALSE);
  totSglBef->GetYaxis()->SetRangeUser(1e-3, 6e-1);
  totSglBef->GetYaxis()->SetTitle("Counts [au]");
  totSglBef->GetYaxis()->SetTitleSize(0.05);
  totSglBef->GetYaxis()->SetTitleOffset(0.9);
  totSglBef->GetYaxis()->SetLabelSize(0.05);
  totSglBef->GetXaxis()->SetRangeUser(2,30);
  totSglBef->GetXaxis()->SetTitle("S [VEM]");
  totSglBef->GetXaxis()->SetTitleSize(0.05);
  totSglBef->GetXaxis()->SetTitleOffset(1.);
  totSglBef->GetXaxis()->SetLabelSize(0.05);
  totSglBef->SetLineColor(kBlue);
  canvasSignal->SetLogy();
  canvasSignal->SetLogx();
  totSglBef->Draw("E1");

  totSglAft->GetFunction("powerLow")->SetLineColor(kRed);
  totSglAft->GetFunction("powerLow")->SetLineWidth(1);
  totSglAft->SetStats(kFALSE);
  totSglAft->SetLineColor(kRed);
  totSglAft->Draw("E1 same");

  doLegendSignal();

  canvasSignal->Print(outputPath+"totSigDist"+stName+".pdf");
}

void plotDiffDist::doLegendDiff() {
  legendDiff = new TLegend(0.6, 0.5, 0.98, 0.95);
  legendDiff->SetHeader("Station "+stName);
  legendDiff->AddEntry(diff12, "PMT_{1} - PMT_{2}","l");
  legendDiff->AddEntry(diff12, Form("Entries: %.f", diff12->GetEntries()), "");
  legendDiff->AddEntry(diff12, Form("Mean: %.2e #pm %.2e", 
        diff12->GetMean(), diff12->GetMeanError()), "");

  legendDiff->AddEntry(diff13, "PMT_{1} - PMT_{3}","l");
  legendDiff->AddEntry(diff13, Form("Entries: %.f", diff13->GetEntries()), "");
  legendDiff->AddEntry(diff13, Form("Mean: %.2e #pm %.2e",
        diff13->GetMean(), diff13->GetMeanError()), "");

  legendDiff->AddEntry(diff23, "PMT_{2} - PMT_{3}","l");
  legendDiff->AddEntry(diff23, Form("Entries: %.f", diff23->GetEntries()), "");
  legendDiff->AddEntry(diff23, Form("Mean: %.2e #pm %.2e",
        diff23->GetMean(), diff23->GetMeanError()), "");

  legendDiff->SetTextSize(0.04);
  legendDiff->SetBorderSize(0);
  legendDiff->SetFillStyle(0);
  legendDiff->Draw();
}

void plotDiffDist::doLegendSignal() {
  legendSignal = new TLegend(0.3, 0.15, 0.68, 0.6);
  legendSignal->SetHeader("Station "+stName);

  legendSignal->AddEntry(totSglBef, "Sep-Nov, Bef","l");
  legendSignal->AddEntry(totSglBef, Form("Entries: %.f", 
        totSglBef->GetEntries()), "");
  legendSignal->AddEntry(totSglBef, Form("Beta: %.2f #pm %.2f", 
        totSglBef->GetFunction("powerLow")->GetParameter(1), 
        totSglBef->GetFunction("powerLow")->GetParError(1)),"");
  legendSignal->AddEntry(totSglBef, Form("#chi^2/ndf: %.2f/%d",
        totSglBef->GetFunction("powerLow")->GetChisquare(), 
        totSglBef->GetFunction("powerLow")->GetNDF() ), "");

  legendSignal->AddEntry(totSglAft, "Feb-Apr, Aft","l");
  legendSignal->AddEntry(totSglAft, Form("Entries: %.f",
        totSglAft->GetEntries()), "");
  legendSignal->AddEntry(totSglAft, Form("Beta: %.2f #pm %.2f", 
        totSglAft->GetFunction("powerLow")->GetParameter(1),
        totSglAft->GetFunction("powerLow")->GetParError(1)), "");
  legendSignal->AddEntry(totSglAft, Form("#chi^2/ndf: %.2f/%d",
        totSglAft->GetFunction("powerLow")->GetChisquare(), 
        totSglAft->GetFunction("powerLow")->GetNDF() ), "");

  legendSignal->SetTextSize(0.04);
  legendSignal->SetBorderSize(0);
  legendSignal->SetFillStyle(0);
  legendSignal->Draw();
}

void plotDiffDist::doLegendTaus(TString pmtId, TH1D *tausBef, TH1D *tausAft) {
  legendTaus = new TLegend(0.1, 0.5, 0.5, 0.95);
  legendTaus->SetHeader("Station "+stName+", "+pmtId);
  legendTaus->AddEntry(tausBef, "Sep-Nov, Bef","l");
  legendTaus->AddEntry(tausBef, Form("Entries: %.f", tausBef->GetEntries()),"");
  legendTaus->AddEntry(tausBef, Form("Mean: %.2e #pm %.2e", 
        tausBef->GetMean(), tausBef->GetMeanError()), "");

  legendTaus->AddEntry(tausAft, "Feb-Apr, Aft","l");
  legendTaus->AddEntry(tausAft, Form("Entries: %.f", tausAft->GetEntries()),"");
  legendTaus->AddEntry(tausAft, Form("Mean: %.2e #pm %.2e",
        tausAft->GetMean(), tausAft->GetMeanError()), "");

  legendTaus->SetTextSize(0.04);
  legendTaus->SetBorderSize(0);
  legendTaus->SetFillStyle(0);
  legendTaus->Draw();
}

void plotDiffDist::writeRootFile() {
  canvasDiff->Write();
  legendDiff->Write();
  diff12->Write();
  diff13->Write();
  diff23->Write();

  canvasSignal->Write();
  legendSignal->Write();
  totSglBef->Write();
  totSglAft->Write();

  canvasTaus->Write();
  legendTaus->Write();
  muonTauBefPmt1->Write();
  muonTauBefPmt2->Write();
  muonTauBefPmt3->Write();

  muonTauAftPmt1->Write();
  muonTauAftPmt2->Write();
  muonTauAftPmt3->Write();

  outputFile->Write();
  outputFile->Close();
}
