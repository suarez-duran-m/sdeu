#ifndef PLOTDIFFDIST_H
#define PLOTDIFFDIST_H

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

class plotDiffDist {
  public:
    plotDiffDist(TString stId, TString printPath, TH1D *pmt12, TH1D *pmt13, 
        TH1D *pmt23, TH1D *totSignalBef, TH1D *totSignalAft);
    void plotTausMuon(TH1D *tauBefPmt1, TH1D *tauBefPmt2, TH1D *tauBefPmt3,
        TH1D *tauAftPmt1, TH1D *tauAftPmt2, TH1D *tauAftPmt3);
    void writeRootFile();

  private:
    TFile *outputFile;
    TH1D *diff12;
    TH1D *diff13;
    TH1D *diff23;
    TH1D *muonTauBefPmt1;
    TH1D *muonTauBefPmt2;
    TH1D *muonTauBefPmt3;
    TH1D *muonTauAftPmt1;
    TH1D *muonTauAftPmt2;
    TH1D *muonTauAftPmt3;

    TH1D *totSglBef;
    TH1D *totSglAft;

    TString stName;
    TString outputPath;
    TCanvas *canvasDiff;
    TLegend *legendDiff;
    TCanvas *canvasSignal;
    TLegend *legendSignal;
    TCanvas *canvasTaus;
    TLegend *legendTaus;

    TCanvas *doCanvas(TString canvasName);
    void doDiffDistPlot();
    void doLegendDiff();
    void doDisTotSglPlot();
    void doLegendSignal();
    void doDisTausPlot(TString pmtId, TH1D *tausBef, TH1D *tausAft);
    void doLegendTaus(TString pmtId, TH1D *tausBef, TH1D *tausAft);
};
#endif
