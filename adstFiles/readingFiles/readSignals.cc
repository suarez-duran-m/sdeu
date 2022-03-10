#include <RecEvent.h>
#include <RecEventFile.h>
#include <DetectorGeometry.h>
#include <Traces.h>

#include <TFile.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 102, 76, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.13);
  canvas->SetTopMargin(0.014);
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderMode(0);
  return canvas;
}

using namespace std;

/********************/
/* Global variables */
/********************/

int main ( int argc, char** argv) {
  if ( argc < 2 ) { cout << endl
      << "Usage: " << argv[0] << " <file_adstfiles> <output> " << endl
      << "<file_adstfiles>: File with list of adst files to read" << endl
      << "<output>: name for adst_output file" << endl;
    exit(0);
  }
  cout << endl << endl;
  // Vectors to store signals values
  vector < double > totSigEvt;
  vector < double > timeSig;
  vector < double > st2read = {846, 860, 1190, 1191, 1220, 1779};
  vector < vector < double > > totSigSt( st2read.size() );
  vector < vector < double > > totSigStErr( st2read.size() );
  vector < vector < double > > accSigSt( st2read.size() );
  vector < vector < double > > timeSigSt( st2read.size() );
  // Charge Vs time, chStPmt[st][pmt][evt]
  vector < vector < vector < double > > > chStPmt( st2read.size() );
  vector < vector < vector < double > > > chStPmtErr( st2read.size() );
  vector < vector < vector < double > > > chStPmtTime( st2read.size() );
  for ( int i=0; i<st2read.size(); i++ ) {
    chStPmt[i].resize(3);
    chStPmtErr[i].resize(3);
    chStPmtTime[i].resize(3);
  }
  // Loaded starting
  const int loadStart = 1324857618;

  double sumSigEvt = 0.;
  double tmpSig = 0.;
  double tmpSigErr = 0.;
  int st_posVec = 0;
  bool readSt = false;
  // Reading ADST files
  for ( int adst_i = 1; adst_i<argc-1; adst_i++ ) {
    RecEventFile inputFile(argv[adst_i]);
    RecEvent *theRecEvent = new RecEvent();
    inputFile.SetBuffers(&theRecEvent);
    cout << "Opened " << argv[adst_i]
      << " with " << inputFile.GetNEvents()
      << " events " << endl;

    // Reading events    
    while ( inputFile.ReadNextEvent() == RecEventFile::eSuccess ) {
      const auto &sdEvent = theRecEvent->GetSDEvent();
      int nCand = sdEvent.GetNumberOfCandidates();
      sumSigEvt = 0.;
      // Reading station in the event
      for ( int st_i=0; st_i<nCand; st_i++ ) {
        readSt = false;
        // Looking if there is a desired station
        for ( int st2r_i=0; st2r_i<st2read.size(); st2r_i++ )
          if ( sdEvent.GetStation(st_i)->GetId() == st2read[st2r_i] ) {
            readSt = true;
            st_posVec = st2r_i;
          }
        if ( !readSt )
          continue;
        // Filling signals values
        tmpSig = sdEvent.GetStation(st_i)->GetTotalSignal();
        tmpSigErr = sdEvent.GetStation(st_i)->GetTotalSignalError();
        sumSigEvt += tmpSig;
        totSigSt[st_posVec].push_back( tmpSig );      
        totSigStErr[st_posVec].push_back( tmpSigErr );
        accSigSt[st_posVec].push_back( tmpSigErr/tmpSig );
        timeSigSt[st_posVec].push_back( sdEvent.GetGPSSecond() );

        // Filling charge values
        for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
          chStPmt[st_posVec][pmt_i].push_back(
              sdEvent.GetStation(st_i)->GetCharge(pmt_i+1) );
          chStPmtErr[st_posVec][pmt_i].push_back(
              sdEvent.GetStation(st_i)->GetChargeError(pmt_i+1) );
          chStPmtTime[st_posVec][pmt_i].push_back( sdEvent.GetGPSSecond() );
        }
      }
      if ( sumSigEvt > 0. ) {
        totSigEvt.push_back( sumSigEvt );
        timeSig.push_back( sdEvent.GetGPSSecond() );
      }
    }
  }

  // Plotting Signal distribution
  double frtbin = 0;
  double lstbin = 2e3;
  int nBins = (int)(lstbin - frtbin)/1.;
   
  TH1D *histTotSig = new TH1D("histTotSig", "", nBins, frtbin, lstbin);
  TH1D *histTotSigBef = new TH1D("histTotSigBef", "", nBins, frtbin, lstbin);
  TH1D *histTotSigAft = new TH1D("histTotSigAft", "", nBins, frtbin, lstbin);
  for ( int sig_i=0; sig_i<totSigEvt.size(); sig_i++ ) {
    if ( timeSig[sig_i] < loadStart )
      histTotSigBef->Fill( totSigEvt[sig_i] );
    else
      histTotSigAft->Fill( totSigEvt[sig_i] );
    histTotSig->Fill( totSigEvt[sig_i] );
  }

  TCanvas *sigDistCanv = canvasStyle("sigDistCanv");
  sigDistCanv->cd();
  sigDistCanv->SetLogy(1);
  sigDistCanv->SetLogx(1);
  
  histTotSig->SetStats(kFALSE);
  histTotSig->GetYaxis()->SetTitle("Counts [au]");
  histTotSig->GetXaxis()->SetTitle("S [VEM]");
  histTotSig->GetXaxis()->SetTitleOffset(1.3);
  histTotSig->SetLineColor(kBlack);
  histTotSig->Draw();

  histTotSigBef->SetLineColor(kBlue);
  histTotSigBef->Draw("same");
  histTotSigAft->SetLineColor(kRed);
  histTotSigAft->Draw("same");

  TLegend *lgnd = new TLegend(0.65, 0.45, 0.98, 0.98);
  lgnd->AddEntry(histTotSig, "St.: 846, 860, 1190", "");
  lgnd->AddEntry(histTotSig, "     1191, 1220, 1779", "");
  lgnd->AddEntry(histTotSig, "Nov to Feb", "l");
  lgnd->AddEntry(histTotSig, Form("Entries %.f",histTotSig->GetEntries()), "");  
  lgnd->AddEntry(histTotSigBef, "Before Loaded", "l");
  lgnd->AddEntry(histTotSigBef, 
      Form("Entries %.f",histTotSigBef->GetEntries()), "");
  lgnd->AddEntry(histTotSigAft, "After Loaded", "l");
  lgnd->AddEntry(histTotSigAft, 
      Form("Entries %.f",histTotSigAft->GetEntries()), "");
  lgnd->SetTextSize(0.04);
  lgnd->SetBorderSize(0);
  lgnd->SetFillStyle(0);
  lgnd->Draw();
  sigDistCanv->Print("totSigDist.pdf");

  TGraphErrors *grpSigSt = new TGraphErrors(timeSigSt[0].size(), 
      &timeSigSt[0].front(), &totSigSt[0].front(), 0, &totSigStErr[0].front() );

  TCanvas *sigStCanv = canvasStyle("sigStCanv");
  sigStCanv->cd();
  sigStCanv->SetLogy(1);
 
  grpSigSt->SetTitle("St. 846");
  grpSigSt->GetXaxis()->SetTimeFormat("%m/%d");
  grpSigSt->GetXaxis()->SetTimeOffset(315964782,"gmt");
  grpSigSt->GetXaxis()->SetTitle("Time [month/day]");
  grpSigSt->GetYaxis()->SetTitle("S [VEM]");
  grpSigSt->SetMarkerStyle(20);
  grpSigSt->SetMarkerSize(1.5);
  grpSigSt->Draw("AP");

  TLine *line = new TLine(loadStart, 0, loadStart, 1e3);
  line->SetLineWidth(1);
  line->SetLineColor(kRed);
  line->Draw();

  sigStCanv->Print("signalStation.pdf");

  TGraph *grpAccSigSt = new TGraph(timeSigSt[0].size(),
      &timeSigSt[0].front(), &accSigSt[0].front() );

  frtbin = 0.;
  lstbin = 10.;
  nBins = (lstbin - frtbin) / 0.01;
  TH1D *accSgDstStBef = new TH1D("accSgDstStBef","", nBins, frtbin, lstbin);
  TH1D *accSgDstStAft = new TH1D("accSgDstStAft","", nBins, frtbin, lstbin);
  for ( int acc_i=0; acc_i<accSigSt[0].size(); acc_i++ ) {
    if ( timeSigSt[0][acc_i] < loadStart )
      accSgDstStBef->Fill( accSigSt[0][acc_i] );
    else
      accSgDstStAft->Fill( accSigSt[0][acc_i] );
  }

  TCanvas *accSigStCanv = canvasStyle("accSigStCanv");
  accSigStCanv->cd();
  accSigStCanv->SetLogy(1);
  grpAccSigSt->SetTitle("St. 846");
  grpAccSigSt->GetXaxis()->SetTimeFormat("%m/%d");
  grpAccSigSt->GetXaxis()->SetTimeOffset(315964782,"gmt");
  grpAccSigSt->GetXaxis()->SetTitle("Time [month/day]");
  grpAccSigSt->GetYaxis()->SetTitle("ErrSignal/Signal [au]");
  grpAccSigSt->SetMarkerColor(kBlue);
  grpAccSigSt->SetMarkerStyle(20);
  grpAccSigSt->SetMarkerSize(1.5);
  grpAccSigSt->Draw("AP");
  
  line = new TLine(loadStart, 0, loadStart, 1.58e2);
  line->SetLineColor(kRed);
  line->SetLineWidth(1);
  line->Draw();

  accSigStCanv->Print("accSigSt.pdf");

  TCanvas *accSigDistStCanv = canvasStyle("accDistSigStCanv");
  accSigDistStCanv->cd();

  accSgDstStBef->SetStats(kFALSE);
  accSgDstStBef->GetXaxis()->SetTitle("ErrSignal/Signal [au]");
  accSgDstStBef->GetXaxis()->SetRangeUser(0, 3);
  accSgDstStBef->GetYaxis()->SetTitle("Counts [au]");
  accSgDstStBef->SetLineColor(kBlue);
  accSgDstStBef->Draw();
  accSgDstStAft->SetLineColor(kRed);
  accSgDstStAft->Draw("same");

  lgnd = new TLegend(0.6, 0.45, 0.98, 0.98);
  lgnd->AddEntry(accSgDstStBef, "St.: 846", "");
  lgnd->AddEntry(accSgDstStBef, "Before Loaded", "l");
  lgnd->AddEntry(accSgDstStBef,
      Form("Entries %.f",accSgDstStBef->GetEntries()), "");
  lgnd->AddEntry(accSgDstStBef, 
      Form("Mean %.2e #pm %.2e", accSgDstStBef->GetMean(), 
        accSgDstStBef->GetMeanError()), "");
  lgnd->AddEntry(accSgDstStAft, "After Loaded", "l");
  lgnd->AddEntry(accSgDstStAft, 
      Form("Entries %.f",accSgDstStAft->GetEntries()), "");
  lgnd->AddEntry(accSgDstStAft,
      Form("Mean %.2e #pm %.2e", accSgDstStAft->GetMean(), 
        accSgDstStAft->GetMeanError()), "");
  lgnd->SetTextSize(0.04);
  lgnd->SetBorderSize(0);
  lgnd->SetFillStyle(0);
  lgnd->Draw();

  accSigDistStCanv->Print("accDistSigSt.pdf");

  TGraphErrors *grpChStPmt1 = new TGraphErrors(chStPmt[0][0].size(), 
      &chStPmtTime[0][0].front(), &chStPmt[0][0].front(), 0, 
      &chStPmtErr[0][0].front());
  TCanvas *chStCanv = canvasStyle("chStCanv");
  chStCanv->cd();

  grpChStPmt1->SetTitle("St. 846 PMT 1");
  grpChStPmt1->GetXaxis()->SetTimeFormat("%m/%d");
  grpChStPmt1->GetXaxis()->SetTimeOffset(315964782,"gmt");
  grpChStPmt1->GetXaxis()->SetTitle("Time [month/day]");
  grpChStPmt1->GetYaxis()->SetTitle("VEM [ADC]");
  grpChStPmt1->SetMarkerStyle(20);
  grpChStPmt1->SetMarkerSize(1.5);
  grpChStPmt1->Draw("AP");

  line = new TLine(loadStart, 1.38e3, loadStart, 1.58e3);
  line->SetLineColor(kRed);
  line->SetLineWidth(1);
  line->Draw();
  chStCanv->Print("chargeStation.pdf");

  return 1;
}
