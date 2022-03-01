#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <IoAuger.h>
#include <Ec.h>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TText.h>

#include "fitcharge.h"

using namespace std;

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

TPaveStats *getPaveStats(TH1D *histo) {
  TPaveStats *ptstats = new TPaveStats(0.55,0.65,0.95,0.95,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);
  ptstats->SetTextSize(0.035);
  ptstats->SetStatFormat(".2e");
  TText *ptstats_LaTex = ptstats->AddText(Form("%s", histo->GetName()));
  ptstats_LaTex->SetTextSize(0.06);
  ptstats->SetOptFit(0);
  return ptstats;
}

void plotZvals(TString strHigh, unsigned int st, vector<vector<double>> zVals, 
    vector<vector<double>> zTime) {
  TGraph *grpZpmt1 = new TGraph( zVals[0].size(), &zTime[0][0], &zVals[0][0] );
  TGraph *grpZpmt2 = new TGraph( zVals[1].size(), &zTime[1][0], &zVals[1][0] );
  TGraph *grpZpmt3 = new TGraph( zVals[2].size(), &zTime[2][0], &zVals[2][0] );
  TH1D *distZpmt1 = new TH1D(Form("St%u PMT1 "+strHigh+" gain",st), "", 100, 0., 10.);
  TH1D *distZpmt2 = new TH1D(Form("St%u PMT2"+strHigh+" gain",st), "", 100, 0., 10.);
  TH1D *distZpmt3 = new TH1D(Form("St%u PMT3"+strHigh+" gain",st), "", 100, 0., 10.);
  for ( auto & i : zVals[0] )
    distZpmt1->Fill( i );
  for ( auto & i : zVals[1] ) 
    distZpmt2->Fill( i );
  for ( auto & i : zVals[2] )
    distZpmt3->Fill( i );

  TCanvas *Zcanvas = canvasStyle("Zcanvas");
  Zcanvas->cd();  
  Zcanvas->Draw();
  TPad *p1 = new TPad("p1","p1",0.01,0.5,0.99,1.);
  p1->Draw();
  p1->Divide(3,1);
  TPad *p11 = (TPad*)p1->cd(1);
  p11->Draw();
  grpZpmt1->SetTitle("");
  grpZpmt1->GetYaxis()->SetTitle("Z [au]");
  grpZpmt1->GetYaxis()->SetRangeUser(0., 5.);
  grpZpmt1->GetXaxis()->SetTimeFormat("%m/%d");
  grpZpmt1->GetXaxis()->SetTimeOffset(0,"gmt");
  grpZpmt1->Draw("AP");
  TPad *p12 = (TPad*)p1->cd(2);
  p12->Draw();
  grpZpmt2->SetTitle("");
  grpZpmt2->GetYaxis()->SetTitle("Z [au]");
  grpZpmt2->GetYaxis()->SetRangeUser(0., 5.);
  grpZpmt2->GetXaxis()->SetTimeFormat("%m/%d");
  grpZpmt2->GetXaxis()->SetTimeOffset(0,"gmt");
  grpZpmt2->Draw("AP");
  TPad *p13 = (TPad*)p1->cd(3);
  p13->Draw();
  grpZpmt3->SetTitle("");
  grpZpmt3->GetYaxis()->SetTitle("Z [au]");
  grpZpmt3->GetYaxis()->SetRangeUser(0., 5.);
  grpZpmt3->GetXaxis()->SetTimeFormat("%m/%d");
  grpZpmt3->GetXaxis()->SetTimeOffset(0,"gmt");
  grpZpmt3->Draw("AP");

  Zcanvas->cd(0);
  gStyle->SetOptStat("neMR");
  TPad *p2 = new TPad("p2","p2",0.01,0.,0.99,0.5);
  p2->Draw();
  p2->Divide(3,1);
  TPad *p21 = (TPad*)p2->cd(1);
  p21->Draw();
  p21->SetLogy(1);
  TPaveStats *ptstats1 = getPaveStats(distZpmt1);
  ptstats1->Draw();
  distZpmt1->GetListOfFunctions()->Add(ptstats1);
  ptstats1->SetParent(distZpmt1);
  distZpmt1->SetTitle("");
  distZpmt1->GetYaxis()->SetTitle("Counts [au]");
  distZpmt1->GetXaxis()->SetTitle("Z [au]");
  distZpmt1->Draw();

  TPad *p22 = (TPad*)p2->cd(2);
  p22->Draw();
  p22->SetLogy(1);
  TPaveStats *ptstats2 = getPaveStats(distZpmt2);
  ptstats2->Draw();
  distZpmt2->GetListOfFunctions()->Add(ptstats2);
  ptstats2->SetParent(distZpmt2);
  distZpmt2->SetTitle(""); 
  distZpmt2->GetYaxis()->SetTitle("Counts [au]");
  distZpmt2->GetXaxis()->SetTitle("Z [au]");
  distZpmt2->Draw();

  TPad *p23 = (TPad*)p2->cd(3);
  p23->Draw();
  p23->SetLogy(1);
  TPaveStats *ptstats3 = getPaveStats(distZpmt3);
  ptstats3->Draw();
  distZpmt3->GetListOfFunctions()->Add(ptstats3);
  ptstats3->SetParent(distZpmt3);
  distZpmt3->SetTitle(""); 
  distZpmt3->GetYaxis()->SetTitle("Counts [au]");
  distZpmt3->GetXaxis()->SetTitle("Z [au]");
  distZpmt3->Draw();

  Zcanvas->Print(Form("../plots2/blZst%u"+strHigh+".pdf",st));
  Zcanvas->Close();
  grpZpmt1->~TGraph();
  grpZpmt2->~TGraph();
  grpZpmt3->~TGraph();
  distZpmt1->~TH1D();
  distZpmt2->~TH1D();
  distZpmt3->~TH1D();
}

void plotQpk(unsigned int st, vector<vector<double>> qpkVals, 
    vector<vector<double>> qpkErr, vector<vector<double>> qpkTime) {
  TGraphErrors *grpQpkpmt1 = new TGraphErrors(qpkVals[0].size(), 
      &qpkTime[0][0], &qpkVals[0][0], 0, &qpkErr[0][0]);
  TGraphErrors *grpQpkpmt2 = new TGraphErrors(qpkVals[1].size(), 
      &qpkTime[1][0], &qpkVals[1][0], 0, &qpkErr[1][0]);
  TGraphErrors *grpQpkpmt3 = new TGraphErrors(qpkVals[2].size(), 
      &qpkTime[2][0], &qpkVals[2][0], 0, &qpkErr[2][0]);
  TH1D *distQpkpmt1 = new TH1D(Form("St%u PMT1",st), "", 3e3, 0., 3e3);
  TH1D *distQpkpmt2 = new TH1D(Form("St%u PMT2",st), "", 3e3, 0., 3e3);
  TH1D *distQpkpmt3 = new TH1D(Form("St%u PMT3",st), "", 3e3, 0., 3e3);
  for ( auto & i : qpkVals[0] )
    distQpkpmt1->Fill( i );
  for ( auto & i : qpkVals[1] ) 
    distQpkpmt2->Fill( i );
  for ( auto & i : qpkVals[2] )
    distQpkpmt3->Fill( i );

  TCanvas *qpkCanvas = canvasStyle("qpkCanvas");
  qpkCanvas->cd();  
  qpkCanvas->Draw();
  TPad *p1 = new TPad("p1","p1",0.005,0.5,1.,1.);
  p1->Draw();
  p1->Divide(3,1);
  TPad *p11 = (TPad*)p1->cd(1);
  p11->Draw();
  grpQpkpmt1->SetTitle("");
  grpQpkpmt1->GetYaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  grpQpkpmt1->GetYaxis()->SetTitleOffset(1.6);
  //grpQpkpmt1->GetYaxis()->SetRangeUser(0., 5.);
  grpQpkpmt1->GetXaxis()->SetTimeFormat("%m/%d");
  grpQpkpmt1->GetXaxis()->SetTimeOffset(0,"gmt");
  grpQpkpmt1->Draw("AP");
  TPad *p12 = (TPad*)p1->cd(2);
  p12->SetLeftMargin(2.);
  p12->Draw();
  grpQpkpmt2->SetTitle("");
  grpQpkpmt2->GetYaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  grpQpkpmt2->GetYaxis()->SetTitleOffset(1.6);
  //grpQpkpmt2->GetYaxis()->SetRangeUser(0., 5.);
  grpQpkpmt2->GetXaxis()->SetTimeFormat("%m/%d");
  grpQpkpmt2->GetXaxis()->SetTimeOffset(0,"gmt");
  grpQpkpmt2->Draw("AP");
  TPad *p13 = (TPad*)p1->cd(3);
  p13->SetLeftMargin(-1.);
  p13->Draw();
  grpQpkpmt3->SetTitle("");
  grpQpkpmt3->GetYaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  grpQpkpmt3->GetYaxis()->SetTitleOffset(1.6);
  //grpQpkpmt3->GetYaxis()->SetRangeUser(0., 5.);
  grpQpkpmt3->GetXaxis()->SetTimeFormat("%m/%d"); grpQpkpmt3->GetXaxis()->SetTimeOffset(0,"gmt");
  grpQpkpmt3->Draw("AP");

  qpkCanvas->cd(0);
  gStyle->SetOptStat("neMR");
  TPad *p2 = new TPad("p2","p2",0.005,0.,1.,0.5);
  p2->Draw();
  p2->Divide(3,1);
  TPad *p21 = (TPad*)p2->cd(1);
  p21->Draw();
  p21->SetLogy(1);
  TPaveStats *ptstats1 = getPaveStats(distQpkpmt1);
  ptstats1->Draw();
  distQpkpmt1->GetListOfFunctions()->Add(ptstats1);
  ptstats1->SetParent(distQpkpmt1);
  distQpkpmt1->SetTitle("");
  distQpkpmt1->GetYaxis()->SetTitle("Counts [au]");
  distQpkpmt1->GetXaxis()->SetTitleOffset(1.2);
  distQpkpmt1->GetXaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  distQpkpmt1->GetXaxis()->SetRangeUser(1e3, 2.5e3);
  distQpkpmt1->Draw();

  TPad *p22 = (TPad*)p2->cd(2);
  p22->Draw();
  p22->SetLogy(1);
  TPaveStats *ptstats2 = getPaveStats(distQpkpmt2);
  ptstats2->Draw();
  distQpkpmt2->GetListOfFunctions()->Add(ptstats2);
  ptstats2->SetParent(distQpkpmt2);
  distQpkpmt2->SetTitle(""); 
  distQpkpmt2->GetYaxis()->SetTitle("Counts [au]");
  distQpkpmt2->GetXaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  distQpkpmt2->GetXaxis()->SetRangeUser(1e3, 2.5e3);
  distQpkpmt2->Draw();

  TPad *p23 = (TPad*)p2->cd(3);
  p23->Draw();
  p23->SetLogy(1);
  TPaveStats *ptstats3 = getPaveStats(distQpkpmt3);
  ptstats3->Draw();
  distQpkpmt3->GetListOfFunctions()->Add(ptstats3);
  ptstats3->SetParent(distQpkpmt3);
  distQpkpmt3->SetTitle(""); 
  distQpkpmt3->GetYaxis()->SetTitle("Counts [au]");
  distQpkpmt3->GetXaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  distQpkpmt3->GetXaxis()->SetRangeUser(1e3, 2.5e3);
  distQpkpmt3->Draw();

  qpkCanvas->Print(Form("../plots2/qpkSt%u.pdf",st));
  qpkCanvas->Close();
}

void plotSgnl(TString strHigh, unsigned int st, vector<vector<double>> sgnlVals, 
    vector<vector<double>> sgnlTime) {
  TGraph *grpSgnlpmt1 = new TGraph( sgnlVals[0].size(), &sgnlTime[0][0], 
      &sgnlVals[0][0] );
  TGraph *grpSgnlpmt2 = new TGraph( sgnlVals[1].size(), &sgnlTime[1][0], 
      &sgnlVals[1][0] );
  TGraph *grpSgnlpmt3 = new TGraph( sgnlVals[2].size(), &sgnlTime[2][0], 
      &sgnlVals[2][0] );
  TH1D *distSgnlpmt1 = new TH1D(Form("St%u PMT1 "+strHigh+" gain",st), 
      "", 1e3, 0., 1e2);
  TH1D *distSgnlpmt2 = new TH1D(Form("St%u PMT2 "+strHigh+" gain",st), 
      "", 1e3, 0., 1e2);
  TH1D *distSgnlpmt3 = new TH1D(Form("St%u PMT3 "+strHigh+" gain",st), 
      "", 1e3, 0., 1e2);
  for ( auto & i : sgnlVals[0] )
    distSgnlpmt1->Fill( i );
  for ( auto & i : sgnlVals[1] ) 
    distSgnlpmt2->Fill( i );
  for ( auto & i : sgnlVals[2] )
    distSgnlpmt3->Fill( i );

  TCanvas *sgnlCanvas = canvasStyle("sgnlCanvas");
  sgnlCanvas->cd();  
  sgnlCanvas->Draw();
  TPad *p1 = new TPad("p1","p1",0.01,0.5,0.99,1.);
  p1->Draw();
  p1->Divide(3,1);
  TPad *p11 = (TPad*)p1->cd(1);
  p11->Draw();
  grpSgnlpmt1->SetTitle("");
  grpSgnlpmt1->GetYaxis()->SetTitle("S [Q^{pk}_{VEM}]");
  //grpSgnlpmt1->GetYaxis()->SetRangeUser(0., 5.);
  grpSgnlpmt1->GetXaxis()->SetTimeFormat("%m/%d");
  grpSgnlpmt1->GetXaxis()->SetTimeOffset(0,"gmt");
  grpSgnlpmt1->Draw("AP");
  TPad *p12 = (TPad*)p1->cd(2);
  p12->Draw();
  grpSgnlpmt2->SetTitle("");
  grpSgnlpmt2->GetYaxis()->SetTitle("S [Q^{pk}_{VEM}]");
  //grpSgnlpmt2->GetYaxis()->SetRangeUser(0., 5.);
  grpSgnlpmt2->GetXaxis()->SetTimeFormat("%m/%d");
  grpSgnlpmt2->GetXaxis()->SetTimeOffset(0,"gmt");
  grpSgnlpmt2->Draw("AP");
  TPad *p13 = (TPad*)p1->cd(3);
  p13->Draw();
  grpSgnlpmt3->SetTitle("");
  grpSgnlpmt3->GetYaxis()->SetTitle("S [Q^{pk}_{VEM}]");
  //grpSgnlpmt3->GetYaxis()->SetRangeUser(0., 5.);
  grpSgnlpmt3->GetXaxis()->SetTimeFormat("%m/%d");
  grpSgnlpmt3->GetXaxis()->SetTimeOffset(0,"gmt");
  grpSgnlpmt3->Draw("AP");

  sgnlCanvas->cd(0);
  gStyle->SetOptStat("neMR");
  TPad *p2 = new TPad("p2","p2",0.01,0.,0.99,0.5);
  p2->Draw();
  p2->Divide(3,1);
  TPad *p21 = (TPad*)p2->cd(1);
  p21->Draw();
  p21->SetLogy(1);
  p21->SetLogx(1);
  TPaveStats *ptstats1 = getPaveStats(distSgnlpmt1);
  ptstats1->Draw();
  distSgnlpmt1->GetListOfFunctions()->Add(ptstats1);
  ptstats1->SetParent(distSgnlpmt1);
  distSgnlpmt1->SetTitle("");
  distSgnlpmt1->GetYaxis()->SetTitle("Counts [au]");
  distSgnlpmt1->GetXaxis()->SetTitle("S [Q^{pk}_{VEM}]");
  distSgnlpmt1->Draw();

  TPad *p22 = (TPad*)p2->cd(2);
  p22->Draw();
  p22->SetLogy(1);
  p22->SetLogx(1);
  TPaveStats *ptstats2 = getPaveStats(distSgnlpmt2);
  ptstats2->Draw();
  distSgnlpmt2->GetListOfFunctions()->Add(ptstats2);
  ptstats2->SetParent(distSgnlpmt2);
  distSgnlpmt2->SetTitle(""); 
  distSgnlpmt2->GetYaxis()->SetTitle("Counts [au]");
  distSgnlpmt2->GetXaxis()->SetTitle("S [Q^{pk}_{VEM}]");
  distSgnlpmt2->Draw();

  TPad *p23 = (TPad*)p2->cd(3);
  p23->Draw();
  p23->SetLogy(1);
  p23->SetLogx(1);
  TPaveStats *ptstats3 = getPaveStats(distSgnlpmt3);
  ptstats3->Draw();
  distSgnlpmt3->GetListOfFunctions()->Add(ptstats3);
  ptstats3->SetParent(distSgnlpmt3);
  distSgnlpmt3->SetTitle(""); 
  distSgnlpmt3->GetYaxis()->SetTitle("Counts [au]");
  distSgnlpmt3->GetXaxis()->SetTitle("S [Q^{pk}_{VEM}]");
  distSgnlpmt3->Draw();

  sgnlCanvas->Print(Form("../plots2/sgnlSt%u"+strHigh+".pdf",st));
  sgnlCanvas->Close();
}

double getmean( vector<double> arr ){
  double mean = 0.;
  for (auto & i : arr)
    mean += i;
  return mean/arr.size();
}

float getrms( int arr[], float meanarr ){
  float rms = 0.;
  for (int i=0; i<100; i++)
    rms += (arr[i] - meanarr)*(arr[i] - meanarr);
  return rms/100.;
}

int main (int argc, char *argv[]) {
   if ( argc < 4 ) {
     cout << endl
         << "Usage: " << argv[0] << " <stationsFile>  <output>  <files>" << endl
         << "  <stationsFile>: file with a list of stations" << endl
         << "  <output>: output file with whatever you want inside" << endl
         << "  <files>: IoSd or IoAuger files to be read" << endl;
    exit(0);
  }
  // Open files with raw data
  const char* stationsFileName = argv[1];
  AugerIoSd input(argc-3, argv+3);
  const unsigned int totalNrEvents = input.NumberOfEvents();  
  ifstream stationsFile(stationsFileName, ios::in);
  if (!stationsFile.is_open()){
    cout << "Could not open file: " << stationsFileName << endl;
    exit(0);
  }
  // Reading stations to study from file
  vector<unsigned int> stationsIds;
  while (stationsFile.good()) {
    unsigned int st = 0;
    stationsFile >> st;
    if (st)
      stationsIds.push_back(st);
  }  
  if (stationsIds.empty()){
    cout << "Please specify the stations ids in the file " << endl;
    exit(0);
  }
  // Auger variable to read data 
  TString nameStati = to_string( stationsIds[0] );
  // Creating root file to save results
  TFile hfile("bl"+nameStati+".root","RECREATE","");

  // Creating vectors to calculate store Z values and base line
  // [st][pmt][evt] = Z
  double tmpZ = 0.;
  vector < vector < vector < double > > > ZlowVals(stationsIds.size());
  vector < vector < vector < double > > > ZhigVals(stationsIds.size());  
  vector < vector < vector < double > > > Ztime(stationsIds.size());  
  int smplBinsBl = 100;
  double tmpFrstErr = 0.;
  double tmpLastErr = 0.;
  // TH1D for first and last bin BL distributions
  TH1D *frstBinsBlLow;
  TH1D *frstBinsBlHig;
  TH1D *lastBinsBlLow;
  TH1D *lastBinsBlHig;
  int stPosVect = 0;
  int tmpPosVect = 0; 
  // For Qpk calculation
  fitcharge fittingQpk;
  TH1F *receivedChHisto = new TH1F();
  TH1F *receivedChCrr = new TH1F(); 
  vector < vector < vector < double > > > qpk(stationsIds.size());
  vector < vector < vector < double > > > qpkErr(stationsIds.size());
  vector < vector < vector < double > > > qpkTime(stationsIds.size());
  double fadc2qpk = 0.;
  // For Signal
  vector < vector < vector < double > > > signalLow(stationsIds.size());
  vector < vector < vector < double > > > signalHig(stationsIds.size());
  vector < vector < vector < double > > > signalTime(stationsIds.size());
  double tmpSigLow = 0.;
  double tmpSigHig = 0.;

  unsigned int previusEvent = 0;
  int sameUtc = 0;
  int nbinsInBl = 0;

  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;
  EventPos pos;
  // Reading Events
  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
    ++nrEventsRead;
    if (nrEventsRead%1000 == 0) {
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }
    bool found = false;
    IoSdEvent event(pos);
    // Searching is the station triggered in this event
    for (unsigned int evt_i = 0 ; evt_i < event.Stations.size(); ++evt_i) {
      found = false;
      tmpPosVect = 0;
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
             iter!= stationsIds.end(); ++iter) {
        if (event.Stations[evt_i].Id == *iter) {
          found = true;
          stPosVect = tmpPosVect;
        }
        tmpPosVect++;
      }
      if (!found)
        continue;
      // Re-size for PMTs
      ZlowVals[stPosVect].resize(3);
      ZhigVals[stPosVect].resize(3);
      Ztime[stPosVect].resize(3);
      qpk[stPosVect].resize(3);
      qpkErr[stPosVect].resize(3);
      qpkTime[stPosVect].resize(3);
      signalLow[stPosVect].resize(3);
      signalHig[stPosVect].resize(3);
      signalTime[stPosVect].resize(3);
      // Asking if it is UUB
      if ( event.Stations[evt_i].IsUUB && event.Id != previusEvent 
          && sameUtc != event.utctime() ) {
        previusEvent = event.Id;
        sameUtc = event.utctime();
        cout << "# Event " << event.Id << " Station " << event.Stations[evt_i].Id 
          << endl;
        IoSdEvent event(pos);
        // Filter for error 
        if ( !(event.Stations[evt_i].Error==256) )
          continue;
        nbinsInBl = event.Stations[evt_i].UFadc->NSample;
        // Reading baseline bins
        for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
          // TH1D for first and last bins of BL
          frstBinsBlLow = new TH1D(Form("frstBinsBlLow%d%d%ld",
                stationsIds[stPosVect], pmt_i, event.utctime()),"", 500, 0, 500);
          frstBinsBlHig = new TH1D(Form("frstBinsBlHig%d%d%ld",
                stationsIds[stPosVect], pmt_i, event.utctime()),"", 500, 0, 500);
          lastBinsBlLow = new TH1D(Form("lastBinsBlLow%d%d%ld",
                stationsIds[stPosVect], pmt_i, event.utctime()),"", 500, 0, 500);
          lastBinsBlHig = new TH1D(Form("lastBinsBlHig%d%d%ld",
                stationsIds[stPosVect], pmt_i, event.utctime()),"", 500, 0, 500);
          for (int k=0; k<smplBinsBl; k++) {
            frstBinsBlLow->Fill( event.Stations[evt_i].UFadc->GetValue(pmt_i,1,k) );
            frstBinsBlHig->Fill( event.Stations[evt_i].UFadc->GetValue(pmt_i,0,k) );
            lastBinsBlLow->Fill( 
                event.Stations[evt_i].UFadc->GetValue(pmt_i,1,nbinsInBl-k) );
            lastBinsBlHig->Fill( 
                event.Stations[evt_i].UFadc->GetValue(pmt_i,0,nbinsInBl-k) );
          }
          // Calculating and storing Z values
          tmpFrstErr = frstBinsBlLow->GetMeanError();
          tmpLastErr = lastBinsBlLow->GetMeanError();
          tmpZ = (frstBinsBlLow->GetMean() - lastBinsBlLow->GetMean()) 
            / (sqrt( tmpFrstErr*tmpFrstErr + tmpLastErr*tmpLastErr ));
          ZlowVals[stPosVect][pmt_i].push_back( tmpZ );
          tmpFrstErr = frstBinsBlHig->GetMeanError();
          tmpLastErr = lastBinsBlHig->GetMeanError();
          tmpZ = (frstBinsBlHig->GetMean() - lastBinsBlHig->GetMean()) 
            / (sqrt( tmpFrstErr*tmpFrstErr + tmpLastErr*tmpLastErr ));
          ZhigVals[stPosVect][pmt_i].push_back( tmpZ );
          Ztime[stPosVect][pmt_i].push_back( event.utctime() );
          // Cut by Z values for high gain
          if ( tmpZ > 2.0 )
            continue;
          receivedChHisto = event.Stations[evt_i].HCharge(pmt_i);
          fittingQpk.setChCrr(*receivedChHisto, 
              event.Stations[evt_i].Histo->Offset[pmt_i+6], 
              Form("%d-%d-%ld", stPosVect, pmt_i, event.utctime()));
          receivedChCrr = fittingQpk.getChCrr();
          fittingQpk.getFitCh(*receivedChCrr, 30 );
          if ( fittingQpk.vemPosCh > 0. ) {
            qpk[stPosVect][pmt_i].push_back( fittingQpk.vemPosCh );
            qpkErr[stPosVect][pmt_i].push_back( fittingQpk.vemPosChErr );
            qpkTime[stPosVect][pmt_i].push_back( event.utctime() );
            fadc2qpk = 1./fittingQpk.vemPosCh;
            tmpSigLow = 0.;
            tmpSigHig = 0.;
            for ( int k=0; k<nbinsInBl; k++ ) {
              tmpSigLow += event.Stations[evt_i].UFadc->GetValue(pmt_i,1,k) 
                - frstBinsBlLow->GetMean();
              tmpSigHig += event.Stations[evt_i].UFadc->GetValue(pmt_i,0,k) 
                - frstBinsBlHig->GetMean();
            }
            signalLow[stPosVect][pmt_i].push_back( tmpSigLow*fadc2qpk );
            signalHig[stPosVect][pmt_i].push_back( tmpSigHig*fadc2qpk );
            signalTime[stPosVect][pmt_i].push_back( event.utctime() );

          }
          frstBinsBlLow->Clear();
          frstBinsBlLow->Delete();
          frstBinsBlHig->Clear();
          frstBinsBlHig->Delete();
          lastBinsBlLow->Clear();
          lastBinsBlLow->Delete();
          lastBinsBlHig->Clear();
          lastBinsBlHig->Delete();
        }
      }
    }
  }
  
  for ( unsigned int st_i=0; st_i<stationsIds.size(); st_i++ ) {
    plotZvals("Low", stationsIds[st_i], ZlowVals[st_i], Ztime[st_i]);
    plotZvals("High", stationsIds[st_i], ZhigVals[st_i], Ztime[st_i]);  
 
    plotQpk(stationsIds[st_i], qpk[st_i], qpkErr[st_i], qpkTime[st_i]);

    plotSgnl("High", stationsIds[st_i], signalHig[st_i], signalTime[st_i]);
  }



  
  hfile.Write();
  hfile.Close();

  return 0;
}
