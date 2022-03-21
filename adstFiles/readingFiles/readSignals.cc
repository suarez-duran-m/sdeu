#include <RecEvent.h>
#include <RecEventFile.h>
#include <DetectorGeometry.h>
#include <Traces.h>
#include <TraceType.h>

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

using namespace std;

TCanvas *canvasStyle(TString name);
void plotSigDist(int st_slct, double frtbin, double lstbin, int nbins, 
    vector<vector<double>> sigValues, vector<vector<double>> time);

/********************/
/* Global variables */
/********************/

// Loaded sPMT starting
const int loadStart = 1324857618;
vector < double > st2read = {846, 860, 1190, 1191, 1220, 1779};

int main ( int argc, char** argv) {
  if ( argc < 2 ) { cout << endl
      << "Usage: " << argv[0] << " <file_adstfiles> <output> " << endl
      << "<file_adstfiles>: File with list of adst files to read" << endl
      << "<output>: name for adst_output file" << endl;
    exit(0);
  }
  cout << endl << endl;
  // Vectors to store signals values
  vector < vector < double > > totSigEvt( st2read.size() );
  vector < vector < double > > timeSig( st2read.size() );
  vector < vector < double > > sigStTime( st2read.size() );
  vector < vector < double > > sigStTimeErr( st2read.size() );
  vector < vector < double > > accSigSt( st2read.size() );
  vector < vector < double > > timeSigSt( st2read.size() );
  // From GAP 2003-030 [st][pmtab][vemBin]
  vector < vector < vector < double > > > deltaVEM( st2read.size() );
  vector < vector < vector < double > > > deltaVEM2( st2read.size() );
  vector < vector < vector < double > > > avePerDel( st2read.size() );
  vector < vector < vector < int > > > nEvtVEMab( st2read.size() );
  // Charge Vs time, chStPmt[st][pmtab][evt]
  vector < vector < vector < double > > > chStPmt( st2read.size() );
  vector < vector < vector < double > > > chStPmtErr( st2read.size() );
  vector < vector < vector < double > > > chStPmtTime( st2read.size() );
  const int nVemBins = 12;
  const double vemBinWidth = 0.25; 
  
  vector < double > vecVemBins(nVemBins);
  for ( int i=0; i<st2read.size(); i++ ) {
    chStPmt[i].resize(3);
    chStPmtErr[i].resize(3);
    chStPmtTime[i].resize(3);
    deltaVEM[i].resize(3);
    deltaVEM2[i].resize(3);
    avePerDel[i].resize(3);
    nEvtVEMab[i].resize(3);
    for ( int j=0; j<3; j++ ) {
      deltaVEM[i][j].resize(nVemBins);
      deltaVEM2[i][j].resize(nVemBins);
      avePerDel[i][j].resize(nVemBins);
      nEvtVEMab[i][j].resize(nVemBins);
      for ( int jj=1; jj<nVemBins+1; jj++ ) {
        deltaVEM[i][j][jj-1] = 0.;
        deltaVEM2[i][j][jj-1] = 0.;
        avePerDel[i][j][jj-1] = 0.;
        nEvtVEMab[i][j][jj-1] = 0;
        vecVemBins[jj-1] = jj*vemBinWidth;
      }
    } 
  }  

  double sumSigEvt = 0.;
  double tmpSig = 0.;
  double tmpSigErr = 0.;
  int st_posVec = 0;
  bool readSt = false;
  double deltaVEMab = 0.;
  double vemPMTa = 0.;
  double vemPMTb = 0.;
  double vemAve = 0.;
  int bin2vem = 0;
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
        totSigEvt[st_posVec].push_back( tmpSig );
        sigStTimeErr[st_posVec].push_back( tmpSigErr );
        accSigSt[st_posVec].push_back( tmpSigErr/tmpSig );
        timeSigSt[st_posVec].push_back( sdEvent.GetGPSSecond() );
        
        // From GAP 2003-030;
        // pmt1-pmt2
        vemPMTa = sdEvent.GetStation(st_i)->GetPMTTraces(eTotalTrace,1).GetVEMSignal();
        vemPMTb = sdEvent.GetStation(st_i)->GetPMTTraces(eTotalTrace,2).GetVEMSignal();
        vemAve = tmpSig; //(vemPMTa + vemPMTb) / 2.;
        if ( vemPMTa > 1. && vemPMTb > 1. && tmpSig > 1. && tmpSig < 1e3 ) {
        deltaVEMab = ( vemPMTa - vemPMTb ) / sqrt(2.);
        bin2vem = int(log10(vemAve)/vemBinWidth);
        deltaVEM[st_posVec][0][bin2vem] += deltaVEMab;
        deltaVEM2[st_posVec][0][bin2vem] += deltaVEMab*deltaVEMab;
        avePerDel[st_posVec][0][bin2vem] += tmpSig;
        nEvtVEMab[st_posVec][0][bin2vem] += 1;
        }        
        // pmt1-pmt3
        vemPMTa = sdEvent.GetStation(st_i)->GetPMTTraces(eTotalTrace,1).GetVEMSignal();
        vemPMTb = sdEvent.GetStation(st_i)->GetPMTTraces(eTotalTrace,3).GetVEMSignal();
        vemAve = tmpSig; //(vemPMTa + vemPMTb) / 2.;
        if ( vemPMTa > 1. && vemPMTb > 1. && tmpSig > 1. && tmpSig < 1e3 ) {
        deltaVEMab = ( vemPMTa - vemPMTb ) / sqrt(2.);
        bin2vem = int(log10(vemAve)/vemBinWidth);
        deltaVEM[st_posVec][1][bin2vem] += deltaVEMab;
        deltaVEM2[st_posVec][1][bin2vem] += deltaVEMab*deltaVEMab;
        avePerDel[st_posVec][1][bin2vem] += vemAve;
        nEvtVEMab[st_posVec][1][bin2vem] += 1;
        }        
        // pmt3-pmt2
        vemPMTa = sdEvent.GetStation(st_i)->GetPMTTraces(eTotalTrace,2).GetVEMSignal();
        vemPMTb = sdEvent.GetStation(st_i)->GetPMTTraces(eTotalTrace,3).GetVEMSignal();
        vemAve = tmpSig; //(vemPMTa + vemPMTb) / 2.; 
        if ( vemPMTa > 1. && vemPMTb > 1. && tmpSig > 1. && tmpSig < 1e3 ) {
        deltaVEMab = ( vemPMTa - vemPMTb ) / sqrt(2.);        
        bin2vem = int(log10(vemAve)/vemBinWidth);
        deltaVEM[st_posVec][2][bin2vem] += deltaVEMab;
        deltaVEM2[st_posVec][2][bin2vem] += deltaVEMab*deltaVEMab;
        avePerDel[st_posVec][2][bin2vem] += vemAve;
        nEvtVEMab[st_posVec][2][bin2vem] += 1;
        }
                        
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
        totSigEvt[st_posVec].push_back( sumSigEvt );
        timeSig[st_posVec].push_back( sdEvent.GetGPSSecond() );
      }
    }
  }

  // Plotting Signal distribution
  double frtbin = 0;
  double lstbin = 3.6;
  int nBins = (int)(lstbin - frtbin)/0.05;
  for ( int st_i=0; st_i<st2read.size(); st_i++ )
    plotSigDist(st_i, frtbin, lstbin, nBins, totSigEvt, timeSig);

  /*
  TGraphErrors *grpSigSt = new TGraphErrors(timeSigSt[st_selected].size(), 
      &timeSigSt[st_selected].front(), &sigStTime[st_selected].front(), 0, 
      &sigStTimeErr[st_selected].front() );

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

  sigStCanv->Print("recSignalStation.pdf");
  //sigStCanv->Print("stdRecSignalStation.pdf");

  TGraph *grpAccSigSt = new TGraph(timeSigSt[st_selected].size(),
      &timeSigSt[st_selected].front(), &accSigSt[st_selected].front() );

  frtbin = 0.;
  lstbin = 10.;
  nBins = (lstbin - frtbin) / 0.01;
  TH1D *accSgDstStBef = new TH1D("accSgDstStBef","", nBins, frtbin, lstbin);
  TH1D *accSgDstStAft = new TH1D("accSgDstStAft","", nBins, frtbin, lstbin);
  for ( int acc_i=0; acc_i<accSigSt[st_selected].size(); acc_i++ ) {
    if ( timeSigSt[st_selected][acc_i] < loadStart )
      accSgDstStBef->Fill( accSigSt[st_selected][acc_i] );
    else
      accSgDstStAft->Fill( accSigSt[st_selected][acc_i] );
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
  //accSigStCanv->Print("stdRecAccSigSt.pdf");

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

  accSigDistStCanv->Print("recAccDistSigSt.pdf");
  //accSigDistStCanv->Print("stdRecAccDistSigSt.pdf");

  TGraphErrors *grpChStPmt1 = new TGraphErrors(chStPmt[st_selected][0].size(), 
      &chStPmtTime[st_selected][0].front(), &chStPmt[st_selected][0].front(), 0, 
      &chStPmtErr[st_selected][0].front());
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
  //chStCanv->Print("stdRecChargeStation.pdf");
  
  
  vector < vector < double > > sgmVEMpmt(3);
  int tmpNevts = 0;
  double sumDltVem = 0.;
  double sumDltVem2 = 0.;
  double aveVem = 0.;
  double frsTerm = 0.;  
  for ( int pmtab=0; pmtab<3; pmtab++ ) {
    sgmVEMpmt[pmtab].resize(12);
    for ( int i=0; i<nEvtVEMab[st_selected][pmtab].size(); i++ ) {
      tmpNevts = nEvtVEMab[st_selected][pmtab][i];      
      if ( tmpNevts < 2 || !(avePerDel[st_selected][pmtab][i] > 0) ) {
        sgmVEMpmt[pmtab][i] = 0;
        continue;
      }
      sumDltVem = deltaVEM[st_selected][pmtab][i] / tmpNevts;
      sumDltVem2 = deltaVEM2[st_selected][pmtab][i] / tmpNevts;
      aveVem = avePerDel[st_selected][pmtab][i] / tmpNevts;
      frsTerm = sumDltVem2 - (sumDltVem*sumDltVem);
      tmpNevts = tmpNevts / (tmpNevts-1);
      sgmVEMpmt[pmtab][i] = (frsTerm*tmpNevts) / (aveVem);
    }
  }

  TGraph *grpSprdVEMpmt12 = new TGraph(vecVemBins.size(),
      &vecVemBins.front(), &sgmVEMpmt[0].front() );
  TGraph *grpSprdVEMpmt13 = new TGraph(vecVemBins.size(),
      &vecVemBins.front(), &sgmVEMpmt[1].front() );
  TGraph *grpSprdVEMpmt23 = new TGraph(vecVemBins.size(),
      &vecVemBins.front(), &sgmVEMpmt[2].front() );

  TCanvas *spreadVEMcanv = canvasStyle("spreadVEMcanv");
  spreadVEMcanv->cd();
  spreadVEMcanv->SetLogy(1);

  grpSprdVEMpmt12->SetMarkerColor(kBlue);
  grpSprdVEMpmt12->SetMarkerStyle(20);
  grpSprdVEMpmt12->SetMarkerSize(1.5);
  grpSprdVEMpmt12->Draw("AP");

  grpSprdVEMpmt13->SetMarkerColor(kGreen+3);
  grpSprdVEMpmt13->SetMarkerStyle(21);
  grpSprdVEMpmt13->SetMarkerSize(1.5);
  grpSprdVEMpmt13->Draw("P same");
  
  grpSprdVEMpmt23->SetMarkerColor(kRed);
  grpSprdVEMpmt23->SetMarkerStyle(22);
  grpSprdVEMpmt23->SetMarkerSize(1.5);
  grpSprdVEMpmt23->Draw("P same");

  spreadVEMcanv->Print("spreadVEMst.pdf");
  //spreadVEMcanv->Print("stdRecSpreadVEMst.pdf");
  
  */  
  return 0;
}


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

void plotSigDist(int st_slct, double frtbin, double lstbin, int nbins, 
    vector<vector<double>> sigValues, vector<vector<double>> time) {
  cerr << "Doing for Station: " << st2read[st_slct] << endl;
  
  TH1D histDistSigBef ("histDistSigBef", "", nbins, frtbin, lstbin);
  TH1D histDistSigAft ("histDistSigAft", "", nbins, frtbin, lstbin);

  for ( int sig_i=0; sig_i<sigValues[st_slct].size(); sig_i++ )
    if ( time[st_slct][sig_i] < loadStart )
      histDistSigBef.Fill( log10(sigValues[st_slct][sig_i]) );
    else 
      histDistSigAft.Fill( log10(sigValues[st_slct][sig_i]) );  

  Double_t factor = 1.;
  histDistSigBef.Scale(factor/histDistSigBef.GetEntries());
  histDistSigAft.Scale(factor/histDistSigAft.GetEntries());

  TCanvas *sigDistCanv = canvasStyle("sigDistCanv");
  sigDistCanv->cd();
  sigDistCanv->SetLogy(1);

  histDistSigBef.SetStats(kFALSE);
  histDistSigBef.GetYaxis()->SetTitle("a.u.");
  histDistSigBef.GetYaxis()->SetTitleOffset(0.8);
  histDistSigBef.GetYaxis()->SetTitleSize(0.06);
  histDistSigBef.GetYaxis()->SetRangeUser(1e-3, 0.21);
  histDistSigBef.GetYaxis()->SetLabelSize(0.05);
  histDistSigBef.GetXaxis()->SetTitle("log_{10}(S/VEM)");
  histDistSigBef.GetXaxis()->SetTitleOffset(1.);
  histDistSigBef.GetXaxis()->SetTitleSize(0.06);
  histDistSigBef.GetXaxis()->SetLabelSize(0.05);
  histDistSigBef.SetLineColor(kBlue);
  histDistSigBef.Draw();
  histDistSigAft.SetLineColor(kRed);
  histDistSigAft.Draw("same");

  TLegend *lgnd = new TLegend(0.67, 0.45, 0.98, 0.98);
  lgnd->AddEntry(&histDistSigBef, Form("St.: %.f",st2read[st_slct]), "");
  lgnd->AddEntry(&histDistSigBef, "Before Loaded", "l");
  lgnd->AddEntry(&histDistSigBef, 
      Form("Entries %.f",histDistSigBef.GetEntries()), "");
  lgnd->AddEntry(&histDistSigAft, "After Loaded", "l");
  lgnd->AddEntry(&histDistSigAft, 
      Form("Entries %.f",histDistSigAft.GetEntries()), "");
  lgnd->SetTextSize(0.05);
  lgnd->SetBorderSize(0);
  lgnd->SetFillStyle(0);
  lgnd->Draw();
  sigDistCanv->Print(Form("totSigDistSt%.f.pdf",st2read[st_slct]));

  histDistSigBef.Delete();
  histDistSigAft.Delete();
  sigDistCanv->Close();
}
