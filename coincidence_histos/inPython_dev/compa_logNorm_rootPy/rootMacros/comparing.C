#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>

#include <iostream>
#include <sstream>

/**************************************************************
 * Author: Mauricio Suárez Durán                              *
 *                                                            *
 * Code to read individual coincidence histograms storing in  *
 * root files. Each of the former, contain the parameters for *
 * the fit of pol1 + LogNormal function performed by          *
 * OffLine/Root and python (Katarina's implementation).       *
 * After read them, a new root file is created containing the *
 * full information of all the histograms read.               *
 *************************************************************/


using namespace std;

TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.11);
  canvas->SetRightMargin(0.03);
  canvas->SetTopMargin(0.02);
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  return canvas; 
 } 

void comparing()
{
  //
  // Charging File IDs
  auto treeFileList = new TTree("treeFileList", "treeFileList");
  TString groupFile = "4";
  TString keysForFiles = "list_keys_rootFiles2_" + groupFile + ".dat";
  treeFileList->ReadFile(keysForFiles);
  int n = treeFileList->Draw("gps:stId:pmtId", "", "goff");
  //
  double *gpsLabel = treeFileList->GetVal(0);
  double *stIdLabel = treeFileList->GetVal(1);
  double *pmtIdLabel = treeFileList->GetVal(2);
  //
  // Creating histos for comparison
  //auto outputRoot = new TFile("kk.root", "recreate");
  auto outputRoot = new TFile("comparison_histos2_" + groupFile + ".root", "recreate");
  //
  double gpsTime = 0.;
  double sdId = 0.;
  double pmtId = 0.;
  //
  double qpk_offLine = 0.;
  double qpkErr_offLine = 0.;
  double qpk_python = 0.;
  double qpkErr_python = 0.;
  //
  double signal = 0.;
  double ldf = 0.;
  double spDist = 0.;
  double energy = 0.;
  double energyErr = 0.;
  double zenith = 0.;
  //
  double cqpk_offLine = 0.;
  double cqpkErr_offLine = 0.;
  double cqpk_python = 0.;
  double cqpkErr_python = 0.;
  double cqpk_pol2 = 0.;
  double cqpkErr_pol2 = 0.;
  //
  outputRoot->cd(); 
  auto outputTree = new TTree("fitValues", "fitValues");
  outputTree->Branch("gpsTime", &gpsTime, "gpsTime/D");
  outputTree->Branch("sdId", &sdId, "sdId/D");
  outputTree->Branch("pmtId", &pmtId, "pmtId/D");
  //
  outputTree->Branch("qpkOffLine", &qpk_offLine, "qpkOffLine/D");
  outputTree->Branch("qpkErrOffLine", &qpkErr_offLine, "qpkErrOffLine/D");
  outputTree->Branch("qpkPython", &qpk_python, "qpkPython/D");
  outputTree->Branch("qpkErrPython", &qpk_python, "qpkErrPython/D");
  //
  outputTree->Branch("signal", &signal, "signal/D");
  outputTree->Branch("ldf", &ldf, "ldf/D");
  outputTree->Branch("spDist", &spDist, "spDist/D");
  outputTree->Branch("energy", &energy, "energy/D");
  outputTree->Branch("energyErr", &energyErr, "energyErr/D");
  outputTree->Branch("zenith", &zenith, "zenith/D");
  //
  outputTree->Branch("cqpkOffLine", &cqpk_offLine, "cqpkOffLine/D");
  outputTree->Branch("cqpkErrOffLine", &cqpkErr_offLine, "cqpkErrOffLine/D");
  outputTree->Branch("cqpkPython", &cqpk_python, "cqpkPython/D");
  outputTree->Branch("cqpkErrPython", &cqpkErr_python, "cqpkErrPython/D");
  outputTree->Branch("cqpkPol2", &cqpk_pol2, "cqpkPol2/D");
  outputTree->Branch("cqpkErrPol2", &cqpkErr_pol2, "cqpkErrPol2/D");
  //
  double startBin = 1e2;
  double lastBin = 3e3;
  int nBins = (int)((lastBin - startBin)/10.);
  //
  double startBinErr = 0;
  double lastBinErr = 5e2;
  int nBinsErr = (int)((lastBinErr - startBinErr)/5.);
  //
  auto compPoLogNormOffPyth = new TH2D("compPoLogNormOffPyth", "", nBins, startBin, lastBin, nBins, startBin, lastBin); 
  auto compPoLogNormOffPythErr = new TH2D("compPoLogNormOffPythErr", "", nBinsErr, startBinErr, lastBinErr, nBinsErr, startBinErr, lastBinErr); 
  auto compPoLogNormPol2 = new TH2D("compPoLogNormPol2", "", nBins, startBin, lastBin, nBins, startBin, lastBin);
  auto compPoLogNormPol2Err = new TH2D("compPoLogNormPol2Err", "", nBinsErr, startBinErr, lastBinErr, nBinsErr, startBinErr, lastBinErr);
  //
  // Reading root files
  //n = 500;
  //TString outputOutLier = "kk_outliers.pdf";
  TString outputOutLier = "outliers2_" + groupFile + ".pdf";
  //TString outputOutLier = "outliers" + groupFile + "_Pol2.pdf";
  auto cPdf = canvasStyle("cPdf");
  cPdf->cd();
  cPdf->Print(outputOutLier + "(");
  //
  for (int hist_i = 0; hist_i < n; hist_i++)
  {
    ostringstream filename;
    filename << "../results/rootFiles2/fittedHisto_delta_"
      << (int)gpsLabel[hist_i] << "_"
      << (int)stIdLabel[hist_i] << "_"
      << (int)pmtIdLabel[hist_i] << ".root";
    //
    // Reading root file
    auto histFile = TFile::Open(filename.str().c_str());
    histFile->cd();
    //
    auto tr = (TTree*)histFile->Get("T");    
    tr->SetBranchAddress("gpsTime", &gpsTime);
    tr->SetBranchAddress("sdId", &sdId);
    tr->SetBranchAddress("pmtId", &pmtId);
    //
    tr->SetBranchAddress("qpkOffLine", &qpk_offLine);
    tr->SetBranchAddress("qpkErrOffLine", &qpkErr_offLine);
    tr->SetBranchAddress("qpk", &qpk_python);
    tr->SetBranchAddress("qpkErr", &qpkErr_python);
    //
    tr->SetBranchAddress("signal", &signal);
    tr->SetBranchAddress("LDF", &ldf);
    tr->SetBranchAddress("spDist", &spDist);
    tr->SetBranchAddress("Energy", &energy);
    tr->SetBranchAddress("EnergyErr", &energyErr);
    tr->SetBranchAddress("Zenith", &zenith);
    //
    tr->SetBranchAddress("cqpkOffLine", &cqpk_offLine);
    tr->SetBranchAddress("cqpkErrOffLine", &cqpkErr_offLine);
    tr->SetBranchAddress("cqpkPoLogNorm", &cqpk_python);
    tr->SetBranchAddress("cqpkErrPoLogNorm", &cqpkErr_python);
    tr->SetBranchAddress("cqpkPy", &cqpk_pol2);
    tr->SetBranchAddress("cqpkErrPy", &cqpkErr_pol2);
    tr->GetEntry(0); 
    /*
    if (cqpk_offLine > 1e3)
    {
      histFile->Close();
      continue;
    }
    if(hist_i > 2000) 
    {
      histFile->Close();
      break; 
    }
    */ 
    //
    outputRoot->cd();
    compPoLogNormOffPyth->Fill(cqpk_offLine, cqpk_python);
    compPoLogNormOffPythErr->Fill(cqpkErr_offLine, cqpkErr_python);
    compPoLogNormPol2->Fill(cqpk_pol2, cqpk_python);
    compPoLogNormPol2Err->Fill(cqpkErr_pol2, cqpkErr_python);
    outputTree->Fill();
    //
    // Checking for outliers
    //if (hist_i < 2000 && cqpk_python < 1000 && cqpk_offLine > 1000 && cqpk_offLine < 1600)
    //if (cqpk_python > 1300 && cqpk_offLine < 1000)
    //if (hist_i < 2000 && cqpk_python < 1000 && cqpk_offLine > 1000 && cqpk_offLine < 1600)
    if (hist_i < 2000 && cqpk_pol2 > 1e3 && cqpk_pol2 < 1.2e3)
    {
      histFile->cd();
      cPdf->cd();
      auto cChisto = (TH1D*)histFile->Get("cch");
      cChisto->SetStats(0);
      cChisto->GetXaxis()->SetTitle("[FADC]");
      cChisto->GetYaxis()->SetTitle("Counts [au]");
      cChisto->Draw();
      //
      auto lOffLine = new TLine(cqpk_offLine, 0, cqpk_offLine, 20);
      lOffLine->SetLineColor(kRed);
      lOffLine->SetLineWidth(2);
      lOffLine->Draw();
      //
      //auto lPython = new TLine(cqpk_python, 0, cqpk_python, 20);
      auto lPython = new TLine(cqpk_python, 0, cqpk_python, 20);
      lPython->SetLineColor(kBlack);
      lPython->SetLineWidth(2);
      lPython->Draw();
      //
      auto lPol2 = new TLine(cqpk_pol2, 0, cqpk_pol2, 20);
      lPol2->SetLineColor(kAzure);
      lPol2->SetLineWidth(2);
      lPol2->Draw();
      //
      auto lgnd = new TLegend(0.4, 0.5, 0.8, 0.9);
      lgnd->AddEntry(lOffLine, Form("CQpk (OffLine): %.2f #\pm %.2f", cqpk_offLine, cqpkErr_offLine), "l");
      lgnd->AddEntry(lPython, Form("CQpk (Python): %.2f #\pm %.2f", cqpk_python, cqpkErr_python), "l");
      lgnd->AddEntry(lPol2, Form("CQpk (Pol2): %2.f #\pm %.2f", cqpk_pol2, cqpkErr_pol2), "l");
      lgnd->SetFillColorAlpha(kWhite, 0.0);
      lgnd->SetBorderSize(0);
      lgnd->SetTextSize(0.05);
      lgnd->Draw();
      //
      cPdf->Print(outputOutLier);
    }
    //
    histFile->Close();
  }
  cPdf->Print(outputOutLier + ")");
  //
  // Writing and closing output root file
  outputRoot->cd();
  outputTree->Write();
  //
  outputRoot->Write();
  outputRoot->Close();
  //
  exit(1);
}
