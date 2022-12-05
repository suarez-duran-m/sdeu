#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

TCanvas *canvasStyle(TString name) 
{
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

void plotForPmt(int doPmt, TH2D *chi2, TH1D *aveChi2, TH2D *logPval, TH1D *aveLogPval,
	TH2D *cQpk, TH1D *aveCQpk, TH2D *cQpkErr, TH1D *aveCQpkErr)
{
  TString pmtStr = to_string(doPmt);
	//
	const int xRangMin = 340;
	const int xRangMax = 350;
  auto c1 = canvasStyle("c1_PMT"+pmtStr);
	c1->SetTitle("PMT"+pmtStr);
  c1->cd();  
  auto padChi2 = new TPad("padChi2", "padChi2", 0.0, 0.5, 0.49, 1.0);
  padChi2->SetLeftMargin(0.15);
  padChi2->Draw();
  padChi2->cd();  
  //
  chi2->GetXaxis()->SetTitle("Ndof [au]");
  chi2->GetYaxis()->SetTitle("#chi^{2}/Ndof [au]");
  chi2->SetMarkerStyle(15);
	chi2->GetXaxis()->SetRangeUser(xRangMin, xRangMax);
  chi2->Draw();
  //
  aveChi2->SetLineColor(kRed);
  aveChi2->SetLineWidth(2);
  aveChi2->SetMarkerStyle(8);
  aveChi2->SetMarkerSize(1);
  aveChi2->SetMarkerColor(kRed);
  //aveChi2->Draw("same");
  //aveChi2->Draw();
	//
	c1->cd();
	auto padPval = new TPad("padPval", "padPval", 0.5, 0.5, 0.98, 1.0);
  padPval->SetLeftMargin(0.15);
  padPval->Draw();
  padPval->cd();
  //
  logPval->GetXaxis()->SetTitle("Ndof [au]");
  logPval->GetYaxis()->SetTitle("Log(Pval) [au]");
  logPval->SetMarkerStyle(15);
	logPval->GetXaxis()->SetRangeUser(xRangMin, xRangMax);
  logPval->Draw();
  //
  aveLogPval->SetLineColor(kRed);
  aveLogPval->SetLineWidth(2);
  aveLogPval->SetMarkerSize(8);
  aveLogPval->SetMarkerSize(1);
  aveLogPval->SetMarkerColor(kRed);
  //aveLogPval->Draw("same");
  //aveLogPval->Draw();
  //
  c1->cd();
  auto padCQpk = new TPad("padCQpk", "padCQpk", 0.0, 0.0, 0.49, 0.49);
  padCQpk->SetLeftMargin(0.15);
  padCQpk->Draw();
  padCQpk->cd();
  //
  cQpk->GetXaxis()->SetTitle("Ndof [au]");
  cQpk->GetYaxis()->SetTitle("CQ [FADC]");
	cQpk->GetXaxis()->SetRangeUser(xRangMin, xRangMax);
  cQpk->SetMarkerStyle(15);
  cQpk->Draw();
  //
  aveCQpk->SetLineColor(kRed);
  aveCQpk->SetLineWidth(2);
  aveCQpk->SetMarkerStyle(8);
  aveCQpk->SetMarkerSize(1);
  aveCQpk->SetMarkerColor(kRed);
  //aveCQpk->Draw("same");
  //aveCQpk->Draw();
  //
  c1->cd();
  auto padCQpkErr = new TPad("padCQpkErr", "padCQpkErr", 0.5, 0.0, 0.98, 0.49);
  padCQpkErr->SetLeftMargin(0.15);	
  padCQpkErr->Draw();
  padCQpkErr->cd();
  //
  cQpkErr->GetXaxis()->SetTitle("Ndof [au]");
  cQpkErr->GetYaxis()->SetTitle("CQerr [FADC]");
  cQpkErr->SetMarkerStyle(15);
	cQpkErr->GetXaxis()->SetRangeUser(xRangMin, xRangMax);
  cQpkErr->Draw();
  //
  aveCQpkErr->SetLineColor(kRed);
  aveCQpkErr->SetLineWidth(2);
  aveCQpkErr->SetMarkerStyle(8);
  aveCQpkErr->SetMarkerSize(1);
  aveCQpkErr->SetMarkerColor(kRed);
  //aveCQpkErr->Draw("same");
  //aveCQpkErr->Draw();
  //
	// Printting canvas
  c1->Print("chi2PvalCQpkCQpkErr_cutParamTrans_x0_pmt"+pmtStr+".pdf");
  //c1->Print("chi2PvalCQpkCQpkErr_Fixed_pmt"+pmtStr+".pdf");
  //c1->Print("chi2PvalCQpkCQpkErr_Low50_Fixed_pmt"+pmtStr+".pdf");
}

using namespace std;

void plotSeekRangeLimits()
{
	auto inFile = TFile::Open("seekRangeLimits_cutParamTrans_x0.root");
	//auto inFile = TFile::Open("seekRangeLimits_cutParamTrans.root");
	//auto inFile = TFile::Open("seekRangeLimits.root");
	//auto inFile = TFile::Open("seekRangeLimits_fixed.root");
	//auto inFile = TFile::Open("seekRangeLimits_CQerrBig50.root");
	//auto inFile = TFile::Open("seekRangeLimits_CQerrLow50.root");
  //
  // Getting histos for PMT1
  auto chi2VsNdofPmt1 = (TH2D*)inFile->Get("chi2VsNdofPmt1");
  auto aveChi2VsNdofPmt1 = (TH1D*)inFile->Get("aveChi2VsNdofPmt1");
  auto logPvalVsNdofPmt1 = (TH2D*)inFile->Get("logPvalVsNdofPmt1");
  auto aveLogPvalVsNdofPmt1 = (TH1D*)inFile->Get("aveLogPvalVsNdofPmt1");
  //
  auto cQpkVsNdofPmt1 = (TH2D*)inFile->Get("cQpkVsNdofPmt1");
  auto aveCQpkVsNdofPmt1 = (TH1D*)inFile->Get("aveCQpkVsNdofPmt1");
  auto cQpkErrVsNdofPmt1 = (TH2D*)inFile->Get("cQpkErrVsNdofPmt1");
  auto aveCQpkErrVsNdofPmt1 = (TH1D*)inFile->Get("aveCQpkErrVsNdofPmt1");
	//
	// Getting histos for PMT2
  auto chi2VsNdofPmt2 = (TH2D*)inFile->Get("chi2VsNdofPmt2");
  auto aveChi2VsNdofPmt2 = (TH1D*)inFile->Get("aveChi2VsNdofPmt2");
  auto logPvalVsNdofPmt2 = (TH2D*)inFile->Get("logPvalVsNdofPmt2");
  auto aveLogPvalVsNdofPmt2 = (TH1D*)inFile->Get("aveLogPvalVsNdofPmt2");
  //
  auto cQpkVsNdofPmt2 = (TH2D*)inFile->Get("cQpkVsNdofPmt2");
  auto aveCQpkVsNdofPmt2 = (TH1D*)inFile->Get("aveCQpkVsNdofPmt2");
  auto cQpkErrVsNdofPmt2 = (TH2D*)inFile->Get("cQpkErrVsNdofPmt2");
  auto aveCQpkErrVsNdofPmt2 = (TH1D*)inFile->Get("aveCQpkErrVsNdofPmt2");
	//
	// Getting histos for PMT3
  auto chi2VsNdofPmt3 = (TH2D*)inFile->Get("chi2VsNdofPmt3");
  auto aveChi2VsNdofPmt3 = (TH1D*)inFile->Get("aveChi2VsNdofPmt3");
  auto logPvalVsNdofPmt3 = (TH2D*)inFile->Get("logPvalVsNdofPmt3");
  auto aveLogPvalVsNdofPmt3 = (TH1D*)inFile->Get("aveLogPvalVsNdofPmt3");
  //
  auto cQpkVsNdofPmt3 = (TH2D*)inFile->Get("cQpkVsNdofPmt3");
  auto aveCQpkVsNdofPmt3 = (TH1D*)inFile->Get("aveCQpkVsNdofPmt3");
  auto cQpkErrVsNdofPmt3 = (TH2D*)inFile->Get("cQpkErrVsNdofPmt3");	
	auto aveCQpkErrVsNdofPmt3 = (TH1D*)inFile->Get("aveCQpkErrVsNdofPmt3");
	//
	// Plotting for PMT1
	plotForPmt(1, chi2VsNdofPmt1, aveChi2VsNdofPmt1, logPvalVsNdofPmt1, aveLogPvalVsNdofPmt1,
		cQpkVsNdofPmt1, aveCQpkVsNdofPmt1, cQpkErrVsNdofPmt1, aveCQpkErrVsNdofPmt1);
	//
	// Plotting for PMT2
	plotForPmt(2, chi2VsNdofPmt2, aveChi2VsNdofPmt2, logPvalVsNdofPmt2, aveLogPvalVsNdofPmt2,
		cQpkVsNdofPmt2, aveCQpkVsNdofPmt2, cQpkErrVsNdofPmt2, aveCQpkErrVsNdofPmt2);
	//
	// Plotting for PMT3
	plotForPmt(3, chi2VsNdofPmt3, aveChi2VsNdofPmt3, logPvalVsNdofPmt3, aveLogPvalVsNdofPmt3,
		cQpkVsNdofPmt3, aveCQpkVsNdofPmt3, cQpkErrVsNdofPmt3, aveCQpkErrVsNdofPmt3);            	
}