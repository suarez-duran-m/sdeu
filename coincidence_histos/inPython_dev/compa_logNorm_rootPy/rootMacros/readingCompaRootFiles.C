#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TString.h>

#include <sstream>
#include <iostream>

/********************************************************/
/* Author: Mauricio Suárez Durán                        */
/* Code to compare the fit of pol1 + LogNormal function */
/* perform by OffLine/Root and python (Katarina's       */
/* implementation).                                     */
/********************************************************/


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


void plotting(TH2D* hist, TString isErr, TString isPython, double frstBin, double lstBin)
{
	frstBin = (isErr == "fCQpkErr" || isErr == "CQpkErr") ? 0 : frstBin;
	lstBin = (isErr == "CQpkErr") ? 500 : lstBin;
	lstBin = (isErr == "fCQpkErr") ? 150 : lstBin;
	//
	auto c = canvasStyle("c");
	c->cd();
  c->SetRightMargin(0.15);
  //
  hist->SetStats(0);
  hist->SetMarkerStyle(8);
  hist->SetMarkerColor(kBlack);
  //
  hist->GetXaxis()->SetTitle("CQpk (OffLine) [FADC]");
  hist->GetYaxis()->SetTitle("CQpk (" + isPython + ") [FADC]");
  hist->GetZaxis()->SetTitle("Counts [au]");
	if(isErr == "fCQpkErr") 
	{
		hist->GetXaxis()->SetRangeUser(0, lstBin);
		hist->GetYaxis()->SetRangeUser(0, lstBin);
	}
	//
  hist->Draw("colz1");
  //
  auto one2one = new TLine(frstBin, frstBin, lstBin, lstBin);
  one2one->SetLineColor(kBlack);
  one2one->Draw();
  //
	TString titleName = (isPython == "Python") ? "Fit by Pol1 + LogNormal" : "Fit by Pol2 (Y-axis), Pol1 + LogNormal (X-axis)";
	//
  auto lgnd = new TLegend(0.15, 0.70, 0.3, 0.95);
	lgnd->AddEntry(hist, titleName, "");
  lgnd->AddEntry(hist, Form("Entries: %.f", hist->GetEntries()), "");
  lgnd->AddEntry(hist, Form("Mean " + isErr + " (OffLine): %.2f #\pm %.2f", hist->GetMean(1), hist->GetMeanError(1)), "");
  lgnd->AddEntry(hist, Form("Mean " + isErr + " (" + isPython + "): %.2f #\pm %.2f", hist->GetMean(2), hist->GetMeanError(2)), "");
  lgnd->SetFillColorAlpha(kWhite, 1.);
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.04);
  lgnd->Draw();
  //
  c->Print(isErr + "_" + isPython + "2.pdf");
}

void plotting(TH1D* hist, TString isErr, TString isPython, double frstBin, double lstBin)
{
	lstBin = (isErr == "CQpkErr") ? 50. : lstBin;
	auto c1 = canvasStyle("c1");
	c1->cd();
  //
  hist->SetStats(0);
  hist->GetXaxis()->SetTitle("CQpk (" + isPython + ") [FADC]");
	hist->GetYaxis()->SetTitle("Counts [au]");
	if(isErr == "CQpkErr")
		hist->GetXaxis()->SetRangeUser(0, lstBin);
	hist->Draw();
  //
  auto lgnd = new TLegend(0.15, 0.70, 0.3, 0.95);
	if(isErr == "CQpkErr")
		lgnd = new TLegend(0.5, 0.70, 0.75, 0.95);
  lgnd->AddEntry(hist, Form("Entries: %.f", hist->GetEntries()), ""); 
  lgnd->AddEntry(hist, Form("Mean " + isErr + " (" + isPython + "): %.2f #\pm %.2f", hist->GetMean(), hist->GetMeanError()), "");
  lgnd->SetFillColorAlpha(kWhite, 1.);
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.04);
  lgnd->Draw();
  //
  c1->Print("dist_" + isErr + "_" + isPython + "2.pdf");
}

void readingCompaRootFiles()
{
	//
	// Variables to get values from TTree
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
	// Histograms for plotting
  double startBin = 1e2;
  double lastBin = 3e3;
  int nBins = (int)((lastBin - startBin) / 10.);
  //
  double startBinErr = 0;
  double lastBinErr = 5e2;
  int nBinsErr = (int)((lastBinErr - startBinErr) / 2.);
  //
	// For Python comparsion
  auto comPoLogNormOffPyth = new TH2D("comPoLogNormOffPyth", "", nBins, startBin, lastBin, nBins, startBin, lastBin); 
  auto comPoLogNormOffPythErr = new TH2D("comPoLogNormOffPythErr", "", nBinsErr, startBinErr, lastBinErr, nBinsErr, startBinErr, lastBinErr); 
	//
	// For Pol2 comparison
	auto comPoLogNormPol2 = new TH2D("comPoLogNormPol2", "", nBins, startBin, lastBin, nBins, startBin, lastBin);
  auto comPoLogNormPol2Err = new TH2D("comPoLogNormPol2Err", "", nBinsErr, startBinErr, lastBinErr, nBinsErr, startBinErr, lastBinErr);
	//
	// For filter
  auto filterOffPyth = new TH2D("filterOffPyth", "", nBins, startBin, lastBin, nBins, startBin, lastBin); 
  auto filterOffPythErr = new TH2D("filterOffPythErr", "", nBinsErr, startBinErr, lastBinErr, nBinsErr, startBinErr, lastBinErr); 
	auto filterPol2 = new TH2D("filterPol2", "", nBins, startBin, lastBin, nBins, startBin, lastBin);
  auto filterPol2Err = new TH2D("filterPol2Err", "", nBinsErr, startBinErr, lastBinErr, nBinsErr, startBinErr, lastBinErr);
	//
	// For OffLine distribution
	auto distOffLine = new TH1D("distOffLine", "", nBins, startBin, lastBin);
	auto distOffLineErr = new TH1D("distOffLineErr", "", nBinsErr, startBinErr, lastBinErr);
	//
	// For Pol2 distribution
	auto distPol2 = new TH1D("distPol2", "", nBins, startBin, lastBin);
	auto distPol2Err = new TH1D("distPol2Err", "", nBinsErr, startBinErr, lastBinErr);
	//
	// Reading root files
	for (int fi = 0; fi < 5; fi++)
	{	
		ostringstream fileName;
		fileName << "comparison_histos2_" << fi << ".root";
		auto file = TFile::Open(fileName.str().c_str());
		auto tr = (TTree*)file->Get("fitValues");
		//
    tr->SetBranchAddress("gpsTime", &gpsTime);
    tr->SetBranchAddress("sdId", &sdId);
    tr->SetBranchAddress("pmtId", &pmtId);
    //
    tr->SetBranchAddress("qpkOffLine", &qpk_offLine);
    tr->SetBranchAddress("qpkErrOffLine", &qpkErr_offLine);
    tr->SetBranchAddress("qpkPython", &qpk_python);
    tr->SetBranchAddress("qpkErrPython", &qpkErr_python);
    //
    tr->SetBranchAddress("signal", &signal);
    tr->SetBranchAddress("ldf", &ldf);
    tr->SetBranchAddress("spDist", &spDist);
    tr->SetBranchAddress("energy", &energy);
    tr->SetBranchAddress("energyErr", &energyErr);
    tr->SetBranchAddress("zenith", &zenith);
    //
    tr->SetBranchAddress("cqpkOffLine", &cqpk_offLine);
    tr->SetBranchAddress("cqpkErrOffLine", &cqpkErr_offLine);
    tr->SetBranchAddress("cqpkPython", &cqpk_python);
    tr->SetBranchAddress("cqpkErrPython", &cqpkErr_python);
    tr->SetBranchAddress("cqpkPol2", &cqpk_pol2);
    tr->SetBranchAddress("cqpkErrPol2", &cqpkErr_pol2);
		//
		// Filling histos
		for (int ntr = 0; ntr < tr->GetEntries(); ntr++)
		{
			tr->GetEntry(ntr);
    	comPoLogNormOffPyth->Fill(cqpk_offLine, cqpk_python);
    	comPoLogNormOffPythErr->Fill(cqpkErr_offLine, cqpkErr_python);
    	comPoLogNormPol2->Fill(cqpk_offLine, cqpk_pol2);
    	comPoLogNormPol2Err->Fill(cqpkErr_offLine, cqpkErr_pol2);
			//
			distOffLine->Fill(cqpk_offLine);
			distOffLineErr->Fill(cqpkErr_offLine);
			distPol2->Fill(cqpk_pol2);
			distPol2Err->Fill(cqpkErr_pol2);
			//
    	if (cqpk_python > 1e3) // && cqpk_offLine > 1e3 && cqpk_offLine  1600)
			{
				filterOffPyth->Fill(cqpk_offLine, cqpk_python);
				filterOffPythErr->Fill(cqpkErr_offLine, cqpkErr_python);
			}
    	if (cqpk_pol2 > 1e3 && cqpk_offLine > 1100)
			{			
				filterPol2->Fill(cqpk_offLine, cqpk_pol2);
				filterPol2Err->Fill(cqpkErr_offLine, cqpkErr_pol2);
			}
		}
	}
	//
	// Ploting for CQpk
	plotting(comPoLogNormOffPyth, "CQpk", "Python", startBin, lastBin);
	plotting(comPoLogNormPol2, "CQpk", "Pol2", startBin, lastBin);
	//
	plotting(filterOffPyth, "fCQpk", "Python", startBin, lastBin);
	plotting(filterPol2, "fCQpk", "Pol2", startBin, lastBin);
	//
  // Plotting for CQpk error
	plotting(comPoLogNormOffPythErr, "CQpkErr", "Python", startBin, lastBin);
	plotting(comPoLogNormPol2Err, "CQpkErr", "Pol2", startBin, lastBin);
	//
	plotting(filterOffPythErr, "fCQpkErr", "Python", startBin, lastBin);
	plotting(filterPol2Err, "fCQpkErr", "Pol2", startBin, lastBin);
	//
	// Plotting distribution
	plotting(distOffLine, "CQpk", "OffLine", startBin, lastBin);
	plotting(distOffLineErr, "CQpkErr", "OffLine", startBin, lastBin);
	plotting(distPol2, "CQpk", "Pol2", startBin, lastBin);
	plotting(distPol2Err, "CQpkErr", "Pol2", startBin, lastBin);
	//
	exit(1);
}
