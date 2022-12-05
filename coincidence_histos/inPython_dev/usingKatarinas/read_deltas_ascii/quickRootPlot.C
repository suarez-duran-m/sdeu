#include <TPad.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveStats.h>

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

void plotVsPar(TString checVal, TH2D *par0, TH2D *par1, TH2D *par2, TH2D *par3, TH2D *par4)
{
	auto c0 = canvasStyle("c0");
	c0->cd();
	//
 	gStyle->SetStatX(0.9);
   	gStyle->SetStatY(0.95);
	gStyle->SetStatFontSize(0.06);
	gStyle->SetPalette(52);
	//
	// Plotting for Par0
	auto padParm0 = new TPad("padParm0", "", 0.0, 0.5, 0.32, 1.0);
	padParm0->SetLeftMargin(0.15);
	padParm0->Draw();
	padParm0->cd();
	//
	par0->SetMarkerStyle(8);
	par0->SetMarkerColor(kBlack);
	par0->GetXaxis()->SetTitle("Par0 [au]"); 
	par0->GetYaxis()->SetTitle(checVal);
	par0->GetZaxis()->SetRangeUser(0, 350);
	par0->Draw("COLZ1");
	//
	// Plotting for Par1
	c0->cd();
	auto padParm1 = new TPad("padParm1", "", 0.33, 0.5, 0.65, 1.0);
	padParm1->SetLeftMargin(0.15);
	padParm1->Draw();
	padParm1->cd();
	//	
	par1->SetMarkerStyle(8);
	par1->SetMarkerColor(kBlack);
	par1->GetXaxis()->SetTitle("Par1 [au]"); 
	par1->GetYaxis()->SetTitle(checVal);
	par1->Draw("colz1");
	TLatex *lat = new TLatex();
	lat->DrawLatexNDC(.4,.95,"PMT 3");
	//
	// Plotting for Par2
	c0->cd();
	auto padParm2 = new TPad("padParm2", "", 0.66, 0.5, 1.0, 1.0);
	padParm2->SetLeftMargin(0.15);
	padParm2->Draw();
	padParm2->cd();
	//
	par2->SetMarkerStyle(8);
	par2->SetMarkerColor(kBlack);
	par2->GetXaxis()->SetTitle("Par2 [au]"); 
	par2->GetYaxis()->SetTitle(checVal);
	par2->Draw("colz1");
	//
	// Plotting for Par3
	c0->cd();
	auto padParm3 = new TPad("padParm3", "", 0.0, 0.0, 0.49, 0.49);
	padParm3->SetLeftMargin(0.15);
	padParm3->Draw();
	padParm3->cd();
	//
	par3->SetMarkerStyle(8);
	par3->SetMarkerColor(kBlack);
	par3->GetXaxis()->SetTitle("Par3 [au]"); 
	par3->GetYaxis()->SetTitle(checVal);
	par3->Draw("colz1");
	//
	// Plotting for Par4
	c0->cd();
	auto padParm4 = new TPad("padParm4 ", "", 0.5, 0.0, 1.0, 0.49);
	padParm4->SetLeftMargin(0.15);
	padParm4->Draw();
	padParm4->cd();
	//
	par4->SetMarkerStyle(8);
	par4->SetMarkerColor(kBlack);
	par4->GetXaxis()->SetTitle("Par4 [au]"); 
	par4->GetYaxis()->SetTitle(checVal);
	par4->Draw("colz1");
	//
	// Printing canvas
	//c0->Print("QpkErrVsParams_"+checVal+"_pmt3.pdf");
	//c0->Print("QpkErrVsParams_Cut_"+checVal+"_pmt3.pdf");
}

void quickRootPlot()
{
	auto fFile = TFile::Open("treeFittedParameters_500.root");
	//auto fFile = TFile::Open("treeFittedParameters2.root");
	auto fFile2 = TFile::Open("treeFittedParameters3.root");
	auto treeAll = (TTree*)fFile2->Get("treeParam");
	auto treeInfo = (TTree*)fFile->Get("treeParam");
	//
	double minVal = 1e2; //100.;//1000.;
	double maxVal = 2.5e3; //2500.;//2100.;
	int nVals = 2.4e3; //2400.;//150;
	//
	int allTreeVals = treeAll->Draw("poLogNormNdf:poLogNormCQpk", "pmtId > 2 && pmtId < 4", "goff");
	//int allTreeVals = treeAll->Draw("poLogNormNdf:poLogNormCQpk", "pmtId > 0 && pmtId < 2", "goff");
	auto CQpkErrVsNdofAll = new TH2D("CQpkErrVsNdofAll", "", 180, 260, 440, nVals, minVal, maxVal);
	//auto CQpkErrVsNdofAll = new TH2D("CQpkErrVsNdofAll", "", 180, 260, 440, nVals, minVal, maxVal);
	double *ndofAll = treeAll->GetVal(0);
	double *CQpkErrAll = treeAll->GetVal(1);
	//
	for(int i = 0; i < allTreeVals; i++)
		CQpkErrVsNdofAll->Fill(ndofAll[i], CQpkErrAll[i]);
	//
	//int selTreeVals = treeInfo->Draw("poLogNormNdf:poLogNormCQpk:poLogNormPar0:poLogNormPar1:poLogNormPar2:poLogNormPar3:poLogNormPar4", "pmtId > 2 && pmtId < 4", "goff");
	int selTreeVals = treeInfo->Draw("poLogNormNdf:poLogNormCQpk:poLogNormPar0:poLogNormPar1:poLogNormPar2:poLogNormPar3:poLogNormPar4:gpsTime:pmtId", "pmtId > 2 && pmtId < 4", "goff");
	//int selTreeVals = treeInfo->Draw("poLogNormNdf:poLogNormCQpk:poLogNormPar0:poLogNormPar1:poLogNormPar2:poLogNormPar3:poLogNormPar4:gpsTime:pmtId", "pmtId > 2 && pmtId < 4 && poLogNormPar0 < 20 && poLogNormPar1 < 0.02 && poLogNormPar2 > 600 && poLogNormPar2 < poLogNormCQpk && poLogNormPar3 > 7. && poLogNormPar4 > 0.2", "goff");
	//auto CQpkErrVsNdofAll = new TH2D("CQpkErrVsNdofAll", "", 180, 260, 440, 500, 0, 500);
	//
	double *ndof = treeInfo->GetVal(0);
	double *CQpkErr = treeInfo->GetVal(1);
	//
	double *par0 = treeInfo->GetVal(2);
	double *par1 = treeInfo->GetVal(3);
	double *par2 = treeInfo->GetVal(4);
	double *par3 = treeInfo->GetVal(5);
	double *par4 = treeInfo->GetVal(6);
	double *gpsTime = treeInfo->GetVal(7);
	double *pmt = treeInfo->GetVal(8);
	//
	auto CQpkErrVsNdof = new TH2D("CQpkErrVsNdof", "", 180, 260, 440, nVals, minVal, maxVal);
	//auto CQpkErrVsNdof = new TH2D("CQpkErrVsNdof", "", 180, 260, 440, 500, 0, 500);
	/*
	auto Par0VsNdof = new TH2D ("Par0VsNdof", "", 12, -38, 40, nVals, minVal, maxVal);
	auto Par1VsNdof = new TH2D ("Par1VsNdof", "", 30, -0.0100, 0.035, nVals, minVal, maxVal);
	auto Par2VsNdof = new TH2D ("Par2VsNdof", "", 30, 500, 1800, nVals, minVal, maxVal);
	auto Par3VsNdof = new TH2D ("Par3VsNdof", "", 5, 4, 9, nVals, minVal, maxVal);
	auto Par4VsNdof = new TH2D ("Par4VsNdof", "", 30, 0.24, 0.55, nVals, minVal, maxVal);
	*/
	auto Par0VsNdof = new TH2D ("Par0VsNdof", "", 12, -8, 17, nVals, minVal, maxVal);
	auto Par1VsNdof = new TH2D ("Par1VsNdof", "", 30, -0.0100, 0.025, nVals, minVal, maxVal);
	auto Par2VsNdof = new TH2D ("Par2VsNdof", "", 30, 650, 1650, nVals, minVal, maxVal);
	auto Par3VsNdof = new TH2D ("Par3VsNdof", "", 3, 6, 9, nVals, minVal, maxVal);
	auto Par4VsNdof = new TH2D ("Par4VsNdof", "", 30, 0.26, 0.45, nVals, minVal, maxVal);
	//
	for(int i = 0; i < selTreeVals; i++)
	{
		CQpkErrVsNdof->Fill(ndof[i], CQpkErr[i]);
		//
		Par0VsNdof->Fill(par0[i], CQpkErr[i]);
		Par1VsNdof->Fill(par1[i], CQpkErr[i]);
		Par2VsNdof->Fill(par2[i], CQpkErr[i]);
		Par3VsNdof->Fill(par3[i], CQpkErr[i]);
		Par4VsNdof->Fill(par4[i], CQpkErr[i]);
		//
		if(par4[i] > 0.4)
			cout << "MSD " << (int)gpsTime[i] << " " << (int)pmt[i] << " " << par4[i] << " " << CQpkErr[i] << endl;
	} 
	//
	// Plotting
	auto c00 = canvasStyle("c00");
	c00->cd();
	//
	CQpkErrVsNdofAll->SetStats(0);
	CQpkErrVsNdofAll->GetXaxis()->SetTitle("Ndof [au]");
	CQpkErrVsNdofAll->GetYaxis()->SetTitle("CQpk");
	CQpkErrVsNdofAll->SetMarkerStyle(8);
	CQpkErrVsNdofAll->SetMarkerColor(kBlack);
	CQpkErrVsNdofAll->GetXaxis()->SetRangeUser(320, 348);
	CQpkErrVsNdofAll->Draw();
	//
	CQpkErrVsNdof->SetStats(0);
	CQpkErrVsNdof->SetMarkerStyle(8);
	CQpkErrVsNdof->SetMarkerColor(kBlue);
	CQpkErrVsNdof->GetXaxis()->SetRangeUser(320, 348);
	CQpkErrVsNdof->Draw("same");
	//
	auto lgnd = new TLegend(0.14, 0.6, 0.35, 0.95);
	lgnd->AddEntry(CQpkErrVsNdofAll, "PMT3", "");
	lgnd->AddEntry(CQpkErrVsNdofAll, Form("ALL, entries: %.f", CQpkErrVsNdofAll->GetEntries()),
	 "p");
	lgnd->AddEntry(CQpkErrVsNdof, Form("Cut, entries: %.f", CQpkErrVsNdof->GetEntries()), "p");
	lgnd->AddEntry(CQpkErrVsNdof, "par0 < 20 && par1 < 0.02", "");
	lgnd->AddEntry(CQpkErrVsNdof, "par2 > 600 && par2 < CQpk", "");
	lgnd->AddEntry(CQpkErrVsNdof, "par3 > 7.0 && par4 > 0.2", "");
	lgnd->SetBorderSize(0);
	lgnd->SetTextSize(0.06); 
	lgnd->Draw(); 
	//
	//c00->Print("cutParTrans_PMT3.pdf");
	//
	// Plotting
	//plotVsPar("CQpk", Par0VsNdof, Par1VsNdof, Par2VsNdof, Par3VsNdof, Par4VsNdof);
}