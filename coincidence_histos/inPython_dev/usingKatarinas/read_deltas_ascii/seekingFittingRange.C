#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TString.h>
#include "Math/PdfFuncMathCore.h"
#include <iostream>
#include <sstream>
#include <TLine.h>


/*********************************************************************/
/* Author: Mauricio Suárez Durán                                     */
/* Code to find the best range to perform the fitting of coincidence */
/* histograms by pol1 + LogNormal function                           */
/*********************************************************************/


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

Double_t linearFunction(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0];
}

Double_t logNormalFunction(Double_t *x, Double_t *par) {
  Double_t num = ROOT::Math::log(x[0]) - par[1];
  Double_t den = 2.0 * par[2] * par[2];
	return par[0] * ROOT::Math::exp( -1.*(num*num) / den);
}

Double_t fitFullFunction(Double_t *x, Double_t *par) {
  Double_t f1Pars[2] = {par[0], par[1]};
  Double_t transition = par[2];    
  Double_t f2Pars[3] = {1, par[3], par[4]};
  Double_t c = linearFunction(&transition, f1Pars) / logNormalFunction(&transition, f2Pars);
  Double_t f3Pars[3] = {c, par[3], par[4]};
  //
  return (x[0] < transition) ? linearFunction(x, f1Pars) : logNormalFunction(x, f3Pars);
}

void fillingHistos(vector<vector<double>> chi2, vector<vector<double>> pVal, 
  vector<vector<double>> cQpk, vector<vector<double>> cQpkErr, 
  TH2D *chi2Hist, TH1D *aveChi2Hist, TH2D *logPvalHist, TH1D *aveLogPvalHist, 
  TH2D *cQpkHist, TH1D *aveCQpkHist, TH2D *cQpkErrHist, TH1D *aveCQpkErrHist)
{
  for (unsigned int i = 0; i < chi2.size(); i++)
  {
    if(chi2[i].size() < 1) 
      continue;
    double tmpAveChi2 = 0.;
    double tmpRmsChi2 = 0.;
    double tmpAvePval = 0.;
    double tmpRmsPval = 0.;
    //
    double tmpAveCQpk = 0.;
    double tmpRmsCQpk = 0.;
    double tmpAveCQpkErr = 0.;
    double tmpRmsCQpkErr = 0.;
    //
    unsigned int tmpSize = chi2[i].size();
    for (unsigned int j = 0; j < tmpSize; j++)
    {
      chi2Hist->Fill(i, chi2[i][j]);
      tmpAveChi2 += chi2[i][j];
      tmpRmsChi2 += pow(chi2[i][j], 2);
      //
      logPvalHist->Fill(i, pVal[i][j]);
      tmpAvePval += pVal[i][j];
      tmpRmsPval += pow(pVal[i][j], 2);
      //
      cQpkHist->Fill(i, cQpk[i][j]);
      tmpAveCQpk += cQpk[i][j];
      tmpRmsCQpk += pow(cQpk[i][j], 2);
      //
      cQpkErrHist->Fill(i, cQpkErr[i][j]);
      tmpAveCQpkErr += cQpkErr[i][j];
      tmpRmsCQpkErr += pow(cQpkErr[i][j], 2);
    }
    tmpAveChi2 /= tmpSize;
    tmpRmsChi2 = tmpRmsChi2/tmpSize - pow(tmpAveChi2, 2);
    aveChi2Hist->Fill(i, tmpAveChi2);
    aveChi2Hist->SetBinError(aveChi2Hist->FindBin(i), ROOT::Math::sqrt(tmpRmsChi2));
    //
    tmpAvePval /= tmpSize;
    tmpRmsPval = tmpRmsPval/tmpSize - pow(tmpAvePval, 2);
    aveLogPvalHist->Fill(i, tmpAvePval);
    aveLogPvalHist->SetBinError(aveLogPvalHist->FindBin(i), ROOT::Math::sqrt(tmpRmsPval));
    //
    tmpAveCQpk /= tmpSize;
    tmpRmsCQpk = tmpRmsCQpk/tmpSize - pow(tmpAveCQpk, 2);
    aveCQpkHist->Fill(i, tmpAveCQpk);
    aveCQpkHist->SetBinError(aveCQpkHist->FindBin(i), ROOT::Math::sqrt(tmpRmsCQpk));
    //
    tmpAveCQpkErr /= tmpSize;
    tmpRmsCQpkErr = tmpRmsCQpkErr/tmpSize - pow(tmpAveCQpkErr, 2);
    aveCQpkErrHist->Fill(i, tmpAveCQpkErr);
    aveCQpkErrHist->SetBinError(aveCQpkErrHist->FindBin(i), ROOT::Math::sqrt(tmpRmsCQpkErr));
  }
}

void seekingFittingRange() {
  //
  // Reading the label list: gsp and st
  auto treeFileList = new TTree("treeFileList", "treeFileList");
  TString fileListName = "../deltas_python2.dat";
  treeFileList->ReadFile(fileListName);
  int nHisto = treeFileList->Draw("gps:stId:pmtId","", "goff");
  //
  // Charging name labels
  double *gpsLabel = treeFileList->GetVal(0);
  double *stIdLabel = treeFileList->GetVal(1);
  double *pmtIdLabel = treeFileList->GetVal(2);
  // 
  // Vectors for chi2 and pval storing
  // chi2Ndof[pmt][Ndof][chi2Ndof]
  vector < vector < vector < double > > > chi2Ndof(3);
  vector < vector < vector < double > > > logPval(3);
  vector < vector < vector < double > > > cQpk(3);
  vector < vector < vector < double > > > cQpkErr(3);
  //
  vector < vector < double > > outlier(3);
  vector < vector < double > > totHistoSt(3);
  //
  for(int i=0; i<3; i++) 
  {
    chi2Ndof[i].resize(500);
    logPval[i].resize(500);
    cQpk[i].resize(500);
    cQpkErr[i].resize(500);
    //
    for(int j = 0; j<2000; j++)
    {
      outlier[i].push_back(0);
      totHistoSt[i].push_back(0);
    }
  }
  //
  int totCoincHisto = 0;
  int cntFitCoincHisto = 0;
  //
  auto c0 = canvasStyle("c0");
  TString outPutPdfHistos = "histos_cutParam_outliers_x0-500_cut_pmt1.pdf";
  c0->Print(outPutPdfHistos+"(");
  //
  auto outputInfo = new TFile("kk.root", "RECREATE");
  //auto outputInfo = new TFile("treeFittedParameters_500.root", "RECREATE");
  //auto outputInfo = new TFile("treeFittedParameters2.root", "RECREATE");
  double timeGps;
  double stid;
  double pmtid;
  double poLogNormCQpk;
  double poLogNormCQpkErr;
  double fitPar0;
  double fitParErr0;
  double fitPar1;
  double fitParErr1;
  double fitPar2;
  double fitParErr2;
  double fitPar3;
  double fitParErr3;
  double fitPar4;
  double fitParErr4;
  double cqpkChi2;
  double cqpkNdf;
  //
  outputInfo->cd();
  auto treeParam = new TTree("treeParam","treeNew");
  treeParam->Branch("gpsTime", &timeGps, "gpsTime/D");
  treeParam->Branch("sdId", &stid, "sdId/D");
  treeParam->Branch("pmtId", &pmtid, "pmtId/D");
  treeParam->Branch("poLogNormCQpk", &poLogNormCQpk, "poLogNormCQpk/D");
  treeParam->Branch("poLogNormCQpkErr", &poLogNormCQpkErr, "poLogNormCQpkErr/D");
  treeParam->Branch("poLogNormPar0", &fitPar0, "poLogNormPar0/D");
  treeParam->Branch("poLogNormParErr0", &fitParErr0, "poLogNormParErr0/D");
  treeParam->Branch("poLogNormPar1", &fitPar1, "poLogNormPar1/D");
  treeParam->Branch("poLogNormParErr1", &fitParErr1, "poLogNormParErr1/D");
  treeParam->Branch("poLogNormPar2", &fitPar2, "poLogNormPar2/D");
  treeParam->Branch("poLogNormParErr2", &fitParErr2, "poLogNormParErr2/D");
  treeParam->Branch("poLogNormPar3", &fitPar3, "poLogNormPar3/D");
  treeParam->Branch("poLogNormParErr3", &fitParErr3, "poLogNormParErr3/D");
  treeParam->Branch("poLogNormPar4", &fitPar4, "poLogNormPar4/D");
  treeParam->Branch("poLogNormParErr4", &fitParErr4, "poLogNormParErr4/D");
  treeParam->Branch("poLogNormChi2", &cqpkChi2, "poLogNormChi2/D");
  treeParam->Branch("poLogNormNdf", &cqpkNdf, "poLogNormNdf/D");
  //
  for (int histo_i=0; histo_i<nHisto; histo_i++) {
    //if((int)gpsLabel[histo_i] > 1342137618) // For Feb.
    //if((int)gpsLabel[histo_i] < 1342137618 || (int)gpsLabel[histo_i] > 1343347218) // For Jul.
    //if((int)gpsLabel[histo_i] < 1343347218 || (int)gpsLabel[histo_i] > 1344124818) // For Aug0109.
    //if((int)gpsLabel[histo_i] < 1344124818  || (int)gpsLabel[histo_i] > 1344729618) // For Aug1016.
    //if((int)gpsLabel[histo_i] < 1344729618) // For Oct.
      //continue;
    //if(totCoincHisto > 10)
      //break;
    //
    // Charging histos from the root file, per PMT
	  ostringstream filename;
	  filename << "../results/plots/fittedHisto_delta_" 
      << (int)gpsLabel[histo_i] << "_"
      << (int)stIdLabel[histo_i] << "_"
      << (int)pmtIdLabel[histo_i] << ".root";
    auto hist_file = TFile::Open(filename.str().c_str());
    //
    // Reading root file with CCH and extra data
    totCoincHisto++;
    hist_file->cd();
    auto cChisto = (TH1D*)hist_file->Get("cch");
    totHistoSt[(int)pmtIdLabel[histo_i]-1][(int)stIdLabel[histo_i]]++;
    //
    // Loop for fitting
		const int x0_o = 400; //600;//48;
		const int x0_f = 1200;
    const int xf_o = 3000;
    const int xf_f = 4000;
		const int delta_x0 = 160;
    //for(int x0 = x0_o; x0 < x0_f; x0 += delta_x0) {
    //for(int xf = xf_o; xf < xf_f; xf += delta_x0) {
			//
      // Creating function for fit
	    //auto fitFcn = new TF1("fitFcn", fitFullFunction, x0, xf, 5);
	    //auto fitFcn = new TF1("fitFcn", fitFullFunction, x0_o, xf, 5);
	    auto fitFcn = new TF1("fitFcn", fitFullFunction, x0_o, xf_f, 5);
      //
			// Fittting final function
	    //fitFcn->SetParNames("a", "b", "t", "m", "s");
      fitFcn->SetParameters(4.6, 0.0017, 1045, 7.32, 0.32);
      cChisto->Fit(fitFcn,"QR+");
      double pkFitFcn = ROOT::Math::exp(fitFcn->GetParameter(3));
      double pkFitFcnErr = pkFitFcn * fitFcn->GetParError(3);
      //
      timeGps = (int)gpsLabel[histo_i];
      stid = (int)stIdLabel[histo_i];
      pmtid = (int)pmtIdLabel[histo_i];
      poLogNormCQpk = pkFitFcn;
      poLogNormCQpkErr = pkFitFcnErr;
      fitPar0 = fitFcn->GetParameter(0);
      fitParErr0 = fitFcn->GetParError(0);
      fitPar1 = fitFcn->GetParameter(1);
      fitParErr1 = fitFcn->GetParError(1);
      fitPar2 = fitFcn->GetParameter(2);
      fitParErr2 = fitFcn->GetParError(2);
      fitPar3 = fitFcn->GetParameter(3);
      fitParErr3 = fitFcn->GetParError(3);
      fitPar4 = fitFcn->GetParameter(4);
      fitParErr4 = fitFcn->GetParError(4);
      cqpkChi2 = fitFcn->GetChisquare();
      cqpkNdf = fitFcn->GetNDF();
      outputInfo->cd();
      treeParam->Fill();
      //
      //if(pkFitFcn < fitFcn->GetParameter(2))
        //continue;
      //
      if(fitFcn->GetParameter(0) > 20 || fitFcn->GetParameter(1) > 0.02 
        || fitFcn->GetParameter(2) < 600 || fitFcn->GetParameter(2) > pkFitFcn 
        || fitFcn->GetParameter(3) < 7.0 || fitFcn->GetParameter(4) < 0.2)      
      //if(fitFcn->GetParameter(0) < 20 && fitFcn->GetParameter(1) < 0.02 
      //  && fitFcn->GetParameter(2) > 600 && fitFcn->GetParameter(2) < pkFitFcn 
      //  && fitFcn->GetParameter(3) > 7.0 && fitFcn->GetParameter(4) > 0.2)
      {        
        //if(fitFcn->GetParameter(4) < 0.4)
          //continue;
        c0->cd();
        cChisto->GetXaxis()->SetRangeUser(100, 4500);
        cChisto->Draw();
        auto fLineCQpk = new TLine(pkFitFcn, 0, pkFitFcn, 22);
        fLineCQpk->SetLineColor(kRed);
        fLineCQpk->SetLineWidth(2);
        fLineCQpk->Draw();
        auto fLinePar = new TLine(fitFcn->GetParameter(2), 0, fitFcn->GetParameter(2), 22);
        fLinePar->SetLineColor(kBlack);
        fLinePar->SetLineWidth(2);
        fLinePar->Draw();
        auto lgnd = new TLegend(0.45, 0.35, 0.65, 0.85);
        lgnd->AddEntry(fLineCQpk, Form("CQpk: %.2f", pkFitFcn), "l");
        lgnd->AddEntry(fLinePar, "ParamTransition", "l");        
        lgnd->AddEntry(fitFcn, Form("Parm0: %.2f",fitFcn->GetParameter(0)), "");
        lgnd->AddEntry(fitFcn, Form("Parm1: %.2f",fitFcn->GetParameter(1)), "");
        lgnd->AddEntry(fitFcn, Form("Parm2: %.2f",fitFcn->GetParameter(2)), "");
        lgnd->AddEntry(fitFcn, Form("Parm3: %.2f",fitFcn->GetParameter(3)), "");
        lgnd->AddEntry(fitFcn, Form("Parm4: %.2f",fitFcn->GetParameter(4)), "");
        lgnd->SetBorderSize(0);
        lgnd->SetTextSize(0.06); 
        lgnd->Draw();
        c0->Print(outPutPdfHistos);
      }
      //
      //if(fitFcn->GetParameter(2) < 600) {
        //continue;
      //}
      //
      /*
      if(pkFitFcnErr > 50) 
      {
        outlier[(int)pmtIdLabel[histo_i]-1][(int)stIdLabel[histo_i]]++;
        if(cntFitCoincHisto < 1500)
          c0->Print(outPutPdfHistos);
        //continue;
      }
      */
      //
      int tmpNdof = fitFcn->GetNDF();
      double tmpChi2 = fitFcn->GetChisquare();
      double tmpLog10Pval = TMath::Log10(TMath::Prob(tmpNdof, tmpChi2));
      tmpChi2 /= tmpNdof;
      //
      //tmpNdof = cChisto->FindBin(xf_f) - cChisto->FindBin(x0_o) - 5;
      //
      chi2Ndof[(int)pmtIdLabel[histo_i] - 1][tmpNdof].push_back(tmpChi2);
      logPval[(int)pmtIdLabel[histo_i] - 1][tmpNdof].push_back(tmpLog10Pval);
      cQpk[(int)pmtIdLabel[histo_i] - 1][tmpNdof].push_back(pkFitFcn);
      cQpkErr[(int)pmtIdLabel[histo_i] - 1][tmpNdof].push_back(pkFitFcnErr);
      //
      if(totCoincHisto > 3000)
        break;
      //
    //}
    cntFitCoincHisto++;
  	hist_file->Close();
  }
  outputInfo->cd();
  treeParam->Write();
  outputInfo->Write();
  outputInfo->Close(); 
  //
  c0->Print(outPutPdfHistos+")");
  //
  cout << endl << "MSD coinc histos fitted: " << endl;
  cout << cntFitCoincHisto << " out of " << totCoincHisto << " " 
      << 100.*(1. - double(cntFitCoincHisto) / double(totCoincHisto)) << endl;
  cout << endl;
  //
  // Creating output root file
  //auto outputRoot = TFile::Open("seekRangeLimits_CQerrLow50.root", "recreate");
  auto outputRoot = TFile::Open("seekRangeLimits_cutParamTrans_x0-500.root", "recreate");
  //
  // Filling the TH2D for PMT1
  auto chi2VsNdofPmt1 = new TH2D ("chi2VsNdofPmt1", "#chi^{2}/Ndof", 180, 260, 440, 100, 0, 0.5);
  auto aveChi2VsNdofPmt1 = new TH1D ("aveChi2VsNdofPmt1", "#chi^{2}/Ndof", 180, 260, 440);
  auto logPvalVsNdofPmt1 = new TH2D ("logPvalVsNdofPmt1", "Log(Pval)", 180, 260, 440, 100, -50, 0);
  auto aveLogPvalVsNdofPmt1 = new TH1D ("aveLogPvalVsNdofPmt1", "#chi^{2}/Ndof", 180, 260, 440);
  //
  auto cQpkVsNdofPmt1 = new TH2D ("cQpkVsNdofPmt1", "", 180, 260, 440, 2400, 100, 2500);
  auto aveCQpkVsNdofPmt1 = new TH1D ("aveCQpkVsNdofPmt1", "", 180, 260, 440);
  auto cQpkErrVsNdofPmt1 = new TH2D ("cQpkErrVsNdofPmt1", "", 180, 260, 440, 500, 0, 500);
  auto aveCQpkErrVsNdofPmt1 = new TH1D ("aveCQpkErrVsNdofPmt1", "", 180, 260, 440);
  //
  int doPmt = 0;
  fillingHistos(chi2Ndof[doPmt], logPval[doPmt], cQpk[doPmt], cQpkErr[doPmt], 
    chi2VsNdofPmt1, aveChi2VsNdofPmt1, logPvalVsNdofPmt1, aveLogPvalVsNdofPmt1, 
    cQpkVsNdofPmt1, aveCQpkVsNdofPmt1, cQpkErrVsNdofPmt1, aveCQpkErrVsNdofPmt1);
  //
  // Filling the TH2D for PMT2
  auto chi2VsNdofPmt2 = new TH2D ("chi2VsNdofPmt2", "#chi^{2}/Ndof", 180, 260, 440, 100, 0, 0.5);
  auto aveChi2VsNdofPmt2 = new TH1D ("aveChi2VsNdofPmt2", "#chi^{2}/Ndof", 180, 260, 440);
  auto logPvalVsNdofPmt2 = new TH2D ("logPvalVsNdofPmt2", "Log(Pval)", 180, 260, 440, 100, -50, 0);
  auto aveLogPvalVsNdofPmt2 = new TH1D ("aveLogPvalVsNdofPmt2", "#chi^{2}/Ndof", 180, 260, 440);
  //
  auto cQpkVsNdofPmt2 = new TH2D ("cQpkVsNdofPmt2", "", 180, 260, 440, 2400, 100, 2500);
  auto aveCQpkVsNdofPmt2 = new TH1D ("aveCQpkVsNdofPmt2", "", 180, 260, 440);
  auto cQpkErrVsNdofPmt2 = new TH2D ("cQpkErrVsNdofPmt2", "", 180, 260, 440, 500, 0, 500);
  auto aveCQpkErrVsNdofPmt2 = new TH1D ("aveCQpkErrVsNdofPmt2", "", 180, 260, 440);
  //
  doPmt = 1;
  fillingHistos(chi2Ndof[doPmt], logPval[doPmt], cQpk[doPmt], cQpkErr[doPmt], 
    chi2VsNdofPmt2, aveChi2VsNdofPmt2, logPvalVsNdofPmt2, aveLogPvalVsNdofPmt2, 
    cQpkVsNdofPmt2, aveCQpkVsNdofPmt2, cQpkErrVsNdofPmt2, aveCQpkErrVsNdofPmt2);
  //
  // Filling the TH2D for PMT3
  auto chi2VsNdofPmt3 = new TH2D ("chi2VsNdofPmt3", "#chi^{2}/Ndof", 180, 260, 440, 100, 0, 0.5);
  auto aveChi2VsNdofPmt3 = new TH1D ("aveChi2VsNdofPmt3", "#chi^{2}/Ndof", 180, 260, 440);
  auto logPvalVsNdofPmt3 = new TH2D ("logPvalVsNdofPmt3", "Log(Pval)", 180, 260, 440, 100, -50, 0);
  auto aveLogPvalVsNdofPmt3 = new TH1D ("aveLogPvalVsNdofPmt3", "#chi^{2}/Ndof", 180, 260, 440);
  //
  auto cQpkVsNdofPmt3 = new TH2D ("cQpkVsNdofPmt3", "", 180, 260, 440, 2400, 100, 2500);
  auto aveCQpkVsNdofPmt3 = new TH1D ("aveCQpkVsNdofPmt3", "", 180, 260, 440);
  auto cQpkErrVsNdofPmt3 = new TH2D ("cQpkErrVsNdofPmt3", "", 180, 260, 440, 500, 0, 500);
  auto aveCQpkErrVsNdofPmt3 = new TH1D ("aveCQpkErrVsNdofPmt3", "", 180, 260, 440);
  //
  doPmt = 2;
  fillingHistos(chi2Ndof[doPmt], logPval[doPmt], cQpk[doPmt], cQpkErr[doPmt], 
    chi2VsNdofPmt3, aveChi2VsNdofPmt3, logPvalVsNdofPmt3, aveLogPvalVsNdofPmt3, 
    cQpkVsNdofPmt3, aveCQpkVsNdofPmt3, cQpkErrVsNdofPmt3, aveCQpkErrVsNdofPmt3);
  //
  // Filling outlier
  auto outlierPmt1 = new TH1D("outlierPmt1", "", 2000, 0, 2000);
  auto totHistStPmt1 = new TH1D("totHistStPmt1", "", 2000, 0, 2000);
  auto outlierPmt2 = new TH1D("outlierPmt2", "", 2000, 0, 2000);
  auto totHistStPmt2 = new TH1D("totHistStPmt2", "", 2000, 0, 2000);
  auto outlierPmt3 = new TH1D("outlierPmt3", "", 2000, 0, 2000);
  auto totHistStPmt3 = new TH1D("totHistStPmt3", "", 2000, 0, 2000);
  //
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 2000; j++)
      if(outlier[i][j] > 1) 
        switch (i)
        {
        case 0:
          outlierPmt1->Fill(j, outlier[i][j]);
          totHistStPmt1->Fill(j, totHistoSt[i][j]);
          break;
        case 1:
          outlierPmt2->Fill(j, outlier[i][j]);
          totHistStPmt2->Fill(j, totHistoSt[i][j]);
          break;
        default:
          outlierPmt3->Fill(j, outlier[i][j]);
          totHistStPmt3->Fill(j, totHistoSt[i][j]);
          break;
        }  
  //
  // Writing and closing output root file
  cout << "MSD, writting and closing" << endl;
  outputRoot->Write();
  outputRoot->Close();
  //
  exit(1);
}
