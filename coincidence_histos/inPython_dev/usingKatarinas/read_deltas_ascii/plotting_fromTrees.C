#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include "Math/PdfFuncMathCore.h"

using namespace std;

Double_t logNormalFunction(Double_t *x, Double_t *par) {
    Double_t num = ROOT::Math::log(x[0]) - par[1];
    Double_t den = 2.0 * par[2] * par[2];
    return par[0] * ROOT::Math::exp( -1.*(num*num) / den);
}

Double_t logNormalFunctionMoved(Double_t *x, Double_t *par) {
    double xnew = x[0] + par[3];
    Double_t num = ROOT::Math::log(xnew) - par[1];
    Double_t den = 2.0 * par[2] * par[2];
    return par[0] * ROOT::Math::exp( -1.*(num*num) / den);
}

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

void fillVectorDeltas(TTree *tr, int isPmt, vector<vector<vector<double>>> &retVec, 
    TH1D *hist) {
    // Setting addresses for input ttree
    double st;
    double pmt;
    double qpk;
    double qpkErr;
    double cqpk;
    double cqpkErr;
    double chi2;
    double ndf;
    tr->SetBranchAddress("sdId", &st);
    tr->SetBranchAddress("pmtId", &pmt);
    tr->SetBranchAddress("qpkPy", &qpk);
    tr->SetBranchAddress("qpkPyErr", &qpkErr);
    tr->SetBranchAddress("poLogNormCQpk", &cqpk);
    tr->SetBranchAddress("poLogNormCQpkErr", &cqpkErr);
    tr->SetBranchAddress("poLogNormChi2", &chi2);
    tr->SetBranchAddress("poLogNormNdf", &ndf);      
    //
    double rangeMax = 0.14 + 0.5*0.013;
    double rangeMin = 0.14 - 0.5*0.013;
    //
    for(int i=0; i<tr->GetEntries(); i++) {
        tr->GetEntry(i);
        if((int)pmt != isPmt)
            continue;
        if(chi2/ndf > rangeMax || chi2/ndf < rangeMin)
            continue;
        //        
        double delta = 100.*(cqpk/qpk - 1.);
        if(delta < -5)
            continue;
        double term1 = (1./qpk)*cqpkErr;
        double term2 = (cqpk/pow(qpk, 2))*qpkErr;
        double deltaErr = 100. * ROOT::Math::sqrt(pow(term1, 2) + pow(term2, 2));        
        retVec[0][(int)st].push_back(delta);
        retVec[1][(int)st].push_back(deltaErr);
        hist->Fill(delta);
    }
    //
    //return retVec;
}

void fillVectorChi2(TTree *tr, int isPmt, TH1D *retHist) {
    double pmt;
    double chi2;
    double ndf;
    tr->SetBranchAddress("pmtId", &pmt);
    tr->SetBranchAddress("poLogNormChi2", &chi2);
    tr->SetBranchAddress("poLogNormNdf", &ndf);
    //
    for(int i=0; i<tr->GetEntries(); i++) {        
        tr->GetEntry(i); 
        if((int)pmt != isPmt)
            continue;
        retHist->Fill(chi2/ndf);
    }
}

TH1D *fillDist(int pmt, vector<vector<vector<double>>> arrayDeltas) {
    TString strPmt = Form("%d", pmt);
    auto retDist = new TH1D("distDeltaStPmt"+strPmt, "", 30, -10, 10);
    for(int st_i=100; st_i<2000; st_i++) {
        double dlt = 0.;
        double dltErr = 0.;
        if (arrayDeltas[0][st_i].size() < 10)
            continue;
        for(int i=0; i<arrayDeltas[0][st_i].size(); i++) {
            dlt += arrayDeltas[0][st_i][i] / pow(arrayDeltas[1][st_i][i], 2);
            dltErr += 1. / pow(arrayDeltas[1][st_i][i], 2);
        }
        retDist->Fill(dlt/dltErr);
    }
    //
    return retDist;
}

void plotting_fromTrees() {
    // Filling vectors with deltas per station
    vector < vector < vector < double > > > deltaStPmt1(2);
    auto deltaDistPmt1 = new TH1D("deltaDistPmt1", "", 50, -10, 20);
    deltaStPmt1[0].resize(2000);
    deltaStPmt1[1].resize(2000);
    vector < vector < vector < double > > > deltaStPmt2(2);
    auto deltaDistPmt2 = new TH1D("deltaDistPmt2", "", 50, -10, 20);
    deltaStPmt2[0].resize(2000);
    deltaStPmt2[1].resize(2000);
    vector < vector < vector < double > > > deltaStPmt3(2);
    auto deltaDistPmt3 = new TH1D("deltaDistPmt3", "", 50, -10, 20);
    deltaStPmt3[0].resize(2000);
    deltaStPmt3[1].resize(2000);
    //
    // TH1D for Chi2/Ndf per PMT
    auto chi2NdfPmt1 = new TH1D("chi2NdfPmt1", "", 1000, 0, 0.5);
    auto chi2NdfPmt2 = new TH1D("chi2NdfPmt2", "", 1000, 0, 0.5);
    auto chi2NdfPmt3 = new TH1D("chi2NdfPmt3", "", 1000, 0, 0.5);
    //
    // Reading file
    auto fIn = TFile::Open("outputFeb.root");
    auto tr = (TTree*)fIn->Get("treeNew");
    //                                           
    fillVectorDeltas(tr, 1, deltaStPmt1, deltaDistPmt1);
    fillVectorDeltas(tr, 2, deltaStPmt2, deltaDistPmt2);
    fillVectorDeltas(tr, 3, deltaStPmt3, deltaDistPmt3);
    tr->Delete();
    tr = (TTree*)fIn->Get("treeNew");
    fillVectorChi2(tr, 1, chi2NdfPmt1);
    fillVectorChi2(tr, 2, chi2NdfPmt2);
    fillVectorChi2(tr, 3, chi2NdfPmt3);
    tr->Delete();
    fIn->Close();
    //    
    fIn = TFile::Open("outputJul.root");
    tr = (TTree*)fIn->Get("treeNew");
    fillVectorDeltas(tr, 1, deltaStPmt1, deltaDistPmt1);
    fillVectorDeltas(tr, 2, deltaStPmt2, deltaDistPmt2);
    fillVectorDeltas(tr, 3, deltaStPmt3, deltaDistPmt3);
    tr->Delete();
    tr = (TTree*)fIn->Get("treeNew");
    fillVectorChi2(tr, 1, chi2NdfPmt1);
    fillVectorChi2(tr, 2, chi2NdfPmt2);
    fillVectorChi2(tr, 3, chi2NdfPmt3);
    tr->Delete();
    fIn->Close();
    //    
    fIn = TFile::Open("outputAug0109.root");
    tr = (TTree*)fIn->Get("treeNew");
    fillVectorDeltas(tr, 1, deltaStPmt1, deltaDistPmt1);
    fillVectorDeltas(tr, 2, deltaStPmt2, deltaDistPmt2);
    fillVectorDeltas(tr, 3, deltaStPmt3, deltaDistPmt3);
    tr->Delete();
    tr = (TTree*)fIn->Get("treeNew");
    fillVectorChi2(tr, 1, chi2NdfPmt1);
    fillVectorChi2(tr, 2, chi2NdfPmt2);
    fillVectorChi2(tr, 3, chi2NdfPmt3);
    tr->Delete();
    fIn->Close();
    //
    fIn = TFile::Open("outputAug1016.root");
    tr = (TTree*)fIn->Get("treeNew");
    fillVectorDeltas(tr, 1, deltaStPmt1, deltaDistPmt1);
    fillVectorDeltas(tr, 2, deltaStPmt2, deltaDistPmt2);
    fillVectorDeltas(tr, 3, deltaStPmt3, deltaDistPmt3);
    tr->Delete();
    tr = (TTree*)fIn->Get("treeNew");    
    fillVectorChi2(tr, 1, chi2NdfPmt1);
    fillVectorChi2(tr, 2, chi2NdfPmt2);
    fillVectorChi2(tr, 3, chi2NdfPmt3);    
    tr->Delete();
    fIn->Close();    
    //
    fIn = TFile::Open("outputOct.root");
    tr = (TTree*)fIn->Get("treeNew");
    fillVectorDeltas(tr, 1, deltaStPmt1, deltaDistPmt1);
    fillVectorDeltas(tr, 2, deltaStPmt2, deltaDistPmt2);
    fillVectorDeltas(tr, 3, deltaStPmt3, deltaDistPmt3);
    tr->Delete();
    tr = (TTree*)fIn->Get("treeNew");    
    fillVectorChi2(tr, 1, chi2NdfPmt1);
    fillVectorChi2(tr, 2, chi2NdfPmt2);
    fillVectorChi2(tr, 3, chi2NdfPmt3);
    tr->Delete();
    fIn->Close();    
    //
    // Filling TH1D with average delta per station
    auto deltaDistStPmt1 = fillDist(1, deltaStPmt1);
    auto deltaDistStPmt2 = fillDist(2, deltaStPmt2);
    auto deltaDistStPmt3 = fillDist(3, deltaStPmt3);
    //
    // Plotting
    auto c0 = canvasStyle("c0");
    c0->cd();
    auto padDistPmt13 = new TPad("padDistPmt13", "padDistPmt13", 0.0, 0.0, 0.49, 1.0);
    padDistPmt13->SetLeftMargin(0.15);
    padDistPmt13->Draw();
    padDistPmt13->cd();
    //
    auto fLogNormDistPmt1 = new TF1("fLogNormDistPmt1", logNormalFunctionMoved, -4, 6, 4);
    fLogNormDistPmt1->SetParameters(2930.70, 1.26, 0.40, 5.40);
    deltaDistPmt1->Fit("fLogNormDistPmt1", "QR");
    double modeDistPmt1 = ROOT::Math::exp(fLogNormDistPmt1->GetParameter(1)) - fLogNormDistPmt1->GetParameter(3);
    double modeDistErrPmt1 = ROOT::Math::sqrt(pow(fLogNormDistPmt1->GetParError(3), 2) 
        + pow(ROOT::Math::exp(fLogNormDistPmt1->GetParameter(1)) * fLogNormDistPmt1->GetParError(1), 2));
    //
    auto fLogNormDistPmt3 = new TF1("fLogNormDistPmt3", logNormalFunctionMoved, -4, 6, 4);    
    fLogNormDistPmt3->SetParameters(2634.84, 1.40, 0.38, 5.35);
    deltaDistPmt3->Fit("fLogNormDistPmt3", "QR");
    double modeDistPmt3 = ROOT::Math::exp(fLogNormDistPmt3->GetParameter(1)) - fLogNormDistPmt3->GetParameter(3);
    double modeDistErrPmt3 = ROOT::Math::sqrt(pow(fLogNormDistPmt3->GetParError(3), 2) 
        + pow(ROOT::Math::exp(fLogNormDistPmt3->GetParameter(1)) * fLogNormDistPmt3->GetParError(1), 2));
    //
    deltaDistPmt1->SetLineColor(kRed);
    deltaDistPmt1->SetLineWidth(2);
    deltaDistPmt1->GetFunction("fLogNormDistPmt1")->SetLineColor(kBlack);
    deltaDistPmt1->GetXaxis()->SetTitle("#Delta [%]");
    deltaDistPmt1->GetYaxis()->SetTitle("Counts [au]");
    deltaDistPmt1->Draw();
    //
    deltaDistPmt3->SetLineColor(kBlue);    
    deltaDistPmt3->SetLineWidth(2);
    deltaDistPmt3->GetFunction("fLogNormDistPmt3")->SetLineColor(kBlack);
    deltaDistPmt3->GetXaxis()->SetTitle("#Delta [%]");
    deltaDistPmt3->GetYaxis()->SetTitle("Counts [au]");
    deltaDistPmt3->Draw("same");
    //
    auto lgnd = new TLegend(0.48, 0.4, 0.89, 0.77);
    lgnd->AddEntry(deltaDistPmt1, "PMT1 Log+Normal fit:", "l");
    lgnd->AddEntry(fLogNormDistPmt1, Form("Mode: %.2f #\pm %.2f [%]", modeDistPmt1, modeDistErrPmt1), "");
    lgnd->AddEntry(fLogNormDistPmt1, Form("#chi^{2}/Ndf: %2.f/%d",
        fLogNormDistPmt1->GetChisquare(), fLogNormDistPmt1->GetNDF()), "");
    lgnd->AddEntry(deltaDistPmt3, "PMT3 Log+Normal fit:", "l");
    lgnd->AddEntry(fLogNormDistPmt3, Form("Mode: %.2f #\pm %.2f [%]", modeDistPmt3, modeDistErrPmt3), "");
    lgnd->AddEntry(fLogNormDistPmt3, Form("#chi^{2}/Ndf: %2.f/%d",
        fLogNormDistPmt3->GetChisquare(), fLogNormDistPmt3->GetNDF()), "");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.03);
    lgnd->Draw();
    //
    c0->cd();
    auto padDistPmt2 = new TPad("padDistPmt2", "padDistPmt2", 0.5, 0.0, 1.0, 1.0);
    padDistPmt2->SetLeftMargin(0.15);
    padDistPmt2->Draw();
    padDistPmt2->cd();
    //
    auto fLogNormDistPmt2 = new TF1("fLogNormDistPmt2", logNormalFunctionMoved, -1, 10, 4);
    fLogNormDistPmt2->SetParameters(2930.70, 1.26, 0.40, 1.0);
    deltaDistPmt2->Fit("fLogNormDistPmt2", "QR");
    double modeDistPmt2 = ROOT::Math::exp(fLogNormDistPmt2->GetParameter(1)) - fLogNormDistPmt2->GetParameter(3);
    double modeDistErrPmt2 = ROOT::Math::sqrt(pow(fLogNormDistPmt2->GetParError(3), 2) 
        + pow(ROOT::Math::exp(fLogNormDistPmt2->GetParameter(1)) * fLogNormDistPmt2->GetParError(1), 2));
    //
    deltaDistPmt2->SetLineColor(kOrange+7);
    deltaDistPmt2->SetLineWidth(2);
    deltaDistPmt2->GetFunction("fLogNormDistPmt2")->SetLineColor(kBlack);
    deltaDistPmt2->GetXaxis()->SetTitle("#Delta [%]");
    deltaDistPmt2->GetYaxis()->SetTitle("Counts [au]");   
    deltaDistPmt2->Draw("same");
    //    
    lgnd = new TLegend(0.5, 0.55, 0.89, 0.77);
    lgnd->AddEntry(deltaDistPmt2, "PMT2 Log+Normal fit:", "l");                                               
    lgnd->AddEntry(fLogNormDistPmt2, Form("Mode: %.2f #\pm %.2f [%]", modeDistPmt2, modeDistErrPmt2), "");
    lgnd->AddEntry(fLogNormDistPmt2, Form("#chi^{2}/Ndf: %2.f/%d",
        fLogNormDistPmt2->GetChisquare(), fLogNormDistPmt2->GetNDF()), "");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.03);
    lgnd->Draw();
    //
    c0->Print("deltaDistPmt.pdf");
    //
    //
    auto c1 = canvasStyle("c1");
    c1->cd();
    //
    // Creating pad for PMT1-3
    auto padPmt13 = new TPad("padPmt13", "padPmt13", 0.01, 0.01, 0.5, 1.);
    padPmt13->Draw();
    padPmt13->cd();    
    // 
    // Plotting for PMT1-3
    auto fLogNormPmt1 = new TF1("fLogNormPmt1", logNormalFunctionMoved, -4, 6, 4);
    fLogNormPmt1->SetParameters(69.08, 1.12, 0.41, 4.96);
    deltaDistStPmt1->Fit("fLogNormPmt1", "QR");
    auto fLogNormPmt3 = new TF1("fLogNormPmt3", logNormalFunctionMoved, -3.5, 6.5, 4);
    fLogNormPmt3->SetParameters(68.07, 0.92, 0.48, 3.89);
    deltaDistStPmt3->Fit("fLogNormPmt3", "QR");    
    //    
    // Calculationg mode
    double modePmt1 = ROOT::Math::exp(fLogNormPmt1->GetParameter(1)) - fLogNormPmt1->GetParameter(3);
    double modeErrPmt1 = ROOT::Math::sqrt(pow(fLogNormPmt1->GetParError(3), 2) 
        + pow(ROOT::Math::exp(fLogNormPmt1->GetParameter(1)) * fLogNormPmt1->GetParError(1), 2));
    double modePmt3 = ROOT::Math::exp(fLogNormPmt3->GetParameter(1)) - fLogNormPmt3->GetParameter(3);
    double modeErrPmt3 = ROOT::Math::sqrt(pow(fLogNormPmt3->GetParError(3), 2) 
        + pow(ROOT::Math::exp(fLogNormPmt3->GetParameter(1)) * fLogNormPmt3->GetParError(1), 2));
    //
    deltaDistStPmt1->GetFunction("fLogNormPmt1")->SetLineColor(kRed);
    deltaDistStPmt1->GetFunction("fLogNormPmt1")->SetLineWidth(2);    
    deltaDistStPmt1->GetXaxis()->SetTitle("#LT #Delta #GT [%]");
    deltaDistStPmt1->GetYaxis()->SetTitle("Counts [au]");
    deltaDistStPmt1->SetLineColor(kRed);
    deltaDistStPmt1->SetLineWidth(2);
    deltaDistStPmt1->Draw();
    //
    deltaDistStPmt3->GetFunction("fLogNormPmt3")->SetLineColor(kBlue);
    deltaDistStPmt3->GetFunction("fLogNormPmt3")->SetLineWidth(2);    
    deltaDistStPmt3->SetLineColor(kBlue);
    deltaDistStPmt3->SetLineWidth(2);    
    deltaDistStPmt3->Draw("same");
    //
    // Doing Legend    
    lgnd = new TLegend(0.52, 0.4, 0.89, 0.77);
    lgnd->AddEntry(deltaDistStPmt1, "PMT1 Log+Normal fit:", "l");
    lgnd->AddEntry(fLogNormPmt1, Form("Mode: %.2f #\pm %.2f [%]", modePmt1, modeErrPmt1), "");
    lgnd->AddEntry(fLogNormPmt1, Form("#chi^{2}/Ndf: %2.f/%d",
        fLogNormPmt1->GetChisquare(), fLogNormPmt1->GetNDF()), "");
    lgnd->AddEntry(deltaDistStPmt3, "PMT3 Log+Normal fit:", "l");
    lgnd->AddEntry(fLogNormPmt3, Form("Mode: %.2f #\pm %.2f [%]", modePmt3, modeErrPmt3), "");
    lgnd->AddEntry(fLogNormPmt3, Form("#chi^{2}/Ndf: %2.f/%d",
        fLogNormPmt3->GetChisquare(), fLogNormPmt3->GetNDF()), "");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.03);
    lgnd->Draw();    
    //
    // Creeating pad for PMT2
    c1->cd();
    auto padPmt2 = new TPad("padPmt2", "padPmt2", 0.5, 0.01, 1., 1.);
    padPmt2->Draw();
    padPmt2->cd();
    //
    auto fLogNormPmt2 = new TF1("fLogNormPmt2", logNormalFunctionMoved, -1.0, 8.0, 4);
    fLogNormPmt2->SetParameters(60.19, 1.97, 0.20, 4.28);
    deltaDistStPmt2->Fit("fLogNormPmt2", "QR");
    double modePmt2 = ROOT::Math::exp(fLogNormPmt2->GetParameter(1)) - fLogNormPmt2->GetParameter(3);
    double modeErrPmt2 = ROOT::Math::sqrt(pow(fLogNormPmt2->GetParError(3), 2));
    //
    fLogNormPmt2->SetLineColor(kBlack);
    deltaDistStPmt2->GetFunction("fLogNormPmt2")->SetLineColor(kBlack);
    deltaDistStPmt2->GetXaxis()->SetTitle("#LT #Delta #GT [%]");
    deltaDistStPmt2->GetYaxis()->SetTitle("Counts [au]");
    deltaDistStPmt2->SetLineColor(kOrange+7);
    deltaDistStPmt2->SetLineWidth(2);
    deltaDistStPmt2->Draw();
    //
    lgnd = new TLegend(0.15, 0.6, 0.4, 0.85);
    lgnd->AddEntry(fLogNormPmt2, "PMT2 Log+Normal fit", "l");
    lgnd->AddEntry(fLogNormPmt2, Form("Mode: %.2f #\pm %.2f [%]", modePmt2, modeErrPmt2), "");
    lgnd->AddEntry(fLogNormPmt2, Form("#chi^{2}/Ndf: %.2f/%d", 
        fLogNormPmt2->GetChisquare(), fLogNormPmt2->GetNDF()), "");
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.03);
    lgnd->Draw();
    //
    // Printing canvas
    c1->Print("deltaStPmt.pdf");
    //
    // Plotting for Chi2/Ndf distribution
    auto c2 = canvasStyle("c2");
    c2->cd();
    auto padChiPmt1 = new TPad("padChiPmt1", "padChiPmt1", 0.01, 0.01, 0.32, 1.0);
    padChiPmt1->Draw();
    padChiPmt1->cd();
    // 
    chi2NdfPmt1->GetXaxis()->SetTitle("#chi^{2}/Ndf [au]");
    chi2NdfPmt1->GetYaxis()->SetTitle("Counts [au]");
    chi2NdfPmt1->Draw();
    //
    c2->cd();
    auto padChiPmt2 = new TPad("padChiPmt2", "padChiPmt2", 0.33, 0.01, 0.65, 1.00);
    padChiPmt2->Draw();
    padChiPmt2->cd();
    // 
    chi2NdfPmt2->GetXaxis()->SetTitle("#chi^{2}/Ndf [au]");
    chi2NdfPmt2->GetYaxis()->SetTitle("Counts [au]");
    chi2NdfPmt2->Draw();
    //
    c2->cd();
    auto padChiPmt3 = new TPad("padChiPmt3", "padChiPmt3", 0.66, 0.01, 1.0, 1.0);
    padChiPmt3->Draw();
    padChiPmt3->cd();
    // 
    chi2NdfPmt3->GetXaxis()->SetTitle("#chi^{2}/Ndf [au]");
    chi2NdfPmt3->GetYaxis()->SetTitle("Counts [au]");
    chi2NdfPmt3->Draw();
    //
    c2->Print("chi2NdfPmt.pdf");
}