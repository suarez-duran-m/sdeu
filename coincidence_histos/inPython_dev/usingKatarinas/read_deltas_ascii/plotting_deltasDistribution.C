#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TVector.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include "Math/PdfFuncMathCore.h"

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

void histoStyle(TH1D *hist) {
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
}

Double_t logNormalFunction(Double_t *x, Double_t *par) {
    Double_t num = ROOT::Math::log(x[0] + 5 ) - par[1];
    Double_t den = 2.0 * par[2] * par[2];
    return par[0] * ROOT::Math::exp( -1.*(num*num) / den);
}

void plotting_deltasDistribution() {
    //auto f = TFile::Open("deltasDistribution.root");
    auto f = TFile::Open("deltasDistribution_notTouch.root");
    auto pmt1Dist = (TH1D*)f->Get("deltaDistPmt1");
    auto pmt2Dist = (TH1D*)f->Get("deltaDistPmt2");
    auto pmt3Dist = (TH1D*)f->Get("deltaDistPmt3");
    auto pmt1WeiDist = (TH1D*)f->Get("weighDeltaDistPmt1");
    auto pmt2WeiDist = (TH1D*)f->Get("weighDeltaDistPmt2");
    auto pmt3WeiDist = (TH1D*)f->Get("weighDeltaDistPmt3");
    //    
    double x0LogNormal = -3.;
    double xfLogNormal = 4.4;
    // 
    auto fLogNormal = new TF1("fLogNormal", logNormalFunction, x0LogNormal, xfLogNormal, 3);
    //fLogNormal->SetParameters(1, 1, 1);
    //fLogNormal->SetParameters(3.52793e+02, 1.41706, 3.26202e-01);
    //pmt1Dist->Fit(fLogNormal, "R");
    //pmt3Dist->Fit(fLogNormal, "R+");
    //
    pmt1Dist->Fit("pol2", "pol2", "QR+", x0LogNormal, xfLogNormal);
    pmt2Dist->Fit("pol2", "pol2", "QR+", 0.2, 4.5);
    pmt3Dist->Fit("pol2", "pol2", "QR", x0LogNormal, xfLogNormal);
    //
    pmt1WeiDist->Rebin(6.0);
    pmt3WeiDist->Rebin(6.0);
    pmt1WeiDist->Fit("pol2", "", "QR", x0LogNormal, xfLogNormal);
    pmt2WeiDist->Fit("pol2", "", "QR", 0.2, 4.5);
    pmt3WeiDist->Fit("pol2", "", "QR", x0LogNormal, xfLogNormal);
    //
    auto c0 = canvasStyle("c0");
    c0->cd();
    // 
    //pmt1Dist->GetFunction("fLogNormal")->SetLineColor(kOrange+1);
    pmt1Dist->GetFunction("pol2")->SetLineColor(kBlue);
    pmt1Dist->GetFunction("pol2")->SetLineWidth(0);
    pmt1Dist->Draw();
    //pmt3Dist->GetFunction("fLogNormal")->SetLineColor(kRed);
    pmt3Dist->GetFunction("pol2")->SetLineColor(kRed);
    pmt3Dist->Draw("same");
    //
    /*
    double pkPmt1 = -1.0 * pmt1Dist->GetFunction("fLogNormal")->GetParameter(1) 
        / (2.0 * pmt1Dist->GetFunction("fLogNormal")->GetParameter(2));
    double pkPmt3 = -1.0 * pmt3Dist->GetFunction("fLogNormal")->GetParameter(1) 
        / (2.0 * pmt1Dist->GetFunction("fLogNormal")->GetParameter(2));
    */    
    double derA = pmt1Dist->GetFunction("pol2")->GetParameter(1)     
        / (2.0 * pow(pmt1Dist->GetFunction("pol2")->GetParameter(2), 2));
    double derB = 1. / (2. * pmt1Dist->GetFunction("pol2")->GetParameter(2));
    double pkPmt1 = -1.0 * pmt1Dist->GetFunction("pol2")->GetParameter(1)
        / (2.0 * pmt1Dist->GetFunction("pol2")->GetParameter(2));
    double pkPmt1Err = sqrt(pow((derA * pmt1Dist->GetFunction("pol2")->GetParError(2)), 2) 
        + pow((derB * pmt1Dist->GetFunction("pol2")->GetParError(1)), 2));
    //
    derA = pmt3Dist->GetFunction("pol2")->GetParameter(1) 
        / (2.0 * pow(pmt3Dist->GetFunction("pol2")->GetParameter(2), 2));
    derB = 1. / (2. * pmt3Dist->GetFunction("pol2")->GetParameter(2));    
    double pkPmt3 = -1.0 * pmt3Dist->GetFunction("pol2")->GetParameter(1)
        / (2.0 * pmt3Dist->GetFunction("pol2")->GetParameter(2));
    double pkPmt3Err = sqrt(pow((derA * pmt3Dist->GetFunction("pol2")->GetParError(2)), 2)
        + pow((derB * pmt3Dist->GetFunction("pol2")->GetParError(1)), 2));
    auto lgnd = new TLegend(0.15, 0.7, 0.4, 0.9);
    lgnd->AddEntry(pmt1Dist, Form("Fit pol2 Peak: %.2f #\pm %.2f", pkPmt1, pkPmt1Err), "l");
    lgnd->AddEntry(pmt3Dist, Form("Fit pol2 Peak: %.2f #\pm %.2f", pkPmt3, pkPmt3Err), "l");
    lgnd->AddEntry(pmt1Dist, Form("Diff.: %.2f [%]", pkPmt3 - pkPmt1), "");
    lgnd->AddEntry(pmt1Dist, Form("Rel. Diff.: %.2f [%]", 
        200. * (pkPmt3 - pkPmt1)/(pkPmt3 + pkPmt1)), "");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.04);
    lgnd->Draw();    
    //
    c0->Print("deltasDistribution.pdf");
    //
    //    
    auto c1 = canvasStyle("c1");
    c1->cd();
    //
    pmt1WeiDist->GetFunction("pol2")->SetLineColor(kBlue);
    pmt1WeiDist->SetLineColor(kBlue);
    pmt1WeiDist->GetXaxis()->SetTitle("Delta [%]");
    pmt1WeiDist->GetYaxis()->SetTitle("Counts [au]");    
    pmt1WeiDist->Draw();
    //
    pmt3WeiDist->SetLineColor(kRed);
    pmt3WeiDist->Draw("same");    
    //
    derA = pmt1WeiDist->GetFunction("pol2")->GetParameter(1) 
        / (2.0 * pow(pmt1WeiDist->GetFunction("pol2")->GetParameter(2), 2));
    derB = 1. / (2. * pmt1WeiDist->GetFunction("pol2")->GetParameter(2));
    pkPmt1 = -1.0 * pmt1WeiDist->GetFunction("pol2")->GetParameter(1)
        / (2.0 * pmt1WeiDist->GetFunction("pol2")->GetParameter(2));
    pkPmt1Err = sqrt(pow((derA * pmt1WeiDist->GetFunction("pol2")->GetParError(2)), 2)
        + pow((derB * pmt1WeiDist->GetFunction("pol2")->GetParError(1)), 2));
    //
    derA = pmt3WeiDist->GetFunction("pol2")->GetParameter(1) 
        / (2.0 * pow(pmt3WeiDist->GetFunction("pol2")->GetParameter(2), 2));
    derB = 1. / (2. * pmt3WeiDist->GetFunction("pol2")->GetParameter(2));
    pkPmt3 = -1.0 * pmt3WeiDist->GetFunction("pol2")->GetParameter(1)
        / (2.0 * pmt3WeiDist->GetFunction("pol2")->GetParameter(2));
    pkPmt3Err = sqrt(pow((derA * pmt3WeiDist->GetFunction("pol2")->GetParError(2)), 2)
        + pow((derB * pmt3WeiDist->GetFunction("pol2")->GetParError(1)), 2));
    //    
    lgnd = new TLegend(0.2, 0.7, 0.4, 0.9);
    lgnd->AddEntry(pmt1WeiDist, Form("Fit Peak: %.2f #\pm %.2f", pkPmt1, pkPmt1Err), "l");
    lgnd->AddEntry(pmt3WeiDist, Form("Fit Peak: %.2f #\pm %.2f", pkPmt3, pkPmt3Err), "l");
    lgnd->AddEntry(pmt1WeiDist, Form("Diff.: %.2f [%]", pkPmt3 - pkPmt1), "");
    lgnd->AddEntry(pmt1WeiDist, Form("Rel. Diff.: %.2f [%]", 
        200. * (pkPmt3 - pkPmt1)/(pkPmt3 + pkPmt1)), "");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.04);
    lgnd->Draw();
    //    
    c1->Print("deltasDistribution_weighed.pdf");
    //    
    // 
    // Plotting for PMT2
    auto c2 = canvasStyle("c2");
    c2->cd();
    //
    pmt2Dist->SetLineColor(kOrange+1);
    pmt2Dist->GetFunction("pol2")->SetLineColor(kOrange+1);
    pmt2Dist->GetFunction("pol2")->SetLineWidth(2);
    pmt2Dist->Draw();
    //
    derA = pmt2Dist->GetFunction("pol2")->GetParameter(1) 
        / (2.0 * pow(pmt2Dist->GetFunction("pol2")->GetParameter(2), 2));                     
    derB = 1. / (2. * pmt2Dist->GetFunction("pol2")->GetParameter(2));
    double pkPmt2 = -1.0 * pmt2Dist->GetFunction("pol2")->GetParameter(1)
        / (2.0 * pmt2Dist->GetFunction("pol2")->GetParameter(2));
    double pkPmt2Err = sqrt(pow((derA * pmt2Dist->GetFunction("pol2")->GetParError(2)), 2)
        + pow((derB * pmt2Dist->GetFunction("pol2")->GetParError(1)), 2));
    //
    derA = pmt2Dist->GetFunction("pol2")->GetParameter(1) 
        / (2.0 * pow(pmt2Dist->GetFunction("pol2")->GetParameter(2), 2));
    derB = 1. / (2. * pmt2Dist->GetFunction("pol2")->GetParameter(2));
    pkPmt2 = -1.0 * pmt2Dist->GetFunction("pol2")->GetParameter(1)
        / (2.0 * pmt2Dist->GetFunction("pol2")->GetParameter(2));
    pkPmt2Err = sqrt(pow((derA * pmt2Dist->GetFunction("pol2")->GetParError(2)), 2)
        + pow((derB * pmt2Dist->GetFunction("pol2")->GetParError(1)), 2));
    //    
    lgnd = new TLegend(0.2, 0.7, 0.4, 0.9);
    lgnd->AddEntry(pmt2Dist, Form("Fit Peak: %.2f #\pm %.2f", pkPmt2, pkPmt2Err), "l");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.04);
    lgnd->Draw();
    //    
    c2->Print("deltasDistribution_pmt2.pdf");
    //
    auto c3 = canvasStyle("c3");
    c3->cd();
    //
    pmt2WeiDist->SetLineColor(kOrange+1);
    pmt2WeiDist->GetFunction("pol2")->SetLineColor(kOrange+1);
    pmt2WeiDist->GetFunction("pol2")->SetLineWidth(2);
    pmt2WeiDist->Draw();

    derA = pmt2WeiDist->GetFunction("pol2")->GetParameter(1) 
        / (2.0 * pow(pmt2WeiDist->GetFunction("pol2")->GetParameter(2), 2));
    derB = 1. / (2. * pmt2WeiDist->GetFunction("pol2")->GetParameter(2));
    pkPmt2 = -1.0 * pmt2WeiDist->GetFunction("pol2")->GetParameter(1)
        / (2.0 * pmt2WeiDist->GetFunction("pol2")->GetParameter(2));
    pkPmt2Err = sqrt(pow((derA * pmt2WeiDist->GetFunction("pol2")->GetParError(2)), 2)
        + pow((derB * pmt2WeiDist->GetFunction("pol2")->GetParError(1)), 2));
    //
    derA = pmt2WeiDist->GetFunction("pol2")->GetParameter(1) 
        / (2.0 * pow(pmt2WeiDist->GetFunction("pol2")->GetParameter(2), 2));
    derB = 1. / (2. * pmt2WeiDist->GetFunction("pol2")->GetParameter(2));
    pkPmt2 = -1.0 * pmt2WeiDist->GetFunction("pol2")->GetParameter(1)
        / (2.0 * pmt2WeiDist->GetFunction("pol2")->GetParameter(2));
    pkPmt2Err = sqrt(pow((derA * pmt2WeiDist->GetFunction("pol2")->GetParError(2)), 2)
        + pow((derB * pmt2WeiDist->GetFunction("pol2")->GetParError(1)), 2));
    //    
    lgnd = new TLegend(0.2, 0.7, 0.4, 0.9);
    lgnd->AddEntry(pmt2WeiDist, Form("Fit Peak: %.2f #\pm %.2f", pkPmt2, pkPmt2Err), "l");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.04);
    lgnd->Draw();
    //    
    c3->Print("deltasDistribution_weighed_pmt2.pdf");
    //
    exit(1);
}