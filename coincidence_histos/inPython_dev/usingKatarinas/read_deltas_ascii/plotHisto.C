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

Double_t linearFunction(Double_t *x, Double_t *par) {
    return par[0] + par[1]*x[0];
}

Double_t logNormalFunction(Double_t *x, Double_t *par) {
    return par[0]*ROOT::Math::lognormal_pdf(x[0], par[1], par[2]);
}

Double_t fitFunction(Double_t *x, Double_t *par) {
    return linearFunction(x, par) + logNormalFunction(x, &par[2]);
}

/*  ==================
        Main code
    ==================
*/
void plotHisto() {
    TString gpsName = "1343048931"; //"1328485447"; //"1343239733"; //"1343048931";
    TString stName = "842";//"840";//"803"; //"842";
    TString pmtName = "3";//"1";//"2"; //"3";
    //auto f = TFile::Open("../results/plots/outlier_delta_1343239733_803_2.root"); 
    auto f = TFile::Open("../results/plots/outlier_delta_"+gpsName+"_"+stName+"_"+pmtName+".root");
    auto cChisto = (TH1D*)f->Get("cch");
    auto rChisto = (TH1D*)f->Get("rch");
    auto qpkVals = (TVectorD*)f->Get("TVectorT<double>;1"); 
    qpkVals->Print();

    double xbins[600];
    for (int bin_i=1; bin_i<cChisto->GetNbinsX()+1; bin_i++) {
        xbins[bin_i] = cChisto->GetBinCenter(bin_i)+4;
    }
    xbins[0] = 0;
    auto cChistoDerivative = new TH1D("cChistoDerivative", "", cChisto->GetNbinsX(),
        xbins);
    double der = 0.;
    int h = 0;
    auto cChistoSmooth = (TH1D*)cChisto->Clone("copy");
    cChistoSmooth->Smooth(10000);

    for (int bin_i=1; bin_i<cChistoSmooth->GetNbinsX()-1; bin_i++) {
        h = cChistoSmooth->GetBinCenter(bin_i+1) - cChistoSmooth->GetBinCenter(bin_i);
        der = cChistoSmooth->GetBinContent(bin_i+1) - cChistoSmooth->GetBinContent(bin_i-1);
        der /= 2.*h;
        cChistoDerivative->SetBinContent(bin_i, der);
    }

    auto c0 = canvasStyle("c0");
    c0->cd();

    auto fLine = new TF1("fLine","[0] + [1]*x", 600, 1200);
    auto f1 = new TF1("f1","[0]*ROOT::Math::lognormal_pdf(x,[1],[2])", 1200, 2900);
    auto fitFcn = new TF1("fitFcn", fitFunction, 500, 3000, 5);

    fLine->SetParameters(1, 1);
    f1->SetParameters(1, 1, 1);

    cChisto->Fit(fLine, "R"); //, "", 350, 1000);
    cout << endl << endl;
    cChisto->Fit(f1, "R+"); //, "", 1250, 3000);
    
    fitFcn->SetParameters(cChisto->GetFunction("fLine")->GetParameter(0), 
        cChisto->GetFunction("fLine")->GetParameter(1), 
        cChisto->GetFunction("f1")->GetParameter(0), 
        cChisto->GetFunction("f1")->GetParameter(1),
        cChisto->GetFunction("f1")->GetParameter(2)
        );

    cChisto->Fit(fitFcn, "R+"); //, "", 600, 3000);
    double modeDist = ROOT::Math::exp( cChisto->GetFunction("f1")->GetParameter(1) 
            - pow(cChisto->GetFunction("f1")->GetParameter(2), 2) );

    Double_t xBins[cChisto->GetNbinsX()];
    Double_t yDeri[cChisto->GetNbinsX()];
    double pkFitFcn = 0.;
    bool maxFound = kFALSE;
    for ( int i=1; i<cChisto->GetNbinsX()+1; i++ ){
        xBins[i-1] = cChisto->GetBinCenter(i);
        yDeri[i-1] = fitFcn->Derivative(cChisto->GetBinCenter(i));
        if ( cChisto->GetBinCenter(i) > 1000 && fitFcn->Derivative(cChisto->GetBinCenter(i)) < 0 )
            if ( !maxFound ) {
                pkFitFcn = cChisto->GetBinCenter(i);
                maxFound = kTRUE;
            }
    }

    cChisto->GetYaxis()->SetTitle("Counts [au]");
    cChisto->GetXaxis()->SetTitle("Charge [FADC]");
    histoStyle(cChisto);
    cChisto->GetFunction("f1")->SetLineColor(kGreen+3);
    cChisto->GetFunction("f1")->SetLineWidth(3);
    cChisto->GetFunction("fLine")->SetLineWidth(0);

    cChisto->GetFunction("fitFcn")->SetLineColor(kRed);
    cChisto->GetFunction("fitFcn")->SetLineWidth(3);
    cChisto->Draw(); 

    auto lCQpk = new TLine((*qpkVals)[2], 0, (*qpkVals)[2], 1.04*cChisto->GetMaximum());
    lCQpk->SetLineColor(kBlack);
    lCQpk->SetLineWidth(3);
    lCQpk->Draw();

    auto lMode = new TLine(modeDist, 0, modeDist, 1.08*cChisto->GetMaximum());    
    lMode->SetLineColor(kGreen+3);
    lMode->SetLineWidth(3);
    lMode->Draw();

    auto lPkFcn = new TLine(pkFitFcn, 0., pkFitFcn, 1.04*cChisto->GetMaximum());
    lPkFcn->SetLineColor(kRed);
    lPkFcn->SetLineWidth(3);
    lPkFcn->Draw();

    auto lgnd = new TLegend(0.57, 0.25, 0.97, 0.75);
    lgnd->SetHeader(""); //Form("Entries: %ld", cChisto->GetEntries()));    
    lgnd->AddEntry(lCQpk, "Fit pol2 (Python)", "l");
    lgnd->AddEntry(lCQpk, Form("CQpk: %.2f #\pm %.2f [FADC]", 
        (*qpkVals)[2], (*qpkVals)[3]),"");
    double deltaErr = ROOT::Math::sqrt(pow((*qpkVals)[3], 2)/(pow((*qpkVals)[0], 2)) 
        + (pow((*qpkVals)[2], 2)/pow((*qpkVals)[0], 4)) * pow((*qpkVals)[1], 2));
    lgnd->AddEntry(lCQpk, Form("Delta: %.2f #\pm %.2f [%]", 100.*((*qpkVals)[2]/(*qpkVals)[0] - 1.), deltaErr), "");
    f1->SetLineColor(kGreen+3);
    f1->SetLineWidth(3);
    lgnd->AddEntry(f1, "Fit: LogNormal", "l");
    double modeDistErr = modeDist * ROOT::Math::sqrt(
            pow(cChisto->GetFunction("f1")->GetParError(1), 2)
            + 4. * pow(cChisto->GetFunction("f1")->GetParameter(2), 2)
            * pow(cChisto->GetFunction("f1")->GetParError(2), 2) );
    lgnd->AddEntry(lMode, Form("CQpk: %.2f #\pm %.2f [FADC]", 
        modeDist, modeDistErr),"");
    deltaErr = ROOT::Math::sqrt(pow(modeDistErr, 2)/(pow((*qpkVals)[0], 2)) 
        + (pow(modeDist, 2)/pow((*qpkVals)[0], 4)) * pow((*qpkVals)[1], 2));
    lgnd->AddEntry(lMode, Form("Delta: %.2f #\pm %.2f [%]", 100.*(modeDist/(*qpkVals)[0] - 1.), deltaErr), "");
    fitFcn->SetLineColor(kRed);
    fitFcn->SetLineWidth(3);
    lgnd->AddEntry(fitFcn, "Fit: pol1+LogNormal", "l");
    lgnd->AddEntry(fitFcn, Form("CQpk: %.2f", pkFitFcn), "");
    lgnd->AddEntry(cChisto->GetFunction("fitFcn"), Form("Delta: %.2f", 100.*(pkFitFcn/(*qpkVals)[0] -1.)), "");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);
    lgnd->SetTextSize(0.04);    
    lgnd->Draw();
    
    c0->Print("histo_"+gpsName+"_"+stName+"_"+pmtName+".pdf");
    
    /*
    cChistoSmooth->SetLineColor(kBlack);
    cChistoSmooth->SetLineWidth(2);
    cChistoSmooth->Draw("same");
    */

    /*
    l = new TLine((*qpkVals)[4], 0, (*qpkVals)[4], 15);
    l->SetLineColor(kBlack);
    l->SetLineWidth(2);
    l->Draw();
    */
    
    auto c1 = new TCanvas("c1", "c1");
    c1->cd();

    auto grphDeri = new TGraph(cChisto->GetNbinsX(), xBins, yDeri);
    grphDeri->Draw("al");
    
    auto lQpk = new TLine((*qpkVals)[0], -0.01, (*qpkVals)[0], 0.015);
    lQpk->SetLineColor(kBlack);
    lQpk->Draw();

    auto lPkLogNorm = new TLine(modeDist, -0.01, modeDist, 0.015);
    lPkLogNorm->SetLineColor(kGreen+3);
    lPkLogNorm->Draw();

    lgnd = new TLegend(0.45, 0.4, 0.95, 0.75);
    lgnd->AddEntry(lQpk, Form("Qpk: %.2f", (*qpkVals)[0]), "l");
    lgnd->SetBorderSize(0);
    lgnd->SetLineWidth(0);      
    lgnd->SetTextSize(0.04);    
    lgnd->Draw();
    
    cout << "Mode: " << modeDist << endl;
    cout << "Delta: " << (*qpkVals)[2]/(*qpkVals)[0] - 1 << " " 
        << modeDist/(*qpkVals)[0] - 1 << endl;
}