#include <TCanvas.h>
#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TGraphErrors.h>
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

Double_t linearFunction(Double_t *x, Double_t *par) {
    return par[0] + par[1]*x[0];
}

Double_t logNormalFunction(Double_t *x, Double_t *par) {
    return par[0]*ROOT::Math::lognormal_pdf(x[0], par[1], par[2]);
}

Double_t fitFunction(Double_t *x, Double_t *par) {
    return linearFunction(x, par) + logNormalFunction(x, &par[2]); //&par[1]); //&par[2]);
}

void fillingDistribution(TH1D *hist, vector<double> deltas) {
    for(int i=0; i<deltas.size(); i++)
        if(deltas[i] > -666.0)
            hist->Fill( deltas[i] );
}

void discartingEmptyStation(vector<double> deltas, vector<double> deltasErr,
    vector<double> vh, vector<double> vhErr, vector<double> &shortDlt, 
    vector<double> &shortDltErr, vector<double> &shortVh, vector<double> 
    &shortVhErr) {
    for(int i=0; i<deltas.size(); i++)    
        if(deltas[i] > -666.0) {
            shortDlt.push_back(deltas[i]);
            shortDltErr.push_back(deltasErr[i]);
            shortVh.push_back(vh[i]);
            shortVhErr.push_back(vhErr[i]);
        }
}

Double_t getErrorFitFunction(double x, double errA[2], double errB[2], double errC[2], 
    double errD[2]) {
    //
    double errX = 0.;
    double piVal = 2.0 * ROOT::Math::Pi();
    double log2Term = pow((log(x) - errD[0]), 2);
    double expTerm = exp( -1.0 * log2Term / ( 2.0 * errC[0] * errC[0]));
    double numTerm = sqrt(errC[0] * piVal);
    double derB = (-1.0 * log2Term * expTerm) / (errC[0] * errC[0] * x * x * numTerm)
        - expTerm / (x * x * numTerm);
    double derC = (errB[0] * piVal * log2Term * expTerm) / (2. * errC[0] * x * x * pow(numTerm, 3))
        + (errB[0] * piVal * expTerm) / (2. * x * x * pow(numTerm, 3))
        - (errB[0] * pow(log2Term, 2) * expTerm) / (pow(errC[0], 5) * x * x * numTerm)
        + (errB[0] * log2Term * expTerm) / (pow(errC[0], 3) * x * x * numTerm);
    double derD = (errB[0] * pow(log2Term, 0.5) * expTerm) / (errC[0] * errC[0] * x * x * numTerm)
        - (errB[0] * pow((log(x) - errD[0]), 3) * expTerm) / (pow(errC[0], 4) * x * x * numTerm);
    //
    errX = sqrt(errA[1]*errA[1] + pow((derB*errB[1]), 2 ) + pow((derC*errC[1]), 2) + pow((derD*errD[1]),2));
    //
    return errX;
}

void fittingCoinciHistos() {
    //
    // Reading the list of label: of gsp and st
    auto treeFileList = new TTree("treeFileList", "treeFileList");
    TString fileListName = "../deltas_python_short.dat";
    //TString fileListName = "../deltas_python.dat";
    //TString fileListName = "list_gpsStId_deltas_short.dat";
    treeFileList->ReadFile(fileListName, "gps/D:stId/I:pmtId/I:cQpk/D:cQpkErr/D:cQpkOff/D:cQpkErrOff/D:Qpk/D:QpkErr/D:QpkOff/D:QpkOffErr/D:delta/D:deltaErr/D:vh/D:vhErr/D");
    //treeFileList->ReadFile(fileListName, "gps/I:stId/I:pmtId/I:cQpk/D");
    int nHisto = treeFileList->Draw("gps:stId:pmtId:cQpk:vh:vhErr", "vh < 0.75", "goff");
    //int nHisto = treeFileList->Draw("gps:stId:pmtId:cQpk", "", "goff");
    //
    // Charging name labels
    double *gpsLabel = treeFileList->GetVal(0);
    double *stIdLabel = treeFileList->GetVal(1);
    double *pmtIdLabel = treeFileList->GetVal(2);
    double *cQpkLabel = treeFileList->GetVal(3);
    double *vhLabel = treeFileList->GetVal(4);
    double *vhLabelErr = treeFileList->GetVal(5);
    // 
    // Vector for Delta storing
    vector < double > chi2FitFcn;
    vector < double > chi2Line;
    vector < double > chi2LogNormal;
    vector < vector < double > > deltasTime(3);
    vector < vector < double > > deltasPmt(3);
    vector < vector < double > > deltasPmtErr(3);
    vector < vector < vector < double > > > deltaPerSt(3);
    vector < vector < vector < double > > > deltaPerStErr(3);
    vector < vector < vector < double > > > vhPerSt(3);
    vector < vector < vector < double > > > vhPerStErr(3);
    for(int i=0; i<3; i++) {
        deltaPerSt[i].resize(2000);
        deltaPerStErr[i].resize(2000);
        vhPerSt[i].resize(2000);
        vhPerStErr[i].resize(2000);
    }
    //
    int firstNhistos = 0;
    auto c0 = canvasStyle("c0");
    TString outPutPdfHistos = "histosFitted.pdf"; //"histosFitted_weirdDelta_2-4.pdf";
    c0->Print(outPutPdfHistos+"(");
    
    //
    // Opening and fitting histos from root files
    for (int histo_i=0; histo_i<nHisto; histo_i++) {
        //
        // Charging histos from the root file, per PMT
        TString filename = "../results/plots/outlier_delta_";
        auto hist_file = TFile::Open(filename+Form("%d_%d_%d.root", 
            (int)gpsLabel[histo_i], (int)stIdLabel[histo_i], 
            (int)pmtIdLabel[histo_i] )
        );
        // 
        // Skipping non-existing files
        if ( hist_file == NULL )
            continue;
        auto cChisto = (TH1D*)hist_file->Get("cch");
        auto qpkVals = (TVectorD*)hist_file->Get("TVectorT<double>;1");
        //
        // Creating function for fit
        int x0Fline = 600; //300;//600;
        int xfFline = 1200; //1000;//1200;
        int x0LogNormal = 1200;
        int xfLogNormal =  9000;//2800;
        int totParameters = 5;
        //auto fLine = new TF1("fLine","[0]", x0Fline, xfFline);
        auto fLine = new TF1("fLine","[0] + [1]*x", x0Fline, xfFline);
        auto fLogNormal = new TF1("fLogNormal","[0]*ROOT::Math::lognormal_pdf(x, [1], [2])", x0LogNormal, xfLogNormal);
        auto fitFcn = new TF1("fitFcn", fitFunction, x0Fline-100, xfLogNormal+100, totParameters);
        //
        // Parameters initialisation
        fLine->SetParameters(1, 1, 1);
        fLogNormal->SetParameters(1, 1, 1);
        //
        // Fitting pol1 and LogNormal independenly
        cChisto->Fit(fLine, "QR");
        cChisto->Fit(fLogNormal, "QR+");
        //
        // Setting parameter for log1+logNormal function
        fitFcn->SetParameters(cChisto->GetFunction("fLine")->GetParameter(0),
            cChisto->GetFunction("fLine")->GetParameter(1),
            cChisto->GetFunction("fLogNormal")->GetParameter(0),
            cChisto->GetFunction("fLogNormal")->GetParameter(1),
            cChisto->GetFunction("fLogNormal")->GetParameter(2)
        );
        //
        // Fitting pol1+logNormal function
        cChisto->Fit(fitFcn, "QR+");
        //
        // Checking for fit quality
        /*
        if ( cChisto->GetFunction("fLine")->GetChisquare() > 18. 
            || cChisto->GetFunction("fLogNormal")->GetChisquare() > 35. 
            || cChisto->GetFunction("fitFcn")->GetChisquare() > 50 ) {
                hist_file->Close();
                continue;
            }
        */
        //
        // Calculating distribution peak from derivative function of fitFcn
        double pkFitFcn = 0.;
        for ( int i=cChisto->GetNbinsX(); i>10; i-- ) {
            if ( fitFcn->Derivative(cChisto->GetBinCenter(i)) > 0 ) {
                pkFitFcn = cChisto->GetBinCenter(i);
                break;
            }
        }
        double deltaFit = 100. * ((pkFitFcn/(*qpkVals)[0] - 1.0 ));
        /*
        if (abs(deltaFit) > 10.)
            continue;
        */
        //
        // Storing Deltas
        double errA[2];
        double errB[2];
        double errC[2];
        double errD[2];
        //
        errA[0] = cChisto->GetFunction("fitFcn")->GetParameter(1);
        errA[1] = cChisto->GetFunction("fitFcn")->GetParError(1);
        errB[0] = cChisto->GetFunction("fitFcn")->GetParameter(2);
        errB[1] = cChisto->GetFunction("fitFcn")->GetParError(2);
        errC[0] = cChisto->GetFunction("fitFcn")->GetParameter(3);
        errC[1] = cChisto->GetFunction("fitFcn")->GetParError(3);
        errD[0] = cChisto->GetFunction("fitFcn")->GetParameter(3); //4);
        errD[1] = cChisto->GetFunction("fitFcn")->GetParError(3); //4);
        double errPk = getErrorFitFunction(pkFitFcn, errA, errB, errC, errD);
        double errDelta = 100. * sqrt(pow((errPk/(*qpkVals)[0]), 2) + pow((pkFitFcn / pow((*qpkVals)[0], 2)), 2));
        deltasTime[(int)pmtIdLabel[histo_i] - 1].push_back( gpsLabel[histo_i] );
        deltasPmt[(int)pmtIdLabel[histo_i] - 1].push_back( deltaFit );
        chi2FitFcn.push_back(cChisto->GetFunction("fitFcn")->GetChisquare());
        chi2Line.push_back(cChisto->GetFunction("fLine")->GetChisquare());
        chi2LogNormal.push_back(cChisto->GetFunction("fLogNormal")->GetChisquare());
        deltaPerSt[(int)pmtIdLabel[histo_i]-1][(int)stIdLabel[histo_i]].push_back(deltaFit);
        deltaPerStErr[(int)pmtIdLabel[histo_i]-1][(int)stIdLabel[histo_i]].push_back(errDelta);
        vhPerSt[(int)pmtIdLabel[histo_i]-1][(int)stIdLabel[histo_i]].push_back(vhLabel[histo_i]);
        vhPerStErr[(int)pmtIdLabel[histo_i]-1][(int)stIdLabel[histo_i]].push_back(vhLabelErr[histo_i]);
        //
        // Plotting outliers
        /*
        if (deltaFit < 2 || deltaFit > 4) {
            hist_file->Close();
            continue;
        }
        */
        firstNhistos++;
        if(firstNhistos > 500)
            continue;
        c0->cd();
        cChisto->SetStats();
        cChisto->GetXaxis()->SetTitle("[FADC]");
        cChisto->GetYaxis()->SetTitle("Counts [au]");
        cChisto->GetFunction("fLine")->SetLineWidth(2);
        cChisto->GetFunction("fLine")->SetLineColor(kGreen+3);
        cChisto->GetFunction("fLogNormal")->SetLineWidth(0);
        cChisto->GetFunction("fitFcn")->SetLineColor(kRed);
        cChisto->GetFunction("fitFcn")->SetLineWidth(3);
        cChisto->Draw();
        //        
        TLegend lgnd(0.6, 0.3, 0.96, 0.7);
        lgnd.AddEntry(cChisto, Form("CQpk from Python: %.2f", cQpkLabel[histo_i]),"");
        lgnd.AddEntry(cChisto, Form("Delta: %.2f", 
            100.*((cQpkLabel[histo_i]/(*qpkVals)[0] - 1.))), "");
        //
        lgnd.AddEntry(cChisto->GetFunction("fitFcn"), Form("CQpk: %.2f", 
            pkFitFcn), "l");
        lgnd.AddEntry(cChisto->GetFunction("fitFcn"), Form("Delta: %.2f", deltaFit, ""), "");
        lgnd.AddEntry(cChisto->GetFunction("fitFcn"), Form("Chi2/NDF: %.2f/%d", 
            cChisto->GetFunction("fitFcn")->GetChisquare(), 
            cChisto->GetFunction("fitFcn")->GetNDF()), "");
        //
        lgnd.AddEntry(cChisto->GetFunction("fLine"), Form("Chi2/NDF: %.2f/%d", 
            cChisto->GetFunction("fLine")->GetChisquare(), 
            cChisto->GetFunction("fLine")->GetNDF()), "");
        //
        lgnd.AddEntry(cChisto->GetFunction("fLogNormal"), Form("Chi2/NDF: %.2f/%d", 
            cChisto->GetFunction("fLogNormal")->GetChisquare(), 
            cChisto->GetFunction("fLogNormal")->GetNDF()), "");
        //
        lgnd.SetBorderSize(0);
        lgnd.SetLineWidth(0);
        lgnd.SetTextSize(0.04);
        lgnd.Draw();
        c0->Print(outPutPdfHistos);
        
        // Closing current root file
        hist_file->Close();
    }
    //
    c0->Print(outPutPdfHistos+")");
    //
    // Creating output root file
    auto outputRoot = TFile::Open("deltasDistribution.root", "recreate");
    //
    // Filling delta distributions
    auto deltaDistPmt1 = new TH1D ("deltaDistPmt1", "", 100, -10, 10);
    auto deltaDistPmt2 = new TH1D ("deltaDistPmt2", "", 100, -10, 10);
    auto deltaDistPmt3 = new TH1D ("deltaDistPmt3", "", 100, -10, 10);
    //
    // Filling
    fillingDistribution(deltaDistPmt1, deltasPmt[0]);
    fillingDistribution(deltaDistPmt2, deltasPmt[1]);
    fillingDistribution(deltaDistPmt3, deltasPmt[2]);
    //
    // Plotting deltas distributions
    auto canvasPMT13 = canvasStyle("pmt13");
    canvasPMT13->cd();
    //
    deltaDistPmt1->GetXaxis()->SetTitle("Delta [%]");
    deltaDistPmt1->GetYaxis()->SetTitle("Counts [au]");    
    deltaDistPmt1->Draw();
    deltaDistPmt3->SetLineColor(kRed);
    deltaDistPmt3->Draw("same");
    //
    // Filling for chi2
    auto distChi2FitFcn = new TH1D ("distChi2FitFcn", "", 3000, 0, 300);
    auto distChi2Line = new TH1D ("distChi2Line", "", 740, 0, 74);
    auto distChi2LogNormal = new TH1D ("distChi2LogNoraml", "", 2100, 0, 210);
    //
    fillingDistribution(distChi2FitFcn, chi2FitFcn);
    fillingDistribution(distChi2Line, chi2Line);
    fillingDistribution(distChi2LogNormal, chi2LogNormal);
    //
    auto canvasChi2 = canvasStyle("canvasChi2");
    canvasChi2->cd();
    distChi2LogNormal->SetLineColor(kRed);
    distChi2LogNormal->Draw();
    distChi2Line->SetLineColor(kGreen+3);
    distChi2Line->Draw("same");
    distChi2LogNormal->SetLineColor(kBlue);
    distChi2LogNormal->Draw("same");
    //
    // Doing weighed deltas
    vector < vector < double > > deltaWeighed(3);
    vector < vector < double > > deltaWeighedErr(3);
    vector < vector < double > > vhWeighed(3);
    vector < vector < double > > vhWeighedErr(3);
    for(int i=0; i<3; i++) {
        deltaWeighed[i].resize(2000);
        deltaWeighedErr[i].resize(2000);
        vhWeighed[i].resize(2000);
        vhWeighedErr[i].resize(2000);
        for(int j=0; j<2000; j++) {
            deltaWeighed[i][j] = -666.0;
            deltaWeighedErr[i][j] = -666.0;
            vhWeighed[i][j] = -666.0;
            vhWeighedErr[i][j] = -666.0;
        }
    }
    //
    for(int pmt_i=0; pmt_i<3; pmt_i++) {
        for(int sd_i=0; sd_i<2000; sd_i++) {
            double numDelta = 0.;
            double denDelta = 0.;
            double numVh = 0.;
            double denVh = 0.;
            if(deltaPerSt[pmt_i][sd_i].size() > 0) {                
                for(int dlt_i=0; dlt_i<deltaPerSt[pmt_i][sd_i].size(); dlt_i++) {
                    numDelta += deltaPerSt[pmt_i][sd_i][dlt_i] / pow(deltaPerStErr[pmt_i][sd_i][dlt_i], 2);
                    denDelta += 1. / pow(deltaPerStErr[pmt_i][sd_i][dlt_i], 2);
                    //
                    numVh += vhPerSt[pmt_i][sd_i][dlt_i] / pow(vhPerStErr[pmt_i][sd_i][dlt_i], 2);
                    denVh += 1. / pow(vhPerStErr[pmt_i][sd_i][dlt_i], 2);
                }
                deltaWeighed[pmt_i][sd_i] = numDelta / denDelta;
                deltaWeighedErr[pmt_i][sd_i] = sqrt(1. / denDelta);
                vhWeighed[pmt_i][sd_i] = numVh / denVh;
                vhWeighedErr[pmt_i][sd_i] = sqrt(1. / denVh);
            }
        }
    }
    //
    // Doing distribution for weighed deltas    
    auto weighDeltaDistPmt1 = new TH1D ("weighDeltaDistPmt1", "", 200, -10, 10);
    auto weighDeltaDistPmt2 = new TH1D ("weighDeltaDistPmt2", "", 200, -10, 10);
    auto weighDeltaDistPmt3 = new TH1D ("weighDeltaDistPmt3", "", 200, -10, 10);
    fillingDistribution(weighDeltaDistPmt1, deltaWeighed[0]);
    fillingDistribution(weighDeltaDistPmt2, deltaWeighed[1]);
    fillingDistribution(weighDeltaDistPmt3, deltaWeighed[2]);
    //
    auto canvasWeigh = canvasStyle("canvasWeigh");
    canvasWeigh->cd();
    weighDeltaDistPmt1->SetLineColor(kBlue);
    weighDeltaDistPmt1->Draw();
    weighDeltaDistPmt3->SetLineColor(kRed);
    weighDeltaDistPmt3->Draw("same");
    //
    //
    // Doing distribution for weighed vh
    auto weighVhDistPmt1 = new TH1D ("weighVhDistPmt1", "", 100, 0, 1);
    auto weighVhDistPmt2 = new TH1D ("weighVhDistPmt2", "", 100, 0, 1);
    auto weighVhDistPmt3 = new TH1D ("weighVhDistPmt3", "", 100, 0, 1);
    fillingDistribution(weighVhDistPmt1, vhWeighed[0]);
    fillingDistribution(weighVhDistPmt2, vhWeighed[1]);
    fillingDistribution(weighVhDistPmt3, vhWeighed[2]);
    //
    auto canvasWeighVh = canvasStyle("canvasWeighVh");
    canvasWeighVh->cd();
    weighVhDistPmt1->SetLineColor(kBlue);
    weighVhDistPmt1->Draw();
    weighVhDistPmt3->SetLineColor(kRed);
    weighVhDistPmt3->Draw("same");
    //
    //
    // Doing Deltas vs VH
    vector < vector < double > > dltWei(3);
    vector < vector < double > > dltWeiErr(3);
    vector < vector < double > > vhWei(3);
    vector < vector < double > > vhWeiErr(3);
    discartingEmptyStation(deltaWeighed[0], deltaWeighedErr[0], vhWeighed[0], vhWeighedErr[0], dltWei[0], dltWeiErr[0], vhWei[0], vhWeiErr[0]);
    discartingEmptyStation(deltaWeighed[1], deltaWeighedErr[1], vhWeighed[1], vhWeighedErr[1], dltWei[1], dltWeiErr[1], vhWei[1], vhWeiErr[1]);
    discartingEmptyStation(deltaWeighed[2], deltaWeighedErr[2], vhWeighed[2], vhWeighedErr[2], dltWei[2], dltWeiErr[2], vhWei[2], vhWeiErr[2]);
    auto deltaVsVhPmt1 = new TGraphErrors (dltWei[0].size(), &vhWei[0].front(), &dltWei[0].front(), &vhWeiErr[0].front(), &dltWeiErr[0].front());
    auto deltaVsVhPmt2 = new TGraphErrors (dltWei[1].size(), &vhWei[1].front(), &dltWei[1].front(), &vhWeiErr[1].front(), &dltWeiErr[1].front());
    auto deltaVsVhPmt3 = new TGraphErrors (dltWei[2].size(), &vhWei[2].front(), &dltWei[2].front(), &vhWeiErr[2].front(), &dltWeiErr[2].front());
    //
    auto canvasDltVh = canvasStyle("canvasDltVh");
    canvasDltVh->cd();
    //    
    deltaVsVhPmt1->Draw("ap");
    //deltaVsVhPmt2->Draw("ap");
    //deltaVsVhPmt3->Draw("ap");
    deltaVsVhPmt1->Write();
    deltaVsVhPmt2->Write();
    deltaVsVhPmt3->Write();
    canvasDltVh->Print("deltasVsVhPmt1.pdf");
    //
    // Doing Deltas Vs time
    auto deltaVsTimePmt1 = new TGraphErrors (nHisto, gpsLabel, cQpkLabel);
    //auto deltaVsTimePmt1 = new TGraphErrors (deltasPmt[0].size(), &deltasTime[0].front(), &deltasPmt[0].front());
    auto deltaVsTimePmt2 = new TGraphErrors (deltasPmt[1].size(), &deltasTime[1].front(), &deltasPmt[1].front());
    auto deltaVsTimePmt3 = new TGraphErrors (deltasPmt[2].size(), &deltasTime[2].front(), &deltasPmt[2].front());
    //
    auto canvasDeltasVsTime = canvasStyle("canvasDeltasVsTime");
    canvasDeltasVsTime->cd();
    //
    deltaVsTimePmt1->SetName("deltaVsTimePmt1");
    deltaVsTimePmt1->SetTitle("");
    deltaVsTimePmt1->GetXaxis()->SetTimeFormat("%m/%d %H");
    deltaVsTimePmt1->GetXaxis()->SetTimeOffset(315964782,"gmt");
    deltaVsTimePmt1->Draw("ap");
    deltaVsTimePmt1->Write();
    canvasDeltasVsTime->Print("deltasVsTimePmt1.pdf");
    //
    deltaVsTimePmt2->SetName("deltaVsTimePmt2");
    deltaVsTimePmt2->SetTitle("");
    deltaVsTimePmt2->GetXaxis()->SetTimeFormat("%m/%d %H");
    deltaVsTimePmt2->GetXaxis()->SetTimeOffset(315964782,"gmt");    
    deltaVsTimePmt2->Draw("ap");
    deltaVsTimePmt2->Write();
    //
    deltaVsTimePmt3->SetName("deltaVsTimePmt3");
    deltaVsTimePmt3->SetTitle("");
    deltaVsTimePmt3->GetXaxis()->SetTimeFormat("%m/%d %H");
    deltaVsTimePmt3->GetXaxis()->SetTimeOffset(315964782,"gmt");    
    deltaVsTimePmt3->Draw("ap");
    deltaVsTimePmt3->Write();
    //
    // Writing and closing output root file
    cout << "MSD, writting and closing" << endl;
    outputRoot->Write();
    outputRoot->Close();
    //
    exit(1);
}