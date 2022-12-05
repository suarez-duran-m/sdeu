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
#include <iostream>
#include <sstream>


/*************************************************************************/
/* Author: Mauricio Suárez Durán                                         */
/* Code to read individual coincidence histograms storing in root files, */
/* fitting them and store the fit information and new root files; one    */
/* file per month.                                                       */
/*************************************************************************/


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

Double_t logNormalFunctionRoot(Double_t *x, Double_t *par) {
    return par[0]*ROOT::Math::lognormal_pdf(x[0], par[1], par[2]);
}

Double_t logNormalFunction(Double_t *x, Double_t *par) {
    Double_t num = ROOT::Math::log(x[0]) - par[1];
    Double_t den = 2.0 * par[2] * par[2];
    return par[0] * ROOT::Math::exp( -1.*(num*num) / den);
}

Double_t fitFunction(Double_t *x, Double_t *par) {
    return linearFunction(x, par) + logNormalFunction(x, &par[2]);
}

Double_t fitFunctionIoana(Double_t *x, Double_t *par) {
    Double_t f1Pars[2] = {par[0], par[1]};
    Double_t transition = par[2];    
    Double_t f2Pars[3] = {1, par[3], par[4]};
    Double_t c = linearFunction(&transition, f1Pars) / logNormalFunction(&transition, f2Pars);
    Double_t f3Pars[3] = {c, par[3], par[4]};
    //
    return (x[0] < transition) ? linearFunction(x, f1Pars) : logNormalFunction(x, f3Pars);
}

void fillingDistribution(TH1D *hist, vector<double> deltas) {
    for(int i=0; i<deltas.size(); i++)        
        hist->Fill( deltas[i] );
}

void fittingCoinciHistos() {
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
    // Vector for Delta storing
    vector < vector < double > > parameters(5);
    vector < vector < double > > deltasTime(3);
    vector < vector < double > > deltasPmt(3);
    vector < vector < double > > deltasPmtErr(3);
    vector < vector < double > > logPval(3);
    vector < vector < double > > chi2FitFcn(3);
    vector < vector < vector < double > > > deltaPerSt(3);
    vector < vector < vector < double > > > deltaPerStErr(3);
    vector < vector < vector < double > > > vhPerSt(3);
    vector < vector < vector < double > > > vhPerStErr(3);
    //
    // Initializing vectors
    for(int i=0; i<3; i++) {
        deltaPerSt[i].resize(2000);
        deltaPerStErr[i].resize(2000);
        vhPerSt[i].resize(2000);
        vhPerStErr[i].resize(2000);
    }
    //
    int totCoincHisto = 0;
    int cntFitCoincHisto = 0;
    //
    auto c0 = canvasStyle("c0");
    TString outPutPdfHistos = "histosFitted.pdf";
    c0->Print(outPutPdfHistos+"(");
    auto outputInfo = new TFile("outputOct.root", "RECREATE");
    //
    // Variables to read and store        
    double timeGps;
    double stid;
    double pmtid;        
    double qpkPy;
    double qpkPyErr;
    double qpkOff;
    double qpkErrOff;
    double cQpkPy;
    double cQpkErrPy;
    double vh;
    double vhErr;
    double signal;
    double ldf;
    double dist;
    double ener;
    double enerErr;
    double angle;
    // 
    // New variables to store
    double poLogNormCQpk;
    double poLogNormCQpkErr;
    double fitPar1;
    double fitParErr1;
    double fitPar2;
    double fitParErr2;
    double fitPar3;
    double fitParErr3;
    double fitPar4;
    double fitParErr4;
    double fitPar5;
    double fitParErr5;
    double cqpkChi2;
    double cqpkNdf;
    //
    // Creating branches for new TTree in output file    
    outputInfo->cd();
    auto treeNew = new TTree("treeNew","treeNew");
    treeNew->Branch("gpsTime", &timeGps, "gpsTime/D");
    treeNew->Branch("sdId", &stid, "sdId/D");
    treeNew->Branch("pmtId", &pmtid, "pmtId/D");
    treeNew->Branch("qpkPy", &qpkPy, "qpkPy/D");
    treeNew->Branch("qpkPyErr", &qpkPyErr, "qpkErrPy/D");
    treeNew->Branch("qpkOffLine", &qpkOff, "qpkOffLine/D");
    treeNew->Branch("qpkErrOffLine", &qpkErrOff, "qpkErrOffLine/D");
    treeNew->Branch("cqpkPy", &cQpkPy, "cqpkPy/D");
    treeNew->Branch("cqpkErrPy", &cQpkErrPy, "cqpkErrPy/D");        
    treeNew->Branch("poLogNormCQpk", &poLogNormCQpk, "poLogNormCQpk/D");
    treeNew->Branch("poLogNormCQpkErr", &poLogNormCQpkErr, "poLogNormCQpkErr/D");
    treeNew->Branch("poLogNormChi2", &cqpkChi2, "poLogNormChi2/D");
    treeNew->Branch("poLogNormNdf", &cqpkNdf, "poLogNormNdf/D");
    treeNew->Branch("poLogNormPar1", &fitPar1, "poLogNormPar1/D");
    treeNew->Branch("poLogNormParErr1", &fitParErr1, "poLogNormParErr1/D");
    treeNew->Branch("poLogNormPar2", &fitPar2, "poLogNormPar2/D");
    treeNew->Branch("poLogNormParErr2", &fitParErr2, "poLogNormParErr2/D");
    treeNew->Branch("poLogNormPar3", &fitPar3, "poLogNormPar3/D");
    treeNew->Branch("poLogNormParErr3", &fitParErr3, "poLogNormParErr3/D");
    treeNew->Branch("poLogNormPar4", &fitPar4, "poLogNormPar4/D");
    treeNew->Branch("poLogNormParErr4", &fitParErr4, "poLogNormParErr4/D");
    treeNew->Branch("poLogNormPar5", &fitPar5, "poLogNormPar5/D");
    treeNew->Branch("poLogNormParErr5", &fitParErr5, "poLogNormParErr5/D");
    treeNew->Branch("vh", &vh, "vh/D");
    treeNew->Branch("vhErr", &vhErr, "vhErr/D");
    treeNew->Branch("signal", &signal, "signal/D");
    treeNew->Branch("LDF", &ldf, "LDF/D");
    treeNew->Branch("Energy", &ener, "Energy/D");
    treeNew->Branch("EnergyErr", &enerErr, "EnergyErr/D");
    treeNew->Branch("spDist", &dist, "spDist/D");
    treeNew->Branch("Zenith", &angle, "Zenith/D");          
    //
    // Opening and fitting histos from root files
    //
    //nHisto = 10;
    for (int histo_i=0; histo_i<nHisto; histo_i++) {
        //if((int)gpsLabel[histo_i] > 1342137618) // For Feb.
        //if((int)gpsLabel[histo_i] < 1342137618 || (int)gpsLabel[histo_i] > 1343347218) // For Jul.
        //if((int)gpsLabel[histo_i] < 1343347218 || (int)gpsLabel[histo_i] > 1344124818) // For Aug0109.
        //if((int)gpsLabel[histo_i] < 1344124818  || (int)gpsLabel[histo_i] > 1344729618) // For Aug1016.
        if((int)gpsLabel[histo_i] < 1344729618) // For Oct.
            continue;
        //
        // Charging histos from the root file, per PMT
	    ostringstream filename;
	    filename << "../results/plots/fittedHisto_delta_" 
            << (int)gpsLabel[histo_i] << "_"
            << (int)stIdLabel[histo_i] << "_"
            << (int)pmtIdLabel[histo_i] << ".root";
        auto hist_file = TFile::Open(filename.str().c_str());
        //
        // Skipping non-existing files
        if ( hist_file == NULL )
            continue;
        //
        // Reading root file with CCH and extra data
        hist_file->cd();
        auto cChisto = (TH1D*)hist_file->Get("cch");
        //
        // Fitting
        //
        // Creating function for fit
        int x0Fline = 600;
        int xfFline = 1200;
        int x0LogNormal = 1200;
        int xfLogNormal = 4000;
        //
        auto fLine = new TF1("fLine", linearFunction, x0Fline, xfFline, 2);
        auto fLogNormal = new TF1("fLogNormal",logNormalFunction, x0LogNormal, xfLogNormal, 3);        
	    auto fitFcn = new TF1("fitFcn", fitFunctionIoana, x0Fline-100, xfLogNormal+100, 5); 
        //
        // Putting colors
        /*
	    fLine->SetLineColor(kGreen+3);
        fLine->SetLineWidth(3);
        fLogNormal->SetLineColor(kBlack);
        fLogNormal->SetLineWidth(3);
	    fitFcn->SetLineColor(kRed);
	    fitFcn->SetLineWidth(3);
        */
        //	
	    // Parameters initialisation
        fLine->SetParameters(10, 1);
        fLogNormal->SetParameters(15.6, 7.3, 0.3);
	    //
        // Fitting pol1 and LogNormal independenly
        cChisto->Fit(fLine, "QR");
        cChisto->Fit(fLogNormal, "QR+");
        //
        // Fittting final function        
	    fitFcn->SetParNames("a", "b", "t", "n", "m", "s");
        fitFcn->SetParameters(4.6, 0.0017, 1045, 7.32, 0.32);
        //
        //
        cChisto->Fit(fitFcn,"QR+");
        double tmpChi2 = 0.;
        int tmpNdf = 0;
        //
        tmpChi2 = fitFcn->GetChisquare();
        tmpNdf = fitFcn->GetNDF();
        double tmpLog10Pval = TMath::Log10(TMath::Prob(tmpNdf, tmpChi2));
        //
        chi2FitFcn[(int)pmtIdLabel[histo_i] - 1].push_back(tmpChi2 / tmpNdf);
        logPval[(int)pmtIdLabel[histo_i] - 1].push_back(tmpLog10Pval);
        //
        // Calculating distribution mode to get CQpk
        double pkFitFcn = ROOT::Math::exp(fitFcn->GetParameter(3));
        double pkFitFcnErr = pkFitFcn * fitFcn->GetParError(3);
        //
        //
        hist_file->cd();
        auto tree = (TTree*)hist_file->Get("T");
        // 
        hist_file->cd();
        tree->SetBranchAddress("gpsTime", &timeGps);
        tree->SetBranchAddress("sdId", &stid);
        tree->SetBranchAddress("pmtId", &pmtid);        
        tree->SetBranchAddress("qpk", &qpkPy);
        tree->SetBranchAddress("qpkErr", &qpkPyErr);
        tree->SetBranchAddress("qpkOffLine", &qpkOff);
        tree->SetBranchAddress("qpkErrOffLine", &qpkErrOff);
        tree->SetBranchAddress("cqpkPy", &cQpkPy);
        tree->SetBranchAddress("cqpkErrPy", &cQpkErrPy);
        tree->SetBranchAddress("vh", &vh);
        tree->SetBranchAddress("vhErr", &vhErr);
        tree->SetBranchAddress("signal", &signal);
        tree->SetBranchAddress("LDF", &ldf);
        tree->SetBranchAddress("Energy", &ener);
        tree->SetBranchAddress("EnergyErr", &enerErr);
        tree->SetBranchAddress("spDist", &dist);
        tree->SetBranchAddress("Zenith", &angle);        
        tree->GetEntry(0);
        //
        poLogNormCQpk = pkFitFcn;
        poLogNormCQpkErr = pkFitFcnErr;
        cqpkChi2 = tmpChi2;
        cqpkNdf = tmpNdf;
        fitPar1 = fitFcn->GetParameter(0);        
        fitParErr1 = fitFcn->GetParError(0);
        fitPar2 = fitFcn->GetParameter(1);
        fitParErr2 = fitFcn->GetParError(1);
        fitPar3 = fitFcn->GetParameter(2);
        fitParErr3 = fitFcn->GetParError(2);
        fitPar4 = fitFcn->GetParameter(3);
        fitParErr4 = fitFcn->GetParError(4);
        fitPar5 = fitFcn->GetParameter(4);
        fitParErr5 = fitFcn->GetParError(4);
        //
        outputInfo->cd();
        treeNew->Fill();
        //
        hist_file->Close();
    }    
    //
    c0->Print(outPutPdfHistos+")");
    cout << endl << "MSD coinc histos fitted: " << endl;
    cout << cntFitCoincHisto << " out of " << totCoincHisto << " " 
        << 100.*(1. - double(cntFitCoincHisto) / double(totCoincHisto)) << endl;
    cout << endl;
    // 
    // Writing and closing output root file
    outputInfo->cd();   
    treeNew->Write();
    outputInfo->Write();
    outputInfo->Close();
    //
    // Computing average for fit parameters:
    cout << endl << "MSD average for fitFcn parameters" << endl;
    double a[5];
    for(int par_i=0; par_i<5; par_i++) {
        a[par_i] = 0;
        for(int j=0; j<parameters[par_i].size(); j++)
            a[par_i] += parameters[par_i][j];   
        a[par_i] /= parameters[par_i].size();
        cout << a[par_i] << ", ";
    }
    cout << endl;
    //
    // Creating output root file
    auto outputRoot = TFile::Open("deltasDistribution.root", "recreate");
    //
    // Filling delta distributions
    auto deltaDistPmt1 = new TH1D ("deltaDistPmt1", "", 200, -20, 20);
    auto deltaDistPmt2 = new TH1D ("deltaDistPmt2", "", 200, -20, 20);
    auto deltaDistPmt3 = new TH1D ("deltaDistPmt3", "", 200, -20, 20);
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
    auto distChi2FitFcnPmt1 = new TH1D ("distChi2FitFcnPmt1", "", 100, 0, 10);
    auto distChi2FitFcnPmt2 = new TH1D ("distChi2FitFcnPmt2", "", 100, 0, 10);
    auto distChi2FitFcnPmt3 = new TH1D ("distChi2FitFcnPmt3", "", 100, 0, 10);
    //
    fillingDistribution(distChi2FitFcnPmt1, chi2FitFcn[0]);
    fillingDistribution(distChi2FitFcnPmt2, chi2FitFcn[1]);
    fillingDistribution(distChi2FitFcnPmt3, chi2FitFcn[2]);
    //
    // Filling for P-Values
    auto distPvalPmt1 = new TH1D ("distPvalPmt1", "", 100, -50, 0);
    auto distPvalPmt2 = new TH1D ("distPvalPmt2", "", 100, -50, 0);
    auto distPvalPmt3 = new TH1D ("distPvalPmt3", "", 100, -50, 0);
    //
    fillingDistribution(distPvalPmt1, logPval[0]);
    fillingDistribution(distPvalPmt2, logPval[1]);
    fillingDistribution(distPvalPmt3, logPval[2]);
    //
    auto canvasPval = canvasStyle("canvasPval");
    canvasPval->cd();
    auto pad01 = new TPad("pad01", "pad01", 0.01, 0.5, 0.32, 1.);
    pad01->Draw();
    pad01->cd();
    distPvalPmt1->SetLineColor(kRed);
    distPvalPmt1->GetXaxis()->SetTitle("Log10(Pval) [au]");
    distPvalPmt1->GetYaxis()->SetTitle("Counts [au]");
    distPvalPmt1->Draw();
    canvasPval->cd();
    auto pad02 = new TPad("pad02", "pad02", 0.33, 0.5, 0.65, 1.);
    pad02->Draw();
    pad02->cd();
    distPvalPmt2->SetLineColor(kOrange+1);
    distPvalPmt2->GetXaxis()->SetTitle("Log10(Pval) [au]");
    distPvalPmt2->GetYaxis()->SetTitle("Counts [au]");         
    distPvalPmt2->Draw();
    canvasPval->cd();
    auto pad03 = new TPad("pad03", "pad03", 0.66, 0.5, 0.99, 1.);
    pad03->Draw();
    pad03->cd();
    distPvalPmt3->SetLineColor(kBlue);
    distPvalPmt3->GetXaxis()->SetTitle("Log10(Pval) [au]");
    distPvalPmt3->GetYaxis()->SetTitle("Counts [au]");
    distPvalPmt3->Draw();
    //
    canvasPval->cd();
    auto pad11 = new TPad("pad11", "pad11", 0.01, 0., 0.32, 0.48);
    pad11->Draw();
    pad11->cd();
    distChi2FitFcnPmt1->SetLineColor(kRed);
    distChi2FitFcnPmt1->GetXaxis()->SetTitle("#chi^{2} Ndf [au]");
    distChi2FitFcnPmt1->GetYaxis()->SetTitle("Counts [au]");
    //distChi2FitFcnPmt1->GetXaxis()->SetRangeUser(0, 1);
    distChi2FitFcnPmt1->Draw();
    canvasPval->cd();
    auto pad12 = new TPad("pad12", "pad12", 0.33, 0., 0.65, 0.48);
    pad12->Draw();
    pad12->cd();
    distChi2FitFcnPmt2->SetLineColor(kOrange+1);
    distChi2FitFcnPmt2->GetXaxis()->SetTitle("#chi^{2} [au]");
    distChi2FitFcnPmt2->GetYaxis()->SetTitle("Counts [au]");    
    //distChi2FitFcnPmt2->GetXaxis()->SetRangeUser(0, 1);
    distChi2FitFcnPmt2->Draw();
    canvasPval->cd();
    auto pad13 = new TPad("pad13", "pad13", 0.66, 0., 0.99, 0.48);
    pad13->Draw();
    pad13->cd();
    distChi2FitFcnPmt3->SetLineColor(kBlue);
    distChi2FitFcnPmt3->GetXaxis()->SetTitle("#chi^{2} [au]");
    distChi2FitFcnPmt3->GetYaxis()->SetTitle("Counts [au]");    
    //distChi2FitFcnPmt3->GetXaxis()->SetRangeUser(0, 1);
    distChi2FitFcnPmt3->Draw();
    //
    canvasPval->cd();    
    canvasPval->Print("distPvalPmts.pdf");
    //
    // Doing weighed deltas
    vector < vector < double > > deltaWeighed(3);
    vector < vector < double > > deltaWeighedErr(3);
    vector < vector < double > > vhWeighed(3);
    vector < vector < double > > vhWeighedErr(3);
    //
    for(int pmt_i=0; pmt_i<3; pmt_i++) {
        for(int sd_i=0; sd_i<2000; sd_i++) {
            double numDelta = 0.;
            double denDelta = 0.;
            double numVh = 0.;
            double denVh = 0.;
            if(deltaPerSt[pmt_i][sd_i].size() < 1)
                continue;
            //
            for(int dlt_i=0; dlt_i<deltaPerSt[pmt_i][sd_i].size(); dlt_i++) {
                numDelta += deltaPerSt[pmt_i][sd_i][dlt_i] / pow(deltaPerStErr[pmt_i][sd_i][dlt_i], 2);
                denDelta += 1. / pow(deltaPerStErr[pmt_i][sd_i][dlt_i], 2);
                //
                numVh += vhPerSt[pmt_i][sd_i][dlt_i] / pow(vhPerStErr[pmt_i][sd_i][dlt_i], 2);
                denVh += 1. / pow(vhPerStErr[pmt_i][sd_i][dlt_i], 2);
            }
            deltaWeighed[pmt_i].push_back(numDelta / denDelta);
            deltaWeighedErr[pmt_i].push_back(sqrt(1. / denDelta));
            vhWeighed[pmt_i].push_back(numVh / denVh);
            vhWeighedErr[pmt_i].push_back(sqrt(1. / denVh));
        }
    }
    //
    // Doing distribution for weighed deltas    
    auto weighDeltaDistPmt1 = new TH1D ("weighDeltaDistPmt1", "", 400, -20, 20);
    auto weighDeltaDistPmt2 = new TH1D ("weighDeltaDistPmt2", "", 400, -20, 20);
    auto weighDeltaDistPmt3 = new TH1D ("weighDeltaDistPmt3", "", 400, -20, 20);
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
    auto deltaVsVhPmt1 = new TGraphErrors (deltaWeighed[0].size(), &vhWeighed[0].front(), 
        &deltaWeighed[0].front(), &vhWeighedErr[0].front(), &deltaWeighedErr[0].front());
    auto deltaVsVhPmt2 = new TGraphErrors (deltaWeighed[1].size(), &vhWeighed[1].front(),
        &deltaWeighed[1].front(), &vhWeighedErr[1].front(), &deltaWeighedErr[1].front());
    auto deltaVsVhPmt3 = new TGraphErrors (deltaWeighed[2].size(), &vhWeighed[2].front(),
        &deltaWeighed[2].front(), &vhWeighedErr[2].front(), &deltaWeighedErr[2].front());
    //
    auto canvasDltVh = canvasStyle("canvasDltVh");
    canvasDltVh->cd();
    pad01 = new TPad("pad01", "pad01", 0.01, 0., 0.32, 1.);
    pad01->Draw();
    pad01->cd();
    //
    deltaVsVhPmt1->SetName("deltaVsVhPmt1");
    deltaVsVhPmt1->SetTitle("PMT1");
    deltaVsVhPmt1->GetXaxis()->SetTitle("v/h [au]");
    deltaVsVhPmt1->GetYaxis()->SetTitle("#Delta [%]");
    deltaVsVhPmt1->SetLineColor(kRed);
    deltaVsVhPmt1->SetMarkerStyle(21);
    deltaVsVhPmt1->Draw("ap");
    deltaVsVhPmt1->Write();
    //
    canvasDltVh->cd();
    pad02 = new TPad("pad02", "pad02", 0.33, 0., 0.65, 1.);
    pad02->Draw();
    pad02->cd();
    deltaVsVhPmt2->SetName("deltaVsVhPmt2");
    deltaVsVhPmt2->SetTitle("PMT2");
    deltaVsVhPmt2->GetXaxis()->SetTitle("v/h [au]");
    deltaVsVhPmt2->GetYaxis()->SetTitle("#Delta [%]");
    deltaVsVhPmt2->SetLineColor(kOrange+1);
    deltaVsVhPmt2->SetMarkerStyle(21);
    deltaVsVhPmt2->Draw("ap");
    deltaVsVhPmt2->Write();
    //
    canvasDltVh->cd();
    pad03 = new TPad("pad03", "pad03", 0.66, 0., 1., 1.);
    pad03->Draw();
    pad03->cd();
    deltaVsVhPmt3->SetName("deltaVsVhPmt3");
    deltaVsVhPmt3->SetTitle("PMT3");
    deltaVsVhPmt3->GetXaxis()->SetTitle("v/h [au]");
    deltaVsVhPmt3->GetYaxis()->SetTitle("#Delta [%]");
    deltaVsVhPmt3->SetLineColor(kBlue);
    deltaVsVhPmt3->SetMarkerStyle(21);
    deltaVsVhPmt3->Draw("ap");
    deltaVsVhPmt3->Write();
    //
    canvasDltVh->Print("deltasVsVhPmt1.pdf");
    //
    // Writing and closing output root file
    cout << "MSD, writting and closing" << endl;
    outputRoot->Write();
    outputRoot->Close();
    //
    exit(1);
}
