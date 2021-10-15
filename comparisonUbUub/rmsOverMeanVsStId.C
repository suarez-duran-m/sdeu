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

void histoStyle(TH1F *hist) {
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}

void histoStyle(TGraphErrors *hist) {
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}

TH1D *getRelQpk( TH1D *qpkVals, vector<double> &qpkAllPmts ) {
  double mean = qpkVals->GetMean();
  TString strMean;
  strMean.Form("%.4f", mean);
  int nbins = 2000;
  double frtBin = 0.;
  double lstBin = 2.;
  double qpk = 0.;
  double relQpk = 0.;
  double relQpk2Bin = 1000.;
  int bin = 0;
  int nCnts = 0;
  TH1D *retRelQpk = new TH1D ("retRelQpk"+strMean, "", nbins, frtBin, lstBin);
  for ( int qpk_i=1; qpk_i<qpkVals->GetNbinsX(); qpk_i++ )
    if ( qpkVals->GetBinContent(qpk_i) > 0 ) {
      qpk = qpkVals->GetBinLowEdge(qpk_i);
      relQpk = qpk/mean;
      bin = relQpk2Bin*(relQpk+5e-4);
      nCnts = qpkVals->GetBinContent(qpk_i);
      qpkAllPmts[bin] += nCnts;
      nCnts += retRelQpk->GetBinContent( bin-1 );
      retRelQpk->SetBinContent( bin-1, nCnts );
    }
  return retRelQpk;
}

void getRelQpk( TH1D *qpkVals, TH1D *retRelQpk ) {
  double qpk = 0.;
  double relQpk = 0.;
  int nCnts = 0;
  for ( int qpk_i=1; qpk_i<qpkVals->GetNbinsX(); qpk_i++ )
    if ( qpkVals->GetBinContent(qpk_i) > 0 ) {
      qpk = qpkVals->GetBinLowEdge(qpk_i);
      nCnts = qpkVals->GetBinContent(qpk_i);
      nCnts += retRelQpk->GetBinContent( qpk_i-1 );
      retRelQpk->SetBinContent( qpk_i-1, nCnts );
    }
}

double getmean( vector<double> arr ) {
  int nb = arr.size();
  int goodVals = 0;
  double mean = 0.;
  for (int bin_i=0; bin_i<nb; bin_i++)
    if ( arr[bin_i] > 0 ) {
      mean += arr[bin_i];
      goodVals++;
    }
  return mean/goodVals;
}

double getrms( vector<double> arr, double meanarr ) {
  int nb = arr.size();
  double rms = 0.;
  int goodVals = 0;
  for (int bin_i=0; bin_i<nb; bin_i++)
    if ( arr[bin_i]> 0 ) {
      rms += (arr[bin_i]-meanarr)*(arr[bin_i]-meanarr);
      goodVals++;
    }

  return sqrt(rms/goodVals);
}

void rmsOverMeanVsStId(bool ifIsUub) {

  TPaveStats *ptstats;
  TLegend *leg;
  TString strEntr;
  TString strMean;
  TString strRms;
 
  TString typeSt = (ifIsUub) ? "StationsUub" : "StationsUb";
  TString strFname = (ifIsUub) ? "qpkValStationsUub.root" : "qpkValStationsUb.root"; 
  TFile *f = TFile::Open(strFname);
  TTree *stIds = (TTree*)f->Get(typeSt);
  TH1D *fetchQpkValsPmt1 = new TH1D();
  TH1D *fetchQpkValsPmt2 = new TH1D();
  TH1D *fetchQpkValsPmt3 = new TH1D();
  
  stIds->SetBranchAddress("qpkValuesPmt1", &fetchQpkValsPmt1);
  stIds->SetBranchAddress("qpkValuesPmt2", &fetchQpkValsPmt2);
  stIds->SetBranchAddress("qpkValuesPmt3", &fetchQpkValsPmt3);

  int nbins = 2000;
  double frtBin = 0.;
  double lstBin = 2.;
  TH1D *distQpkAveQpkPmt1Uub = new TH1D(); 
  TH1D *distQpkAveQpkPmt2Uub = new TH1D(); 
  TH1D *distQpkAveQpkPmt3Uub = new TH1D(); 
  vector < double > qpkAveQpkUubPmts;
  for ( int i=0; i<nbins; i++ )
    qpkAveQpkUubPmts.push_back(0.);
  TH1D *distQpkAveQpkUubPmts = new TH1D("distQpkAveQpkUubPmts", "", 
      nbins, frtBin, lstBin);
  TH1D *distQpkAveQpkUub = new TH1D("distQpkAveQpkUub", "",
      nbins, frtBin, lstBin);

  TH1D *distRelRMSqpkUub = new TH1D("distRelRMSqpkUub", "", 20, 0., 20.);

  double rmsQpkRelAllPmt = 0.;
  double rmsQpkRelPmt = 0.;
  double farFromSigma = 1.1;

  for ( int etry=0; etry<stIds->GetEntries(); etry++ ) {
    stIds->GetEntry( etry );
    for ( int i=0; i<nbins; i++ )
      qpkAveQpkUubPmts.push_back(0.);

    distQpkAveQpkPmt1Uub = getRelQpk( fetchQpkValsPmt1, qpkAveQpkUubPmts );
    distQpkAveQpkPmt2Uub = getRelQpk( fetchQpkValsPmt2, qpkAveQpkUubPmts );
    distQpkAveQpkPmt3Uub = getRelQpk( fetchQpkValsPmt3, qpkAveQpkUubPmts );

    for ( int bin=1; bin<qpkAveQpkUubPmts.size(); bin++ )
      if ( qpkAveQpkUubPmts[bin] > 0 ) {
        distQpkAveQpkUubPmts->SetBinContent( bin-1, qpkAveQpkUubPmts[bin] );
      }

    rmsQpkRelAllPmt = farFromSigma*distQpkAveQpkUubPmts->GetRMS();
    rmsQpkRelPmt = distQpkAveQpkPmt1Uub->GetRMS();
    if ( rmsQpkRelPmt < rmsQpkRelAllPmt )
      getRelQpk(distQpkAveQpkPmt1Uub, distQpkAveQpkUub);
    rmsQpkRelPmt = distQpkAveQpkPmt2Uub->GetRMS();
    if ( rmsQpkRelPmt < rmsQpkRelAllPmt )
      getRelQpk(distQpkAveQpkPmt2Uub, distQpkAveQpkUub);
    rmsQpkRelPmt = distQpkAveQpkPmt3Uub->GetRMS();
    if ( rmsQpkRelPmt < rmsQpkRelAllPmt )
      getRelQpk(distQpkAveQpkPmt3Uub, distQpkAveQpkUub);

    rmsQpkRelAllPmt = 100.*distQpkAveQpkUub->GetRMS()/distQpkAveQpkUub->GetMean();
    distRelRMSqpkUub->Fill( rmsQpkRelAllPmt );

    distQpkAveQpkPmt1Uub->Reset();
    distQpkAveQpkPmt2Uub->Reset();
    distQpkAveQpkPmt3Uub->Reset();
    distQpkAveQpkUubPmts->Reset();
    distQpkAveQpkUub->Reset();

    qpkAveQpkUubPmts.clear();
  }

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  distRelRMSqpkUub->SetStats(kFALSE);
  distRelRMSqpkUub->SetLineColor(kRed);
  distRelRMSqpkUub->SetLineWidth(2);
  distRelRMSqpkUub->GetYaxis()->SetTitle("Counts [au]");
  distRelRMSqpkUub->GetXaxis()->SetTitle("#LT RMS/#LT Q^{pk}_{VEM} #GT #GT [%]");
  distRelRMSqpkUub->Draw();
/*
  distRelRMSqpkUb->SetLineColor(kBlack);
  distRelRMSqpkUb->SetLineWidth(2);
  distRelRMSqpkUb->Draw("same");  
*/
  leg = new TLegend(0.7,0.65,0.9,0.95);
  /*
  strMean.Form("%.3f", distRelRMSqpkUb->GetMean());
  strRms.Form("%.3f", distRelRMSqpkUb->GetRMS());
  leg->AddEntry(distRelRMSqpkUb, "UB", "l" );
  leg->AddEntry(distRelRMSqpkUb, "MEAN: "+strMean, "");
  leg->AddEntry(distRelRMSqpkUb, "RMS: "+strRms, "");
  */
  strMean.Form("%.3f", distRelRMSqpkUub->GetMean());
  strRms.Form("%.3f", distRelRMSqpkUub->GetRMS());
  leg->AddEntry(distRelRMSqpkUub, "UUB", "l" );
  leg->AddEntry(distRelRMSqpkUub, "MEAN: "+strMean, "");
  leg->AddEntry(distRelRMSqpkUub, "RMS: "+strRms, "");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print("kk.pdf");
  //c1->Print("../plots/accuracyQpksFitsUbUubAllStAllPmt.pdf");

  //exit(0);
} 
