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

  return sqrt(rms/goodVals)/meanarr;
}

TH1D * getQpksDist( TString bname, int StId, int pmt, bool ifIsUub ) {

  TString monthUub[2] = {"Aug", "Sep"};
  TString pmtId;
  TString strStId;
  TString strYear[3] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;

  TString strChargeData = "ChargeData";
  pmtId.Form("%d", pmt);

  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int fetchEvId = 0;
  int fetchTime = 0;

  int nBins = (ifIsUub) ? 2500 : 250;
  int firstBin = (ifIsUub) ? 500 : 50;
  int lastBin = (ifIsUub) ? 3000 : 300;

  TString strTHname = (ifIsUub) ? "uub" : "ub";
  TH1D *qpkDist = new TH1D("qpkDist"+strTHname, "", nBins, firstBin, lastBin);

  for ( int year=stYear; year<=lstYear; year++ )
    for ( int month_i=0; month_i<sizeof(monthUub)/sizeof(*monthUub); month_i++ ) {
      for ( int pmt_i=1; pmt_i<4; pmt_i++ ) {
      strStId.Form("St%d", StId);
      pmtId.Form("%d", pmt_i);
      fname = bname + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
     
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);
      chargeInfo->SetBranchAddress("eventId", &fetchEvId);
      chargeInfo->SetBranchAddress("timeEvnt", &fetchTime);
      
      fetchQpkVals = 0.;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
        chargeInfo->GetEntry(etry);
        if ( fetchQpkVals == 0 || fetchQpkVals < 0 )
          continue;
        qpkDist->Fill( fetchQpkVals );
      }
    }
    }
  chargeInfo->Delete();
  f->Delete();
  return qpkDist;
}

void rmsOverMeanPerStIdPmt(int st, double rmsQpk, bool ifIsUub=true) {

  TString bnCdas = (ifIsUub) ? 
    "~/2021/sdeu/underHisto/results/uubChPkPMT" :
    "~/2021/sdeu/nouub/underHistos/results/ubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strSt;
  strSt.Form("%d", st);
  TString struub;

  TH1D *distQpkStPmt1;
  TH1D *distQpkStPmt2;
  TH1D *distQpkStPmt3;

  distQpkStPmt1 = getQpksDist(bnCdas, st, 1, ifIsUub);
  //distQpkStPmt2 = getQpksDist(bnCdas, st, 2, ifIsUub);
  //distQpkStPmt3 = getQpksDist(bnCdas, st, 3, ifIsUub);

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  
  distQpkStPmt1->SetStats(kFALSE);
  distQpkStPmt1->SetTitle("");
  distQpkStPmt1->GetYaxis()->SetTitle("Counts [au]");
  distQpkStPmt1->GetXaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  distQpkStPmt1->GetXaxis()->SetTitleOffset(1.5);
  distQpkStPmt1->SetLineColor(kGreen+2);//kRed);
  distQpkStPmt1->Draw();

  //distQpkStPmt2->SetLineColor(kBlack);
  //distQpkStPmt2->Draw("same");

  //distQpkStPmt3->SetLineColor(kGreen+2);
  //distQpkStPmt3->Draw("same");

  leg = new TLegend(0.6,0.5,0.9,0.95);
  TString strAve;
  strAve.Form("%.1f", rmsQpk);
  struub = (ifIsUub) ? "UUB" : "UB";
  leg->SetHeader(struub+" Station "+strSt+"; #LT RMS/#LT Q^{pk}_{VEM} #GT #GT: "+strAve+" [%]");
  //leg->AddEntry(distQpkStPmt1, "PMT1", "l");
  
  strAve.Form("%.3f", distQpkStPmt1->GetMean());
  leg->AddEntry(distQpkStPmt1, "MEAN: "+strAve+" [FADC]", "");
  strAve.Form("%.3f", distQpkStPmt1->GetRMS());
  leg->AddEntry(distQpkStPmt1, "RMS: "+strAve+" [FADC]", "");
/*
  leg->AddEntry(distQpkStPmt2, "PMT2", "l");
  strAve.Form("%.3f", distQpkStPmt2->GetMean());
  leg->AddEntry(distQpkStPmt2, "MEAN: "+strAve+" [FADC]", "");
  strAve.Form("%.3f", distQpkStPmt2->GetRMS());
  leg->AddEntry(distQpkStPmt2, "RMS: "+strAve+" [FADC]", "");

  leg->AddEntry(distQpkStPmt3, "PMT3", "l");
  strAve.Form("%.3f", distQpkStPmt3->GetMean());
  leg->AddEntry(distQpkStPmt3, "MEAN: "+strAve+" [FADC]", "");
  strAve.Form("%.3f", distQpkStPmt3->GetRMS());
  leg->AddEntry(distQpkStPmt3, "RMS: "+strAve+" [FADC]", "");
*/
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  struub = (ifIsUub) ? "Uub" : "Ub";  
  //c1->Print("../plots/qpksFitsCdas"+struub+"St"+strSt+".pdf");

  //exit(0);
} 
