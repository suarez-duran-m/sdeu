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

TH1D *doRelRmsDist( TString bname, 
    vector<double> listSt, bool ifIsUub ) {

  TString monthUub[2] = {"Aug", "Sep"};
  TString pmtId;
  TString strStId;
  TString strYear[3] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;

  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int nBins = (ifIsUub) ? 2500 : 250;
  int firstBin = (ifIsUub) ? 500 : 50;
  int lastBin = (ifIsUub) ? 3000 : 300;

  TString strTHname = (ifIsUub) ? "uub" : "ub";

  TH1D *qpkDist = new TH1D("qpkDist"+strTHname, "", nBins, firstBin, lastBin);
  TH1D *retRmsQpkDist = new TH1D("retQpkDist"+strTHname, "", 30, 0., 30.); // 0% to 20%
  vector < double > cumuRelRmsQpk;
  for (int i=0; i<3; i++ ) 
    cumuRelRmsQpk.push_back(0.);
  int nPmtWithQpk = 0;

  for ( int stId_i=0; stId_i<listSt.size(); stId_i++ ) {
    for ( int pmt=1; pmt<4; pmt++ ) {
      for ( int year=stYear; year<=lstYear; year++ ) {
        for ( int month_i=0; month_i<sizeof(monthUub)/sizeof(*monthUub); month_i++ ) {
          pmtId.Form("%d", pmt);
          strStId.Form("St%d", (int)listSt[stId_i]);
          fname = bname + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
  
          f = TFile::Open(fname+".root");
          chargeInfo = (TTree*)f->Get(strChargeData);
          chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);

          fetchQpkVals = 0.;
          for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
            chargeInfo->GetEntry(etry);
            if ( fetchQpkVals == 0 || fetchQpkVals < 0 )
              continue;
            qpkDist->Fill( fetchQpkVals );
          }
        }
        f->Clear();
        f->Close();
      }
      if ( qpkDist->GetMean() > 0 ) { 
        cumuRelRmsQpk[pmt-1] += 100.*qpkDist->GetRMS()/qpkDist->GetMean();
        nPmtWithQpk++;   
      }
      qpkDist->Reset();
    }
    double tmp = 0.;
    for ( int i=0; i<3; i++ ) {
      tmp += cumuRelRmsQpk[i];
      cumuRelRmsQpk[i] = 0.;
    }
    retRmsQpkDist->Fill( tmp/nPmtWithQpk );
    cout << ifIsUub << " " << listSt[stId_i] << " " << tmp << " " << tmp/nPmtWithQpk 
      << endl;
    cumuRelRmsQpk.clear();
    nPmtWithQpk = 0;
  }
  chargeInfo->Delete();
  delete f;
  return retRmsQpkDist;
}

TH1D *getQpksDist( TString bname, int StId, bool ifIsUub ) {

  TString monthUub[2] = {"Aug", "Sep"};
  TString pmtId;
  TString strStId;
  TString strYear[3] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;

  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int nBins = (ifIsUub) ? 2500 : 250;
  int firstBin = (ifIsUub) ? 500 : 50;
  int lastBin = (ifIsUub) ? 3000 : 300;

  TString strTHname = (ifIsUub) ? "uub" : "ub";

  TH1D *qpkDist = new TH1D("qpkDist"+strTHname, "", nBins, firstBin, lastBin);

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( int month_i=0; month_i<sizeof(monthUub)/sizeof(*monthUub); month_i++ ) {
      for ( int pmt=1; pmt<4; pmt++ ) {
        pmtId.Form("%d", pmt);
        strStId.Form("St%d", StId);
        fname = bname + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];

        f = TFile::Open(fname+".root");
        chargeInfo = (TTree*)f->Get(strChargeData);
        chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);
        
        fetchQpkVals = 0.;
        for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
          chargeInfo->GetEntry(etry);
          if ( fetchQpkVals == 0 || fetchQpkVals < 0 )
            continue;
          qpkDist->Fill( fetchQpkVals );
        }
      }
    }
    f->Clear(); 
    f->Close();
  }
  chargeInfo->Delete();
  delete f;
  return qpkDist;
}

void rmsOverMeanVsStId() {

  //string stListPath = "/home/msd/2021/sdeu/listStationsShort.txt";
  string stListPath = "/home/msd/2021/sdeu/fullUubStationsListVert.txt";
  TString bnOffl = "~/2021/sdeu/offline/forUb/results/offlineUb";
  TString bnUubOffl = "~/2021/sdeu/offline/forUub/results/offlineUub";
  TString bnCdas = "~/2021/sdeu/nouub/underHistos/results/ubChPkPMT";
  TString bnUubCdas = "~/2021/sdeu/underHisto/results/uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strEntr;
  TString strMean;
  TString strRms;
  bool ifIsUub = false;

  ifstream fileStList;
  double St_i;
  vector < double > stListId;
  vector< double > stList;
  fileStList.open(stListPath);
  int cnt = 1;
  while( fileStList.good() ) {
    fileStList >> St_i;
    stListId.push_back(St_i);
    stList.push_back(cnt);
    cnt++;
  }
  stListId.pop_back();
  stList.pop_back();
  fileStList.close();

  int nStations = stList.size();
  TH1D *distQpkStUb;
  TH1D *distQpkStUub;

  ifIsUub = false;
  distQpkStUb = doRelRmsDist(bnCdas, stListId, ifIsUub);

  ifIsUub = true;
  distQpkStUub = doRelRmsDist(bnUubCdas, stListId, ifIsUub);

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  distQpkStUub->SetStats(kFALSE);
  distQpkStUub->SetLineColor(kRed);
  distQpkStUub->SetLineWidth(2);
  distQpkStUub->GetYaxis()->SetTitle("Counts [au]");
  distQpkStUub->GetXaxis()->SetTitle("#LT RMS/#LT Q^{pk}_{VEM} #GT #GT [%]");
  distQpkStUub->Draw();
  distQpkStUb->SetLineColor(kBlack);
  distQpkStUb->SetLineWidth(2);
  distQpkStUb->Draw("same");

  leg = new TLegend(0.7,0.65,0.9,0.95);
  strMean.Form("%.3f", distQpkStUb->GetMean());
  strRms.Form("%.3f", distQpkStUb->GetRMS());
  leg->AddEntry(distQpkStUb, "UB", "l" );
  leg->AddEntry(distQpkStUb, "MEAN: "+strMean, "");
  leg->AddEntry(distQpkStUb, "RMS: "+strRms, "");
  strMean.Form("%.3f", distQpkStUub->GetMean());
  strRms.Form("%.3f", distQpkStUub->GetRMS());
  leg->AddEntry(distQpkStUub, "UUB", "l" );
  leg->AddEntry(distQpkStUub, "MEAN: "+strMean, "");
  leg->AddEntry(distQpkStUub, "RMS: "+strRms, "");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print("../plots/accuracyQpksFitsUbUubAllStAllPmt.pdf");

  exit(0);
} 
