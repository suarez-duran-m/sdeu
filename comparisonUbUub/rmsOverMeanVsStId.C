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

  return sqrt(rms/goodVals);
}

vector < double > getQpkValues( TString bname, double StId, int pmt, 
    bool ifIsUub ) {

  TString monthUub[2] = {"Aug", "Sep"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString strStId;
  strStId.Form("St%d", (int)StId);
  TString strYear[3] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;
  TString strTHname = (ifIsUub) ? "uub" : "ub";

  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  vector < double > retQpkDist;

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( int month_i=0; month_i<sizeof(monthUub)/sizeof(*monthUub); month_i++ ) {
      fname = bname + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);
      
      fetchQpkVals = 0.;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
        chargeInfo->GetEntry(etry);
        if ( fetchQpkVals == 0 || fetchQpkVals < 0 )
          continue;
        retQpkDist.push_back( fetchQpkVals );
      }
      f->Clear();
      f->Close();
    }
  }
  chargeInfo->Delete();
  delete f;
  return retQpkDist;
}

void rmsOverMeanVsStId() {

  //string stListPath = "/home/msd/2021/sdeu/listStationsShort.txt";
  string stListPath = "/home/msd/2021/sdeu/fullUubStationsListVert.txt";
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
  fileStList.open(stListPath);
  while( fileStList.good() ) {
    fileStList >> St_i;
    stListId.push_back(St_i);
  }
  stListId.pop_back();
  fileStList.close();

  vector < vector < double > > distQpkUb;
  vector < vector < double > > distQpkUub;
  distQpkUb.resize(3);
  distQpkUub.resize(3);
  vector < double > aveQpkUb;
  vector < double > aveQpkUub;
  vector < double > rmsQpkUub;
  TH1D *distQpk_AveQpk3PmtUb = new TH1D("distQpkOverAveQpkUb", "", 20, 0., 2.);
  TH1D *distQpk_AveQpk3PmtUub = new TH1D("distQpkOverAveQpkUub", "", 2000, 0., 2.);
  TH1D *distRelRMSqpkUb = new TH1D("distRelRMSqpkUb", "", 20, 0., 20.); 
  TH1D *distRelRMSqpkUub = new TH1D("distRelRMSqpkUub", "", 20, 0., 20.);
  double fecthMean = 0.;
  TCanvas *c0 = canvasStyle("c0");
  c0->cd();
  TString strStId;

  for ( auto & st_i : stListId ) {
    aveQpkUb.resize(3);
    aveQpkUub.resize(3);
    rmsQpkUub.resize(3);

    for ( int pmt=1; pmt<4; pmt++ ) {
      distQpkUb[pmt-1] = getQpkValues(bnCdas, st_i, pmt, false);
      distQpkUub[pmt-1] = getQpkValues(bnUubCdas, st_i, pmt, true);
      fecthMean = getmean(distQpkUb[pmt-1]);
      aveQpkUb[pmt-1] = ( fecthMean > 100 && fecthMean < 200 ) ? fecthMean : 0.;
      fecthMean = getmean(distQpkUub[pmt-1]);  
      aveQpkUub[pmt-1] = ( fecthMean > 1000 && fecthMean < 2000 ) ? fecthMean : 0.;
    }
    for ( int pmt=0; pmt<3; pmt++ ) {
      if ( aveQpkUb[pmt] > 0 )
        for ( auto & qpk_i : distQpkUb[pmt] )
          if ( qpk_i > 100 && qpk_i < 200 )
            distQpk_AveQpk3PmtUb->Fill( qpk_i/aveQpkUb[pmt] );
      if ( aveQpkUub[pmt] > 0 )
        for ( auto & qpk_i : distQpkUub[pmt] )
          if ( qpk_i > 1000 && qpk_i < 2000 )
            distQpk_AveQpk3PmtUub->Fill( qpk_i/aveQpkUub[pmt]);
      distQpkUub[pmt].clear();
      aveQpkUub.clear();
    }
    double tmpRel = 0.;
    tmpRel = 100.*distQpk_AveQpk3PmtUb->GetRMS()/distQpk_AveQpk3PmtUb->GetMean(); 
    distRelRMSqpkUb->Fill( tmpRel );

    strStId.Form("%d",(int)st_i);
    distQpk_AveQpk3PmtUb->GetYaxis()->SetTitle("Counts [au]");
    distQpk_AveQpk3PmtUb->GetXaxis()->SetTitle("Q^{pk}_{VEM}/#LT Q^{pk}_{VEM} #GT_{pmt}");
    distQpk_AveQpk3PmtUb->GetXaxis()->SetTitleOffset(1.3);
    distQpk_AveQpk3PmtUb->Draw();
    c0->Print("distQpk_AveQpk3PmtUb"+strStId+".pdf");
    c0->Clear();
    distQpk_AveQpk3PmtUub->GetYaxis()->SetTitle("Counts [au]");
    distQpk_AveQpk3PmtUub->GetXaxis()->SetTitle("Q^{pk}_{VEM}/#LT Q^{pk}_{VEM} #GT_{pmt}");
    distQpk_AveQpk3PmtUub->GetXaxis()->SetTitleOffset(1.3);
    distQpk_AveQpk3PmtUub->Draw();
    c0->Print("distQpk_AveQpk3PmtUub"+strStId+".pdf");
  
    //if ( tmpRel > 6 )
      cout << "UB " << st_i << " " 
        << 100.*distQpk_AveQpk3PmtUb->GetRMS()/distQpk_AveQpk3PmtUb->GetMean()
        << endl;
    
    tmpRel = 100.*distQpk_AveQpk3PmtUub->GetRMS()/distQpk_AveQpk3PmtUub->GetMean();
    distRelRMSqpkUub->Fill( tmpRel );

    //if ( tmpRel > 4 && tmpRel < 6 )
      cout << "UUB " << st_i << " " 
        << 100.*distQpk_AveQpk3PmtUub->GetRMS()/distQpk_AveQpk3PmtUub->GetMean()
        << endl;
    distQpk_AveQpk3PmtUb->Reset();
    distQpk_AveQpk3PmtUub->Reset();
  }
 
  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  distRelRMSqpkUub->SetStats(kFALSE);
  distRelRMSqpkUub->SetLineColor(kRed);
  distRelRMSqpkUub->SetLineWidth(2);
  distRelRMSqpkUub->GetYaxis()->SetTitle("Counts [au]");
  distRelRMSqpkUub->GetXaxis()->SetTitle("#LT RMS/#LT Q^{pk}_{VEM} #GT #GT [%]");
  distRelRMSqpkUub->Draw();

  distRelRMSqpkUb->SetLineColor(kBlack);
  distRelRMSqpkUb->SetLineWidth(2);
  distRelRMSqpkUb->Draw("same");  

  leg = new TLegend(0.7,0.65,0.9,0.95);
  strMean.Form("%.3f", distRelRMSqpkUb->GetMean());
  strRms.Form("%.3f", distRelRMSqpkUb->GetRMS());
  leg->AddEntry(distRelRMSqpkUb, "UB", "l" );
  leg->AddEntry(distRelRMSqpkUb, "MEAN: "+strMean, "");
  leg->AddEntry(distRelRMSqpkUb, "RMS: "+strRms, "");
  strMean.Form("%.3f", distRelRMSqpkUub->GetMean());
  strRms.Form("%.3f", distRelRMSqpkUub->GetRMS());
  leg->AddEntry(distRelRMSqpkUub, "UUB", "l" );
  leg->AddEntry(distRelRMSqpkUub, "MEAN: "+strMean, "");
  leg->AddEntry(distRelRMSqpkUub, "RMS: "+strRms, "");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print("../plots/accuracyQpksFitsUbUubAllStAllPmt.pdf");
  
  exit(0);
} 
