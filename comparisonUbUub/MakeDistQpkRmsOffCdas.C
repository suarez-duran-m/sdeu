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
  double mean = 0.;
  int nb = arr.size();
  int ngoodb = 0;
    for (auto & elem : arr)
      if ( elem > 0 ) { 
        mean += elem;
        ngoodb++;
      }
  return mean/ngoodb;
}

double getrms( vector<double> arr, double meanarr ) {
  double rms = 0.;
  int nb = arr.size();
  int ngoodb = 0;
  for (auto & elem : arr)
    if ( elem > 0 ) {
      rms += (elem - meanarr)*(elem - meanarr);
      ngoodb++;
    }
  return sqrt(rms/ngoodb)/meanarr;
}

TH1F *getFitVals( TString bname, TString st, int pmt, bool ifIsUub, bool ifoff ) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strYear[] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;
  // Sept. 1st. for 2019, 2020, and 2021
  int beforeSept [3] = {1251331218, 1282953618, 1314489618};
  TString strChargeData = (ifoff) ? "charge" : "ChargeData";
  pmtId.Form("%d", pmt);

  TFile *f;
  TTree *chargeInfo;
  double getQpkVals = 0.;
  int nBins = (ifIsUub) ? 2000 : 200;
  double stBin = 0.; // (ifIsUub) ? 1000. : 100.;
  double endBin = (ifIsUub) ? 2000. : 200.; 
  TH1F *retQpkDist = new TH1F("retQpkDist", "", nBins, stBin, endBin);
  int nMonths = 0;

  int evttime = 0;
  int prevTime = 0;
  int EvtId = 0;
  //vector < double > tmpQpk;

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( int month=nMonths; month<nMonths+1; month++ ) {
      fname = (ifoff) ? 
        (bname + monthUub[month] + strYear[year] + st + "Pmt" + pmtId)
        : (bname + pmtId + st + "lrb35" + monthUub[month] + strYear[year]);
 
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &getQpkVals); 
      if ( ifoff ) {
        chargeInfo->SetBranchAddress("GpsTime", &evttime);
        chargeInfo->SetBranchAddress("evtId", &EvtId);
      }
      else {
        chargeInfo->SetBranchAddress("timeEvnt", &evttime);
        chargeInfo->SetBranchAddress("eventId", &EvtId);
      }

      prevTime = 0;
      getQpkVals = 0.;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {

        chargeInfo->GetEntry(etry);
        if ( evttime < beforeSept[year] ) { // Events before Sept. 1st, 2021
          // To avoid read the same event twice from Offline
          if ( ifoff ) {
            if ( prevTime != evttime ) {
              retQpkDist->Fill( getQpkVals );
              prevTime = evttime;
              //if ( ifIsUub )
                //tmpQpk.push_back( getQpkVals );
            }
          }
          else
            retQpkDist->Fill( getQpkVals );
        }
      }
    }
  }
/*
  double tmpAve2 = 0.;
  double tmpRelRms = 0.;
  if ( ifIsUub && ifoff ) {
    tmpAve2 = getmean( tmpQpk );
    tmpRelRms = getrms( tmpQpk, tmpAve2 );
    cout << "tmpAve2 " << tmpAve2 << " relRms " << tmpRelRms << endl;
  }
  */
  chargeInfo->Delete();
  f->Delete();
  return retQpkDist;
  retQpkDist->Delete();
}

vector < double > getFailsVals( TString bname, TString st, int pmt, int year, bool ifoff) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  // Still without program properly.
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strSt;
  TString fname;
  TString strYear;
  TString strChargeData;
  pmtId.Form("%d", pmt);
  strYear.Form("%d", year);

  TFile *f;
  TTree *chargeInfo;
  double getQpk = 0.;
  vector < double > retQpk;
  int nMonths = 0;

  int evttime = 0;
  int prevTime = 0;

  if ( !ifoff )
    strChargeData = "ChargeData";
  else
    strChargeData = "charge";

  for ( int month=nMonths; month<nMonths+1; month++ ) {
    if ( !ifoff )
      fname = bname + pmtId + st + "lrb35" + monthUub[month] + strYear;
    else if ( ifoff )
      fname = bname + monthUub[month] + strYear + st + "Pmt" + pmtId;

    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);

    getQpk = 0.;
    chargeInfo->SetBranchAddress("chargeVal", &getQpk);
    if ( ifoff && year==2021 )
      chargeInfo->SetBranchAddress("GpsTime", &evttime);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
      chargeInfo->GetEntry(etry);

      if ( evttime <= 0 || int(getQpk) <= 0 )  
        continue;

      if ( ifoff && year==2021 ) {
        if ( prevTime != evttime ) {
          if ( getQpk == 0 )
            retQpk.push_back( getQpk );
          prevTime = evttime;
        }
      }
      else {
       if ( getQpk == 0 )
         retQpk.push_back( getQpk );
      }
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return retQpk;
}

TGraph *getGraph( vector<int> time, vector<double> values) {
  int nPoints = time.size();
  double xtime[ nPoints ];
  double yvalues[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtime[i] = time[i];
    yvalues[i] = values[i];
  }

  TGraph *grp = new TGraph(nPoints,xtime,yvalues);
  return grp;
}


// =================================
// *** *** *** MAIN CODE *** *** ***

void MakeDistQpkRmsOffCdas (int st, int pmt) {

  TString statId;
  statId.Form("St%d", st);

  TString bnOffl = "~/2021/sdeu/offline/forUb/Aug/offlineUb";
  TString bnUubOffl = "~/2021/sdeu/offline/forUub/AugResults/offlineUub";
  TString bnCdas = "~/2021/sdeu/nouub/underHistos/AugResults/ubChPkPMT";
  TString bnUubCdas = "~/2021/sdeu/underHisto/AugResults/uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strMean;
  TString strRms;
  TString strFails;
  TString strGoods;
  TString strSt;
  TString strPmt;
  strSt.Form("%d", st);
  strPmt.Form("%d", pmt);
  bool ifoff = false;
  bool ifIsUub = false;
  TString whinfo = "chargeVal";
  double uub2ub = 11.4;

  TH1F *qpkUbOff;  
  TH1F *qpkUubOff; 
  TH1F *qpkUbCdas; 
  TH1F *qpkUubCdas;

  ifoff = true;
  ifIsUub = false; 
  qpkUbOff = getFitVals(bnOffl, statId, pmt, ifIsUub, ifoff);
  ifIsUub = true;
  qpkUubOff = getFitVals(bnUubOffl, statId, pmt, ifIsUub, ifoff);

  ifoff = false;
  ifIsUub = false;
  qpkUbCdas = getFitVals(bnCdas, statId, pmt, ifIsUub, ifoff);
  ifIsUub = true;
  qpkUubCdas = getFitVals(bnUubCdas, statId, pmt, ifIsUub, ifoff);

  // Doing plots
  TCanvas *c1 = canvasStyle("c1");
  c1->cd(0);
  qpkUbOff->SetStats(kFALSE);
  qpkUbOff->SetLineColor(kGreen+3);
  qpkUbOff->SetLineWidth(2);
  qpkUbOff->SetFillStyle(3001);
  qpkUbOff->SetFillColor(kGreen+3);
  qpkUbOff->GetXaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  qpkUbOff->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(qpkUbOff);
  qpkUbOff->Draw();
  qpkUbCdas->SetLineColor(kMagenta+2);
  qpkUbCdas->SetLineWidth(2);
  qpkUbCdas->SetFillStyle(3001);
  qpkUbCdas->SetFillColor(kMagenta+2);
  qpkUbCdas->Draw("same");

  leg = new TLegend(0.55, 0.6, 0.9, 0.95);
  leg->SetHeader("UB: St "+strSt+", PMT"+strPmt);
  leg->AddEntry(qpkUbOff, "From SdCalibrator");
  strMean.Form("%.2f", qpkUbOff->GetMean());
  strRms.Form("%.4f", qpkUbOff->GetRMS()/qpkUbOff->GetMean());
  leg->AddEntry(qpkUbOff, "#mu="+strMean+"; RMS/#mu="+strRms,"");
  leg->AddEntry(qpkUbOff, "", "");
  leg->AddEntry(qpkUbCdas, "Local Implementation");
  strMean.Form("%.2f", qpkUbCdas->GetMean());
  strRms.Form("%.4f", qpkUbCdas->GetRMS()/qpkUbCdas->GetMean()); 
  leg->AddEntry(qpkUbCdas, "#mu="+strMean+"; RMS/#mu="+strRms,"");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("qpkDistributionUb.pdf");

  TCanvas *c2 = canvasStyle("c2");
  c2->cd(0);
  qpkUubOff->SetStats(kFALSE);
  qpkUubOff->SetLineColor(kGreen+3);
  qpkUubOff->SetLineWidth(2);
  qpkUubOff->SetFillStyle(3001);
  qpkUubOff->SetFillColor(kGreen+3);
  qpkUubOff->GetXaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  qpkUubOff->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(qpkUubOff);
  qpkUubOff->Draw();
  qpkUubCdas->SetLineColor(kMagenta+2);
  qpkUubCdas->SetLineWidth(2);
  qpkUubCdas->SetFillStyle(3001);
  qpkUubCdas->SetFillColor(kMagenta+2);
  qpkUubCdas->Draw("same");

  leg = new TLegend(0.55, 0.6, 0.9, 0.95);
  leg->SetHeader("UUB: St "+strSt+", PMT"+strPmt);
  leg->AddEntry(qpkUubOff, "From SdCalibrator");
  strMean.Form("%.2f", qpkUubOff->GetMean());
  strRms.Form("%.4f", qpkUubOff->GetRMS()/qpkUubOff->GetMean());
  leg->AddEntry(qpkUubOff, "#mu="+strMean+"; RMS/#mu="+strRms,"");
  leg->AddEntry(qpkUubOff, "", "");
  leg->AddEntry(qpkUubCdas, "Local Implementation");
  strMean.Form("%.2f", qpkUubCdas->GetMean());
  strRms.Form("%.4f", qpkUubCdas->GetRMS()/qpkUubCdas->GetMean());
  leg->AddEntry(qpkUubCdas, "#mu="+strMean+"; RMS/#mu="+strRms,"");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c2->Print("qpkDistributionUub.pdf");
  
  //exit(0);
}
