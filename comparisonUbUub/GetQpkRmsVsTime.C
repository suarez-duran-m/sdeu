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

void tgraphErrorStyle(TGraphErrors *hist) {
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}

void tgraphStyle(TGraph *hist) {
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
  return sqrt(rms/ngoodb);
}

double getRelRms( vector<double> arr, double meanarr ) {
  double rms = 0.;
  int nb = arr.size();
  int ngoodb = 0;
  for (auto & elem : arr)
    if ( elem > 0 ) {
      rms += (elem - meanarr)*(elem - meanarr);
      ngoodb++;
    }
  //cout << "inRelRms rms " << sqrt(rms/ngoodb) << " mean " << meanarr << endl;
  return sqrt(rms/ngoodb)/meanarr;
}

vector < vector < double > > getFitQpkTime( TString bname, TString st, int pmt, int year, bool ifoff ) {
  TString monthUub[] = {"Aug", "Sep"};
  TString pmtId;
  TString strYear;
  TString fname;
  TString strChargeData = (ifoff) ? "charge" : "ChargeData";
  pmtId.Form("%d", pmt);
  strYear.Form("%d", year);

  TFile *f;
  TTree *chargeInfo;
  double fetchQpk = 0.;
  vector < vector < double > > retQpkTime(61); // Number of days for August+Sept.
  int nMonths = 0;

  int evttime = 0;
  int prevTime = 0;
  int EvtId = 0;
  int crrDay = 0;

  int aug1st = 0;
  double uub2ub = (year==2021) ? 11.4 : 1.;
  if ( year == 2019 )
    aug1st = 1248652818;
  else if ( year == 2020 )
    aug1st = 1280275218;
  else 
    aug1st = 1311811218;

  //vector < double > tmpQpk;

  for ( int month=nMonths; month<nMonths+1; month++ ) {
    fname = (ifoff) ? (fname = bname + monthUub[month] + strYear + st + "Pmt" + pmtId)
      : (bname + pmtId + st + "lrb35" + monthUub[month] + strYear);

    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);

    chargeInfo->SetBranchAddress("chargeVal", &fetchQpk);
    // To avoid read the same event twice from Offline
    if ( ifoff ) {
      chargeInfo->SetBranchAddress("GpsTime", &evttime);
      chargeInfo->SetBranchAddress("evtId", &EvtId);
    }
    else {
      chargeInfo->SetBranchAddress("timeEvnt", &evttime);
      chargeInfo->SetBranchAddress("eventId", &EvtId);
    }

    prevTime = 0;
    fetchQpk = 0.;
    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {  
      chargeInfo->GetEntry(etry);
      if ( fetchQpk/uub2ub < 100 )
        continue;
      
      // To avoid read the same event twice from Offline
      crrDay = int((evttime-aug1st)/86400);
      if ( crrDay < 61 ) {
        if ( ifoff ) {
          if ( prevTime != evttime ) {
            retQpkTime[crrDay].push_back( fetchQpk/uub2ub );
            prevTime = evttime;
            //tmpQpk.push_back(fetchQpk);
          }
        }
        else
          retQpkTime[crrDay].push_back( fetchQpk/uub2ub );
      }
    }
  }
/*
  double tmpAve2 = getmean( tmpQpk );
  cout << "From GetVals" << endl;
  double tmpRms = getRelRms( tmpQpk, tmpAve2 );
  cout << "For all tmpAve2 " << tmpAve2 << " tmpRms " << tmpRms << endl;
  tmpQpk.clear();
  */
  chargeInfo->Delete();
  f->Delete();
  return retQpkTime;
}

vector < double > getFailsVals( TString bname, TString st, int pmt, int year, TString whichInfo, bool ifoff) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
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
  double tmpVals = 0.;
  vector < double > returnVals;
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

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals);
    if ( ifoff && year==2021 )
      chargeInfo->SetBranchAddress("GpsTime", &evttime);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
      chargeInfo->GetEntry(etry);
      if ( ifoff && year==2021 ) {
        if ( prevTime != evttime ) {
          if ( tmpVals == 0 )
            returnVals.push_back( tmpVals );
          prevTime = evttime;
        }
      }
      else {
       if ( tmpVals == 0 )
         returnVals.push_back( tmpVals );
      }
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}


TGraphErrors *bldGraphErrorQpk( vector< vector < double > > timeQpk ) {
  int nDays = timeQpk.size();
  vector < double > AugDays;
  vector < double > aveQpkDay;
  vector < double > rmsQpkDay;
  double tmpAve = 0.;
  for ( int day=0; day<nDays; day++ ) {
    AugDays.push_back(day+1);
    tmpAve = getmean(timeQpk[day]);
    if ( tmpAve > 0 ) {
      aveQpkDay.push_back( tmpAve );
      rmsQpkDay.push_back( getrms(timeQpk[day], tmpAve));
    }
    else {
      aveQpkDay.push_back( 0. );
      rmsQpkDay.push_back( 0. );
    }
  }
  
  TGraphErrors *retGrpErr = new TGraphErrors(nDays,&AugDays[0], &aveQpkDay[0], 0, &rmsQpkDay[0]);
  return retGrpErr;
  AugDays.clear();
  aveQpkDay.clear();
  rmsQpkDay.clear(); 
}

TGraph *bldGraphRmsQpk( vector< vector < double > > timeQpk ) {
  int nWeeks = 6;
  int do5days = 5;
  vector < vector < double > > qpk5days (nWeeks);
  int daysAug = 31;
  int crr5day = 0;

  //vector < double > tmpQpk;

  for ( int day=0; day<daysAug; day++ ) {
    for ( int cQpk=0; cQpk<timeQpk[day].size(); cQpk++ ) {
      qpk5days[ crr5day ].push_back( timeQpk[day][cQpk] );
      //tmpQpk.push_back( timeQpk[day][cQpk] );
    }
    crr5day = int(day/do5days);
  }
  
  /*
  double tmpAve2 = getmean( tmpQpk );
  cout << "HERE" << endl;
  double tmpRms = getRelRms( tmpQpk, tmpAve2 );
  cout << "sizeTmpQpk " << tmpQpk.size() << " tmpAve " << tmpAve2 << " tmpRms " << tmpRms << endl;
  cout << "END" << endl;
  */

  int n5Days = qpk5days.size();
  vector < double > Aug5Days;
  vector < double > relRmsQpkDay;
  double tmpAve = 0.;
  for ( int day=0; day<n5Days; day++ ) {
    Aug5Days.push_back(do5days*(day+1));
    tmpAve = getmean(qpk5days[day]);
    if ( tmpAve > 0 )
      relRmsQpkDay.push_back( getRelRms(qpk5days[day], tmpAve) );
    else
      relRmsQpkDay.push_back( 0. );
    //cout << "For day " << day << " Nevt " << qpk5days[day].size() << " " << getRelRms(qpk5days[day], tmpAve) << endl;
  }

  TGraph *retGraphRmsQpk = new TGraph (n5Days, &Aug5Days[0], &relRmsQpkDay[0]);

  qpk5days.clear();
  Aug5Days.clear();
  relRmsQpkDay.clear();
  return retGraphRmsQpk;
}

double getMax(TGraphErrors *grpErr1, TGraphErrors *grpErr2, TGraphErrors *grpErr3, TGraphErrors *grpErr4, TGraphErrors *grpErr5, TGraphErrors *grpErr6) {
  double retMax = 0.;
  double max1 = TMath::MaxElement(grpErr1->GetN(), grpErr1->GetY());
  double max2 = TMath::MaxElement(grpErr2->GetN(), grpErr2->GetY());
  retMax = (max1 > max2) ? max1 : max2;
  max1 = TMath::MaxElement(grpErr3->GetN(), grpErr3->GetY());
  retMax = (max1 > retMax) ? max1 : retMax;
  max1 = TMath::MaxElement(grpErr4->GetN(), grpErr4->GetY());
  retMax = (max1 > retMax) ? max1 : retMax;
  max1 = TMath::MaxElement(grpErr5->GetN(), grpErr5->GetY());
  retMax = (max1 > retMax) ? max1 : retMax;
  max1 = TMath::MaxElement(grpErr6->GetN(), grpErr6->GetY());
  retMax = (max1 > retMax) ? max1 : retMax;
  return retMax;
}

double getMax(TGraph *grp1, TGraph *grp2, TGraph *grp3, TGraph *grp4, TGraph *grp5, TGraph *grp6) {
  double retMax = 0.;
  double max1 = TMath::MaxElement(grp1->GetN(), grp1->GetY());
  double max2 = TMath::MaxElement(grp2->GetN(), grp2->GetY());
  retMax = (max1 > max2) ? max1 : max2;
  max1 = TMath::MaxElement(grp3->GetN(), grp3->GetY());
  retMax = (max1 > retMax) ? max1 : retMax;
  max1 = TMath::MaxElement(grp4->GetN(), grp4->GetY());
  retMax = (max1 > retMax) ? max1 : retMax;
  max1 = TMath::MaxElement(grp5->GetN(), grp5->GetY());
  retMax = (max1 > retMax) ? max1 : retMax;  
  max1 = TMath::MaxElement(grp6->GetN(), grp6->GetY());
  retMax = (max1 > retMax) ? max1 : retMax;
  return retMax;
}


// =================================
// *** *** *** MAIN CODE *** *** ***

void GetQpkRmsVsTime(int st, int pmt) {

  TString statId;
  statId.Form("St%d", st);

  TString bnOffl = "~/2021/sdeu/offline/forUb/results/offlineUb";
  TString bnUubOffl = "~/2021/sdeu/offline/forUub/AugResults/offlineUub";
  TString bnCdas = "~/2021/sdeu/nouub/underHistos/results/ubChPkPMT";
  TString bnUubCdas = "~/2021/sdeu/underHisto/results/uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strMean;
  TString strRms;
  TString strFails;
  TString strGoods;
  TString strSt;
  TString strPmt;
  strPmt.Form("%d", pmt);

  bool ifoff = false;

  // Getting the values
  ifoff = true;
  vector < vector <double> > fitQpkOffUb2019
    = getFitQpkTime(bnOffl, statId, pmt, 2019, ifoff);
  vector < vector <double> > fitQpkOffUb2020
    = getFitQpkTime(bnOffl, statId, pmt, 2020, ifoff);
  vector < vector <double> > fitQpkOffUub2021 
    = getFitQpkTime(bnUubOffl, statId, pmt, 2021, ifoff);

  ifoff = false;
  vector < vector <double> > fitQpkCdasUb2019
    = getFitQpkTime(bnCdas, statId, pmt, 2019, ifoff);
  vector < vector <double> > fitQpkCdasUb2020
    = getFitQpkTime(bnCdas, statId, pmt, 2020, ifoff);
  vector < vector <double> > fitQpkCdasUub2021 
    = getFitQpkTime(bnUubCdas, statId, pmt, 2021, ifoff);

  // Constructing the TGraphErrors
  TGraphErrors *grpQpkDayOffUb2019 = bldGraphErrorQpk( fitQpkOffUb2019 ); 
  TGraphErrors *grpQpkDayOffUb2020 = bldGraphErrorQpk( fitQpkOffUb2020 ); 
  TGraphErrors *grpQpkDayOffUub2021 = bldGraphErrorQpk( fitQpkOffUub2021 ); 
  TGraphErrors *grpQpkDayCdasUb2019 = bldGraphErrorQpk( fitQpkCdasUb2019 ); 
  TGraphErrors *grpQpkDayCdasUb2020 = bldGraphErrorQpk( fitQpkCdasUb2020 ); 
  TGraphErrors *grpQpkDayCdasUub2021 = bldGraphErrorQpk( fitQpkCdasUub2021 ); 

  TGraph *grpRms5DayOffUb2019 = bldGraphRmsQpk( fitQpkOffUb2019 ); 
  TGraph *grpRms5DayOffUb2020 = bldGraphRmsQpk( fitQpkOffUb2020 ); 
  TGraph *grpRms5DayOffUub2021 = bldGraphRmsQpk( fitQpkOffUub2021 ); 
  TGraph *grpRms5DayCdasUb2019 = bldGraphRmsQpk( fitQpkCdasUb2019 ); 
  TGraph *grpRms5DayCdasUb2020 = bldGraphRmsQpk( fitQpkCdasUb2020 ); 
  TGraph *grpRms5DayCdasUub2021 = bldGraphRmsQpk( fitQpkCdasUub2021 ); 

  double yUpLimFactor = 1.15;
  double yUpLimQpk = getMax(grpQpkDayOffUb2019, grpQpkDayOffUb2020, grpQpkDayOffUub2021,
      grpRms5DayCdasUb2019, grpRms5DayCdasUb2020, grpRms5DayCdasUub2021);
  yUpLimQpk *= yUpLimFactor;

  double yUpLimRms = getMax(grpRms5DayOffUb2019, grpRms5DayOffUb2020, grpRms5DayOffUub2021,
      grpRms5DayCdasUb2019, grpRms5DayCdasUb2020, grpRms5DayCdasUub2021);
  yUpLimRms *= yUpLimFactor;
    
  TCanvas *c1 = canvasStyle("c1");
  c1->Divide(1,2,0,0);
  c1->cd(1);
  grpQpkDayOffUb2019->SetTitle("");
  grpQpkDayOffUb2019->GetYaxis()->SetTitle("<Q^{pk}_{VEM}/1.01>_{day} [FADC]");
  grpQpkDayOffUb2019->GetYaxis()->SetRangeUser(106, yUpLimQpk);
  //grpQpkDayOffUb2019->GetXaxis()->SetRangeUser(0, 32);
  grpQpkDayOffUb2019->GetXaxis()->SetTitle("August [days]");
  grpQpkDayOffUb2019->SetMarkerStyle(72);
  grpQpkDayOffUb2019->SetMarkerSize(1.5);
  grpQpkDayOffUb2019->SetMarkerColor(kGray);
  tgraphErrorStyle(grpQpkDayOffUb2019);
  grpQpkDayOffUb2019->GetYaxis()->SetTitleSize(.07);
  grpQpkDayOffUb2019->GetYaxis()->SetLabelSize(0.07);
  grpQpkDayOffUb2019->GetYaxis()->SetTitleOffset(.55);
  grpQpkDayOffUb2019->GetXaxis()->SetTickLength(0.06);
  grpQpkDayOffUb2019->Draw("AP"); 

  grpQpkDayOffUb2020->SetMarkerStyle(71);
  grpQpkDayOffUb2020->SetMarkerSize(1.5);
  grpQpkDayOffUb2020->SetMarkerColor(kGray+2);
  grpQpkDayOffUb2020->Draw("P");

  grpQpkDayOffUub2021->SetMarkerStyle(73);
  grpQpkDayOffUub2021->SetMarkerColor(kBlack);
  grpQpkDayOffUub2021->SetMarkerSize(1.5);
  grpQpkDayOffUub2021->Draw("P");

  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.01);
  gPad->SetRightMargin(0.01);

  leg = new TLegend(0.3,0.7,0.8,0.95);
  leg->SetHeader("St. "+statId+", PMT"+strPmt);
  leg->AddEntry(grpQpkDayOffUb2019,"2019","p");
  leg->AddEntry(grpQpkDayOffUb2020,"2020","p");
  leg->AddEntry(grpQpkDayOffUub2021,"2021 (Q^{pk}_{VEM}/11.4)","p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  leg->SetColumnSeparation(.03);
  leg->Draw();

  c1->cd(2);
  grpRms5DayOffUb2019->SetTitle("");
  grpRms5DayOffUb2019->GetYaxis()->SetTitle("RMS/#mu_{5-days} [au]");
  grpRms5DayOffUb2019->GetYaxis()->SetRangeUser(0., yUpLimRms);
  grpRms5DayOffUb2019->GetXaxis()->SetTitle("August [days]");
  //grpRms5DayOffUb2019->GetXaxis()->SetRangeUser(0, 32);
  grpRms5DayOffUb2019->SetMarkerStyle(72);
  grpRms5DayOffUb2019->SetMarkerSize(1.5);
  grpRms5DayOffUb2019->SetMarkerColor(kGray);
  tgraphStyle(grpRms5DayOffUb2019);
  grpRms5DayOffUb2019->GetYaxis()->SetLabelSize(0.05);
  grpRms5DayOffUb2019->GetYaxis()->SetTitleSize(.06);
  grpRms5DayOffUb2019->GetYaxis()->SetTitleOffset(.6);
  grpRms5DayOffUb2019->GetXaxis()->SetTickLength(0.06);
  grpRms5DayOffUb2019->Draw("AP");

  grpRms5DayOffUb2020->SetMarkerStyle(71);
  grpRms5DayOffUb2020->SetMarkerSize(1.5);
  grpRms5DayOffUb2020->SetMarkerColor(kGray+2);
  grpRms5DayOffUb2020->Draw("P");

  grpRms5DayOffUub2021->SetMarkerStyle(73);
  grpRms5DayOffUub2021->SetMarkerSize(1.5);  
  grpRms5DayOffUub2021->SetMarkerColor(kBlack);
  grpRms5DayOffUub2021->Draw("P");

  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.01);
  //c1->Print("../plots/qpkRmsOffSt"+statId+"Pmt"+strPmt+".pdf");


  TCanvas *c2 = canvasStyle("c2");
  c2->Divide(1,2,0,0);
  c2->cd(1);
  grpQpkDayCdasUb2019->SetTitle("");
  grpQpkDayCdasUb2019->GetYaxis()->SetTitle("<Q^{pk}_{VEM}/1.01>_{day} [FADC]");
  grpQpkDayCdasUb2019->GetYaxis()->SetRangeUser(106, yUpLimQpk);
  grpQpkDayCdasUb2019->GetXaxis()->SetRangeUser(0, 32);
  grpQpkDayCdasUb2019->GetXaxis()->SetTitle("August [days]");
  grpQpkDayCdasUb2019->SetMarkerStyle(72);
  grpQpkDayCdasUb2019->SetMarkerSize(1.5);
  grpQpkDayCdasUb2019->SetMarkerColor(kGray);
  tgraphErrorStyle(grpQpkDayCdasUb2019);
  grpQpkDayCdasUb2019->GetYaxis()->SetTitleSize(.07);
  grpQpkDayCdasUb2019->GetYaxis()->SetLabelSize(0.07);
  grpQpkDayCdasUb2019->GetYaxis()->SetTitleOffset(.55);
  grpQpkDayCdasUb2019->GetXaxis()->SetTickLength(0.06);
  grpQpkDayCdasUb2019->Draw("AP"); 

  grpQpkDayCdasUb2020->SetMarkerStyle(71);
  grpQpkDayCdasUb2020->SetMarkerSize(1.5);
  grpQpkDayCdasUb2020->SetMarkerColor(kGray+2);
  grpQpkDayCdasUb2020->Draw("P");

  grpQpkDayCdasUub2021->SetMarkerStyle(73);
  grpQpkDayCdasUub2021->SetMarkerColor(kBlack);
  grpQpkDayCdasUub2021->SetMarkerSize(1.5);
  grpQpkDayCdasUub2021->Draw("P");

  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.01);
  gPad->SetRightMargin(0.01);

  leg = new TLegend(0.3,0.7,0.8,0.95);
  leg->SetHeader("St. "+statId+", PMT"+strPmt);
  leg->AddEntry(grpQpkDayCdasUb2019,"2019","p");
  leg->AddEntry(grpQpkDayCdasUb2020,"2020","p");
  leg->AddEntry(grpQpkDayCdasUub2021,"2021 (Q^{pk}_{VEM}/11.4)","p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  leg->SetColumnSeparation(.03);
  leg->Draw();

  c2->cd(2);
  grpRms5DayCdasUb2019->SetTitle("");
  grpRms5DayCdasUb2019->GetYaxis()->SetTitle("RMS/#mu_{5-days} [au]");
  grpRms5DayCdasUb2019->GetYaxis()->SetRangeUser(0., yUpLimRms);
  grpRms5DayCdasUb2019->GetXaxis()->SetTitle("August [days]");
  grpRms5DayCdasUb2019->GetXaxis()->SetRangeUser(0, 32);
  grpRms5DayCdasUb2019->SetMarkerStyle(72);
  grpRms5DayCdasUb2019->SetMarkerSize(1.5);
  grpRms5DayCdasUb2019->SetMarkerColor(kGray);
  tgraphStyle(grpRms5DayCdasUb2019);
  grpRms5DayCdasUb2019->GetYaxis()->SetLabelSize(0.05);
  grpRms5DayCdasUb2019->GetYaxis()->SetTitleSize(.06);
  grpRms5DayCdasUb2019->GetYaxis()->SetTitleOffset(.6);
  grpRms5DayCdasUb2019->GetXaxis()->SetTickLength(0.06);
  grpRms5DayCdasUb2019->Draw("AP");

  grpRms5DayCdasUb2020->SetMarkerStyle(71);
  grpRms5DayCdasUb2020->SetMarkerSize(1.5);
  grpRms5DayCdasUb2020->SetMarkerColor(kGray+2);
  grpRms5DayCdasUb2020->Draw("P");

  grpRms5DayCdasUub2021->SetMarkerStyle(73);
  grpRms5DayCdasUub2021->SetMarkerSize(1.5);  
  grpRms5DayCdasUub2021->SetMarkerColor(kBlack);
  grpRms5DayCdasUub2021->Draw("P");

  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.01);
  //c2->Print("../plots/qpkRmsCdasSt"+statId+"Pmt"+strPmt+".pdf");

  //exit(0);
}
