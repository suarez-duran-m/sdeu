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


void histoStyle(TGraph *hist) {
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.06);
}


vector < double > getFitVals( TString bname, TString st, int pmt, int year, TString whichInfo, bool ifoff ) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strYear;
  TString fname;
  TString strChargeData;
  pmtId.Form("%d", pmt);
  strYear.Form("%d", year);

  TFile *f;
  TTree *chargeInfo;
  double tmpVals = 0.;
  vector < double > returnVals;
  int nMonths = 0;

  if ( !ifoff )
    strChargeData = "ChargeData";
  else
    strChargeData = "charge";

  int evttime = 0;
  int prevTime = 0;
  int EvtId = 0;

  for ( int month=nMonths; month<nMonths+1; month++ ) {
    if ( !ifoff )
      fname = bname + pmtId + st + "lrb35" + monthUub[month] + strYear;
    else if ( ifoff )
      fname = bname + monthUub[month] + strYear + st + "Pmt" + pmtId;

    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals);
    // To avoid read the same event twice from Offline
    if ( ifoff ) {
      chargeInfo->SetBranchAddress("GpsTime", &evttime);
      chargeInfo->SetBranchAddress("evtId", &EvtId);
    }
    else
      chargeInfo->SetBranchAddress("eventId", &EvtId);

    prevTime = 0;
    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {  
      chargeInfo->GetEntry(etry);
      //if ( year==2021 && EvtId==64693194)
        //cout << "ifoff " << ifoff << " EvtId " << EvtId << " tmpVals " << tmpVals << endl;
      // To avoid read the same event twice from Offline
      if ( ifoff ) {
        if ( prevTime != evttime ) {
          returnVals.push_back( tmpVals );
          prevTime = evttime;
        }
      }
      else
        returnVals.push_back( tmpVals );
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}


vector < int > getFitVals( TString bname, TString st, int pmt, int year, bool ifoff) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strYear;
  TString fname;
  TString strChargeData;
  TString strTime;
  pmtId.Form("%d", pmt);
  strYear.Form("%d", year);
  
  TFile *f;
  TTree *chargeInfo;
  int tmpVals = 0;
  vector < int > returnVals;
  int tmpEvtId = 0;
  int nMonths = 0;
  int prevTime = 0;

  if ( !ifoff ) {
    strChargeData = "ChargeData";
    strTime = "timeEvnt";
  }
  else {
    strChargeData = "charge";
    strTime = "GpsTime";
  }

  for ( int month=nMonths; month<nMonths+1; month++ ) {
    if ( !ifoff )
      fname = bname + pmtId + st + "lrb35" + monthUub[month] + strYear;
    else if ( ifoff )
      fname = bname + monthUub[month] + strYear + st + "Pmt" + pmtId;

    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);

    tmpVals = 0.;  
    chargeInfo->SetBranchAddress(strTime, &tmpVals);
    
    prevTime = 0;    
    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
      chargeInfo->GetEntry(etry);
      // To avoid read the same event twice from Offline
      if ( ifoff ) {
        if ( prevTime != tmpVals ) {
          returnVals.push_back( tmpVals );
          prevTime = tmpVals;
        }
      }
      else
        returnVals.push_back( tmpVals );
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
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
  //cout << sqrt(rms/ngoodb) << endl;
  return sqrt(rms/ngoodb)/meanarr;
}

vector <double> getQpkPerDay( vector < int >time, vector < double >qpkVals, int totdays, int year) {
  vector<double> qpkPerDay;
  int aug1st = 0;
  double uub2ub = 1.;
  if ( year == 2019 )
    aug1st = 1248652818;
  else if ( year == 2020 )
    aug1st = 1280275218;
  else {
    aug1st = 1311811218;
    uub2ub = 11.4;
  }
  
  vector < vector < double > > tmpQpk (totdays);
  int tmpDay = 0;

  for ( int i=0; i<totdays; i++ )
    qpkPerDay.push_back(-10.);

  for ( int i=0; i<time.size(); i++ ) {
    if ( time[i] <= 0 || int(qpkVals[i]) <= 0 )
      continue;

    // 86400 for Average per 1 day.
    tmpDay = int((time[i]-aug1st)/86400);
    if ( tmpDay < 31 )
      tmpQpk[tmpDay].push_back( qpkVals[i] );
  }

  for ( int dys=0; dys<totdays; dys++ ) {
    if ( tmpQpk[dys].size() > 0 )
      qpkPerDay[dys] = getmean(tmpQpk[dys] ) / uub2ub;
  }
  tmpQpk.clear();
  return qpkPerDay;
}

vector < double > getRmsQpkPerWeek( vector < int >time, vector < double >qpkVals, int totweeks, int year ) {
  vector<double> rmsQpkPerWeek;
  int aug1st = 0;
  if ( year == 2019 )
    aug1st = 1248652818;
  else if ( year == 2020 )
    aug1st = 1280275218;
  else
    aug1st = 1311811218;

  vector < vector < double > > qpk5days (totweeks);
  int crrnt5days = 0;

  vector < double > allQpk;

  for ( int i=0; i<=totweeks; i++ )
    rmsQpkPerWeek.push_back(-10.);

  for ( int i=0; i<time.size(); i++ ) {
    if ( time[i] <= 0 || int(qpkVals[i]) <= 0 )
      continue;

    // 432000 for Average per 5 days.
    crrnt5days = int((time[i]-aug1st)/432000);
    if ( crrnt5days < 31 )
      qpk5days[ crrnt5days ].push_back( qpkVals[i] );
  }

  for ( int wk=0; wk<totweeks; wk++ ) 
    if ( qpk5days[wk].size() > 0 ) {
      rmsQpkPerWeek[wk+1] = getrms( qpk5days[wk], getmean(qpk5days[wk]) );
      /*if ( year==2021 ) {
        cout << "5Week " << wk << " Nelem " << qpk5days[wk].size() << " getmean " << getmean(qpk5days[wk]) << " getrms " << getrms( qpk5days[wk], getmean(qpk5days[wk]) ) << endl;
        //for ( auto & elem : qpk5days[wk] )
          //allQpk.push_back(elem);
      }
      */
    }
  /*
  double ave = 0.;
  if ( year==2021 ) {
    ave = getmean(allQpk);
    cout << "ave " << ave << endl;
    ave = getrms(allQpk, ave);
    cout << "rms " << ave << endl;
  }
  */

  qpk5days.clear();
  return rmsQpkPerWeek;
}

// =================================
// *** *** *** MAIN CODE *** *** ***

void readingChOffCdas(int st, int pmt) {

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
  TString strPmt;
  strPmt.Form("%d", pmt);
  bool gettime = true;
  bool ifoff = false;
  TString whinfo = "chargeVal";
  double vemFactor = 11.4;

  vector < int > time2019off;
  vector < int > time2020off;
  vector < int > time2021off;
  vector < double > qpk2019off;
  vector < double > qpk2020off;
  vector < double > qpk2021off;
  vector < double > qpkFail2019off;
  vector < double > qpkFail2020off;
  vector < double > qpkFail2021off;

  vector < int > time2019cdas;
  vector < int > time2020cdas;
  vector < int > time2021cdas;
  vector < double > qpk2019cdas;
  vector < double > qpk2020cdas;
  vector < double > qpk2021cdas;
  vector < double > qpkFail2019cdas;
  vector < double > qpkFail2020cdas;
  vector < double > qpkFail2021cdas;

  ifoff = true;
  time2019off = getFitVals( bnOffl, statId, pmt, 2019, ifoff);
  time2020off = getFitVals( bnOffl, statId, pmt, 2020, ifoff);
  time2021off = getFitVals( bnUubOffl, statId, pmt, 2021, ifoff);
  ifoff = false;
  time2019cdas = getFitVals( bnCdas, statId, pmt, 2019, ifoff);
  time2020cdas = getFitVals( bnCdas, statId, pmt, 2020, ifoff);
  time2021cdas = getFitVals( bnUubCdas, statId, pmt, 2021, ifoff);

  whinfo = "chargeVal";
  ifoff = true;
  qpk2019off = getFitVals( bnOffl, statId, pmt, 2019, whinfo, ifoff);
  qpk2020off = getFitVals( bnOffl, statId, pmt, 2020, whinfo, ifoff);
  qpk2021off = getFitVals( bnUubOffl, statId, pmt, 2021, whinfo, ifoff);
  qpkFail2019off = getFailsVals( bnOffl, statId, pmt, 2019, whinfo, ifoff);
  qpkFail2020off = getFailsVals( bnOffl, statId, pmt, 2020, whinfo, ifoff);
  qpkFail2021off = getFailsVals( bnUubOffl, statId, pmt, 2021, whinfo, ifoff);

  ifoff = false;  
  qpk2019cdas = getFitVals( bnCdas, statId, pmt, 2019, whinfo, ifoff);
  qpk2020cdas = getFitVals( bnCdas, statId, pmt, 2020, whinfo, ifoff);
  qpk2021cdas = getFitVals( bnUubCdas, statId, pmt, 2021, whinfo, ifoff);
  qpkFail2019cdas = getFailsVals( bnCdas, statId, pmt, 2019, whinfo, ifoff);
  qpkFail2020cdas = getFailsVals( bnCdas, statId, pmt, 2020, whinfo, ifoff);
  qpkFail2021cdas = getFailsVals( bnUubCdas, statId, pmt, 2021, whinfo, ifoff);

  int nQpkFails2019off = qpkFail2019off.size();
  int nQpkFails2020off = qpkFail2020off.size();
  int nQpkFails2021off = qpkFail2021off.size();

  int nQpkFails2019cdas = qpkFail2019cdas.size();
  int nQpkFails2020cdas = qpkFail2020cdas.size();
  int nQpkFails2021cdas = qpkFail2021cdas.size();

  double meanQpk2019off = getmean( qpk2019off );
  double rmsQpk2019off = getrms( qpk2019off, meanQpk2019off );
  double meanQpk2020off = getmean( qpk2020off );
  double rmsQpk2020off = getrms( qpk2020off, meanQpk2020off );
  double meanQpk2021off = getmean( qpk2021off );
  double rmsQpk2021off = getrms( qpk2021off, meanQpk2021off );

  double meanQpk2019cdas = getmean( qpk2019cdas );
  double rmsQpk2019cdas = getrms( qpk2019cdas, meanQpk2019cdas );
  double meanQpk2020cdas = getmean( qpk2020cdas );
  double rmsQpk2020cdas = getrms( qpk2020cdas, meanQpk2020cdas );
  double meanQpk2021cdas = getmean( qpk2021cdas ); 
  double rmsQpk2021cdas = getrms( qpk2021cdas, meanQpk2021cdas ); 
 
  int nPointsOff = time2019off.size()+time2020off.size()+time2021off.size();
  cout << "nPointsOff : " << nPointsOff << " 2019: " << time2019off.size() << " 2020: " << time2020off.size() << " 2021: " << time2021off.size() << endl;

  int nPointsCdas = time2019cdas.size()+time2020cdas.size()+time2021cdas.size();
  cout << "nPointsCdas: " << nPointsCdas << " 2019: " << time2019cdas.size() << " 2020: " << time2020cdas.size() << " 2021: " << time2021cdas.size() << endl;

  int augDays = 31;
  vector <double> xAugDays;
  for ( int i=0; i<augDays; i++ )
    xAugDays.push_back(i+1.);

  vector <double> augQpk2019off = getQpkPerDay(time2019off, qpk2019off, augDays, 2019);
  vector <double> augQpk2020off = getQpkPerDay(time2020off, qpk2020off, augDays, 2020);
  vector <double> augQpk2021off = getQpkPerDay(time2021off, qpk2021off, augDays, 2021);
  vector <double> augQpk2019cdas = getQpkPerDay(time2019cdas, qpk2019cdas, augDays, 2019);
  vector <double> augQpk2020cdas = getQpkPerDay(time2020cdas, qpk2020cdas, augDays, 2020);
  vector <double> augQpk2021cdas = getQpkPerDay(time2021cdas, qpk2021cdas, augDays, 2021);

  // Doing RMS
  int nWeeks = 7;
  vector <double> xAugDaysRms;
  for ( int i=0; i<nWeeks-1; i++ )
    xAugDaysRms.push_back(i*5);
  xAugDaysRms.push_back(31.);
  
  vector <double> augRmsQpk2019off = getRmsQpkPerWeek(time2019off, qpk2019off, nWeeks, 2019); 
  vector <double> augRmsQpk2020off = getRmsQpkPerWeek(time2020off, qpk2020off, nWeeks, 2020);
  vector <double> augRmsQpk2021off = getRmsQpkPerWeek(time2021off, qpk2021off, nWeeks, 2021);
  cout << "For CDAS: " << endl;
  vector <double> augRmsQpk2019cdas = getRmsQpkPerWeek(time2019cdas, qpk2019cdas, nWeeks, 2019);
  vector <double> augRmsQpk2020cdas = getRmsQpkPerWeek(time2020cdas, qpk2020cdas, nWeeks, 2020);
  vector <double> augRmsQpk2021cdas = getRmsQpkPerWeek(time2021cdas, qpk2021cdas, nWeeks, 2021);

  /*
  double tmp = 0.;
  for ( auto & elem : augRmsQpk2021cdas ) {
    tmp += elem;
    cout << elem << endl;
  }
  cout << "ave: " << (tmp+10)/(augRmsQpk2021cdas.size()-1) << endl;
  */
 
  // Doing TGraphs
  TGraph *grp2019off = new TGraph(augDays, &xAugDays[0], &augQpk2019off[0]);
  TGraph *grp2020off = new TGraph(augDays, &xAugDays[0], &augQpk2020off[0]);
  TGraph *grp2021off = new TGraph(augDays, &xAugDays[0], &augQpk2021off[0]);
  TGraph *grp2019cdas = new TGraph(augDays, &xAugDays[0], &augQpk2019cdas[0]);
  TGraph *grp2020cdas = new TGraph(augDays, &xAugDays[0], &augQpk2020cdas[0]);
  TGraph *grp2021cdas = new TGraph(augDays, &xAugDays[0], &augQpk2021cdas[0]);

  TGraph *grp2019RmsOff = new TGraph(nWeeks, &xAugDaysRms[0], &augRmsQpk2019off[0]);
  TGraph *grp2020RmsOff = new TGraph(nWeeks, &xAugDaysRms[0], &augRmsQpk2020off[0]);
  TGraph *grp2021RmsOff = new TGraph(nWeeks, &xAugDaysRms[0], &augRmsQpk2021off[0]);
  TGraph *grp2019RmsCdas = new TGraph(nWeeks, &xAugDaysRms[0], &augRmsQpk2019cdas[0]);
  TGraph *grp2020RmsCdas = new TGraph(nWeeks, &xAugDaysRms[0], &augRmsQpk2020cdas[0]);
  TGraph *grp2021RmsCdas = new TGraph(nWeeks, &xAugDaysRms[0], &augRmsQpk2021cdas[0]); 

  vector <double> tmpMax;
  double yrangeFactor = 1.4;

  tmpMax.push_back( *max_element(augQpk2019off.begin(), augQpk2019off.end()) );
  tmpMax.push_back( *max_element(augQpk2020off.begin(), augQpk2020off.end()) );
  tmpMax.push_back( *max_element(augQpk2021off.begin(), augQpk2021off.end()) );
  double yUpperLimOff = *max_element(tmpMax.begin(), tmpMax.end());
  if ( yUpperLimOff > 0 )
    yUpperLimOff *= yrangeFactor;
  else 
    yUpperLimOff = 200.;
  tmpMax.clear();

  tmpMax.push_back( *max_element(augQpk2019cdas.begin(), augQpk2019cdas.end()) );
  tmpMax.push_back( *max_element(augQpk2020cdas.begin(), augQpk2020cdas.end()) );
  tmpMax.push_back( *max_element(augQpk2021cdas.begin(), augQpk2021cdas.end()) );
  double yUpperLimCdas = *max_element(tmpMax.begin(), tmpMax.end());
  if ( yUpperLimCdas > 0 )
    yUpperLimCdas *= yrangeFactor;
  else 
    yUpperLimCdas = 200.;
  tmpMax.clear();

  tmpMax.push_back( *max_element(augRmsQpk2019off.begin(), augRmsQpk2019off.end()) );
  tmpMax.push_back( *max_element(augRmsQpk2020off.begin(), augRmsQpk2020off.end()) );
  tmpMax.push_back( *max_element(augRmsQpk2021off.begin(), augRmsQpk2021off.end()) );
  double yUpperLimRmsOff = *max_element(tmpMax.begin(), tmpMax.end());
   if ( yUpperLimOff > 0 )         
     yUpperLimRmsOff *= yrangeFactor; 
   else                            
     yUpperLimRmsOff = 200.;          
  tmpMax.clear();

  tmpMax.push_back( *max_element(augRmsQpk2019cdas.begin(), augRmsQpk2019cdas.end()) );
  tmpMax.push_back( *max_element(augRmsQpk2020cdas.begin(), augRmsQpk2020cdas.end()) );
  tmpMax.push_back( *max_element(augRmsQpk2021cdas.begin(), augRmsQpk2021cdas.end()) );
  double yUpperLimRmsCdas = *max_element(tmpMax.begin(), tmpMax.end());
  if ( yUpperLimRmsCdas > 0 )
    yUpperLimRmsCdas *= yrangeFactor;
  else
    yUpperLimRmsCdas = 200.;          

  // Doing plots
  TCanvas *c1 = canvasStyle("c1");
  c1->Divide(1,2,0,0);
  c1->cd(1);

  grp2019off->SetTitle("");
  grp2019off->GetYaxis()->SetRangeUser(100, yUpperLimOff);
  grp2019off->SetMarkerStyle(72);
  grp2019off->SetMarkerSize(2);
  grp2019off->SetMarkerColor(kGray);
  grp2019off->GetYaxis()->SetTitle("<Q_{VEM}^{pk} /1.01>_{day} [FADC]");
  grp2019off->GetXaxis()->SetTitle("August [days]");
  grp2019off->GetXaxis()->SetRangeUser(0, 35);
  histoStyle(grp2019off);
  grp2019off->GetYaxis()->SetTitleOffset(.45);
  grp2019off->GetYaxis()->SetTitleSize(.075);
  grp2019off->GetXaxis()->SetTickLength(0.06);
  grp2019off->Draw("AP same");
  
  grp2020off->SetMarkerStyle(71);
  grp2020off->SetMarkerSize(2);
  grp2020off->SetMarkerColor(kGray+2);
  grp2020off->Draw("P same");

  grp2021off->SetMarkerStyle(73);
  grp2021off->SetMarkerSize(2);
  grp2021off->SetMarkerColor(kBlack);
  grp2021off->Draw("P same");
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.01);
  gPad->SetRightMargin(0.01);

  leg = new TLegend(0.82,0.45,0.95,0.85);
  leg->SetHeader("St. "+statId+", PMT"+strPmt);
  leg->AddEntry(grp2019off,"(OffLine)","");
  leg->AddEntry(grp2019off,"2019","p");
  leg->AddEntry(grp2020off,"2020","p");
  leg->AddEntry(grp2021off,"2021 (Q^{pk}_{VEM}/11.4)","p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  c1->cd(2);
  grp2019RmsOff->SetTitle("");
  grp2019RmsOff->SetMarkerStyle(72);
  grp2019RmsOff->SetMarkerSize(2);
  grp2019RmsOff->SetMarkerColor(kGray);
  grp2019RmsOff->GetYaxis()->SetTitle("RMS/#mu_{5days} [FADC]");
  grp2019RmsOff->GetXaxis()->SetTitle("August [days]");
  grp2019RmsOff->GetYaxis()->SetRangeUser(0, yUpperLimRmsOff);
  grp2019RmsOff->GetXaxis()->SetRangeUser(0, 35);
  histoStyle(grp2019RmsOff);
  grp2019RmsOff->GetYaxis()->SetTitleSize(.06);
  grp2019RmsOff->GetYaxis()->SetTitleOffset(.55);
  grp2019RmsOff->GetXaxis()->SetTickLength(0.06);
  grp2019RmsOff->Draw("AP same");
  
  grp2020RmsOff->SetMarkerStyle(71);
  grp2020RmsOff->SetMarkerSize(2);
  grp2020RmsOff->SetMarkerColor(kGray+2);
  grp2020RmsOff->Draw("P same");

  grp2021RmsOff->SetMarkerStyle(73);
  grp2021RmsOff->SetMarkerSize(2);
  grp2021RmsOff->SetMarkerColor(kBlack);
  grp2021RmsOff->Draw("P same");
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.01);

  leg = new TLegend(0.86,0.45,0.95,0.85);
  leg->SetHeader("St. "+statId+", PMT"+strPmt);
  leg->AddEntry(grp2019off,"(OffLine)","");
  leg->AddEntry(grp2019RmsOff,"2019","p");
  leg->AddEntry(grp2020RmsOff,"2020","p");
  leg->AddEntry(grp2021RmsOff,"2021","p");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print("../plots/qpkRmsOffSt"+statId+"Pmt"+strPmt+".pdf");


  TCanvas *c2 = canvasStyle("c2");
  c2->Divide(1,2,0,0);
  c2->cd(1); 

  grp2019cdas->SetTitle("");
  grp2019cdas->GetYaxis()->SetRangeUser(100,yUpperLimCdas);
  grp2019cdas->SetMarkerStyle(72);
  grp2019cdas->SetMarkerSize(2);
  grp2019cdas->SetMarkerColor(kGray);
  grp2019cdas->GetYaxis()->SetTitle("<Q_{VEM}^{pk}>_{day} [FADC]");
  grp2019cdas->GetXaxis()->SetTitle("August [days]");
  grp2019cdas->GetXaxis()->SetRangeUser(0, 35);
  histoStyle(grp2019cdas);
  grp2019cdas->GetYaxis()->SetTitleSize(.08);
  grp2019cdas->GetYaxis()->SetTitleOffset(.45);
  grp2019cdas->GetXaxis()->SetTickLength(0.06);
  grp2019cdas->Draw("AP same");
  
  grp2020cdas->SetMarkerStyle(71);
  grp2020cdas->SetMarkerSize(2);
  grp2020cdas->SetMarkerColor(kGray+2);
  grp2020cdas->Draw("P same");

  grp2021cdas->SetMarkerStyle(73);
  grp2021cdas->SetMarkerSize(2);
  grp2021cdas->SetMarkerColor(kBlack);
  grp2021cdas->Draw("P same");

  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.01);
  gPad->SetRightMargin(0.01);

  leg = new TLegend(0.82,0.55,0.95,0.85);
  leg->SetHeader("St. "+statId+", PMT"+strPmt);
  leg->AddEntry(grp2019cdas,"2019","p");
  leg->AddEntry(grp2020cdas,"2020","p");
  leg->AddEntry(grp2021cdas,"2021 (Q^{pk}_{VEM}/11.4)","p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  c2->cd(2);
  grp2019RmsCdas->SetTitle("");
  grp2019RmsCdas->SetMarkerStyle(72);
  grp2019RmsCdas->SetMarkerSize(2);
  grp2019RmsCdas->SetMarkerColor(kGray);
  grp2019RmsCdas->GetYaxis()->SetTitle("RMS/#mu_{5days} [FADC]");
  grp2019RmsCdas->GetXaxis()->SetTitle("August [days]");
  grp2019RmsCdas->GetYaxis()->SetRangeUser(0, yUpperLimRmsCdas);
  grp2019RmsCdas->GetXaxis()->SetRangeUser(0, 35);
  histoStyle(grp2019RmsCdas);
  grp2019RmsCdas->GetYaxis()->SetTitleOffset(.55);
  grp2019RmsCdas->GetXaxis()->SetTickLength(0.06);
  grp2019RmsCdas->Draw("AP same");
  
  grp2020RmsCdas->SetMarkerStyle(71);
  grp2020RmsCdas->SetMarkerSize(2);
  grp2020RmsCdas->SetMarkerColor(kGray+2);
  grp2020RmsCdas->Draw("P same");

  grp2021RmsCdas->SetMarkerStyle(73);
  grp2021RmsCdas->SetMarkerSize(2);
  grp2021RmsCdas->SetMarkerColor(kBlack);
  grp2021RmsCdas->Draw("P same");
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.01);

  leg = new TLegend(0.86,0.55,0.95,0.85);
  leg->SetHeader("St. "+statId+", PMT"+strPmt);
  leg->AddEntry(grp2019RmsCdas,"2019","p");
  leg->AddEntry(grp2020RmsCdas,"2020","p");
  leg->AddEntry(grp2021RmsCdas,"2021","p");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c2->Modified();
  c2->SetSelected(c2);
  c2->Print("../plots/qpkRmsCdasSt"+statId+"Pmt"+strPmt+".pdf");
  exit(0);
}
