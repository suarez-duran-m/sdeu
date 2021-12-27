TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 102, 76, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.09); 
  canvas->SetRightMargin(0.14);
  canvas->SetTopMargin(0.04); 
  canvas->SetFrameBorderMode(0);
  return canvas;
}

double getMean(vector<double> vecValues) {
  double mean = 0;
  for( auto & val_i : vecValues )
    mean += val_i;

  return mean/vecValues.size();
}

void fillingDist( int pmt, int st_i, vector<vector<vector<double>>> vecVals,
    TH1D *retDist) {
  double average = 0.;
  if ( vecVals[pmt].size() > 0 && vecVals[pmt][st_i].size() > 0 ) {
    average = getMean( vecVals[pmt][st_i] );
    for ( auto & qpk_i : vecVals[pmt][st_i] )
      retDist->Fill( qpk_i/average );
  }
}


void fillingQpksFullDays(vector<double> inVals, 
    vector<vector<double>> &retVals) {
  retVals[0].push_back( inVals[0] );
  retVals[1].push_back( inVals[1] );
}

bool doMovingWindow(vector<double> qpksDay0, vector<double> qpksDay1,
    vector<double> qpksDay2, vector<double> qpksDay3, vector<double> qpksDay4, 
    vector<double> qpksDay5, int pmt, bool ifUUB, double retSlope) {

  bool criterium = false;
  double chi2CutUB = 18.;//20.;
  double chi2CutUUB = 12.; //15.;
  vector < vector < double > > qpksFullDays;
  vector < double > xAxis{1., 2., 3., 4., 5., 6.};

  qpksFullDays.resize(2);

  fillingQpksFullDays(qpksDay0, qpksFullDays);
  fillingQpksFullDays(qpksDay1, qpksFullDays);
  fillingQpksFullDays(qpksDay2, qpksFullDays);
  fillingQpksFullDays(qpksDay3, qpksFullDays);
  fillingQpksFullDays(qpksDay4, qpksFullDays);
  fillingQpksFullDays(qpksDay5, qpksFullDays);

  TGraphErrors *distQpks = new TGraphErrors(
      qpksFullDays[0].size(), &xAxis[0], &qpksFullDays[0][0], 0, &qpksFullDays[1][0]);

  bool ifFitOk = distQpks->Fit("pol1","0Q");

  if ( ifFitOk==0 ) {
    double slp = distQpks->GetFunction("pol1")->GetParameter(1);
    double chi2 = distQpks->GetFunction("pol1")->GetChisquare();

    if ( fabs(slp) < 0.5 && ifUUB ) {
      switch ( pmt ) {
        case 1:
          if ( chi2 < chi2CutUUB ) {
            criterium = true;
            retSlope = slp;
          }
          break;
        case 2:
          if ( chi2 < chi2CutUUB ) {
            criterium = true;
            retSlope = slp;
          }
          break;
        case 3:
          if ( chi2 < chi2CutUUB ) {
            criterium = true;
            retSlope = slp;
          }        
          break;
      }
    }
    if ( fabs(slp) < 0.5 && !ifUUB ) {
      switch ( pmt ) {
        case 1:
          if ( chi2 < chi2CutUB ) {
            criterium = true;
            retSlope = slp;
          }
          break;
        case 2:
          if ( chi2 < chi2CutUB ) {
            criterium = true;
            retSlope = slp;
          }
          break;
        case 3:
          if ( chi2 < chi2CutUB ) {
            criterium = true;
            retSlope = slp;
          }
          break;
      }
    }
  }

  if ( criterium )
    return criterium;
  else
    return false;
}


void fillingChi2VsSlop( vector<double>chi2, vector<double>slp, TH2D *hist ) {
  for ( int slp_i=0; slp_i<slp.size(); slp_i++ )
    hist->Fill( slp[slp_i], chi2[slp_i] );
}


void fillQpkTimeVals(bool ifIsUub, int pmt, int st_id, 
    vector<double> &retQpkVect, vector<double> &retTimeVect) {

  TString bnCdas = (ifIsUub) ?
    "~/2021/sdeu/underHisto/results/uubChPkPMT" :
    "~/2021/sdeu/nouub/underHistos/results/ubChPkPMT"; 

  TString pmtId;
  TString strStId;
  TString fname;
  pmtId.Form("%d", pmt);
  strStId.Form("St%d", (int)st_id);

  TString monthUub[4] = {"Aug", "Sep", "Oct", "Nov"};
  int nMonths = sizeof(monthUub)/sizeof(*monthUub);

  int nYears = 5;
  TString strYear[5] = {"2016", "2018", "2019", "2020", "2021"};
  int stYear = (ifIsUub) ? nYears-1 : 1;
  int lstYear = (ifIsUub) ? nYears-1 : 1;
  
  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int fetchTime = 0;

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( int month_i=0; month_i<nMonths; month_i++ ) {
      fname = bnCdas + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);
      chargeInfo->SetBranchAddress("timeEvnt", &fetchTime);
      fetchQpkVals = 0.;
      fetchTime = 0;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
        chargeInfo->GetEntry(etry);
        if ( fetchQpkVals <= 0 )
          continue;
        retQpkVect.push_back( fetchQpkVals );
        retTimeVect.push_back( fetchTime );      
      }
      f->Clear();
      f->Close();
    }
  }
}

bool doAvePerTime(bool ifIsUub, vector<double> qpkVect, 
    vector<double> timeVect, int pmt, vector<double> &retDistSelec) {
  int currentDay = (ifIsUub) ? 1311811218 : 1217116818;//1248652818; //1217116818; // August 1st
  const int oneDay = 86400;
  double aveDay = 0.;
  double qpk2 = 0.;
  double rms = 0.;
  int qpkInDay = 0;
  int timeDiff = 0;

  vector < double > distQpkSelect;
  vector < vector < double > > distQpkTemp;
  double currSlope;
  
  const int daysForMw = 6;
  int crrDayForMw = 0;
  distQpkTemp.resize(daysForMw);

  vector < double > qpksDay0(2);
  vector < double > qpksDay1(2);
  vector < double > qpksDay2(2);
  vector < double > qpksDay3(2);
  vector < double > qpksDay4(2);
  vector < double > qpksDay5(2);

  bool selectionOk = false;
  double minSlp = 10.;
 
  for ( int qpk_i=0; qpk_i < qpkVect.size(); qpk_i++ ) {
    timeDiff = timeVect[qpk_i] - currentDay;  
    if ( timeDiff > oneDay ) {
      if ( aveDay > 0 ) {
        aveDay /= qpkInDay;
        if ( qpkInDay > 1 ) {
          rms = sqrt(qpk2/qpkInDay - aveDay*aveDay);
          // Using the error of the mean
          //sqrt(qpk2/qpkInDay - aveDay*aveDay) );
          //retDistErrMean->Fill( rms/sqrt(qpkInDay) );
          switch ( crrDayForMw ) {
            case 0 :
              qpksDay0[0] = aveDay;
              qpksDay0[1] = rms/sqrt(qpkInDay);
              break;
            case 1 :
              qpksDay1[0] = aveDay;
              qpksDay1[1] = rms/sqrt(qpkInDay);
              break;
            case 2 :
              qpksDay2[0] = aveDay;
              qpksDay2[1] = rms/sqrt(qpkInDay);
              break;
            case 3 :
              qpksDay3[0] = aveDay;
              qpksDay3[1] = rms/sqrt(qpkInDay);
              break;
            case 4 :
              qpksDay4[0] = aveDay;
              qpksDay4[1] = rms/sqrt(qpkInDay);
              break;
            case 5 :
              qpksDay5[0] = aveDay;
              qpksDay5[1] = rms/sqrt(qpkInDay);
              break;
          }
        }
      }
      if ( aveDay > 0 )
        crrDayForMw++;
      else {
        crrDayForMw = 0;
        for ( int day_i=0; day_i<daysForMw; day_i++ )
          distQpkTemp[day_i].clear();
        for (int j=0; j<2; j++ ) {
          qpksDay0[j] = 0.;
          qpksDay1[j] = 0.;
          qpksDay2[j] = 0.;
          qpksDay3[j] = 0.;
          qpksDay4[j] = 0.;
          qpksDay5[j] = 0.;
        }
      }
      aveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;
      currentDay += oneDay;
      bool isDistOk = false;
      if ( crrDayForMw == daysForMw ) {
        crrDayForMw = 5;
        //do MW retChi2 retSlope
        isDistOk = doMovingWindow(qpksDay0, qpksDay1, qpksDay2, qpksDay3, 
            qpksDay4, qpksDay5, pmt, ifIsUub, currSlope);
        qpksDay0 = qpksDay1;
        qpksDay1 = qpksDay2;
        qpksDay2 = qpksDay3;
        qpksDay3 = qpksDay4;
        qpksDay4 = qpksDay5;

        if ( isDistOk ) {
          if ( minSlp > currSlope ) {
            minSlp = currSlope;
            distQpkSelect.clear();
            for ( int day_i=0; day_i<daysForMw; day_i++ )
              for ( auto & qpks : distQpkTemp[day_i] )
                distQpkSelect.push_back( qpks );              
            selectionOk = true;
          }
        }
        distQpkTemp[crrDayForMw].clear();
      }
    }
/*
    if ( pmt == 1 && ifIsUub )
      cout << crrDayForMw << " " 
        << (int)timeVect[qpk_i] << " " 
        << distQpkTemp[crrDayForMw].size() << endl;
  */    
    if ( qpkVect[qpk_i] > 0. ) {
      aveDay += qpkVect[qpk_i];
      qpk2 += qpkVect[qpk_i]*qpkVect[qpk_i];
      qpkInDay++;
      distQpkTemp[crrDayForMw].push_back( qpkVect[qpk_i] );
    }
    if ( timeDiff > 2*oneDay ) {
      currentDay += oneDay*(timeDiff/oneDay)-oneDay;
      crrDayForMw = 0;
      aveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;
      bool isDistOk = false;
      for (int j=0; j<2; j++ ) {
        qpksDay0[j] = 0.;
        qpksDay1[j] = 0.;
        qpksDay2[j] = 0.;
        qpksDay3[j] = 0.;
        qpksDay4[j] = 0.;
        qpksDay5[j] = 0.;
      }
      for ( int day_i=0; day_i<daysForMw; day_i++ )
        distQpkTemp[day_i].clear();
    }
  }

  if ( selectionOk ) {
    for ( auto & qpks : distQpkSelect )
      retDistSelec.push_back( qpks );
  }
  
  distQpkTemp.clear();
  distQpkSelect.clear();

  return selectionOk;
}


void makingQpkAccuracy() {

  // Getting the working stations
  TFile *outFile = new TFile("makingQpkAccuracy.root","RECREATE");
  
  ifstream stationSelected;
  TString strFileSelecSt = "../fullUubStationsListVert.txt"; //"selectedStations.dat";
  int readStId = 0;
  vector < int > stId;
  stationSelected.open(strFileSelecSt);
  while( stationSelected.good() ) {
    stationSelected >> readStId;
    stId.push_back( readStId );
  }
  stId.pop_back();
  stationSelected.close();

  // Vector for store Qpk and time data
  vector < vector < double > > qpkUb;
  vector < vector < double > > timeUb;
  vector < vector < double > > qpkUub;
  vector < vector < double > > timeUub;

  qpkUb.resize(3);
  timeUb.resize(3);
  qpkUub.resize(3);
  timeUub.resize(3);

  // Vector for Qpks selected distributions
  vector < double > qpkSelUbTemp;
  vector < double > qpkSelUubTemp;
  vector < vector < vector < double > > > qpkSelUb;
  vector < vector < vector < double > > > qpkSelUub;

  qpkSelUb.resize(3);
  qpkSelUub.resize(3);

  bool isSelectionUBOk = false;
  bool isSelectionUUBOk = false;

  vector < vector < int > > stationOk;
  stationOk.resize(3);
  stationOk[0].resize( stId.size() );
  stationOk[1].resize( stId.size() );
  stationOk[2].resize( stId.size() );

  for ( int pmt_j=0; pmt_j<3; pmt_j++ ) {
    qpkSelUb[pmt_j].resize( stId.size() );
    qpkSelUub[pmt_j].resize( stId.size() );
    for ( int st_j=0; st_j<stId.size(); st_j++ )
      stationOk[pmt_j][st_j] = 0;
  }
   
  int countStations = 0;

  for ( auto & st_i : stId ) {
    cout << "Moving window for station " << st_i << endl;
    //if ( st_i != 864 )
      //continue;
    
    // Fetching Qpk and time values
    for ( int pmt_i=1; pmt_i<4; pmt_i++ ) {
      isSelectionUBOk = false;
      isSelectionUUBOk = false;

      fillQpkTimeVals(false, pmt_i, st_i,
          qpkUb[pmt_i-1], timeUb[pmt_i-1]);
      fillQpkTimeVals(true, pmt_i, st_i,
          qpkUub[pmt_i-1], timeUub[pmt_i-1]);
      
      // Fitting and choosing applying moving-window       
      isSelectionUBOk = doAvePerTime(false, qpkUb[pmt_i-1], 
          timeUb[pmt_i-1], pmt_i, qpkSelUbTemp);
      isSelectionUUBOk = doAvePerTime(true, qpkUub[pmt_i-1], 
          timeUub[pmt_i-1], pmt_i, qpkSelUubTemp);
  
      if ( isSelectionUBOk==true && isSelectionUUBOk==true ) {
        for ( auto qpks : qpkSelUbTemp )
          qpkSelUb[pmt_i-1][countStations].push_back( qpks );
        for ( auto qpks : qpkSelUubTemp )
          qpkSelUub[pmt_i-1][countStations].push_back( qpks );
        stationOk[pmt_i-1][countStations] = 1;
        cout << "St " << st_i 
          << " PMT " << pmt_i 
          << " " << countStations << endl;
      }
      qpkSelUbTemp.clear();
      qpkSelUubTemp.clear();
    }

    countStations++;

    qpkUb.clear();
    qpkUub.clear();
    timeUb.clear();
    timeUub.clear();

    qpkUb.resize(3);
    qpkUub.resize(3);
    timeUb.resize(3);
    timeUub.resize(3);

    //if ( st_i == 846 )
      //break;
  }

  int nbins = 200;
  double frsBin = 0.;
  double lstBin = 2.;
  TH1D *uniqRelQpkDistUB = new TH1D("uniqRelQpkDistUB", "", nbins, frsBin, lstBin);
  TH1D *uniqRelQpkDistUUB = new TH1D("uniqRelQpkDistUUB", "", nbins, frsBin, lstBin);

  lstBin = 20.;
  TH1D *distAccUB = new TH1D("distAccUB", "", nbins, frsBin, lstBin);
  TH1D *distAccUUB = new TH1D("distAccUUB", "", nbins, frsBin, lstBin); 

  double mu = 0.;
  double sgm = 0.;
  Int_t status = -1;
  bool doFit = false;
 
  cout << endl;
  for ( int st_i=0; st_i<stationOk[0].size(); st_i++ ) {
    cout << "Doing accuracy for Station " << st_i << endl;
    doFit = false;
    
    for ( int pmt_i=0; pmt_i<3; pmt_i++ )
      if ( stationOk[pmt_i][st_i] == 1 ) {       
        fillingDist(pmt_i, st_i, qpkSelUb, uniqRelQpkDistUB);
        fillingDist(pmt_i, st_i, qpkSelUub, uniqRelQpkDistUUB);
        doFit = true;
      }
      
    if ( doFit ) {
      status = uniqRelQpkDistUB->Fit("gaus","0Q","",0,2);
      if ( status == 0 ) {
        mu = uniqRelQpkDistUB->GetFunction("gaus")->GetParameter(1);
        sgm = uniqRelQpkDistUB->GetFunction("gaus")->GetParameter(2);
        distAccUB->Fill( 100.*sgm/mu );
      }
 
      status = uniqRelQpkDistUUB->Fit("gaus","0Q","",0,2);
      if ( status == 0 ) {
        mu = uniqRelQpkDistUUB->GetFunction("gaus")->GetParameter(1);
        sgm = uniqRelQpkDistUUB->GetFunction("gaus")->GetParameter(2);
        distAccUUB->Fill( 100.*sgm/mu );
      }
    }
    uniqRelQpkDistUB->Reset();
    uniqRelQpkDistUUB->Reset();
  }

  TString strMean;
  TString strRms;
  TString strMeanErr;
  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  distAccUB->SetStats(kFALSE);
  distAccUB->SetLineColor(kBlue);
  distAccUB->GetXaxis()->SetRangeUser(0.6, 3.8);
  distAccUB->GetXaxis()->SetTitle("#sigma/#mu [%]");
  distAccUB->GetYaxis()->SetTitle("Counts [au]");
  distAccUB->SetLineWidth(2);
  distAccUB->SetFillStyle(3345);
  distAccUB->SetFillColor(kBlue);
  distAccUB->Draw();

  distAccUUB->SetLineColor(kRed);
  distAccUUB->SetLineWidth(2);
  distAccUUB->SetFillStyle(3354);
  distAccUUB->SetFillColor(kRed);
  distAccUUB->Draw("same");

  TLegend *leg = new TLegend(0.45,0.6,0.9,0.95);
  strMean.Form("%.2f", distAccUB->GetMean());
  strRms.Form("%.2f", distAccUB->GetRMS());
  strMeanErr.Form("%.2f", distAccUB->GetMeanError());
  leg->AddEntry(distAccUB, "UB", "l" );
  leg->AddEntry(distAccUB, "Mean: "+strMean+" % #pm "+strMeanErr+" %", "");
  leg->AddEntry(distAccUB, "RMS: "+strRms,"");
  
  strMean.Form("%.2f", distAccUUB->GetMean());
  strRms.Form("%.2f", distAccUUB->GetRMS());
  strMeanErr.Form("%.2f", distAccUUB->GetMeanError());
  leg->AddEntry(distAccUUB, "UUB", "l" );
  leg->AddEntry(distAccUUB, "Mean: "+strMean+" % #pm "+strMeanErr+" %", "");
  leg->AddEntry(distAccUUB, "RMS: "+strRms,"");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print("../plots2/accQpkFitUbUubAllStAllPmt_Stats.pdf");

  cout << "For UB" << endl;
  cout << distAccUB->GetMean() << " " 
    << distAccUB->GetMeanError() << endl;

  cout << "For UUB" << endl;
  cout << distAccUUB->GetMean() << " " 
    << distAccUUB->GetMeanError() << endl;

  outFile->cd();
  uniqRelQpkDistUB->Write();
  uniqRelQpkDistUUB->Write();

  distAccUB->Write();
  distAccUUB->Write();

  outFile->Write();
  outFile->Close();
  //exit(0);
}
