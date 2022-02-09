TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 102, 76, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.09);
  canvas->SetRightMargin(0.05);
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


void resetQpksDays(vector<double> &qpksDay0, vector<double> &qpksDay1, 
    vector<double> &qpksDay2, vector<double> &qpksDay3, vector<double> &qpksDay4, 
    vector<double> &qpksDay5, vector<double> &qpksDay6 ) {
  for (int j=0; j<2; j++ ) {
    qpksDay0[j] = 0.;
    qpksDay1[j] = 0.;
    qpksDay2[j] = 0.;
    qpksDay3[j] = 0.;
    qpksDay4[j] = 0.;
    qpksDay5[j] = 0.;
    qpksDay6[j] = 0.;
  }
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


int fillingQpksFullDays(vector<double> inVals, 
    vector<vector<double>> &retVals) {
  if ( inVals[0] > 0. ) {
    retVals[0].push_back( inVals[0] );
    retVals[1].push_back( inVals[1] );
    return 1;
  }
  else
    return 0;
}

bool doMovingWindow(vector<double> qpksDay0, vector<double> qpksDay1,
    vector<double> qpksDay2, vector<double> qpksDay3, vector<double> qpksDay4, 
    vector<double> qpksDay5, vector<double> qpksDay6, int pmt, bool ifUUB, 
    double retSlope) {
  // Vector to store <Qpk>-day time series
  vector < vector < double > > qpksFullDays;
  qpksFullDays.resize(2);
                                                        
  // Check for holes into time window
  int daysOk = 0;
  daysOk += fillingQpksFullDays(qpksDay0, qpksFullDays);
  daysOk += fillingQpksFullDays(qpksDay1, qpksFullDays);
  daysOk += fillingQpksFullDays(qpksDay2, qpksFullDays);
  daysOk += fillingQpksFullDays(qpksDay3, qpksFullDays);
  daysOk += fillingQpksFullDays(qpksDay4, qpksFullDays);
  daysOk += fillingQpksFullDays(qpksDay5, qpksFullDays);
  daysOk += fillingQpksFullDays(qpksDay6, qpksFullDays);
  if ( daysOk < 7 )
    return false;

  // Bool to check if the current slope is Ok
  bool criterium = false;
  double cutPval = -10.;
  double cutSlpMn = 0.;
  double cutSlpMx = 0.;

  // TGraph to apply fit and get Slope and Chi2
  vector < double > xAxis{1., 2., 3., 4., 5., 6., 7.};
  TGraphErrors *distQpks = new TGraphErrors(
      qpksFullDays[0].size(), &xAxis[0], &qpksFullDays[0][0], 0, &qpksFullDays[1][0]);

  // Applying linear fit to <Qpk>-per-day during the time-window
  // and checking if the fit was successful
  bool ifFitOk = distQpks->Fit("pol1","Q");
  if ( distQpks->GetFunction("pol1")->GetNDF() != 5 )
    cout << "Wrong NDF" << endl;
  if ( ifFitOk==0 ) {
    double slp = distQpks->GetFunction("pol1")->GetParameter(1);
    double chi2 = distQpks->GetFunction("pol1")->GetChisquare();
    chi2 = log10(TMath::Prob(chi2, 5));

    // Normalizing Slope by <qpk> during time-window
    double meanTimeWidow = 0.;
    for ( int i=0; i<qpksFullDays.size(); i++ )
      meanTimeWidow += qpksFullDays[0][i];
    meanTimeWidow /= qpksFullDays.size();
    slp /= meanTimeWidow;
    // Asking if the slope is according with 1-sigma of
    // slope distribution (makeStatsFitMovingWindow.C)
    // From here, the current slope is returned 
    if ( ifUUB ) {
      switch ( pmt ) {
        case 1:
          cutSlpMn = -7.97e-04 - 2.36e-03;
          cutSlpMx = -7.97e-04 + 2.36e-03;
          if ( slp > cutSlpMn && slp < cutSlpMx && chi2 > cutPval ) {
            criterium = true;
            retSlope = slp;
          }
          break;
        case 2:
          cutSlpMn = -8.43e-04 - 2.76e-03;
          cutSlpMx = -8.434e-04 + 2.76e-03;
          if ( slp > cutSlpMn && slp < cutSlpMx && chi2 > cutPval ) {
            criterium = true;
            retSlope = slp;
          }
          break;
        case 3:
          cutSlpMn = -6.91e-04 - 2.69e-03;
          cutSlpMx = -6.91e-04 + 2.69e-03;
          if ( slp > cutSlpMn && slp < cutSlpMx && chi2 > cutPval ) {
            criterium = true;
            retSlope = slp;
          }        
          break;
      }
    }
    if ( !ifUUB ) {
      switch ( pmt ) {
        case 1:
          cutSlpMn = -6.07e-04 - 1.75e-03;
          cutSlpMx = -6.07e-04 + 1.75e-03;
          if ( slp > cutSlpMn && slp < cutSlpMx && chi2 > cutPval ) {
            criterium = true;
            retSlope = slp;
          }
          break;
        case 2:
          cutSlpMn = -7.07e-04 - 2.15e-03;
          cutSlpMx = -7.07e-04 + 2.15e-03;
          if ( slp > cutSlpMn && slp < cutSlpMx && chi2 > cutPval ) {
            criterium = true;
            retSlope = slp;
          }
          break;
        case 3:
          cutSlpMn = -7.17e-04 - 2.34e-03;
          cutSlpMx = -7.17e-04 + 2.34e-03;
          if ( slp > cutSlpMn && slp < cutSlpMx && chi2 > cutPval ) {
            criterium = true;
            retSlope = slp;
          }
          break;
      }
    }
  }
  // Returning if the 7-day is ok. 
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

  // String variable to read root-files with Qpk values
  TString bnCdas = (ifIsUub) ?
    "~/postdoc/sdeu/underHisto/results/uubChPkPMT" :
    "~/postdoc/sdeu/nouub/underHistos/results/ubChPkPMT"; 

  TString pmtId;
  TString strStId;
  TString fname;
  pmtId.Form("%d", pmt);
  strStId.Form("St%d", (int)st_id);

  // Interval time for Qpk accuracy calculation
  TString monthUub[4] = {"Aug", "Sep", "Oct", "Nov"};
  int nMonths = sizeof(monthUub)/sizeof(*monthUub);

  int nYears = 5;
  TString strYear[5] = {"2016", "2018", "2019", "2020", "2021"};
  int stYear = (ifIsUub) ? nYears-1 : 1;
  int lstYear = (ifIsUub) ? nYears-1 : 1;

  // Reading TTree with Qpk values 
  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int fetchTime = 0;

  // Fetching Qpk as function of time
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

  // Starting date for average calculation
  int currentDay = (ifIsUub) ? 1311811218 : 1217116818;// August 1st
  const int oneDay = 86400;
  double aveDay = 0.;
  double qpk2 = 0.;
  double rms = 0.;
  int qpkInDay = 0;
  int timeDiff = 0;
  // Number of days for the time-window
  const int daysForMw = 7;
  int crrDayForMw = 0;
  bool fitWasOk = false;
  
  // vectors to storage the average of Qpk per day.
  vector < double > qpksDay0(2);
  vector < double > qpksDay1(2);
  vector < double > qpksDay2(2);
  vector < double > qpksDay3(2);
  vector < double > qpksDay4(2);
  vector < double > qpksDay5(2);
  vector < double > qpksDay6(2);
  // Vector to store the selected Qpk distribution
  vector < double > distQpkSelect;
  vector < vector < double > > distQpkTemp;
  distQpkTemp.resize(daysForMw);

  double currSlope;
  bool selectionOk = false;
  bool isDistOk = false;
  // Arbitrary slope to find the smaller one.
  double minSlp = 10.;
  // Doing Qpk average per day and applying MW
  for ( int qpk_i=0; qpk_i < qpkVect.size(); qpk_i++ ) {
    // Time difference to identify if current day has gone
    timeDiff = timeVect[qpk_i] - currentDay;
    // Checking for time discontinuity in Qpk series.
    if ( timeDiff > 2*oneDay ) {
      currentDay += oneDay*(timeDiff/oneDay)-oneDay;
      crrDayForMw = 0;
      aveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;
      isDistOk = false;
      resetQpksDays(qpksDay0, qpksDay1, qpksDay2, qpksDay3, qpksDay4,
              qpksDay5, qpksDay6);
      for ( int day_i=0; day_i<daysForMw; day_i++ )
        distQpkTemp[day_i].clear();
      continue;
    }
    // Checking if current day has gone
    if ( timeDiff > oneDay ) {
      if ( aveDay > 0 ) {
        aveDay /= qpkInDay;
        if ( qpkInDay > 1 ) {
          // Using the error of the mean
          // retDistErrMean->Fill( rms/sqrt(qpkInDay) );
          rms = sqrt(qpk2/qpkInDay - aveDay*aveDay);
          // Filling the vector corresponding with current Qpk-average-day
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
            case 6 :
              qpksDay6[0] = aveDay;
              qpksDay6[1] = rms/sqrt(qpkInDay);
              break;
          }
          crrDayForMw++;         
        }
        // Checking if the number of days has gone = time-window
        if ( crrDayForMw == daysForMw ) {
          crrDayForMw = daysForMw - 1;
          // Asking if current 7-days is Ok
          isDistOk = doMovingWindow(qpksDay0, qpksDay1, qpksDay2, qpksDay3, 
              qpksDay4, qpksDay5, qpksDay6, pmt, ifIsUub, currSlope);
          if ( isDistOk ) {
            if ( minSlp > fabs(currSlope) ) {
              minSlp = fabs(currSlope);
              distQpkSelect.clear();
              for ( int day_i=0; day_i<daysForMw; day_i++ )
                for ( auto & qpks : distQpkTemp[day_i] )
                  distQpkSelect.push_back( qpks );
              selectionOk = true;
            }
            qpksDay0 = qpksDay1;
            qpksDay1 = qpksDay2;
            qpksDay2 = qpksDay3;
            qpksDay3 = qpksDay4;
            qpksDay4 = qpksDay5;
            qpksDay5 = qpksDay6;
          }
          else {
            resetQpksDays(qpksDay0, qpksDay1, qpksDay2, qpksDay3, qpksDay4,
                qpksDay5, qpksDay6);
            crrDayForMw = 0;
          }
          // Cleaning vector with current <Qpk>-day
          distQpkTemp[crrDayForMw].clear();
        }
      }
      else {
        resetQpksDays(qpksDay0, qpksDay1, qpksDay2, qpksDay3, qpksDay4,
            qpksDay5, qpksDay6);
        crrDayForMw = 0;
      }
      aveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;
      currentDay += oneDay;
      isDistOk = false;
      // Checking if the number of days has gone = time-window

    }
    // Doing Qpk day average
    if ( !(qpkVect[qpk_i] > 0.) )
      continue;
    // Doing day average
    aveDay += qpkVect[qpk_i];
    qpk2 += qpkVect[qpk_i]*qpkVect[qpk_i];
    qpkInDay++;
    distQpkTemp[crrDayForMw].push_back( qpkVect[qpk_i] );
  }
  // Returning the Qpk-7days distribution with
  // the smaller Slope
  if ( selectionOk ) {
    for ( auto & qpks : distQpkSelect )
      retDistSelec.push_back( qpks );
  }
  
  distQpkTemp.clear();
  distQpkSelect.clear();

  return selectionOk;
}


void makingQpkAccuracy() {

  // output file storing results 
  TFile *outFile = new TFile("makingQpkAccuracy.root","RECREATE");
 
  // Reading stations from list
  ifstream stationSelected;
  TString strFileSelecSt = "../fullUubStationsListVert.txt";
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

  // Vector for store Accuracy
  vector < double > qpkSelUbTemp;
  vector < double > qpkSelUubTemp;
  vector < vector < vector < double > > > qpkSelUb;
  vector < vector < vector < double > > > qpkSelUub;

  qpkSelUb.resize(3);
  qpkSelUub.resize(3);

  // Vector to stores Station IDs.
  vector < double > vecStID;

  bool isSelectionUBOk = false;
  bool isSelectionUUBOk = false;

  // Vector for store station, and PMT, with Ok slope
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
   
  // Applying MW and calculating Accuracy
  int countStations = 0;
  // flag to store station with proper Qpk distribution
  bool stOk = false;
  for ( auto & st_i : stId ) {
    cout << "Moving window for station " << st_i << endl;
    //if ( st_i != 1798 )
      //continue;
    
    // Fetching Qpk and time values
    for ( int pmt_i=1; pmt_i<4; pmt_i++ ) {
      isSelectionUBOk = false;
      isSelectionUUBOk = false;

      fillQpkTimeVals(false, pmt_i, st_i, qpkUb[pmt_i-1], 
          timeUb[pmt_i-1]);
      fillQpkTimeVals(true, pmt_i, st_i, qpkUub[pmt_i-1], 
          timeUub[pmt_i-1]);
      // Applying moving-window and choosing the Qpk-6days
      isSelectionUBOk = doAvePerTime(false, qpkUb[pmt_i-1], 
          timeUb[pmt_i-1], pmt_i, qpkSelUbTemp);
      isSelectionUUBOk = doAvePerTime(true, qpkUub[pmt_i-1], 
          timeUub[pmt_i-1], pmt_i, qpkSelUubTemp);
      // Asking if the selection was Ok for two PMT (ub and uub)
      // and then storing it in qpkSelUb/qpkSelUub
      if ( isSelectionUBOk==true && isSelectionUUBOk==true ) {
        for ( auto qpks : qpkSelUbTemp )
          qpkSelUb[pmt_i-1][countStations].push_back( qpks );        
        for ( auto qpks : qpkSelUubTemp )
          qpkSelUub[pmt_i-1][countStations].push_back( qpks );
        stationOk[pmt_i-1][countStations] = 1;
        stOk = true;
      }
      qpkSelUbTemp.clear();
      qpkSelUubTemp.clear();
    }
    if ( stOk )
      vecStID.push_back( st_i );
    stOk = false;
    countStations++;

    qpkUb.clear();
    qpkUub.clear();
    timeUb.clear();
    timeUub.clear();

    qpkUb.resize(3);
    qpkUub.resize(3);
    timeUb.resize(3);
    timeUub.resize(3);

    //if ( st_i == 830 )
      //break;
  }
  // TH1D to fit Qpk distribution per station 
  int nbins = 200;
  double frsBin = 0.;
  double lstBin = 2.;
  TH1D *uniqRelQpkDistUB = new TH1D("uniqRelQpkDistUB", "", nbins, frsBin, lstBin);
  TH1D *uniqRelQpkDistUUB = new TH1D("uniqRelQpkDistUUB", "", nbins, frsBin, lstBin);

  // TH1D to store the diff between sigma/mu ub and uub per ID
  TH1 *diffAccPerID = new TH1D("diffAcc", "", countStations, 0, countStations);
  double tmpAccUub = 0.;
  double tmpAccUb = 0.;
  vector < double > aveAccUub(2);
  vector < double > aveAccUb(2);
  for (int i=0; i<2; i++) {
    aveAccUb[i] = 0.;
    aveAccUub[i] = 0.;
  }

  // TH1D to store dist. diff between sigma/mu ub and uub
  TH1D *distDiffAcc = new TH1D("distDiffAcc","", 240, -40., 40.);

  // TH1D to store sigma/mu per Station ID
  TH1D *accPerIDub = new TH1D("accPerIDub", "", countStations, 0, countStations);
  TH1D *accPerIDuub = new TH1D("accPerIDuub", "", countStations, 0, countStations);
  double accErrUb = 0.;
  double accErrUub = 0.;
  double accErrMuTerm = 0.;
  double accErrSgTerm = 0.;
  TString strStName;

  // TH1D to plot and calculate the accuracy
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
    // Filling the relative Qpk  
    for ( int pmt_i=0; pmt_i<3; pmt_i++ )
      // Checking if selection was Ok for current PMT
      if ( stationOk[pmt_i][st_i] == 1 ) {       
        fillingDist(pmt_i, st_i, qpkSelUb, uniqRelQpkDistUB);
        fillingDist(pmt_i, st_i, qpkSelUub, uniqRelQpkDistUUB);
        doFit = true;
      }
    // Filling accuracy distribution histo
    if ( doFit ) {
      tmpAccUb = 0.;
      tmpAccUub = 0.;
      strStName.Form("%d", (int)stId[st_i]);
      status = uniqRelQpkDistUB->Fit("gaus","Q","",0,2);
      /*
      TCanvas *c0 = canvasStyle("c0");
      c0->cd();
      uniqRelQpkDistUB->GetXaxis()->SetTitle("#sigma/#mu [%]");
      uniqRelQpkDistUB->GetYaxis()->SetTitle("Counts [au]");
      uniqRelQpkDistUB->Draw();
      TLegend *leg = new TLegend(0.,0.8,0.6,0.9);
      leg->AddEntry(uniqRelQpkDistUB,"St. 1798 UB","");
      leg->SetTextSize(0.05);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->Draw();
      c0->Print("../plots2/accDistSt1798ub.pdf");
      */
      if ( status == 0 ) {
        mu = uniqRelQpkDistUB->GetFunction("gaus")->GetParameter(1);
        sgm = uniqRelQpkDistUB->GetFunction("gaus")->GetParameter(2);
        tmpAccUb = 100.*sgm/mu;
        if ( tmpAccUb < 5.0 )
          distAccUB->Fill( tmpAccUb );

        accPerIDub->SetBinContent(st_i, tmpAccUb);
        if ( tmpAccUb < 5.0 ) {
          aveAccUb[0] += tmpAccUb;
          aveAccUb[1]++;
        }
        accErrMuTerm = (1./mu)*uniqRelQpkDistUB->GetRMSError(1);
        accErrSgTerm = (sgm/(mu*mu))*uniqRelQpkDistUB->GetMeanError(1);
        accErrUb = 100.*( accErrSgTerm*accErrSgTerm + accErrMuTerm*accErrMuTerm );
        accPerIDub->SetBinError(st_i, accErrUb);
        accPerIDub->GetXaxis()->SetBinLabel(st_i, strStName);
      } 
      status = uniqRelQpkDistUUB->Fit("gaus","Q","",0,2);
      /*
      uniqRelQpkDistUUB->GetXaxis()->SetTitle("#sigma/#mu [%]");
      uniqRelQpkDistUUB->GetYaxis()->SetTitle("Counts [au]");
      uniqRelQpkDistUUB->Draw();
      leg = new TLegend(0.05,0.8,0.6,0.9);
      leg->AddEntry(uniqRelQpkDistUB,"St. 1798 UUB","");
      leg->SetTextSize(0.05);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->Draw();
      c0->Print("../plots2/accDistSt1798uub.pdf");
      */
      if ( status == 0 ) {
        mu = uniqRelQpkDistUUB->GetFunction("gaus")->GetParameter(1);
        sgm = uniqRelQpkDistUUB->GetFunction("gaus")->GetParameter(2);
        tmpAccUub = 100.*sgm/mu;
        if ( tmpAccUub < 5.0 )
          distAccUUB->Fill( tmpAccUub );

        accPerIDuub->SetBinContent(st_i, tmpAccUub);
         if ( tmpAccUub < 5.0 ) {
          aveAccUub[0] += tmpAccUub;
          aveAccUub[1]++;
        }
        accErrMuTerm = (1./mu)*uniqRelQpkDistUUB->GetRMSError(1);
        accErrSgTerm = (sgm/(mu*mu))*uniqRelQpkDistUUB->GetMeanError(1);
        accErrUub = 100.*( accErrSgTerm*accErrSgTerm + accErrMuTerm*accErrMuTerm );
        accPerIDuub->SetBinError(st_i, accErrUub);
        accPerIDuub->GetXaxis()->SetBinLabel(st_i, strStName);
      }
      if ( tmpAccUb > 0 && tmpAccUub > 0 ) {
        if ( (tmpAccUb-tmpAccUub) > -2.0 ) { 
          diffAccPerID->SetBinContent(st_i, tmpAccUb-tmpAccUub);
          diffAccPerID->SetBinError(st_i, sqrt( pow(accErrUb,2) + pow(accErrUub,2)));
          diffAccPerID->GetXaxis()->SetBinLabel(st_i, strStName);
          distDiffAcc->Fill( tmpAccUb-tmpAccUub );
        }
        else
          cout << strStName << " " 
            << tmpAccUb << " " 
            << tmpAccUub << " "
            << tmpAccUb-tmpAccUub << endl;
      }
    }
    uniqRelQpkDistUB->Reset();
    uniqRelQpkDistUUB->Reset();
  }

  // Plotting the accuracy results
  TString strMean;
  TString strRms;
  TString strMeanErr;
  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  distAccUB->SetStats(kFALSE);
  distAccUB->SetLineColor(kBlue);
  distAccUB->GetXaxis()->SetRangeUser(0., 5.);
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

  TLegend *leg = new TLegend(0.59,0.6,0.9,0.95);
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

  TCanvas *c2 = canvasStyle("c2");
  TLine *lineAveAccUb;
  TLine *lineAveAccUub;
  c2->cd();

  accPerIDuub->SetTitle("");
  accPerIDuub->SetStats(kFALSE);
  accPerIDuub->GetXaxis()->SetTitle("Station ID.");
  accPerIDuub->GetXaxis()->SetTitleOffset(1.5); 
  //accPerIDuub->GetYaxis()->SetRangeUser(0.2, 2.4);
  accPerIDuub->GetYaxis()->SetTitle("#sigma/#mu [%]");
  accPerIDuub->SetMarkerStyle(71);
  accPerIDuub->SetMarkerColor(kRed);
  accPerIDuub->SetLineColor(kRed);
  accPerIDuub->SetMarkerSize(1.2);
  accPerIDuub->Draw("E1");

  lineAveAccUub = new TLine(0.2, aveAccUub[0]/aveAccUub[1], 74.5, aveAccUub[0]/aveAccUub[1]);
  lineAveAccUub->SetLineWidth(2);
  lineAveAccUub->SetLineColor(kRed);
  lineAveAccUub->Draw();

  accPerIDub->SetMarkerStyle(73);
  accPerIDub->SetMarkerColor(kBlue);
  accPerIDub->SetLineColor(kBlue);
  accPerIDub->SetMarkerSize(1.2);
  accPerIDub->Draw("E1 same");

  lineAveAccUb = new TLine(0.2, aveAccUb[0]/aveAccUb[1], 74.5, aveAccUb[0]/aveAccUb[1]);
  lineAveAccUb->SetLineWidth(2);
  lineAveAccUb->SetLineColor(kBlue);
  lineAveAccUb->Draw();

  leg = new TLegend(0.6,0.7,0.95,0.95);
  strMean.Form("%.2f", aveAccUub[0]/aveAccUub[1]);
  leg->AddEntry(accPerIDuub,"UUB Ave.: "+strMean+" %", "p");
  strMean.Form("%.2f", aveAccUb[0]/aveAccUb[1]);
  leg->AddEntry(accPerIDub,"UB Ave.: "+strMean+" %", "p");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c2->Print("../plots2/accQpkFitUbUubPerSt_Stats.pdf");

  accPerIDub->Write();
  accPerIDuub->Write();
  
  TCanvas *c3 = canvasStyle("c3");
  c3->cd();

  diffAccPerID->SetTitle("");
  diffAccPerID->SetStats(kFALSE);                       
  diffAccPerID->GetXaxis()->SetTitle("Station ID.");
  diffAccPerID->GetXaxis()->SetTitleOffset(1.5);
  diffAccPerID->GetYaxis()->SetTitle("#left(#sigma/#mu#right)_{UB} - #left(#sigma/#mu#right)_{UUB} [%]");
  diffAccPerID->SetMarkerStyle(71);
  diffAccPerID->SetMarkerColor(kGreen+3);
  diffAccPerID->SetLineColor(kGreen+3);
  diffAccPerID->SetMarkerSize(1.2);
  diffAccPerID->Draw("P");
 
  lineAveAccUb = new TLine(0.5,0.,74.5,0.);
  lineAveAccUb->SetLineWidth(2);
  lineAveAccUb->SetLineColor(kGray);
  lineAveAccUb->SetLineStyle(2);
  lineAveAccUb->Draw();

  c3->Print("../plots2/accQpkFitUbUubDiffPerSt_Stats.pdf");

  diffAccPerID->Write();

  TCanvas *c4 = canvasStyle("c4");
  c4->cd();

  distDiffAcc->SetTitle("");
  distDiffAcc->SetStats(kFALSE);
  distDiffAcc->GetXaxis()->SetTitle("#left(#sigma/#mu#right)_{UB} - #left(#sigma/#mu#right)_{UUB} [%]");
  distDiffAcc->GetXaxis()->SetTitleOffset(1.3);
  distDiffAcc->GetXaxis()->SetRangeUser(-6., 6.);
  distDiffAcc->GetYaxis()->SetTitle("Counts [au]");
  distDiffAcc->GetYaxis()->SetRangeUser(-2., 2.);
  strMean.Form("%.3f", distDiffAcc->GetMean());
  strRms.Form("%.3f", distDiffAcc->GetMeanError());
  distDiffAcc->GetYaxis()->SetRangeUser(0., 26.);
  distDiffAcc->SetLineColor(kGreen+3);
  distDiffAcc->Draw();

  leg = new TLegend(0.5,0.78,0.98,0.9);
  leg->AddEntry(distDiffAcc, "MEAN: "+strMean+" % #pm "+strRms+"%", "");
  strRms.Form("%.3f", distDiffAcc->GetRMS());
  leg->AddEntry(distDiffAcc, "RMS:     "+strRms+" %", "");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  lineAveAccUb = new TLine(0.,0.,0.,26.);
  lineAveAccUb->SetLineWidth(2);
  lineAveAccUb->SetLineColor(kGray);
  lineAveAccUb->SetLineStyle(2);
  lineAveAccUb->Draw();
  c4->Print("../plots2/accQpkFitUbUubDistDiff_Stats.pdf");
  distDiffAcc->Write();


  outFile->Write();
  outFile->Close();
  //exit(0);
}
