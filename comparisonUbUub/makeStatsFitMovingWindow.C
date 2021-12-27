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

void plottingAndSaving(TString canvasName, TString outputName, TH2D *histo) {
  TCanvas *canvas = canvasStyle(canvasName); 
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.95);
  gStyle->SetPalette(56);

  histo->GetYaxis()->SetTitle("#chi^{2}");
  histo->GetXaxis()->SetTitle("Slope [FADC/day]");
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetZaxis()->SetTitle("Counts [au]");
  histo->Draw("COLZ");
  canvas->Print("../plots2/"+outputName+".pdf");
}

void plottingProj(TString canvasName, TString outputName, TH2D *histo, bool ifXi) {
  TCanvas *canvas = canvasStyle(canvasName);

  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetRightMargin(0.1245307);
  canvas->SetLeftMargin(0.1170213);
  canvas->SetTopMargin(0.01541096);
  canvas->SetBottomMargin(0.125);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderMode(0);

  //canvas->SetTopMargin(0.);
  //canvas->SetBottomMargin(-0.8);

  //canvas->Divide(1,2);
  //canvas->cd();
  TH1D *hc = new TH1D("hc", "", 200, 0., 100.);
  int maxCounts = 0;
  double sum = 0.;

  if ( ifXi ) {
    canvas->cd();
    histo->ProjectionY()->GetYaxis()->SetTitle("Counts [au]");
    histo->ProjectionY()->GetYaxis()->SetTitleSize(0.06);
    histo->ProjectionY()->GetYaxis()->SetLabelSize(0.05);
    histo->ProjectionY()->GetYaxis()->SetTitleOffset(1.0);
    histo->GetYaxis()->SetTitleSize(0.06);
    histo->GetYaxis()->SetLabelSize(0.05);
    histo->ProjectionY()->Draw();
    canvas->Update();
  
    maxCounts = histo->ProjectionY()->GetEntries();
    for ( int i=0; i<histo->ProjectionY()->GetNbinsX()+1; i++ ) {
      sum += histo->ProjectionY()->GetBinContent(i); ///maxCounts;
      hc->SetBinContent( histo->ProjectionY()->GetBin(i+0.5), sum );
      hc->SetBinError( histo->ProjectionY()->GetBin(i+0.5), 0 );
    }
    hc->Scale(1./sum);
    Float_t rightmax = 1.1*hc->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    hc->Scale(scale);
    hc->SetLineColor(kRed);
    hc->Draw("same"); 

    TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), 
        gPad->GetUxmax(), gPad->GetUymax(), 0., rightmax, 510, "+L");
    axis->SetLabelSize(0.05);
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetTitle("Cumulative [au]");
    axis->SetTitleColor(kRed);
    axis->SetTitleSize(0.05);
    axis->SetTitleOffset(1.0);
    axis->Draw();
    
    /*
    canvas->cd(2);
    hc->SetStats(kFALSE);
    hc->GetYaxis()->SetTitle("[%]");
    hc->GetXaxis()->SetTitle("#chi^{2}");
    hc->GetYaxis()->SetLabelSize(0.08);
    hc->GetYaxis()->SetTitleOffset(0.5);
    hc->Draw();
    .canvas->Update();
    */
  }
  else {
    histo->ProjectionX()->GetYaxis()->SetTitle("Counts [au]");
    histo->ProjectionX()->GetYaxis()->SetTitleSize(0.06);
    histo->ProjectionX()->GetYaxis()->SetTitleOffset(1.0);
    histo->ProjectionX()->GetYaxis()->SetLabelSize(0.05);
    histo->ProjectionX()->GetYaxis()->SetTitleOffset(0.7);
    //histo->ProjectionX()->GetXaxis()->SetTitleSize(.5);
    histo->GetXaxis()->SetTitleSize(0.06);
    histo->GetYaxis()->SetLabelSize(0.05);
    histo->ProjectionX()->Draw();
  }
  canvas->Print("../plots2/"+outputName+".pdf");
}


void fillingQpksFullDays(vector<double> inVals, 
    vector<vector<double>> &retVals) {
  retVals[0].push_back( inVals[0] );
  retVals[1].push_back( inVals[1] );
}

void doMovingWindow(vector<double> qpksDay0, vector<double> qpksDay1,
    vector<double> qpksDay2, vector<double> qpksDay3, vector<double> qpksDay4, 
    vector<double> qpksDay5, vector<double> &retChi2, vector<double> &retSlope) {

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

  distQpks->Fit("pol1","0Q");
  double slp = distQpks->GetFunction("pol1")->GetParameter(1);

  retChi2.push_back( distQpks->GetFunction("pol1")->GetChisquare() );
  retSlope.push_back( slp );
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

void doAvePerTime(bool ifIsUub, vector<double> qpkVect, vector<double> timeVect,
    vector<double> &retChi2, vector<double> &retSlope) {
  int currentDay = (ifIsUub) ? 1311811218 : 1217116818;//1248652818; //1217116818; // August 1st
  const int oneDay = 86400;
  double aveDay = 0.;
  double qpk2 = 0.;
  double rms = 0.;
  int qpkInDay = 0;
  int timeDiff = 0;
  
  const int daysForMw = 6;
  int crrDayForMw = 0;
  vector < double > qpksDay0(2);
  vector < double > qpksDay1(2);
  vector < double > qpksDay2(2);
  vector < double > qpksDay3(2);
  vector < double > qpksDay4(2);
  vector < double > qpksDay5(2);
 
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
      aveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;
      currentDay += oneDay;
      crrDayForMw++;
      if ( crrDayForMw == daysForMw ) {
        crrDayForMw = 5;
        //do MW retChi2 retSlope
        doMovingWindow(qpksDay0, qpksDay1, qpksDay2, qpksDay3, qpksDay4, qpksDay5,
            retChi2, retSlope);
        qpksDay0 = qpksDay1;
        qpksDay1 = qpksDay2;
        qpksDay2 = qpksDay3;
        qpksDay3 = qpksDay4;
        qpksDay4 = qpksDay5;
      }
    }
    if ( qpkVect[qpk_i] > 0. ) {
      aveDay += qpkVect[qpk_i];
      qpk2 += qpkVect[qpk_i]*qpkVect[qpk_i];
      qpkInDay++;
    }
    if ( timeDiff > 2*oneDay ) {
      currentDay += oneDay*(timeDiff/oneDay)-oneDay;
      crrDayForMw = 0;
    }
  }
}


void makeStatsFitMovingWindow() {

  // Getting the working stations
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

  // Vectors for store chi2 and slope
  vector < vector < double > > distChi2Ub;
  vector < vector < double > > distSlopUb;
  vector < vector < double > > distChi2Uub;
  vector < vector < double > > distSlopUub;

  distChi2Ub.resize(3);
  distSlopUb.resize(3);
  distChi2Uub.resize(3);
  distSlopUub.resize(3);

  int nBinsX = 80; // X axis for Slope
  double xLow = -20.;
  double xUp = 20.;
  int nBinsY = 200; // Y axis for Chi2
  double yLow = 0.;
  double yUp = 100;
  TH2D *chi2VsSlopPmt1Ub = new TH2D("chi2VsSlopPmt1Ub", "", 
      nBinsX, xLow, xUp, nBinsY, yLow, yUp);
  TH2D *chi2VsSlopPmt2Ub = new TH2D("chi2VsSlopPmt2Ub", "", 
      nBinsX, xLow, xUp, nBinsY, yLow, yUp);
  TH2D *chi2VsSlopPmt3Ub = new TH2D("chi2VsSlopPmt3Ub", "", 
      nBinsX, xLow, xUp, nBinsY, yLow, yUp);

  TH2D *chi2VsSlopPmt1Uub = new TH2D("chi2VsSlopPmt1Uub", "", 
      nBinsX, xLow, xUp, nBinsY, yLow, yUp);
  TH2D *chi2VsSlopPmt2Uub = new TH2D("chi2VsSlopPmt2Uub", "", 
      nBinsX, xLow, xUp, nBinsY, yLow, yUp);
  TH2D *chi2VsSlopPmt3Uub = new TH2D("chi2VsSlopPmt3Uub", "", 
      nBinsX, xLow, xUp, nBinsY, yLow, yUp);

  for ( auto & st_i : stId ) {
    cout << "Doing for station " << st_i << endl;
    //if ( st_i != 863 )
      //continue;
    // Fetching Qpk and time values
    for ( int pmt_i=1; pmt_i<4; pmt_i++ ) {
      cout << "Doing for PMT " << pmt_i << endl;
      fillQpkTimeVals(false, pmt_i, st_i,
          qpkUb[pmt_i-1], timeUb[pmt_i-1]);
      fillQpkTimeVals(true, pmt_i, st_i,
          qpkUub[pmt_i-1], timeUub[pmt_i-1]);
      
      // Fitting applying moving-window 
      doAvePerTime(false, qpkUb[pmt_i-1], timeUb[pmt_i-1],
          distChi2Ub[pmt_i-1], distSlopUb[pmt_i-1]);
      doAvePerTime(true, qpkUub[pmt_i-1], timeUub[pmt_i-1],
          distChi2Uub[pmt_i-1], distSlopUub[pmt_i-1]);
    }

    qpkUb.clear();
    qpkUb.resize(3);
    timeUb.clear();
    timeUb.resize(3);
    qpkUub.clear();
    qpkUub.resize(3);
    timeUub.clear();
    timeUub.resize(3);

    fillingChi2VsSlop( distChi2Ub[0], distSlopUb[0], chi2VsSlopPmt1Ub );
    fillingChi2VsSlop( distChi2Ub[1], distSlopUb[1], chi2VsSlopPmt2Ub );
    fillingChi2VsSlop( distChi2Ub[2], distSlopUb[2], chi2VsSlopPmt3Ub );

    fillingChi2VsSlop( distChi2Uub[0], distSlopUub[0], chi2VsSlopPmt1Uub );
    fillingChi2VsSlop( distChi2Uub[1], distSlopUub[1], chi2VsSlopPmt2Uub );
    fillingChi2VsSlop( distChi2Uub[2], distSlopUub[2], chi2VsSlopPmt3Uub );

    //if ( st_i == 545 )
      //break;
  }

  // Plotting and saving
  
  plottingAndSaving("Pmt1Ub", "chi2VsSlopPmt1Ub", chi2VsSlopPmt1Ub);
  plottingAndSaving("Pmt2Ub", "chi2VsSlopPmt2Ub", chi2VsSlopPmt2Ub);
  plottingAndSaving("Pmt3Ub", "chi2VsSlopPmt3Ub", chi2VsSlopPmt3Ub);

  plottingAndSaving("Pmt1Uub", "chi2VsSlopPmt1Uub", chi2VsSlopPmt1Uub);
  plottingAndSaving("Pmt2Uub", "chi2VsSlopPmt2Uub", chi2VsSlopPmt2Uub);
  plottingAndSaving("Pmt3Uub", "chi2VsSlopPmt3Uub", chi2VsSlopPmt3Uub);
  

  plottingProj("SlopDistUB", "chi2VsSlopPmt1UbProjSlop", chi2VsSlopPmt1Ub, 0);
  plottingProj("SlopDistUB", "chi2VsSlopPmt2UbProjSlop", chi2VsSlopPmt2Ub, 0);
  plottingProj("SlopDistUB", "chi2VsSlopPmt3UbProjSlop", chi2VsSlopPmt3Ub, 0);

  plottingProj("Chi2DistUB", "chi2VsSlopPmt1UbProjChi2", chi2VsSlopPmt1Ub, 1);
  plottingProj("Chi2DistUB", "chi2VsSlopPmt2UbProjChi2", chi2VsSlopPmt2Ub, 1);
  plottingProj("Chi2DistUB", "chi2VsSlopPmt3UbProjChi2", chi2VsSlopPmt3Ub, 1);

  plottingProj("SlopDistUUB", "chi2VsSlopPmt1UubProjSlop", chi2VsSlopPmt1Uub, 0);
  plottingProj("SlopDistUUB", "chi2VsSlopPmt2UubProjSlop", chi2VsSlopPmt2Uub, 0);
  plottingProj("SlopDistUUB", "chi2VsSlopPmt3UubProjSlop", chi2VsSlopPmt3Uub, 0);

  plottingProj("Chi2DistUUB", "chi2VsSlopPmt1UubProjChi2", chi2VsSlopPmt1Uub, 1);
  plottingProj("Chi2DistUUB", "chi2VsSlopPmt2UubProjChi2", chi2VsSlopPmt2Uub, 1);
  plottingProj("Chi2DistUUB", "chi2VsSlopPmt3UubProjChi2", chi2VsSlopPmt3Uub, 1);

  //exit(0);
}
