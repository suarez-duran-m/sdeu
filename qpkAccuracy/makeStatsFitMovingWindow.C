TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 102, 76, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12); 
  canvas->SetRightMargin(0.14);
  canvas->SetTopMargin(0.015); 
  canvas->SetFrameBorderMode(0);
  return canvas;
}

void plottingAndSaving(TString canvasName, TString outputName, TH2D *histo) {
  TCanvas *canvas = canvasStyle(canvasName); 
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.95);
  gStyle->SetPalette(56);

  histo->GetYaxis()->SetTitle("#chi^{2}");
  histo->GetXaxis()->SetTitle("Slope [FADC/day/#LTFADC#GT_{7days}]");
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetZaxis()->SetTitle("Counts [au]");
  histo->Draw("COLZ");
  canvas->Print("../plots2/"+outputName+".pdf");
  canvas->Close();
}


void plottingAndSaving(bool isLog, TString canvasName, TString outputName, 
    TH1D *histo) {
  TCanvas *canvas = canvasStyle(canvasName); 
  histo->GetYaxis()->SetTitle("Counts [au]");
  if ( !isLog ) {
    histo->GetXaxis()->SetTitle("P-Value [au]");
    histo->GetYaxis()->SetRangeUser(0, 1e3);
  }
  else
    histo->GetXaxis()->SetTitle("Log10(Pval) [au]");
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->Draw();
  canvas->Print("../plots2/"+outputName+".pdf");
}


void plottingProj(TString canvasName, TString outputName, TH2D *histo, bool ifXi) {
  TCanvas *canvas = canvasStyle(canvasName);

  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.1270213);
  canvas->SetTopMargin(0.01541096);
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderMode(0);

  //TH1D *hc = new TH1D("hc", "", 200, 0., 100.);
  //int maxCounts = 0;
  //double sum = 0.;

  if ( ifXi ) {
    canvas->cd();
    histo->ProjectionY()->GetYaxis()->SetTitle("Counts [au]");
    histo->ProjectionY()->GetYaxis()->SetTitleSize(0.06);
    histo->ProjectionY()->GetYaxis()->SetLabelSize(0.05);
    histo->ProjectionY()->GetYaxis()->SetTitleOffset(1.0);
    histo->GetYaxis()->SetTitleSize(0.06);
    histo->GetYaxis()->SetLabelSize(0.05);
    histo->Scale( 1./histo->ProjectionY()->GetEntries() );
    histo->ProjectionY()->Draw();
    canvas->Update();
    /*
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
    */
    
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
    TH1D *proj_x = histo->ProjectionX();
    proj_x->Fit("gaus","Q");
    proj_x->GetYaxis()->SetTitle("Counts [au]");
    proj_x->GetYaxis()->SetTitleSize(0.06);
    proj_x->GetYaxis()->SetTitleOffset(1.1);
    proj_x->GetYaxis()->SetLabelSize(0.05);

    proj_x->GetXaxis()->SetTitleSize(0.06);
    proj_x->GetXaxis()->SetLabelSize(0.05);
    proj_x->GetXaxis()->SetTitleOffset(1.1);
    proj_x->GetXaxis()->SetRangeUser(-0.06, 0.06);
    proj_x->Draw();
    
    TString strVal;
    TLegend *leg = new TLegend(0.641, 0.52, 0.85, 0.7);
    leg->AddEntry(proj_x->GetFunction("gaus"), "Gaus fitted parameters:", "l");
    strVal.Form("%.2e", proj_x->GetFunction("gaus")->GetParameter(1));
    leg->AddEntry(histo, "#mu = "+strVal,"");
    strVal.Form("%.2e", proj_x->GetFunction("gaus")->GetParameter(2));
    leg->AddEntry(histo, "#sigma = "+strVal,"");    
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
  }
  canvas->Print("../plots2/"+outputName+".pdf");
  canvas->Close();
}

void plottingCrossCheck(TString stStr, TString ifUubStr, TGraphErrors *qpk, 
    TGraphErrors *slope, TGraph *chi2) {  
  TCanvas *canvCross = new TCanvas("canvCross", ifUubStr);
  canvCross->SetTopMargin(0.03125);
  canvCross->SetBottomMargin(0.03333334);

  canvCross->cd();
  
  TPad *canvCrossQpks = new TPad("canvCrossQpks", "", 0,0.66,0.97,0.97);
  canvCrossQpks->Draw();
  canvCrossQpks->cd();
  canvCrossQpks->SetLeftMargin(0.07662624);
  canvCrossQpks->SetRightMargin(0.002205072);
  canvCrossQpks->SetTopMargin(0.003236246);
  canvCrossQpks->SetBottomMargin(0.05501618);

  qpk->SetTitle("");
  qpk->GetYaxis()->SetTitle("#LT Qpk #GT_{1-Day}");
  qpk->GetXaxis()->SetTimeFormat("%y/%m/%d");
  qpk->GetXaxis()->SetTimeOffset(315964782,"gmt");
  if ( ifUubStr == "UB" )
    qpk->GetXaxis()->SetRange(5, 100);
  else
    qpk->GetXaxis()->SetRange(6, 122);
  qpk->SetMarkerStyle(70);
  qpk->SetMarkerSize(.5);
  qpk->SetMarkerColor(kBlue);
  qpk->SetLineColor(kBlue);
  qpk->GetYaxis()->SetLabelSize(0.05);
  qpk->GetYaxis()->SetTitleSize(0.08);
  qpk->GetYaxis()->SetTitleOffset(0.45);
  qpk->GetYaxis()->SetTitleFont(42);
  qpk->GetYaxis()->SetLabelSize(0.07);
  qpk->GetXaxis()->SetTickLength(0.06);
  qpk->Draw("AP");

  TLegend *legend = new TLegend(0.85,0.75,0.98,0.98);
  legend->AddEntry(qpk, stStr+" "+ifUubStr, "");
  legend->SetTextSize(0.1);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw();

  canvCrossQpks->Modified();

  canvCross->cd();
  TPad *canvCrossSlop = new TPad("canvCrossSlop", "", 0,0.34,0.97,0.66);
  canvCrossSlop->Draw();
  canvCrossSlop->cd();
  canvCrossSlop->SetLeftMargin(0.07662624);
  canvCrossSlop->SetRightMargin(0.002205072);
  canvCrossSlop->SetTopMargin(0.003236246);
  canvCrossSlop->SetBottomMargin(0.05501618);

  slope->SetTitle("");
  slope->GetYaxis()->SetTitle("Slope [FADC/day/#LTFADC#GT_{7Day}]");
  slope->GetXaxis()->SetTimeFormat("%y/%m/%d");
  slope->GetXaxis()->SetTimeOffset(315964782,"gmt");
  slope->SetMarkerStyle(70);
  slope->SetMarkerSize(.5);
  slope->SetMarkerColor(kBlue);
  slope->SetLineColor(kBlue);
  slope->GetYaxis()->SetLabelSize(0.05);
  slope->GetYaxis()->SetTitleSize(.06);
  slope->GetYaxis()->SetTitleOffset(0.55);
  slope->GetXaxis()->SetTickLength(0.06);
  slope->GetYaxis()->SetLabelSize(0.05);
  slope->GetYaxis()->SetTitleSize(0.07);
  slope->GetYaxis()->SetTitleFont(42);
  slope->GetYaxis()->SetLabelSize(0.07);
  slope->Draw("AP");

  canvCross->cd();
  TPad *canvCrossChi2 = new TPad("canvCrossChi2", "", 0,0.0,0.97,0.34);
  canvCrossChi2->Draw();
  canvCrossChi2->cd();
  canvCrossChi2->SetLeftMargin(0.07662624);
  canvCrossChi2->SetRightMargin(0.002205072);
  canvCrossChi2->SetTopMargin(0.003236246);
  canvCrossChi2->SetBottomMargin(1.05501618);

  chi2->SetTitle("");
  chi2->GetYaxis()->SetRangeUser(-10,100);
  chi2->GetYaxis()->SetTitle("#chi^{2}");
  chi2->GetXaxis()->SetTimeFormat("%y/%m/%d");
  chi2->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chi2->SetMarkerStyle(70);
  chi2->SetMarkerSize(.5);
  chi2->SetMarkerColor(kBlue);
  chi2->GetYaxis()->SetRangeUser(-10, 98);
  chi2->GetYaxis()->SetTitleSize(.08);
  chi2->GetYaxis()->SetTitleOffset(.4);
  chi2->GetYaxis()->SetTitleFont(42);
  chi2->GetYaxis()->SetLabelSize(0.07);

  chi2->GetXaxis()->SetLabelSize(0.08);
  chi2->GetXaxis()->SetTickLength(0.06);
  chi2->Draw("AP");

  canvCrossChi2->Modified();   
  canvCross->cd();
  canvCross->Modified();
  canvCross->cd();
  canvCross->SetSelected(canvCross);
  canvCross->Print("../plots2/crossCheck"+stStr+ifUubStr+".pdf");
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

void doMovingWindow(vector<double> qpksDay0, vector<double> qpksDay1,
    vector<double> qpksDay2, vector<double> qpksDay3, vector<double> qpksDay4, 
    vector<double> qpksDay5, vector<double> qpksDay6, vector<double> &retChi2, 
    vector<double> &retPval, vector<double> &retSlope, vector<double> &retSlopErr) {
  // Vector to store <Qpk>-day time series and its respective error
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
    return;

  // TGraph to apply fit and get Slope and Chi2
  vector < double > xAxis{1., 2., 3., 4., 5., 6., 7.};
  TGraphErrors *distQpks = new TGraphErrors(
      qpksFullDays[0].size(), &xAxis[0], &qpksFullDays[0][0], 0, &qpksFullDays[1][0]);

  // Applying linear fit to <Qpk>-per-day during the time-window
  distQpks->Fit("pol1","0Q");
  if ( distQpks->GetFunction("pol1")->GetNDF() != 5 )
    cout << "Wrong NDF" << endl;
 
  double slp = distQpks->GetFunction("pol1")->GetParameter(1);
  double chi2 = distQpks->GetFunction("pol1")->GetChisquare();

  // Storing Chi2 and Slope
  retChi2.push_back( chi2 );
  retPval.push_back( TMath::Prob(chi2, 5) );

  // Normalizing Slope by <qpk> during time-window 
  double meanTimeWidow = 0.;
  for ( int i=0; i<qpksFullDays.size(); i++ )
    meanTimeWidow += qpksFullDays[0][i];
  meanTimeWidow /= qpksFullDays.size();
  retSlope.push_back( slp / meanTimeWidow );
  retSlopErr.push_back( distQpks->GetFunction("pol1")->GetParError(1) / meanTimeWidow );
  qpksFullDays.clear();
}


void fillingChi2VsSlop( vector<double>chi2, vector<double>slp, TH2D *hist ) {
  for ( int slp_i=0; slp_i<slp.size(); slp_i++ )
    hist->Fill( slp[slp_i], chi2[slp_i] );
}


void fillingPval( vector<double>pVal, TH1D *hist, TH1D *histlog ) {
  for ( int pVal_i=0; pVal_i<pVal.size(); pVal_i++ ) {
    hist->Fill( pVal[pVal_i] );
    histlog->Fill( log10(pVal[pVal_i]) );
  }
}


void fillQpkTimeVals(bool ifIsUub, int pmt, int st_id, 
    vector<double> &retQpkVect, vector<double> &retTimeVect) {
  // Path to fitting-root-files
  TString bnCdas = (ifIsUub) ?
    "~/postdoc/sdeu/underHisto/results/uubChPkPMT" :
    "~/postdoc/sdeu/nouub/underHistos/results/ubChPkPMT"; 

  // String for fitting-root-files names
  TString pmtId;
  TString strStId;
  TString fname;
  pmtId.Form("%d", pmt);
  strStId.Form("St%d", (int)st_id);

  // Strings for months and years to read
  TString monthUub[4] = {"Aug", "Sep", "Oct", "Nov"};
  int nMonths = sizeof(monthUub)/sizeof(*monthUub);

  int nYears = 5;
  TString strYear[5] = {"2016", "2018", "2019", "2020", "2021"};
  int stYear = (ifIsUub) ? nYears-1 : 1;
  int lstYear = (ifIsUub) ? nYears-1 : 1;
  
  // Constants to read branch with Qpk values
  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int fetchTime = 0;

  // Reading and fetching fitting-root-files
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
    vector<double> &retChi2, vector<double> &retPval, vector<double> &retSlope, 
    vector<double> &retSlopErr, vector<double> &retTime, vector<double> &retQpkDays,
    vector<double> &retQpkErrDays, vector<double> &retTimeQpk) {
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
  // vectors to storage the average and error-of-mean of Qpk per day.
  vector < double > qpksDay0(2);
  vector < double > qpksDay1(2);
  vector < double > qpksDay2(2);
  vector < double > qpksDay3(2);
  vector < double > qpksDay4(2);
  vector < double > qpksDay5(2);
  vector < double > qpksDay6(2);

  // Doing Qpk average per day
  for ( int qpk_i=0; qpk_i < qpkVect.size(); qpk_i++ ) {
    // Time difference to identify if current day has gone
    timeDiff = timeVect[qpk_i] - currentDay;

    // Checking for time discontinuity in Qpk series.
    if ( timeDiff > 2*oneDay ) {
      currentDay += oneDay*(timeDiff/oneDay)-oneDay;
      aveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;
      crrDayForMw = 0;
      for (int j=0; j<2; j++ ) {
        qpksDay0[j] = 0.;
        qpksDay1[j] = 0.;
        qpksDay2[j] = 0.;
        qpksDay3[j] = 0.;
        qpksDay4[j] = 0.;
        qpksDay5[j] = 0.;
      }
      continue;
    }
    // Checking if current day has gone
    if ( timeDiff > oneDay ) {
      if ( aveDay > 0 ) {
        aveDay /= qpkInDay;
        if ( qpkInDay > 1 ) {
          rms = sqrt(qpk2/qpkInDay - aveDay*aveDay);
          retQpkDays.push_back(aveDay);
          retQpkErrDays.push_back(rms/sqrt(qpkInDay));
          retTimeQpk.push_back(currentDay);
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
        }
      }
      aveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;
      currentDay += oneDay;
      crrDayForMw++;
      // Checking if the number of days has gone = time-window
      if ( crrDayForMw == daysForMw ) {
        crrDayForMw = daysForMw - 1;
        // Doing MW and returning Chi2 and Slope
        doMovingWindow(qpksDay0, qpksDay1, qpksDay2, qpksDay3, qpksDay4, qpksDay5,
            qpksDay6, retChi2, retPval, retSlope, retSlopErr);
        qpksDay0 = qpksDay1;
        qpksDay1 = qpksDay2;
        qpksDay2 = qpksDay3;
        qpksDay3 = qpksDay4;
        qpksDay4 = qpksDay5;
        qpksDay5 = qpksDay6;
        retTime.push_back(currentDay);
      }
    }
    // Doing day average
    if ( qpkVect[qpk_i] > 0. ) {
      aveDay += qpkVect[qpk_i];
      qpk2 += qpkVect[qpk_i]*qpkVect[qpk_i];
      qpkInDay++;
    }
  }
}


void makeStatsFitMovingWindow() {
  // Reading working stations' ID
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

  // Vectors for store chi2 and slope
  vector < vector < double > > distChi2Ub;
  vector < vector < double > > distPvalUb;
  vector < vector < double > > distSlopUb;
  vector < vector < double > > distSlopErrUb;
  vector < vector < double > > qpkDayUb;
  vector < vector < double > > qpkErrDayUb;
  vector < vector < double > > qpkTimeUb;
  vector < vector < double > > timeChi2SlopeUb;

  vector < vector < double > > distChi2Uub;
  vector < vector < double > > distPvalUub;
  vector < vector < double > > distSlopUub;
  vector < vector < double > > distSlopErrUub;
  vector < vector < double > > qpkDayUub;
  vector < vector < double > > qpkErrDayUub;
  vector < vector < double > > qpkTimeUub;
  vector < vector < double > > timeChi2SlopeUub;

  distChi2Ub.resize(3);
  distPvalUb.resize(3);
  distSlopUb.resize(3);
  distSlopErrUb.resize(3);
  qpkDayUb.resize(3);
  qpkErrDayUb.resize(3);
  qpkTimeUb.resize(3);
  timeChi2SlopeUb.resize(3);

  distChi2Uub.resize(3);
  distPvalUub.resize(3);
  distSlopUub.resize(3);
  distSlopErrUub.resize(3);
  qpkDayUub.resize(3);
  qpkErrDayUub.resize(3);
  qpkTimeUub.resize(3);
  timeChi2SlopeUub.resize(3);

  // Constant to initialize TH2 to plot chi2 and slope
  int nBinsX = 200; // bins for Slope's X axis
  double xLow = -.1;
  double xUp = .1;
  int nBinsY = 200; // bins for Chi2's Y axis
  double yLow = 0.;
  double yUp = 100;

  // TH2 to plot chi2 and slope
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

  nBinsX = 200;
  xLow = 0;
  xUp = 1;

  TH1D *distPvalPmt1Ub = new TH1D("distPvalPmt1Ub", "", nBinsX, xLow, xUp);
  TH1D *distPvalPmt2Ub = new TH1D("distPvalPmt2Ub", "", nBinsX, xLow, xUp);
  TH1D *distPvalPmt3Ub = new TH1D("distPvalPmt3Ub", "", nBinsX, xLow, xUp);

  TH1D *distPvalPmt1Uub = new TH1D("distPvalPmt1Uub", "", nBinsX, xLow, xUp);
  TH1D *distPvalPmt2Uub = new TH1D("distPvalPmt2Uub", "", nBinsX, xLow, xUp);
  TH1D *distPvalPmt3Uub = new TH1D("distPvalPmt3Uub", "", nBinsX, xLow, xUp);

  nBinsX = 202;
  xLow = -200;
  xUp = 2;

  TH1D *distLogPvalPmt1Ub = new TH1D("distLogPvalPmt1Ub", "", nBinsX, xLow, xUp);
  TH1D *distLogPvalPmt2Ub = new TH1D("distLogPvalPmt2Ub", "", nBinsX, xLow, xUp);
  TH1D *distLogPvalPmt3Ub = new TH1D("distLogPvalPmt3Ub", "", nBinsX, xLow, xUp);

  TH1D *distLogPvalPmt1Uub = new TH1D("distLogPvalPmt1Uub", "", nBinsX, xLow, xUp);
  TH1D *distLogPvalPmt2Uub = new TH1D("distLogPvalPmt2Uub", "", nBinsX, xLow, xUp);
  TH1D *distLogPvalPmt3Uub = new TH1D("distLogPvalPmt3Uub", "", nBinsX, xLow, xUp);

  // Applying moving-window
  for ( auto & st_i : stId ) {
    //int chosenSt = 864;
    //if ( st_i != chosenSt )
      //continue;
    cout << "Doing for station " << st_i << endl;
    // Fetching Qpk and time values from fitting-root-files
    for ( int pmt_i=1; pmt_i<4; pmt_i++ ) {
      cout << "Doing for PMT " << pmt_i << endl;
      fillQpkTimeVals(false, pmt_i, st_i,
          qpkUb[pmt_i-1], timeUb[pmt_i-1]);
      fillQpkTimeVals(true, pmt_i, st_i,
          qpkUub[pmt_i-1], timeUub[pmt_i-1]);
      
      // Applying moving-window and returning Chi2 and Slope
      doAvePerTime(false, qpkUb[pmt_i-1], timeUb[pmt_i-1],
          distChi2Ub[pmt_i-1], distPvalUb[pmt_i-1], distSlopUb[pmt_i-1], 
          distSlopErrUb[pmt_i-1], timeChi2SlopeUb[pmt_i-1], qpkDayUb[pmt_i-1], 
          qpkErrDayUb[pmt_i-1], qpkTimeUb[pmt_i-1]);
      doAvePerTime(true, qpkUub[pmt_i-1], timeUub[pmt_i-1],
          distChi2Uub[pmt_i-1], distPvalUub[pmt_i-1], distSlopUub[pmt_i-1], 
          distSlopErrUub[pmt_i-1], timeChi2SlopeUub[pmt_i-1], qpkDayUub[pmt_i-1], 
          qpkErrDayUub[pmt_i-1], qpkTimeUub[pmt_i-1]);
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

    fillingPval( distPvalUb[0], distPvalPmt1Ub, distLogPvalPmt1Ub );
    fillingPval( distPvalUb[1], distPvalPmt2Ub, distLogPvalPmt2Ub );
    fillingPval( distPvalUb[2], distPvalPmt3Ub, distLogPvalPmt3Ub );

    fillingPval( distPvalUub[0], distPvalPmt1Uub, distLogPvalPmt1Uub );
    fillingPval( distPvalUub[1], distPvalPmt2Uub, distLogPvalPmt2Uub );
    fillingPval( distPvalUub[2], distPvalPmt3Uub, distLogPvalPmt3Uub );

    // TH1 to plot Slope and Chi2 as function of time
    /*
    TGraphErrors *qpkVstimeUb = new TGraphErrors(qpkTimeUb[0].size(),
        &qpkTimeUb[0][0], &qpkDayUb[0][0], 0, &qpkErrDayUb[0][0]);
    TGraphErrors *slopVstimePmt1Ub = new TGraphErrors(timeChi2SlopeUb[0].size(), 
        &timeChi2SlopeUb[0][0], &distSlopUb[0][0], 0, &distSlopErrUb[0][0]);
    TGraph *chi2VstimePmt1Ub = new TGraph(timeChi2SlopeUb[0].size(), 
        &timeChi2SlopeUb[0][0], &distChi2Ub[0][0]);

    TGraphErrors *qpkVstimeUub = new TGraphErrors(qpkTimeUub[0].size(),
        &qpkTimeUub[0][0], &qpkDayUub[0][0], 0, &qpkErrDayUub[0][0]);
    TGraphErrors *slopVstimePmt1Uub = new TGraphErrors(timeChi2SlopeUub[0].size(), 
        &timeChi2SlopeUub[0][0], &distSlopUub[0][0], 0, &distSlopErrUub[0][0]);
    TGraph *chi2VstimePmt1Uub = new TGraph(timeChi2SlopeUub[0].size(), 
        &timeChi2SlopeUub[0][0], &distChi2Uub[0][0]);
    
    plottingCrossCheck("846", "UB", qpkVstimeUb, slopVstimePmt1Ub, chi2VstimePmt1Ub);
    plottingCrossCheck("846", "UUB", qpkVstimeUub, slopVstimePmt1Uub, 
        chi2VstimePmt1Uub);
    */
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

  bool plotLogPval = false;
  plottingAndSaving(plotLogPval, "Pmt1Ub", "pValDistPmt1Ub", distPvalPmt1Ub);
  plottingAndSaving(plotLogPval, "Pmt2Ub", "pValDistPmt2Ub", distPvalPmt2Ub);
  plottingAndSaving(plotLogPval, "Pmt3Ub", "pValDistPmt3Ub", distPvalPmt3Ub);

  plottingAndSaving(plotLogPval, "Pmt1Uub", "pValDistPmt1Uub", distPvalPmt1Uub);
  plottingAndSaving(plotLogPval, "Pmt2Uub", "pValDistPmt2Uub", distPvalPmt2Uub);
  plottingAndSaving(plotLogPval, "Pmt3Uub", "pValDistPmt3Uub", distPvalPmt3Uub);

  plotLogPval = true;
  plottingAndSaving(plotLogPval, "Pmt1Ub", "logPvalDistPmt1Ub", distLogPvalPmt1Ub);
  plottingAndSaving(plotLogPval, "Pmt2Ub", "logPvalDistPmt2Ub", distLogPvalPmt2Ub);
  plottingAndSaving(plotLogPval, "Pmt3Ub", "logPvalDistPmt3Ub", distLogPvalPmt3Ub);

  plottingAndSaving(plotLogPval, "Pmt1Uub", "logPvalDistPmt1Uub", distLogPvalPmt1Uub);
  plottingAndSaving(plotLogPval, "Pmt2Uub", "logPvalDistPmt2Uub", distLogPvalPmt2Uub);
  plottingAndSaving(plotLogPval, "Pmt3Uub", "logPvalDistPmt3Uub", distLogPvalPmt3Uub);

  //exit(0);
}
