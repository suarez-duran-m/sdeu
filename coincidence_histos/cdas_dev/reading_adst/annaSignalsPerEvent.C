void setCanvasStyle(TCanvas& canvas);
void fillSigVec(vector<vector<double>> signal, vector<vector<int>> evtIds,
    vector<vector<vector<double>>> &slctSigPmt, vector<double> &slctSigTot,
    vector<int> &slctEvtId);
void plottingSigPmt(TGraphErrors grp1, TGraphErrors grp2, TString pmtName, 
    TString pmtId);

void annaSignalsPerEvent() {
  ifstream dataVEM("signalsVEM.dat");
  ifstream dataVEMCoinci("signalsVEMCoinci.dat");

  int tmpEvtId = 0;
  int tmpSt = 0;
  int tmpPmt = 0;
  double tmpSignal = 0.;

  // To storage signals 
  // signals[pmt_i][signal_i]
  vector < vector < double > > signalsVEM(3);
  vector < vector < double > > signalsVEMCoinci(3);
  vector < int > stIdsVEM;
  vector < int > stIdsVEMCoinci;
  vector < vector < int > > evtIdsVEM(3);
  vector < vector < int > > evtIdsVEMCoinci(3);

  // Reading data from ASCII files
  //
  while ( dataVEM >> tmpEvtId >> tmpSt >> tmpPmt >> tmpSignal ) {
    evtIdsVEM[tmpPmt-1].push_back( tmpEvtId );
    stIdsVEM.push_back( tmpSt );
    signalsVEM[tmpPmt-1].push_back( tmpSignal );
  }
  dataVEM.close();

  while ( dataVEMCoinci >> tmpEvtId >> tmpSt >> tmpPmt >> tmpSignal ) {
    evtIdsVEMCoinci[tmpPmt-1].push_back( tmpEvtId );
    stIdsVEMCoinci.push_back( tmpSt );
    signalsVEMCoinci[tmpPmt-1].push_back( tmpSignal );
  }
  dataVEMCoinci.close();

  // Charging data into vectors
  // sig[evt][pmt] ; signal per PMT
  // sigTot[evt] ; totals signal

  vector < vector < vector < double > > > slctSigPmt(3);
  vector < double > slctSigTot;
  vector < int > slctEvtId;

  vector < vector < vector < double > > > slctSigCoinciPmt(3);
  vector < double > slctSigCoinciTot;
  vector < int > slctEvtCoinciId;

  for (int i=0; i<3; i++ ) {
    slctSigPmt[i].resize(1100);
    slctSigCoinciPmt[i].resize(1100);
  }
  
  fillSigVec(signalsVEM, evtIdsVEM, slctSigPmt, slctSigTot, slctEvtId);
  cout << "MSD doing for coincidence" << endl;
  fillSigVec(signalsVEMCoinci, evtIdsVEMCoinci, slctSigCoinciPmt, 
      slctSigCoinciTot, slctEvtCoinciId);

  // Selecting signal for same events
  //
  // sigPmt[pmt][binEner]
  vector < vector < double > > sigPmt(3);
  vector < vector < double > > sigErrPmt(3);
  vector < double > totCntEnerBin;

  vector < vector < double > > sigCoinciPmt(3);
  vector < vector < double > > sigCoinciErrPmt(3);
  vector < double > totCntEnerCoinciBin;

  vector < double > binTotEner;

  // Doing histograms (TGraph) with n-bins per decade
  // bin = int(nbinsDecade*log(x))
  // TotBins = int( maxVal*nbinsDecade )
  int nbinsDecade = 4;
  int totBins = int( log10(1e4)*nbinsDecade );

  for ( int i=0; i<totBins; i++ ) {
    totCntEnerBin.push_back( 0 );
    totCntEnerCoinciBin.push_back( 0 );
    binTotEner.push_back( pow(10., i/(1.*nbinsDecade)) );
    for ( int j=0; j<3; j++ ) {
      sigPmt[j].push_back( 0 );
      sigErrPmt[j].push_back( 0 );
      sigCoinciPmt[j].push_back( 0 );
      sigCoinciErrPmt[j].push_back( 0 );
    }
  }

  for ( int pmt_i = 0; pmt_i<3; pmt_i++ ) {
    for ( int evt_i = 0; evt_i<slctEvtId.size(); evt_i++ ) {
      for ( int evtCo_i = 0; evtCo_i<slctEvtCoinciId.size(); evtCo_i++ ) {
        if ( !(slctEvtCoinciId[evtCo_i] == slctEvtId[evt_i]) )
          continue;
        for ( int sig_i=0; sig_i<slctSigPmt[pmt_i][evt_i].size(); sig_i++ ) {
          if ( slctSigTot[evt_i] < 1 )
            continue;
          sigPmt[pmt_i][ int( nbinsDecade*log10(slctSigTot[evt_i]) ) ]
            += slctSigPmt[pmt_i][evt_i][sig_i];
          sigErrPmt[pmt_i][ int( nbinsDecade*log10(slctSigTot[evt_i]) ) ]
            += slctSigPmt[pmt_i][evt_i][sig_i] * slctSigPmt[pmt_i][evt_i][sig_i];
          totCntEnerBin[ int( nbinsDecade*log10(slctSigTot[evt_i]) ) ]++;
        }
        for ( int sig_i=0; sig_i<slctSigCoinciPmt[pmt_i][evtCo_i].size(); sig_i++ ) {
          if ( slctSigCoinciTot[evtCo_i] < 1 )
            continue;
          sigCoinciPmt[pmt_i][ int( nbinsDecade*log10(slctSigCoinciTot[evtCo_i]) ) ] 
            += slctSigCoinciPmt[pmt_i][evtCo_i][sig_i];
          sigCoinciErrPmt[pmt_i][ int( nbinsDecade*log10(slctSigCoinciTot[evtCo_i]) ) ]
            += slctSigCoinciPmt[pmt_i][evtCo_i][sig_i] 
            * slctSigCoinciPmt[pmt_i][evtCo_i][sig_i];
          totCntEnerCoinciBin[int( nbinsDecade*log10(slctSigCoinciTot[evtCo_i]) )]++;
        }
      }
    }
    double tmp = 0.;
    for ( int i=0; i<sigPmt[pmt_i].size(); i++ ) {
      if ( totCntEnerBin[i] < 1 ) {
        sigPmt[pmt_i][i] = 12345;
        continue;
      }
      sigPmt[pmt_i][i] /= totCntEnerBin[i];
      tmp = sigErrPmt[pmt_i][i] / totCntEnerBin[i] - sigPmt[pmt_i][i]*sigPmt[pmt_i][i];
      tmp = sqrt(tmp);
      sigErrPmt[pmt_i][i] = tmp / sqrt(totCntEnerBin[i]);
    }    
    for ( int i=0; i<sigCoinciPmt[pmt_i].size(); i++ ) {
      if ( totCntEnerCoinciBin[i] < 1 ) {
        sigCoinciPmt[pmt_i][i] = 12345;
        continue;
      }
      sigCoinciPmt[pmt_i][i] /= totCntEnerCoinciBin[i];
      tmp = sigCoinciErrPmt[pmt_i][i] / totCntEnerCoinciBin[i] 
        - sigCoinciPmt[pmt_i][i]*sigCoinciPmt[pmt_i][i];
      tmp = sqrt(tmp);
      sigCoinciErrPmt[pmt_i][i] = tmp / sqrt(totCntEnerCoinciBin[i]);
    } 
  }

  // Plotting
  //
  TGraphErrors grphPmt1( sigPmt[0].size(), &binTotEner.front(), &sigPmt[0].front(), 0, 
       &sigErrPmt[0].front() );
  TGraphErrors grphCoinciPmt1( sigCoinciPmt[0].size(), &binTotEner.front(), 
      &sigCoinciPmt[0].front(), 0, &sigCoinciErrPmt[0].front() );

  TGraphErrors grphPmt2( sigPmt[1].size(), &binTotEner.front(), &sigPmt[1].front(), 0, 
       &sigErrPmt[1].front() );
  TGraphErrors grphCoinciPmt2( sigCoinciPmt[1].size(), &binTotEner.front(), 
      &sigCoinciPmt[1].front(), 0, &sigCoinciErrPmt[1].front() );

  TGraphErrors grphPmt3( sigPmt[2].size(), &binTotEner.front(), &sigPmt[2].front(), 0, 
       &sigErrPmt[2].front() );
  TGraphErrors grphCoinciPmt3( sigCoinciPmt[2].size(), &binTotEner.front(), 
      &sigCoinciPmt[2].front(), 0, &sigCoinciErrPmt[2].front() );

  plottingSigPmt(grphPmt1, grphCoinciPmt1, "S_{1}", "1");
  plottingSigPmt(grphPmt2, grphCoinciPmt2, "S_{2}", "2");
  plottingSigPmt(grphPmt3, grphCoinciPmt3, "S_{3}", "3");


  // Doing for resolution 
  //


  exit(0);
} 

void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetBottomMargin(0.155);
  canvas.SetLeftMargin(0.09);
  canvas.SetRightMargin(0.02);
}

void fillSigVec(vector<vector<double>> signal,  vector<vector<int>> evtIds,
    vector<vector<vector<double>>> &slctSigPmt, vector<double> &slctSigTot, 
    vector<int> &slctEvtId) { 

  double tmpSum = 0.;
  int crrEvt = evtIds[0][0];
  int binCrrEvt = 0;

  for ( int evt_i=0; evt_i<evtIds[0].size(); evt_i++ ) {
    if ( !(signal[0][evt_i] > 0) || !(signal[1][evt_i] > 0)
        || !(signal[2][evt_i] > 0) )
      continue;
    tmpSum = signal[0][evt_i] + signal[1][evt_i] + signal[2][evt_i];
    tmpSum /= 3.0;
    for ( int i=0; i<3; i++ )
      slctSigPmt[i][binCrrEvt].push_back( 100.*signal[i][evt_i] / tmpSum - 100. );
    slctSigTot.push_back( tmpSum );    
    if ( crrEvt != evtIds[0][evt_i] ) {
      slctEvtId.push_back( evtIds[0][evt_i] );
      crrEvt = evtIds[0][evt_i];
      binCrrEvt++;
    }
    if ( (signal[0][evt_i] / (signal[1][evt_i] + signal[2][evt_i])/2.0) < 0.1 )
      cout << "MSD " << evtIds[0][evt_i] << " " 
        << signal[0][evt_i] << " " << signal[1][evt_i] << " " << signal[2][evt_i] << endl;
  }
}

void plottingSigPmt(TGraphErrors grp1, TGraphErrors grp2, TString pmtName, 
    TString pmtId) {
  TCanvas cvns("cvns", "", 1.6e3, 9e2);
  cvns.cd();
  cvns.SetLogx();

  grp1.SetTitle("");
  grp1.GetXaxis()->SetTitle("#LT S_{(1,2,3)} #GT [VEM]");
  grp1.GetYaxis()->SetTitle(pmtName+" / #LT S_{(1,2,3)} #GT - 100 [%]");
  grp1.SetLineColor(kBlue);
  grp1.SetMarkerColor(kBlue);
  grp1.SetMarkerStyle(54);
  grp1.SetMarkerSize(2);
  grp1.GetXaxis()->SetRangeUser(1, 2e3);
  grp1.GetYaxis()->SetRangeUser(-20, 20);
  grp1.Draw("AP");

  grp2.SetLineColor(kRed);
  grp2.SetMarkerColor(kRed);
  grp2.SetMarkerStyle(53);
  grp2.SetMarkerSize(2);
  grp2.Draw("P same");

  TLegend lgnd(0.7, 0.65, 0.88, 0.88);
  lgnd.AddEntry(&grp1, "Using Q^{pk}", "ep");
  lgnd.AddEntry(&grp2, "Using CQ^{pk}", "ep");

  lgnd.SetBorderSize(0);
  lgnd.SetTextSize(0.055); 
  lgnd.Draw();

  TLine line(1, 0, 2e3, 0);
  line.SetLineStyle(5);
  line.SetLineColor(kGray+2);
  line.Draw();

  cvns.Print("results/signalPMT"+pmtId+"OverTot.pdf");
}
