void setCanvasStyle(TCanvas& canvas);
void plotCoinciFactor(TH1D &dist, TString pmt);

void fromQpk2CQpk() {
  ifstream data("qpksForCoinciFactor.dat");

  // Defining vectors to store Qpk, CQpk and their errors
  vector < vector < vector < double > > > qpkSt(1900);
  vector < vector < vector < double > > > qpkErrSt(1900);
  vector < vector < vector < double > > > cqpkSt(1900);
  vector < vector < vector < double > > > cqpkErrSt(1900);
  vector < vector < vector < int > > > timeSt(1900);

  // Re-sizing per PMT
  for ( int i=0; i<1900; i++ ) {
    qpkSt[i].resize(3);
    qpkErrSt[i].resize(3);
    cqpkSt[i].resize(3);
    cqpkErrSt[i].resize(3);
    timeSt[i].resize(3);
  }

  // Temporal variables to fetch data
  double tmpStId = 0.;
  double tmpPmtId = 0.;
  double tmpQpk = 0.;
  double tmpQpkErr = 0.;
  double tmpCQpk = 0.;
  double tmpCQpkErr = 0.;
  double tmpTime = 0.;

  // Reading and storing data
  while ( data >> tmpStId >> tmpPmtId >> tmpTime
      >> tmpQpk >> tmpQpkErr >> tmpCQpk >> tmpCQpkErr ) {
    qpkSt[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpQpk );
    qpkErrSt[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpQpkErr );
    cqpkSt[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpCQpk );
    cqpkErrSt[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpCQpkErr );
    timeSt[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpTime );
  }
  data.close();

  int nBins = 6e2;
  TH1D coinciFactorPmt1("coinciFactorPmt1", "", nBins, -20, 40);
  TH1D coinciFactorPmt2("coinciFactorPmt2", "", nBins, -20, 40);
  TH1D coinciFactorPmt3("coinciFactorPmt3", "", nBins, -20, 40);

  // Building RelErr vs Station
  double tmpDelta = 0.;
  for ( int st_i=0; st_i<1900; st_i++ ) {
    for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
      if ( timeSt[st_i][pmt_i].size() < 1 )
        continue;
      for ( int utc = 0; utc < timeSt[st_i][pmt_i].size(); utc++ ) {
        tmpDelta = 100.*(cqpkSt[st_i][pmt_i][utc] - qpkSt[st_i][pmt_i][utc]);
        tmpDelta /= qpkSt[st_i][pmt_i][utc];
        switch( pmt_i ) {
          case 0:
            coinciFactorPmt1.Fill( tmpDelta );
            break;
          case 1:
            coinciFactorPmt2.Fill( tmpDelta );
            break;
          case 2:
            coinciFactorPmt3.Fill( tmpDelta );
            break;
        }
      }
    }
  }

  plotCoinciFactor(coinciFactorPmt1, "1");
  plotCoinciFactor(coinciFactorPmt2, "2");
  plotCoinciFactor(coinciFactorPmt3, "3");

  exit(0);
}

void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetLeftMargin(0.08);
  canvas.SetRightMargin(0.02);
}

void plotCoinciFactor(TH1D &dist, TString pmt) {
  TCanvas cvnsFactorPmt("cvnsFactorPmt"+pmt,"", 1.6e3, 9e2);
  setCanvasStyle(cvnsFactorPmt);
  cvnsFactorPmt.SetLogy();

  dist.Fit("gaus", "Q");
  dist.GetFunction("gaus")->SetLineStyle(9);
  dist.GetFunction("gaus")->SetLineColor(kBlack);

  dist.SetStats(kFALSE);
  dist.SetLineColor(kBlue);
  dist.SetLineWidth(2);
  dist.GetXaxis()->SetTitle("#Delta [%]");
  dist.GetYaxis()->SetTitle("Counts [au]");
  dist.SetFillColorAlpha(kBlue, 0.15);
  dist.Draw();

  TLegend *lgnd = new TLegend(0.65, 0.5, 0.9, 0.96);
  lgnd->AddEntry(&dist, "PMT "+pmt, "h");
  lgnd->AddEntry(&dist, Form( "Mean: %.2f [%] #pm %.2f [%]",
        dist.GetMean(), dist.GetMeanError()), "l");
  lgnd->AddEntry(&dist, Form( "RMS: %.2f [%] #pm %.2f [%]",
        dist.GetRMS(), dist.GetRMSError()), "l");
  lgnd->AddEntry(dist.GetFunction("gaus"), "Gaussian fit", "l");
  lgnd->AddEntry(dist.GetFunction("gaus"), Form(
        "#mu: %.2f [%] #pm %.2f", 
        dist.GetFunction("gaus")->GetParameter(1),
        dist.GetFunction("gaus")->GetParError(1)),"");
  lgnd->AddEntry(dist.GetFunction("gaus"), Form( 
        "#sigma: %.2f [%] #pm %.2f", 
        dist.GetFunction("gaus")->GetParameter(2),
        dist.GetFunction("gaus")->GetParError(2)),"");
  lgnd->AddEntry(&dist, "", "h");

  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.04);
  lgnd->Draw();


  cvnsFactorPmt.Print("results/CQpkOverQpkPmt"+pmt+".pdf");
}
