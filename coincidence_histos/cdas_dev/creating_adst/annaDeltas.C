void setCanvasStyle(TCanvas& canvas);
void setTGraphStyle(TH1D& graph);

void annaDeltas() {
  double tmpStId = 0;
  double tmpPmtId = 0;
  double tmpQpk = 0.;
  double tmpQpkErr = 0.;
  double tmpCQpk = 0.;
  double tmpCQpkErr = 0.;
  double tmpDelta = 0.;

  TH1D delta1("delta1", "", 35, -5, 20);
  TH1D delta2("delta2", "", 35, -5, 20);
  TH1D delta3("delta3", "", 35, -5, 20);
  TH1D delta13("delta13", "", 35, -5, 20);

  TH1D distMsdQpk("distMsdQpk", "", 2e3, 1e3, 3e3);
  TH1D distMsdQpkErr("distMsdQpkErr", "", 100, 0, 100);
  TH1D distMsdCQpk("distMsdCQpk", "", 2e3, 1e3, 3e3);
  TH1D distMsdCQpkErr("distMsdCQpkErr", "", 2e3, 1e3, 3e3);

  vector < vector < double > > deltaPmt;
  vector < vector < double > > countDeltas;
  deltaPmt.resize(3);
  countDeltas.resize(3);
  for ( int i=0; i<3; i++ )
    for ( int j=0; j<2e3; j++ ) {
      deltaPmt[i].push_back( 0. );
      countDeltas[i].push_back( 0. );
  }

  ifstream datafile("qpksForDeltas.dat");

  while( datafile >> tmpStId >> tmpPmtId >> tmpQpk >> tmpQpkErr >> tmpCQpk >> tmpCQpkErr ) {
    tmpDelta = 100.*(tmpCQpk - tmpQpk)/tmpQpk;
    distMsdQpk.Fill( tmpQpk );
    distMsdQpkErr.Fill( tmpQpkErr );
    distMsdCQpk.Fill( tmpCQpk );
    distMsdCQpkErr.Fill( tmpCQpkErr );
    switch ( tmpPmtId ) {
      case 1:
        delta1.Fill( tmpDelta );
        delta13.Fill( tmpDelta );
        deltaPmt[0][int(tmpStId)] += tmpDelta;
        countDeltas[0][int(tmpStId)]++;
        break;
      case 2:
        delta2.Fill( tmpDelta );
        deltaPmt[1][int(tmpStId)] += tmpDelta;
        countDeltas[1][int(tmpStId)]++;
        break;
      case 3:        
        delta3.Fill( tmpDelta );
        delta13.Fill( tmpDelta );
        deltaPmt[2][int(tmpStId)] += tmpDelta;
        countDeltas[2][int(tmpStId)]++;
        break;
    }
  }
  datafile.close();

  double tmpTime = 0.;
  double tmpQpkErr = 0.;
  double tmpCQpkErr = 0.;

  ifstream fileKataFits("cch_fits.csv");
  TH1D distKatQpk("distKatQpk", "", 2e3, 1e3, 3e3);
  TH1D distKatQpkErr("distKatQpkErr", "", 100, 0, 100);
  TH1D distKatCQpk("distKatCQpk", "", 2e3, 1e3, 3e3); 
  TH1D distKatCQpkErr("distKatCQpkErr", "", 2e3, 1e3, 3e3); 

  while( fileKataFits >> tmpStId >> tmpPmtId >> tmpTime >> tmpCQpk >> tmpCQpkErr >> tmpQpk >> tmpQpkErr ) { 
    distKatQpk.Fill( tmpQpk );
    distKatQpkErr.Fill( tmpQpkErr );
    distKatCQpk.Fill( tmpCQpk );
    distKatCQpkErr.Fill( tmpCQpkErr );
  }
  
  for ( int i=0; i<3; i++ )
    for ( int j=0; j<2e3; j++ )
      if ( countDeltas[i][j] > 0 )
        deltaPmt[i][j] /= countDeltas[i][j];

  TCanvas cvnsDeltaDist13("cvnsDeltaDist13","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDeltaDist13);
  cvnsDeltaDist13.cd();
  
  delta13.SetStats(kFALSE); 
  delta13.SetLineColor(kPink+6);
  delta13.GetXaxis()->SetTitle("#Delta [%]");
  delta13.GetYaxis()->SetTitle("Counts [au]");
  delta13.Fit("gaus","Q");
  delta13.GetFunction("gaus")->SetLineColor(kBlack);
  delta13.GetFunction("gaus")->SetLineStyle(9);
  delta13.SetFillColorAlpha(kPink+6, 0.15);
  delta13.Draw();

  delta1.SetLineColor(kRed);
  delta1.Draw("same");

  delta3.SetLineColor(kBlue);
  delta3.Draw("same");


  TLegend *lgnd = new TLegend(0.7, 0.55, 0.9, 0.94);
  lgnd->AddEntry(&delta13, Form("Entries for #Delta_{13}: %.f",
        delta13.GetEntries()), "h");
  lgnd->AddEntry(&delta1, "PMT1", "l"); 
  lgnd->AddEntry(&delta3, "PMT3", "l");
  lgnd->AddEntry(&delta13, "PMT1 and PMT3", "f");
  lgnd->AddEntry(&delta13.GetFunction("gaus"), 
      "Gaussian Fit", "l");
  lgnd->AddEntry(&delta13.GetFunction("gaus"), 
      Form("#mu = %.2f #pm %.2f", 
        delta13.GetFunction("gaus")->GetParameter(1), 
        delta13.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(&delta13.GetFunction("gaus"), 
      Form("#sigma = %.2f #pm %.2f", 
        delta13.GetFunction("gaus")->GetParameter(2), 
        delta13.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDeltaDist13.Print("results/deltaDist13.pdf");


  TCanvas cvnsDeltaDist2("cvnsDeltaDist2","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDeltaDist2);
  cvnsDeltaDist2.cd();

  delta2.SetStats(kFALSE);
  delta2.SetLineColor(kPink+10);
  delta2.SetLineColor(kOrange);
  delta2.GetXaxis()->SetTitle("#Delta [%]");
  delta2.GetYaxis()->SetTitle("Counts [au]");
  delta2.Fit("gaus","Q");
  delta2.GetFunction("gaus")->SetLineColor(kBlack);
  delta2.GetFunction("gaus")->SetLineStyle(9);
  delta2.SetFillColorAlpha(kOrange+10, 0.15);
  delta2.Draw();

  lgnd = new TLegend(0.7, 0.6, 0.9, 0.94);
  lgnd->AddEntry(&delta2, Form("Entries: %.f",
        delta2.GetEntries()), "h");
  lgnd->AddEntry(&delta2, "PMT2", "f");
  lgnd->AddEntry(&delta2.GetFunction("gaus"), 
      "Gaussian Fit", "l");
  lgnd->AddEntry(&delta2.GetFunction("gaus"), 
      Form("#mu = %.2f #pm %.2f", 
        delta2.GetFunction("gaus")->GetParameter(1), 
        delta2.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(&delta13.GetFunction("gaus"), 
      Form("#sigma = %.2f #pm %.2f", 
        delta2.GetFunction("gaus")->GetParameter(2), 
        delta2.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDeltaDist2.Print("results/deltaDist2.pdf");


  TH1D deltaPmt1("deltaPmt1", "", 35, -5, 20);
  TH1D deltaPmt2("deltaPmt2", "", 35, -5, 20);
  TH1D deltaPmt3("deltaPmt3", "", 35, -5, 20);
  TH1D deltaPmt13("deltaPmt13", "", 35, -5, 20);

  for ( int i=0; i<3; i++ )
    for ( int j=0; j<2e3; j++ )
      if ( countDeltas[i][j] > 0 )
        switch ( i ) {
          case 0:
            deltaPmt1.Fill( deltaPmt[i][j] );
            deltaPmt13.Fill( deltaPmt[i][j] );
            break;
          case 1:
            deltaPmt2.Fill( deltaPmt[i][j] );
            break;
          case 2:
            deltaPmt3.Fill( deltaPmt[i][j] );
            deltaPmt13.Fill( deltaPmt[i][j] );
            break;
        }

  TCanvas cvnsDeltaPmt13("cvnsDeltaPmt13","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDeltaPmt13);
  cvnsDeltaPmt13.cd();
  
  deltaPmt13.SetStats(kFALSE); 
  deltaPmt13.SetLineColor(kPink+6);
  deltaPmt13.GetXaxis()->SetTitle("#Delta [%]");
  deltaPmt13.GetYaxis()->SetTitle("Counts [au]");
  deltaPmt13.Fit("gaus","Q","R", 0, 5);
  deltaPmt13.GetFunction("gaus")->SetLineColor(kBlack);
  deltaPmt13.GetFunction("gaus")->SetLineStyle(9);
  deltaPmt13.SetFillColorAlpha(kPink+6, 0.15);
  deltaPmt13.Draw();

  deltaPmt1.SetLineColor(kRed);
  deltaPmt1.Draw("same");

  deltaPmt3.SetLineColor(kBlue);
  deltaPmt3.Draw("same");


  TLegend *lgnd = new TLegend(0.7, 0.55, 0.9, 0.94);
  lgnd->AddEntry(&deltaPmt13, Form("Entries for #Delta_{13}: %.f",
        deltaPmt13.GetEntries()), "h");
  lgnd->AddEntry(&deltaPmt1, "PMT1", "l"); 
  lgnd->AddEntry(&deltaPmt3, "PMT3", "l");
  lgnd->AddEntry(&deltaPmt13, "PMT1 and PMT3", "f");
  lgnd->AddEntry(&deltaPmt13.GetFunction("gaus"), 
      "Gaussian Fit", "l");
  lgnd->AddEntry(&deltaPmt13.GetFunction("gaus"), 
      Form("#mu = %.2f #pm %.2f", 
        deltaPmt13.GetFunction("gaus")->GetParameter(1), 
        deltaPmt13.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(&deltaPmt13.GetFunction("gaus"), 
      Form("#sigma = %.2f #pm %.2f", 
        deltaPmt13.GetFunction("gaus")->GetParameter(2), 
        deltaPmt13.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDeltaPmt13.Print("results/deltaPmt13.pdf");


  TCanvas cvnsDeltaPmt2("cvnsDeltaPmt2","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDeltaPmt2);
  cvnsDeltaPmt2.cd();

  deltaPmt2.SetStats(kFALSE);
  deltaPmt2.SetLineColor(kPink+10);
  deltaPmt2.SetLineColor(kOrange);
  deltaPmt2.GetXaxis()->SetTitle("#Delta [%]");
  deltaPmt2.GetYaxis()->SetTitle("Counts [au]");
  deltaPmt2.Fit("gaus","Q", "R", 3., 10);
  deltaPmt2.GetFunction("gaus")->SetLineColor(kBlack);
  deltaPmt2.GetFunction("gaus")->SetLineStyle(9);
  deltaPmt2.SetFillColorAlpha(kOrange+10, 0.15);
  deltaPmt2.Draw();

  lgnd = new TLegend(0.7, 0.6, 0.9, 0.94);
  lgnd->AddEntry(&deltaPmt2, Form("Entries: %.f",
        deltaPmt2.GetEntries()), "h");
  lgnd->AddEntry(&deltaPmt2, "PMT2", "f");
  lgnd->AddEntry(&deltaPmt2.GetFunction("gaus"), 
      "Gaussian Fit", "l");
  lgnd->AddEntry(&deltaPmt2.GetFunction("gaus"), 
      Form("#mu = %.2f #pm %.2f", 
        deltaPmt2.GetFunction("gaus")->GetParameter(1), 
        deltaPmt2.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(&delta13.GetFunction("gaus"), 
      Form("#sigma = %.2f #pm %.2f", 
        deltaPmt2.GetFunction("gaus")->GetParameter(2), 
        deltaPmt2.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDeltaPmt2.Print("results/deltaPmt2.pdf");


  TCanvas cvnsKatQpkVsMsdQpk("cvnsKatQpkVsMsdQpk","", 1.6e3, 9e2);
  setCanvasStyle(cvnsKatQpkVsMsdQpk);
  cvnsKatQpkVsMsdQpk.cd();

  distKatQpk.SetStats(kFALSE);
  distKatQpk.SetLineColor(kRed);
  distKatQpk.GetXaxis()->SetTitle("Qpk [FADC]");
  distKatQpk.GetYaxis()->SetTitle("Counts [au]");
  distKatQpk.Draw();

  distMsdQpk.SetLineColor(kBlue);
  distMsdQpk.Draw("same");

  lgnd = new TLegend(0.7, 0.4, 0.9, 0.94);
  lgnd->AddEntry(&distKatQpk, "Independent:","l");
  lgnd->AddEntry(&distKatQpk, Form("Entries: %.f",
        distKatQpk.GetEntries()), "h");
  lgnd->AddEntry(&distKatQpk, Form("#mu: %.f #pm %.f",
        distKatQpk.GetMean(), distKatQpk.GetMeanError()), "h");
  lgnd->AddEntry(&distKatQpk, Form("RMS: %.f #pm %.f",
        distKatQpk.GetRMS(), distKatQpk.GetRMSError()), "h");
  lgnd->AddEntry(&distMsdQpk, "OffLine:","l");
  lgnd->AddEntry(&distMsdQpk, Form("Entries: %.f",
        distMsdQpk.GetEntries()), "h");
  lgnd->AddEntry(&distMsdQpk, Form("#mu: %.f #pm %.f",
        distMsdQpk.GetMean(), distMsdQpk.GetMeanError()), "h");
  lgnd->AddEntry(&distMsdQpk, Form("RMS: %.f #pm %.f",
        distMsdQpk.GetRMS(), distMsdQpk.GetRMSError()), "h");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsKatQpkVsMsdQpk.Print("results/KatQpkVsMsdQpk.pdf");


  TCanvas cvnsKatQpkErrVsMsdQpkErr("cvnsKatQpkErrVsMsdQpkErr","", 1.6e3, 9e2);
  setCanvasStyle(cvnsKatQpkErrVsMsdQpkErr);
  cvnsKatQpkErrVsMsdQpkErr.cd();

  cvnsKatQpkErrVsMsdQpkErr.SetLogy();

  distKatQpkErr.SetStats(kFALSE);
  distKatQpkErr.SetLineColor(kRed);
  distKatQpkErr.GetXaxis()->SetTitle("Qpk [FADC]");
  distKatQpkErr.GetYaxis()->SetTitle("Counts [au]");
  distKatQpkErr.GetXaxis()->SetRangeUser(0, 40);
  distKatQpkErr.Draw();

  distMsdQpkErr.SetLineColor(kBlue);
  distMsdQpkErr.Draw("same");

  lgnd = new TLegend(0.7, 0.4, 0.9, 0.94);
  lgnd->AddEntry(&distKatQpkErr, "Independent:","l");
  lgnd->AddEntry(&distKatQpkErr, Form("Entries: %.f",
        distKatQpkErr.GetEntries()), "h");
  lgnd->AddEntry(&distKatQpkErr, Form("#mu: %.f #pm %.1f",
        distKatQpkErr.GetMean(), distKatQpkErr.GetMeanError()), "h");
  lgnd->AddEntry(&distKatQpkErr, Form("RMS: %.f #pm %.1f",
        distKatQpkErr.GetRMS(), distKatQpkErr.GetRMSError()), "h");
  lgnd->AddEntry(&distMsdQpkErr, "OffLine:","l");
  lgnd->AddEntry(&distMsdQpkErr, Form("Entries: %.f",
        distMsdQpkErr.GetEntries()), "h");
  lgnd->AddEntry(&distMsdQpkErr, Form("#mu: %.f #pm %.1f",
        distMsdQpkErr.GetMean(), distMsdQpkErr.GetMeanError()), "h");
  lgnd->AddEntry(&distMsdQpkErr, Form("RMS: %.f #pm %.1f",
        distMsdQpkErr.GetRMS(), distMsdQpkErr.GetRMSError()), "h");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsKatQpkErrVsMsdQpkErr.Print("results/KatQpkErrVsMsdQpkErr.pdf");



  exit(0);
}

void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetLeftMargin(0.08);
  canvas.SetRightMargin(0.02);
}

void setTGraphStyle(TH1D& graph) {
  graph.GetYaxis()->SetLabelSize(0.05);
  graph.GetYaxis()->SetTitleSize(0.05);
  graph.GetXaxis()->SetLabelSize(0.05);
  graph.GetXaxis()->SetTitleSize(0.05);
}
