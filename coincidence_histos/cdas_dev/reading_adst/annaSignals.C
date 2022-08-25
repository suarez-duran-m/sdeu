void setCanvasStyle(TCanvas& canvas);
TGraph fillingSiVsStot(vector<vector<double>> signals, int pmt);
void plottingSiVsStot(TString cnvsName, TGraph &graph, TString outName);
TGraph fillingAsymmetry(vector<vector<double>> signals, int pmt);
void plottingAsymmetry(TString cnvsName, TGraph &graph, TString outName);
void fillingAsymmetryDist(TH1D &dist, vector<vector<double>> signals, int pmt);
void plottingAsymmetryDist(TString cnvsName, TH1D &dist, TString outName);
TGraph fillingDiffAsymmetry(vector<vector<double>> signals, 
    vector<vector<double>> signalsCoinci, int pmt);
void plottingDiffAsymmetry(TString cnvsName, TGraph &graph, TString outName);
void fillingDiff(TH1D &hist, vector<double> pmtA, vector<double> pmtB);
void plottingDiff(TString cnvsName, TH1D &hist, TString outName);

void annaSignals() {
  ifstream dataVEM("signalsVEM.dat");
  ifstream dataVEMCoinci("signalsVEMCoinci.dat");

  double tmpEvtId = 0.;
  double tmpSt = 0.;
  double tmpPmt = 0.;
  double tmpSignal = 0.;

  // To storage signals 
  // signals[pmt_i][signal_i]
  vector < vector < double > > signalsVEM(3);
  vector < vector < double > > signalsVEMCoinci(3);
  vector < int > stIdsVEM;
  vector < int > stIdsVEMCoinci;
  vector < int > evtIdsVEM;
  vector < int > evtIdsVEMCoinci;

  while ( dataVEM >> tmpEvtId >> tmpSt >> tmpPmt >> tmpSignal ) {
    evtIdsVEM.push_back( tmpEvtId );
    stIdsVEM.push_back( tmpSt );
    signalsVEM[int(tmpPmt-1)].push_back( tmpSignal );
  }
  dataVEM.close();

  while ( dataVEMCoinci >> tmpEvtId >> tmpSt >> tmpPmt >> tmpSignal ) {
    evtIdsVEMCoinci.push_back( tmpEvtId );
    stIdsVEMCoinci.push_back( tmpSt );
    signalsVEMCoinci[int(tmpPmt-1)].push_back( tmpSignal );
  }
  dataVEMCoinci.close();

  // Doing for relative signal magnitude
  //
  TGraph relative_S1 = fillingSiVsStot(signalsVEM, 1);
  relative_S1.SetTitle("; log #LT S_{1,2,3} #GT; S_{1} / #LT S_{2,3} #GT");
  plottingSiVsStot("cvnsRelative_S1", relative_S1, "relative_S1");

  TGraph relative_S2 = fillingSiVsStot(signalsVEM, 2);
  relative_S2.SetTitle("; log #LT S_{1,2,3} #GT ; S_{2} / #LT S_{1,3} #GT ");
  plottingSiVsStot("cvnsRelative_S2", relative_S2, "relative_S2");

  TGraph relative_S3 = fillingSiVsStot(signalsVEM, 3);
  relative_S3.SetTitle("; log #LT S_{1,2,3} #GT ; S_{3} / #LT S_{1,2} #GT ");
  plottingSiVsStot("cvnsRelative_S3", relative_S3, "relative_S3");

  TGraph relativeCoinci_S1 = fillingSiVsStot(signalsVEMCoinci, 1);
  relativeCoinci_S1.SetTitle("; log #LT S_{1,2,3} #GT ; S_{1} / #LT S_{2,3} #GT ");
  plottingSiVsStot("cvnsRelativeCoinci_S1", relativeCoinci_S1, "relativeCoinci_S1");

  TGraph relativeCoinci_S2 = fillingSiVsStot(signalsVEMCoinci, 2);
  relativeCoinci_S2.SetTitle("; log #LT S_{1,2,3} #GT ; S_{2} / #LT S_{1,3} #GT ");
  plottingSiVsStot("cvnsRelativeCoinci_S2", relativeCoinci_S2, "relativeCoinci_S2");

  TGraph relativeCoinci_S3 = fillingSiVsStot(signalsVEMCoinci, 3);
  relativeCoinci_S3.SetTitle("; log #LT S_{1,2,3} #GT ; S_{3} / #LT S_{1,2} #GT ");
  plottingSiVsStot("cvnsRelativeCoinci_S3", relativeCoinci_S3, "relativeCoinci_S3");

  // Doing for asymmetry of signal
  //  
  TGraph asymmetry_S1 = fillingAsymmetry(signalsVEM, 1);
  asymmetry_S1.SetTitle("; log #LT S_{1,2,3} #GT; #left( S_{1} - S_{2} #right) / #LT S_{1,2} #GT");
  plottingAsymmetry("cvnsRelative_S1", asymmetry_S1, "asymmetry_S1");

  TGraph asymmetry_S2 = fillingAsymmetry(signalsVEM, 2);
  asymmetry_S2.SetTitle("; log #LT S_{1,2,3} #GT ; #left( S_{2} - S_{3} #right) / #LT S_{2,3} #GT ");
  plottingAsymmetry("cvnsRelative_S2", asymmetry_S2, "asymmetry_S2");

  TGraph asymmetry_S3 = fillingAsymmetry(signalsVEM, 3);
  asymmetry_S3.SetTitle("; log #LT S_{1,2,3} #GT ; #left( S_{3} - S_{1} #right) / #LT S_{1,3} #GT ");
  plottingAsymmetry("cvnsRelative_S3", asymmetry_S3, "asymmetry_S3");

  TGraph asymmetryCoinci_S1 = fillingAsymmetry(signalsVEMCoinci, 1);
  asymmetryCoinci_S1.SetTitle("; log #LT S_{1,2,3} #GT ; #left( S_{1} - S_{2} #right) / #LT S_{1,2} #GT ");
  plottingAsymmetry("cvnsRelativeCoinci_S1", asymmetryCoinci_S1, "asymmetryCoinci_S1");

  TGraph asymmetryCoinci_S2 = fillingAsymmetry(signalsVEMCoinci, 2);
  asymmetryCoinci_S2.SetTitle("; log #LT S_{1,2,3} #GT ; #left( S_{2} - S_{1} #right) / #LT S_{2,1} #GT ");
  plottingAsymmetry("cvnsRelativeCoinci_S2", asymmetryCoinci_S2, "asymmetryCoinci_S2");

  TGraph asymmetryCoinci_S3 = fillingAsymmetry(signalsVEMCoinci, 3);
  asymmetryCoinci_S3.SetTitle("; log #LT S_{1,2,3} #GT ; #left( S_{3} - S_{1} #right) / #LT S_{3,1} #GT ");
  plottingAsymmetry("cvnsRelativeCoinci_S3", asymmetryCoinci_S3, "asymmetryCoinci_S3");


  // Doing distribution for asymmetry 
  //
  int nBins = 100;

  TH1D distri_S1("distri_S1", "", nBins, -2, 2);
  distri_S1.SetTitle(";  #left( S_{1} - S_{2} #right) /#LT S_{1,2} #GT; Counts [au]");
  fillingAsymmetryDist(distri_S1, signalsVEM, 1);
  plottingAsymmetryDist("cvnsAsyDist_S1", distri_S1, "asymmetryDist_S1");

  TH1D distri_S2("distri_S2", "", nBins, -2, 2);
  distri_S2.SetTitle(";  #left( S_{2} - S_{3} #right) /#LT S_{2,3} #GT; Counts [au]");
  fillingAsymmetryDist(distri_S2, signalsVEM, 2);
  plottingAsymmetryDist("cvnsAsyDist_S2", distri_S2, "asymmetryDist_S2");

  TH1D distri_S3("distri_S3", "", nBins, -2, 2);
  distri_S3.SetTitle(";  #left( S_{3} - S_{1} #right) /#LT S_{3,1} #GT; Counts [au]");
  fillingAsymmetryDist(distri_S3, signalsVEM, 3);
  plottingAsymmetryDist("cvnsAsyDist_S3", distri_S3, "asymmetryDist_S3");

  TH1D distriCoinci_S1("distriCoinci_S1", "", nBins, -2, 2);
  distriCoinci_S1.SetTitle(";  #left( S_{1} - S_{2} #right) /#LT S_{1,2} #GT; Counts [au]");
  fillingAsymmetryDist(distriCoinci_S1, signalsVEMCoinci, 1);
  plottingAsymmetryDist("cvnsAsyDistCoinci_S1", distriCoinci_S1, "asymmetryDistCoinci_S1");

  TH1D distriCoinci_S2("distriCoinci_S2", "", nBins, -2, 2);
  distriCoinci_S2.SetTitle(";  #left( S_{2} - S_{3} #right) /#LT S_{2,3} #GT; Counts [au]");
  fillingAsymmetryDist(distriCoinci_S2, signalsVEMCoinci, 2);
  plottingAsymmetryDist("cvnsAsyDistCoinci_S2", distriCoinci_S2, "asymmetryDistCoinci_S2");

  TH1D distriCoinci_S3("distriCoinci_S3", "", nBins, -2, 2);
  distriCoinci_S3.SetTitle(";  #left( S_{3} - S_{1} #right) /#LT S_{3,1} #GT; Counts [au]");
  fillingAsymmetryDist(distriCoinci_S3, signalsVEMCoinci, 3);
  plottingAsymmetryDist("cvnsAsyDistCoinci_S3", distriCoinci_S3, "asymmetryDistCoinci_S3");


  // Doing for symmetry comparison S_i - S_j 
  //
  TGraph diffAsymmetry_S1 = fillingDiffAsymmetry(signalsVEM, signalsVEMCoinci, 1);
  diffAsymmetry_S1.SetTitle("; log #LT S_{1,2,3} #GT; RelS_{1,2} [%]");
  plottingDiffAsymmetry("cvnsDiffAsymmetry_S1", diffAsymmetry_S1, "diffAsymmetry_S1");

  TGraph diffAsymmetry_S2 = fillingDiffAsymmetry(signalsVEM, signalsVEMCoinci, 2);
  diffAsymmetry_S2.SetTitle("; log #LT S_{1,2,3} #GT ; RelS_{2,3} [%]");
  plottingDiffAsymmetry("cvnsDiffAsymmetry_S2", diffAsymmetry_S2, "diffAsymmetry_S2");

  TGraph diffAsymmetry_S3 = fillingDiffAsymmetry(signalsVEM, signalsVEMCoinci, 3);
  diffAsymmetry_S3.SetTitle("; log #LT S_{1,2,3} #GT ; RelS_{3,1} [%]");
  plottingDiffAsymmetry("cvnsDiffAsymmetry_S3", diffAsymmetry_S3, "diffAsymmetry_S3");


  nBins = 100;
  TH1D diffS_12("diffS_12","", nBins, -1, 1); 
  diffS_12.SetTitle("; (S_{1} - S_{2}) / #LT S_{1,2} #GT; Counts [au]");
  TH1D diffS_23("diffS_23","", nBins, -1, 1);
  diffS_23.SetTitle("; (S_{2} - S_{3}) / #LT S_{2,3} #GT; Counts [au]");
  TH1D diffS_31("diffS_31","", nBins, -1, 1);
  diffS_31.SetTitle("; (S_{3} - S_{1}) / #LT S_{3,1} #GT; Counts [au]");

  fillingDiff( diffS_12, signalsVEM[0], signalsVEM[1] );
  fillingDiff( diffS_23, signalsVEM[1], signalsVEM[2] );
  fillingDiff( diffS_31, signalsVEM[2], signalsVEM[0] );

  plottingDiff("cvnsDiffS_12", diffS_12, "diffS_12");
  plottingDiff("cvnsDiffS_23", diffS_23, "diffS_23");
  plottingDiff("cvnsDiffS_31", diffS_31, "diffS_31");
 

  TH1D diffSCoinci_12("diffSCoinci_12","", nBins, -1, 1);
  diffSCoinci_12.SetTitle("; (S_{1} - S_{2}) / #LT S_{1,2} #GT; Counts[au]");
  TH1D diffSCoinci_23("diffSCoinci_23","", nBins, -1, 1);
  diffSCoinci_23.SetTitle("; (S_{2} - S_{3}) / #LT S_{2,3} #GT; Counts[au]");
  TH1D diffSCoinci_31("diffSCoinci_31","", nBins, -1, 1);
  diffSCoinci_31.SetTitle("; (S_{3} - S_{1}) / #LT S_{3,1} #GT; Counts[au]");

  fillingDiff( diffSCoinci_12, signalsVEM[0], signalsVEM[1] );
  fillingDiff( diffSCoinci_23, signalsVEM[1], signalsVEM[2] );
  fillingDiff( diffSCoinci_31, signalsVEM[2], signalsVEM[0] );

  plottingDiff("cvnsDiffSCoinci_12", diffSCoinci_12, "diffSCoinci_12");
  plottingDiff("cvnsDiffSCoinci_23", diffSCoinci_23, "diffSCoinci_23");
  plottingDiff("cvnsDiffSCoinci_31", diffSCoinci_31, "diffSCoinci_31");

  exit(0);
}

void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetBottomMargin(0.155);
  canvas.SetLeftMargin(0.09);
  canvas.SetRightMargin(0.02);
}

TGraph fillingSiVsStot(vector<vector<double>> signals, int pmt) {
  vector < double > Sa;
  vector < double > Sabc;
  int a, b, c = 0;
  switch( pmt ) {
    case 1:
    a = 0;
    b = 1;
    c = 2;
    break;
    case 2:
    a = 1;
    b = 2;
    c = 0;
    break;
    case 3:
    a = 2;
    b = 0;
    c = 1;
    break;
  }
  double tmpVal = 0.;
  for ( int i=0; i<signals[a].size(); i++ )
    if (signals[a][i] > 0 && signals[b][i] > 0 && signals[c][i] > 0 ) {
      tmpVal = signals[b][i] + signals[c][i];
      Sa.push_back( signals[a][i]/(tmpVal/2.0) );
      Sabc.push_back( (tmpVal+signals[a][i])/3.0 );
    }
  TGraph retGraph( Sa.size(), &Sabc.front(), &Sa.front() );
  Sa.clear();
  Sabc.clear();
  return retGraph;
}

void plottingSiVsStot(TString cnvsName, TGraph &graph, TString outName) {
  TCanvas cvnsRelaS(cnvsName, "", 1.6e3, 9e2);
  setCanvasStyle(cvnsRelaS);
  cvnsRelaS.SetLogy();
  cvnsRelaS.SetLogx();

  graph.SetMarkerStyle(3);
  graph.SetMarkerSize(2);
  graph.SetMarkerColor(kBlue);
  
  graph.GetXaxis()->SetLabelSize(0.05);
  graph.GetXaxis()->SetTitleSize(0.06);
  graph.GetXaxis()->SetTitleOffset(1.2);
  graph.GetXaxis()->SetRangeUser(1e-1,1e4);

  graph.GetYaxis()->SetLabelSize(0.05);
  graph.GetYaxis()->SetTitleSize(0.06);
  graph.GetYaxis()->SetTitleOffset(0.7);
  graph.GetYaxis()->SetRangeUser(1e-2,1e2);

  graph.Draw("AP");
  cvnsRelaS.Print("results/"+outName+".pdf");
}

TGraph fillingAsymmetry(vector<vector<double>> signals, int pmt) {
  vector < double > Sa;
  vector < double > Sabc;
  int a, b, c = 0;
  switch( pmt ) {
    case 1:
    a = 0;
    b = 1;
    c = 2;
    break;
    case 2:
    a = 1;
    b = 2;
    c = 0;
    break;
    case 3:
    a = 2;
    b = 0;
    c = 1;
    break;
  }
  double tmpVal = 0.;
  double tmpDiff = 0.;
  for ( int i=0; i<signals[a].size(); i++ )
    if (signals[a][i] > 0 && signals[b][i] > 0 && signals[c][i] > 0 ) {
      tmpVal = signals[a][i] + signals[b][i];
      tmpDiff = signals[a][i]-signals[b][i];
      Sa.push_back( ( tmpDiff / (tmpVal/2.0) ) );
      Sabc.push_back( (tmpVal+signals[c][i])/3.0 );
    }
  TGraph retGraph( Sa.size(), &Sabc.front(), &Sa.front() );
  Sa.clear();
  Sabc.clear();
  return retGraph;
}

void plottingAsymmetry(TString cnvsName, TGraph &graph, TString outName) {
  TCanvas cvnsRelaS(cnvsName, "", 1.6e3, 9e2);
  setCanvasStyle(cvnsRelaS);
  cvnsRelaS.SetLogx();

  graph.SetMarkerStyle(3);
  graph.SetMarkerSize(2);
  graph.SetMarkerColor(kBlue);
  
  graph.GetXaxis()->SetLabelSize(0.05);
  graph.GetXaxis()->SetTitleSize(0.06);
  graph.GetXaxis()->SetTitleOffset(1.2);
  graph.GetXaxis()->SetRangeUser(1e-1,1e4);

  graph.GetYaxis()->SetLabelSize(0.05);
  graph.GetYaxis()->SetTitleSize(0.06);
  graph.GetYaxis()->SetTitleOffset(0.7);
  graph.GetYaxis()->SetRangeUser(-2,2);

  graph.Draw("AP");
  cvnsRelaS.Print("results/"+outName+".pdf");
}

void fillingAsymmetryDist(TH1D &dist, vector<vector<double>> signals, int pmt) {
  int a, b, c = 0;
  switch( pmt ) {
    case 1:
    a = 0;
    b = 1;
    c = 2;
    break;
    case 2:
    a = 1;
    b = 2;
    c = 0;
    break;
    case 3:
    a = 2;
    b = 0;
    c = 1;
    break;
  }
  double tmpVal = 0.;
  double tmpDiff = 0.;  
  for ( int i=0; i<signals[a].size(); i++ )
    if (signals[a][i] > 0 && signals[b][i] > 0 && signals[c][i] > 0 ) {
      tmpVal = signals[a][i] + signals[b][i];
      tmpDiff = signals[a][i]-signals[b][i];
      dist.Fill( tmpDiff / (tmpVal/2.0) ) ;      
    }
}

void plottingAsymmetryDist(TString cnvsName, TH1D &dist, TString outName) {
  TCanvas cvnsAsyDist(cnvsName, "", 1.6e3, 9e2);
  setCanvasStyle(cvnsAsyDist);
  cvnsAsyDist.SetLogy();

  dist.SetStats(kFALSE);
  dist.SetLineColor(kBlue);
  dist.SetLineWidth(2);
  
  dist.GetXaxis()->SetLabelSize(0.05);
  dist.GetXaxis()->SetTitleSize(0.06);
  dist.GetXaxis()->SetTitleOffset(1.2);

  dist.GetYaxis()->SetLabelSize(0.05);
  dist.GetYaxis()->SetTitleSize(0.06);
  dist.GetYaxis()->SetTitleOffset(0.7);

  dist.Fit("gaus", "Q");
  dist.Draw();

  TLegend lgnd(0.62, 0.8, 0.9, 0.96);
  lgnd.AddEntry(dist.GetFunction("gaus"),
      Form("#mu = %.1e #pm %.1e",
        dist.GetFunction("gaus")->GetParameter(1), 
        dist.GetFunction("gaus")->GetParError(1)), "l");
  lgnd.AddEntry(dist.GetFunction("gaus"),
      Form("#sigma = %.1e #pm %.1e",
        dist.GetFunction("gaus")->GetParameter(2), 
        dist.GetFunction("gaus")->GetParError(2)), "l");
  lgnd.SetBorderSize(0);
  lgnd.SetTextSize(0.055);
  lgnd.Draw();
  cvnsAsyDist.Print("results/"+outName+".pdf");
}


TGraph fillingDiffAsymmetry(vector<vector<double>> signals, 
  vector<vector<double>> signalsCoinci, int pmt) {
  vector < double > Sa;
  vector < double > Sabc;
  int a, b, c = 0;
  switch( pmt ) {
    case 1:
    a = 0;
    b = 1;
    c = 2;
    break;
    case 2:
    a = 1;
    b = 2;
    c = 0;
    break;
    case 3:
    a = 2;
    b = 0;
    c = 1;
    break;
  }
  double tmpVal = 0.;
  double tmpValCoinci = 0.;
  double tmpDiff = 0.;
  double tmpDiffCoinci = 0.;
  for ( int i=0; i<signals[a].size(); i++ ) {
    if ( !(signals[a][i] > 0) || !(signals[b][i] > 0) || !(signals[c][i] > 0) )
      continue;
    if ( !(signalsCoinci[a][i] > 0) || !(signalsCoinci[b][i] > 0) || !(signalsCoinci[c][i] > 0) )
      continue;
    tmpVal = signals[a][i] + signals[b][i];
    tmpDiff = signals[a][i] - signals[b][i];
    tmpDiff /= tmpVal / 2.0;
    tmpValCoinci = signalsCoinci[a][i] + signalsCoinci[b][i];
    tmpDiffCoinci = signalsCoinci[a][i] - signalsCoinci[b][i];
    tmpDiffCoinci /= tmpValCoinci / 2.0;

    Sa.push_back( 100.*(tmpDiffCoinci/tmpDiff - 1) );
    Sabc.push_back( (tmpVal+signals[c][i])/3.0 );
  }
  TGraph retGraph( Sa.size(), &Sabc.front(), &Sa.front() );
  Sa.clear();
  Sabc.clear();
  return retGraph;
}

void plottingDiffAsymmetry(TString cnvsName, TGraph &graph, TString outName) {
  TCanvas cvnsRelaS(cnvsName, "", 1.6e3, 9e2);
  setCanvasStyle(cvnsRelaS);
  cvnsRelaS.SetLogx();

  graph.SetMarkerStyle(3);
  graph.SetMarkerSize(2);
  graph.SetMarkerColor(kBlue);
  
  graph.GetXaxis()->SetLabelSize(0.05);
  graph.GetXaxis()->SetTitleSize(0.06);
  graph.GetXaxis()->SetTitleOffset(1.2);
  graph.GetXaxis()->SetRangeUser(1e1,1e4);

  graph.GetYaxis()->SetLabelSize(0.05);
  graph.GetYaxis()->SetTitleSize(0.06);
  graph.GetYaxis()->SetTitleOffset(0.7);
  graph.GetYaxis()->SetRangeUser(-1e2,1e2);

  graph.Draw("AP");

  TLine ceroLine(1e1, 0, 1e4, 0);
  ceroLine.Draw();
  cvnsRelaS.Print("results/"+outName+".pdf");
}


void fillingDiff(TH1D &hist, vector<double> pmtA, vector<double> pmtB) {
  double average = 0.;
  int nEvts = 0;
  for ( int i=0; i<pmtA.size(); i++ )
    if ( pmtA[i] > 0 && pmtB[i] > 0 ) {
      average += pmtA[i] + pmtB[i];
      nEvts++;
    }
  average /= nEvts;

  for ( int i=0; i<pmtA.size(); i++ ) 
    if ( pmtA[i] > 0 && pmtB[i] > 0 )
      hist.Fill( (pmtA[i] - pmtB[i]) / average );
}

void plottingDiff(TString cnvsName, TH1D &hist, TString outName) {
  TCanvas cvnsDiffS(cnvsName, "", 1.6e3, 9e2);
  setCanvasStyle(cvnsDiffS);

  cvnsDiffS.cd();
  cvnsDiffS.SetLogy();

  hist.Draw();
  cvnsDiffS.Print("results/"+outName+".pdf");  
}
