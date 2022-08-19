void fillingDiff(TH1D &hist, vector<double> pmtA, vector<double> pmtB);
void setCanvasStyle(TCanvas& canvas);
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

  int nBins = 100;
  TH1D diffS_12("diffS_12","", nBins, -1, 1); 
  diffS_12.SetTitle("; (S_1 - S_2) / <S_{1,2}>; Counts [au]");
  TH1D diffS_23("diffS_23","", nBins, -1, 1);
  diffS_23.SetTitle("; (S_2 - S_3) / <S_{2,3}>; Counts [au]");
  TH1D diffS_31("diffS_31","", nBins, -1, 1);
  diffS_31.SetTitle("; (S_3 - S_1) / <S_{3,1}>; Counts [au]");

  fillingDiff( diffS_12, signalsVEM[0], signalsVEM[1] );
  fillingDiff( diffS_23, signalsVEM[1], signalsVEM[2] );
  fillingDiff( diffS_31, signalsVEM[2], signalsVEM[0] );

  plottingDiff("cvnsDiffS_12", diffS_12, "diffS_12");
  plottingDiff("cvnsDiffS_23", diffS_23, "diffS_23");
  plottingDiff("cvnsDiffS_31", diffS_31, "diffS_31");
 

  TH1D diffSCoinci_12("diffSCoinci_12","", nBins, -1, 1);
  diffSCoinci_12.SetTitle("; (S_1 - S_2) / <S_{1,2}>; Counts[au]");
  TH1D diffSCoinci_23("diffSCoinci_23","", nBins, -1, 1);
  diffSCoinci_23.SetTitle("; (S_2 - S_3) / <S_{2,3}>; Counts[au]");
  TH1D diffSCoinci_31("diffSCoinci_31","", nBins, -1, 1);
  diffSCoinci_31.SetTitle("; (S_3 - S_1) / <S_{3,1}>; Counts[au]");

  fillingDiff( diffSCoinci_12, signalsVEM[0], signalsVEM[1] );
  fillingDiff( diffSCoinci_23, signalsVEM[1], signalsVEM[2] );
  fillingDiff( diffSCoinci_31, signalsVEM[2], signalsVEM[0] );

  plottingDiff("cvnsDiffSCoinci_12", diffSCoinci_12, "diffSCoinci_12");
  plottingDiff("cvnsDiffSCoinci_23", diffSCoinci_23, "diffSCoinci_23");
  plottingDiff("cvnsDiffSCoinci_31", diffSCoinci_31, "diffSCoinci_31");

  exit(0);
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

void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetLeftMargin(0.08);
  canvas.SetRightMargin(0.02);
}


void plottingDiff(TString cnvsName, TH1D &hist, TString outName) {
  TCanvas cvnsDiffS(cnvsName, "", 1.6e3, 9e2);
  setCanvasStyle(cvnsDiffS);

  cvnsDiffS.cd();
  cvnsDiffS.SetLogy();

  hist.Draw();
  cvnsDiffS.Print("results/"+outName+".pdf");
  
}
