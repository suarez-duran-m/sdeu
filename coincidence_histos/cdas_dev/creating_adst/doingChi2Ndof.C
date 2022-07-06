void doingChi2Ndof() {
  ifstream dataFile("ndfVsRedChi2_Ipk.dat", ios::in);
  //ifstream dataFile("ndfVsRedChi2_Qpk.dat", ios::in);

  int nOfBins = 400;
  vector < vector < double > > lowHigChi2(nOfBins);
  vector < vector < double > > cntsLowHig(nOfBins);
  vector < vector < double > > arrErrIpk(nOfBins);
  vector < vector < double > > arrPval(nOfBins);

  for ( int i=0; i<nOfBins; i++ ) {
    lowHigChi2[i].resize(nOfBins);
    arrErrIpk[i].resize(nOfBins);
    cntsLowHig[i].resize(nOfBins);
    arrPval[i].resize(nOfBins);

    for ( int j=0; j<nOfBins; j++ ) {
      lowHigChi2[i][j] = 0.;
      arrErrIpk[i][j] = 0.;
      cntsLowHig[i][j] = 0.;
      arrPval[i][j] = 0.;
    }
  }

  double ndof = 0.;
  double chi2ndof = 0.;
  double ipkErr = 0.;
  double lowBin = 0.;
  double higBin = 0.;

  TGraph2D *lowHigChi2Ndof = new TGraph2D();
  TGraph2D *lowHigErrIpk = new TGraph2D();
  TGraph2D *lowHigPval = new TGraph2D();
  lowHigChi2Ndof->SetTitle("; binLow [FADC]; binHigh [FADC]; #chi^{2}/Ndof [FADC]"); 
  lowHigPval->SetTitle("; binLow [FADC]; binHigh [FADC]; Log(Pval)");
  lowHigErrIpk->SetTitle("; binLow [FADC]; binHigh [FADC]; ErrIpk [%]"); 
  //lowHigErrIpk->SetTitle("; binLow [FADC]; binHigh [FADC]; ErrQpk [%]"); 

  while (dataFile.good()) {
    dataFile >> ndof >> chi2ndof >> ipkErr >> lowBin >> higBin;
      lowHigChi2[lowBin][higBin] += chi2ndof;
      arrErrIpk[lowBin][higBin] += ipkErr;
      arrPval[lowBin][higBin] += log10( TMath::Prob(chi2ndof*ndof, ndof) );
      cntsLowHig[lowBin][higBin]++;
  }
 
  int nPoint = 0;
  int nPointPval = 0;

  for ( int i=0; i<nOfBins; i++ )
    for ( int j=0; j<nOfBins; j++ )
      if ( cntsLowHig[i][j] > 0 ) {
        lowHigChi2[i][j] /= cntsLowHig[i][j];
        lowHigChi2Ndof->SetPoint( nPoint, i, j, lowHigChi2[i][j] );
        //lowHigChi2Ndof->SetPoint( nPoint, i*8, j*8, lowHigChi2[i][j] );
        arrErrIpk[i][j] /= cntsLowHig[i][j];
        lowHigErrIpk->SetPoint( nPoint, i, j, 100.*arrErrIpk[i][j] );
        //lowHigErrIpk->SetPoint( nPoint, i*8, j*8, 100.*arrErrIpk[i][j] );
        arrPval[i][j] /= cntsLowHig[i][j];
        if ( arrPval[i][j] > -10 && arrPval[i][j] < -0.001 ) {
          lowHigPval->SetPoint( nPointPval, i, j, arrPval[i][j] );
          //lowHigPval->SetPoint( nPointPval, i*8, j*8, arrPval[i][j] );
          nPointPval++;
        }
        nPoint++;
      }
  int usePalette = 56;

  TCanvas c2dChi2("c2dChi2","",1.6e3, 9e2);
  c2dChi2.cd();
  c2dChi2.SetTopMargin(0.03);
  c2dChi2.SetLeftMargin(0.08);
  c2dChi2.SetRightMargin(0.13);
  gStyle->SetPalette(usePalette);

  lowHigChi2Ndof->Draw("COLZ CONT");
  c2dChi2.Print("lowHigChi2NdofIpk.pdf");
  //c2dChi2.Print("lowHigChi2NdofQpk.pdf");

  TCanvas c2dErrIpk("c2dErrIpk","",1.6e3, 9e2);
  c2dErrIpk.cd();
  c2dErrIpk.SetTopMargin(0.03);
  c2dErrIpk.SetLeftMargin(0.08);
  c2dErrIpk.SetRightMargin(0.13);
  gStyle->SetPalette(usePalette);

  lowHigErrIpk->Draw("COLZ CONT");
  c2dErrIpk.Print("lowHigErrIpk.pdf");
  //c2dErrIpk.Print("lowHigErrQpk.pdf");


  usePalette = 53;
  TCanvas c2dPval("c2dPval","",1.6e3, 9e2);
  c2dPval.cd();
  c2dPval.SetTopMargin(0.03);
  c2dPval.SetLeftMargin(0.08);
  c2dPval.SetRightMargin(0.13);
  gStyle->SetPalette(usePalette);
 
  lowHigPval->Draw("COLZ CONT");
  c2dPval.Print("lowHigPvalIpk.pdf");
  //c2dPval.Print("lowHigPvalQpk.pdf");

  exit(0);
}
