void doingChi2Ndof() {
  ifstream dataFile("ndfVsRedChi2_Ipk.dat", ios::in);
  //ifstream dataFile("ndfVsRedChi2_Qpk.dat", ios::in);

  int nBinLow = 400;
  int nNdof = 100;
  // redChi2[binlow][ndof] = redChi2
  vector < vector < double > > redChi2(nBinLow);
  vector < vector < double > > errIQpk(nBinLow);
  vector < vector < double > > pval(nBinLow);
  vector < vector < double > > cntsPerBin(nBinLow);
  vector < double > xAxisNdof(nNdof);

  for ( int i=0; i<nBinLow; i++ ) {
    redChi2[i].resize( nNdof );
    errIQpk[i].resize( nNdof );
    pval[i].resize( nNdof );
    cntsPerBin[i].resize( nNdof );

    for ( int j=0; j<nNdof; j++ ) {
      redChi2[i][j] = 0.;
      errIQpk[i][j] = 0.;
      pval[i][j] = 0.;
      cntsPerBin[i][j] = 0.;
      xAxisNdof[j] = j;
    }
  }

  double tmpNdof = 0.;
  double tmpRedChi2 = 0.;
  double tmpErrIQpk = 0.;
  double tmpLowBin = 0.;
  double tmpHigBin = 0.;

  while (dataFile.good()) {
    dataFile >> tmpNdof >> tmpRedChi2 >> tmpErrIQpk >> tmpLowBin >> tmpHigBin;
    redChi2[tmpLowBin][tmpNdof] += tmpRedChi2;
    errIQpk[tmpLowBin][tmpNdof] += tmpErrIQpk;
    pval[tmpLowBin][tmpNdof] += log10( TMath::Prob(tmpRedChi2*tmpNdof, tmpNdof) );    
    cntsPerBin[tmpLowBin][tmpNdof]++;
  }

  for ( int i=0; i<nBinLow; i++ )
    for ( int j=0; j<nNdof; j++ )
      if ( cntsPerBin[i][j] > 0 ) {
        redChi2[i][j] /= cntsPerBin[i][j];
        errIQpk[i][j] /= cntsPerBin[i][j];
        errIQpk[i][j] *= 100.;
      }

  TGraph redChi2Grph100(xAxisNdof.size(), &xAxisNdof.front(), &redChi2[100].front());
  redChi2Grph100.SetTitle("; Ndof [au]; #chi^{2}/Ndof [FADC]");
  
  TGraph redChi2Grph125(xAxisNdof.size(), &xAxisNdof.front(), &redChi2[124].front());
  TGraph redChi2Grph150(xAxisNdof.size(), &xAxisNdof.front(), &redChi2[152].front());
  TGraph redChi2Grph175(xAxisNdof.size(), &xAxisNdof.front(), &redChi2[176].front());
  TGraph redChi2Grph200(xAxisNdof.size(), &xAxisNdof.front(), &redChi2[200].front());

  TCanvas cvnsRedChi2("cvnsRedChi2","",1.6e3, 9e2);
  cvnsRedChi2.cd();
  cvnsRedChi2.SetTopMargin(0.03);
  cvnsRedChi2.SetLeftMargin(0.08);
  cvnsRedChi2.SetRightMargin(0.02);

  redChi2Grph100.GetXaxis()->SetRangeUser(0, 74);
  redChi2Grph100.GetYaxis()->SetRangeUser(0.8, 4.5);
  redChi2Grph100.SetMarkerStyle(20);
  redChi2Grph100.SetMarkerSize(2);
  redChi2Grph100.SetMarkerColor(kBlue);
  redChi2Grph100.Draw("ap");

  redChi2Grph125.SetMarkerStyle(21);
  redChi2Grph125.SetMarkerSize(2);
  redChi2Grph125.SetMarkerColor(kYellow+2);
  redChi2Grph125.Draw("p same");

  redChi2Grph150.SetMarkerStyle(22);
  redChi2Grph150.SetMarkerSize(2.5);
  redChi2Grph150.SetMarkerColor(kGreen+2);
  redChi2Grph150.Draw("p same");

  redChi2Grph175.SetMarkerStyle(23);
  redChi2Grph175.SetMarkerSize(2.5);
  redChi2Grph175.SetMarkerColor(kRed);
  redChi2Grph175.Draw("p same");

  redChi2Grph200.SetMarkerStyle(33);
  redChi2Grph200.SetMarkerSize(2.5);
  redChi2Grph200.SetMarkerColor(kGray+3);
  redChi2Grph200.Draw("p same");

  TLegend *lgn = new TLegend(0.12, 0.75, 0.35, 0.92);
  lgn->AddEntry(&redChi2Grph100,"binLow = 100 [FADC]", "p");
  lgn->AddEntry(&redChi2Grph125,"binLow = 125 [FADC]", "p");
  lgn->AddEntry(&redChi2Grph150,"binLow = 150 [FADC]", "p");
  lgn->AddEntry(&redChi2Grph175,"binLow = 175 [FADC]", "p");
  lgn->AddEntry(&redChi2Grph200,"binLow = 200 [FADC]", "p");
  lgn->SetBorderSize(0);
  lgn->SetTextSize(0.04);
  lgn->Draw();
  
  cvnsRedChi2.Print("redChi2Ipk.pdf");


  // *********************************
  // =========== ErrIQpk =============

  TGraph errIQpkGrph100(xAxisNdof.size(), &xAxisNdof.front(), &errIQpk[100].front());
  errIQpkGrph100.SetTitle("; Ndof [au]; #sigma/Ipk [%]");
  
  TGraph errIQpkGrph125(xAxisNdof.size(), &xAxisNdof.front(), &errIQpk[124].front());
  TGraph errIQpkGrph150(xAxisNdof.size(), &xAxisNdof.front(), &errIQpk[152].front());
  TGraph errIQpkGrph175(xAxisNdof.size(), &xAxisNdof.front(), &errIQpk[176].front());
  TGraph errIQpkGrph200(xAxisNdof.size(), &xAxisNdof.front(), &errIQpk[200].front());

  TCanvas cvnsErrIQpk("cvnsErrIQpk","",1.6e3, 9e2);
  cvnsErrIQpk.cd();
  cvnsErrIQpk.SetTopMargin(0.03);
  cvnsErrIQpk.SetLeftMargin(0.08);
  cvnsErrIQpk.SetRightMargin(0.02);

  errIQpkGrph100.GetXaxis()->SetRangeUser(0, 74);
  errIQpkGrph100.GetYaxis()->SetRangeUser(0.5, 12);
  errIQpkGrph100.SetMarkerStyle(20);
  errIQpkGrph100.SetMarkerSize(2);
  errIQpkGrph100.SetMarkerColor(kBlue);
  errIQpkGrph100.Draw("ap");

  errIQpkGrph125.SetMarkerStyle(21);
  errIQpkGrph125.SetMarkerSize(2);
  errIQpkGrph125.SetMarkerColor(kYellow+2);
  errIQpkGrph125.Draw("p same");

  errIQpkGrph150.SetMarkerStyle(22);
  errIQpkGrph150.SetMarkerSize(2.5);
  errIQpkGrph150.SetMarkerColor(kGreen+2);
  errIQpkGrph150.Draw("p same");

  errIQpkGrph175.SetMarkerStyle(23);
  errIQpkGrph175.SetMarkerSize(2.5);
  errIQpkGrph175.SetMarkerColor(kRed);
  errIQpkGrph175.Draw("p same");

  errIQpkGrph200.SetMarkerStyle(33);
  errIQpkGrph200.SetMarkerSize(2.5);
  errIQpkGrph200.SetMarkerColor(kGray+3);
  errIQpkGrph200.Draw("p same");

  lgn = new TLegend(0.65, 0.75, 0.88, 0.92);
  lgn->AddEntry(&errIQpkGrph100,"binLow = 100 [FADC]", "p");
  lgn->AddEntry(&errIQpkGrph125,"binLow = 125 [FADC]", "p");
  lgn->AddEntry(&errIQpkGrph150,"binLow = 150 [FADC]", "p");
  lgn->AddEntry(&errIQpkGrph175,"binLow = 175 [FADC]", "p");
  lgn->AddEntry(&errIQpkGrph200,"binLow = 200 [FADC]", "p");
  lgn->SetBorderSize(0);
  lgn->SetTextSize(0.04);
  lgn->Draw();
  
  cvnsErrIQpk.Print("errIQpkIpk.pdf");

  exit(0);
}
