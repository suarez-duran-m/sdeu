void setCanvasStyle(TCanvas& canvas);

void annaBinsLRQpk() {
  ifstream dataFile("binsLefltRightQpk.dat");

  const int totalHistos = 2464;
  int nbins = 3000;
  vector < double > xAxisFadcLR(nbins);
  vector < double > arrErrIQpk(nbins);
  vector < double > arrErrIQpkDev(nbins);
  vector < double > arrChi2Ndf(nbins);
  vector < double > arrChi2NdfDev(nbins);
  vector < double > xAxisNdof(nbins);
  vector < double > arrPval(nbins);
  vector < double > arrPvalDev(nbins);
  vector < double > nCntsInFadc(nbins);

  for ( int i=0; i<nbins; i++ ) {
    xAxisFadcLR[i] = 0.;
    arrErrIQpk[i] = 0.;
    arrErrIQpkDev[i] = 0.;
    arrChi2Ndf[i] = 0.;
    arrChi2NdfDev[i] = 0.;
    xAxisNdof[i] = 0.;
    arrPval[i] = 0.;
    nCntsInFadc[i] = 0;
  }

  double tmpFadc = 0.;
  double tmpIQpk = 0.;
  double tmpErrIQpk = 0.;
  double tmpChi2Ndf = 0.;
  double tmpNdof = 0.; 
  double tmpPval = 0.;

  double tmp = 0.;

  vector < double > valsChi2;
  vector < double > valsErrIpk;

  while ( dataFile >> tmpFadc >> tmpIQpk >> tmpErrIQpk >> tmpChi2Ndf >> tmpNdof ) {
    if ( !tmpNdof || !tmpChi2Ndf )
      continue;
    if ( tmpFadc < 9 )
      continue;
    
    tmpErrIQpk *= 100.;

    valsChi2.push_back( tmpChi2Ndf );
    valsErrIpk.push_back( tmpErrIQpk );

    arrErrIQpk[int(tmpFadc)] += tmpErrIQpk;
    arrErrIQpkDev[int(tmpFadc)] += tmpErrIQpk*tmpErrIQpk;

    arrChi2Ndf[int(tmpFadc)] += tmpChi2Ndf;
    arrChi2NdfDev[int(tmpFadc)] += tmpChi2Ndf*tmpChi2Ndf;

    tmpPval = log10( TMath::Prob(tmpChi2Ndf*tmpNdof, tmpNdof) );
    arrPval[int(tmpFadc)] += tmpPval;
    arrPvalDev[int(tmpFadc)] += tmpPval*tmpPval;

    xAxisFadcLR[int(tmpFadc)] = int(tmpFadc);
    xAxisNdof[int(tmpFadc)] = int(tmpNdof);
    nCntsInFadc[int(tmpFadc)]++;
  }
  dataFile.close();

  ifstream dataFileLeftRight("distLefltRightQpk.dat");
  TH1D distCntsLeftRight("distCntsLeftRight", "", 1.25e2, -1., 1.);
  distCntsLeftRight.SetTitle("; Counts (Left-Right)/Total counts [au]; counts [au]");
  double tmpCntLeftRight = 0;
  
  while ( dataFileLeftRight >> tmpCntLeftRight >> tmpCntLeftRight )
    distCntsLeftRight.Fill( tmpCntLeftRight );
  dataFileLeftRight.close();

  /*
  ifstream dataFileAfterFit("afterFitLefltRight.dat");
  TH1D distIpk("distIpk", "", 1500, 750, 3000);
  distIpk.SetTitle("; Ipk [FADC]; counts [au]");
  TH1D distChi2("distChi2", "", 65, 0, 5);
  distChi2.SetTitle("; #chi^{2}/Ndof; counts [au]");
  TH1D distErrIpk("distErrIpk", "", 100, 0, 10);
  distErrIpk.SetTitle("; #sigma/Ipk [%]; counts [au]");
  TH1D distLogPval("distLogPval", "", 100, 0, 10);
  distLogPval.SetTitle("; -1.Log(Pval); counts [au]");

  while ( dataFileAfterFit >> tmpIQpk >> tmpErrIQpk >> tmpChi2Ndf >> tmpNdof ) {
    distIpk.Fill( tmpIQpk );
    distChi2.Fill( tmpChi2Ndf );
    distErrIpk.Fill( 100.*tmpErrIQpk );
    distLogPval.Fill( -1.*log10( TMath::Prob(tmpChi2Ndf*tmpNdof, tmpNdof) ) );
  }
  dataFileAfterFit.close();
  */
  

  int nMinus1 = 0;
  for ( int i=0; i<nbins; i++ ) {
    if ( nCntsInFadc[i] == 0 )
      continue;

    nMinus1 = nCntsInFadc[i]-1;
    
    arrErrIQpk[i] = arrErrIQpk[i]/nCntsInFadc[i];
    arrErrIQpkDev[i] = sqrt(arrErrIQpkDev[i]/nMinus1 - nCntsInFadc[i]*arrErrIQpk[i]*arrErrIQpk[i]/nMinus1);
    arrErrIQpkDev[i] /= sqrt(nCntsInFadc[i]);

    arrChi2Ndf[i] = arrChi2Ndf[i]/nCntsInFadc[i];
    arrChi2NdfDev[i] = sqrt(arrChi2NdfDev[i]/nMinus1 - nCntsInFadc[i]*arrChi2Ndf[i]*arrChi2Ndf[i]/nMinus1);
    arrChi2NdfDev[i] /= sqrt(nCntsInFadc[i]);

    arrPval[i] = arrPval[i]/nCntsInFadc[i];
    arrPvalDev[i] = sqrt(arrPvalDev[i]/nMinus1 - nCntsInFadc[i]*arrPval[i]*arrPval[i]/nMinus1);
    arrPvalDev[i] /= sqrt(nCntsInFadc[i]);
    arrPval[i] *= -1.;

    nCntsInFadc[i] = 100.*(nCntsInFadc[i]/totalHistos);
  }

  double minNdof = 5.;
  double maxNdof = 250.;

  TGraphErrors errIQpkgrph(xAxisNdof.size(), &xAxisNdof.front(), &arrErrIQpk.front(), 0, &arrErrIQpkDev.front());
  errIQpkgrph.SetTitle("; Ndof [au]; #sigma/Qpk [%]");
  
  TCanvas cvnsErrIQpk("cvnsErrIQpk","",1.6e3, 9e2);
  setCanvasStyle(cvnsErrIQpk);
  cvnsErrIQpk.cd();
  cvnsErrIQpk.SetLogy(); 

  errIQpkgrph.GetXaxis()->SetRangeUser(minNdof, maxNdof);
  errIQpkgrph.GetYaxis()->SetRangeUser(0.3, 5);
  errIQpkgrph.SetMarkerStyle(20);
  errIQpkgrph.SetMarkerSize(2);
  errIQpkgrph.SetMarkerColor(kBlue);
  errIQpkgrph.SetLineColor(kBlue);
  errIQpkgrph.Draw("ap");

  TLine *line = new TLine(minNdof,3, maxNdof,3);
  line->SetLineStyle(2);
  line->SetLineColor(kGray+2);
  line->Draw();

  cvnsErrIQpk.Print("results/resLeftRight_ErrQpk.pdf");


  TGraphErrors redChi2grph(xAxisNdof.size(), &xAxisNdof.front(), &arrChi2Ndf.front(), 0, &arrChi2NdfDev.front());
  redChi2grph.SetTitle("; Ndof [au]; #chi^{2}/Ndof [FADC]");
  
  TCanvas cvnsRedChi2("cvnsRedChi2","",1.6e3, 9e2);
  setCanvasStyle(cvnsRedChi2);
  cvnsRedChi2.cd();

  redChi2grph.GetXaxis()->SetRangeUser(minNdof, maxNdof);
  redChi2grph.GetYaxis()->SetRangeUser(0.8, 2.8);
  redChi2grph.SetMarkerStyle(20);
  redChi2grph.SetMarkerSize(2);
  redChi2grph.SetMarkerColor(kBlue);
  redChi2grph.SetLineColor(kBlue);
  redChi2grph.Draw("ap");

  line = new TLine(minNdof,2, maxNdof,2);
  line->SetLineStyle(2);
  line->SetLineColor(kGray+2);
  line->Draw();

  cvnsRedChi2.Print("results/resLeftRight_redChi2Qpk.pdf");


  TGraphErrors logPvalgrph(xAxisNdof.size(), &xAxisNdof.front(), &arrPval.front(), 0, &arrPvalDev.front());
  logPvalgrph.SetTitle("; Ndof [au]; (-1)*Log(Pval) [au]");
  
  TCanvas cvnsLogPval("cvnsLogPval","",1.6e3, 9e2);
  setCanvasStyle(cvnsLogPval);
  cvnsLogPval.cd();
  cvnsLogPval.SetLogy();

  logPvalgrph.GetXaxis()->SetRangeUser(minNdof, maxNdof);
  logPvalgrph.GetYaxis()->SetRangeUser(.3, 5);
  logPvalgrph.SetMarkerStyle(20);
  logPvalgrph.SetMarkerSize(2);
  logPvalgrph.SetMarkerColor(kBlue);
  logPvalgrph.SetLineColor(kBlue);
  logPvalgrph.Draw("ap");

  line = new TLine(minNdof,1, maxNdof,1);
  line->SetLineStyle(2);
  line->SetLineColor(kGray+2);
  line->Draw();

  cvnsLogPval.Print("results/resLeftRight_logPvalQpk.pdf");


  TGraph succesHistos(xAxisNdof.size(), &xAxisNdof.front(), &nCntsInFadc.front());
  succesHistos.SetTitle("; Ndof [au]; Fit_{Success}/TotalHistos [%]");

  TCanvas cvnsSuccHistos("cvnsSuccHistos","",1.6e3, 9e2);
  setCanvasStyle(cvnsSuccHistos);
  cvnsSuccHistos.cd();

  succesHistos.GetXaxis()->SetRangeUser(minNdof, maxNdof);
  //succesHistos.GetYaxis()->SetRangeUser(76, 96);
  succesHistos.SetMarkerStyle(20);
  succesHistos.SetMarkerSize(2);
  succesHistos.SetMarkerColor(kBlue);
  succesHistos.SetLineColor(kBlue);
  succesHistos.Draw("ap");

  cvnsSuccHistos.Print("resLeftRight_succesHistsoQpk.pdf");

/*
  TGraph errIpkVsChi2(valsChi2.size(), &valsChi2.front(), &valsErrIpk.front());  
  TGraph errIpkVsChi2Ave(arrErrIQpk.size(), &arrChi2Ndf.front(), &arrErrIQpk.front()); 
  errIpkVsChi2.SetTitle("; #chi^{2}/Ndof [FADC]; #sigma/Qpk [%]");
  
  TCanvas cvnsErrIpkVsChi2("cvnsErrIpkVsChi2","",1.6e3, 9e2);
  setCanvasStyle(cvnsErrIpkVsChi2);
  cvnsErrIpkVsChi2.cd();
  cvnsErrIpkVsChi2.SetLogy();

  errIpkVsChi2.SetMarkerStyle(20);
  errIpkVsChi2.SetMarkerSize(1);
  errIpkVsChi2.SetMarkerColor(kBlue);
  errIpkVsChi2.Draw("ap");

  errIpkVsChi2Ave.SetMarkerStyle(20);
  errIpkVsChi2Ave.SetMarkerSize(2);
  errIpkVsChi2Ave.SetMarkerColor(kRed);
  errIpkVsChi2Ave.Draw("p same");

  cvnsErrIpkVsChi2.Print("resLeftRight_ErrIpkVsChi2Qpk.pdf");
*/

  TCanvas cvnsDistCntLeft("cvnsDistCntLeft","",1.6e3, 9e2);
  setCanvasStyle(cvnsDistCntLeft);
  cvnsDistCntLeft.cd();
  cvnsDistCntLeft.SetLogy();

  distCntsLeftRight.SetStats(kFALSE);
  distCntsLeftRight.SetLineColor(kBlue);
  distCntsLeftRight.Draw();

  TLegend *lgnd = new TLegend(0.7, 0.75, 0.9, 0.94);
  lgnd->AddEntry(&distCntsLeftRight,Form("Mean: %.2f #pm %.2f",
      distCntsLeftRight.GetMean(), distCntsLeftRight.GetMeanError()), "l");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.04);
  lgnd->Draw();

  cvnsDistCntLeft.Print("results/resLeftRight_DistLeftRightQpk.pdf");

/*
  TCanvas cvnsAfterDistIpk("cvnsAfterDistIpk","",1.6e3, 9e2);
  setCanvasStyle(cvnsAfterDistIpk);
  cvnsAfterDistIpk.cd();
  cvnsAfterDistIpk.SetLogy();

  distIpk.Draw();
  cvnsAfterDistIpk.Print("results/resLeftRight_DistQpk.pdf");


  TCanvas cvnsAfterDistErrIpk("cvnsAfterDistErrIpk","",1.6e3, 9e2);
  setCanvasStyle(cvnsAfterDistErrIpk);
  cvnsAfterDistErrIpk.cd();
  cvnsAfterDistErrIpk.SetLogy();

  distErrIpk.Draw();
  cvnsAfterDistErrIpk.Print("results/resLeftRight_DistErrQpk.pdf");
  
  
  TCanvas cvnsAfterDistChi2("cvnsAfterDistChi2","",1.6e3, 9e2);
  setCanvasStyle(cvnsAfterDistChi2);
  cvnsAfterDistChi2.cd();
  cvnsAfterDistChi2.SetLogy();
  
  distChi2.Draw();
  cvnsAfterDistChi2.Print("results/resLeftRight_DistChi2Qpk.pdf");

  
  TCanvas cvnsAfterDistPval("cvnsAfterDistPval","",1.6e3, 9e2);
  setCanvasStyle(cvnsAfterDistPval);
  cvnsAfterDistPval.cd();
  cvnsAfterDistPval.SetLogy();
  
  distLogPval.Draw();
  cvnsAfterDistPval.Print("results/resLeftRight_DisPvallQpk.pdf");
*/
  exit(0);
}


void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetLeftMargin(0.08);
  canvas.SetRightMargin(0.02);
}
