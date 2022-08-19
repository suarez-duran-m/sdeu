void setCanvasStyle(TCanvas& canvas);
void setTGraphStyle(TGraphErrors& graph);
void setTGraphStyle(TGraph& graph);
void setTGraphStyle(TH1D& graph);
void setTGraphStyle(TH2D& graph);

void annaBinsLRIpk() {
  ifstream dataFile("binsLefltRightIpk.dat");

  const int totalHistos = 2198; //2464;
  //int nbins = 152;
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

  double tmpNfadc = 0.;
  double tmpIpk = 0.;
  double tmpErrIpk = 0.;
  double tmpChi2Ndf = 0.;
  double tmpNdof = 0.; 
  double tmpPval = 0.;

  vector < double > valsChi2;
  vector < double > valsErrIpk;

  while ( dataFile >> tmpNfadc >> tmpIpk >> tmpErrIpk >> tmpChi2Ndf >> tmpNdof ) {
    if ( !tmpNdof || !tmpChi2Ndf )
      continue;
    if ( tmpNfadc < 9 )
      continue;
    
    tmpErrIpk *= 100.;

    valsChi2.push_back( tmpChi2Ndf );
    valsErrIpk.push_back( tmpErrIpk );

    arrErrIQpk[int(tmpNfadc)] += tmpErrIpk;
    arrErrIQpkDev[int(tmpNfadc)] += tmpErrIpk*tmpErrIpk;

    arrChi2Ndf[int(tmpNfadc)] += tmpChi2Ndf;
    arrChi2NdfDev[int(tmpNfadc)] += tmpChi2Ndf*tmpChi2Ndf;

    tmpPval = log10( TMath::Prob(tmpChi2Ndf*tmpNdof, tmpNdof) );
    arrPval[int(tmpNfadc)] += tmpPval;
    arrPvalDev[int(tmpNfadc)] += tmpPval*tmpPval;

    xAxisFadcLR[int(tmpNfadc)] = int(tmpNfadc);
    xAxisNdof[int(tmpNfadc)] = int(tmpNdof);
    nCntsInFadc[int(tmpNfadc)]++;
  }
  dataFile.close();  

  ifstream dataFileLeftRight("distLefltRightIpk.dat");
  TH1D distCntsLeftRight("distCntsLeftRight", "", 1.25e2, -1., 1.);
  distCntsLeftRight.SetTitle("; Counts (Left-Right)/Total [au]; counts [au]");
  double tmpCntLeftRight = 0;
  
  while ( dataFileLeftRight >> tmpCntLeftRight )
    distCntsLeftRight.Fill( tmpCntLeftRight );
  dataFileLeftRight.close();

  ifstream dataFileAfterFit("IpkVsCoinciIpk.dat");
  TH2D distIpk("distIpk", "", 75, 100, 400, 75, 100, 400);
  distIpk.SetTitle("; Ipk [FADC]; CIpk [FADC]; Counts [au]");
  TH2D distErrIpk("distErrIpk", "", 500, 0, 100, 500, 0, 100);
  distErrIpk.SetTitle("; #sigma/Ipk [%]; #sigma/CIpk [%]; Counts [au]");

  double tmpCIpk = 0.;
  double tmpErrCIpk = 0.;

  while ( dataFileAfterFit >> tmpIpk >> tmpErrIpk >> tmpCIpk >> tmpErrCIpk ) {
    if ( tmpIpk < 1 || tmpCIpk < 1 )
      continue;
    distIpk.Fill( tmpIpk, tmpCIpk );
    distErrIpk.Fill( 100.*tmpErrIpk/tmpIpk, 100.*tmpErrCIpk/tmpCIpk );
  }
  dataFileAfterFit.close();  

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

  double minNdof = 1.;
  double maxNdof = 48.;

  TGraphErrors errIQpkgrph(xAxisNdof.size(), &xAxisNdof.front(), &arrErrIQpk.front(), 0, &arrErrIQpkDev.front());
  errIQpkgrph.SetTitle("; Ndof [au]; #sigma/Ipk [%]");
  
  TCanvas cvnsErrIQpk("cvnsErrIQpk","",1.6e3, 9e2);
  setCanvasStyle(cvnsErrIQpk);
  cvnsErrIQpk.cd();
  cvnsErrIQpk.SetLogy();

  errIQpkgrph.GetXaxis()->SetRangeUser(minNdof, maxNdof);
  errIQpkgrph.GetYaxis()->SetRangeUser(0.5, 6);
  errIQpkgrph.SetMarkerStyle(20);
  errIQpkgrph.SetMarkerSize(2);
  errIQpkgrph.SetMarkerColor(kBlue);
  errIQpkgrph.SetLineColor(kBlue);
  setTGraphStyle(errIQpkgrph);
  errIQpkgrph.GetYaxis()->SetTitleOffset(0.7);
  errIQpkgrph.Draw("ap");

  cvnsErrIQpk.Print("results/resLeftRight_ErrIpk.pdf");


  TGraphErrors redChi2grph(xAxisNdof.size(), &xAxisNdof.front(), &arrChi2Ndf.front(), 0, &arrChi2NdfDev.front());
  redChi2grph.SetTitle("; Ndof [au]; #chi^{2}/Ndof [FADC]");
  
  TCanvas cvnsRedChi2("cvnsRedChi2","",1.6e3, 9e2);
  setCanvasStyle(cvnsRedChi2);
  cvnsRedChi2.cd();

  redChi2grph.GetXaxis()->SetRangeUser(minNdof, maxNdof);
  redChi2grph.GetYaxis()->SetRangeUser(1.05, 1.77);
  redChi2grph.SetMarkerStyle(20);
  redChi2grph.SetMarkerSize(2);
  redChi2grph.SetMarkerColor(kBlue);
  redChi2grph.SetLineColor(kBlue);
  setTGraphStyle(redChi2grph);
  redChi2grph.GetYaxis()->SetTitleOffset(0.8);
  redChi2grph.Draw("ap");

  cvnsRedChi2.Print("results/resLeftRight_redChi2Ipk.pdf");


  TGraphErrors logPvalgrph(xAxisNdof.size(), &xAxisNdof.front(), &arrPval.front(), 0, &arrPvalDev.front());
  logPvalgrph.SetTitle("; Ndof [au]; (-1)*Log(Pval) [au]");
  
  TCanvas cvnsLogPval("cvnsLogPval","",1.6e3, 9e2);
  setCanvasStyle(cvnsLogPval);
  cvnsLogPval.cd();
  cvnsLogPval.SetLogy();

  logPvalgrph.GetXaxis()->SetRangeUser(minNdof, maxNdof);
  logPvalgrph.GetYaxis()->SetRangeUser(.5, 3.7);
  logPvalgrph.SetMarkerStyle(20);
  logPvalgrph.SetMarkerSize(2);
  logPvalgrph.SetMarkerColor(kBlue);
  logPvalgrph.SetLineColor(kBlue);
  setTGraphStyle(logPvalgrph);
  logPvalgrph.GetYaxis()->SetTitleOffset(0.7);
  logPvalgrph.Draw("ap");

  cvnsLogPval.Print("results/resLeftRight_logPvalIpk.pdf");

  TGraph succesHistos(xAxisNdof.size(), &xAxisNdof.front(), &nCntsInFadc.front());
  succesHistos.SetTitle("; Ndof [au]; Fit_{Success}/TotalHistos [%]");

  TCanvas cvnsSuccHistos("cvnsSuccHistos","",1.6e3, 9e2);
  setCanvasStyle(cvnsSuccHistos);
  cvnsSuccHistos.cd();

  succesHistos.GetXaxis()->SetRangeUser(minNdof, maxNdof);
  succesHistos.GetYaxis()->SetRangeUser(82, 100);
  succesHistos.SetMarkerStyle(20);
  succesHistos.SetMarkerSize(2);
  succesHistos.SetMarkerColor(kBlue);
  succesHistos.SetLineColor(kBlue);
  setTGraphStyle(succesHistos);
  succesHistos.GetYaxis()->SetTitleOffset(0.7);
  succesHistos.Draw("ap");

  cvnsSuccHistos.Print("results/resLeftRight_succesHistsoIpk.pdf");


  TCanvas cvnsDistCntLeft("cvnsDistCntLeft","",1.6e3, 9e2);
  setCanvasStyle(cvnsDistCntLeft);
  cvnsDistCntLeft.cd();
  cvnsDistCntLeft.SetLogy();

  distCntsLeftRight.SetStats(kFALSE);
  distCntsLeftRight.SetLineColor(kBlue);
  distCntsLeftRight.Fit("gaus","QR");
  setTGraphStyle(distCntsLeftRight);
  distCntsLeftRight.GetYaxis()->SetTitleOffset(0.8);
  distCntsLeftRight.Draw();

  double cut = distCntsLeftRight.GetFunction("gaus")->GetParameter(1) 
    + 3*distCntsLeftRight.GetFunction("gaus")->GetParameter(2);

  line = new TLine(cut, 0, cut, 1e2);
  line->SetLineWidth(2);
  line->Draw();

  TLegend *lgnd = new TLegend(0.7, 0.75, 0.9, 0.94);
  lgnd->AddEntry(&distCntsLeftRight,Form("Mean: %.3f #pm %.3f",
      distCntsLeftRight.GetMean(), distCntsLeftRight.GetMeanError()), "l");
  lgnd->AddEntry(&distCntsLeftRight.GetFunction("gaus"), Form("#mu = %.3f; #sigma = %.3f",
        distCntsLeftRight.GetFunction("gaus")->GetParameter(1),
        distCntsLeftRight.GetFunction("gaus")->GetParameter(2)),"l");
  lgnd->AddEntry(line, Form("Cut at #mu+3#sigma: %.3f", cut), "l");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDistCntLeft.Print("results/resLeftRight_DistLeftRightIpk.pdf");


  TCanvas cvnsIpkVsCIpk("cvnsIpkVsCIpk","",1.6e3, 9e2);
  cvnsIpkVsCIpk.SetTopMargin(0.03);
  cvnsIpkVsCIpk.SetBottomMargin(0.12);
  cvnsIpkVsCIpk.SetLeftMargin(0.09);
  cvnsIpkVsCIpk.SetRightMargin(0.12);
  cvnsIpkVsCIpk.cd();

  distIpk.Fit("pol1", "Q", "R", 125, 260);
  distIpk.SetStats(kFALSE);
  distIpk.GetYaxis()->SetTitleOffset(0.85);
  distIpk.GetZaxis()->SetTitleOffset(0.65);
  distIpk.GetYaxis()->SetRangeUser(120,320);
  distIpk.GetXaxis()->SetRangeUser(120,320);
  distIpk.SetMarkerStyle(20);
  distIpk.SetMarkerColor(kBlue);
  distIpk.SetMarkerSize(2);
  distIpk.GetFunction("pol1")->SetLineWidth(4);
  setTGraphStyle(distIpk);
  distIpk.Draw("colz");

  lgnd = new TLegend(0.6, 0.16, 0.8, 0.7);
  lgnd->AddEntry(&distIpk, "Ipk", "");
  lgnd->AddEntry(&distIpk, Form("Total entries: %.f",
        distIpk.GetEntries()), ""); 
  lgnd->AddEntry(&distIpk, Form("Mean: %.2f #pm %.2f",
        distIpk.ProfileX()->GetMean(), distIpk.ProfileX()->GetMeanError()), "");
  lgnd->AddEntry(&distIpk, Form("RMS: %.2f #pm %.2f",
        distIpk.ProfileX()->GetRMS(), distIpk.ProfileX()->GetRMSError()), "");
  lgnd->AddEntry(&distIpk, "CIpk", "");
  lgnd->AddEntry(&distIpk, Form("Mean: %.2f #pm %.2f",
        distIpk.ProfileY()->GetMean(), distIpk.ProfileY()->GetMeanError()), "");
  lgnd->AddEntry(&distIpk, Form("RMS: %.2f #pm %.2f",
        distIpk.ProfileY()->GetRMS(), distIpk.ProfileY()->GetRMSError()), "");
  lgnd->AddEntry(&distIpk.GetFunction("pol1"), "Fit", "");
  lgnd->AddEntry(&distIpk.GetFunction("pol1"), Form("Slope = %.2f #pm %.2f",
        distIpk.GetFunction("pol1")->GetParameter(1), 
        distIpk.GetFunction("pol1")->GetParError(1)), "");
  lgnd->AddEntry(&distIpk.GetFunction("pol1"), Form("b = %.2f #pm %.2f",
        distIpk.GetFunction("pol1")->GetParameter(0),
        distIpk.GetFunction("pol1")->GetParError(0)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();
/*
  lgnd = new TLegend(0.7, 0.15, 0.8, 0.3);
  lgnd->AddEntry(&distIpk.GetFunction("pol1"), Form("slope = %.3f", 
        distIpk.GetFunction("pol1")->GetParameter(1)), "l");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();
*/
  cvnsIpkVsCIpk.Print("results/resLeftRight_IpkVsCIpk.pdf");


  TCanvas cvnsErrIpkVsCIpk("cvnsErrIpkVsCIpk","",1.6e3, 9e2);
  cvnsErrIpkVsCIpk.SetTopMargin(0.03);
  cvnsErrIpkVsCIpk.SetBottomMargin(0.12);
  cvnsErrIpkVsCIpk.SetLeftMargin(0.07);
  cvnsErrIpkVsCIpk.SetRightMargin(0.14);
  cvnsErrIpkVsCIpk.cd();

  distErrIpk.SetStats(kFALSE);
  distErrIpk.GetYaxis()->SetTitleOffset(0.6);
  distErrIpk.GetZaxis()->SetTitleOffset(0.85);
  distErrIpk.GetYaxis()->SetRangeUser(0.5,6);
  distErrIpk.GetXaxis()->SetRangeUser(0.5,100);
  distErrIpk.SetMarkerStyle(20);
  distErrIpk.SetMarkerColor(kBlue);
  distErrIpk.SetMarkerSize(2);
  setTGraphStyle(distErrIpk);
  distErrIpk.Draw("colz");

  lgnd = new TLegend(0.6, 0.55, 0.8, 0.95);
  lgnd->AddEntry(&distErrIpk, "ErrIpk", "");
  lgnd->AddEntry(&distErrIpk, Form("Total entries: %.f",
        distErrIpk.GetEntries()), "");
  lgnd->AddEntry(&distErrIpk, Form("Mean: %.2f #pm %.2f",
        distErrIpk.ProfileX()->GetMean(), distErrIpk.ProfileX()->GetMeanError()), "");
  lgnd->AddEntry(&distErrIpk, Form("RMS: %.2f #pm %.2f",
        distErrIpk.ProfileX()->GetRMS(), distErrIpk.ProfileX()->GetRMSError()), "");
  lgnd->AddEntry(&distErrIpk, "ErrCIpk", "");
  lgnd->AddEntry(&distErrIpk, Form("Mean: %.2f #pm %.2f",
        distErrIpk.ProfileY()->GetMean(), distErrIpk.ProfileY()->GetMeanError()), "");
  lgnd->AddEntry(&distErrIpk, Form("RMS: %.2f #pm %.2f",
        distErrIpk.ProfileY()->GetRMS(), distErrIpk.ProfileY()->GetRMSError()), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsErrIpkVsCIpk.Print("results/resLeftRight_ErrIpkVsCIpk.pdf");


  TCanvas cvnsZErrIpkVsCIpk("cvnsZErrIpkVsCIpk","",1.6e3, 9e2);
  cvnsZErrIpkVsCIpk.SetTopMargin(0.03);
  cvnsZErrIpkVsCIpk.SetBottomMargin(0.12);
  cvnsZErrIpkVsCIpk.SetLeftMargin(0.075);
  cvnsZErrIpkVsCIpk.SetRightMargin(0.14);
  cvnsZErrIpkVsCIpk.cd();

  distErrIpk.SetStats(kFALSE);
  distErrIpk.GetYaxis()->SetTitleOffset(0.7);
  distErrIpk.GetZaxis()->SetTitleOffset(0.85);
  distErrIpk.GetYaxis()->SetRangeUser(0.5,5);
  distErrIpk.GetXaxis()->SetRangeUser(0.5,5);
  distErrIpk.SetMarkerStyle(20);
  distErrIpk.SetMarkerColor(kBlue);
  distErrIpk.SetMarkerSize(2);
  setTGraphStyle(distErrIpk);
  distErrIpk.Draw("colz");

  lgnd = new TLegend(0.6, 0.55, 0.8, 0.95);
  lgnd->AddEntry(&distErrIpk, "ErrIpk", "");
  lgnd->AddEntry(&distErrIpk, Form("Total entries: %.f",
        distErrIpk.GetEntries()), "");
  lgnd->AddEntry(&distErrIpk, Form("Mean: %.2f #pm %.2f",
        distErrIpk.ProfileX()->GetMean(), distErrIpk.ProfileX()->GetMeanError()), "");
  lgnd->AddEntry(&distErrIpk, Form("RMS: %.2f #pm %.2f",
        distErrIpk.ProfileX()->GetRMS(), distErrIpk.ProfileX()->GetRMSError()), "");
  lgnd->AddEntry(&distErrIpk, "ErrCIpk", "");
  lgnd->AddEntry(&distErrIpk, Form("Mean: %.2f #pm %.2f",
        distErrIpk.ProfileY()->GetMean(), distErrIpk.ProfileY()->GetMeanError()), "");
  lgnd->AddEntry(&distErrIpk, Form("RMS: %.2f #pm %.2f",
        distErrIpk.ProfileY()->GetRMS(), distErrIpk.ProfileY()->GetRMSError()), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsZErrIpkVsCIpk.Print("results/resLeftRight_ZErrIpkVsCIpk.pdf");

  exit(0);
}


void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetLeftMargin(0.08);
  canvas.SetRightMargin(0.02);
}


void setTGraphStyle(TGraphErrors& graph) {
  graph.GetYaxis()->SetLabelSize(0.05);
  graph.GetYaxis()->SetTitleSize(0.05);
  graph.GetXaxis()->SetLabelSize(0.05);
  graph.GetXaxis()->SetTitleSize(0.05);
}
void setTGraphStyle(TGraph& graph) {
  graph.GetYaxis()->SetLabelSize(0.05);
  graph.GetYaxis()->SetTitleSize(0.05);
  graph.GetXaxis()->SetLabelSize(0.05);
  graph.GetXaxis()->SetTitleSize(0.05);
}
void setTGraphStyle(TH1D& graph) {
  graph.GetYaxis()->SetLabelSize(0.05);
  graph.GetYaxis()->SetTitleSize(0.05);
  graph.GetXaxis()->SetLabelSize(0.05);
  graph.GetXaxis()->SetTitleSize(0.05);
}
void setTGraphStyle(TH2D& graph) {
  graph.GetYaxis()->SetLabelSize(0.05);
  graph.GetYaxis()->SetTitleSize(0.05);
  graph.GetXaxis()->SetLabelSize(0.05);
  graph.GetXaxis()->SetTitleSize(0.05);
  graph.GetZaxis()->SetLabelSize(0.05);
  graph.GetZaxis()->SetTitleSize(0.05);
}
