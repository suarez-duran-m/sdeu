void setCanvasStyle(TCanvas& canvas);
void setCanvasStyle2D(TCanvas& canvas);
void setTGraphStyle(TGraphErrors& graph);
void setTGraphStyle(TGraph& graph);
void setTGraphStyle(TH1D& graph);
void setTGraphStyle(TH2D& graph);

void annaBinsLRQpk() {
  ifstream dataFile("binsLefltRightQpk_FebJul.dat");

  const int totalHistos = 102290;
  int nbins = 1000;
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
  double tmpQpk = 0.;
  double tmpErrQpk = 0.;
  double tmpChi2Ndf = 0.;
  double tmpNdof = 0.; 
  double tmpPval = 0.;

  double checkPval = 0.;

  while ( dataFile >> tmpNfadc >> tmpQpk >> tmpErrQpk >> tmpChi2Ndf >> tmpNdof ) {
    checkPval = TMath::Prob(tmpChi2Ndf*tmpNdof, tmpNdof);    
    if ( checkPval < 1e-10 )
      continue;
    // Error coming as ErrQpk/Qpk
    tmpErrQpk *= 100.;

    arrErrIQpk[int(tmpNfadc)] += tmpErrQpk;
    arrErrIQpkDev[int(tmpNfadc)] += tmpErrQpk*tmpErrQpk;

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

  ifstream dataFileLeftRight("distLefltRightQpk.dat");
  TH1D distCntsLeftRight("distCntsLeftRight", "", 1.25e2, -1., 1.);
  distCntsLeftRight.SetTitle("; Counts (Left-Right)/Total [au]; counts [au]");
  double tmpCntLeftRight = 0;
  
  while ( dataFileLeftRight >> tmpCntLeftRight )
    distCntsLeftRight.Fill( tmpCntLeftRight );
  dataFileLeftRight.close();

  ifstream dataFileAfterFit("QpkVsCoinciQpk_FebJul2.dat");
  TH2D distQpk("distQpk", "", 50, 1e3, 3e3, 50, 1e3, 3e3);
  distQpk.SetTitle("; Qpk [FADC]; CQpk [FADC]; Counts [au]");
  TH2D distErrQpk("distErrQpk", "", 2e3, 0, 100, 2e3, 0, 100);
  distErrQpk.SetTitle("; ErrQ^{pk}/Q^{pk} [%]; ErrCQ^{pk}/CQ^{pk} [%]; Counts [au]");
  TH2D distChi2("distChi2", "", 90, 0, 5, 90, 0, 5);
  distChi2.SetTitle("; #chi^{2}/Ndof; #chi^{2}_{Coinci}/Ndof; Counts [au]");
  TH2D distEntropy("distEntropy", "", 200, 0, 10, 2e3, 0, 100);
  distEntropy.SetTitle("; Entropy [au]; ErrCQ^{pk}/CQ^{pk} [%]; Counts [au]");

  double tmpCQpk = 0.;
  double tmpErrCQpk = 0.;
  double tmpChi2NdfQpk = 0.;
  double tmpChi2NdfCQpk = 0.;
  double tmpEntropy = 0.;
  double cqpkOverqpk = 0.;
  int n = 0;

  while ( dataFileAfterFit >> tmpQpk >> tmpErrQpk >> tmpCQpk >> tmpErrCQpk >> tmpChi2NdfQpk >> tmpChi2NdfCQpk >> tmpEntropy) {
    if ( tmpQpk < 1 || tmpCQpk < 1 )
      continue;
    distQpk.Fill( tmpQpk, tmpCQpk );
    distErrQpk.Fill( 100.*tmpErrQpk, 100.*tmpErrCQpk );
    distChi2.Fill( tmpChi2NdfQpk, tmpChi2NdfCQpk );
    distEntropy.Fill( tmpEntropy, 100.*tmpErrCQpk );
    if ( tmpQpk > 1e3 && tmpQpk < 1.9e3 ) {
      cqpkOverqpk += tmpCQpk / tmpQpk;
      n++;
    }
  }
  dataFileAfterFit.close();
  cout << "MSD " << cqpkOverqpk/n << " " << n << endl;

  int nMinus1 = 0;
  for ( int i=0; i<nbins; i++ ) {
    if ( nCntsInFadc[i] == 0 )
      continue;
   
    arrErrIQpk[i] = arrErrIQpk[i]/nCntsInFadc[i];
    arrErrIQpkDev[i] 
      = sqrt(arrErrIQpkDev[i]/nCntsInFadc[i] - arrErrIQpk[i]*arrErrIQpk[i]);
    arrErrIQpkDev[i] /= sqrt(nCntsInFadc[i]);

    arrChi2Ndf[i] = arrChi2Ndf[i]/nCntsInFadc[i];
    arrChi2NdfDev[i] 
      = sqrt(arrChi2NdfDev[i]/nCntsInFadc[i] - arrChi2Ndf[i]*arrChi2Ndf[i]);
    arrChi2NdfDev[i] /= sqrt(nCntsInFadc[i]);

    arrPval[i] = arrPval[i]/nCntsInFadc[i];
    arrPvalDev[i] = sqrt(arrPvalDev[i]/nCntsInFadc[i] - arrPval[i]*arrPval[i]);
    arrPvalDev[i] /= sqrt(nCntsInFadc[i]);
    arrPval[i] *= -1.;

    nCntsInFadc[i] = 100.*(nCntsInFadc[i]/totalHistos);
  }

  TGraphErrors errIQpkgrph(xAxisNdof.size(), &xAxisNdof.front(), &arrErrIQpk.front(), 
      0, &arrErrIQpkDev.front());
  errIQpkgrph.SetTitle("; Ndof [au]; #LT ErrQpk/Qpk #GT [%]");
  
  TCanvas cvnsErrIQpk("cvnsErrIQpk","",1.6e3, 9e2);
  setCanvasStyle(cvnsErrIQpk);
  cvnsErrIQpk.cd();
  cvnsErrIQpk.SetLogy();

  errIQpkgrph.GetXaxis()->SetRangeUser(1, 250);
  errIQpkgrph.GetYaxis()->SetRangeUser(2e-1, 3e3);
  errIQpkgrph.SetMarkerStyle(20);
  errIQpkgrph.SetMarkerSize(2);
  errIQpkgrph.SetMarkerColor(kBlue);
  errIQpkgrph.SetLineColor(kBlue);
  setTGraphStyle(errIQpkgrph);
  errIQpkgrph.GetYaxis()->SetTitleOffset(0.7);
  errIQpkgrph.Draw("ap");

  cvnsErrIQpk.Print("results/resLeftRight_ErrQpk.pdf");


  TGraphErrors redChi2grph(xAxisNdof.size(), &xAxisNdof.front(), &arrChi2Ndf.front(), 
      0, &arrChi2NdfDev.front());
  redChi2grph.SetTitle("; Ndof [au]; #chi^{2}/Ndof [FADC]");
  
  TCanvas cvnsRedChi2("cvnsRedChi2","",1.6e3, 9e2);
  setCanvasStyle(cvnsRedChi2);
  cvnsRedChi2.cd();

  redChi2grph.GetXaxis()->SetRangeUser(1, 250);  
  redChi2grph.GetYaxis()->SetRangeUser(0.9, 3);
  redChi2grph.SetMarkerStyle(20);
  redChi2grph.SetMarkerSize(2);
  redChi2grph.SetMarkerColor(kBlue);
  redChi2grph.SetLineColor(kBlue);
  setTGraphStyle(redChi2grph);
  redChi2grph.GetYaxis()->SetTitleOffset(0.8);
  redChi2grph.Draw("ap");

  cvnsRedChi2.Print("results/resLeftRight_redChi2Qpk.pdf");


  TGraphErrors logPvalgrph(xAxisNdof.size(), &xAxisNdof.front(), &arrPval.front(), 
      0, &arrPvalDev.front());
  logPvalgrph.SetTitle("; Ndof [au]; (-1)*Log(Pval) [au]");
  
  TCanvas cvnsLogPval("cvnsLogPval","",1.6e3, 9e2);
  setCanvasStyle(cvnsLogPval);
  cvnsLogPval.cd();
  cvnsLogPval.SetLogy();

  logPvalgrph.GetXaxis()->SetRangeUser(1, 250);
  logPvalgrph.GetYaxis()->SetRangeUser(0.3, 10.);
  logPvalgrph.SetMarkerStyle(20);
  logPvalgrph.SetMarkerSize(2);
  logPvalgrph.SetMarkerColor(kBlue);
  logPvalgrph.SetLineColor(kBlue);
  setTGraphStyle(logPvalgrph);
  logPvalgrph.GetYaxis()->SetTitleOffset(0.7);
  logPvalgrph.Draw("ap");

  cvnsLogPval.Print("results/resLeftRight_logPvalQpk.pdf");

  TGraph succesHistos(xAxisNdof.size(), &xAxisNdof.front(), &nCntsInFadc.front());
  succesHistos.SetTitle("; Ndof [au]; Fit_{Success}/TotalHistos [%]");

  TCanvas cvnsSuccHistos("cvnsSuccHistos","",1.6e3, 9e2);
  setCanvasStyle(cvnsSuccHistos);
  cvnsSuccHistos.cd();

  succesHistos.GetXaxis()->SetRangeUser(1, 250);
  //succesHistos.GetYaxis()->SetRangeUser(78, 100);
  succesHistos.SetMarkerStyle(20);
  succesHistos.SetMarkerSize(2);
  succesHistos.SetMarkerColor(kBlue);
  succesHistos.SetLineColor(kBlue);
  setTGraphStyle(succesHistos);
  succesHistos.GetYaxis()->SetTitleOffset(0.7);
  succesHistos.Draw("ap");

  cvnsSuccHistos.Print("results/resLeftRight_succesHistsoQpk.pdf");


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

  TLegend *lgnd = new TLegend(0.7, 0.75, 0.9, 0.94);
  lgnd->AddEntry(&distCntsLeftRight,Form("Mean: %.3f #pm %.3f",
      distCntsLeftRight.GetMean(), distCntsLeftRight.GetMeanError()), "l");
  lgnd->AddEntry(distCntsLeftRight.GetFunction("gaus"), Form("#mu = %.3f; #sigma = %.3f",
        distCntsLeftRight.GetFunction("gaus")->GetParameter(1),
        distCntsLeftRight.GetFunction("gaus")->GetParameter(2)),"l");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDistCntLeft.Print("results/resLeftRight_DistLeftRightQpk.pdf");


  TCanvas cvnsQpkVsCQpk("cvnsQpkVsCQpk","",1.6e3, 9e2);
  setCanvasStyle2D(cvnsQpkVsCQpk);
  cvnsQpkVsCQpk.cd();

  distQpk.Fit("pol1", "Q", "R", 1.e3, 1.9e3);
  distQpk.SetStats(kFALSE);
  distQpk.GetYaxis()->SetRangeUser(1.1e3,2.7e3);
  distQpk.GetXaxis()->SetRangeUser(1.1e3,2.7e3);
  distQpk.SetMarkerStyle(20);
  distQpk.SetMarkerSize(2);
  distQpk.GetFunction("pol1")->SetLineWidth(4);
  setTGraphStyle(distQpk);
  distQpk.Draw("colz");

  lgnd = new TLegend(0.15, 0.75, 0.23, 0.96);
  lgnd->AddEntry(&distQpk, "Qpk:", "");
  lgnd->AddEntry(&distQpk, Form("Etrs/TotCHist: %.2f",
        distQpk.GetEntries()/totalHistos), ""); 
  lgnd->AddEntry(&distQpk, Form("Mean: %.2f #pm %.2f",
        distQpk.ProfileX()->GetMean(), distQpk.ProfileX()->GetMeanError()), "");
  lgnd->AddEntry(&distQpk, Form("RMS: %.2f #pm %.2f",
        distQpk.ProfileX()->GetRMS(), distQpk.ProfileX()->GetRMSError()), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();
  
  lgnd = new TLegend(0.4, 0.8, 0.48, 0.96);
  lgnd->AddEntry(&distQpk, "CQpk:", "");
  lgnd->AddEntry(&distQpk, Form("Mean: %.2f #pm %.2f",
        distQpk.ProfileY()->GetMean(), distQpk.ProfileY()->GetMeanError()), "");
  lgnd->AddEntry(&distQpk, Form("RMS: %.2f #pm %.2f",
        distQpk.ProfileY()->GetRMS(), distQpk.ProfileY()->GetRMSError()), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  lgnd = new TLegend(0.52, 0.25, 0.72, 0.42);
  lgnd->AddEntry(distQpk.GetFunction("pol1"), "Fit:", "l");
  lgnd->AddEntry(distQpk.GetFunction("pol1"), Form("Slope = %.2f #pm %.2f",
        distQpk.GetFunction("pol1")->GetParameter(1), 
        distQpk.GetFunction("pol1")->GetParError(1)), "");
  lgnd->AddEntry(distQpk.GetFunction("pol1"), Form("b = %.2f #pm %.2f",
        distQpk.GetFunction("pol1")->GetParameter(0),
        distQpk.GetFunction("pol1")->GetParError(0)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.05);
  lgnd->Draw();

  cvnsQpkVsCQpk.Print("results/resLeftRight_QpkVsCQpk.pdf");


  TCanvas cvnsErrQpkVsCQpk("cvnsErrQpkVsCQpk","",1.6e3, 9e2);
  setCanvasStyle2D(cvnsErrQpkVsCQpk);
  cvnsErrQpkVsCQpk.SetLeftMargin(0.09);
  cvnsErrQpkVsCQpk.SetRightMargin(0.16);
  cvnsErrQpkVsCQpk.cd();

  distErrQpk.SetStats(kFALSE);
  distErrQpk.GetYaxis()->SetRangeUser(0.25,10.5);
  distErrQpk.GetXaxis()->SetRangeUser(0.25,2);
  distErrQpk.SetMarkerStyle(20);
  distErrQpk.SetMarkerColor(kBlue);
  distErrQpk.SetMarkerSize(2);
  setTGraphStyle(distErrQpk);
  distErrQpk.GetYaxis()->SetTitleOffset(.55);
  distErrQpk.GetZaxis()->SetTitleOffset(.75);
  distErrQpk.Draw("colz");

  lgnd = new TLegend(0.55, 0.55, 0.75, 0.95);
  lgnd->AddEntry(&distErrQpk, "ErrQpk", "");
  lgnd->AddEntry(&distErrQpk, Form("Total entries: %.f",
        distErrQpk.GetEntries()), "");
  lgnd->AddEntry(&distErrQpk, Form("Mean: %.2f #pm %.2f",
        distErrQpk.ProfileX()->GetMean(), distErrQpk.ProfileX()->GetMeanError()), "");
  lgnd->AddEntry(&distErrQpk, Form("RMS: %.2f #pm %.2f",
        distErrQpk.ProfileX()->GetRMS(), distErrQpk.ProfileX()->GetRMSError()), "");
  lgnd->AddEntry(&distErrQpk, "ErrCQpk", "");
  lgnd->AddEntry(&distErrQpk, Form("Mean: %.2f #pm %.2f",
        distErrQpk.ProfileY()->GetMean(), distErrQpk.ProfileY()->GetMeanError()), "");
  lgnd->AddEntry(&distErrQpk, Form("RMS: %.2f #pm %.2f",
        distErrQpk.ProfileY()->GetRMS(), distErrQpk.ProfileY()->GetRMSError()), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsErrQpkVsCQpk.Print("results/resLeftRight_ErrQpkVsCQpk.pdf");


  TCanvas cvnsChi2("cvnsChi2","",1.6e3, 9e2);
  setCanvasStyle2D(cvnsChi2);
  cvnsChi2.SetLeftMargin(0.11);
  cvnsChi2.SetRightMargin(0.17);
  cvnsChi2.cd();

  distChi2.SetStats(kFALSE);
  distChi2.GetYaxis()->SetRangeUser(0.3, 2.6);
  distChi2.GetXaxis()->SetRangeUser(0.3, 2.6);
  distChi2.SetMarkerStyle(20);
  distChi2.SetMarkerSize(2);
  setTGraphStyle(distChi2);
  distChi2.GetYaxis()->SetTitleOffset(.7);
  distChi2.GetZaxis()->SetTitleOffset(0.9);
  distChi2.Draw("colz");

  lgnd = new TLegend(0.55, 0.8, 0.65, 0.96);
  lgnd->AddEntry(&distChi2, "Regular histograms:", "");
  lgnd->AddEntry(&distChi2, Form("Mean: %.3f #pm %.3f",
        distChi2.ProfileX()->GetMean(), distChi2.ProfileX()->GetMeanError()), "");
  lgnd->AddEntry(&distChi2, Form("RMS: %.3f #pm %.3f",
        distChi2.ProfileX()->GetRMS(), distChi2.ProfileX()->GetRMSError()), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();
  
  lgnd = new TLegend(0.55, 0.59, 0.65, 0.75);
  lgnd->AddEntry(&distChi2, "Coincidence histogram:", "");
  lgnd->AddEntry(&distChi2, Form("Mean: %.3f #pm %.3f",
        distChi2.ProfileY()->GetMean(), distChi2.ProfileY()->GetMeanError()), "");
  lgnd->AddEntry(&distChi2, Form("RMS: %.3f #pm %.3f",
        distChi2.ProfileY()->GetRMS(), distChi2.ProfileY()->GetRMSError()), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsChi2.Print("results/resLeftRight_Chi2.pdf");

  TCanvas cvnsEntropy("cvnsEntropy","",1.6e3, 9e2);
  setCanvasStyle2D(cvnsEntropy);
  cvnsEntropy.cd();

  distEntropy.SetStats(kFALSE);
  setTGraphStyle(distEntropy);
  distEntropy.GetXaxis()->SetRangeUser(5.1, 6.3);
  distEntropy.GetYaxis()->SetRangeUser(0., 10);
  distEntropy.Draw("colz");

  cvnsEntropy.Print("results/entropyCQ.pdf");

  exit(0);
}


void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetBottomMargin(0.17);
  canvas.SetLeftMargin(0.11);
  canvas.SetRightMargin(0.07);
}

void setCanvasStyle2D(TCanvas& canvas) {
  canvas.SetTopMargin(0.02);
  canvas.SetBottomMargin(0.14);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.18);
}


void setTGraphStyle(TGraphErrors& graph) {
  graph.GetYaxis()->SetTitleSize(0.07);
  graph.GetYaxis()->SetLabelSize(0.06);
  graph.GetYaxis()->SetTitleOffset(0.1);
  graph.GetXaxis()->SetTitleSize(0.07);
  graph.GetXaxis()->SetLabelSize(0.06);
  graph.GetXaxis()->SetTitleOffset(1.);
}
void setTGraphStyle(TGraph& graph) {
  graph.GetYaxis()->SetTitleSize(0.07);
  graph.GetYaxis()->SetLabelSize(0.06);
  graph.GetYaxis()->SetTitleOffset(0.1);
  graph.GetXaxis()->SetTitleSize(0.07);
  graph.GetXaxis()->SetLabelSize(0.06);
  graph.GetXaxis()->SetTitleOffset(1.);
}
void setTGraphStyle(TH1D& graph) {
  graph.GetYaxis()->SetTitleSize(0.07);
  graph.GetYaxis()->SetLabelSize(0.06);
  graph.GetYaxis()->SetTitleOffset(0.1);
  graph.GetXaxis()->SetTitleSize(0.07);
  graph.GetXaxis()->SetLabelSize(0.06);
  graph.GetXaxis()->SetTitleOffset(1.);
}
void setTGraphStyle(TH2D& graph) {
  graph.GetYaxis()->SetTitleSize(0.07);
  graph.GetYaxis()->SetLabelSize(0.06);
  graph.GetYaxis()->SetTitleOffset(.88);
  graph.GetXaxis()->SetTitleSize(0.07);
  graph.GetXaxis()->SetLabelSize(0.06);
  graph.GetXaxis()->SetTitleOffset(.9);
  graph.GetZaxis()->SetTitleSize(0.07);
  graph.GetZaxis()->SetLabelSize(0.06);
  graph.GetZaxis()->SetTitleOffset(0.92);
}
