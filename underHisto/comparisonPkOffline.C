TCanvas *canvasStyle(TString name)
{
  TCanvas *canvas = new TCanvas(name, name, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.11); 
  canvas->SetRightMargin(0.03);
  canvas->SetTopMargin(0.02); 
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  return canvas;
}

void histoStyle(TH1F *hist)
{
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


void histoStyle(TGraphErrors *hist)
{
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


void histoStyle(TGraph *hist)
{
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


void getResiduals( TGraphErrors *grphErr, TF1 *func,
    double rangMin, double rangMax,
    vector < double > &x, vector < double > &y, vector < double > &err )
{
  Double_t *xpnts = grphErr->GetX();
  Double_t *ypnts = grphErr->GetY();

  int nbins = grphErr->GetXaxis()->GetNbins();
  double tmp = 0.;
  for ( int kbin=1; kbin<nbins; kbin++ )
    if ( xpnts[kbin] >= rangMin && xpnts[kbin] <= rangMax ) 
    {
      x.push_back( xpnts[kbin] );
      tmp = func->Eval( xpnts[kbin] ) - ypnts[kbin];
      y.push_back( tmp / sqrt( ypnts[kbin] ) );
      err.push_back( 
          sqrt( pow(sqrt( ypnts[kbin] ),2)
            + pow(sqrt( sqrt(func->Eval( xpnts[kbin] ) ) ),2)
            ) / sqrt( ypnts[kbin] )
          );
    }
}

// =================================
// *** *** *** MAIN CODE *** *** ***

void comparisonPkOffline()
{
  TString offlineName = "/home/msd/2021/offlineCourse/Offline2020/practice/uub/offlineUubDecSt863Pmt1.root";
  TString cdasmeName = "uubAoPPMT1St863Mthdec.root";

  TPaveStats *ptstats;
  
  TFile *offlineF = TFile::Open(offlineName);
  TFile *cdasmeF = TFile::Open(cdasmeName);

  TTree *offlineInfo = (TTree*)offlineF->Get("peak");
  TTree *cdasmeInfo = (TTree*)cdasmeF->Get("PeakData");

  int pkOffEvId = 0;
  int pkOffTime = 0;
  double pkOffVem = 0.;
  double pkOffChi2 = 0.;
  int pkOffNdof = 0;
  double pkOffLow = 0.;
  double pkOffHigh = 0.;
  double pkOffPar0 = 0.;
  double pkOffPar1 = 0.;
  double pkOffPar2 = 0.;

  offlineInfo->SetBranchAddress("evtId", &pkOffEvId);
  offlineInfo->SetBranchAddress("GpsTime", &pkOffTime);
  offlineInfo->SetBranchAddress("peakVal", &pkOffVem);
  offlineInfo->SetBranchAddress("peakChi2", &pkOffChi2);
  offlineInfo->SetBranchAddress("peakNdofl", &pkOffNdof);
  offlineInfo->SetBranchAddress("peakLow", &pkOffLow);
  offlineInfo->SetBranchAddress("peakHigh", &pkOffHigh);
  offlineInfo->SetBranchAddress("peakP0", &pkOffPar0);
  offlineInfo->SetBranchAddress("peakP1", &pkOffPar1);
  offlineInfo->SetBranchAddress("peakP2", &pkOffPar2);

  int NntrsOff = offlineInfo->GetEntries();
  double pkoffId[NntrsOff];
  double pkoffTime[NntrsOff];
  double pkoffVem[NntrsOff];
  double pkoffChi2[NntrsOff];
  double pkoffNdf[NntrsOff];
  double pkoffLow[NntrsOff];
  double pkoffHigh[NntrsOff];
  double pkoffch2Ndof[NntrsOff];
  
  int nPkOks = 0;
  double avePkOff = 0.;
  double rmsPkOff = 0.;
  double aveChi2NdosOf = 0.;
  double rmsChi2NdosOf = 0.;
  int nFailsOff = 0;
  vector < int > idFailsOff;

  for ( int ntry=0; ntry<NntrsOff; ntry++ )
  {
    offlineInfo->GetEntry(ntry);
    pkoffId[ntry] = pkOffEvId;
    pkoffTime[ntry] = pkOffTime;
    pkoffVem[ntry] = pkOffVem;
    pkoffChi2[ntry] = pkOffChi2;
    pkoffNdf[ntry] = pkOffNdof;
    pkoffLow[ntry] = pkOffLow;
    pkoffHigh[ntry] = pkOffHigh;
    if ( pkOffVem > 0 )
    {
      avePkOff += pkOffVem;
      aveChi2NdosOf += pkOffChi2/pkOffNdof;
      pkoffch2Ndof[ntry] = pkOffChi2/pkOffNdof;
      nPkOks++;
      /*
      cerr << ntry << " " 
        << pkOffLow << " "
        << pkOffHigh << " "
        << pkOffVem << " "
        << pkOffPar0 << " " 
        << pkOffPar1 << " " 
        << pkOffPar2 
        << endl;
        */
    }
    else
    {
      nFailsOff++;
      idFailsOff.push_back( pkOffEvId );
      /*
      cerr << "FailOff: " << 
        ntry << " " 
        << pkOffEvId
        << endl;
        */
    }
  }
  avePkOff /= nPkOks;
  aveChi2NdosOf /= nPkOks;
  for ( int i=0; i<NntrsOff; i++ )
    if ( pkoffVem[i] > 0 )
    {
      rmsPkOff += (pkoffVem[i] - avePkOff)*(pkoffVem[i] - avePkOff);
      rmsChi2NdosOf += (pkoffch2Ndof[i] - aveChi2NdosOf)*(pkoffch2Ndof[i] - aveChi2NdosOf);
    }
  rmsPkOff = sqrt( rmsPkOff/nPkOks );
  rmsChi2NdosOf = sqrt( rmsChi2NdosOf/nPkOks );

  int pkCdasEvId = 0;
  int pkCdasTime = 0;
  double pkCdasVem = 0.;
  double pkCdasVemDer = 0.;
  double pkCdasChi2 = 0.;
  int pkCdasNdof = 0;
  double pkCdasLow = 0.;
  double pkCdasHigh = 0.;
  double pkCdasPar0 = 0.;
  double pkCdasPar1 = 0.;
  double pkCdasPar2 = 0.;

  cdasmeInfo->SetBranchAddress("eventId", &pkCdasEvId);
  cdasmeInfo->SetBranchAddress("timeEvnt", &pkCdasTime);
  cdasmeInfo->SetBranchAddress("peakVal", &pkCdasVem);
  cdasmeInfo->SetBranchAddress("peakValDer", &pkCdasVemDer);
  cdasmeInfo->SetBranchAddress("chi2", &pkCdasChi2);
  cdasmeInfo->SetBranchAddress("ndf", &pkCdasNdof);
  cdasmeInfo->SetBranchAddress("low", &pkCdasLow);
  cdasmeInfo->SetBranchAddress("high", &pkCdasHigh);
  cdasmeInfo->SetBranchAddress("pkPar0", &pkCdasPar0);
  cdasmeInfo->SetBranchAddress("pkPar1", &pkCdasPar1);
  cdasmeInfo->SetBranchAddress("pkPar2", &pkCdasPar2); 
 
  int NntrsCdas = cdasmeInfo->GetEntries();
  double pkcdasId[NntrsCdas];
  double pkcdasTime[NntrsCdas];
  double pkcdasVem[NntrsCdas];
  double pkcdasVemDer[NntrsCdas];
  double pkcdasChi2[NntrsCdas];
  double pkcdasNdf[NntrsCdas];
  double pkcdasLow[NntrsCdas];
  double pkcdasHigh[NntrsCdas];
  double pkcdasP0[NntrsCdas];
  double pkcdasP1[NntrsCdas];
  double pkcdasP2[NntrsCdas];
  double pkcdasch2Ndof[NntrsCdas];

  double avePkCdas = 0.;
  double avePkCdasDer = 0.;
  double rmsPkCdas = 0.;
  double rmsPkCdasDer = 0.;
  double aveChi2NdosCd = 0.;
  double rmsChi2NdosCd = 0.;
  int nFailsCdas = 0;
  nPkOks = 0;
  for ( int ntry=0; ntry<NntrsCdas; ntry++ )
  {
    cdasmeInfo->GetEntry(ntry);
    pkcdasId[ntry] = pkCdasEvId;
    pkcdasTime[ntry] = pkCdasTime;
    pkcdasVem[ntry] = pkCdasVem;
    pkcdasVemDer[ntry] = pkCdasVemDer;
    pkcdasChi2[ntry] = pkCdasChi2;
    pkcdasNdf[ntry] = pkCdasNdof;
    pkcdasLow[ntry] = pkCdasLow;
    pkcdasHigh[ntry] = pkCdasHigh;
    pkcdasP0[ntry] = pkCdasPar0;
    pkcdasP1[ntry] = pkCdasPar1;
    pkcdasP2[ntry] = pkCdasPar2;
    pkcdasch2Ndof[ntry] = pkCdasChi2/pkCdasNdof;
    if ( pkCdasVem > 0 )
    {
      avePkCdas += pkCdasVem;
      avePkCdasDer += pkCdasVemDer;
      aveChi2NdosCd += pkCdasChi2/pkCdasNdof;
      nPkOks++;
    }
    else
    {
      nFailsCdas++;
      //cerr << "Failed: " << pkCdasEvId << " " << ntry << " " << pkCdasTime << endl;
    }
    //if ( pkCdasVem > 0 && pkCdasVem < 130.)
      //cerr << pkCdasEvId << " " << ntry << " " << pkCdasTime << endl;
  }
  avePkCdas /= nPkOks;
  avePkCdasDer /= nPkOks;
  aveChi2NdosCd /= nPkOks;
  rmsPkCdas = 0.;
  for ( int i=0; i<NntrsCdas; i++ )
    if ( pkcdasVem[i] > 0 )
    {
      rmsPkCdas += (pkcdasVem[i] - avePkCdas)*(pkcdasVem[i] - avePkCdas);
      rmsPkCdasDer += (pkcdasVemDer[i] - avePkCdasDer)*(pkcdasVemDer[i] - avePkCdasDer);
      rmsChi2NdosCd += (pkcdasch2Ndof[i] - aveChi2NdosCd)*(pkcdasch2Ndof[i] - aveChi2NdosCd);
    }
  rmsPkCdas = sqrt( rmsPkCdas/nPkOks );
  rmsPkCdasDer = sqrt( rmsPkCdasDer/nPkOks );
  rmsChi2NdosCd = sqrt( rmsChi2NdosCd/nPkOks );

  TLegend *leg;
  TString strPkMean;
  TString strPkRms;
  TString strFails;
  TString strNtrs;
  TString strPerFails;
  strPkMean.Form("%.2f", avePkOff);
  strPkRms.Form("%.2f", rmsPkOff);
  strFails.Form("%d", nFailsOff);
  strNtrs.Form("%d", NntrsOff);
  strPerFails.Form("%.2f", (double)nFailsOff/NntrsOff);

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  TGraph *gr1 = new TGraph (NntrsOff, pkoffTime, pkoffVem);
  gr1->SetTitle("");
  gr1->GetYaxis()->SetTitle("VEM-Peak [FADC/8.33 ns]");
  gr1->GetYaxis()->SetRangeUser(0, 180);
  gr1->GetXaxis()->SetTitle("Time since December 1st, 2020 [day/month]");
  gr1->GetXaxis()->SetTimeFormat("%m/%d");
  gr1->GetXaxis()->SetTimeOffset(315964782,"gmt"); 
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerSize(2);
  gr1->SetLineWidth(0);
  histoStyle(gr1);
  gr1->Draw("AP");

  TGraph *gr2 = new TGraph (NntrsCdas, pkcdasTime, pkcdasVem);
  gr2->SetMarkerStyle(32);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerSize(2);
  gr2->Draw("P same");

  leg = new TLegend(0.52,0.22,0.82,0.62);
  leg->AddEntry(gr1, "Average Peak-OffLine: "+strPkMean,"p");
  leg->AddEntry(gr1, "RMS:  "+strPkRms,"");
  leg->AddEntry(gr1, "Fail Fits: "+strFails+"/"+strNtrs+" = "+strPerFails,"");

  strPkMean.Form("%.2f", avePkCdas);
  strPkRms.Form("%.2f", rmsPkCdas);
  strFails.Form("%d", nFailsCdas);
  strNtrs.Form("%d", NntrsCdas);
  strPerFails.Form("%.2f", (double)nFailsCdas/NntrsCdas);

  leg->AddEntry(gr2, "Average Peak-BXL: "+strPkMean,"p");
  leg->AddEntry(gr2, "RMS:  "+strPkRms,"");
  leg->AddEntry(gr2, "Fail Fits: "+strFails+"/"+strNtrs+" = "+strPerFails,"");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("../plots/offlineVEMpkSt863Pmt1.pdf");

  TCanvas *c21 = canvasStyle("c21");
  c21->cd();
  TGraph *gr21 = new TGraph (NntrsCdas, pkcdasTime, pkcdasVemDer);
  gr21->SetTitle(""); 
  gr21->GetYaxis()->SetTitle("VEMpk-BXL [FADC/8.33 ns]");
  gr21->GetYaxis()->SetRangeUser(0, 180);
  gr21->GetXaxis()->SetTitle("Time since December 1st, 2020 [day/month]");
  gr21->GetXaxis()->SetTimeFormat("%m/%d");
  gr21->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr21->SetMarkerStyle(32);
  gr21->SetMarkerColor(kRed);
  gr21->SetMarkerSize(1.5);
  gr21->SetLineWidth(0);
  histoStyle(gr21);
  gr21->Draw("AP");

  strPkMean.Form("%.2f", avePkCdasDer);
  strPkRms.Form("%.2f", rmsPkCdasDer);

  leg = new TLegend(0.12,0.31,0.42,0.5);
  leg->AddEntry(gr2, "Average Peak from Derivative: "+strPkMean,"f");
  leg->AddEntry(gr2, "RMS:  "+strPkRms,"f");
  leg->SetTextSize(0.065);
  leg->SetBorderSize(0);
  leg->Draw();


  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  TGraph *gr3 = new TGraph (NntrsCdas, pkoffTime, pkoffNdf);
  gr3->SetTitle(""); 
  gr3->GetYaxis()->SetTitle("NDF-OffLine");
  gr3->GetYaxis()->SetRangeUser(0, 33);
  gr3->GetXaxis()->SetTitle("Time since December 1st, 2020 [day/month]");
  gr3->GetXaxis()->SetTimeFormat("%m/%d");
  gr3->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr3->SetMarkerStyle(8);
  gr3->SetMarkerColor(kBlue);
  gr3->SetMarkerSize(1.5);
  gr3->SetLineWidth(0);
  histoStyle(gr3);
  gr3->Draw("AP");
  c3->Print("../plots/offlineNDFpkSt863Pmt1.pdf");

  TCanvas *c4 = canvasStyle("c4");
  c4->cd();
  TGraph *gr4 = new TGraph (NntrsCdas, pkcdasTime, pkcdasNdf);
  gr4->SetTitle(""); 
  gr4->GetYaxis()->SetTitle("NDF-BXL");
  gr4->GetYaxis()->SetRangeUser(0, 33);
  gr4->GetXaxis()->SetTitle("Time since December 1st, 2020 [day/month]");
  gr4->GetXaxis()->SetTimeFormat("%m/%d");
  gr4->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr4->SetMarkerStyle(32);
  gr4->SetMarkerColor(kRed);
  gr4->SetMarkerSize(1.5);
  gr4->SetLineWidth(0);
  histoStyle(gr4);
  gr4->Draw("AP");
  c4->Print("../plots/cdasNDFpkSt863Pmt1.pdf");


  TCanvas *c5 = canvasStyle("c5");
  c5->cd();
  TGraph *gr5 = new TGraph (NntrsCdas, pkoffTime, pkoffChi2);
  gr5->SetTitle(""); 
  gr5->GetYaxis()->SetTitle("Chi2-OffLine");
  gr5->GetYaxis()->SetRangeUser(0, 62);
  gr5->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr5->GetXaxis()->SetTimeFormat("%m/%d");
  gr5->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr5->SetMarkerStyle(8);
  gr5->SetMarkerColor(kBlue);
  gr5->SetMarkerSize(1.5);
  gr5->SetLineWidth(0);
  histoStyle(gr5);
  gr5->Draw("AP");
  c5->Print("../plots/offlineChi2pkSt863Pmt1.pdf");


  TCanvas *c6 = canvasStyle("c6");
  c6->cd();
  TGraph *gr6 = new TGraph (NntrsCdas, pkcdasTime, pkcdasChi2);
  gr6->SetTitle(""); 
  gr6->GetYaxis()->SetTitle("Chi2-BXL");
  gr6->GetYaxis()->SetRangeUser(0, 62);
  gr6->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr6->GetXaxis()->SetTimeFormat("%m/%d");
  gr6->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr6->SetMarkerStyle(32);
  gr6->SetMarkerColor(kRed);
  gr6->SetMarkerSize(1.5);
  gr6->SetLineWidth(0);
  histoStyle(gr6);
  gr6->Draw("AP");
  c6->Print("../plots/cdasChi2pkSt863Pmt1.pdf");


  TCanvas *c7 = canvasStyle("c7");
  c7->cd();
  TGraph *gr7 = new TGraph (NntrsCdas, pkoffLow, pkoffHigh);
  gr7->SetTitle(""); 
  gr7->GetYaxis()->SetTitle("RangeHigh-OffLine [FADC]");
  gr7->GetYaxis()->SetRangeUser(0, 270);
  gr7->GetXaxis()->SetTitle("RangeLow-OffLine [FADC]");
  gr7->SetMarkerStyle(8);
  gr7->SetMarkerColor(kBlue);
  gr7->SetMarkerSize(1.5);
  gr7->SetLineWidth(0);
  histoStyle(gr7);
  gr7->Draw("AP");
  c7->Print("../plots/offlineLowHihgpkSt863Pmt1.pdf");


  TCanvas *c8 = canvasStyle("c8");
  c8->cd();
  TGraph *gr8 = new TGraph (NntrsCdas, pkcdasLow, pkcdasHigh);
  gr8->SetTitle(""); 
  gr8->GetYaxis()->SetTitle("RangeHigh-BXL [FADC]");
  gr8->GetYaxis()->SetRangeUser(0, 270);
  gr8->GetXaxis()->SetTitle("RangeLow-BXL [FADC]");
  gr8->SetMarkerStyle(32);
  gr8->SetMarkerColor(kRed);
  gr8->SetMarkerSize(1.5);
  gr8->SetLineWidth(0);
  histoStyle(gr8);
  gr8->Draw("AP");
  c8->Print("../plots/cdasLowHighpkSt863Pmt1.pdf");


  TString strAveChi2Ndf;
  TString strRmsChi2Ndf;
  strAveChi2Ndf.Form("%.2f", aveChi2NdosOf);
  strRmsChi2Ndf.Form("%.2f", rmsChi2NdosOf);

  TCanvas *c9 = canvasStyle("c9");
  c9->cd();
  TGraph *gr9 = new TGraph (NntrsOff, pkoffTime, pkoffch2Ndof);
  gr9->SetTitle("");
  gr9->GetYaxis()->SetTitle("Chi2/NDF");
  gr9->GetYaxis()->SetRangeUser(0, 4.2);
  gr9->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr9->GetXaxis()->SetTimeFormat("%m/%d");
  gr9->GetXaxis()->SetRangeUser(pkcdasTime[0]-1e5, pkcdasTime[NntrsCdas-1]+1e5);
  gr9->GetXaxis()->SetTimeOffset(315964782,"gmt");

  gr9->SetMarkerStyle(8);
  gr9->SetMarkerColor(kBlue);
  gr9->SetMarkerSize(1.5);
  gr9->SetLineWidth(0);
  histoStyle(gr9);
  gr9->Draw("AP");

  TLine *line = new TLine(pkcdasTime[0]-1e5,3.5,pkcdasTime[NntrsCdas-1]+1e5,3.5);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  leg = new TLegend(0.12,0.2,0.42,0.4);
  leg->AddEntry(gr9, "Average Chi2/NDF: "+strAveChi2Ndf,"f");
  leg->AddEntry(gr9, "RMS:  "+strRmsChi2Ndf,"f");
  leg->SetTextSize(0.065);
  leg->SetBorderSize(0);
  leg->Draw();

  c9->Print("../plots/offlineChi2NdfPkSt863Pmt1.pdf");


  strAveChi2Ndf.Form("%.2f", aveChi2NdosCd);
  strRmsChi2Ndf.Form("%.2f", rmsChi2NdosCd);

  TCanvas *c10 = canvasStyle("c10");
  c10->cd();
  TGraph *gr10 = new TGraph (NntrsCdas, pkcdasTime, pkcdasch2Ndof);
  gr10->SetTitle("");
  gr10->GetYaxis()->SetTitle("Chi2/NDF");
  gr10->GetYaxis()->SetRangeUser(0, 4.2);
  gr10->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr10->GetXaxis()->SetTimeFormat("%m/%d");
  gr10->GetXaxis()->SetRangeUser(pkcdasTime[0]-1e5, pkcdasTime[NntrsCdas-1]+1e5);
  gr10->GetXaxis()->SetTimeOffset(315964782,"gmt");

  gr10->SetMarkerStyle(32);
  gr10->SetMarkerColor(kRed);
  gr10->SetMarkerSize(1.5);
  gr10->SetLineWidth(0);
  histoStyle(gr10);
  gr10->Draw("AP");

  line = new TLine(pkcdasTime[0]-1e5,3.5,pkcdasTime[NntrsCdas-1]+1e5,3.5);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  leg = new TLegend(0.12,0.6,0.42,0.8);
  leg->AddEntry(gr10, "Average Chi2/NDF: "+strAveChi2Ndf,"f");
  leg->AddEntry(gr10, "RMS:  "+strRmsChi2Ndf,"f");
  leg->SetTextSize(0.065);
  leg->SetBorderSize(0);
  leg->Draw();

  c10->Print("../plots/cdasChi2NdfPkSt863Pmt1.pdf");


  // =======================
  // *** Doing Residuals ***

  TGraphErrors *graphFitted = new TGraphErrors();

  TCanvas *c11 = canvasStyle("c11");
  cdasmeInfo->SetBranchAddress("graph", &graphFitted);
  c11->cd();
  cdasmeInfo->GetEntry(0); // For event 61219267

  TF1 *offFunct = new TF1("offFunct","0.0242158*(x-151.596)*(x-151.596) + 577.861",98,178);
  TF1 *poly = graphFitted->GetFunction("fitFcn");

  graphFitted->SetTitle("");
  graphFitted->GetXaxis()->SetTitle("[FADC/8.33 ns]");
  graphFitted->GetYaxis()->SetTitle("Counts/FADC");
  graphFitted->GetXaxis()->SetRangeUser(0, 400);
  graphFitted->GetYaxis()->SetRangeUser(1e2, 1e3);
  graphFitted->SetMarkerStyle(89);
  graphFitted->SetMarkerColor(kBlack);
  graphFitted->SetLineColor(kGray);
  histoStyle(graphFitted);
  graphFitted->Draw("ap");

  poly->SetLineColor(kRed);
  poly->SetLineWidth(4);
  poly->Draw("same");

  offFunct->SetLineColor(kBlue);
  offFunct->SetLineWidth(4);
  offFunct->Draw("same");

  TString vemFit;
  TString vemDer;
  TString vemOff;
  vemFit.Form("%.2f", pkCdasVem);
  vemDer.Form("%.2f", pkCdasVemDer);
  vemOff.Form("%.2f", 151.596);

  leg = new TLegend(0.5,.66,0.8,0.86);
  leg->AddEntry(offFunct, "VEM from OffLine: "+vemOff,"p");
  leg->AddEntry(poly, "VEM from Poly2: "+vemFit,"p");
  leg->AddEntry(poly, "VEM from BXL-method: "+vemDer,"p");
  leg->AddEntry(poly, "(Green line)","");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  line = new TLine(pkCdasVemDer, 1e2, pkCdasVemDer, 1e3);
  line->SetLineColor(kGreen+3);
  line->SetLineStyle(4);
  line->SetLineWidth(3);
  line->Draw();
  c11->Print("../plots/offlinePkCompaSt863Pmt1.pdf");


  vector < double > xRsdCdas;
  vector < double > yRsdCdas;
  vector < double > errRsdCdas;
  vector < double > xRsdOff;
  vector < double > yRsdOff;
  vector < double > errRsdOff;

  getResiduals( graphFitted, offFunct, 98, 178, xRsdOff, yRsdOff, errRsdOff );
  getResiduals( graphFitted, poly, 98, pkCdasHigh, xRsdCdas, yRsdCdas, errRsdCdas );

  TGraphErrors* rsdGrphOff = new TGraphErrors( xRsdOff.size(), &xRsdOff.front(),
      &yRsdOff.front(), 0, &errRsdOff.front() );

  TGraphErrors* rsdGrphCdas = new TGraphErrors( xRsdCdas.size(), &xRsdCdas.front(),
      &yRsdCdas.front(), 0, &errRsdCdas.front() );

  TCanvas *c12 = canvasStyle("c12");

  rsdGrphCdas->SetTitle("");
  rsdGrphCdas->GetXaxis()->SetTitle("[FADC]");
  rsdGrphCdas->GetYaxis()->SetTitle("Residuals");
  rsdGrphCdas->GetXaxis()->SetRangeUser(96, pkCdasHigh+2);
  rsdGrphCdas->GetYaxis()->SetRangeUser(-6, 5);
  rsdGrphCdas->SetLineColor(kRed);
  rsdGrphCdas->SetLineWidth(2);
  rsdGrphCdas->SetMarkerStyle(20);
  rsdGrphCdas->SetMarkerColor(kRed);
  rsdGrphCdas->SetMarkerSize(2);
  histoStyle(rsdGrphCdas);
  rsdGrphCdas->Draw("AP same");

  rsdGrphOff->SetTitle("");
  rsdGrphOff->SetLineColor(kBlue);
  rsdGrphOff->SetMarkerColor(kBlue);
  rsdGrphOff->SetLineWidth(2);
  rsdGrphOff->SetMarkerStyle(20);
  rsdGrphOff->SetMarkerSize(2);
  histoStyle(rsdGrphOff);
  rsdGrphOff->Draw("P same");

  line = new TLine(pkCdasLow, -6, pkCdasLow, 5);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->Draw();

  line = new TLine(96, 0, pkCdasHigh+2, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();  
  
  line = new TLine(96, 1, pkCdasHigh+2, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(96, -1, pkCdasHigh+2, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  TString strchi2;
  TString strndf;
  TString strchi2Ndf;
  strchi2.Form("%.2f", 43.7414);
  strndf.Form("%d", 17);
  strchi2Ndf.Form("%.2f", 43.7414/17.);

  leg = new TLegend(0.50,0.19,0.76,0.33);
  leg->AddEntry(rsdGrphOff, "OffLine, #chi^{2}/ndf = "+strchi2+"/"+strndf+" = "+strchi2Ndf,"ep");
  strchi2.Form("%.2f", pkCdasChi2);
  strndf.Form("%d", pkCdasNdof);
  strchi2Ndf.Form("%.2f", pkCdasChi2/pkCdasNdof);
  leg->AddEntry(rsdGrphCdas, "Poly2, #chi^{2}/ndf = "+strchi2+"/"+strndf+" = "+strchi2Ndf,"ep");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  c12->Print("../plots/offlineResidPkSt863Pmt1.pdf");

  int chosenEvt = 0;
  for ( int i=0; i<cdasmeInfo->GetEntries(); i++ )
  {
    cdasmeInfo->GetEntry(i);
    if ( pkCdasEvId == idFailsOff[15] )
    {
      chosenEvt = i; 
      break;
    }
  }

  cdasmeInfo->GetEntry(chosenEvt);

  TCanvas *c13 = canvasStyle("c13");
  c13->cd();

  poly = graphFitted->GetFunction("fitFcn");

  graphFitted->SetTitle("");
  graphFitted->GetXaxis()->SetTitle("[FADC/8.33 ns]");
  graphFitted->GetYaxis()->SetTitle("Counts/FADC");
  graphFitted->GetXaxis()->SetRangeUser(0, 400);
  graphFitted->GetYaxis()->SetRangeUser(1e2, 1e3);
  graphFitted->SetMarkerStyle(89);
  graphFitted->SetMarkerColor(kBlack);
  graphFitted->SetLineColor(kGray);
  histoStyle(graphFitted);
  graphFitted->Draw("ap");

  poly->SetLineColor(kRed);
  poly->SetLineWidth(4);
  poly->Draw("same");

  vemFit.Form("%.2f", pkCdasVem);

  line = new TLine(pkCdasVemDer,1e2,pkCdasVemDer,1e3);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(pkCdasVem,1e2,pkCdasVem,1e3);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  leg = new TLegend(0.39,0.66,0.69,0.86);
  leg->AddEntry(graphFitted, "Failed fit for OffLine", "f");
  leg->AddEntry(poly, "Fit from BXL, VEM-Peak: "+vemDer,"f");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  TString strEvt;
  strEvt.Form("%d", pkCdasEvId);
  cerr << pkCdasTime << endl;
  c13->Print("../plots/offlineFailedPeakSt863PMT1Evt"+strEvt+".pdf");
}
