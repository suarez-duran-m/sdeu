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

void comparisonChOffline()
{
  TString offlineName = "/home/msd/2021/offlineCourse/Offline2020/practice/uub/offlineUubDecSt863Pmt1.root";
  //TString offlineName = "/home/msd/2021/offlineCourse/Offline2020/practice/uub/Dec/offlineUubDecSt863Pmt1.root";
  TString cdasmeName = "uubAoPPMT1St863Mthdec.root";

  TPaveStats *ptstats;
  
  TFile *offlineF = TFile::Open(offlineName);
  TFile *cdasmeF = TFile::Open(cdasmeName);

  TTree *offlineInfo = (TTree*)offlineF->Get("charge");
  TTree *cdasmeInfo = (TTree*)cdasmeF->Get("ChargeData");

  int chOffEvId = 0;
  int chOffTime = 0;
  double chOffVem = 0.;
  double chOffChi2 = 0.;
  int chOffNdof = 0;
  double chOffLow = 0.;
  double chOffHigh = 0.;
  double chOffP0 = 0.;
  double chOffP1 = 0.;
  double chOffP2 = 0.;

  offlineInfo->SetBranchAddress("evtId", &chOffEvId);
  offlineInfo->SetBranchAddress("GpsTime", &chOffTime);
  offlineInfo->SetBranchAddress("chargeVal", &chOffVem);
  offlineInfo->SetBranchAddress("chargeChi2", &chOffChi2);
  offlineInfo->SetBranchAddress("chargeNdofl", &chOffNdof);
  offlineInfo->SetBranchAddress("chargeLow", &chOffLow);
  offlineInfo->SetBranchAddress("chargeHigh", &chOffHigh);
  offlineInfo->SetBranchAddress("chargeP0", &chOffP0);
  offlineInfo->SetBranchAddress("chargeP1", &chOffP1);
  offlineInfo->SetBranchAddress("chargeP2", &chOffP2);

  int NntrsOff = offlineInfo->GetEntries();
  double choffId[NntrsOff];
  double choffTime[NntrsOff];
  double choffVem[NntrsOff];
  double choffChi2[NntrsOff];
  double choffNdf[NntrsOff];
  double choffLow[NntrsOff];
  double choffHigh[NntrsOff];
  double choffch2Ndof[NntrsOff];
  
  int nChOks = 0;
  double aveChOff = 0.;
  double rmsChOff = 0.;
  double aveChi2NdosOf = 0.;
  double rmsChi2NdosOf = 0.;
  int nFailsOff = 0;
  vector < int > idFailsOff;
  vector < int > etryFailsOff;

  for ( int ntry=0; ntry<NntrsOff; ntry++ )
  {
    offlineInfo->GetEntry(ntry);
    choffId[ntry] = chOffEvId;
    choffTime[ntry] = chOffTime;
    choffVem[ntry] = chOffVem;
    choffChi2[ntry] = chOffChi2;
    choffNdf[ntry] = chOffNdof;
    choffLow[ntry] = chOffLow;
    choffHigh[ntry] = chOffHigh;
    choffch2Ndof[ntry] = chOffChi2/chOffNdof;
    if ( chOffVem > 1e3 )
    {
      aveChOff += chOffVem;
      aveChi2NdosOf += chOffChi2/chOffNdof;
      nChOks++;
    }
    else
    {
      nFailsOff++;
      idFailsOff.push_back( chOffEvId );
      etryFailsOff.push_back( ntry );
    }
  }
  aveChOff /= nChOks;
  aveChi2NdosOf /= nChOks;
  for ( int i=0; i<NntrsOff; i++ )
    if ( choffVem[i] > 1e3 )
    {
      rmsChOff += (choffVem[i] - aveChOff)*(choffVem[i] - aveChOff);
      rmsChi2NdosOf += (choffch2Ndof[i] - aveChi2NdosOf)*(choffch2Ndof[i] - aveChi2NdosOf);
    }
  rmsChOff = sqrt( rmsChOff/nChOks );
  rmsChi2NdosOf = sqrt( rmsChi2NdosOf/nChOks );

  int chCdasEvId = 0;
  int chCdasTime = 0;
  double chCdasVem = 0.;
  double chCdasVemDer = 0.;
  double chCdasChi2 = 0.;
  int chCdasNdof = 0;
  double chCdasLow = 0.;
  double chCdasHigh = 0.;
  double chCdasPar0 = 0.;
  double chCdasPar1 = 0.;
  double chCdasPar2 = 0.;

  cdasmeInfo->SetBranchAddress("eventId", &chCdasEvId);
  cdasmeInfo->SetBranchAddress("timeEvnt", &chCdasTime);
  cdasmeInfo->SetBranchAddress("chargeVal", &chCdasVem);
  cdasmeInfo->SetBranchAddress("chargeValDer", &chCdasVemDer);
  cdasmeInfo->SetBranchAddress("chi2", &chCdasChi2);
  cdasmeInfo->SetBranchAddress("ndf", &chCdasNdof);
  cdasmeInfo->SetBranchAddress("low", &chCdasLow);
  cdasmeInfo->SetBranchAddress("high", &chCdasHigh);
  cdasmeInfo->SetBranchAddress("chPar0", &chCdasPar0);
  cdasmeInfo->SetBranchAddress("chPar1", &chCdasPar1);
  cdasmeInfo->SetBranchAddress("chPar1", &chCdasPar2); 
 
  int NntrsCdas = cdasmeInfo->GetEntries();
  double chcdasId[NntrsCdas];
  double chcdasTime[NntrsCdas];
  double chcdasVem[NntrsCdas];
  double chcdasVemDer[NntrsCdas];
  double chcdasChi2[NntrsCdas];
  double chcdasNdf[NntrsCdas];
  double chcdasLow[NntrsCdas];
  double chcdasHigh[NntrsCdas];
  double chcdasch2Ndof[NntrsCdas];

  double aveChCdas = 0.;
  double aveChCdasDer = 0.;
  double rmsChCdas = 0.;
  double rmsChCdasDer = 0.;
  double aveChi2NdosCd = 0.;
  double rmsChi2NdosCd = 0.;
  int nFailsCdas = 0;
  nChOks = 0;
  for ( int ntry=0; ntry<NntrsCdas; ntry++ )
  {
    cdasmeInfo->GetEntry(ntry);
    chcdasId[ntry] = chCdasEvId;
    chcdasTime[ntry] = chCdasTime;
    chcdasVem[ntry] = chCdasVem;
    chcdasVemDer[ntry] = chCdasVemDer;
    chcdasChi2[ntry] = chCdasChi2;
    chcdasNdf[ntry] = chCdasNdof;
    chcdasLow[ntry] = chCdasLow;
    chcdasHigh[ntry] = chCdasHigh;
    chcdasch2Ndof[ntry] = chCdasChi2/chCdasNdof;
    if ( chCdasVem > 1e3 )
    {
      aveChCdas += chCdasVem;
      aveChCdasDer += chCdasVemDer;
      aveChi2NdosCd += chCdasChi2/chCdasNdof;
      nChOks++;
    }
    else
      nFailsCdas++;
  }
  aveChCdas /= nChOks;
  aveChCdasDer /= nChOks;
  aveChi2NdosCd /= nChOks;
  rmsChCdas = 0.;
  for ( int i=0; i<NntrsCdas; i++ )
    if ( chcdasVem[i] > 1e3 )
    {
      rmsChCdas += (chcdasVem[i] - aveChCdas)*(chcdasVem[i] - aveChCdas);
      rmsChCdasDer += (chcdasVemDer[i] - aveChCdasDer)*(chcdasVemDer[i] - aveChCdasDer);
      rmsChi2NdosCd += (chcdasch2Ndof[i] - aveChi2NdosCd)*(chcdasch2Ndof[i] - aveChi2NdosCd);
    }
  rmsChCdas = sqrt( rmsChCdas/nChOks );
  rmsChCdasDer = sqrt( rmsChCdasDer/nChOks );
  rmsChi2NdosCd = sqrt( rmsChi2NdosCd/nChOks );

  TLegend *leg;
  TLine *line;
  TString strChMean;
  TString strChRms;
  TString strFails;
  TString strNtrs;
  TString strPerFails;
  strChMean.Form("%.2f", aveChOff);
  strChRms.Form("%.2f", rmsChOff);
  strFails.Form("%d", nFailsOff);
  strNtrs.Form("%d", NntrsOff);
  strPerFails.Form("%.2f", (double)nFailsOff/NntrsOff);

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  TGraph *gr1 = new TGraph (NntrsOff, choffTime, choffVem);
  gr1->SetTitle("");
  gr1->GetYaxis()->SetTitle("VEM Charge [FADC]");
  gr1->GetYaxis()->SetRangeUser(200, 1400);
  gr1->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr1->GetXaxis()->SetTimeFormat("%m/%d");
  gr1->GetXaxis()->SetTimeOffset(315964782,"gmt");
  
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerSize(2);
  gr1->SetLineWidth(0);
  histoStyle(gr1);
  gr1->Draw("AP");

  TGraph *gr2 = new TGraph (NntrsCdas, chcdasTime, chcdasVem);
  gr2->SetMarkerStyle(32);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerSize(2);
  gr2->Draw("P same");

  leg = new TLegend(0.52,0.3,0.82,0.72);
  leg->SetHeader("        PMT1");
  leg->AddEntry(gr1, "Average Charge: "+strChMean,"p");
  leg->AddEntry(gr1, "RMS:  "+strChRms,"f");
  leg->AddEntry(gr1, "Fail Fits: "+strFails+"/"+strNtrs+" = "+strPerFails,"f");

  strChMean.Form("%.2f", aveChCdas);
  strChRms.Form("%.2f", rmsChCdas);
  strFails.Form("%d", nFailsCdas);
  strNtrs.Form("%d", NntrsCdas);
  strPerFails.Form("%.2f", (double)nFailsCdas/NntrsCdas);

  leg->SetHeader("        PMT1");
  leg->AddEntry(gr2, "Average Charge: "+strChMean,"p");
  leg->AddEntry(gr2, "RMS:  "+strChRms,"");
  leg->AddEntry(gr2, "Fail Fits: "+strFails+"/"+strNtrs+" = "+strPerFails,"");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("../plots/offlineVEMchSt863Pmt1.pdf");


  TCanvas *c21 = canvasStyle("c21");
  c21->cd();
  TGraph *gr21 = new TGraph (NntrsCdas, chcdasTime, chcdasVemDer);
  gr21->SetTitle(""); 
  gr21->GetYaxis()->SetTitle("VEMch-BXL [FADC]");
  gr21->GetYaxis()->SetRangeUser(200, 1400);
  gr21->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr21->GetXaxis()->SetTimeFormat("%m/%d");
  gr21->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr21->SetMarkerStyle(32);
  gr21->SetMarkerColor(kRed);
  gr21->SetMarkerSize(1.5);
  gr21->SetLineWidth(0);
  histoStyle(gr21);
  gr21->Draw("AP");

  strChMean.Form("%.2f", aveChCdasDer);
  strChRms.Form("%.2f", rmsChCdasDer);

  leg = new TLegend(0.12,0.31,0.42,0.55);
  leg->SetHeader("        PMT1");
  leg->AddEntry(gr2, "Average Charge from Derivative: "+strChMean,"f");
  leg->AddEntry(gr2, "RMS:  "+strChRms,"f");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();


  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  TGraph *gr3 = new TGraph (NntrsCdas, choffTime, choffNdf);
  gr3->SetTitle(""); 
  gr3->GetYaxis()->SetTitle("NDF-OffLine");
  gr3->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr3->GetXaxis()->SetTimeFormat("%m/%d");
  gr3->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr3->SetMarkerStyle(8);
  gr3->SetMarkerColor(kBlue);
  gr3->SetMarkerSize(1.5);
  gr3->SetLineWidth(0);
  histoStyle(gr3);
  gr3->Draw("AP");
  c3->Print("../plots/offlineNDFchSt863Pmt1.pdf");

  TCanvas *c4 = canvasStyle("c4");
  c4->cd();
  TGraph *gr4 = new TGraph (NntrsCdas, chcdasTime, chcdasNdf);
  gr4->SetTitle(""); 
  gr4->GetYaxis()->SetTitle("NDF-BXL");
  gr4->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr4->GetXaxis()->SetTimeFormat("%m/%d");
  gr4->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr4->SetMarkerStyle(32);
  gr4->SetMarkerColor(kRed);
  gr4->SetMarkerSize(1.5);
  gr4->SetLineWidth(0);
  histoStyle(gr4);
  gr4->Draw("AP");
  c4->Print("../plots/cdasNDFchSt863Pmt1.pdf");

  TCanvas *c5 = canvasStyle("c5");
  c5->cd();
  TGraph *gr5 = new TGraph (NntrsCdas, choffTime, choffChi2);
  gr5->SetTitle(""); 
  gr5->GetYaxis()->SetTitle("Chi2-OffLine");
  gr5->GetYaxis()->SetRangeUser(0, 220);
  gr5->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr5->GetXaxis()->SetTimeFormat("%m/%d");
  gr5->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr5->SetMarkerStyle(8);
  gr5->SetMarkerColor(kBlue);
  gr5->SetMarkerSize(1.5);
  gr5->SetLineWidth(0);
  histoStyle(gr5);
  gr5->Draw("AP");
  c5->Print("../plots/offlineChi2chSt863Pmt1.pdf");
                                                                          
  TCanvas *c6 = canvasStyle("c6");
  c6->cd();
  TGraph *gr6 = new TGraph (NntrsCdas, chcdasTime, chcdasChi2);
  gr6->SetTitle(""); 
  gr6->GetYaxis()->SetTitle("Chi2-BXL");
  gr6->GetYaxis()->SetRangeUser(0, 220);
  gr6->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr6->GetXaxis()->SetTimeFormat("%m/%d");
  gr6->GetXaxis()->SetTimeOffset(315964782,"gmt");
  gr6->SetMarkerStyle(32);
  gr6->SetMarkerColor(kRed);
  gr6->SetMarkerSize(1.5);
  gr6->SetLineWidth(0);
  histoStyle(gr6);
  gr6->Draw("AP");
  c6->Print("../plots/cdasChi2chSt863Pmt1.pdf");

  TCanvas *c7 = canvasStyle("c7");
  c7->cd();
  TGraph *gr7 = new TGraph (NntrsCdas, choffLow, choffHigh);
  gr7->SetTitle(""); 
  gr7->GetYaxis()->SetTitle("RangeHigh-OffLine [FADC]");
  gr7->GetYaxis()->SetRangeUser(350,2100);
  gr7->GetXaxis()->SetTitle("RangeLow-OffLine [FADC]");
  gr7->SetMarkerStyle(8);
  gr7->SetMarkerColor(kBlue);
  gr7->SetMarkerSize(1.5);
  gr7->SetLineWidth(0);
  histoStyle(gr7);
  gr7->Draw("AP");
  c7->Print("../plots/offlineLowHihgchSt863Pmt1.pdf");

  TCanvas *c8 = canvasStyle("c8");
  c8->cd();
  TGraph *gr8 = new TGraph (NntrsCdas, chcdasLow, chcdasHigh);
  gr8->SetTitle(""); 
  gr8->GetYaxis()->SetTitle("RangeHigh-BXL [FADC]");
  gr8->GetYaxis()->SetRangeUser(350,2100);
  gr8->GetXaxis()->SetTitle("RangeLow-BXL [FADC]");
  gr8->SetMarkerStyle(32);
  gr8->SetMarkerColor(kRed);
  gr8->SetMarkerSize(1.5);
  gr8->SetLineWidth(0);
  histoStyle(gr8);
  gr8->Draw("AP");
  c8->Print("../plots/cdasLowHighchSt863Pmt1.pdf");

  TString strAveChi2Ndf;
  TString strRmsChi2Ndf;
  strAveChi2Ndf.Form("%.2f", aveChi2NdosOf);
  strRmsChi2Ndf.Form("%.2f", rmsChi2NdosOf);

  TCanvas *c9 = canvasStyle("c9");
  c9->cd();
  TGraph *gr9 = new TGraph (NntrsOff, choffTime, choffch2Ndof);
  gr9->SetTitle("");
  gr9->GetYaxis()->SetTitle("Chi2/NDF");
  gr9->GetYaxis()->SetRangeUser(0, 2);
  gr9->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr9->GetXaxis()->SetTimeFormat("%m/%d");
  gr9->GetXaxis()->SetTimeOffset(315964782,"gmt");

  gr9->SetMarkerStyle(8);
  gr9->SetMarkerColor(kBlue);
  gr9->SetMarkerSize(1.5);
  gr9->SetLineWidth(0);
  histoStyle(gr9);
  gr9->Draw("AP");

  leg = new TLegend(0.12,0.2,0.42,0.45);
  leg->SetHeader("        PMT1");
  leg->AddEntry(gr9, "Average Chi2/NDF: "+strAveChi2Ndf,"f");
  leg->AddEntry(gr9, "RMS:  "+strRmsChi2Ndf,"f");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();

  c9->Print("../plots/offlineChi2NdfChSt863Pmt1.pdf");


  strAveChi2Ndf.Form("%.2f", aveChi2NdosCd);
  strRmsChi2Ndf.Form("%.2f", rmsChi2NdosCd);

  TCanvas *c10 = canvasStyle("c10");
  c10->cd();
  TGraph *gr10 = new TGraph (NntrsCdas, chcdasTime, chcdasch2Ndof);
  gr10->SetTitle("");
  gr10->GetYaxis()->SetTitle("Chi2/NDF");
  gr10->GetYaxis()->SetRangeUser(0, 2);
  gr10->GetXaxis()->SetTitle("Time since December 1st, 2020 [month/day]");
  gr10->GetXaxis()->SetTimeFormat("%m/%d");
  gr10->GetXaxis()->SetTimeOffset(315964782,"gmt");

  gr10->SetMarkerStyle(32);
  gr10->SetMarkerColor(kRed);
  gr10->SetMarkerSize(1.5);
  gr10->SetLineWidth(0);
  histoStyle(gr10);
  gr10->Draw("AP");

  leg = new TLegend(0.12,0.4,0.42,0.65);
  leg->SetHeader("        PMT1");
  leg->AddEntry(gr10, "Average Chi2/NDF: "+strAveChi2Ndf,"f");
  leg->AddEntry(gr10, "RMS:  "+strRmsChi2Ndf,"f");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();

  c10->Print("../plots/cdasChi2NdfChSt863Pmt1.pdf");


  TGraphErrors *graphFitted = new TGraphErrors();

  TCanvas *c11 = canvasStyle("c11");
  cdasmeInfo->SetBranchAddress("graph", &graphFitted);
  c11->cd();
  cdasmeInfo->GetEntry(0);

  TF1 *offFunct = new TF1("offFunct","-0.000164817*(x-1238.37)*(x-1238.37) + 89.2575",916,1556);
  TF1 *poly22 = new TF1("poly22","[0]*x*x + [1]*x + [2]",916,1556);
  TF1 *poly = graphFitted->GetFunction("poly2");

  graphFitted->SetTitle("");
  graphFitted->GetXaxis()->SetTitle("[FADC]");
  graphFitted->GetYaxis()->SetTitle("Counts [FADC]");
  graphFitted->GetXaxis()->SetRangeUser(0, 2500);
  graphFitted->GetYaxis()->SetRangeUser(0, 170);
  graphFitted->SetMarkerStyle(89);
  graphFitted->SetMarkerColor(kBlack);
  graphFitted->SetLineColor(kGray);
  histoStyle(graphFitted);
  graphFitted->Draw("ap");

  
  poly->SetLineColor(kRed);
  poly->SetLineWidth(0);
  poly->Draw("same");
  

  offFunct->SetLineColor(kBlue);
  offFunct->SetLineWidth(4);
  offFunct->Draw("same");

  TString vemFit;
  TString vemDer;
  TString vemOff;
  vemFit.Form("%.2f", chCdasVem);
  vemDer.Form("%.2f", chCdasVemDer);
  vemOff.Form("%.2f", 1238.37);

  leg = new TLegend(0.5,.66,0.8,0.95);
  //leg->SetHeader("        PMT1");
  //leg->AddEntry(offFunct, "VEM from OffLine: "+vemOff,"p");
  //leg->AddEntry(poly, "VEM from Poly2: "+vemFit,"p");
  //leg->AddEntry(poly, "VEM from BXL-method: "+vemDer,"p");
  //leg->AddEntry(poly, "(Green line)","");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();
/*
  line = new TLine(chCdasVemDer, 0, chCdasVemDer, 170);
  line->SetLineColor(kGreen+3);
  line->SetLineStyle(4);
  line->SetLineWidth(3);
  line->Draw();
  */
  c11->Print("../plots/offlineChCompaSt863Pmt1.pdf");

  vector < double > xRsdCdas;
  vector < double > yRsdCdas;
  vector < double > errRsdCdas;
  vector < double > xRsdOff;
  vector < double > yRsdOff;
  vector < double > errRsdOff;

  TGraphErrors *grphFitt2 = (TGraphErrors*)graphFitted->Clone();
  grphFitt2->Fit("poly22","R");
  cerr << poly22->GetChisquare() << " " << poly22->GetNDF() << " "
    << poly22->GetChisquare() / poly22->GetNDF() << " "
    << -poly22->GetParameter(1) / (2.*poly22->GetParameter(0)) << " "
    << -poly->GetParameter(1) / (2.*poly->GetParameter(0)) 
    << endl;

  getResiduals( graphFitted, offFunct, 916, 1556, xRsdOff, yRsdOff, errRsdOff );
  //getResiduals( grphFitt2, poly22, 916, 1556, xRsdCdas, yRsdCdas, errRsdCdas );
  getResiduals( graphFitted, poly, chCdasLow, chCdasHigh, xRsdCdas, yRsdCdas, errRsdCdas );

  TGraphErrors* rsdGrphOff = new TGraphErrors( xRsdOff.size(), &xRsdOff.front(),
      &yRsdOff.front(), 0, &errRsdOff.front() );

  TGraphErrors* rsdGrphCdas = new TGraphErrors( xRsdCdas.size(), &xRsdCdas.front(),
      &yRsdCdas.front(), 0, &errRsdCdas.front() );


  TCanvas *c12 = canvasStyle("c12");
  c12->cd();

  rsdGrphCdas->SetTitle("");
  rsdGrphCdas->GetXaxis()->SetTitle("[FADC]");
  rsdGrphCdas->GetYaxis()->SetTitle("Residuals");
  rsdGrphCdas->GetXaxis()->SetRangeUser(chCdasLow-2, chCdasHigh+2);
  rsdGrphCdas->GetYaxis()->SetRangeUser(-5, 2.5);
  rsdGrphCdas->SetLineColor(kRed);
  rsdGrphCdas->SetLineWidth(2);
  rsdGrphCdas->SetMarkerStyle(20);//23);
  rsdGrphCdas->SetMarkerColor(kRed);
  rsdGrphCdas->SetMarkerSize(2);//3);
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

  leg = new TLegend(0.12,0.16,0.75,0.48);
  leg->SetHeader("        PMT1");
  TString strchi2;
  TString strndf;
  TString strchi2Ndf;
  
  strchi2.Form("%.2f", 71.1869);
  strndf.Form("%d", 77);
  strchi2Ndf.Form("%.2f", 71.1869/77.);
  leg->AddEntry(rsdGrphOff,"OffLine, #chi^{2}/ndf = "+strchi2+"/"+strndf+" = "+strchi2Ndf,"ep");
  strchi2.Form("%.2f", chCdasChi2 );
  strndf.Form("%d", chCdasNdof );
  strchi2Ndf.Form("%.2f", chCdasChi2/chCdasNdof);
  leg->AddEntry(rsdGrphCdas,"Poly2, #chi^{2}/ndf = "+strchi2+"/"+strndf+" = "+strchi2Ndf,"ep"); 
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  line = new TLine(chCdasLow-2, 0, chCdasHigh+2, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();  
  
  line = new TLine(chCdasLow-2, 1, chCdasHigh+2, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(chCdasLow-2, -1, chCdasHigh+2, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  c12->Print("../plots/offlineResidChSt863Pmt1.pdf");

  int chosenEvt = 0;
  for ( int i=0; i<cdasmeInfo->GetEntries(); i++ )
  {
    cdasmeInfo->GetEntry(i);
    if ( chCdasEvId == idFailsOff[15] )
    {
      chosenEvt = i; 
      break;
    }
  }
  
  cdasmeInfo->GetEntry(chosenEvt);
  
  TCanvas *c13 = canvasStyle("c13");
  c13->cd();

  poly = graphFitted->GetFunction("poly2");
  //offFunct = new TF1("offFunct","-0.00126428*(x-138.443)*(x-138.443) + 318.509",chOffLow,chOffHigh);
  offFunct = new TF1("offFunct","-0.00131403*(x-321.677)*(x-321.677) + 140.319",190,450);
  offlineInfo->GetEntry( etryFailsOff[15] );
  cout << "Ch: " << chOffVem << endl;
  cout << "P0: " << chOffP0 << endl;
  cout << "P1: " << chOffP1 << endl;
  cout << "P2: " << chOffP2 << endl;
  cout << "Charge Low: " << chOffLow << endl;
  cout << "Charge High: " << chOffHigh << endl;

  graphFitted->SetTitle("");
  graphFitted->GetXaxis()->SetTitle("[FADC/8.33 ns]");
  graphFitted->GetYaxis()->SetTitle("Counts/FADC");
  graphFitted->GetXaxis()->SetRangeUser(0, 3500);
  graphFitted->GetYaxis()->SetRangeUser(0, 170);
  graphFitted->SetMarkerStyle(89);
  graphFitted->SetMarkerColor(kBlack);
  graphFitted->SetLineColor(kGray);
  histoStyle(graphFitted);
  graphFitted->Draw("ap");

  poly->SetLineColor(kRed);
  poly->SetLineWidth(0);
  poly->Draw("same");

  offFunct->SetLineColor(kBlue);
  offFunct->SetLineWidth(4);
  offFunct->Draw("same");

  //vemCdas.Form("%.2f", chCdasVem);
  vemOff.Form("%.2f", 281.677);
  TString strEvt;
  strEvt.Form("%d", chCdasEvId);

  leg = new TLegend(0.39,0.66,0.69,0.86);
  leg->AddEntry(offFunct, "Failed fit for OffLine, VEM-Charge: "+vemOff, "l");
  leg->AddEntry(offFunct, "(Event "+strEvt+")", "");
  //leg->AddEntry(poly, "Fit from BXL, VEM-Charge: "+vemCdas,"l");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  cerr << chCdasTime << endl;
  cerr << strEvt << endl;
  c13->Print("../plots/offlineFailedChargeSt863PMT1Evt"+strEvt+".pdf");
/*
  line = new TLine(chCdasVemDer, 0, chCdasVemDer, 170);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(chCdasVem, 0, chCdasVem, 170);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  xRsdCdas.clear();
  yRsdCdas.clear();
  errRsdCdas.clear();
  xRsdOff.clear();
  yRsdOff.clear();
  errRsdOff.clear();

  //getResiduals( graphFitted, offFunct, 98, 178, xRsdOff, yRsdOff, errRsdOff );
  getResiduals( graphFitted, poly, chCdasLow, chCdasHigh, xRsdCdas, yRsdCdas, errRsdCdas );

  //TGraphErrors* rsdGrphOff = new TGraphErrors( xRsdOff.size(), &xRsdOff.front(),
      //&yRsdOff.front(), 0, &errRsdOff.front() );

  TGraphErrors *rsdGrphCdas2 = new TGraphErrors( xRsdCdas.size(), &xRsdCdas.front(),
      &yRsdCdas.front(), 0, &errRsdCdas.front() );

  TCanvas *c14 = canvasStyle("c14");

  rsdGrphCdas2->SetTitle("");
  rsdGrphCdas2->GetXaxis()->SetTitle("[FADC]");
  rsdGrphCdas2->GetYaxis()->SetTitle("Residuals");
  rsdGrphCdas2->GetXaxis()->SetRangeUser(chCdasLow-2, chCdasHigh+2);
  rsdGrphCdas2->GetYaxis()->SetRangeUser(-6, 5);
  rsdGrphCdas2->SetLineColor(kRed);
  rsdGrphCdas2->SetLineWidth(2);
  rsdGrphCdas2->SetMarkerStyle(20);
  rsdGrphCdas2->SetMarkerColor(kRed);
  rsdGrphCdas2->SetMarkerSize(2);
  histoStyle(rsdGrphCdas2);
  rsdGrphCdas2->Draw("AP same");
  */
/*
  rsdGrphOff->SetTitle("");
  rsdGrphOff->SetLineColor(kBlue);
  rsdGrphOff->SetMarkerColor(kBlue);
  rsdGrphOff->SetLineWidth(2);
  rsdGrphOff->SetMarkerStyle(20);
  rsdGrphOff->SetMarkerSize(2);
  histoStyle(rsdGrphOff);
  rsdGrphOff->Draw("P same");
*/
/*
  line = new TLine(chCdasLow-2, 0, chCdasHigh+2, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();  
  
  line = new TLine(chCdasLow-2, 1, chCdasHigh+2, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(chCdasLow-2, -1, chCdasHigh+2, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  strchi2.Form("%.2f", chCdasChi2);
  strndf.Form("%d", chCdasNdof);
  strchi2Ndf.Form("%.2f", chCdasChi2/chCdasNdof);

  leg = new TLegend(0.50,0.19,0.76,0.33);
  leg->AddEntry(rsdGrphCdas2, "BXL, #chi^{2}/ndf = "+strchi2+"/"+strndf+" = "+strchi2Ndf,"f");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  strEvt.Form("%d", chCdasEvId);
  c14->Print("../plots/offlineFailedChResidSt863PMT1Evt"+strEvt+".pdf");
  */
}
