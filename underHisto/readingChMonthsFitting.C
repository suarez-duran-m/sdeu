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


vector < double > fillingCh( TString bname, TString st, int pmt, TString whichInfo)
{
  TString monthUub[] = {"dec", "jan", "feb", "mar", "abr", "may", "jun", "jul"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname = bname + pmtId+st+"Mth";

  TFile *f;
  TTree *chargeInfo;
  double tmpVals = 0.;
  double tmpDer = 0.;
  vector < double > returnVals;
  int evtId = 0;
  TGraphErrors *gr = new TGraphErrors();
  TString ngr;
  TString vemDer;
  TString vemFit;
  TString pmtStr;

  for ( int month=0; month<8; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    chargeInfo = (TTree*)f->Get("ChargeData");

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals);
    if ( whichInfo=="chargeVal" )
      chargeInfo->SetBranchAddress("chargeValDer", &tmpDer);
    chargeInfo->SetBranchAddress("eventId", &evtId);
    chargeInfo->SetBranchAddress("graph", &gr);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      returnVals.push_back( tmpVals );
      //if ( pmt==2 && whichInfo=="chargeValDer" && tmpVals < 80 && tmpVals > 60 )
        //cerr << etry << endl;
     
      //if ( pmt==1 && whichInfo=="chargeVal" && tmpVals < 1000 && month == 5 )
      //if ( pmt==2 && whichInfo=="chargeVal" && tmpVals < 1200 && month == 4 )
      if ( pmt==3 && whichInfo=="chargeVal" && tmpVals < 1000 && month == 3 )
      {
        TCanvas *c0 = canvasStyle("c0");
        c0->ResetDrawn();
        c0->cd();
        ngr.Form("%d",etry);
        vemDer.Form("%.2f", tmpDer);
        vemFit.Form("%.2f", tmpVals);
        pmtStr.Form("%d", pmt);
        gr->SetTitle("");
        gr->GetXaxis()->SetTitle("[FADC]");
        gr->GetXaxis()->SetRangeUser(0, 4e3);
        gr->GetYaxis()->SetRangeUser(0, 165);
        gr->GetYaxis()->SetTitle("Counts [au]");
        histoStyle(gr);
        gr->Draw();

        TF1 *poly;
        poly = gr->GetFunction("poly2");
        poly->SetLineColor(kRed);
        poly->SetLineWidth(2);
        poly->Draw("same");
        cout << poly->GetParameter(0) << " " << poly->GetParameter(1) << " " << poly->GetParameter(2) << endl;

        TLegend *leg = new TLegend(0.5,0.5,0.76,0.8);
        leg->AddEntry(gr,"St 863 PMT "+pmtStr,"");
        leg->AddEntry(gr,"VEM from Fit: "+vemFit,"");
        leg->AddEntry(gr,"VEM from Der.: "+vemDer,"");
        leg->SetLineWidth(0);
        leg->SetTextSize(0.06);
        leg->Draw();

        TLine *line = new TLine(tmpDer, 0, tmpDer, 165);
        line->SetLineColor(kGreen+3);
        line->SetLineStyle(4);
        line->SetLineWidth(3);
        line->Draw();
        //c0->Print("kk.C");
        c0->Print("../plots/sampleChHistoDerVem"+ngr+"Pmt"+pmt+".pdf");
        //break;
      }
     
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}

vector < int > fillingCh( TString bname, TString st, int pmt, bool whichInfo)
{
  TString monthUub[] = {"dec", "jan", "feb", "mar", "abr", "may", "jun", "jul"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname = bname + pmtId+st+"Mth";

  TFile *f;
  TTree *peakInfo;
  int tmpVals = 0.;
  vector < int > returnVals;

  TString getinfo;
  if ( whichInfo )
    getinfo = "timeEvnt";
  else
   getinfo = "eventId";

  for ( int month=0; month<8; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    peakInfo = (TTree*)f->Get("PeakData");

    tmpVals = 0.;
    peakInfo->SetBranchAddress(getinfo, &tmpVals);

    for( int etry=0; etry<peakInfo->GetEntries(); etry++)
    {
      peakInfo->GetEntry(etry);
      returnVals.push_back( tmpVals );
    }
  }
  peakInfo->Delete();
  f->Delete();
  return returnVals;
}


double getmean( vector<double> arr )
{
  double mean = 0.;
  int nb = arr.size();
  int ngoodb = 0;
    for (int i=0; i<nb; i++)
      if ( arr[i] > 100 )
      {
        mean += arr[i];
        ngoodb++;
      }
  return mean/ngoodb;
}

double getrms( vector<double> arr, double meanarr )
{
  double rms = 0.;
  int nb = arr.size();
  int ngoodb = 0;
  for (int i=0; i<nb; i++)
    if ( arr[i] > 100 )
    {
      rms += (arr[i] - meanarr)*(arr[i] - meanarr);
      ngoodb++;
    }
  return sqrt(rms/ngoodb);
}


// =================================
// *** *** *** MAIN CODE *** *** ***

void readingChMonthsFitting(int st)
{
  TString statId;
  statId.Form("St%d", st);
  TString basename = "uubAoPPMT";

  TPaveStats *ptstats;
  TLegend *leg;

  vector < int > timePmt1;
  vector < int > timePmt2;
  vector < int > timePmt3;
  vector < double > chargePmt1;
  vector < double > chargePmt2;
  vector < double > chargePmt3;
  vector < double > chargeDerPmt1;
  vector < double > chargeDerPmt2;
  vector < double > chargeDerPmt3;

  bool gettime = true;
  timePmt1 = fillingCh( basename, statId, 1, gettime);
  timePmt2 = fillingCh( basename, statId, 2, gettime);
  timePmt3 = fillingCh( basename, statId, 3, gettime);
  TString whinfo = "chargeVal";
  chargePmt1 = fillingCh( basename, statId, 1, whinfo);
  chargePmt2 = fillingCh( basename, statId, 2, whinfo);
  chargePmt3 = fillingCh( basename, statId, 3, whinfo);
  whinfo = "chargeValDer";
  chargeDerPmt1 = fillingCh( basename, statId, 1, whinfo);
  chargeDerPmt2 = fillingCh( basename, statId, 2, whinfo);
  chargeDerPmt3 = fillingCh( basename, statId, 3, whinfo);

  int nPoints = timePmt1.size();
  double xtimePmt1[ nPoints ];
  double yChPmt1[ nPoints ];
  double yChDerPmt1[ nPoints ];
  double yPkMeanPmt1 = 0.;
  double yPkRmsPmt1 = 0.;
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt1[i] = timePmt1[i];
    yChPmt1[i] = chargePmt1[i];
    yChDerPmt1[i] = chargeDerPmt1[i];
  }
  double meanPmt1 = getmean( chargePmt1 );
  double rmsPmt1 = getrms( chargePmt1, meanPmt1 );
  double meanDerPmt1 = getmean( chargeDerPmt1 );
  double rmsDerPmt1 = getrms( chargeDerPmt1, meanDerPmt1 );
  TGraph *chPmt1 = new TGraph(nPoints,xtimePmt1,yChPmt1);
  TGraph *chDerPmt1 = new TGraph(nPoints,xtimePmt1,yChDerPmt1);

  nPoints = timePmt2.size();
  double xtimePmt2[ nPoints ];
  double yChPmt2[ nPoints ];
  double yChDerPmt2[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt2[i] = timePmt2[i];
    yChPmt2[i] = chargePmt2[i];
    yChDerPmt2[i] = chargeDerPmt2[i];
  }
  double meanPmt2 = getmean( chargePmt2 );
  double rmsPmt2 = getrms( chargePmt2, meanPmt2 );
  double meanDerPmt2 = getmean( chargeDerPmt2 );
  double rmsDerPmt2 = getrms( chargeDerPmt2, meanDerPmt2 );
  TGraph *chPmt2 = new TGraph(nPoints,xtimePmt2,yChPmt2);
  TGraph *chDerPmt2 = new TGraph(nPoints,xtimePmt2,yChDerPmt2);

  nPoints = timePmt3.size();
  double xtimePmt3[ nPoints ];
  double yChPmt3[ nPoints ];
  double yChDerPmt3[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt3[i] = timePmt3[i];
    yChPmt3[i] = chargePmt3[i];
    yChDerPmt3[i] = chargeDerPmt3[i];
  }
  double meanPmt3 = getmean( chargePmt3 );
  double rmsPmt3 = getrms( chargePmt3, meanPmt3 );
  double meanDerPmt3 = getmean( chargeDerPmt3 );
  double rmsDerPmt3 = getrms( chargeDerPmt3, meanDerPmt3 );
  TGraph *chPmt3 = new TGraph(nPoints,xtimePmt3,yChPmt3);
  TGraph *chDerPmt3 = new TGraph(nPoints,xtimePmt3,yChDerPmt3);

  TString avePkStr;
  TString rmsPkStr;

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt1->SetTitle("");
  chPmt1->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt1->GetXaxis()->SetTitle("Time since Dec. 2020 [month/day]");
  chPmt1->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chPmt1->GetYaxis()->SetTitle("Charge [FADC]");
  chPmt1->GetYaxis()->SetRangeUser(0, 1800);
  chPmt1->SetMarkerStyle(25);
  chPmt1->SetMarkerSize(2);
  chPmt1->SetMarkerColor(kAzure+10);
  chPmt1->SetLineColor(kAzure+10);
  chPmt1->SetLineWidth(2);
  histoStyle(chPmt1);
  chPmt1->Draw("AP same");

  chDerPmt1->SetMarkerStyle(32);
  chDerPmt1->SetMarkerSize(2);
  chDerPmt1->SetMarkerColor(kOrange+10);
  chDerPmt1->SetLineColor(kOrange+10);
  chDerPmt1->SetLineWidth(2);
  chDerPmt1->Draw("P same");

  avePkStr.Form("%.2f", meanPmt1);
  rmsPkStr.Form("%.2f", rmsPmt1);
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT1");
  leg->AddEntry(chPmt1, "Average Charge-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  avePkStr.Form("%.2f", meanDerPmt1);
  rmsPkStr.Form("%.2f", rmsDerPmt1);
  leg->AddEntry(chDerPmt1, "Average Charge-Der.: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("../plots/uubChargeFromDerSt"+statId+"pmt1.pdf");


  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt2->SetTitle("");
  chPmt2->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt2->GetXaxis()->SetTitle("Time since Dec. 2020 [month/day]");
  chPmt2->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chPmt2->GetYaxis()->SetTitle("Charge [FADC]");
  chPmt2->GetYaxis()->SetRangeUser(0, 1800);
  chPmt2->SetMarkerStyle(25);
  chPmt2->SetMarkerSize(2);
  chPmt2->SetMarkerColor(kAzure+10);
  chPmt2->SetLineColor(kAzure+10);
  chPmt2->SetLineWidth(2);
  histoStyle(chPmt2);
  chPmt2->Draw("AP same");

  chDerPmt2->SetMarkerStyle(32);
  chDerPmt2->SetMarkerSize(2);
  chDerPmt2->SetMarkerColor(kOrange+10);
  chDerPmt2->SetLineColor(kOrange+10);
  chDerPmt2->SetLineWidth(2);
  chDerPmt2->Draw("P same");

  avePkStr.Form("%.2f", meanPmt2);
  rmsPkStr.Form("%.2f", rmsPmt2);
  leg = new TLegend(0.15,0.21,0.52,0.4);
  leg->SetHeader("PMT2");
  leg->AddEntry(chPmt2, "Average Charge-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  avePkStr.Form("%.2f", meanDerPmt2);
  rmsPkStr.Form("%.2f", rmsDerPmt2);
  leg->AddEntry(chDerPmt2, "Average Charge-Der.: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c2->Print("../plots/uubChargeFromDerSt"+statId+"pmt2.pdf");


  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt3->SetTitle("");
  chPmt3->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt3->GetXaxis()->SetTitle("Time since Dec. 2020 [month/day]");
  chPmt3->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chPmt3->GetYaxis()->SetTitle("Charge [FADC]");
  chPmt3->GetYaxis()->SetRangeUser(0, 1800);
  chPmt3->SetMarkerStyle(25);
  chPmt3->SetMarkerSize(2);
  chPmt3->SetMarkerColor(kAzure+10);
  chPmt3->SetLineColor(kAzure+10);
  chPmt3->SetLineWidth(2);
  histoStyle(chPmt3);
  chPmt3->Draw("AP same");

  chDerPmt3->SetMarkerStyle(32);
  chDerPmt3->SetMarkerSize(2);
  chDerPmt3->SetMarkerColor(kOrange+10);
  chDerPmt3->SetLineColor(kOrange+10);
  chDerPmt3->SetLineWidth(2);
  chDerPmt3->Draw("P same");

  avePkStr.Form("%.2f", meanPmt3);
  rmsPkStr.Form("%.2f", rmsPmt3);
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT3");
  leg->AddEntry(chPmt3, "Average Charge-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  avePkStr.Form("%.2f", meanDerPmt3);
  rmsPkStr.Form("%.2f", rmsDerPmt3);
  leg->AddEntry(chDerPmt3, "Average Charge-Der.: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c3->Print("../plots/uubChargeFromDerSt"+statId+"pmt3.pdf");
}
