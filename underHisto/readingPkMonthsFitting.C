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


vector < double > fillingPk( TString bname, TString st, int pmt, TString whichInfo)
{
  TString monthUub[] = {"dec", "jan", "feb", "mar", "abr", "may", "jun", "jul"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname = bname + pmtId+st+"Mth";

  TFile *f;
  TTree *peakInfo;
  double tmpVals = 0.;
  double tmpDer = 0.;
  vector < double > returnVals;
  int evtId = 0;
  int time = 0;
  TGraphErrors *gr = new TGraphErrors();
  TString ngr;
  TString vemFit;
  TString vemDer;
  TString pmtStr;

  for ( int month=0; month<8; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    peakInfo = (TTree*)f->Get("PeakData");

    tmpVals = 0.;
    peakInfo->SetBranchAddress(whichInfo, &tmpVals);
    
    if ( whichInfo=="peakVal" )
    {
      peakInfo->SetBranchAddress("peakValDer", &tmpDer);
      peakInfo->SetBranchAddress("graph", &gr);
    }
    peakInfo->SetBranchAddress("eventId", &evtId);
    peakInfo->SetBranchAddress("timeEvnt", &time);
    
    for( int etry=0; etry<peakInfo->GetEntries(); etry++)
    {
      peakInfo->GetEntry(etry);
      returnVals.push_back( tmpVals );
      //if ( tmpVals == 0 && tmpDer==0 && pmt==3)
        //cerr << "time: " << time << " " << evtId << " " << tmpDer << endl;
      //if ( pmt==2 && whichInfo=="peakValDer" && tmpVals < 80 && tmpVals > 60 )
        //cerr << etry << endl;
      
      //if ( pmt==1 && whichInfo=="peakVal" && evtId==62955574 ) // PMT1 vemFit = 0; vemDer = 130
      //if ( pmt==2 && whichInfo=="peakVal" && evtId==63819183 ) // PMT2 vemFit = 0; vemDer = 134
      //if ( pmt==3 && whichInfo=="peakVal" && evtId==61422934 ) // PMT2 vemFit = 0; vemDer = 138
      //if ( pmt==3 && whichInfo=="peakVal" && evtId==64188933 ) // PMT3 July vemFit = 0; vemDer = 0
      if ( pmt==3 && whichInfo=="peakVal" && evtId==62175266 ) // PMT3 Feb vemFit = 0; vemDer = 138
      {
        cerr << time << endl;
        TCanvas *c0 = canvasStyle("c0");
        c0->ResetDrawn();
        c0->cd();
        ngr.Form("%d",etry);
        vemFit.Form("%.2f", tmpVals);
        vemDer.Form("%.2f", tmpDer);
        pmtStr.Form("%d", pmt);
        gr->SetTitle("");
        //gr->SetTitle("St 863 PMT "+pmtStr+" VEM from Fit: "+vemFit+" VEM from Der.: "+vemDer);
        gr->GetXaxis()->SetTitle("[FADC/8.33 ns]");
        gr->GetXaxis()->SetRangeUser(0, 600);
        gr->GetYaxis()->SetTitle("Counts [au]");
        gr->Draw();
      
        TLegend *leg = new TLegend(0.5,0.5,0.76,0.8);
        leg->AddEntry(gr,"St 863 PMT "+pmtStr,"");
        leg->AddEntry(gr,"VEM from Fit: "+vemFit,"");
        leg->AddEntry(gr,"VEM from Der.: "+vemDer,"");
        leg->SetLineWidth(0);
        leg->SetTextSize(0.06);
        leg->Draw();

        TLine *line = new TLine(tmpDer, 0, tmpDer, 900);
        line->SetLineColor(kGreen+3);
        line->SetLineStyle(4);
        line->SetLineWidth(3);
        line->Draw();
        //c0->Print("kk.pdf");
        //c0->Print("../plots/samplePkHistoDerVem"+ngr+"Pmt"+pmt+".pdf");
        c0->Print("../plots/samplePkHistoDerVem"+ngr+"Pmt"+pmt+".pdf");
        //break;
      }
    }
  }
  peakInfo->Delete();
  f->Delete();
  return returnVals;
}

vector < int > fillingPk( TString bname, TString st, int pmt, bool whichInfo)
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

void readingPkMonthsFitting(int st)
{
  TString statId;
  statId.Form("St%d", st);
  TString basename = "uubAoPPMT";

  TPaveStats *ptstats;
  TLegend *leg;

  vector < int > timePmt1;
  vector < int > timePmt2;
  vector < int > timePmt3;
  vector < double > peakPmt1;
  vector < double > peakPmt2;
  vector < double > peakPmt3;
  vector < double > peakDerPmt1;
  vector < double > peakDerPmt2;
  vector < double > peakDerPmt3;

  bool gettime = true;
  timePmt1 = fillingPk( basename, statId, 1, gettime);
  timePmt2 = fillingPk( basename, statId, 2, gettime);
  timePmt3 = fillingPk( basename, statId, 3, gettime);
  TString whinfo = "peakVal";
  peakPmt1 = fillingPk( basename, statId, 1, whinfo);
  peakPmt2 = fillingPk( basename, statId, 2, whinfo);
  peakPmt3 = fillingPk( basename, statId, 3, whinfo);
  whinfo = "peakValDer";
  peakDerPmt1 = fillingPk( basename, statId, 1, whinfo);
  peakDerPmt2 = fillingPk( basename, statId, 2, whinfo);
  peakDerPmt3 = fillingPk( basename, statId, 3, whinfo);

  int nPoints = 0;

  nPoints = timePmt1.size();
  double xtimePmt1[ nPoints ];
  double yPkPmt1[ nPoints ];
  double yPkDerPmt1[ nPoints ];
  double yPkMeanPmt1 = 0.;
  double yPkRmsPmt1 = 0.;
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt1[i] = timePmt1[i];
    yPkPmt1[i] = peakPmt1[i];
    yPkDerPmt1[i] = peakDerPmt1[i];
  }
  double meanPmt1 = getmean( peakPmt1 );
  double rmsPmt1 = getrms( peakPmt1, meanPmt1 );
  double meanDerPmt1 = getmean( peakDerPmt1 );
  double rmsDerPmt1 = getrms( peakDerPmt1, meanDerPmt1 );
  TGraph *pkPmt1 = new TGraph(nPoints,xtimePmt1,yPkPmt1);
  TGraph *pkDerPmt1 = new TGraph(nPoints,xtimePmt1,yPkDerPmt1);

  nPoints = timePmt2.size();
  double xtimePmt2[ nPoints ];
  double yPkPmt2[ nPoints ];
  double yPkDerPmt2[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt2[i] = timePmt2[i];
    yPkPmt2[i] = peakPmt2[i];
    yPkDerPmt2[i] = peakDerPmt2[i];
  }
  double meanPmt2 = getmean( peakPmt2 );
  double rmsPmt2 = getrms( peakPmt2, meanPmt2 );
  double meanDerPmt2 = getmean( peakDerPmt2 );
  double rmsDerPmt2 = getrms( peakDerPmt2, meanDerPmt2 );
  TGraph *pkPmt2 = new TGraph(nPoints,xtimePmt2,yPkPmt2);
  TGraph *pkDerPmt2 = new TGraph(nPoints,xtimePmt2,yPkDerPmt2);

  nPoints = timePmt3.size();
  double xtimePmt3[ nPoints ];
  double yPkPmt3[ nPoints ];
  double yPkDerPmt3[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt3[i] = timePmt3[i];
    yPkPmt3[i] = peakPmt3[i];
    yPkDerPmt3[i] = peakDerPmt3[i];
  }
  double meanPmt3 = getmean( peakPmt3 );
  double rmsPmt3 = getrms( peakPmt3, meanPmt3 );
  double meanDerPmt3 = getmean( peakDerPmt3 );
  double rmsDerPmt3 = getrms( peakDerPmt3, meanDerPmt3 );
  TGraph *pkPmt3 = new TGraph(nPoints,xtimePmt3,yPkPmt3);
  TGraph *pkDerPmt3 = new TGraph(nPoints,xtimePmt3,yPkDerPmt3);

  TString avePkStr;
  TString rmsPkStr;

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  statId = "";
  statId.Form("%d", st);

  pkPmt1->SetTitle("");
  pkPmt1->GetXaxis()->SetTimeFormat("%m/%d");
  pkPmt1->GetXaxis()->SetTitle("Time since Dec. 2020 [month/day]");
  pkPmt1->GetXaxis()->SetTimeOffset(315964782,"gmt");
  pkPmt1->GetYaxis()->SetTitle("VEM Peak [FADC/8.33 ns]");
  pkPmt1->GetYaxis()->SetRangeUser(0, 230);
  pkPmt1->SetMarkerStyle(25);
  pkPmt1->SetMarkerSize(2);
  pkPmt1->SetMarkerColor(kAzure+10);
  pkPmt1->SetLineColor(kAzure+10);
  pkPmt1->SetLineWidth(2);
  histoStyle(pkPmt1);
  pkPmt1->Draw("AP same");

  pkDerPmt1->SetMarkerStyle(32);
  pkDerPmt1->SetMarkerSize(2);
  pkDerPmt1->SetMarkerColor(kOrange+10);
  pkDerPmt1->SetLineColor(kOrange+10);
  pkDerPmt1->SetLineWidth(2);
  pkDerPmt1->Draw("P same");

  avePkStr.Form("%.2f", meanPmt1);
  rmsPkStr.Form("%.2f", rmsPmt1);
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT1");
  leg->AddEntry(pkPmt1, "Average Peak-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  avePkStr.Form("%.2f", meanDerPmt1);
  rmsPkStr.Form("%.2f", rmsDerPmt1);
  leg->AddEntry(pkDerPmt1, "Average Peak-Der.: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("../plots/uubPeakFromDerSt"+statId+"pmt1.pdf");


  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  statId = "";
  statId.Form("%d", st);

  pkPmt2->SetTitle("");
  pkPmt2->GetXaxis()->SetTimeFormat("%m/%d");
  pkPmt2->GetXaxis()->SetTitle("Time since Dec. 2020 [month/day]");
  pkPmt2->GetXaxis()->SetTimeOffset(315964782,"gmt");
  pkPmt2->GetYaxis()->SetTitle("VEM Peak [FADC/8.33 ns]");
  pkPmt2->GetYaxis()->SetRangeUser(0, 230);
  pkPmt2->SetMarkerStyle(25);
  pkPmt2->SetMarkerSize(2);
  pkPmt2->SetMarkerColor(kAzure+10);
  pkPmt2->SetLineColor(kAzure+10);
  pkPmt2->SetLineWidth(2);
  histoStyle(pkPmt2);
  pkPmt2->Draw("AP same");

  pkDerPmt2->SetMarkerStyle(32);
  pkDerPmt2->SetMarkerSize(2);
  pkDerPmt2->SetMarkerColor(kOrange+10);
  pkDerPmt2->SetLineColor(kOrange+10);
  pkDerPmt2->SetLineWidth(2);
  pkDerPmt2->Draw("P same");

  avePkStr.Form("%.2f", meanPmt2);
  rmsPkStr.Form("%.2f", rmsPmt2);
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT2");
  leg->AddEntry(pkPmt2, "Average Peak-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  avePkStr.Form("%.2f", meanDerPmt2);
  rmsPkStr.Form("%.2f", rmsDerPmt2);
  leg->AddEntry(pkDerPmt2, "Average Peak-Der.: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c2->Print("../plots/uubPeakFromDerSt"+statId+"pmt2.pdf");


  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  statId = "";
  statId.Form("%d", st);

  pkPmt3->SetTitle("");
  pkPmt3->GetXaxis()->SetTimeFormat("%m/%d");
  pkPmt3->GetXaxis()->SetTitle("Time since Dec. 2020 [month/day]");
  pkPmt3->GetXaxis()->SetTimeOffset(315964782,"gmt");
  pkPmt3->GetYaxis()->SetTitle("VEM Peak [FADC/8.33 ns]");
  pkPmt3->GetYaxis()->SetRangeUser(0, 230);
  pkPmt3->SetMarkerStyle(25);
  pkPmt3->SetMarkerSize(2);
  pkPmt3->SetMarkerColor(kAzure+10);
  pkPmt3->SetLineColor(kAzure+10);
  pkPmt3->SetLineWidth(2);
  histoStyle(pkPmt3);
  pkPmt3->Draw("AP same");

  pkDerPmt3->SetMarkerStyle(32);
  pkDerPmt3->SetMarkerSize(2);
  pkDerPmt3->SetMarkerColor(kOrange+10);
  pkDerPmt3->SetLineColor(kOrange+10);
  pkDerPmt3->SetLineWidth(2);
  pkDerPmt3->Draw("P same");

  avePkStr.Form("%.2f", meanPmt3);
  rmsPkStr.Form("%.2f", rmsPmt3);
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT3");
  leg->AddEntry(pkPmt3, "Average Peak-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  avePkStr.Form("%.2f", meanDerPmt3);
  rmsPkStr.Form("%.2f", rmsDerPmt3);
  leg->AddEntry(pkDerPmt3, "Average Peak-Der.: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c3->Print("../plots/uubPeakFromDerSt"+statId+"pmt3.pdf");
}
