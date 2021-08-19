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
  TString monthUub[] = {"Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname;

  TFile *f;
  TTree *peakInfo;
  double tmpVals = 0.;
  vector < double > returnVals;
  int evtId = 0;
  TString ngr;
  TString vem;
  TString pmtStr;

  for ( int month=0; month<8; month++ )
  {
    fname = bname+monthUub[month]+"/offlineUub"+monthUub[month]+st+"Pmt"+pmtId+".root";
    f = TFile::Open(fname);
    peakInfo = (TTree*)f->Get("peak");

    tmpVals = 0.;
    peakInfo->SetBranchAddress(whichInfo, &tmpVals); 
    peakInfo->SetBranchAddress("evtId", &evtId);

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

vector < int > fillingPk( TString bname, TString st, int pmt, bool whichInfo)
{
  TString monthUub[] = {"Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname;

  TFile *f;
  TTree *peakInfo;
  int tmpVals = 0.;
  vector < int > returnVals;

  TString getinfo;
  if ( whichInfo )
    getinfo = "GpsTime";
  else
   getinfo = "evtId";

  for ( int month=0; month<8; month++ )
  {
    fname = bname+monthUub[month]+"/offlineUub"+monthUub[month]+st+"Pmt"+pmtId+".root";
    f = TFile::Open(fname);
    peakInfo = (TTree*)f->Get("peak");

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

void readingPkMonthsFittingOff(int st)
{
  TString statId;
  statId.Form("St%d", st);
  TString basename = "/home/msd/2021/offlineCourse/Offline2020/practice/uub/";

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
  int nPoints = 0;

  nPoints = timePmt1.size();
  double xtimePmt1[ nPoints ];
  double yPkPmt1[ nPoints ];
  double yPkMeanPmt1 = 0.;
  double yPkRmsPmt1 = 0.;
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt1[i] = timePmt1[i];
    yPkPmt1[i] = peakPmt1[i];
  }
  double meanPmt1 = getmean( peakPmt1 );
  double rmsPmt1 = getrms( peakPmt1, meanPmt1 );
  TGraph *pkPmt1 = new TGraph(nPoints,xtimePmt1,yPkPmt1);

  nPoints = timePmt2.size();
  double xtimePmt2[ nPoints ];
  double yPkPmt2[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt2[i] = timePmt2[i];
    yPkPmt2[i] = peakPmt2[i];
  }
  double meanPmt2 = getmean( peakPmt2 );
  double rmsPmt2 = getrms( peakPmt2, meanPmt2 );
  TGraph *pkPmt2 = new TGraph(nPoints,xtimePmt2,yPkPmt2);

  nPoints = timePmt3.size();
  double xtimePmt3[ nPoints ];
  double yPkPmt3[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt3[i] = timePmt3[i];
    yPkPmt3[i] = peakPmt3[i];
  }
  double meanPmt3 = getmean( peakPmt3 );
  double rmsPmt3 = getrms( peakPmt3, meanPmt3 );
  TGraph *pkPmt3 = new TGraph(nPoints,xtimePmt3,yPkPmt3);

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
  pkPmt1->GetYaxis()->SetTitle("Peak-OffLine [FADC/8.33 ns]");
  pkPmt1->GetYaxis()->SetRangeUser(0, 300);
  pkPmt1->SetMarkerStyle(25);
  pkPmt1->SetMarkerSize(2);
  pkPmt1->SetMarkerColor(kAzure+10);
  pkPmt1->SetLineColor(kAzure+10);
  pkPmt1->SetLineWidth(2);
  histoStyle(pkPmt1);
  pkPmt1->Draw("AP same");

  avePkStr.Form("%.2f", meanPmt1);
  rmsPkStr.Form("%.2f", rmsPmt1);
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT1");
  leg->AddEntry(pkPmt1, "Average Peak-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("../plots/uubPeakFromOffSt"+statId+"pmt1.pdf");

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  statId = "";
  statId.Form("%d", st);

  pkPmt2->SetTitle("");
  pkPmt2->GetXaxis()->SetTimeFormat("%m/%d");
  pkPmt2->GetXaxis()->SetTitle("Time since Dec. 2020 [month/day]");
  pkPmt2->GetXaxis()->SetTimeOffset(315964782,"gmt");
  pkPmt2->GetYaxis()->SetTitle("Peak-OffLine [FADC/8.33 ns]");
  pkPmt2->GetYaxis()->SetRangeUser(0, 300);
  pkPmt2->SetMarkerStyle(25);
  pkPmt2->SetMarkerSize(2);
  pkPmt2->SetMarkerColor(kAzure+10);
  pkPmt2->SetLineColor(kAzure+10);
  pkPmt2->SetLineWidth(2);
  histoStyle(pkPmt2);
  pkPmt2->Draw("AP same");

  avePkStr.Form("%.2f", meanPmt2);
  rmsPkStr.Form("%.2f", rmsPmt2);
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT2");
  leg->AddEntry(pkPmt2, "Average Peak-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c2->Print("../plots/uubPeakFromOffSt"+statId+"pmt2.pdf");

  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  statId = "";
  statId.Form("%d", st);

  pkPmt3->SetTitle("");
  pkPmt3->GetXaxis()->SetTimeFormat("%m/%d");
  pkPmt3->GetXaxis()->SetTitle("Time since Dec. 2020 [month/day]");
  pkPmt3->GetXaxis()->SetTimeOffset(315964782,"gmt");
  pkPmt3->GetYaxis()->SetTitle("Peak-OffLine [FADC/8.33 ns]");
  pkPmt3->GetYaxis()->SetRangeUser(0, 300);
  pkPmt3->SetMarkerStyle(25);
  pkPmt3->SetMarkerSize(2);
  pkPmt3->SetMarkerColor(kAzure+10);
  pkPmt3->SetLineColor(kAzure+10);
  pkPmt3->SetLineWidth(2);
  histoStyle(pkPmt3);
  pkPmt3->Draw("AP same");

  avePkStr.Form("%.2f", meanPmt3);
  rmsPkStr.Form("%.2f", rmsPmt3);
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT3");
  leg->AddEntry(pkPmt3, "Average Peak-Fit: "+avePkStr+"; RMS: "+rmsPkStr,"p");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c3->Print("../plots/uubPeakFromOffSt"+statId+"pmt3.pdf");
}
