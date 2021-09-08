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
  TString monthUub[] = {"Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname;

  TFile *f;
  TTree *chargeInfo;
  double tmpVals = 0.;
  vector < double > returnVals;
  int evtId = 0;
  TString ngr;
  TString vem;
  TString pmtStr;
  int nMonths = 8;

  for ( int month=nMonths; month<nMonths+1; month++ )
  {
    fname = bname+monthUub[month]+"/offlineUub"+monthUub[month]+st+"Pmt"+pmtId+".root";

    f = TFile::Open(fname);
    chargeInfo = (TTree*)f->Get("charge");

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals); 
    chargeInfo->SetBranchAddress("evtId", &evtId);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      returnVals.push_back( tmpVals );
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}

vector < int > fillingCh( TString bname, TString st, int pmt, bool whichInfo)
{
  TString monthUub[] = {"Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname;
  TFile *f;
  TTree *chargeInfo;
  
  int tmpVals = 0;
  vector < int > returnVals;
  int nMonths = 8;

  TString getinfo;
  if ( whichInfo )
    getinfo = "GpsTime";
  else
   getinfo = "evtId";

  for ( int month=nMonths; month<nMonths+1; month++ )
  {
    fname = bname+monthUub[month]+"/offlineUub"+monthUub[month]+st+"Pmt"+pmtId+".root";

    f = TFile::Open(fname);
    chargeInfo = (TTree*)f->Get("charge");

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(getinfo, &tmpVals);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      returnVals.push_back( tmpVals );
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}

vector < double > failingCh( TString bname, int st, int pmt, TString whichInfo)
{
  TString monthUub[] = {"Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"};
  TString pmtId;
  TString strSt;
  strSt.Form("St%d", st);
  pmtId.Form("%d", pmt);
  TString fname;

  TFile *f;
  TTree *chargeInfo;
  double tmpVals = 0.;
  vector < double > returnVals;
  int evtId = 0;
  int evttime = 0;
  TString ngr;
  TString vem;
  TString pmtStr;
  int nMonths = 8;

  for ( int month=nMonths; month<nMonths+1; month++ )
  {
    fname = bname+monthUub[month]+"/offlineUub"+monthUub[month]+strSt+"Pmt"+pmtId+".root";

    f = TFile::Open(fname);
    chargeInfo = (TTree*)f->Get("charge");

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals); 
    chargeInfo->SetBranchAddress("evtId", &evtId);
    chargeInfo->SetBranchAddress("GpsTime", &evttime);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      if ( st!=1729 )
      {
        if ( evttime < 1308441618 || evttime > 1314057618 )
          if ( tmpVals < 1e3 || tmpVals > 2e3 )
            returnVals.push_back( tmpVals );
      }
      else
      {
        if ( tmpVals < 1e3 || tmpVals > 2e3 )
          returnVals.push_back( tmpVals );
      }
    }
  }
  chargeInfo->Delete();
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
        if ( arr[i] > 1e3 && arr[i] < 2e3 )
        {
          mean += arr[i];
          ngoodb++;
        }
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
      if ( arr[i] > 1e3 && arr[i] < 2e3 )
      {
        rms += (arr[i] - meanarr)*(arr[i] - meanarr);
        ngoodb++;
      }
    }
  return sqrt(rms/ngoodb);
}


// =================================
// *** *** *** MAIN CODE *** *** ***

void readingChMonthsFittingOff(int st)
{
  TString statId;
  statId.Form("St%d", st);
  TString basename = "/home/msd/2021/offlineCourse/Offline2020/practice/uub/";

  TPaveStats *ptstats;
  TLegend *leg;

  vector < int > timePmt1;
  vector < int > timePmt2;
  vector < int > timePmt3;
  vector < double > chargePmt1;
  vector < double > chargePmt2;
  vector < double > chargePmt3;
  vector < double > chfailpmt1;
  vector < double > chfailpmt2;
  vector < double > chfailpmt3;
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
  whinfo = "chargeVal";
  chfailpmt1 = failingCh( basename, st, 1, whinfo);
  chfailpmt2 = failingCh( basename, st, 2, whinfo);
  chfailpmt3 = failingCh( basename, st, 3, whinfo);

  int nPoints = 0;

  nPoints = timePmt1.size();
  double xtimePmt1[ nPoints ];
  double yChPmt1[ nPoints ];
  double yPkMeanPmt1 = 0.;
  double yPkRmsPmt1 = 0.;
  int pmt1fails = chfailpmt1.size();
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt1[i] = timePmt1[i];
    yChPmt1[i] = chargePmt1[i];
  }
  double meanPmt1 = getmean( chargePmt1 );
  double rmsPmt1 = getrms( chargePmt1, meanPmt1 );
  TGraph *chPmt1 = new TGraph(nPoints,xtimePmt1,yChPmt1);

  nPoints = timePmt2.size();
  double xtimePmt2[ nPoints ];
  double yChPmt2[ nPoints ];
  int pmt2fails = chfailpmt2.size();
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt2[i] = timePmt2[i];
    yChPmt2[i] = chargePmt2[i];
  }
  double meanPmt2 = getmean( chargePmt2 );
  double rmsPmt2 = getrms( chargePmt2, meanPmt2 );
  TGraph *chPmt2 = new TGraph(nPoints,xtimePmt2,yChPmt2);

  nPoints = timePmt3.size();
  double xtimePmt3[ nPoints ];
  double yChPmt3[ nPoints ];
  int pmt3fails = chfailpmt3.size();
  for ( int i=0; i<nPoints; i++ )
  {
    xtimePmt3[i] = timePmt3[i];
    yChPmt3[i] = chargePmt3[i];
  }
  double meanPmt3 = getmean( chargePmt3 );
  double rmsPmt3 = getrms( chargePmt3, meanPmt3 );
  TGraph *chPmt3 = new TGraph(nPoints,xtimePmt3,yChPmt3);

  TString aveChStr;
  TString rmsChStr;
  TString strFails;
  TString strTotEvt;

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt1->SetTitle("");
  chPmt1->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt1->GetXaxis()->SetTitle("Time [month/day]");
  chPmt1->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chPmt1->GetYaxis()->SetTitle("VEM-Charge [FADC]");
  chPmt1->GetYaxis()->SetRangeUser(0, 2200);
  chPmt1->SetMarkerStyle(25);
  chPmt1->SetMarkerSize(2);
  chPmt1->SetMarkerColor(kAzure+10);
  chPmt1->SetLineColor(kAzure+10);
  chPmt1->SetLineWidth(2);
  histoStyle(chPmt1);
  chPmt1->Draw("AP same");

  aveChStr.Form("%.2f", meanPmt1);
  rmsChStr.Form("%.2f", rmsPmt1);
  strFails.Form("%d", pmt1fails);
  strTotEvt.Form("%d", (int)chargePmt1.size());
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT1");
  leg->AddEntry(chPmt1, "From OffLine (Ave.: "+aveChStr+", RMS: "+rmsChStr+")","p");
  leg->AddEntry(chPmt1, "Fails: "+strFails+"/"+strTotEvt,"");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  //c1->Print("kk.pdf");
  
  c1->Print("../plots/uubChargeFromOffSt"+statId+"pmt1.pdf");

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt2->SetTitle("");
  chPmt2->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt2->GetXaxis()->SetTitle("Time [month/day]");
  chPmt2->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chPmt2->GetYaxis()->SetTitle("VEM-Charge [FADC]");
  chPmt2->GetYaxis()->SetRangeUser(0, 2200);
  chPmt2->SetMarkerStyle(25);
  chPmt2->SetMarkerSize(2);
  chPmt2->SetMarkerColor(kAzure+10);
  chPmt2->SetLineColor(kAzure+10);
  chPmt2->SetLineWidth(2);
  histoStyle(chPmt2);
  chPmt2->Draw("AP same");

  aveChStr.Form("%.2f", meanPmt2);
  rmsChStr.Form("%.2f", rmsPmt2);
  strFails.Form("%d", pmt2fails);
  strTotEvt.Form("%d", (int)chargePmt2.size());
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT2");
  leg->AddEntry(chPmt2, "From OffLine (Ave.: "+aveChStr+", RMS: "+rmsChStr+")","p");
  leg->AddEntry(chPmt2, "Fails: "+strFails+"/"+strTotEvt,"");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c2->Print("../plots/uubChargeFromOffSt"+statId+"pmt2.pdf");

  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt3->SetTitle("");
  chPmt3->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt3->GetXaxis()->SetTitle("Time [month/day]");
  chPmt3->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chPmt3->GetYaxis()->SetTitle("VEM-Charge [FADC]");
  chPmt3->GetYaxis()->SetRangeUser(0, 2200);
  chPmt3->SetMarkerStyle(25);
  chPmt3->SetMarkerSize(2);
  chPmt3->SetMarkerColor(kAzure+10);
  chPmt3->SetLineColor(kAzure+10);
  chPmt3->SetLineWidth(2);
  histoStyle(chPmt3);
  chPmt3->Draw("AP same");

  aveChStr.Form("%.2f", meanPmt3);
  rmsChStr.Form("%.2f", rmsPmt3);
  strFails.Form("%d", pmt3fails);
  strTotEvt.Form("%d", (int)chargePmt3.size());
  leg = new TLegend(0.15,0.31,0.52,0.5);
  leg->SetHeader("PMT3");
  leg->AddEntry(chPmt3, "From OffLine (Ave.: "+aveChStr+", RMS: "+rmsChStr+")","p");
  leg->AddEntry(chPmt3, "Fails: "+strFails+"/"+strTotEvt,"");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  c3->Print("../plots/uubChargeFromOffSt"+statId+"pmt3.pdf");
}
