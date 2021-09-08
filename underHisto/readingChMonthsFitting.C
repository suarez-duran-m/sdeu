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
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "abr", "may", "jun", "jul"};
  TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  //TString fname = bname + pmtId+st+"Mth";
  TString fname = bname + pmtId+st+"lrb35";

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
  int nMonths = 8;

  for ( int month=nMonths; month<nMonths+1; month++ )
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
      //if ( month==4 && pmt==3 && tmpVals>1680 && tmpVals<1720 )
        //cout << evtId << " " << etry << endl;
      //if ( whichInfo=="chargeVal" && pmt==2 && tmpVals>1750 && tmpVals<1850 )
      //if ( whichInfo=="chargeVal" && month==8 && pmt==1 && tmpVals>10 && tmpVals<1750 )
        //cout << monthUub[month] << " " << etry << " " << evtId << endl;
      //if ( pmt==2 && whichInfo=="chargeValDer" && tmpVals < 80 && tmpVals > 60 )
        //cerr << etry << endl;
     
      //if ( whichInfo=="chargeVal" && pmt==1 && tmpVals > 40 && tmpVals < 60 )
      //if ( pmt==2 && whichInfo=="chargeVal" && tmpVals < 1200 && month == 4 )
      //if ( pmt==3 && whichInfo=="chargeVal" && tmpVals < 1000 && month == 3 )
      //if ( whichInfo=="chargeVal" && pmt==1 && month == 5 && etry == 206 )
      //if ( pmt==3 && whichInfo=="chargeVal" && evtId==63051877 )
     /* 
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
        c0->Print("kk.C");
        cout << "MSD" << endl;
        //c0->Print("../plots/sampleChHistoDerVem"+ngr+"Pmt"+pmt+".pdf");
        break;
      } 
     */
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}


vector < int > fillingCh( TString bname, TString st, int pmt, bool whichInfo)
{
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "abr", "may", "jun", "jul"};
  TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  //TString fname = bname + pmtId+st+"Mth";
  TString fname = bname + pmtId+st+"lrb35";

  TFile *f;
  TTree *chargeInfo;
  int tmpVals = 0;
  vector < int > returnVals;
  int tmpEvtId = 0;
  int nMonths = 8;

  TString getinfo;
  if ( whichInfo )
    getinfo = "timeEvnt";
  else
    getinfo = "eventId";

  for ( int month=nMonths; month<nMonths+1; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    chargeInfo = (TTree*)f->Get("ChargeData");

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(getinfo, &tmpVals);
    chargeInfo->SetBranchAddress("eventId", &tmpEvtId);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);      
      returnVals.push_back( tmpVals );
      //if ( pmt==1 && tmpVals > 1308268818 && tmpVals < 1308441618 )
        //cout << "MSD " << etry << " " << tmpEvtId << endl;
      //if ( pmt==2 && month == 4 && etry == 277 )
        //cout << "MSD: " << tmpEvtId << endl;
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}


vector < double > failingCh( TString bname, int st, int pmt, TString whichInfo)
{
  TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString pmtId;
  TString strSt;
  strSt.Form("St%d", st);
  pmtId.Form("%d", pmt);
  TString fname = bname + pmtId+strSt+"lrb35";

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
    f = TFile::Open(fname+monthUub[month]+".root");
    chargeInfo = (TTree*)f->Get("ChargeData");

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals); 
    chargeInfo->SetBranchAddress("eventId", &evtId);
    chargeInfo->SetBranchAddress("timeEvnt", &evttime);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      if ( st==863 )
      {
        if ( pmt==1 || pmt==3 )
          if ( evttime < 1308441618 || evttime > 1314057618 )
            if ( tmpVals < 1e3 || tmpVals > 2e3 )
              returnVals.push_back( tmpVals );
      }
      else if ( st==1740 )
      {
        if ( pmt==2 || pmt==3 )
          if ( evttime < 1308441618 || evttime > 1314057618 )
            if ( tmpVals < 1e3 || tmpVals > 2e3 )
              returnVals.push_back( tmpVals );
      }
      else if ( st==1729 )
      {
        if ( tmpVals < 1e3 || tmpVals > 2.1e3 )
          returnVals.push_back( tmpVals );
      }
      else
        if ( tmpVals < 1e3 || tmpVals > 2.1e3 )
          returnVals.push_back( tmpVals );
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}

vector < double > fillingDayCh( TString bname, TString st, int pmt, TString whichInfo)
{
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "abr", "may", "jun", "jul"};
  TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  //TString fname = bname + pmtId+st+"Mth";
  TString fname = bname + pmtId+st+"lrb35";

  TFile *f;
  TTree *chargeInfo;
  double tmpVals = 0.;
  double aveChDay = 0.;
  int tmpTime = 0;
  int nvals = 0;
  int startTime = 1290816018; //Dec 1st, 2020. //1290902418; //Dec 2nd, 2020.
  int diff = 0;
  vector < double > returnVals;
  int nMonths = 8;

  for ( int month=0; month<nMonths; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    chargeInfo = (TTree*)f->Get("ChargeData");

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals);
    chargeInfo->SetBranchAddress("timeEvnt", &tmpTime);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      if ( etry==0 )
      {
        diff = startTime - tmpTime;
        while ( diff < 0 )
        {
          startTime += 86400;
          diff = startTime - tmpTime; // Set time for first event
        }
      }
      if ( tmpTime >= startTime-86400 && tmpTime < startTime )
      {
        if ( tmpVals > 0 )
        {
          aveChDay += tmpVals;
          nvals++;
        }
      }
      if ( tmpTime >= startTime )
      {
        if ( nvals > 0 )
          returnVals.push_back( aveChDay/nvals );
        else
          returnVals.push_back(0);
        if ( tmpVals > 0 )
          nvals = 1;
        else
          nvals = 0;
        aveChDay = tmpVals;
        startTime += 86400;
      }
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}

vector < int > fillingDayCh( TString bname, TString st, int pmt)
{
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "abr", "may", "jun", "jul"};
  TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  //TString fname = bname + pmtId+st+"Mth";
  TString fname = bname + pmtId+st+"lrb35";

  TFile *f;
  TTree *chargeInfo;
  double aveChDay = 0.;
  int tmpTime = 0;
  int nvals = 0;
  int startTime = 1290816018; //Dec 1st, 2020. //1290902418; //Dec 2nd, 2020.
  int diff = 0;
  vector < int > returnVals;
  int nMonths = 9;

  for ( int month=0; month<nMonths; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    chargeInfo = (TTree*)f->Get("ChargeData");
    chargeInfo->SetBranchAddress("timeEvnt", &tmpTime);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      if ( etry==0 )
      {
        diff = startTime - tmpTime;
        while ( diff < 0 )
        {
          startTime += 86400;
          diff = startTime - tmpTime; // Set time for first event
        }
      }
      if ( tmpTime >= startTime )
      {
        returnVals.push_back( startTime-43200 );
        startTime += 86400;
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

void readingChMonthsFitting(int st)
{
  TString statId;
  statId.Form("St%d", st);
  //TString basename = "uubAoPPMT";
  TString basename = "uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  bool gettime = true;
  TString whinfo = "chargeVal";
  int nPoints = 0;

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

  gettime = true;
  timePmt1 = fillingCh( basename, statId, 1, gettime);
  timePmt2 = fillingCh( basename, statId, 2, gettime);
  timePmt3 = fillingCh( basename, statId, 3, gettime);
  whinfo = "chargeVal";
  chargePmt1 = fillingCh( basename, statId, 1, whinfo);
  chargePmt2 = fillingCh( basename, statId, 2, whinfo);
  chargePmt3 = fillingCh( basename, statId, 3, whinfo);
  whinfo = "chargeValDer";
  chargeDerPmt1 = fillingCh( basename, statId, 1, whinfo);
  chargeDerPmt2 = fillingCh( basename, statId, 2, whinfo);
  chargeDerPmt3 = fillingCh( basename, statId, 3, whinfo);
  whinfo = "chargeVal";
  chfailpmt1 = failingCh( basename, st, 1, whinfo);
  chfailpmt2 = failingCh( basename, st, 2, whinfo);
  chfailpmt3 = failingCh( basename, st, 3, whinfo);

  nPoints = timePmt1.size();
  double xtimePmt1[ nPoints ];
  double yChPmt1[ nPoints ];
  double yChDerPmt1[ nPoints ];
  double yPkMeanPmt1 = 0.;
  double yPkRmsPmt1 = 0.;
  int pmt1fails = chfailpmt1.size();
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
  int pmt2fails = chfailpmt2.size();
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
  int pmt3fails = chfailpmt3.size();
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
  TString strFails;
  TString strTotEvt;

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt1->SetTitle("");
  chPmt1->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt1->GetXaxis()->SetTitle("Time [month/day]");
  chPmt1->GetXaxis()->SetTimeOffset(315964800,"gmt");
  chPmt1->GetYaxis()->SetTitle("VEM-Charge [FADC]");
  chPmt1->GetYaxis()->SetRangeUser(0, 2200); //2200 for 1729; 1800 others;
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
  strFails.Form("%d", pmt1fails);
  strTotEvt.Form("%d", (int)chargePmt1.size());
  leg = new TLegend(0.15,0.25,0.52,0.5);
  leg->SetHeader("PMT1");
  leg->AddEntry(chPmt1, "From Fit (Ave.: "+avePkStr+", RMS: "+rmsPkStr+")","p");  
  avePkStr.Form("%.2f", meanDerPmt1);
  rmsPkStr.Form("%.2f", rmsDerPmt1);
  leg->AddEntry(chDerPmt1, "From Der (Ave.: "+avePkStr+", RMS: "+rmsPkStr+")","p");
  leg->AddEntry(chPmt1, "Fails: "+strFails+"/"+strTotEvt,"");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  //c1->Print("../plots/uubChargeFromDerSt"+statId+"pmt1.pdf");
  c1->Print("../plots/uubChargeFromBlrSt"+statId+"pmt1.pdf");


  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt2->SetTitle("");
  chPmt2->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt2->GetXaxis()->SetTitle("Time [month/day]");
  chPmt2->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chPmt2->GetYaxis()->SetTitle("VEM-Charge [FADC]");
  chPmt2->GetYaxis()->SetRangeUser(0, 2200); //2200 for 1729; 1800 others
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
  leg->AddEntry(chPmt2, "From Fit (Ave.: "+avePkStr+", RMS: "+rmsPkStr+")","p");
  avePkStr.Form("%.2f", meanDerPmt2);
  rmsPkStr.Form("%.2f", rmsDerPmt2);
  strFails.Form("%d", pmt2fails);
  strTotEvt.Form("%d", (int)chargePmt2.size());
  leg->AddEntry(chDerPmt2, "From Der. (Ave.: "+avePkStr+", RMS: "+rmsPkStr+")","p");
  leg->AddEntry(chPmt2, "Fails: "+strFails+"/"+strTotEvt,"");  
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  //c2->Print("../plots/uubChargeFromDerSt"+statId+"pmt2.pdf");
  c2->Print("../plots/uubChargeFromBlrSt"+statId+"pmt2.pdf");


  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  statId = "";
  statId.Form("%d", st);

  chPmt3->SetTitle("");
  chPmt3->GetXaxis()->SetTimeFormat("%m/%d");
  chPmt3->GetXaxis()->SetTitle("Time [month/day]");
  chPmt3->GetXaxis()->SetTimeOffset(315964782,"gmt");
  chPmt3->GetYaxis()->SetTitle("VEM-Charge [FADC]");
  chPmt3->GetYaxis()->SetRangeUser(0, 2200); //2200 for 1729; 1800 others
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
  leg->AddEntry(chPmt3, "From Fit (Ave.: "+avePkStr+", RMS: "+rmsPkStr+")","p");
  avePkStr.Form("%.2f", meanDerPmt3);
  rmsPkStr.Form("%.2f", rmsDerPmt3);
  strFails.Form("%d", pmt3fails);
  strTotEvt.Form("%d", (int)chargePmt3.size());
  leg->AddEntry(chDerPmt3, "From Der. (Ave.: "+avePkStr+", RMS: "+rmsPkStr+")","p");
  leg->AddEntry(chPmt3, "Fails: "+strFails+"/"+strTotEvt,"");
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->Draw();
  //c3->Print("../plots/uubChargeFromDerSt"+statId+"pmt3.pdf");
  c3->Print("../plots/uubChargeFromBlrSt"+statId+"pmt3.pdf");

/*
  TCanvas *c4 = canvasStyle("c4");
  c4->cd();
  vector < double > chargeDayPmt1;
  vector < int > timeChDayPmt1;
  statId.Form("St%d", st);

  chargeDayPmt1 = fillingDayCh( basename, statId, 1, whinfo);
  timeChDayPmt1 = fillingDayCh( basename, statId, 1 );
  nPoints = chargeDayPmt1.size();
  double xChDayPmt1[ nPoints ];
  double yChDayPmt1[ nPoints ];
  //cout << timeChDayPmt1[0] << endl;
  for ( int i=0; i<nPoints; i++ )
  {
    xChDayPmt1[i] = timeChDayPmt1[i];
    yChDayPmt1[i] = chargeDayPmt1[i];
  }
  TGraph *chDayPmt1 = new TGraph(nPoints,xChDayPmt1,yChDayPmt1);
  chDayPmt1->SetTitle("");
  chDayPmt1->GetXaxis()->SetTimeFormat("%m/%d");
  chDayPmt1->GetXaxis()->SetTitle("Time [month/day]");
  chDayPmt1->GetXaxis()->SetTimeOffset(315964800,"gmt"); //315980000,"gmt");
  chDayPmt1->GetYaxis()->SetTitle("VEM-Charge [FADC]");
  chDayPmt1->GetYaxis()->SetRangeUser(0, 2200);
  chDayPmt1->SetMarkerStyle(25);
  chDayPmt1->SetMarkerSize(2);
  chDayPmt1->SetMarkerColor(kAzure+10);
  chDayPmt1->SetLineColor(kAzure+10);
  chDayPmt1->SetLineWidth(2);
  histoStyle(chDayPmt1);
  chDayPmt1->Draw("AP"); 
  */
}
