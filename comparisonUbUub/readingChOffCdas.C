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


vector < double > getFitVals( TString bname, TString st, int pmt, int year, TString whichInfo, bool ifoff ) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strYear;
  TString fname;
  TString strChargeData;
  pmtId.Form("%d", pmt);
  strYear.Form("%d", year);

  TFile *f;
  TTree *chargeInfo;
  double tmpVals = 0.;
  vector < double > returnVals;
  int nMonths = 0;

  if ( !ifoff )
    strChargeData = "ChargeData";
  else
    strChargeData = "charge";

  int evttime = 0;
  int prevTime = 0;

  for ( int month=nMonths; month<nMonths+1; month++ ) {
    if ( !ifoff && year!=2021 )
      fname = bname + pmtId + st + "lrb35" + monthUub[month] + strYear;
    else if ( ifoff && year!=2021 )
      fname = bname + monthUub[month] + strYear + st + "Pmt" + pmtId;
    else if ( !ifoff && year==2021 )
      fname = bname + pmtId + st + "lrb35" + monthUub[month];
    else if ( ifoff && year==2021 )
      fname = bname + monthUub[month] + st + "Pmt" + pmtId;

    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals);
    if ( ifoff && year==2021 )
      chargeInfo->SetBranchAddress("GpsTime", &evttime);

    prevTime = 0;
    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {  
      chargeInfo->GetEntry(etry);
      
      if ( ifoff && year==2021 ) {
        if ( prevTime != evttime ) {
          returnVals.push_back( tmpVals );
          prevTime = evttime;
        }
      }
      else 
        returnVals.push_back( tmpVals );
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}


vector < int > getFitVals( TString bname, TString st, int pmt, int year, bool ifoff) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strYear;
  TString fname;
  TString strChargeData;
  TString strTime;
  pmtId.Form("%d", pmt);
  strYear.Form("%d", year);
  
  TFile *f;
  TTree *chargeInfo;
  int tmpVals = 0;
  vector < int > returnVals;
  int tmpEvtId = 0;
  int nMonths = 0;
  int prevTime = 0;

  if ( !ifoff ) {
    strChargeData = "ChargeData";
    strTime = "timeEvnt";
  }
  else {
    strChargeData = "charge";
    strTime = "GpsTime";
  }

  for ( int month=nMonths; month<nMonths+1; month++ ) {
    if ( !ifoff && year!=2021 )
      fname = bname + pmtId + st + "lrb35" + monthUub[month] + strYear;
    else if ( ifoff && year!=2021 )
      fname = bname + monthUub[month] + strYear + st + "Pmt" + pmtId;
    else if ( !ifoff && year==2021 )
      fname = bname + pmtId + st + "lrb35" + monthUub[month];
    else if ( ifoff && year==2021 )
      fname = bname + monthUub[month] + st + "Pmt" + pmtId;

    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);

    tmpVals = 0.;  
    chargeInfo->SetBranchAddress(strTime, &tmpVals);
    
    prevTime = 0;    
    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
      chargeInfo->GetEntry(etry);
    
      if ( ifoff && year==2021 ) {
        if ( prevTime != tmpVals ) {
          returnVals.push_back( tmpVals );
          prevTime = tmpVals;
        }
      }
      else
        returnVals.push_back( tmpVals );
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}


vector < double > getFailsVals( TString bname, TString st, int pmt, int year, TString whichInfo, bool ifoff) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strSt;
  TString fname;
  TString strYear;
  TString strChargeData;
  pmtId.Form("%d", pmt);
  strYear.Form("%d", year);

  TFile *f;
  TTree *chargeInfo;
  double tmpVals = 0.;
  vector < double > returnVals;
  int nMonths = 0;

  int evttime = 0;
  int prevTime = 0;

  if ( !ifoff )
    strChargeData = "ChargeData";
  else
    strChargeData = "charge";

  for ( int month=nMonths; month<nMonths+1; month++ ) {
    if ( !ifoff && year!=2021 )
      fname = bname + pmtId + st + "lrb35" + monthUub[month] + strYear;
    else if ( ifoff && year!=2021 )
      fname = bname + monthUub[month] + strYear + st + "Pmt" + pmtId;
    else if ( !ifoff && year==2021 )
      fname = bname + pmtId + st + "lrb35" + monthUub[month];
    else if ( ifoff && year==2021 )
      fname = bname + monthUub[month] + st + "Pmt" + pmtId;

    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals);
    if ( ifoff && year==2021 )
      chargeInfo->SetBranchAddress("GpsTime", &evttime);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
      chargeInfo->GetEntry(etry);
      if ( ifoff && year==2021 ) {
        if ( prevTime != evttime ) {
          if ( tmpVals == 0 )
            returnVals.push_back( tmpVals );
          prevTime = evttime;
        }
      }
      else {
       if ( tmpVals == 0 )
         returnVals.push_back( tmpVals );
      }
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}

TGraph *getGraph( vector<int> time, vector<double> values) {
  int nPoints = time.size();
  double xtime[ nPoints ];
  double yvalues[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtime[i] = time[i];
    yvalues[i] = values[i];
  }

  TGraph *grp = new TGraph(nPoints,xtime,yvalues);
  return grp;
}

double getmean( vector<double> arr ) {
  double mean = 0.;
  int nb = arr.size();
  int ngoodb = 0;
    for (int i=0; i<nb; i++)
      if ( arr[i] > 0 ) { 
        mean += arr[i];
        ngoodb++;
      }
  return mean/ngoodb;
}

double getrms( vector<double> arr, double meanarr ) {
  double rms = 0.;
  int nb = arr.size();
  int ngoodb = 0;
  for (int i=0; i<nb; i++)
    if ( arr[i] > 0 ) {
      rms += (arr[i] - meanarr)*(arr[i] - meanarr);
      ngoodb++;
    }
  return sqrt(rms/ngoodb)/meanarr;
}

// =================================
// *** *** *** MAIN CODE *** *** ***

void readingChOffCdas(int st, int pmt)
{
  TString statId;
  statId.Form("St%d", st);
  TString bnCdas = "~/2021/sdeu/nouub/underHistos/ubChPkPMT";
  TString bnOffl = "~/2021/offlineCourse/Offline2020/practice/ub/aug2020/offlineUb";
  TString bnUubOffl = "~/2021/offlineCourse/Offline2020/practice/uub/Aug/offlineUub";
  TString bnUubCdas = "~/2021/sdeu/underHisto/uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strMean;
  TString strRms;
  TString strFails;
  TString strGoods;
  TString strPmt;
  bool gettime = true;
  bool ifoff = false;
  TString whinfo = "chargeVal";
  double vemFactor = 11.4;

  vector < int > time2019off;
  vector < int > time2020off;
  vector < int > time2021off;
  vector < double > qpk2019off;
  vector < double > qpk2020off;
  vector < double > qpk2021off;
  vector < double > qpkFail2019off;
  vector < double > qpkFail2020off;
  vector < double > qpkFail2021off;

  vector < int > time2019cdas;
  vector < int > time2020cdas;
  vector < int > time2021cdas;
  vector < double > qpk2019cdas;
  vector < double > qpk2020cdas;
  vector < double > qpk2021cdas;
  vector < double > qpkFail2019cdas;
  vector < double > qpkFail2020cdas;
  vector < double > qpkFail2021cdas;

  ifoff = true;
  time2019off = getFitVals( bnOffl, statId, pmt, 2019, ifoff);
  time2020off = getFitVals( bnOffl, statId, pmt, 2020, ifoff);
  time2021off = getFitVals( bnUubOffl, statId, pmt, 2021, ifoff);
  ifoff = false;
  time2019cdas = getFitVals( bnCdas, statId, pmt, 2019, ifoff);
  time2020cdas = getFitVals( bnCdas, statId, pmt, 2020, ifoff);
  time2021cdas = getFitVals( bnUubCdas, statId, pmt, 2021, ifoff);

  whinfo = "chargeVal";
  ifoff = true;
  qpk2019off = getFitVals( bnOffl, statId, pmt, 2019, whinfo, ifoff);
  qpk2020off = getFitVals( bnOffl, statId, pmt, 2020, whinfo, ifoff);
  qpk2021off = getFitVals( bnUubOffl, statId, pmt, 2021, whinfo, ifoff);
  qpkFail2019off = getFailsVals( bnOffl, statId, pmt, 2019, whinfo, ifoff);
  qpkFail2020off = getFailsVals( bnOffl, statId, pmt, 2020, whinfo, ifoff);
  qpkFail2021off = getFailsVals( bnUubOffl, statId, pmt, 2021, whinfo, ifoff);

  ifoff = false;  
  qpk2019cdas = getFitVals( bnCdas, statId, pmt, 2019, whinfo, ifoff);
  qpk2020cdas = getFitVals( bnCdas, statId, pmt, 2020, whinfo, ifoff);
  qpk2021cdas = getFitVals( bnUubCdas, statId, pmt, 2021, whinfo, ifoff);
  qpkFail2019cdas = getFailsVals( bnCdas, statId, pmt, 2019, whinfo, ifoff);
  qpkFail2020cdas = getFailsVals( bnCdas, statId, pmt, 2020, whinfo, ifoff);
  qpkFail2021cdas = getFailsVals( bnUubCdas, statId, pmt, 2021, whinfo, ifoff);

  int nQpkFails2019off = qpkFail2019off.size();
  int nQpkFails2020off = qpkFail2020off.size();
  int nQpkFails2021off = qpkFail2021off.size();

  int nQpkFails2019cdas = qpkFail2019cdas.size();
  int nQpkFails2020cdas = qpkFail2020cdas.size();
  int nQpkFails2021cdas = qpkFail2021cdas.size();

  double meanQpk2019off = getmean( qpk2019off );
  double rmsQpk2019off = getrms( qpk2019off, meanQpk2019off );
  double meanQpk2020off = getmean( qpk2020off );
  double rmsQpk2020off = getrms( qpk2020off, meanQpk2020off );
  double meanQpk2021off = getmean( qpk2021off );
  double rmsQpk2021off = getrms( qpk2021off, meanQpk2021off );

  double meanQpk2019cdas = getmean( qpk2019cdas );
  double rmsQpk2019cdas = getrms( qpk2019cdas, meanQpk2019cdas );
  double meanQpk2020cdas = getmean( qpk2020cdas );
  double rmsQpk2020cdas = getrms( qpk2020cdas, meanQpk2020cdas );
  double meanQpk2021cdas = getmean( qpk2021cdas ); 
  double rmsQpk2021cdas = getrms( qpk2021cdas, meanQpk2021cdas );
 
  TGraph *GrpQpk2019off = getGraph( time2019off, qpk2019off);
  TGraph *GrpQpk2020off = getGraph( time2020off, qpk2020off);
  TGraph *GrpQpk2021cdas = getGraph( time2021cdas, qpk2021cdas);
/*
  TCanvas *c0 = canvasStyle("c0");
  c0->cd();
  GrpQpk2019off->GetXaxis()->SetTimeFormat("%m/%d");
  GrpQpk2019off->SetMarkerStyle(25);
  GrpQpk2019off->SetMarkerSize(2);
  GrpQpk2019off->GetXaxis()->SetTimeOffset(315964800,"gmt");
  GrpQpk2019off->Draw("AP");
*/
  
  int nPointsOff = time2019off.size()+time2020off.size()+time2021off.size();
  cout << "nPointsOff OffLine: " << nPointsOff << " " << time2019off.size() << " " << time2020off.size() << " " << time2021off.size() << endl;
  double xtime[ nPointsOff ];
  double y2019off[ nPointsOff ];
  double y2020off[ nPointsOff ];
  double y2021off[ nPointsOff ];
  for ( int i=0; i<nPointsOff; i++ ) {
    if ( i < time2019off.size() ) {
      xtime[i] = i;
      y2019off[i] = qpk2019off[i];
      y2020off[i] = -10.;
      y2021off[i] = -10.;
    }
    else if ( i >=time2019off.size() && i<(time2019off.size()+time2020off.size()) ) {
      xtime[i] = i;
      y2019off[i] = -10.;
      y2020off[i] = qpk2020off[i-time2019off.size()];
      y2021off[i] = -10.;
    }
    else {
      xtime[i] = i;
      y2019off[i] = -10.;
      y2020off[i] = -10.;
      y2021off[i] = qpk2021off[i-time2019off.size()-time2020off.size()]/vemFactor;
    }
  }

  int nPointsCdas = time2019cdas.size()+time2020cdas.size()+time2021cdas.size();
  cout << "nPointsCdas cdas: " << nPointsCdas << " " << time2019cdas.size() << " " << time2020cdas.size() << " " << time2021cdas.size() << endl;
  double xtimeCdas[ nPointsCdas ];
  double y2019cdas[ nPointsCdas ];
  double y2020cdas[ nPointsCdas ];
  double y2021cdas[ nPointsCdas ];
  for ( int i=0; i<nPointsCdas; i++ ) {
    if ( i < time2019cdas.size() ) {
      xtimeCdas[i] = i;
      y2019cdas[i] = qpk2019cdas[i];
      y2020cdas[i] = -10.;
      y2021cdas[i] = -10.;
    }
    else if ( i >=time2019cdas.size() && i<(time2019cdas.size()+time2020cdas.size()) ) {
      xtimeCdas[i] = i;
      y2019cdas[i] = -10.;
      y2020cdas[i] = qpk2020cdas[i-time2019cdas.size()];
      y2021cdas[i] = -10.;
    }
    else {
      xtimeCdas[i] = i;
      y2019cdas[i] = -10.;
      y2020cdas[i] = -10.;
      y2021cdas[i] = qpk2021cdas[i-time2019cdas.size()-time2020cdas.size()]/vemFactor;
    }
  }

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  TGraph *grp2019off = new TGraph(nPointsOff,xtime,y2019off);
  TGraph *grp2020off = new TGraph(nPointsOff,xtime,y2020off);
  TGraph *grp2021off = new TGraph(nPointsOff,xtime,y2021off);

  grp2019off->SetTitle("");
  grp2019off->GetYaxis()->SetRangeUser(100,200);
  grp2019off->SetMarkerStyle(25);
  grp2019off->SetMarkerSize(2);
  grp2019off->SetMarkerColor(kAzure+10);
  grp2019off->GetYaxis()->SetTitle("OffLine (Q_{VEM}^{pk} /1.01) [FADC]");
  histoStyle(grp2019off);
  grp2019off->Draw("AP same");

  grp2020off->SetMarkerStyle(8);
  grp2020off->SetMarkerSize(2);
  grp2020off->SetMarkerColor(kAzure+10);
  grp2020off->Draw("P same");

  grp2021off->SetMarkerStyle(3);
  grp2021off->SetMarkerSize(2);
  grp2021off->SetMarkerColor(kAzure+10);
  grp2021off->Draw("P same");

  leg = new TLegend(0.55,0.6,0.92,0.9);
  strPmt.Form("%d", pmt);
  statId.Form("%d", st);
  leg->SetHeader("Station "+statId+", PMT"+strPmt);
  strMean.Form("%.2f", meanQpk2019off);
  strRms.Form("%.2f", rmsQpk2019off*100.);
  strFails.Form("%d", nQpkFails2019off);
  strGoods.Form("%ld", qpk2019off.size()); 
  leg->AddEntry(grp2019off, "2019: #mu: "+strMean+", (RMS/#mu)%: "+strRms,"p");  
  leg->AddEntry(grp2019off, "Fails: "+strFails+"/"+strGoods,"");
  strMean.Form("%.2f", meanQpk2020off);
  strRms.Form("%.2f", rmsQpk2020off*100.);
  strFails.Form("%d", nQpkFails2020off);
  strGoods.Form("%ld", qpk2020off.size()); 
  leg->AddEntry(grp2020off, "2020: #mu: "+strMean+", (RMS/#mu)%: "+strRms,"p");
  leg->AddEntry(grp2020off, "Fails: "+strFails+"/"+strGoods,"");
  strMean.Form("%.2f", meanQpk2021off);
  strRms.Form("%.2f", rmsQpk2021off*100.);
  strFails.Form("%d", nQpkFails2021off);
  strGoods.Form("%ld", qpk2021off.size());
  leg->AddEntry(grp2021off, "2021: #mu: "+strMean+", (RMS/#mu)%: "+strRms,"p");
  leg->AddEntry(grp2021off, "Fails: "+strFails+"/"+strGoods,"");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->Draw();

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  TGraph *grp2019cdas = new TGraph(nPointsCdas,xtimeCdas,y2019cdas);
  TGraph *grp2020cdas = new TGraph(nPointsCdas,xtimeCdas,y2020cdas);
  TGraph *grp2021cdas = new TGraph(nPointsCdas,xtimeCdas,y2021cdas);

  grp2019cdas->SetTitle("");
  grp2019cdas->GetYaxis()->SetRangeUser(100, 200);
  grp2019cdas->SetMarkerStyle(25);
  grp2019cdas->SetMarkerSize(2);
  grp2019cdas->SetMarkerColor(kAzure+10);
  grp2019cdas->GetYaxis()->SetTitle("Q_{VEM}^{pk} [FADC]");
  histoStyle(grp2019cdas);
  grp2019cdas->Draw("AP same");

  grp2020cdas->SetMarkerStyle(8);
  grp2020cdas->SetMarkerSize(2);
  grp2020cdas->SetMarkerColor(kAzure+10);
  grp2020cdas->Draw("P same");

  grp2021cdas->SetMarkerStyle(3);
  grp2021cdas->SetMarkerSize(2);
  grp2021cdas->SetMarkerColor(kAzure+10);
  grp2021cdas->Draw("P same");

  leg = new TLegend(0.55,0.6,0.92,0.9);
  strPmt.Form("%d", pmt);
  statId.Form("%d", st);
  leg->SetHeader("Station "+statId+", PMT"+strPmt);
  strMean.Form("%.2f", meanQpk2019cdas);
  strRms.Form("%.2f", rmsQpk2019cdas*100.);
  strFails.Form("%d", nQpkFails2019cdas);
  strGoods.Form("%ld", qpk2019cdas.size()); 
  leg->AddEntry(grp2019cdas, 
      "2019: #mu: "+strMean+", (RMS/#mu)%: "+strRms,"p");  
  leg->AddEntry(grp2019cdas, "Fails: "+strFails+"/"+strGoods,"");
  strMean.Form("%.2f", meanQpk2020cdas);
  strRms.Form("%.2f", rmsQpk2020cdas*100.);
  strFails.Form("%d", nQpkFails2020cdas);
  strGoods.Form("%ld", qpk2020cdas.size()); 
  leg->AddEntry(grp2020cdas, 
      "2020: #mu: "+strMean+", (RMS/#mu)%: "+strRms,"p");
  leg->AddEntry(grp2020cdas, "Fails: "+strFails+"/"+strGoods,"");
  strMean.Form("%.2f", meanQpk2021cdas);
  strRms.Form("%.2f", rmsQpk2021cdas*100.);
  strFails.Form("%d", nQpkFails2021cdas);
  strGoods.Form("%ld", qpk2021cdas.size());
  leg->AddEntry(grp2021cdas, 
      "2021: #mu: "+strMean+", (RMS/#mu)%: "+strRms,"p");
  leg->AddEntry(grp2021cdas, "Fails: "+strFails+"/"+strGoods,"");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->Draw();
}
