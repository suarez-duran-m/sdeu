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

void fillingPk( TString bname, TString st, int pmt, TH1F *hist )
{
  TString monthUub[] = {"aug", "sep", "oct", "nov"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname = bname + pmtId+st+"Mth";

  TFile *f;
  TTree *peakInfo; 
  double tmppeak = 0.;

  for ( int month=0; month<4; month++ )
  {
    cerr << fname+monthUub[month]+".root" << endl;
    f = TFile::Open(fname+monthUub[month]+".root");
    peakInfo = (TTree*)f->Get("PeakData");

    tmppeak = 0.;
    peakInfo->SetBranchAddress("peakVal", &tmppeak);

    for( int etry=0; etry<peakInfo->GetEntries(); etry++)
    {
      peakInfo->GetEntry(etry);
      hist->Fill( tmppeak );
    }
  }
  peakInfo->Delete();
  f->Delete();
}


void fillingCh( TString bname, TString st, int pmt, TH1F *hist )
{
  TString monthUub[] = {"aug", "sep", "oct", "nov"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname = bname + pmtId+st+"Mth";

  TFile *f;
  TTree *chargeInfo; 
  double tmpcharge = 0.;

  for ( int month=0; month<4; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    chargeInfo = (TTree*)f->Get("ChargeData");

    tmpcharge = 0.;
    chargeInfo->SetBranchAddress("chargeVal", &tmpcharge);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      hist->Fill( tmpcharge );
    }
  }
  chargeInfo->Delete();
  f->Delete();
}

TH1F *fillingAoP( TString bname, TString st, int pmt )
{
  TString monthUub[] = {"aug", "sep", "oct", "nov"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname = bname + pmtId+st+"Mth";

  TFile *f;
  TTree *peakInfo;
  TTree *chargeInfo; 
  double tmppk = 0.;
  double tmpch = 0.;
  int tmptmpk = 0;
  int tmptmch = 0;
  vector < double > tmppeak;
  vector < double > tmpcharge;
  vector < int > tmptimePk;
  vector < int > tmptimeCh;

  for ( int month=0; month<4; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    peakInfo = (TTree*)f->Get("PeakData");
    chargeInfo = (TTree*)f->Get("ChargeData");

    peakInfo->SetBranchAddress("peakVal", &tmppk);
    peakInfo->SetBranchAddress("timeEvnt", &tmptmpk);
    chargeInfo->SetBranchAddress("chargeVal", &tmpch);
    chargeInfo->SetBranchAddress("timeEvnt", &tmptmch);

    for( int etry=0; etry<peakInfo->GetEntries(); etry++)
    {
      peakInfo->GetEntry(etry);
      tmppeak.push_back( tmppk );
      tmptimePk.push_back( tmptmpk );
    }
    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      tmpcharge.push_back( tmpch );
      tmptimeCh.push_back( tmpch );
    }
  }
  
  double tmpaop = 0.;
  int nevts = 0;
  int cday = 1596326400; // August 2nd, 2020, 00h:00m:00s
  int nday = 0;
  const int dday = 86400;
  const int ndays = 122;
  double xtime[ndays];
  double yaop[ndays];
  double yerr[ndays]; 

  for ( int kk=0; kk<tmptimePk.size(); kk++ )
  {
    if ( tmppeak[kk] > 0 )
      tmpaop += tmpcharge[kk]  / tmppeak[kk];  
    nevts++;

    if ( tmptimePk[kk] > cday )
    {
      xtime[nday] = cday;
      if ( cday > 1603797703 && cday < 1605395135 )
        yaop[nday] = 0;
      if ( nevts==1 )
        yaop[nday] = 0;
      else
        yaop[nday] = tmpaop / nevts;

      tmpaop = 0.;
      nevts = 0;
      nday++;
      cday += dday;
    }
  }
  
  cday = 1596326400;
  nday = 0;
  nevts = 0;
  tmpaop = 0.;
  for ( int kk=0; kk<tmptimePk.size(); kk++ )
  {
    if ( tmppeak[kk] > 0 )
      tmpaop += ( (tmpcharge[kk]  / tmppeak[kk]) - yaop[nday] )*( (tmpcharge[kk]  / tmppeak[kk]) - yaop[nday] );
    nevts++;
    if ( tmptimePk[kk] > cday )
    {
      if ( cday > 1603797703 && cday < 1605395135 )
        yerr[nday] = 0;
      if ( nevts==1 )
        yerr[nday] = 0;
      yerr[nday] = sqrt(tmpaop / nevts);

      tmpaop = 0.;
      nevts = 0;
      nday++;
      cday += dday;
    }
  }
  

  TH1F *aop = new TH1F ("aop", "", ndays-1, xtime);
  double tmp = 0.;
  for ( int kk=0; kk<ndays-1; kk++ )
    if ( yaop[kk] > 1 )
    {
      aop->SetBinContent( kk+1, yaop[kk] ); 
      aop->SetBinError( kk+1, yerr[kk] );
      cerr << pmt << " " << xtime[kk] << " " << yaop[kk] << " " << yerr[kk] << endl;
    }
    else
      cerr << pmt << " " << xtime[kk] << " " << 0 << " " << 0 << endl;

  return aop;
  peakInfo->Delete();
  chargeInfo->Delete();
  f->Delete();
}


Double_t doMean(TH1F &ave)
{
  Double_t mean = 0.;
  int tmpcnt = 0;
  for ( int k=1; k<100; k++ )
  {
    mean += ave.GetBinCenter(k)*ave.GetBinContent(k);
    tmpcnt += ave.GetBinContent(k);
  }
  return mean/tmpcnt;
}

Double_t doMean(vector< double > &ave)
{
  Double_t mean = 0.;
  int tmpcnt = 0;
  for ( int k=0; k<ave.size(); k++ )
    mean += ave[k];
  
  return mean/ave.size();
}

Double_t doRms(TH1F &ave, double mean)
{
  Double_t rms = 0.;
  int tmpcnt = 0;
  for ( int k=0; k<100; k++ )
  {
    rms += ave.GetBinContent(k)*( ave.GetBinCenter(k) - mean )*( ave.GetBinCenter(k) - mean );
    tmpcnt += ave.GetBinContent(k);
  }
  return sqrt(rms/tmpcnt);
}


// =================================
// *** *** *** MAIN CODE *** *** ***

void plottingAoP(int st)
{
  TString statId;
  statId.Form("St%d", st);
  TString basename = "ubAoPPMT";

  TPaveStats *ptstats;
  TLegend *leg;

  // ===============================================
  // *** *** *** Reading Peak and Charge *** *** ***

  // ==============================
  // *** *** Doing for Peak *** ***
  
  int nbins = 250;
  TH1F *hPkpmt1 = new TH1F("hPkpmt1", "", nbins/2, 0, nbins);
  TH1F *hPkpmt2 = new TH1F("hPkpmt2", "", nbins/2, 0, nbins);
  TH1F *hPkpmt3 = new TH1F("hPkpmt3", "", nbins/2, 0, nbins);

  fillingPk( basename, statId, 1, hPkpmt1 );
  fillingPk( basename, statId, 2, hPkpmt2 );
  fillingPk( basename, statId, 3, hPkpmt3 );

  // =========================
  // *** Plotting for Peak ***

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  statId = "";
  statId.Form("%d", st);

  hPkpmt1->GetXaxis()->SetTitle("Peak [FADC]");
  hPkpmt1->GetYaxis()->SetTitle("Counts [au]");
  hPkpmt1->GetYaxis()->SetRangeUser(0, 2200); // 760); for 863
  hPkpmt1->GetXaxis()->SetRangeUser(25, 55); // 20, 60); for 863
  hPkpmt1->SetLineColor(kBlue);
  hPkpmt1->SetLineWidth(2);
  hPkpmt1->SetFillColor(kBlue);
  hPkpmt1->SetFillStyle(3001);
  histoStyle(hPkpmt1);
  hPkpmt1->GetYaxis()->SetTitleOffset(1.1);
  hPkpmt1->Draw();

  ptstats = new TPaveStats(0.73, 0.77, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  hPkpmt1->SetName("Station "+statId+" PMT1");
  hPkpmt1->GetListOfFunctions()->Add(ptstats);

  hPkpmt2->SetLineColor(kGreen+2);
  hPkpmt2->SetLineWidth(2);
  hPkpmt2->SetFillColor(kGreen+2);
  hPkpmt2->SetFillStyle(3001);
  hPkpmt2->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.55, 0.96, 0.75,"brNDC");
  ptstats->SetTextColor(kGreen+2);
  hPkpmt2->SetName("Station "+statId+" PMT2");
  hPkpmt2->GetListOfFunctions()->Add(ptstats);

  hPkpmt3->SetLineColor(kRed+1);
  hPkpmt3->SetLineWidth(2);
  hPkpmt3->SetFillColor(kRed+1);
  hPkpmt3->SetFillStyle(3001);
  hPkpmt3->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.33, 0.96, 0.53,"brNDC");
  ptstats->SetTextColor(kRed+1);
  hPkpmt3->SetName("Station "+statId+" PMT3");
  hPkpmt3->GetListOfFunctions()->Add(ptstats);

  c1->Print("../../plots/ubPeakDistPmts"+statId+".pdf");

  // ================================
  // *** *** Doing for Charge *** ***

  nbins = 2000;
  TH1F *hChpmt1 = new TH1F("hChpmt1", "", nbins/2, 0, nbins);
  TH1F *hChpmt2 = new TH1F("hChpmt2", "", nbins/2, 0, nbins);
  TH1F *hChpmt3 = new TH1F("hChpmt3", "", nbins/2, 0, nbins);

  statId = "";
  statId.Form("St%d", st);

  fillingCh( basename, statId, 1, hChpmt1 );
  fillingCh( basename, statId, 2, hChpmt2 );
  fillingCh( basename, statId, 3, hChpmt3 );

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  statId = "";
  statId.Form("%d", st);

  hChpmt1->GetXaxis()->SetTitle("Charge [FADC*8.33 ns]");
  hChpmt1->GetYaxis()->SetTitle("Counts [au]");
  hChpmt1->GetYaxis()->SetRangeUser(0, 280); // 98); for 863
  hChpmt1->GetXaxis()->SetRangeUser(80, 200); // 80, 190); for 863
  hChpmt1->SetLineColor(kBlue);
  hChpmt1->SetLineWidth(2);
  hChpmt1->SetFillColor(kBlue);
  hChpmt1->SetFillStyle(3001);
  histoStyle(hChpmt1);
  hChpmt1->Draw();

  ptstats = new TPaveStats(0.73, 0.77, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  hChpmt1->SetName("Station "+statId+" PMT1");
  hChpmt1->GetListOfFunctions()->Add(ptstats);

  hChpmt2->SetLineColor(kGreen+2);
  hChpmt2->SetLineWidth(2);
  hChpmt2->SetFillColor(kGreen+2);
  hChpmt2->SetFillStyle(3001);
  hChpmt2->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.55, 0.96, 0.75,"brNDC");
  ptstats->SetTextColor(kGreen+2);
  hChpmt2->SetName("Station "+statId+" PMT2");
  hChpmt2->GetListOfFunctions()->Add(ptstats);

  hChpmt3->SetLineColor(kRed+1);
  hChpmt3->SetLineWidth(2);
  hChpmt3->SetFillColor(kRed+1);
  hChpmt3->SetFillStyle(3001);
  hChpmt3->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.33, 0.96, 0.53,"brNDC");
  ptstats->SetTextColor(kRed+1);
  hChpmt3->SetName("Station "+statId+" PMT3");
  hChpmt3->GetListOfFunctions()->Add(ptstats);

  c2->Print("../../plots/ubChargeDistPmts"+statId+".pdf");

  // =====================================
  // *** *** *** Doing for AoP *** *** ***

  TCanvas *c3 = canvasStyle("c3");
  statId = "";
  statId.Form("St%d", st);
                                            
  TH1F *aophistpmt1 = fillingAoP( basename, statId, 1 );
  TH1F *aophistpmt2 = fillingAoP( basename, statId, 2 );
  TH1F *aophistpmt3 = fillingAoP( basename, statId, 3 );

  aophistpmt1->SetStats(0);
  aophistpmt1->SetMarkerColor(kBlue);
  aophistpmt1->SetLineColor(kBlue);
  aophistpmt1->SetMarkerSize(1.5);
  aophistpmt1->SetMarkerStyle(8);
  aophistpmt1->GetXaxis()->SetTimeDisplay(1);
  aophistpmt1->GetXaxis()->SetTimeFormat("%d/%m");
  aophistpmt1->GetXaxis()->SetTitle("Time since August 1st, 2020 [day/month]");
  aophistpmt1->GetYaxis()->SetTitle("AoP [8.33 ns]");
  aophistpmt1->GetYaxis()->SetRangeUser(2., 4.4);
  histoStyle(aophistpmt1); 
  aophistpmt1->Draw("P");

  aophistpmt2->SetStats(0);
  aophistpmt2->SetMarkerColor(kGreen+2);
  aophistpmt2->SetLineColor(kGreen+2);
  aophistpmt2->SetMarkerSize(1.5);
  aophistpmt2->SetMarkerStyle(8);
  aophistpmt2->Draw("P sames");

  aophistpmt3->SetStats(0);
  aophistpmt3->SetMarkerColor(kRed+1);
  aophistpmt3->SetLineColor(kRed+1);
  aophistpmt3->SetMarkerSize(1.5);
  aophistpmt3->SetMarkerStyle(8);
  aophistpmt3->Draw("P sames");

  leg = new TLegend(0.74, 0.16, 0.96, 0.44);
  leg->SetHeader("Station "+statId);
  leg->SetTextSize(0.06);
  leg->AddEntry(aophistpmt1, "PMT1");
  leg->AddEntry(aophistpmt2, "PMT2");
  leg->AddEntry(aophistpmt3, "PMT3");
  leg->SetTextAlign(22);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->Draw();

  c3->Print("../../plots/ubAoP"+statId+"Pmts.pdf");

}
