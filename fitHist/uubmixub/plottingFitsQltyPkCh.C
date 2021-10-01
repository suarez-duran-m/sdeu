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


void plottingFitsQltyPkCh(int st)
{
  TString stat;
  TString basename = "../aoptimeTraceSelec/uubAoPtimePMT";
  TString basenameUb = "../../nouub/fitHist/aoptimeTraceSelec/ubAoPtimePMT";
  stat.Form("St%d", st);
  TString pmtId = "1";
  TString month = "Mthdec";

  TString monthUub[] = {"dec", "jan", "feb", "mar", "abr"};
  TString monthUb[] = {"aug", "sep", "oct", "nov"};

  TFile *f; 

  TTree *histos; // For averages
  TTree *histosUb; // For averages
  TTree *histosPk; // For Peak
  TTree *histosUbPk; // For Peak
  TTree *histosCh; // For Charge
  TTree *histosUbCh; // For Charge

  double tmpChiPk = 0.;
  double tmpChiCh = 0.;
  double tmpPeak = 0.;
  double tmpCharge = 0.;
 
  TCanvas *c1 = new TCanvas("c1","c1", 1600, 900); // For Peak
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  
  c1->cd();
  TPad *padA1 = new TPad("padA1","", 0 ,0.65, .5, 1.);
  padA1->SetFillColor(0);
  padA1->SetBottomMargin(0.2);
  padA1->Draw();

  c1->cd();
  TPad *padA11 = new TPad("padA11","", 0.5 ,0.65, 1., 1.);
  padA11->SetFillColor(0);
  padA11->SetBottomMargin(0.2);
  padA11->Draw();


  c1->cd();
  TPad *padA2 = new TPad("padA2","", 0 ,0.35 , .5, 0.65);
  padA2->SetFillColor(0);
  padA2->SetBottomMargin(0.2);
  padA2->Draw();

  c1->cd();
  TPad *padA21 = new TPad("padA21","", 0.5 ,0.35 , 1., 0.65);
  padA21->SetFillColor(0);
  padA21->SetBottomMargin(0.2);
  padA21->Draw();

  c1->cd();
  TPad *padA3 = new TPad("padA3","", 0 ,0.01, .5, 0.35);
  padA3->SetFillColor(0);
  padA3->SetBottomMargin(0.2);
  padA3->Draw();

  c1->cd();
  TPad *padA31 = new TPad("padA31","", 0.5 ,0.01, 1., 0.35);
  padA31->SetFillColor(0);
  padA31->SetBottomMargin(0.2);
  padA31->Draw();

  TCanvas *c2 = new TCanvas("c2","c2", 1600, 900); // For Charge 
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  
  c2->cd();
  TPad *padB1 = new TPad("padB1","", 0 ,0.65 ,.5 ,1);
  padB1->SetFillColor(0);
  padB1->SetBottomMargin(0.2);
  padB1->Draw();

  c2->cd();
  TPad *padB11 = new TPad("padB11","", 0.5 ,0.65 ,1. ,1);
  padB11->SetFillColor(0);
  padB11->SetBottomMargin(0.2);
  padB11->Draw();

  c2->cd();
  TPad *padB2 = new TPad("padB2","", 0 ,0.35, .5 ,0.65);
  padB2->SetFillColor(0);
  padB2->SetBottomMargin(0.2);
  padB2->Draw();

  c2->cd();
  TPad *padB21 = new TPad("padB21","", .5 ,0.35, 1. ,0.65);
  padB21->SetFillColor(0);
  padB21->SetBottomMargin(0.2);
  padB21->Draw();

  c2->cd();
  TPad *padB3 = new TPad("padB3","", 0 ,0.01, .5 ,0.35);
  padB3->SetFillColor(0);
  padB3->SetBottomMargin(0.2);
  padB3->Draw();

  c2->cd();
  TPad *padB31 = new TPad("padB31","", 0.5 ,0.01, 1. ,0.35);
  padB31->SetFillColor(0);
  padB31->SetBottomMargin(0.2);
  padB31->Draw();

  TCanvas *c3 = new TCanvas("c3","c3", 1600, 900); // For AoP
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  
  c3->cd();
  TPad *padC1 = new TPad("padC1","", 0. ,0.65 ,.5 , 1);
  padC1->SetFillColor(0);
  padC1->SetBottomMargin(0.2);
  padC1->Draw();

  c3->cd();
  TPad *padC11 = new TPad("padC11","", 0.5 ,0.65 ,1. ,1);
  padC11->SetFillColor(0);
  padC11->SetBottomMargin(0.2);
  padC11->Draw();

  c3->cd();
  TPad *padC2 = new TPad("padC2","", 0., 0.35, 0.5, 0.65);
  padC2->SetFillColor(0);
  padC2->SetBottomMargin(0.2);
  padC2->Draw();

  TPad *padC21 = new TPad("padC21","", 0.5 ,0.35, 1., 0.65);
  padC21->SetFillColor(0);
  padC21->SetBottomMargin(0.2);
  padC21->Draw();

  c3->cd();
  TPad *padC3 = new TPad("padC3","", 0 ,0.01 ,0.5, 0.35);
  padC3->SetFillColor(0);
  padC3->SetBottomMargin(0.2);
  padC3->Draw();

  c3->cd();
  TPad *padC31 = new TPad("padc31","", 0.5 ,0.01 ,1., 0.35);
  padC31->SetFillColor(0);
  padC31->SetBottomMargin(0.2);
  padC31->Draw();

  TH1F *chi2ndfAvePkPmt1 = new TH1F("chi2ndfAvePkPmt1", "", 500, 0, 50);
  TH1F *chi2ndfAvePkPmt2 = new TH1F("chi2ndfAvePkPmt2", "", 500, 0, 50);
  TH1F *chi2ndfAvePkPmt3 = new TH1F("chi2ndfAvePkPmt3", "", 500, 0, 50);

  TH1F *chi2ndfAveChPmt1 = new TH1F("chi2ndfAveChPmt1", "", 500, 0, 50);
  TH1F *chi2ndfAveChPmt2 = new TH1F("chi2ndfAveChPmt2", "", 500, 0, 50);
  TH1F *chi2ndfAveChPmt3 = new TH1F("chi2ndfAveChPmt3", "", 500, 0, 50);

  TH1F *chi2UbAvePkPmt1 = new TH1F("chi2UbAvePkPmt1", "", 500, 0, 50);
  TH1F *chi2UbAvePkPmt2 = new TH1F("chi2UbAvePkPmt2", "", 500, 0, 50);
  TH1F *chi2UbAvePkPmt3 = new TH1F("chi2UbAvePkPmt3", "", 500, 0, 50);

  TH1F *chi2UbAveChPmt1 = new TH1F("chi2UbAveChPmt1", "", 500, 0, 50);
  TH1F *chi2UbAveChPmt2 = new TH1F("chi2UbAveChPmt2", "", 500, 0, 50);
  TH1F *chi2UbAveChPmt3 = new TH1F("chi2UbAveChPmt3", "", 500, 0, 50);

  TPaveStats *ptstats;

  // ===============================================
  // *** *** *** Chis Average all Months *** *** *** 
  TString tmp;

  // =======================
  // *** *** FOR UUB *** *** 
  for ( int month=0; month<5; month++ )
  {
    for ( int pmt=1; pmt<4; pmt++ )
    {
      tmp.Form("%d", pmt);
      f = TFile::Open(basename+tmp+stat+"Mth"+monthUub[month]+"chpk.root");

      histos = (TTree*)f->Get("HistForChi2");

      tmpChiPk = 0.;
      tmpChiCh = 0.;

      histos->SetBranchAddress("pkChi2", &tmpChiPk);
      histos->SetBranchAddress("peak", &tmpPeak);
      histos->SetBranchAddress("chChi2", &tmpChiCh);
      histos->SetBranchAddress("charge", &tmpCharge);

      for( int etry=0; etry<histos->GetEntries(); etry++)
      {
        histos->GetEntry(etry);

        if ( tmpChiPk > 10. || tmpPeak < 100. )
          tmpChiPk = 10.;
        if ( tmpChiCh > 10. || tmpCharge < 1000. )
          tmpChiCh = 10.;

        if ( pmt==1 )
        {
          chi2ndfAvePkPmt1->Fill( tmpChiPk, 1);
          chi2ndfAveChPmt1->Fill( tmpChiCh, 1);
        }
        else if ( pmt==2 )
        {
          chi2ndfAvePkPmt2->Fill( tmpChiPk, 1);
          chi2ndfAveChPmt2->Fill( tmpChiCh, 1);
        }
        else if ( pmt==3 )
        {
          chi2ndfAvePkPmt3->Fill( tmpChiPk, 1);
          chi2ndfAveChPmt3->Fill( tmpChiCh, 1);
        }
      }
    }
  }

  // ======================
  // *** *** FOR UB *** ***

  for ( int month=0; month<4; month++ )
  {
    for ( int pmt=1; pmt<4; pmt++ )
    {
      tmp.Form("%d", pmt);
      f = TFile::Open(basenameUb+tmp+stat+"Mth"+monthUb[month]+".root");

      histosUb = (TTree*)f->Get("HistForChi2");

      tmpChiPk = 0.;
      tmpChiCh = 0.;

      histosUb->SetBranchAddress("pkChi2", &tmpChiPk);
      histosUb->SetBranchAddress("peak", &tmpPeak);
      histosUb->SetBranchAddress("chChi2", &tmpChiCh);
      histosUb->SetBranchAddress("charge", &tmpCharge);

      for( int etry=0; etry<histos->GetEntries(); etry++)
      {
        histosUb->GetEntry(etry);

        if ( tmpChiPk > 10. ) //|| tmpPeak < 100. )
          tmpChiPk = 10.;
        if ( tmpChiCh > 10. ) //|| tmpCharge < 1000. )
          tmpChiCh = 10.;

        if ( pmt==1 )
        {
          chi2UbAvePkPmt1->Fill( tmpChiPk, 1);
          chi2UbAveChPmt1->Fill( tmpChiCh, 1);
        }
        else if ( pmt==2 )
        {
          chi2UbAvePkPmt2->Fill( tmpChiPk, 1);
          chi2UbAveChPmt2->Fill( tmpChiCh, 1);
        }
        else if ( pmt==3 )
        {
          chi2UbAvePkPmt3->Fill( tmpChiPk, 1);
          chi2UbAveChPmt3->Fill( tmpChiCh, 1);
        }
      }
    }
  }
  // =========================================================
  // *** *** *** Ploting Peak and Charge over time *** *** ***

  double tmpChiAvePk = 0.;
  double tmpChiRmsPk = 0.;

  double tmpChiAveCh = 0.;
  double tmpChiRmsCh = 0.;

  // ========================
  // *** *** FOR UUB *** ***

  int timevt = 0;
  vector < int > xtimepk[3];
  vector < double > ypeak[3];
  vector < int > xtimech[3];
  vector < double > ycharge[3];
  vector < int > xtimeaop[3];
  vector < double > yaop[3];

  for ( int month=0; month<5; month++ )
  {
    for ( int pmt=1; pmt<4; pmt++ )
    {
      tmp.Form("%d", pmt);
      f = TFile::Open(basename+tmp+stat+"Mth"+monthUub[month]+"chpk.root");

      histosPk = (TTree*)f->Get("HistForChi2");

      histosPk->SetBranchAddress("pkChi2", &tmpChiPk);
      histosPk->SetBranchAddress("peak", &tmpPeak);
      histosPk->SetBranchAddress("chChi2", &tmpChiCh);
      histosPk->SetBranchAddress("charge", &tmpCharge);
      histosPk->SetBranchAddress("evtTime", &timevt);

      if ( pmt == 1 )
      {
        tmpChiAvePk = doMean(*chi2ndfAvePkPmt1); //->GetMean();
        tmpChiRmsPk = doRms(*chi2ndfAvePkPmt1, tmpChiAvePk); //->GetRMS();
        tmpChiAveCh = doMean(*chi2ndfAveChPmt1); //->GetMean();
        tmpChiRmsCh = doRms(*chi2ndfAveChPmt1, tmpChiAveCh); //->GetRMS();
      }
      else if ( pmt == 2 )
      {
        tmpChiAvePk = doMean(*chi2ndfAvePkPmt2); //->GetMean();
        tmpChiRmsPk = doRms(*chi2ndfAvePkPmt2, tmpChiAvePk); //->GetRMS();
        tmpChiAveCh = doMean(*chi2ndfAveChPmt2); //->GetMean();
        tmpChiRmsCh = doRms(*chi2ndfAveChPmt2, tmpChiAveCh); //->GetRMS();
      }
      else if ( pmt == 3 )
      {
        tmpChiAvePk = doMean(*chi2ndfAvePkPmt3); //->GetMean();
        tmpChiRmsPk = doRms(*chi2ndfAvePkPmt3, tmpChiAvePk); //->GetRMS();
        tmpChiAveCh = doMean(*chi2ndfAveChPmt3); //->GetMean();
        tmpChiRmsCh = doRms(*chi2ndfAveChPmt3, tmpChiAveCh); //->GetRMS();
      }

      for( int etry=0; etry<histos->GetEntries(); etry++)
      {
        histosPk->GetEntry(etry);

        if ( tmpChiPk > tmpChiAvePk-tmpChiRmsPk && tmpChiPk < tmpChiAvePk+tmpChiRmsPk )
          if ( tmpPeak > 0 )
          {
            xtimepk[pmt-1].push_back( timevt );
            ypeak[pmt-1].push_back( tmpPeak );
          }

        if ( tmpChiCh > tmpChiAveCh-tmpChiRmsCh && tmpChiCh < tmpChiAveCh+tmpChiRmsCh )
          if ( tmpCharge > 0 )
          {
            xtimech[pmt-1].push_back( timevt );
            ycharge[pmt-1].push_back( tmpCharge );
          }

        if ( tmpPeak > 0 && tmpCharge > 0 )
          if ( tmpChiPk > tmpChiAvePk-tmpChiRmsPk && tmpChiPk < tmpChiAvePk+tmpChiRmsPk
            && tmpChiCh > tmpChiAveCh-tmpChiRmsCh && tmpChiCh < tmpChiAveCh+tmpChiRmsCh )
          {
            xtimeaop[pmt-1].push_back( timevt );
            yaop[pmt-1].push_back( tmpCharge/tmpPeak );
          }
      }
    }
  }

  // ======================
  // *** *** FOR UB *** ***

  timevt = 0;
  vector < int > xtimeUbpk[3];
  vector < double > yUbpeak[3];
  vector < int > xtimeUbch[3];
  vector < double > yUbcharge[3];
  vector < int > xtimeUbaop[3];
  vector < double > yUbaop[3];

  for ( int month=0; month<4; month++ )
  {
    for ( int pmt=1; pmt<4; pmt++ )
    {
      tmp.Form("%d", pmt);
      f = TFile::Open(basenameUb+tmp+stat+"Mth"+monthUb[month]+".root");

      histosUbPk = (TTree*)f->Get("HistForChi2");

      histosUbPk->SetBranchAddress("pkChi2", &tmpChiPk);
      histosUbPk->SetBranchAddress("peak", &tmpPeak);
      histosUbPk->SetBranchAddress("chChi2", &tmpChiCh);
      histosUbPk->SetBranchAddress("charge", &tmpCharge);
      histosUbPk->SetBranchAddress("evtTime", &timevt);

      if ( pmt == 1 )
      {
        tmpChiAvePk = doMean(*chi2UbAvePkPmt1); //->GetMean();
        tmpChiRmsPk = doRms(*chi2UbAvePkPmt1, tmpChiAvePk); //->GetRMS();
        tmpChiAveCh = doMean(*chi2UbAveChPmt1); //->GetMean();
        tmpChiRmsCh = doRms(*chi2UbAveChPmt1, tmpChiAveCh); //->GetRMS();
      }
      else if ( pmt == 2 )
      {
        tmpChiAvePk = doMean(*chi2UbAvePkPmt2); //->GetMean();
        tmpChiRmsPk = doRms(*chi2UbAvePkPmt2, tmpChiAvePk); //->GetRMS();
        tmpChiAveCh = doMean(*chi2UbAveChPmt2); //->GetMean();
        tmpChiRmsCh = doRms(*chi2UbAveChPmt2, tmpChiAveCh); //->GetRMS();
      }
      else if ( pmt == 3 )
      {
        tmpChiAvePk = doMean(*chi2UbAvePkPmt3); //->GetMean();
        tmpChiRmsPk = doRms(*chi2UbAvePkPmt3, tmpChiAvePk); //->GetRMS();
        tmpChiAveCh = doMean(*chi2UbAveChPmt3); //->GetMean();
        tmpChiRmsCh = doRms(*chi2UbAveChPmt3, tmpChiAveCh); //->GetRMS();
      }

      for( int etry=0; etry<histos->GetEntries(); etry++)
      {
        histosUbPk->GetEntry(etry);

        if ( tmpChiPk > tmpChiAvePk-tmpChiRmsPk && tmpChiPk < tmpChiAvePk+tmpChiRmsPk )
          if ( tmpPeak > 0 && timevt > 1596239999 && tmpPeak < 100 )
          {
            xtimeUbpk[pmt-1].push_back( timevt );
            yUbpeak[pmt-1].push_back( tmpPeak );
          }

        if ( tmpChiCh > tmpChiAveCh-tmpChiRmsCh && tmpChiCh < tmpChiAveCh+tmpChiRmsCh )
          if ( tmpCharge > 0 && timevt > 1596239999 && tmpPeak < 100 )
          {
            xtimeUbch[pmt-1].push_back( timevt );
            yUbcharge[pmt-1].push_back( tmpCharge );
          }

        if ( tmpPeak > 0 && tmpCharge > 0 && timevt > 1596239999 && tmpPeak < 100 )
          if ( tmpChiPk > tmpChiAvePk-tmpChiRmsPk && tmpChiPk < tmpChiAvePk+tmpChiRmsPk
            && tmpChiCh > tmpChiAveCh-tmpChiRmsCh && tmpChiCh < tmpChiAveCh+tmpChiRmsCh )
          {
            xtimeUbaop[pmt-1].push_back( timevt );
            yUbaop[pmt-1].push_back( tmpCharge/tmpPeak );
          }
      }
    }
  }

  // ===============================================
  // *** *** *** Plotting Peak over Time *** *** *** 

  Double_t tmpXpkpmt1[xtimepk[0].size()];
  Double_t tmpYpkpmt1[xtimepk[0].size()];
  Double_t tmpXpkpmt2[xtimepk[1].size()];
  Double_t tmpYpkpmt2[xtimepk[1].size()];
  Double_t tmpXpkpmt3[xtimepk[2].size()];
  Double_t tmpYpkpmt3[xtimepk[2].size()];

  Double_t tmpXUbpkpmt1[xtimeUbpk[0].size()];
  Double_t tmpYUbpkpmt1[xtimeUbpk[0].size()];
  Double_t tmpXUbpkpmt2[xtimeUbpk[1].size()];
  Double_t tmpYUbpkpmt2[xtimeUbpk[1].size()];
  Double_t tmpXUbpkpmt3[xtimeUbpk[2].size()];
  Double_t tmpYUbpkpmt3[xtimeUbpk[2].size()];

  int n = 0;
  int nUb = 0;
  for ( int pmt=0; pmt<3; pmt++ )
  {
    n = xtimepk[pmt].size();
    nUb = xtimeUbpk[pmt].size();

    if ( pmt==0 )
    {
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbpkpmt1[k] = xtimeUbpk[pmt][k];
        tmpYUbpkpmt1[k] = yUbpeak[pmt][k];
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXpkpmt1[k] = xtimepk[pmt][k];
        tmpYpkpmt1[k] = ypeak[pmt][k];
      }

      padA1->cd();
      TGraph *pltUbPeak1 = new TGraph(nUb, tmpXUbpkpmt1, tmpYUbpkpmt1);
      pltUbPeak1->SetTitle("");
      pltUbPeak1->SetLineColor(kBlue-3);
      pltUbPeak1->SetMarkerStyle(20);
      pltUbPeak1->SetMarkerColor(kBlue-3);
      pltUbPeak1->GetXaxis()->SetTimeDisplay(1);
      pltUbPeak1->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbPeak1->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbPeak1->GetXaxis()->SetTitleSize(0.05);
      pltUbPeak1->GetXaxis()->SetLabelSize(0.05);
      pltUbPeak1->GetYaxis()->SetTitle("Peak [FADC]");
      pltUbPeak1->GetYaxis()->SetTitleSize(0.05);
      pltUbPeak1->GetYaxis()->SetLabelSize(0.05);
      pltUbPeak1->Draw("AP");
      padA1->Update();

      padA11->cd();
      TGraph *pltPeak1 = new TGraph(n, tmpXpkpmt1, tmpYpkpmt1);
      pltPeak1->SetTitle("");
      pltPeak1->SetLineColor(kOrange+8);
      pltPeak1->SetMarkerStyle(20);
      pltPeak1->SetMarkerColor(kOrange+8);
      pltPeak1->GetXaxis()->SetTimeDisplay(1);
      pltPeak1->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltPeak1->GetXaxis()->SetTitle("Time since (month/day)");
      pltPeak1->GetXaxis()->SetTitleSize(0.05);
      pltPeak1->GetXaxis()->SetLabelSize(0.05);
      //pltPeak1->GetYaxis()->SetRangeUser(0, 157);
      pltPeak1->GetYaxis()->SetTitle("Peak [FADC]");
      pltPeak1->GetYaxis()->SetTitleSize(0.05);
      pltPeak1->GetYaxis()->SetLabelSize(0.05);
      pltPeak1->Draw("AP");
      padA11->Update();
    }
    else if ( pmt==1 )
    {
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbpkpmt2[k] = xtimeUbpk[pmt][k];
        tmpYUbpkpmt2[k] = yUbpeak[pmt][k];
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXpkpmt2[k] = xtimepk[pmt][k];
        tmpYpkpmt2[k] = ypeak[pmt][k];
      }

      padA2->cd();
      TGraph *pltUbPeak2 = new TGraph(nUb, tmpXUbpkpmt2, tmpYUbpkpmt2);
      pltUbPeak2->SetTitle("");
      pltUbPeak2->SetLineColor(kBlue-3);
      pltUbPeak2->SetMarkerStyle(20);
      pltUbPeak2->SetMarkerColor(kBlue-3);
      pltUbPeak2->GetXaxis()->SetTimeDisplay(1);
      pltUbPeak2->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbPeak2->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbPeak2->GetXaxis()->SetTitleSize(0.05);
      pltUbPeak2->GetXaxis()->SetLabelSize(0.05);
      pltUbPeak2->GetYaxis()->SetTitle("Peak [FADC]");
      pltUbPeak2->GetYaxis()->SetTitleSize(0.05);
      pltUbPeak2->GetYaxis()->SetLabelSize(0.05);
      pltUbPeak2->Draw("AP");
      padA2->Update();

      padA21->cd();
      TGraph *pltPeak2 = new TGraph(n, tmpXpkpmt2, tmpYpkpmt2);
      pltPeak2->SetTitle("");
      pltPeak2->SetLineColor(kOrange+8);
      pltPeak2->SetMarkerStyle(20);
      pltPeak2->SetMarkerColor(kOrange+8);
      pltPeak2->GetXaxis()->SetTimeDisplay(1);
      pltPeak2->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltPeak2->GetXaxis()->SetTitle("Time since (month/day)");
      pltPeak2->GetXaxis()->SetTitleSize(0.05);
      pltPeak2->GetXaxis()->SetLabelSize(0.05);
      pltPeak2->GetYaxis()->SetTitle("Peak [FADC]");
      pltPeak2->GetYaxis()->SetTitleSize(0.05);
      pltPeak2->GetYaxis()->SetLabelSize(0.05);
      pltPeak2->Draw("AP");
      padA21->Update();
    }
    else if ( pmt==2 )
    {
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbpkpmt3[k] = xtimeUbpk[pmt][k];
        tmpYUbpkpmt3[k] = yUbpeak[pmt][k];
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXpkpmt3[k] = xtimepk[pmt][k];
        tmpYpkpmt3[k] = ypeak[pmt][k];
      }
      padA3->cd();
      TGraph *pltUbPeak3 = new TGraph(nUb, tmpXUbpkpmt3, tmpYUbpkpmt3);
      pltUbPeak3->SetTitle("");
      pltUbPeak3->SetLineColor(kBlue-3);
      pltUbPeak3->SetMarkerStyle(20);
      pltUbPeak3->SetMarkerColor(kBlue-3);
      pltUbPeak3->GetXaxis()->SetTimeDisplay(1);
      pltUbPeak3->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbPeak3->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbPeak3->GetXaxis()->SetTitleSize(0.05);
      pltUbPeak3->GetXaxis()->SetLabelSize(0.05);
      pltUbPeak3->GetYaxis()->SetTitle("Peak [FADC]");
      pltUbPeak3->GetYaxis()->SetTitleSize(0.05);
      pltUbPeak3->GetYaxis()->SetLabelSize(0.05);
      pltUbPeak3->Draw("AP");
      padA3->Update();

      padA31->cd();
      TGraph *pltPeak3 = new TGraph(n, tmpXpkpmt3, tmpYpkpmt3);
      pltPeak3->SetTitle("");
      pltPeak3->SetLineColor(kOrange+8);
      pltPeak3->SetMarkerStyle(20);
      pltPeak3->SetMarkerColor(kOrange+8);
      pltPeak3->GetXaxis()->SetTimeDisplay(1);
      pltPeak3->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltPeak3->GetXaxis()->SetTitle("Time since (month/day)");
      pltPeak3->GetXaxis()->SetTitleSize(0.05);
      pltPeak3->GetXaxis()->SetLabelSize(0.05);
      pltPeak3->GetYaxis()->SetTitle("Peak [FADC]");
      pltPeak3->GetYaxis()->SetTitleSize(0.05);
      pltPeak3->GetYaxis()->SetLabelSize(0.05);
      pltPeak3->Draw("AP");
      padA31->Update();
    }
  }
  c1->Update();

  // =================================================
  // *** *** *** Plotting Charge over Time *** *** *** 

  Double_t tmpXchpmt1[xtimech[0].size()];
  Double_t tmpYchpmt1[xtimech[0].size()];
  Double_t tmpXchpmt2[xtimech[1].size()];
  Double_t tmpYchpmt2[xtimech[1].size()];
  Double_t tmpXchpmt3[xtimech[2].size()];
  Double_t tmpYchpmt3[xtimech[2].size()];

  Double_t tmpXUbchpmt1[xtimeUbch[0].size()];
  Double_t tmpYUbchpmt1[xtimeUbch[0].size()];
  Double_t tmpXUbchpmt2[xtimeUbch[1].size()];
  Double_t tmpYUbchpmt2[xtimeUbch[1].size()];
  Double_t tmpXUbchpmt3[xtimeUbch[2].size()];
  Double_t tmpYUbchpmt3[xtimeUbch[2].size()];

  for ( int pmt=0; pmt<3; pmt++ )
  {
    n = xtimech[pmt].size();
    nUb = xtimeUbch[pmt].size();

    if ( pmt==0 )
    {
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbchpmt1[k] = xtimeUbch[pmt][k];
        tmpYUbchpmt1[k] = yUbcharge[pmt][k];
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXchpmt1[k] = xtimech[pmt][k];
        tmpYchpmt1[k] = ycharge[pmt][k];
      }

      padB1->cd();
      TGraph *pltUbCharge1 = new TGraph(nUb, tmpXUbchpmt1, tmpYUbchpmt1);
      pltUbCharge1->SetTitle("");
      pltUbCharge1->SetLineColor(kBlue-3);
      pltUbCharge1->SetMarkerStyle(20);
      pltUbCharge1->SetMarkerColor(kBlue-3);
      pltUbCharge1->GetXaxis()->SetTimeDisplay(1);
      pltUbCharge1->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbCharge1->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbCharge1->GetXaxis()->SetTitleSize(0.05);
      pltUbCharge1->GetXaxis()->SetLabelSize(0.05);
      pltUbCharge1->GetYaxis()->SetTitle("Charge [FADC]");
      pltUbCharge1->GetYaxis()->SetTitleSize(0.05);
      pltUbCharge1->GetYaxis()->SetLabelSize(0.05);
      pltUbCharge1->Draw("AP");
      padB1->Update();

      padB11->cd();
      TGraph *pltCharge1 = new TGraph(n, tmpXchpmt1, tmpYchpmt1);
      pltCharge1->SetTitle("");
      pltCharge1->SetLineColor(kOrange+8);
      pltCharge1->SetMarkerStyle(20);
      pltCharge1->SetMarkerColor(kOrange+8);
      pltCharge1->GetXaxis()->SetTimeDisplay(1);
      pltCharge1->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltCharge1->GetXaxis()->SetTitle("Time since (month/day)");
      pltCharge1->GetYaxis()->SetTitle("Charge [FADC]");
      pltCharge1->Draw("AP");
      padB11->Update();
    }
    else if ( pmt==1 )
    {
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbchpmt2[k] = xtimeUbch[pmt][k];
        tmpYUbchpmt2[k] = yUbcharge[pmt][k];
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXchpmt2[k] = xtimech[pmt][k];
        tmpYchpmt2[k] = ycharge[pmt][k];
      }
      padB2->cd();
      TGraph *pltUbCharge2 = new TGraph(nUb, tmpXUbchpmt2, tmpYUbchpmt2);
      pltUbCharge2->SetTitle("");
      pltUbCharge2->SetLineColor(kBlue-3);
      pltUbCharge2->SetMarkerStyle(20);
      pltUbCharge2->SetMarkerColor(kBlue-3);
      pltUbCharge2->GetXaxis()->SetTimeDisplay(1);
      pltUbCharge2->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbCharge2->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbCharge2->GetXaxis()->SetTitleSize(0.05);
      pltUbCharge2->GetXaxis()->SetLabelSize(0.05);
      pltUbCharge2->GetYaxis()->SetTitle("Charge [FADC]");
      pltUbCharge2->GetYaxis()->SetTitleSize(0.05);
      pltUbCharge2->GetYaxis()->SetLabelSize(0.05);
      pltUbCharge2->Draw("AP");
      padB2->Update();

      padB21->cd();
      TGraph *pltCharge2 = new TGraph(n, tmpXchpmt2, tmpYchpmt2);
      pltCharge2->SetTitle("");
      pltCharge2->SetLineColor(kOrange+8);
      pltCharge2->SetMarkerStyle(20);
      pltCharge2->SetMarkerColor(kOrange+8);
      pltCharge2->GetXaxis()->SetTimeDisplay(1);
      pltCharge2->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltCharge2->GetXaxis()->SetTitle("Time since (month/day)");
      pltCharge2->GetXaxis()->SetTitleSize(0.05);
      pltCharge2->GetXaxis()->SetLabelSize(0.05);
      pltCharge2->GetYaxis()->SetTitle("Charge [FADC]");
      pltCharge2->GetYaxis()->SetTitleSize(0.05);
      pltCharge2->GetYaxis()->SetLabelSize(0.05);
      pltCharge2->Draw("AP");
      padB21->Update();
    }
    else if ( pmt==2 )
    {
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbchpmt3[k] = xtimeUbch[pmt][k];
        tmpYUbchpmt3[k] = yUbcharge[pmt][k];
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXchpmt3[k] = xtimech[pmt][k];
        tmpYchpmt3[k] = ycharge[pmt][k];
      }
      padB3->cd();
      TGraph *pltUbCharge3 = new TGraph(nUb, tmpXUbchpmt3, tmpYUbchpmt3);
      pltUbCharge3->SetTitle("");
      pltUbCharge3->SetLineColor(kBlue-3);
      pltUbCharge3->SetMarkerStyle(20);
      pltUbCharge3->SetMarkerColor(kBlue-3);
      pltUbCharge3->GetXaxis()->SetTimeDisplay(1);
      pltUbCharge3->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbCharge3->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbCharge3->GetXaxis()->SetTitleSize(0.05);
      pltUbCharge3->GetXaxis()->SetLabelSize(0.05);
      pltUbCharge3->GetYaxis()->SetTitle("Charge [FADC]");
      pltUbCharge3->GetYaxis()->SetTitleSize(0.05);
      pltUbCharge3->GetYaxis()->SetLabelSize(0.05);
      pltUbCharge3->Draw("AP");
      padB3->Update();

      padB31->cd();
      TGraph *pltCharge3 = new TGraph(n, tmpXchpmt3, tmpYchpmt3);
      pltCharge3->SetTitle("");
      pltCharge3->SetLineColor(kOrange+8);
      pltCharge3->SetMarkerStyle(20);
      pltCharge3->SetMarkerColor(kOrange+8);
      pltCharge3->GetXaxis()->SetTimeDisplay(1);
      pltCharge3->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltCharge3->GetXaxis()->SetTitle("Time since (month/day)");
      pltCharge3->GetXaxis()->SetTitleSize(0.05);
      pltCharge3->GetXaxis()->SetLabelSize(0.05);
      pltCharge3->GetYaxis()->SetTitle("Charge [FADC]");
      pltCharge3->GetYaxis()->SetTitleSize(0.05);
      pltCharge3->GetYaxis()->SetLabelSize(0.05);
      pltCharge3->Draw("AP");
      padB31->Update();
    }
  }
  c2->Update();


  Double_t tmpXaoppmt1[xtimeaop[0].size()];
  Double_t tmpYaoppmt1[xtimeaop[0].size()];
  Double_t tmpXaoppmt2[xtimeaop[1].size()];
  Double_t tmpYaoppmt2[xtimeaop[1].size()];
  Double_t tmpXaoppmt3[xtimeaop[2].size()];
  Double_t tmpYaoppmt3[xtimeaop[2].size()];

  Double_t tmpXUbaoppmt1[xtimeaop[0].size()];
  Double_t tmpYUbaoppmt1[xtimeaop[0].size()];
  Double_t tmpXUbaoppmt2[xtimeaop[1].size()];
  Double_t tmpYUbaoppmt2[xtimeaop[1].size()];
  Double_t tmpXUbaoppmt3[xtimeaop[2].size()];
  Double_t tmpYUbaoppmt3[xtimeaop[2].size()];

  for ( int pmt=0; pmt<3; pmt++ )
  {
    n = xtimeaop[pmt].size();
    nUb = xtimeUbaop[pmt].size();

    if ( pmt==0 )
    {
      /*
      int cday = 1596239999;
      int tmpcntday = 0;
      vector < double > tmpSum;
      double tmp = 0.;
      */
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbaoppmt1[k] = xtimeUbaop[pmt][k];
        tmpYUbaoppmt1[k] = yUbaop[pmt][k];
        /*
        if ( xtimeUbaop[pmt][k] > cday && xtimeUbaop[pmt][k] < (cday+86400) )
          tmpSum.push_back( yUbaop[pmt][k] );
        else
        {
          tmpXUbaoppmt1[tmpcntday] = (1000000000.*cday)/1000000000;
          tmp = doMean( tmpSum );
          //cerr << tmp << endl;
          if ( tmp > 0 )
          {
            tmpYUbaoppmt1[tmpcntday] = tmp;
            cerr << int(tmpXUbaoppmt1[tmpcntday]) << " " << tmp << endl;
          }
          else
            tmpYUbaoppmt1[tmpcntday] = 0.;
          
          tmpSum.clear();
          cday += 86400;
          tmpcntday++;
        }
        */
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXaoppmt1[k] = xtimeaop[pmt][k];
        tmpYaoppmt1[k] = yaop[pmt][k];
      }

      padC1->cd();
      TGraph *pltUbAop1 = new TGraph(nUb, tmpXUbaoppmt1, tmpYUbaoppmt1);
      pltUbAop1->SetTitle("");
      pltUbAop1->SetLineColor(kBlue-3);
      pltUbAop1->SetMarkerStyle(20);
      pltUbAop1->SetMarkerColor(kBlue-3);
      pltUbAop1->GetXaxis()->SetTimeDisplay(1);
      pltUbAop1->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbAop1->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbAop1->GetXaxis()->SetTitleSize(0.05);
      pltUbAop1->GetXaxis()->SetLabelSize(0.05);
      pltUbAop1->GetYaxis()->SetTitle("AoP [25 ns]");
      pltUbAop1->GetYaxis()->SetTitleSize(0.05);
      pltUbAop1->GetYaxis()->SetLabelSize(0.05);
      pltUbAop1->Draw("AP");
      padC1->Update();

      padC11->cd();
      TGraph *pltAop1 = new TGraph(n, tmpXaoppmt1, tmpYaoppmt1);
      pltAop1->SetTitle("");
      pltAop1->SetLineColor(kOrange+8);
      pltAop1->SetMarkerStyle(20);
      pltAop1->SetMarkerColor(kOrange+8);
      pltAop1->GetXaxis()->SetTimeDisplay(1);
      pltAop1->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltAop1->GetXaxis()->SetTitle("Time since (month/day)");
      pltAop1->GetXaxis()->SetTitleSize(0.05);
      pltAop1->GetXaxis()->SetLabelSize(0.05);
      pltAop1->GetYaxis()->SetTitle("AoP [8.33 ns]");
      pltAop1->GetYaxis()->SetTitleSize(0.05);
      pltAop1->GetYaxis()->SetLabelSize(0.05);
      pltAop1->Draw("AP");
      padC11->Update();
    
    }
    
    else if ( pmt==1 )
    {
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbaoppmt2[k] = xtimeUbaop[pmt][k];
        tmpYUbaoppmt2[k] = yUbaop[pmt][k];
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXaoppmt2[k] = xtimeaop[pmt][k];
        tmpYaoppmt2[k] = yaop[pmt][k];
      }
      padC2->cd();
      TGraph *pltUbAop2 = new TGraph(nUb, tmpXUbaoppmt2, tmpYUbaoppmt2);
      pltUbAop2->SetTitle("");
      pltUbAop2->SetLineColor(kBlue-3);
      pltUbAop2->SetMarkerStyle(20);
      pltUbAop2->SetMarkerColor(kBlue-3);
      pltUbAop2->GetXaxis()->SetTimeDisplay(1);
      pltUbAop2->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbAop2->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbAop2->GetXaxis()->SetTitleSize(0.05);
      pltUbAop2->GetXaxis()->SetLabelSize(0.05);
      pltUbAop2->GetYaxis()->SetTitle("AoP [25 ns]");
      pltUbAop2->GetYaxis()->SetTitleSize(0.05);
      pltUbAop2->GetYaxis()->SetLabelSize(0.05);
      pltUbAop2->Draw("AP");
      padC2->Update();
      
      padC21->cd();
      TGraph *pltAop2 = new TGraph(n, tmpXaoppmt2, tmpYaoppmt2);
      pltAop2->SetTitle("");
      pltAop2->SetLineColor(kOrange+8);
      pltAop2->SetMarkerStyle(20);
      pltAop2->SetMarkerColor(kOrange+8);
      pltAop2->GetXaxis()->SetTimeDisplay(1);
      pltAop2->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltAop2->GetXaxis()->SetTitle("Time since (month/day)");
      pltAop2->GetXaxis()->SetTitleSize(0.05);
      pltAop2->GetXaxis()->SetLabelSize(0.05);
      pltAop2->GetYaxis()->SetTitle("AoP [8.33 ns]");
      pltAop2->GetYaxis()->SetTitleSize(0.05);
      pltAop2->GetYaxis()->SetLabelSize(0.05);
      pltAop2->Draw("AP");
      padC21->Update();
    }
    else if ( pmt==2 )
    {
      for ( int k=0; k<nUb; k++ )
      {
        tmpXUbaoppmt3[k] = xtimeUbaop[pmt][k];
        tmpYUbaoppmt3[k] = yUbaop[pmt][k];
      }
      for ( int k=0; k<n; k++ )
      {
        tmpXaoppmt3[k] = xtimeaop[pmt][k];
        tmpYaoppmt3[k] = yaop[pmt][k];
      }
      padC3->cd();
      TGraph *pltUbAop3 = new TGraph(nUb, tmpXUbaoppmt3, tmpYUbaoppmt3);
      pltUbAop3->SetTitle("");
      pltUbAop3->SetLineColor(kBlue-3);
      pltUbAop3->SetMarkerStyle(20);
      pltUbAop3->SetMarkerColor(kBlue-3);
      pltUbAop3->GetXaxis()->SetTimeDisplay(1);
      pltUbAop3->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltUbAop3->GetXaxis()->SetTitle("Time since (month/day)");
      pltUbAop3->GetXaxis()->SetTitleSize(0.05);
      pltUbAop3->GetXaxis()->SetLabelSize(0.05);
      pltUbAop3->GetYaxis()->SetTitle("AoP [25 ns]");
      pltUbAop3->GetYaxis()->SetRangeUser(2,5);
      pltUbAop3->GetYaxis()->SetTitleSize(0.05);
      pltUbAop3->GetYaxis()->SetLabelSize(0.05);
      pltUbAop3->Draw("AP");
      padC3->Update();

      padC31->cd();
      TGraph *pltAop3 = new TGraph(n, tmpXaoppmt3, tmpYaoppmt3);
      pltAop3->SetTitle("");
      pltAop3->SetLineColor(kOrange+8);
      pltAop3->SetMarkerStyle(20);
      pltAop3->SetMarkerColor(kOrange+8);
      pltAop3->GetXaxis()->SetTimeDisplay(1);
      pltAop3->GetXaxis()->SetTimeFormat("%m/%d"); // %F 1970-01-01 00:00:00");
      pltAop3->GetXaxis()->SetTitle("Time since (month/day)");
      pltAop3->GetXaxis()->SetTitleSize(0.05);
      pltAop3->GetXaxis()->SetLabelSize(0.05);
      pltAop3->GetYaxis()->SetTitle("AoP [8.33 ns]");
      pltAop3->GetYaxis()->SetTitleSize(0.05);
      pltAop3->GetYaxis()->SetLabelSize(0.05);
      pltAop3->Draw("AP");
      padC31->Update();
    }
  }
  c3->Update();

}
