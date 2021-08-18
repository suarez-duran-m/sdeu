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


void plottingFitsQualities(int st)
{
  TString stat;
  TString basename = "../aoptimeTraceSelec/uubAoPtimePMT";
  stat.Form("St%d", st);
  TString pmtId = "1";
  TString month = "Mthdec";

  TString monthUub[] = {"dec", "jan", "feb", "mar", "abr"};

  TFile *f; 

  TTree *histos;
  TTree *histos2;
  TTree *histos3;
  TGraphErrors *chiGrpPk = new TGraphErrors();
  TGraphErrors *chiGrpCh = new TGraphErrors();

  double tmpChiPk = 0.;
  double tmpChiCh = 0.;
  double tmpPeak = 0.;
  double tmpCharge = 0.;
 
  TCanvas *c1 = new TCanvas("c1","c1", 1600, 900);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  
  c1->cd();
  TPad *padA1 = new TPad("padA1","", 0 ,0.65, 1., 1.);
  padA1->SetFillColor(0);
  padA1->SetBottomMargin(0.2);
  padA1->Draw();

  c1->cd();
  TPad *padA2 = new TPad("padA2","", 0 ,0.35 , 1., 0.65);
  padA2->SetFillColor(0);
  padA2->SetBottomMargin(0.2);
  padA2->Draw();

  c1->cd();
  TPad *padA3 = new TPad("padA3","", 0 ,0.01, 1., 0.35);
  padA3->SetFillColor(0);
  padA3->SetBottomMargin(0.2);
  padA3->Draw();


  TCanvas *c2 = new TCanvas("c2","c2", 1600, 900); 
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  
  c2->cd();
  TPad *pad1 = new TPad("pad1","", 0 ,0.66 ,0.5 ,1);
  pad1->SetFillColor(0);
  pad1->SetBottomMargin(0.2);
  pad1->Draw();

  c2->cd();
  TPad *pad11 = new TPad("pad11","", 0.5 ,0.66 ,1 ,1.);
  pad11->SetFillColor(0);
  pad11->SetBottomMargin(0.2);
  pad11->Draw();

  c2->cd();
  TPad *pad2 = new TPad("pad2","", 0 ,0.33 ,0.5 ,0.66);
  pad2->SetFillColor(0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();

  c2->cd();
  TPad *pad21 = new TPad("pad21","", 0.5 ,0.33 ,1 ,0.66);
  pad21->SetFillColor(0);
  pad21->SetBottomMargin(0.2);
  pad21->Draw();

  c2->cd();
  TPad *pad3 = new TPad("pad3","", 0 ,0.0 ,0.5 ,0.33);
  pad3->SetFillColor(0);
  pad3->SetBottomMargin(0.2);
  pad3->Draw();

  c2->cd();
  TPad *pad31 = new TPad("pad31","", 0.5 ,0.0, 1 ,0.33);
  pad31->SetFillColor(0);
  pad31->SetBottomMargin(0.2);
  pad31->Draw();


  TCanvas *c3 = new TCanvas("c3","c3", 1600, 900); 
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  
  c3->cd();
  TPad *padc1 = new TPad("padc1","", 0 ,0.66 ,0.5 ,1);
  padc1->SetFillColor(0);
  padc1->SetBottomMargin(0.2);
  padc1->Draw();

  c3->cd();
  TPad *padc11 = new TPad("padc11","", 0.5 ,0.66 ,1 ,1.);
  padc11->SetFillColor(0);
  padc11->SetBottomMargin(0.2);
  padc11->Draw();

  c3->cd();
  TPad *padc2 = new TPad("padc2","", 0 ,0.33 ,0.5 ,0.66);
  padc2->SetFillColor(0);
  padc2->SetBottomMargin(0.2);
  padc2->Draw();

  c3->cd();
  TPad *padc21 = new TPad("padc21","", 0.5 ,0.33 ,1 ,0.66);
  padc21->SetFillColor(0);
  padc21->SetBottomMargin(0.2);
  padc21->Draw();

  c3->cd();
  TPad *padc3 = new TPad("padc3","", 0 ,0.0 ,0.5 ,0.33);
  padc3->SetFillColor(0);
  padc3->SetBottomMargin(0.2);
  padc3->Draw();

  c3->cd();
  TPad *padc31 = new TPad("padc31","", 0.5 ,0.0, 1 ,0.33);
  padc31->SetFillColor(0);
  padc31->SetBottomMargin(0.2);
  padc31->Draw();

  double x, y;
  int nPoints = chiGrpPk->GetN();
  int nbins = (chiGrpPk->GetPointX(nPoints-1) - chiGrpPk->GetPointX(0) ) / 4.; // Binwidth for Peak histograms

  int nPoints2 = 0;
  int nbins2 = 0;

  TH1F *chi2ndfAvePkPmt1 = new TH1F("chi2ndfAvePkPmt1", "", 500, 0, 50);
  TH1F *chi2ndfAvePkPmt2 = new TH1F("chi2ndfAvePkPmt2", "", 500, 0, 50);
  TH1F *chi2ndfAvePkPmt3 = new TH1F("chi2ndfAvePkPmt3", "", 500, 0, 50);

  TH1F *chi2ndfAveChPmt1 = new TH1F("chi2ndfAveChPmt1", "", 500, 0, 50);
  TH1F *chi2ndfAveChPmt2 = new TH1F("chi2ndfAveChPmt2", "", 500, 0, 50);
  TH1F *chi2ndfAveChPmt3 = new TH1F("chi2ndfAveChPmt3", "", 500, 0, 50);

  TH1F *chiProbPkPmt1 = new  TH1F("chiProbPkPmt1", "", 10000, 0, 1);
  TH1F *chiProbPkPmt2 = new  TH1F("chiProbPkPmt2", "", 10000, 0, 1);
  TH1F *chiProbPkPmt3 = new  TH1F("chiProbPkPmt3", "", 10000, 0, 1);

  TH1F *chiProbChPmt1 = new  TH1F("chiProbChPmt1", "", 10000, 0, 1);
  TH1F *chiProbChPmt2 = new  TH1F("chiProbChPmt2", "", 10000, 0, 1);
  TH1F *chiProbChPmt3 = new  TH1F("chiProbChPmt3", "", 10000, 0, 1);

  TH1F *residMean; 
  TH1F *residMean1;

  TH1F *residMean2; 
  TH1F *residMean21;

  TRatioPlot *rpmean;
  TRatioPlot *rpmean1;

  TRatioPlot *rpmean2;
  TRatioPlot *rpmean21;

  int ok2040 = 0;

  TRandom3 *rdmMonth = new TRandom3();
  TDatime d;
  rdmMonth->SetSeed(d.Convert());
  TPaveStats *ptstats;

  // ==========================
  // *** Average all Months ***
  TString tmp;

  for ( int month=0; month<5; month++ )
  {
    for ( int pmt=1; pmt<4; pmt++ )
    {
      tmp.Form("%d", pmt);
      f = TFile::Open(basename+tmp+stat+"Mth"+monthUub[month]+"chpk.root");

      histos = (TTree*)f->Get("HistForChi2");
      chiGrpPk = new TGraphErrors();

      tmpChiPk = 0.;
      tmpChiCh = 0.;

      histos->SetBranchAddress("pkHistFit", &chiGrpPk);
      histos->SetBranchAddress("pkChi2", &tmpChiPk);
      histos->SetBranchAddress("peak", &tmpPeak);
      histos->SetBranchAddress("chChi2", &tmpChiCh);
      histos->SetBranchAddress("charge", &tmpCharge);

      histos->GetEntry(0);
      nPoints = chiGrpPk->GetN();
      nbins = (chiGrpPk->GetPointX(nPoints-1) - chiGrpPk->GetPointX(0) ) / 4.; // Binwidth for Peak histograms
      residMean = new TH1F("residMean", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));
      residMean1 = new TH1F("residMean1", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));

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
    
      if ( pmt==1 )
      {
        padA1->cd();
        chi2ndfAvePkPmt1->SetLineColor(kBlue);
        chi2ndfAvePkPmt1->SetLineWidth(1);
        chi2ndfAvePkPmt1->SetFillColor(kBlue);
        chi2ndfAvePkPmt1->SetFillStyle(3001);
        chi2ndfAvePkPmt1->SetLabelSize(0.06,"X");
        chi2ndfAvePkPmt1->SetLabelSize(0.06,"Y");
        chi2ndfAvePkPmt1->GetXaxis()->SetTitle("#chi^{2}/ndf");
        chi2ndfAvePkPmt1->GetXaxis()->SetTitleSize(0.06);
        chi2ndfAvePkPmt1->GetXaxis()->SetRangeUser(0, 12);
        chi2ndfAvePkPmt1->GetYaxis()->SetTitle("Counts [au]");
        chi2ndfAvePkPmt1->GetYaxis()->SetTitleOffset(0.4);
        chi2ndfAvePkPmt1->GetYaxis()->SetTitleSize(0.06);
        chi2ndfAvePkPmt1->Draw("SAMES HIST"); 

        chi2ndfAveChPmt1->SetLineColor(kOrange+10);
        chi2ndfAveChPmt1->SetLineWidth(1);
        chi2ndfAveChPmt1->SetFillColor(kOrange+10);
        chi2ndfAveChPmt1->SetFillStyle(3001);
        chi2ndfAveChPmt1->Draw("SAMES HIST");

        ptstats = new TPaveStats(0.4, 0.4, 0.5, 0.8,"brNDC");
        ptstats->SetTextColor(kBlue);
        chi2ndfAvePkPmt1->SetName("PMT1 Peak Fit");
        chi2ndfAvePkPmt1->GetListOfFunctions()->Add(ptstats);

        ptstats = new TPaveStats(0.6,0.4, 0.7,0.8,"brNDC");
        ptstats->SetTextColor(kOrange+10);
        chi2ndfAveChPmt1->SetName("PMT1 Charge Fit");
        chi2ndfAveChPmt1->GetListOfFunctions()->Add(ptstats); 

        padA1->Update();
      }
      else if ( pmt==2 )
      {
        padA2->cd();
        chi2ndfAvePkPmt2->SetLineColor(kBlue);
        chi2ndfAvePkPmt2->SetLineWidth(1);
        chi2ndfAvePkPmt2->SetFillColor(kBlue);
        chi2ndfAvePkPmt2->SetFillStyle(3001);
        chi2ndfAvePkPmt2->SetLabelSize(0.06,"X");
        chi2ndfAvePkPmt2->SetLabelSize(0.06,"Y");
        chi2ndfAvePkPmt2->GetXaxis()->SetTitle("#chi^{2}/ndf");
        chi2ndfAvePkPmt2->GetXaxis()->SetTitleSize(0.06);
        chi2ndfAvePkPmt2->GetXaxis()->SetRangeUser(0, 12);
        chi2ndfAvePkPmt2->GetYaxis()->SetTitle("Counts [au]");
        chi2ndfAvePkPmt2->GetYaxis()->SetTitleOffset(0.4);
        chi2ndfAvePkPmt2->GetYaxis()->SetTitleSize(0.06);
        chi2ndfAvePkPmt2->Draw("SAMES HIST");

        chi2ndfAveChPmt2->SetLineColor(kOrange+10);
        chi2ndfAveChPmt2->SetLineWidth(1);
        chi2ndfAveChPmt2->SetFillColor(kOrange+10);
        chi2ndfAveChPmt2->SetFillStyle(3001);
        chi2ndfAveChPmt2->Draw("SAMES HIST");

        ptstats = new TPaveStats(0.4, 0.4, 0.5, 0.8,"brNDC");
        ptstats->SetTextColor(kBlue);
        chi2ndfAvePkPmt2->SetName("PMT2 Peak Fit");
        chi2ndfAvePkPmt2->GetListOfFunctions()->Add(ptstats);

        ptstats = new TPaveStats(0.6,0.4, 0.7,0.8,"brNDC");
        ptstats->SetTextColor(kOrange+10);
        chi2ndfAveChPmt2->SetName("PMT2 Charge Fit");
        chi2ndfAveChPmt2->GetListOfFunctions()->Add(ptstats); 

        padA2->Update();
      }
      else if ( pmt==3 )
      {
        padA3->cd();
        chi2ndfAvePkPmt3->SetLineColor(kBlue);
        chi2ndfAvePkPmt3->SetLineWidth(1);
        chi2ndfAvePkPmt3->SetFillColor(kBlue);
        chi2ndfAvePkPmt3->SetFillStyle(3001);
        chi2ndfAvePkPmt3->SetLabelSize(0.06,"X");
        chi2ndfAvePkPmt3->SetLabelSize(0.06,"Y");
        chi2ndfAvePkPmt3->GetXaxis()->SetTitle("#chi^{2}/ndf");
        chi2ndfAvePkPmt3->GetXaxis()->SetTitleSize(0.06);
        chi2ndfAvePkPmt3->GetXaxis()->SetRangeUser(0, 12);
        chi2ndfAvePkPmt3->GetYaxis()->SetTitle("Counts [au]");
        chi2ndfAvePkPmt3->GetYaxis()->SetTitleOffset(0.4);
        chi2ndfAvePkPmt3->GetYaxis()->SetTitleSize(0.06);
        chi2ndfAvePkPmt3->Draw("SAMES HIST");

        chi2ndfAveChPmt3->SetLineColor(kOrange+10);
        chi2ndfAveChPmt3->SetLineWidth(1);
        chi2ndfAveChPmt3->SetFillColor(kOrange+10);
        chi2ndfAveChPmt3->SetFillStyle(3001);
        chi2ndfAveChPmt3->Draw("SAMES HIST");

        ptstats = new TPaveStats(0.4,0.4, 0.5, 0.8,"brNDC");
        ptstats->SetTextColor(kBlue);
        chi2ndfAvePkPmt3->SetName("PMT3 Peak Fit");
        chi2ndfAvePkPmt3->GetListOfFunctions()->Add(ptstats);

        ptstats = new TPaveStats(0.6, 0.4, 0.7,0.8,"brNDC");
        ptstats->SetTextColor(kOrange+10);
        chi2ndfAveChPmt3->SetName("PMT3 Charge Fit");
        chi2ndfAveChPmt3->GetListOfFunctions()->Add(ptstats); 
        
        padA3->Update();
      } 
    }
  }

  c1->Print("test.pdf");

  // =========================
  // *** Plotting Examples ***
  double tmpAvepk = 0.;
  double tmpRmspk = 0.;

  int tmpMonth = 2; //rdmMonth->Integer(5);
  vector < int > tmpInMeanpk;
  vector < int > tmpOutMeanpk;

  for ( int pmt=1; pmt<4; pmt++ )
  {
    tmp.Form("%d", pmt);
    while ( tmpMonth==0 )
      tmpMonth = rdmMonth->Integer(5);
    f = TFile::Open(basename+tmp+stat+"Mth"+monthUub[tmpMonth]+"chpk.root");

    histos2 = (TTree*)f->Get("HistForChi2");
    chiGrpPk = new TGraphErrors();

    tmpChiPk = 0.;
    
    histos2->SetBranchAddress("pkHistFit", &chiGrpPk);
    histos2->SetBranchAddress("pkChi2", &tmpChiPk);
    histos2->SetBranchAddress("peak", &tmpPeak);

    histos2->GetEntry(0);
    nPoints = chiGrpPk->GetN();
    nbins = (chiGrpPk->GetPointX(nPoints-1) - chiGrpPk->GetPointX(0) ) / 4.; // Binwidth for Peak histograms
    residMean = new TH1F("residMean", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));
    residMean1 = new TH1F("residMean1", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));

    if ( pmt==1 )
    {
      tmpAvepk = doMean(*chi2ndfAvePkPmt1); //chi2ndfAvePkPmt1->GetMean();
      tmpRmspk = doRms(*chi2ndfAvePkPmt1, tmpAvepk); //chi2ndfAvePkPmt1->GetRMS();  
    }
    else if ( pmt==2 )
    {
      tmpAvepk = doMean(*chi2ndfAvePkPmt2); //chi2ndfAvePkPmt2->GetMean();
      tmpRmspk = doRms(*chi2ndfAvePkPmt2, tmpAvepk); //chi2ndfAvePkPmt2->GetRMS();
    }
    else if ( pmt==3 )
    {
      tmpAvepk = doMean(*chi2ndfAvePkPmt3); //chi2ndfAvePkPmt3->GetMean();
      tmpRmspk = doRms(*chi2ndfAvePkPmt3, tmpAvepk); //chi2ndfAvePkPmt3->GetRMS();
    }

    tmpInMeanpk.clear();
    tmpOutMeanpk.clear();

    for ( int ntry=0; ntry<histos2->GetEntries(); ntry++ )
    {
      histos2->GetEntry( ntry );
      if ( tmpChiPk > tmpAvepk-0.5*tmpRmspk && tmpChiPk < tmpAvepk+0.5*tmpRmspk && tmpPeak > 100. )
        tmpInMeanpk.push_back( ntry );
      if ( tmpChiPk > tmpAvepk + 2.*tmpRmspk )
        tmpOutMeanpk.push_back( ntry );
    }

    // ================
    // *** For Peak ***

    histos2->GetEntry( tmpInMeanpk[ rdmMonth->Integer(tmpInMeanpk.size() ) ] );

    for(int i=0; i < nPoints; ++i)
    {
      chiGrpPk->GetPoint(i, x, y);
      residMean->SetBinContent(i, y);
    } 
    if ( pmt==1 )
      pad1->cd();
    else if ( pmt==2 )
      pad2->cd();
    else if ( pmt==3 )
      pad3->cd();

    TF1 *fitmean = chiGrpPk->GetFunction("fitFcn");
    residMean->Fit(fitmean,"QR");
    gStyle->SetOptFit(1110);
    residMean->GetXaxis()->SetTitleSize(0.06);
    residMean->GetXaxis()->SetTitle("Peak [FADC]");
    residMean->GetXaxis()->SetRangeUser(0, 400);
    residMean->GetYaxis()->SetTitle("Counts [au]");
    chiProbChPmt1->Fill(fitmean->GetProb());

    rpmean = new TRatioPlot(residMean);
    rpmean->SetGraphDrawOpt("L");
    rpmean->Draw();
    rpmean->GetLowerRefYaxis()->SetRangeUser(-3, 3);
    rpmean->GetLowerRefYaxis()->SetTitle("Residual");

    if ( pmt==1 )
      pad1->Update();
    else if ( pmt==2 )
      pad2->Update();
    else if ( pmt==3 )
      pad3->Update();
    
    histos2->GetEntry( tmpOutMeanpk[ rdmMonth->Integer(tmpOutMeanpk.size() ) ] );
    for(int i=0; i < nPoints; ++i)
    {
      chiGrpPk->GetPoint(i, x, y);
      residMean1->SetBinContent(i, y);
    }
    if ( pmt==1 )
      pad11->cd();
    else if ( pmt==2 )
      pad21->cd();
    else if ( pmt==3 )
      pad31->cd();
    TF1 *fitmean1 = chiGrpPk->GetFunction("fitFcn");
    residMean1->Fit(fitmean1,"QR");
    residMean1->GetXaxis()->SetRangeUser(0, 400);
    residMean1->GetYaxis()->SetTitle("Counts [au]");
    residMean1->GetXaxis()->SetTitleSize(0.06);
    residMean1->GetXaxis()->SetTitle("Peak [FADC]");
    rpmean1 = new TRatioPlot(residMean1);
    rpmean1->SetGraphDrawOpt("L");
    rpmean1->Draw();
    rpmean1->GetLowerRefYaxis()->SetRangeUser(-3, 3);
    rpmean1->GetLowerRefYaxis()->SetTitle("Residual");
    if ( pmt==1 )
      pad11->Update();
    else if ( pmt==2 )
      pad21->Update();
    else if ( pmt==3 )
      pad31->Update();
  }

  c2->Update();
  c2->Print("test2.pdf");

  // =================================
  // *** *** ** For Charge *** *** ***

  double tmpAvech = 0.;
  double tmpRmsch = 0.;

  vector < int > tmpInMeanch;
  vector < int > tmpOutMeanch;

  tmpMonth = 2; //rdmMonth->Integer(5);

  for ( int pmt=1; pmt<4; pmt++ )
  {
    tmp.Form("%d", pmt);
    while ( tmpMonth==0 )
      tmpMonth = rdmMonth->Integer(5);
    f = TFile::Open(basename+tmp+stat+"Mth"+monthUub[tmpMonth]+"chpk.root");

    histos3 = (TTree*)f->Get("HistForChi2");
    chiGrpCh = new TGraphErrors();

    tmpChiCh = 0.;

    histos3->SetBranchAddress("chHistFit", &chiGrpCh);
    histos3->SetBranchAddress("chChi2", &tmpChiCh);
    histos3->SetBranchAddress("charge", &tmpCharge);

    histos3->GetEntry(0);

    nPoints2 = chiGrpCh->GetN();
    nbins2 = (chiGrpCh->GetPointX(nPoints2-1) - chiGrpCh->GetPointX(0)) / 8.;

    residMean2 = new TH1F("residMean2", "", nbins2, chiGrpCh->GetPointX(0), chiGrpCh->GetPointX(nPoints2-1));
    residMean21 = new TH1F("residMean21", "", nbins2, chiGrpCh->GetPointX(0), chiGrpCh->GetPointX(nPoints2-1));

    if ( pmt==1 )
    {
      tmpAvech = doMean(*chi2ndfAveChPmt1); //->GetMean();
      tmpRmsch = doRms(*chi2ndfAveChPmt1, tmpAvech); //->GetRMS();
    }
    else if ( pmt==2 )
    {
      tmpAvech = doMean(*chi2ndfAveChPmt2); //->GetMean();
      tmpRmsch = doRms(*chi2ndfAveChPmt2, tmpAvech); //->GetRMS();
    }
    else if ( pmt==3 )
    {    
      tmpAvech = doMean(*chi2ndfAveChPmt3); //->GetMean();
      tmpRmsch = doRms(*chi2ndfAveChPmt3, tmpAvech); //chi2ndfAveChPmt3->GetRMS();
    }

    tmpInMeanch.clear();
    tmpOutMeanch.clear();

    for ( int ntry=0; ntry<histos3->GetEntries(); ntry++ )
    {
      histos3->GetEntry( ntry );
      if ( tmpChiCh > tmpAvech-0.5*tmpRmsch && tmpChiCh < tmpAvech+0.5*tmpRmsch && tmpCharge > 1000. )
        tmpInMeanch.push_back( ntry );
      if ( tmpChiCh > tmpAvech + 2.*tmpRmsch )
        tmpOutMeanch.push_back( ntry );
    }
    
    histos3->GetEntry( tmpInMeanch[ rdmMonth->Integer(tmpInMeanch.size() ) ] );

    for(int i=0; i < nPoints2; ++i)
    {
      chiGrpCh->GetPoint(i, x, y);
      residMean2->SetBinContent(i, y);
    } 
    if ( pmt==1 )
      padc1->cd();
    else if ( pmt==2 )
      padc2->cd();
    else if ( pmt==3 )
      padc3->cd();

    TF1 *fitmean2 = chiGrpCh->GetFunction("fitFcn");
    residMean2->Fit(fitmean2,"QR");
    gStyle->SetOptFit(1110);
    residMean2->GetXaxis()->SetTitleSize(0.06);
    residMean2->GetXaxis()->SetTitle("Charge [FADC]");
    residMean2->GetXaxis()->SetRangeUser(0, 3000);
    residMean2->GetYaxis()->SetTitle("Counts [au]");
    rpmean2 = new TRatioPlot(residMean2);
    rpmean2->SetGraphDrawOpt("L");
    rpmean2->Draw();
    rpmean2->GetLowerRefYaxis()->SetRangeUser(-3, 3);
    rpmean2->GetLowerRefYaxis()->SetTitle("Residual");
    if ( pmt==1 )
      padc1->Update();
    else if ( pmt==2 )
      padc2->Update();
    else if ( pmt==3 )
      padc3->Update();
    
    histos3->GetEntry( tmpOutMeanch[ rdmMonth->Integer(tmpOutMeanch.size() ) ] );
    for(int i=0; i < nPoints2; ++i)
    {
      chiGrpCh->GetPoint(i, x, y);
      residMean21->SetBinContent(i, y);
    }
    if ( pmt==1 )
      padc11->cd();
    else if ( pmt==2 )
      padc21->cd();
    else if ( pmt==3 )
      padc31->cd();
    TF1 *fitmean21 = chiGrpCh->GetFunction("fitFcn");
    residMean21->Fit(fitmean21,"QR");
    residMean21->GetXaxis()->SetRangeUser(0, 3000);
    residMean21->GetYaxis()->SetTitle("Counts [au]");
    residMean21->GetXaxis()->SetTitleSize(0.06);
    residMean21->GetXaxis()->SetTitle("Charge [FADC]");
    rpmean21 = new TRatioPlot(residMean21);
    rpmean21->SetGraphDrawOpt("L");
    rpmean21->Draw();
    rpmean21->GetLowerRefYaxis()->SetRangeUser(-3, 3);
    rpmean21->GetLowerRefYaxis()->SetTitle("Residual");
    if ( pmt==1 )
      padc11->Update();
    else if ( pmt==2 )
      padc21->Update();
    else if ( pmt==3 )
      padc31->Update();
  }
  c3->Update();
}
