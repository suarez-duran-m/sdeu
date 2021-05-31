void plottingFitsQualities(int st)
{
  TString stat;
  TString basename = "../aoptimeTraceSelec/uubAoPtimePMT";
  stat.Form("St%d", st);
  TString pmtId = "1";
  TString month = "Mthdec";

  TFile *f = TFile::Open(basename+pmtId+stat+month+"chpk.root");

  TTree *histos = (TTree*)f->Get("HistForChi2");
  TGraphErrors *chiGrpPk = new TGraphErrors();

  vector < double > chiPk;
  vector < double > chiCh;

  double tmpChiPk = 0.;
  double tmpChiCh = 0.;

  histos->SetBranchAddress("pkHistFit", &chiGrpPk);
  histos->SetBranchAddress("pkChi2", &tmpChiPk);
  histos->SetBranchAddress("chChi2", &tmpChiCh);

  TCanvas *c1 = new TCanvas("c1","c1", 1600, 900); 
  gStyle->SetOptStat(0);
  c1->cd();
  TPad *pad1 = new TPad("pad1","", 0 ,0.5 ,0.5 ,1);
  pad1->SetFillColor(0);
  pad1->SetBottomMargin(0.2);
  pad1->Draw();

  c1->cd();
  TPad *pad11 = new TPad("pad11","", 0.5 ,0.5 ,1 ,1);
  pad11->SetFillColor(0);
  pad11->SetBottomMargin(0.2);
  pad11->Draw();

  c1->cd();
  TPad *pad2 = new TPad("pad2","", 0 ,0.05 ,0.5 ,0.5);
  pad2->SetFillColor(0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();

  c1->cd();
  TPad *pad21 = new TPad("pad21","", 0.5 ,0.05 ,1 ,0.5);
  pad21->SetFillColor(0);
  pad21->SetBottomMargin(0.2);
  pad21->Draw();

  double x, y;
  histos->GetEntry(0);
  int nPoints = chiGrpPk->GetN();
  int nbins = (chiGrpPk->GetPointX(nPoints-1) - chiGrpPk->GetPointX(0) ) / 4.; // Binwidth for Peak histograms
  TH1F *residMean = new TH1F("residMean", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));
  TH1F *residMean1 = new TH1F("residMean1", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));
  TH1F *resid2040 = new TH1F("resid2040", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));
  TH1F *resid20401 = new TH1F("resid20401", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));
  TRatioPlot *rpmean;
  TRatioPlot *rpmean1;
  TRatioPlot *rp2040;
  TRatioPlot *rp20401;
  int okmean = 0;
  int ok2040 = 0;

  for( int i=0; i<histos->GetEntries(); i++)
  {
    histos->GetEntry(i);
    if ( tmpChiPk > 2. && tmpChiPk < 8. && okmean<2 )
    {
      if ( okmean==0 )
      {
        for(int i=0; i < nPoints; ++i)
        {
          chiGrpPk->GetPoint(i, x, y);
          residMean->SetBinContent(i, y);
        } 
        pad1->cd();
        TF1 *fitmean = chiGrpPk->GetFunction("fitFcn");
        residMean->Fit(fitmean,"QR");
        residMean->GetXaxis()->SetRangeUser(0, 400);
        residMean->GetYaxis()->SetTitle("Counts [au]");
        residMean->GetXaxis()->SetTitle("Peak [FADC]");
        rpmean = new TRatioPlot(residMean);
        rpmean->SetGraphDrawOpt("L");
        rpmean->Draw();
        rpmean->GetLowerRefYaxis()->SetRangeUser(-3, 3);
        pad1->Update();
      }
      else
      {
        for(int i=0; i < nPoints; ++i)
        {
          chiGrpPk->GetPoint(i, x, y);
          residMean1->SetBinContent(i, y);
        } 
        pad11->cd();
        TF1 *fitmean1 = chiGrpPk->GetFunction("fitFcn");
        residMean1->Fit(fitmean1,"QR");
        residMean1->GetXaxis()->SetRangeUser(0, 400);
        residMean1->GetYaxis()->SetTitle("Counts [au]");
        residMean1->GetXaxis()->SetTitle("Peak [FADC]");
        rpmean1 = new TRatioPlot(residMean1);
        rpmean1->SetGraphDrawOpt("L");
        rpmean1->Draw();
        rpmean1->GetLowerRefYaxis()->SetRangeUser(-3, 3);
        pad11->Update();
      }
      okmean++;
    }
    if ( tmpChiPk > 20. && tmpChiPk < 37. && ok2040<2 )
    {
      if ( ok2040==0 )
      {
        for(int i=0; i < nPoints; ++i)
        {
          chiGrpPk->GetPoint(i, x, y);
          resid2040->SetBinContent(i, y);
        }
        pad2->cd();
        TF1 *fit2040 = chiGrpPk->GetFunction("fitFcn");
        resid2040->Fit(fit2040,"QR");
        resid2040->GetXaxis()->SetRangeUser(0, 400);
        resid2040->GetYaxis()->SetTitle("Counts [au]");
        resid2040->GetXaxis()->SetTitle("Peak [FADC]");
        rp2040 = new TRatioPlot(resid2040);
        rp2040->SetGraphDrawOpt("L");
        rp2040->Draw();
        rp2040->GetLowerRefYaxis()->SetRangeUser(-3, 3);
        pad2->Update();
      }
      else
      {
        for(int i=0; i < nPoints; ++i)
        {
          chiGrpPk->GetPoint(i, x, y);
          resid20401->SetBinContent(i, y);
        }
        pad21->cd();
        TF1 *fit20401 = chiGrpPk->GetFunction("fitFcn");
        resid20401->Fit(fit20401,"QR");
        resid20401->GetXaxis()->SetRangeUser(0, 400);
        resid20401->GetYaxis()->SetTitle("Counts [au]");
        resid20401->GetXaxis()->SetTitle("Peak [FADC]");
        rp20401 = new TRatioPlot(resid20401);
        rp20401->SetGraphDrawOpt("L");
        rp20401->Draw();
        rp20401->GetLowerRefYaxis()->SetRangeUser(-3, 3);
        pad21->Update();
      }
      ok2040++;
    }
  }
  c1->Update();

  // =============================
  // *** Calculationg Averages ***
  double aveChiPk = 0.;
  int tmpcnt = 0;
  TF1 *tmpfit;
  TH1F *tmpresid; //= new TH1F("resid20401", "", nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));
  TRatioPlot *tmpRp;
  for ( int ety=0; ety<histos->GetEntries(); ety++ )
  {
    histos->GetEntry( ety );
    if ( tmpChiPk < 20. )
    {
      tmpresid =  new TH1F("tmpresid"+to_string(ety), nbins, chiGrpPk->GetPointX(0), chiGrpPk->GetPointX(nPoints-1));
      for(int i=0; i < nPoints; ++i)
      {
        chiGrpPk->GetPoint(i, x, y);
        tmpresid->SetBinContent(i, y);
      }
      tmpfit = chiGrpPk->GetFunction("fitFcn");
      tmpresid->Fit(fit20401,"QR");
      tmpRp = new TRatioPlot(tmpresid);


      aveChiPk += tmpChiPk;
      tmpcnt++;
    }
  }
  cerr << aveChiPk/tmpcnt << endl;
}
