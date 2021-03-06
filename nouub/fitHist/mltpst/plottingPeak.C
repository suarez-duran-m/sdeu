void plottingPeak(int pmt) {

	TString pmtId;
	pmtId.Form("PMT%d", pmt);

	TFile *f = TFile::Open("ubCalibHist"+pmtId+".root");
	cout << "Reading file: " << "ubCalibHist"+pmtId+".root" <<endl;

 	TTree *histos = (TTree*)f->Get("Histograms"); 
	/*
	TH1F *hcharge = new TH1F();
	TH1F *hpeak = new TH1F();
	TH1F *peak = new TH1F(); 
	TH1F *charge = new TH1F();
	TH1F *peakCalibBl = new TH1F();
	TH1F *chargeCalibBl = new TH1F();
 	TGraphErrors *fittedPk = new TGraphErrors(); 
 	TGraphErrors *fittedCh = new TGraphErrors(); 
	*/
	TH1F *firstCntBinPk = new TH1F();
	TH1F *evtSt= new TH1F();
	int evTime;
 /*	
	histos->SetBranchAddress("receCh", &hcharge);	
	histos->SetBranchAddress("recePk", &hpeak);	
	histos->SetBranchAddress("corrPk",&peak); 	
	histos->SetBranchAddress("corrCh",&charge);
	histos->SetBranchAddress("fllHPk", &fittedPk); 
	histos->SetBranchAddress("fllHCh", &fittedCh);
	histos->SetBranchAddress("baselineHbase", &blHbase);
	histos->SetBranchAddress("baselineCalib", &blCalib);
	histos->SetBranchAddress("offSetPk", &offsetpk);
	histos->SetBranchAddress("offSetCh", &offsetch);
	//histos->SetBranchAddress("ap", &chpk);
	histos->SetBranchAddress("firstBinCntPk", &firstCntBinPk);
	histos->SetBranchAddress("firstBinCntCh", &firstCntBinCh);
  */
  histos->SetBranchAddress("firstBinCntPk", &firstCntBinPk);
  histos->SetBranchAddress("eventStat", &evtSt);
	histos->SetBranchAddress("evtTime", &evTime);

	TCanvas *c0 = new TCanvas("c0","C0",500,10,1200,800);
  TGraph *gr;

  int nx = histos->GetEntries();
  double totetrySt = 0.;
  double tmpave = 0.;
  double averageStat[19];
  double xst[19];

  for ( int st=0; st<19; st++ )
  {
    tmpave = 0.;
    totetrySt = 0.;
    for ( int etry=0; etry<nx; etry++ )
    {
      histos->GetEntry( etry );
      totetrySt = evtSt->GetBinContent( st+1 );
      tmpave = firstCntBinPk->GetBinContent(st);
    }
    xst[st] = st;
    averageStat[st] = tmpave / totetrySt;
  }

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UB Average of counts in first bin for Peak Histograms "+pmtId);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	//c1->Print("../plots/ubPeakFirstBinHBase"+pmtId+stat+".pdf");







/*
	int nx = histos->GetEntries();
	double xl[nx];
	double yl[nx];
	TString nameGraph;
	TGraph* gr[totSt];
	TMultiGraph *mg = new TMultiGraph();

	gStyle->SetOptTitle(kFALSE);
  gStyle->SetPalette(kSolar);

	for ( int st=0; st<totSt; st++ )
	{
		for	( int etry=0; etry<nx; etry++ )
		{
			histos->GetEntry( etry );
			xl[etry] = evtTime;
			yl[etry] = firstCntBinPk->GetBinContent(st);
		}
		nameGraph.Form("Station %d\n",st);
		gr[st] = new TGraph(nx,xl,yl);
		gr[st]->SetName(nameGraph);
		gr[st]->SetTitle(nameGraph);
		//gr[st]->SetTitle("UUB Offset for Peak Histograms "+pmtId);
		gr[st]->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
		gr[st]->GetXaxis()->SetTimeDisplay(1);
		gr[st]->GetXaxis()->SetTimeFormat("%m/%d");
		gr[st]->GetYaxis()->SetTitle("1. / FADC");
		gr[st]->GetYaxis()->SetTitleOffset(1.3);
		gr[st]->SetFillStyle(1000);
 	 	gr[st]->SetMarkerStyle(20);
		gr[st]->SetMarkerColor(30+st);
		gr[st]->SetMarkerSize(0.8);	
		//gr[st]->Draw("AP PLC PFC");
		gr[st]->SetDrawOption("AP");
		mg->Add( gr[st] );
	}
	
	mg->GetXaxis()->SetTimeDisplay(1);
	mg->GetXaxis()->SetTimeFormat("%m/%d");
	mg->Draw("ALP");
	gPad->BuildLegend();
	//c0->BuildLegend();
	


		//gr->Draw("AP");
*/
/*

	int nx = histos->GetEntries();
	TGraph* gr[totSt];
	double minlimHb = 1e4;
	double maxlimHb = 0.;
	double minlimCl = 1e4;
	double maxlimCl = 0.;

	TString p0nm;
	TString p1nm;
	TPaveStats *ptstats;
	TText *ptstats_LaTex;
	histos->GetEntry(0);
	c0->cd();

	// ========================================
	// *** For counts in first bin for Peak ***
	vector < vector < int > > xx; // For each St.
	vector < vector < double > > yy; // For each St.
	xx.resize( totSt );
	yy.resize( totSt );
	
	for ( Int_t etry=0; etry<nx; etry++ )
	{
		histos->GetEntry( etry );
		for ( int st=0; st<totSt; st++ )
			if ( firstCntBinPk[st] != -1 )
			{
				xx[st].push_back( evtTime[st] );
				yy[st].push_back( firstCntBinPk[st] );
			}
	}
	double xbl[nx];
	double ybl[nx];
*/
/*
	c0->cd();
	TMultiGraph *mg = new TMultiGraph();
	for ( int st=1; st<totSt-17; st++ )
	{
		for ( int tt=0; tt<xx[st].size(); tt++ )
		{
			xbl[tt] = xx[st][tt];
			ybl[tt] = yy[st][tt];
		}
		gr[st] = new TGraph(xx[st].size(),xbl,ybl);
		//gr[st]->SetTitle("UUB Offset for Peak Histograms "+pmtId);
		//gr[st]->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
		//gr[st]->GetXaxis()->SetRangeUser(xbl[0], xbl[tmplim]);
		//gr[st]->GetXaxis()->SetTimeDisplay(1);
		//gr[st]->GetXaxis()->SetTimeFormat("%m/%d");
		//gr[st]->GetYaxis()->SetTitle("1. / FADC");
		//gr[st]->GetYaxis()->SetTitleOffset(1.3);
		//gr[st]->SetFillStyle(21);
	  gr[st]->SetMarkerStyle(21);
 		gr[st]->SetDrawOption("AP");
 		gr[st]->SetLineColor(2+st);
 		gr[st]->SetLineWidth(4);
 		gr[st]->SetFillStyle(0);
		
		gr[st]->SetMarkerStyle(20);
		gr[st]->SetMarkerSize(0.8);
		gr[st]->SetMarkerColor(4+st);
		
		mg->Add(gr[st]);//, "AP");
	}
	//mg->GetXaxis()->SetRangeUser(xx[0][0], xx[0][int(xx[0].size())-1]);
	mg->GetXaxis()->SetTimeDisplay(1);
	mg->GetXaxis()->SetTimeFormat("%m/%d");
	mg->Draw("AP");
*/
	//c0->Print("../plots/uubOffsetPk"+pmtId+stat+".pdf");


/*	
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = offsetch;
	}

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Offset for Charge Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	if ( pmt==3 )
		gr->GetYaxis()->SetTitleOffset(1.5);
	else
		gr->GetYaxis()->SetTitleOffset(1.3);
	if ( pmt!=3 )
		gr->GetYaxis()->SetRangeUser(-1,1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.8);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubOffsetCh"+pmtId+stat+".pdf");
*/



	// ====================
	// *** For Baseline ***

/*
	TH1F *histBlHbase = new TH1F("histBlHbase", "UUB Mean Baseline HBase "+pmtId+" "+stat, 100, 200, 300);
	TH1F *histBlCalib = new TH1F("histBlCalib", "UUB Baseline Calib "+pmtId+" "+stat, 100, 200, 300);

	TF1 *linef = new TF1("linef", "[1]*x + [0]");

	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = blHbase;
		histBlHbase->Fill( blHbase );
		if ( blHbase < minlimHb )
			minlimHb = blHbase;
		if ( blHbase > maxlimHb )
			maxlimHb = blHbase;
	}
	minlimHb -= 1.;
	maxlimHb += 2.;

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Mean Baseline HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("Mean Baseline HBase / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Fit(linef, "Q");

	p0nm.Form("%.2f", linef->GetParameters()[0]);
	p1nm.Form("%.2e", linef->GetParameters()[1]);
	ptstats = new TPaveStats(0.74,0.775,0.98,0.935,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);	
  ptstats_LaTex = ptstats->AddText("p1*x + p0");
  ptstats_LaTex = ptstats->AddText("p0 = "+p0nm);
  ptstats_LaTex = ptstats->AddText("p1 = "+p1nm);
  ptstats->SetOptStat(0);
  ptstats->Draw();
  gr->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(gr->GetListOfFunctions());

	gr->Draw("AP");
	c1->Print("../plots/uubBlHbase"+pmtId+stat+".pdf");


	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );		
		xbl[tmp] = evtTime;
		ybl[tmp] = blCalib;
		histBlCalib->Fill( blCalib );
		if ( blCalib < minlimCl )
			minlimCl = blCalib;
		if ( blCalib > maxlimCl )
			maxlimCl = blCalib;
	}

	minlimCl -= 1.;
	maxlimCl += 2.;

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Baseline Calib "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("Calib.Base / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Fit(linef, "Q");

	p0nm.Form("%.2f", linef->GetParameters()[0]);
	p1nm.Form("%.2e", linef->GetParameters()[1]);
	ptstats = new TPaveStats(0.74,0.775,0.98,0.935,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);	
  ptstats_LaTex = ptstats->AddText("p1*x + p0");
  ptstats_LaTex = ptstats->AddText("p0 = "+p0nm);
  ptstats_LaTex = ptstats->AddText("p1 = "+p1nm);
  ptstats->SetOptStat(0);
  ptstats->Draw();
  gr->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(gr->GetListOfFunctions());

	gr->Draw("AP");
	c1->Print("../plots/ubBlCalib"+pmtId+stat+".pdf");

	// ===========================
	// *** For Baseline histos ***

	c0->cd();
	histBlHbase->SetLineWidth(2.0);
	histBlHbase->GetXaxis()->SetRangeUser((int)minlimHb, (int)maxlimHb);	
	histBlHbase->SetStats(kTRUE);
	histBlHbase->GetYaxis()->SetTitle("Counts / au");
	histBlHbase->GetXaxis()->SetTitle("1 / FADC");
	histBlHbase->Draw();
	c0->Print("../plots/ubBlHistHbase"+pmtId+stat+".pdf");

	c0->cd();
	histBlCalib->SetLineWidth(2.0);
	histBlCalib->GetXaxis()->SetRangeUser((int)minlimCl, (int)maxlimCl);
	histBlCalib->SetStats(kTRUE);
	histBlCalib->GetYaxis()->SetTitle("Counts / au");
	histBlCalib->GetXaxis()->SetTitle("1 / FADC");
	histBlCalib->Draw();
	c0->Print("../plots/ubBlHistCalib"+pmtId+stat+".pdf");
*/


	// =================================
	// *** For Corrected Peak Histos ***
/*
	histos->GetEntry(0); // Set for first event
	int cnt = 0;

	for (int i=0; i<peak->GetXaxis()->GetNbins(); i++)
		cnt += peak->GetBinContent(i);

	p0nm.Form("%.2d",cnt);

	cerr << blHbase << " " << blCalib << endl;
	cerr << peak->GetBinCenter(1) << endl;
	*/

/*
	c1->cd();
	peak->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 283");
	peak->SetStats(0);
	peak->SetLineColor(kBlue);
	peak->GetXaxis()->SetTitle("1. / FADC");
	peak->GetYaxis()->SetTitle("Counts / au");
	peak->GetYaxis()->SetTitleOffset(1.3);
	peak->Draw();
	c1->Print("../plots/ubPeakCoorHBase"+pmtId+stat+".pdf");
*/
/*
	c1->cd();	
	peak->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 283");
	peak->SetStats(0);
	peak->SetLineColor(kBlue);
	peak->GetXaxis()->SetTitle("1. / FADC");
	peak->GetXaxis()->SetRangeUser(-60, 60);
	peak->GetYaxis()->SetTitle("Counts / au");
	peak->GetYaxis()->SetTitleOffset(1.3);
	peak->Draw();
	c1->Print("../plots/ubPeakCoorHBaseZoom"+pmtId+stat+".pdf");
*/

	// =====================================
	// *** First Bin for Peak Histograms ***
/*
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = peak->GetBinCenter(1);
	}
	*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Peak histogram First Bin Center HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/ubPeakFirstBinHBase"+pmtId+stat+".pdf");
*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Peak histogram First Bin Center Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/ubPeakFirstBinCalibBase"+pmtId+stat+".pdf");

*/
















	// =============================
	// *** For Charge Histograms ***
/*
	histos->GetEntry(0);
	c1->cd();
	charge->SetTitle("UUB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 283");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../plots/ubChargeCoorHBase"+pmtId+stat+".pdf");

	c1->cd();	
	charge->SetTitle("UUB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 283");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetXaxis()->SetRangeUser(charge->GetXaxis()->GetBinLowEdge(1), charge->GetXaxis()->GetBinLowEdge(30));
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../plots/ubChargeCoorHBaseZoom"+pmtId+stat+".pdf");
*/

	/*
	c1->cd();
	charge->SetTitle("UUB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 233");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../plots/ubChargeCoorCalibBase"+pmtId+stat+".pdf");

	c1->cd();	
	charge->SetTitle("UUB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 233");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetXaxis()->SetRangeUser(charge->GetXaxis()->GetBinLowEdge(1), charge->GetXaxis()->GetBinLowEdge(30));
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../plots/ubChargeCoorCalibBaseZoom"+pmtId+stat+".pdf");
*/


	// =====================================
	// *** First Bin for Peak Histograms ***
/*
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = peak->GetBinCenter(1);
	}
	*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Peak histogram First Bin Center HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/ubPeakFirstBinHBase"+pmtId+stat+".pdf");
*/


	// =======================================
	// *** First Bin for Charge Histograms ***
/*
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = charge->GetBinCenter(1);
	}
	*/

/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Charge histogram First Bin Center HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/ubChargeFirstBinHBase"+pmtId+stat+".pdf");
	*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Charge histogram First Bin Center Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/ubChargeFirstBinCalibBase"+pmtId+stat+".pdf");
	*/


	// =============================
	// *** Area/Peak Calculation ***
/*
	nx = histos->GetEntries();
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = chpk;
	}
	*/

	/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB A/P "+pmtId+" "+stat+" HBase-Offset");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / 8.33*ns");
	gr->GetYaxis()->SetRangeUser(0, 12);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/ubApOffset"+pmtId+stat+".pdf");
	*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB A/P "+pmtId+" "+stat+" Calib.Base-Offset");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / 8.33*ns");
	gr->GetYaxis()->SetRangeUser(0, 12);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	//c1->Print("../plots/ubApOffsetCalib"+pmtId+stat+".pdf");
*/

/*
	c1->cd();
	peak->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 233");
	peak->SetStats(0);
	peak->SetLineColor(kBlue);
	peak->GetXaxis()->SetTitle("1. / FADC");
	peak->GetYaxis()->SetTitle("Counts / au");
	peak->GetYaxis()->SetTitleOffset(1.3);
	peak->Draw();
	c1->Print("../plots/ubPeakCoorCalibBase"+pmtId+stat+".pdf");

	peak->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 233");
	peak->SetStats(0);
	peak->SetLineColor(kBlue);
	peak->GetXaxis()->SetTitle("1. / FADC");
	peak->GetXaxis()->SetRangeUser(40, 60);
	peak->GetYaxis()->SetTitle("Counts / au");
	peak->GetYaxis()->SetTitleOffset(1.3);
	peak->Draw();
	c1->Print("../plots/ubPeakCoorCalibBaseZoom"+pmtId+stat+".pdf");
*/


/*	
	TLatex latex(9e2,3.4e3,"Entries "+p0nm);
	latex.SetTextSize(0.04);
	latex.DrawClone("Same");
*/

	
/*
	int nx = histos->GetEntries();
	double xbl[nx];
	double ybl[nx];

	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) {
		histos->GetEntry(tmp);
		xbl[tmp] = tmp;
		ybl[tmp] = chpk;
	}

	TGraph* gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB First bin for Peak Histogram BL from Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020 until March 31st, 2021");
	gr->GetYaxis()->SetTitle("1 / FADCp");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d:%H");
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
*/

	/*
	// ================================
	// *** For Baseline Differences ***

	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xbl[tmp] = tmp;
		ybl[tmp] = blHbase - blCalib;
	}

	TF1 *linef = new TF1("linef", "[1]*x + [0]");

	gr = new TGraph(nx,xbl,ybl);
	c0->cd();
	gr->SetTitle("UUB Baseline Difference HBase - Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020 until March 31st, 2021");
	gr->GetYaxis()->SetTitle("(HBase - Calib.Base) / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Fit(linef, "Q");

	TString p0nm;
	TString p1nm;
	p0nm.Form("%f", linef->GetParameters()[0]);
	p1nm.Form("%f", linef->GetParameters()[1]);
	TPaveStats *ptstats = new TPaveStats(0.62,0.775,0.98,0.935,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);
  TText *ptstats_LaTex = ptstats->AddText("p1*x + p0");
  ptstats_LaTex = ptstats->AddText("p0 = "+p0nm);
  ptstats_LaTex = ptstats->AddText("p1 = "+p1nm);
  ptstats->SetOptStat(0);
  ptstats->Draw();
  gr->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(gr->GetListOfFunctions());
	gr->Draw("AP");
	c0->Print("../plots/ubDiffBlHbaseCalib"+pmtId+stat+".pdf");


	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xbl[tmp] = tmp;
		ybl[tmp] = offsetpk - blCalib;
	}

	gr = new TGraph(nx,xbl,ybl);
	c0->cd();
	gr->SetTitle("UUB Offset-Peak - Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020 until March 31st, 2021");
	gr->GetYaxis()->SetTitle("(OffsetPeak - Calib.Base) / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);

	gr->Fit(linef, "Q");
	p0nm.Form("%f", linef->GetParameters()[0]);
	p1nm.Form("%f", linef->GetParameters()[1]);
	ptstats = new TPaveStats(0.62,0.775,0.98,0.935,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);
  ptstats_LaTex = ptstats->AddText("p1*x + p0");
  ptstats_LaTex = ptstats->AddText("p0 = "+p0nm);
  ptstats_LaTex = ptstats->AddText("p1 = "+p1nm);
  ptstats->SetOptStat(0);
  ptstats->Draw();
  gr->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(gr->GetListOfFunctions());

	gr->Draw("AP");
	c0->Print("../plots/ubDiffBlOffsetPkCalib"+pmtId+stat+".pdf");
*/

	// =================================
	// *** Plotting Fitted histograms ***


	//histos->GetEntry(0);
/*
	c1->cd();
	peakCalibBl->SetTitle("UUB Peak histogram "+pmtId+" "+stat+" Corr. for Calib.Base") ;
	peakCalibBl->SetStats(kFALSE);
	peakCalibBl->GetXaxis()->SetTitle("1 / FADCp");
	peakCalibBl->GetYaxis()->SetTitle("Counts / au");
	peakCalibBl->GetYaxis()->SetTitleOffset(1.3);
	peakCalibBl->SetLineColor(kBlue);
	peakCalibBl->Draw();
	c1->Print("../plots/ubPeakCorrCalib"+pmtId+stat+".pdf");
	*/
/*
	c1->cd();
	peakCalibBl->SetTitle("UUB Peak histogram "+pmtId+" "+stat+" Corr. for Calib.Base") ;
	peakCalibBl->SetStats(kFALSE);
	peakCalibBl->GetXaxis()->SetRange(-50,40);
	peakCalibBl->GetXaxis()->SetTitle("1 / FADCp");
	peakCalibBl->GetYaxis()->SetTitle("Counts / au");
	peakCalibBl->GetYaxis()->SetTitleOffset(1.3);
	peakCalibBl->SetLineColor(kBlue);
	peakCalibBl->Draw();
	c1->Print("../plots/ubPeakCorrCalibZoom"+pmtId+stat+".pdf");
*/

	// ==============================================
	// *** Checking for first bin Peak Histograms ***
	
	/*
	int nx = histos->GetEntries();
	double xbl[nx];
	double ybl[nx];

	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) {
		histos->GetEntry(tmp);
		xbl[tmp] = tmp;
		ybl[tmp] = peakCalibBl->GetBinLowEdge(1);
	}

	TGraph* gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB First bin for Peak Histogram BL from Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020 until March 31st, 2021");
	gr->GetYaxis()->SetTitle("1 / FADCp");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	//c1->Print("../plots/uubFirstBinPeakCrrCalib"+pmtId+stat+".pdf");

*/
	//chargeCalibBl->Draw();
	//hcharge->Draw();

	//c1->Print("../plots/uubFitPk"+pmtId+stat+".pdf");

/*
	TCanvas *c2 = new TCanvas("c2","C2",500,10,1200,800);

	histos->GetEntry(0);
	c2->cd();
	//charge->SetLineColor(kBlack);
	//charge->Draw("SAME");
	fittedCh->GetYaxis()->SetTitle("Counts / au");
	fittedCh->GetYaxis()->SetTitleOffset(1.3);
	fittedCh->GetXaxis()->SetTitle("1 / FADCq");
	fittedCh->Draw("AL");
	c2->Print("../plots/uubFitCh"+pmtId+stat+".pdf");

	int nx = histos->GetEntries();
	double xb[nx];
	double yb[nx];

	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xb[tmp] = tmp;
		yb[tmp] = chpk * (8.33/25);
	}

	TGraph* gr = new TGraph(nx,xb,yb);
	c1->cd();
	gr->SetTitle("UUB A/P "+pmtId+" "+stat+" Here A/P*(8.33/25)");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("Events");
	gr->GetYaxis()->SetTitle("Counts / au");
	gr->GetYaxis()->SetRangeUser(0, 5);
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubAp"+pmtId+stat+".pdf");
	*/
}
